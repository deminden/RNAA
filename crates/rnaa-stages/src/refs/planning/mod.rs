//! Bounded decoded-source observations and structural build/publication promises.
mod resources;
#[cfg(test)]
mod tests;
mod unique;
use super::annotation::Patterns;
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use rnaa_core::{
    ManifestArtifact,
    config::RefsConfig,
    resources::{DiskRequest, ExecutionResources, ResourceRequest},
};
use sha2::{Digest, Sha256};
use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};
use unique::UniqueWords;
const MIB: u64 = 1024 * 1024;

#[derive(Clone, Copy)]
pub(super) struct Limits {
    pub line: usize,
    pub decoded: Option<u64>,
}
impl TryFrom<&RefsConfig> for Limits {
    type Error = anyhow::Error;
    fn try_from(config: &RefsConfig) -> Result<Self> {
        config.validate()?;
        Ok(Self {
            line: usize::try_from(config.max_source_line_mb * MIB)?,
            decoded: if config.max_decoded_mb == 0 {
                None
            } else {
                Some(config.max_decoded_mb * MIB)
            },
        })
    }
}
pub(super) fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    Ok(if super::source::is_gzip(path)? {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}
pub(super) struct Lines<'a> {
    reader: Box<dyn BufRead>,
    pub bytes: u64,
    pub rows: u64,
    pub max_line: usize,
    limits: Limits,
    hash: Sha256,
    next_check: u64,
    check: &'a dyn Fn() -> Result<()>,
}
impl<'a> Lines<'a> {
    pub fn new(path: &Path, limits: Limits, check: &'a dyn Fn() -> Result<()>) -> Result<Self> {
        Ok(Self {
            reader: open(path)?,
            bytes: 0,
            rows: 0,
            max_line: 0,
            limits,
            hash: Sha256::new(),
            next_check: 0,
            check,
        })
    }
    pub fn next(&mut self, line: &mut Vec<u8>) -> Result<bool> {
        if self.bytes >= self.next_check {
            (self.check)()?;
            self.next_check = self.bytes.saturating_add(8 * MIB);
        }
        line.clear();
        let allowed = (self.limits.line as u64).min(
            self.limits
                .decoded
                .map_or(u64::MAX, |n| n.saturating_sub(self.bytes)),
        );
        let count = (&mut self.reader)
            .take(allowed + 1)
            .read_until(b'\n', line)?;
        if count == 0 {
            return Ok(false);
        }
        self.bytes = self
            .bytes
            .checked_add(count as u64)
            .context("reference decoded byte count overflows")?;
        if count > self.limits.line || self.limits.decoded.is_some_and(|n| self.bytes > n) {
            bail!("reference source exceeds its decoded line/group byte limit");
        }
        self.rows = self
            .rows
            .checked_add(1)
            .context("reference row counter overflows")?;
        self.max_line = self.max_line.max(count);
        self.hash.update(&*line);
        Ok(true)
    }
    pub fn digest(&self) -> String {
        hex::encode(self.hash.clone().finalize())
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize)]
pub(super) struct FastaPlan {
    pub decoded_bytes: u64,
    pub decoded_sha256: String,
    pub max_line: usize,
    pub transcripts: u64,
    pub bases: u64,
    pub max_sequence: u64,
    pub name_bytes: u64,
    pub kmer_positions: u64,
    pub clean_unique_kmers: u64,
    pub uncertain_kmers: u64,
    pub unique_kmers_upper: u64,
}
#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize)]
pub(super) struct GtfPlan {
    pub decoded_bytes: u64,
    pub decoded_sha256: String,
    pub max_line: usize,
    pub rows: u64,
    pub string_bytes: u64,
}
#[derive(Clone, Debug, serde::Serialize)]
pub(super) struct Plan {
    pub model: &'static str,
    pub fasta: FastaPlan,
    pub gtf: GtfPlan,
    pub memory_bytes: u64,
    pub copy_bytes: u64,
    pub staging_bytes: u64,
    pub index_bytes: u64,
    pub table_bytes: u64,
    pub metadata_bytes: u64,
    pub threads: usize,
}
fn add(value: &mut u64, n: u64) -> Result<()> {
    *value = value
        .checked_add(n)
        .context("reference observation overflows")?;
    Ok(())
}
fn scan_fasta(
    path: &Path,
    limits: Limits,
    k: usize,
    words: &mut UniqueWords<'_>,
    check: &dyn Fn() -> Result<()>,
) -> Result<FastaPlan> {
    let mut reader = Lines::new(path, limits, check)?;
    let mut line = Vec::new();
    let mut plan = FastaPlan::default();
    let (mut length, mut clean, mut fwd, mut rev) = (0_u64, 0_usize, 0_u64, 0_u64);
    let mask = (1_u64 << (2 * k)) - 1;
    while reader.next(&mut line)? {
        let text = std::str::from_utf8(&line)?.trim_end_matches(['\n', '\r']);
        if let Some(name) = text.strip_prefix('>') {
            plan.max_sequence = plan.max_sequence.max(length);
            length = 0;
            clean = 0;
            fwd = 0;
            rev = 0;
            add(&mut plan.transcripts, 1)?;
            let name = name.split_whitespace().next().unwrap_or("");
            add(
                &mut plan.name_bytes,
                if name.is_empty() {
                    20
                } else {
                    name.len() as u64
                },
            )?;
        } else {
            if plan.transcripts == 0 && !text.is_empty() {
                bail!("reference FASTA sequence precedes its first header");
            }
            if !text.is_ascii() {
                bail!("reference FASTA sequence must contain ASCII symbols");
            }
            for byte in text.bytes() {
                length = length
                    .checked_add(1)
                    .context("reference sequence length overflows")?;
                add(&mut plan.bases, 1)?;
                let code = match byte.to_ascii_uppercase() {
                    b'A' => Some(0),
                    b'C' => Some(1),
                    b'G' => Some(2),
                    b'T' | b'U' => Some(3),
                    _ => None,
                };
                if let Some(code) = code {
                    fwd = ((fwd << 2) | code) & mask;
                    rev = (rev >> 2) | ((3 - code) << (2 * (k - 1)));
                    clean = clean.saturating_add(1);
                } else {
                    clean = 0;
                    fwd = 0;
                    rev = 0;
                }
                if length >= k as u64 {
                    add(&mut plan.kmer_positions, 1)?;
                    if clean >= k {
                        words.push(fwd.min(rev))?;
                    } else {
                        add(&mut plan.uncertain_kmers, 1)?;
                    }
                }
            }
            if length > i32::MAX as u64 {
                bail!("reference transcript exceeds the engine's positional width");
            }
        }
    }
    plan.max_sequence = plan.max_sequence.max(length);
    if plan.transcripts == 0 || plan.bases == 0 || plan.transcripts > i32::MAX as u64 {
        bail!("reference FASTA has no usable transcript collection");
    }
    plan.decoded_bytes = reader.bytes;
    plan.max_line = reader.max_line;
    plan.decoded_sha256 = reader.digest();
    Ok(plan)
}
pub(super) fn scan_gtf(
    path: &Path,
    limits: Limits,
    check: &dyn Fn() -> Result<()>,
) -> Result<GtfPlan> {
    let patterns = Patterns::new()?;
    let mut reader = Lines::new(path, limits, check)?;
    let mut line = Vec::new();
    let mut result = GtfPlan::default();
    while reader.next(&mut line)? {
        let text = std::str::from_utf8(&line)?.trim_end_matches(['\n', '\r']);
        if text.trim().is_empty() || text.starts_with('#') {
            continue;
        }
        if text.split('\t').count() != 9 {
            bail!(
                "GTF line {} must have exactly nine tab-separated columns",
                reader.rows
            );
        }
        add(&mut result.rows, 1)?;
        add(
            &mut result.string_bytes,
            patterns.string_bytes(text.split('\t').nth(8).unwrap_or_default())?,
        )?;
    }
    if result.rows == 0 {
        bail!("reference GTF has no annotation rows");
    }
    result.decoded_bytes = reader.bytes;
    result.decoded_sha256 = reader.digest();
    result.max_line = reader.max_line;
    Ok(result)
}
