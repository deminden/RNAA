//! GTF annotation derivation and normalized transcript mapping validation.
use anyhow::{Context, Result, anyhow, bail};
use regex::Regex;
use rnaa_core::model::ReferenceBundle;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub(super) struct Patterns {
    transcript: Regex,
    gene: Regex,
    name: Regex,
    biotype: Regex,
}
impl Patterns {
    pub fn new() -> Result<Self> {
        Ok(Self {
            transcript: Regex::new(r#"transcript_id "([^"]+)""#)?,
            gene: Regex::new(r#"gene_id "([^"]+)""#)?,
            name: Regex::new(r#"gene_name "([^"]+)""#)?,
            biotype: Regex::new(r#"(gene_biotype|gene_type) "([^"]+)""#)?,
        })
    }
    pub fn string_bytes(&self, attributes: &str) -> Result<u64> {
        [
            (&self.transcript, 1),
            (&self.gene, 1),
            (&self.name, 1),
            (&self.biotype, 2),
        ]
        .iter()
        .try_fold(0_u64, |sum, (pattern, group)| {
            let bytes = pattern
                .captures(attributes)
                .and_then(|c| c.get(*group))
                .map_or(0, |m| m.len() as u64);
            sum.checked_add(bytes)
                .context("GTF string accounting overflows")
        })
    }
}

pub(super) fn validate_index_mapping(bundle: &ReferenceBundle) -> Result<serde_json::Value> {
    let index = kallistors::index::Index::load(&bundle.kallisto_index_path)
        .map_err(|err| anyhow!("new reference index cannot be loaded: {err}"))?;
    let mut table = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&bundle.tx2gene_path)?;
    let mapped = table
        .records()
        .map(|row| Ok(row?.get(0).unwrap_or_default().to_string()))
        .collect::<Result<BTreeSet<_>>>()?;
    let mut unique = BTreeSet::new();
    let mut missing = Vec::new();
    let mut matched = 0;
    for target in &index.transcripts {
        let id = normalize_stable_id(&target.name);
        if id.is_empty() || !unique.insert(id.clone()) {
            bail!("ambiguous normalized reference transcript ID: {id}");
        }
        if mapped.contains(&id) {
            matched += 1;
        } else {
            missing.push(id);
        }
    }
    if matched == 0 {
        bail!("reference transcriptome and GTF have no matching transcript IDs");
    }
    Ok(
        json!({"indexed_transcripts":index.transcripts.len(),"mapped_transcripts":matched,"unmapped_transcripts":missing.len(),"unmapped_examples":missing.into_iter().take(20).collect::<Vec<_>>()}),
    )
}

pub(super) fn build_annotation_tables(
    gtf_path: &Path,
    tx2gene_path: &Path,
    gene_annotation_path: &Path,
    expected: &super::planning::GtfPlan,
    check: &impl Fn() -> Result<()>,
) -> Result<()> {
    let mut reader = super::planning::Lines::new(
        gtf_path,
        super::planning::Limits {
            line: expected.max_line,
            decoded: Some(expected.decoded_bytes),
        },
        check,
    )?;
    let mut raw = Vec::new();
    let patterns = Patterns::new()?;

    let mut tx2gene: BTreeMap<String, String> = BTreeMap::new();
    let mut gene_annotation: BTreeMap<String, (String, String)> = BTreeMap::new();

    let (mut rows, mut strings) = (0_u64, 0_u64);
    while reader.next(&mut raw)? {
        let number = reader.rows - 1;
        let line = std::str::from_utf8(&raw)?.trim_end_matches(['\r', '\n']);
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        if line.split('\t').count() != 9 {
            bail!(
                "GTF line {} must have exactly nine tab-separated columns",
                number + 1
            );
        }
        let mut fields = line.split('\t');
        let _seqname = fields.next();
        let _source = fields.next();
        let feature = fields.next().unwrap_or_default();
        let _start = fields.next();
        let _end = fields.next();
        let _score = fields.next();
        let _strand = fields.next();
        let _frame = fields.next();
        let attributes = fields.next().unwrap_or_default();
        rows += 1;
        strings = strings
            .checked_add(patterns.string_bytes(attributes)?)
            .context("GTF string count overflows")?;
        if rows > expected.rows || strings > expected.string_bytes {
            bail!("GTF exceeds its planned row/string promise");
        }

        let Some(gene_id) = patterns.gene.captures(attributes).and_then(|captures| {
            captures
                .get(1)
                .map(|value| normalize_stable_id(value.as_str()))
        }) else {
            bail!("GTF line {} is missing gene_id", number + 1);
        };
        if gene_id.is_empty() {
            bail!("empty normalized gene ID at GTF line {}", number + 1);
        }
        if let Some(transcript_id) = patterns
            .transcript
            .captures(attributes)
            .and_then(|captures| {
                captures
                    .get(1)
                    .map(|value| normalize_stable_id(value.as_str()))
            })
        {
            if transcript_id.is_empty() {
                bail!("empty normalized transcript ID at GTF line {}", number + 1);
            }
            if let Some(previous) = tx2gene.insert(transcript_id.clone(), gene_id.clone())
                && previous != gene_id
            {
                bail!("conflicting transcript-to-gene mapping for {transcript_id}");
            }
        }

        if feature == "gene" || !gene_annotation.contains_key(&gene_id) {
            let gene_name = patterns
                .name
                .captures(attributes)
                .and_then(|captures| captures.get(1).map(|value| value.as_str().to_string()))
                .unwrap_or_default();
            let biotype = patterns
                .biotype
                .captures(attributes)
                .and_then(|captures| captures.get(2).map(|value| value.as_str().to_string()))
                .unwrap_or_default();
            gene_annotation.insert(gene_id, (gene_name, biotype));
        }
    }

    if reader.bytes != expected.decoded_bytes || reader.digest() != expected.decoded_sha256 {
        bail!("GTF changed after reference planning");
    }
    check()?;

    if tx2gene.is_empty() || gene_annotation.is_empty() {
        bail!("GTF contains no transcript-to-gene mappings");
    }
    let mut tx_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(tx2gene_path)
        .with_context(|| format!("failed to create {}", tx2gene_path.display()))?;
    tx_writer.write_record(["transcript_id", "gene_id"])?;
    let mut written = 0_u64;
    let mut next_check = 0_u64;
    for (transcript_id, gene_id) in tx2gene {
        written = written.saturating_add((transcript_id.len() + gene_id.len()) as u64);
        if written >= next_check {
            check()?;
            next_check = written.saturating_add(8 * 1024 * 1024);
        }
        tx_writer.write_record([transcript_id, gene_id])?;
    }
    tx_writer.flush()?;

    let mut gene_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(gene_annotation_path)
        .with_context(|| format!("failed to create {}", gene_annotation_path.display()))?;
    gene_writer.write_record(["gene_id", "gene_name", "gene_biotype"])?;
    for (gene_id, (gene_name, biotype)) in gene_annotation {
        written = written.saturating_add((gene_id.len() + gene_name.len() + biotype.len()) as u64);
        if written >= next_check {
            check()?;
            next_check = written.saturating_add(8 * 1024 * 1024);
        }
        gene_writer.write_record([gene_id, gene_name, biotype])?;
    }
    gene_writer.flush()?;

    Ok(())
}

pub(super) fn normalize_stable_id(value: &str) -> String {
    rnaa_core::util::normalize_feature_id(value)
}
