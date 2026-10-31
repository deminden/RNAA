use super::{FastqLimits, FastqSummary};
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
    time::{Duration, Instant},
};

pub(super) struct ScanControl<'a, F> {
    check: &'a F,
    next_bytes: u64,
    next_records: u64,
    last: Instant,
}
impl<'a, F: Fn() -> Result<()>> ScanControl<'a, F> {
    pub fn new(check: &'a F) -> Self {
        Self {
            check,
            next_bytes: 0,
            next_records: 0,
            last: Instant::now(),
        }
    }
    fn fence(&mut self, summary: &FastqSummary) -> Result<()> {
        if summary.decoded_bytes >= self.next_bytes
            || summary.read_records >= self.next_records
            || (summary.read_records.is_multiple_of(256)
                && self.last.elapsed() >= Duration::from_millis(100))
        {
            (self.check)()?;
            self.next_bytes = summary.decoded_bytes.saturating_add(8 * 1024 * 1024);
            self.next_records = summary.read_records.saturating_add(4096);
            self.last = Instant::now();
        }
        Ok(())
    }
}

pub(crate) fn open(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let mut file = BufReader::new(File::open(path)?);
    if file.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(file))
    }
}
pub(super) struct Record {
    pub lines: [Vec<u8>; 4],
    pub normalized_bytes: u64,
}
pub(super) fn record(
    reader: &mut dyn BufRead,
    limits: FastqLimits,
    summary: &mut FastqSummary,
    control: &mut ScanControl<'_, impl Fn() -> Result<()>>,
) -> Result<Option<Record>> {
    let mut lines: [Vec<u8>; 4] = Default::default();
    let mut wire = 0_u64;
    let mut maximum = 0;
    for (index, line) in lines.iter_mut().enumerate() {
        control.fence(summary)?;
        let decoded_remaining = limits
            .max_decoded_bytes
            .map(|limit| limit.saturating_sub(summary.decoded_bytes));
        let allowed = (limits.max_line_bytes as u64).min(decoded_remaining.unwrap_or(u64::MAX));
        let read = (&mut *reader).take(allowed + 1).read_until(b'\n', line)?;
        if read == 0 {
            if index == 0 {
                return Ok(None);
            }
            bail!("truncated FASTQ record");
        }
        summary.decoded_bytes = summary
            .decoded_bytes
            .checked_add(read as u64)
            .context("decoded FASTQ byte count overflows")?;
        if limits
            .max_decoded_bytes
            .is_some_and(|limit| summary.decoded_bytes > limit)
        {
            bail!(
                "decoded FASTQ group exceeds its byte limit {}",
                limits.max_decoded_bytes.unwrap()
            );
        }
        if read > limits.max_line_bytes {
            bail!(
                "FASTQ line exceeds max_line_bytes={}",
                limits.max_line_bytes
            );
        }
        wire += read as u64;
        maximum = maximum.max(read);
        if line.last() == Some(&b'\n') {
            line.pop();
        } else {
            summary.needs_staging = true;
        }
        if line.last() == Some(&b'\r') {
            line.pop();
            summary.needs_staging = true;
        }
        line.push(b'\n');
    }
    let seq = &lines[1][..lines[1].len() - 1];
    let quality = &lines[3][..lines[3].len() - 1];
    if !lines[0].starts_with(b"@")
        || !lines[2].starts_with(b"+")
        || seq.is_empty()
        || seq.len() != quality.len()
        || !seq.is_ascii()
        || quality.iter().any(|byte| !(33..=126).contains(byte))
    {
        bail!("malformed FASTQ record");
    }
    // Headers and plus lines keep the public UTF-8 input contract.
    for line in [&lines[0], &lines[2]] {
        std::str::from_utf8(line).context("FASTQ record is not UTF-8")?;
    }
    let normalized = lines.iter().try_fold(0_u64, |total, line| {
        total
            .checked_add(line.len() as u64)
            .context("normalized FASTQ byte count overflows")
    })?;
    summary.max_sequence_bytes = summary.max_sequence_bytes.max(seq.len());
    summary.min_sequence_bytes = if summary.min_sequence_bytes == 0 {
        seq.len()
    } else {
        summary.min_sequence_bytes.min(seq.len())
    };
    summary.max_record_bytes = summary.max_record_bytes.max(wire);
    summary.max_line_bytes = summary.max_line_bytes.max(maximum);
    summary.min_record_bytes = summary.min_record_bytes.min(wire);
    Ok(Some(Record {
        lines,
        normalized_bytes: normalized,
    }))
}
pub(super) fn read_id(header: &[u8]) -> &[u8] {
    let id = header
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    id.strip_suffix(b"/1")
        .or_else(|| id.strip_suffix(b"/2"))
        .unwrap_or(id)
}
