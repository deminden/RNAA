//! Enforce the decoded shape before the engine allocates sequence batches.
use crate::fastq::FastqSummary;
use anyhow::Result;
use std::{
    io::{BufRead, Read},
    path::Path,
};

struct Bounded<R> {
    inner: R,
    maximum_bytes: u64,
    maximum_records: u64,
    maximum_record: u64,
    maximum_line: usize,
    maximum_sequence: usize,
    minimum_sequence: usize,
    bytes: u64,
    records: u64,
    record_bytes: u64,
    line: usize,
    line_bytes: usize,
}
pub(super) fn open(
    path: &Path,
    summary: &FastqSummary,
    mate: usize,
) -> Result<impl kallistors::io::ReadSource + kallistors::io::PackedSeqBatchSource + Send> {
    let input = crate::fastq::open_reader(path)?;
    Ok(kallistors::io::FastqReader::new(Bounded {
        inner: input,
        maximum_bytes: if mate == 0 {
            summary.r1_normalized_bytes
        } else {
            summary.r2_normalized_bytes
        },
        maximum_records: summary.fragments,
        maximum_record: summary.max_record_bytes.saturating_add(4),
        maximum_line: summary.max_line_bytes.saturating_add(1),
        maximum_sequence: summary.max_sequence_bytes,
        minimum_sequence: summary.min_sequence_bytes,
        bytes: 0,
        records: 0,
        record_bytes: 0,
        line: 0,
        line_bytes: 0,
    }))
}
impl<R: BufRead> Read for Bounded<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        let input = self.fill_buf()?;
        let n = out.len().min(input.len());
        out[..n].copy_from_slice(&input[..n]);
        self.consume(n);
        Ok(n)
    }
}
impl<R: BufRead> BufRead for Bounded<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        let input = self.inner.fill_buf()?;
        let (mut total, mut record, mut line, mut bytes) =
            (self.bytes, self.record_bytes, self.line, self.line_bytes);
        for (i, &byte) in input.iter().enumerate() {
            let reason = if total >= self.maximum_bytes {
                Some("decoded bytes")
            } else if self.records >= self.maximum_records {
                Some("record count")
            } else if record >= self.maximum_record {
                Some("record bytes")
            } else if bytes >= self.maximum_line {
                Some("line bytes")
            } else if line == 1 && byte != b'\n' && bytes >= self.maximum_sequence {
                Some("sequence length")
            } else if line == 1 && byte == b'\n' && bytes < self.minimum_sequence {
                Some("minimum sequence length")
            } else {
                None
            };
            if let Some(reason) = reason {
                if i > 0 {
                    return Ok(&input[..i]);
                }
                return Err(std::io::Error::other(format!(
                    "quant input exceeds its preflight {reason} promise"
                )));
            }
            total += 1;
            record += 1;
            bytes += 1;
            if byte == b'\n' {
                line += 1;
                bytes = 0;
                if line == 4 {
                    return Ok(&input[..=i]);
                }
            }
        }
        Ok(input)
    }
    fn consume(&mut self, n: usize) {
        if let Ok(bytes) = self.inner.fill_buf() {
            for &byte in bytes.iter().take(n) {
                self.bytes += 1;
                self.record_bytes += 1;
                self.line_bytes += 1;
                if byte == b'\n' {
                    self.line += 1;
                    self.line_bytes = 0;
                    if self.line == 4 {
                        self.line = 0;
                        self.record_bytes = 0;
                        self.records += 1;
                    }
                }
            }
        }
        self.inner.consume(n);
    }
}

#[cfg(test)]
mod tests;
