use super::{NativeProcessOptions, window::BatchWindow};
use crate::pipeline::{Batch, PairedBatch};
use anyhow::{Context, Result, bail};
use crossbeam_channel::Sender;
use std::{
    io::{BufRead, Read},
    sync::Arc,
};

/// Append one normalized complete record without an unbounded read_line.
fn record(
    reader: &mut dyn BufRead,
    buffer: &mut Vec<u8>,
    limit: usize,
) -> Result<Option<[usize; 4]>> {
    let mut starts = [0; 4];
    let mut wire_bytes = 0_usize;
    for (line, start) in starts.iter_mut().enumerate() {
        *start = buffer.len();
        let read = (&mut *reader)
            .take((limit - wire_bytes) as u64 + 1)
            .read_until(b'\n', buffer)?;
        if read == 0 {
            if line == 0 {
                return Ok(None);
            }
            bail!("truncated native FASTQ record");
        }
        wire_bytes = wire_bytes
            .checked_add(read)
            .context("FASTQ record length overflow")?;
        if wire_bytes > limit {
            bail!("native FASTQ record exceeds max_record_bytes={limit}");
        }
        if buffer.last() == Some(&b'\n') {
            buffer.pop();
        }
        if buffer.last() == Some(&b'\r') {
            buffer.pop();
        }
        buffer.push(b'\n');
    }
    let seq = &buffer[starts[1]..starts[2] - 1];
    let qual = &buffer[starts[3]..buffer.len() - 1];
    if buffer[starts[0]] != b'@'
        || buffer[starts[2]] != b'+'
        || seq.is_empty()
        || seq.len() != qual.len()
        || !seq.is_ascii()
        || qual.iter().any(|byte| !(33..=126).contains(byte))
    {
        bail!("malformed native FASTQ record");
    }
    Ok(Some(starts))
}
fn read_id(header: &[u8]) -> &[u8] {
    let id = header
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    id.strip_suffix(b"/1")
        .or_else(|| id.strip_suffix(b"/2"))
        .unwrap_or(id)
}

pub(super) fn single(
    input: &str,
    sender: Sender<Option<Batch>>,
    options: NativeProcessOptions,
    window: &Arc<BatchWindow>,
) -> Result<()> {
    let mut reader = self::input(input, options, 0)?;
    let mut id = 0_u64;
    loop {
        let permit = window.acquire()?;
        let mut bytes = Vec::new();
        let mut records = Vec::new();
        let mut eof = false;
        while bytes.len() < options.batch_bytes {
            match record(&mut reader, &mut bytes, options.max_record_bytes)? {
                Some(record) => records.push(record),
                None => {
                    eof = true;
                    break;
                }
            }
        }
        if !records.is_empty() {
            sender
                .send(Some(Batch {
                    id,
                    buf: Arc::new(bytes),
                    recs: records,
                    permit: Some(permit),
                }))
                .map_err(|_| anyhow::anyhow!("native worker input closed"))?;
            id = id
                .checked_add(1)
                .context("native batch identity overflow")?;
        }
        if eof {
            break;
        }
    }
    let _ = sender.send(None);
    Ok(())
}

pub(super) fn paired(
    input1: &str,
    input2: &str,
    sender: Sender<Option<PairedBatch>>,
    options: NativeProcessOptions,
    window: &Arc<BatchWindow>,
) -> Result<()> {
    let mut left = input(input1, options, 0)?;
    let mut right = input(input2, options, 1)?;
    let mut id = 0_u64;
    loop {
        let permit = window.acquire()?;
        let (mut a, mut b, mut ra, mut rb) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut eof = false;
        while a.len().max(b.len()) < options.batch_bytes {
            match (
                record(&mut left, &mut a, options.max_record_bytes)?,
                record(&mut right, &mut b, options.max_record_bytes)?,
            ) {
                (Some(r1), Some(r2)) => {
                    if read_id(&a[r1[0]..r1[1] - 1]) != read_id(&b[r2[0]..r2[1] - 1]) {
                        bail!("native paired FASTQ read identifiers differ");
                    }
                    ra.push(r1);
                    rb.push(r2);
                }
                (None, None) => {
                    eof = true;
                    break;
                }
                _ => bail!("native paired FASTQs contain unequal record counts"),
            }
        }
        if !ra.is_empty() {
            sender
                .send(Some(PairedBatch {
                    id,
                    buf_r1: Arc::new(a),
                    buf_r2: Arc::new(b),
                    recs_r1: ra,
                    recs_r2: rb,
                    permit: Some(permit),
                }))
                .map_err(|_| anyhow::anyhow!("native paired worker input closed"))?;
            id = id
                .checked_add(1)
                .context("native batch identity overflow")?;
        }
        if eof {
            break;
        }
    }
    let _ = sender.send(None);
    Ok(())
}

/// Bound the existing serial reader before its record parser allocates strings.
pub(super) fn input(
    path: &str,
    options: NativeProcessOptions,
    mate: usize,
) -> Result<RecordCap<InputCap<Box<dyn BufRead + Send>>>> {
    Ok(RecordCap::new(
        InputCap {
            inner: crate::io::open_input(path)?,
            remaining: options.input_byte_limits.map(|v| v[mate]),
        },
        options,
    ))
}
pub(super) struct InputCap<R> {
    inner: R,
    remaining: Option<u64>,
}
impl<R: BufRead> Read for InputCap<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        let bytes = self.fill_buf()?;
        let count = bytes.len().min(out.len());
        out[..count].copy_from_slice(&bytes[..count]);
        self.consume(count);
        Ok(count)
    }
}
impl<R: BufRead> BufRead for InputCap<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        let bytes = self.inner.fill_buf()?;
        if let Some(remaining) = self.remaining {
            if remaining == 0 && !bytes.is_empty() {
                return Err(std::io::Error::other(
                    "native input exceeds its preflight decoded byte promise",
                ));
            }
            return Ok(&bytes[..bytes
                .len()
                .min(usize::try_from(remaining).unwrap_or(usize::MAX))]);
        }
        Ok(bytes)
    }
    fn consume(&mut self, amount: usize) {
        if let Some(remaining) = &mut self.remaining {
            *remaining = remaining.saturating_sub(amount as u64);
        }
        self.inner.consume(amount);
    }
}

pub(super) struct RecordCap<R> {
    inner: R,
    limit: usize,
    bytes: usize,
    lines: usize,
    records: u64,
    minimum: usize,
    maximum_records: Option<u64>,
}
impl<R> RecordCap<R> {
    pub fn new(inner: R, options: NativeProcessOptions) -> Self {
        Self {
            inner,
            limit: options.max_record_bytes,
            bytes: 0,
            lines: 0,
            records: 0,
            minimum: options.min_record_bytes,
            maximum_records: options.max_input_records,
        }
    }
}
impl<R: BufRead> Read for RecordCap<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        let bytes = self.fill_buf()?;
        let count = bytes.len().min(out.len());
        out[..count].copy_from_slice(&bytes[..count]);
        self.consume(count);
        Ok(count)
    }
}
impl<R: BufRead> BufRead for RecordCap<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        let available = self.inner.fill_buf()?;
        let mut bytes = self.bytes;
        let mut lines = self.lines;
        let records = self.records;
        for (index, &byte) in available.iter().enumerate() {
            let shape_error = if lines == 0
                && bytes == 0
                && self.maximum_records.is_some_and(|max| records >= max)
            {
                Some("native input exceeds its preflight record count")
            } else if byte == b'\n' && lines == 3 && bytes + 1 < self.minimum {
                Some("native input record is shorter than its preflight minimum")
            } else {
                None
            };
            if let Some(message) = shape_error {
                if index == 0 {
                    return Err(std::io::Error::other(message));
                }
                return Ok(&available[..index]);
            }
            if bytes == self.limit {
                if index == 0 {
                    return Err(std::io::Error::other(format!(
                        "native FASTQ record exceeds max_record_bytes={}",
                        self.limit
                    )));
                }
                return Ok(&available[..index]);
            }
            bytes += 1;
            if byte == b'\n' {
                lines += 1;
                if lines == 4 {
                    // Return at this record boundary. Rescanning the whole
                    // decoder buffer for each read_until would be quadratic in
                    // its size for short reads, especially on the parallel path.
                    return Ok(&available[..=index]);
                }
            }
        }
        Ok(available)
    }
    fn consume(&mut self, amount: usize) {
        if let Ok(bytes) = self.inner.fill_buf() {
            for &byte in bytes.iter().take(amount) {
                self.bytes += 1;
                if byte == b'\n' {
                    self.lines += 1;
                    if self.lines == 4 {
                        self.lines = 0;
                        self.bytes = 0;
                        self.records = self.records.saturating_add(1);
                    }
                }
            }
        }
        self.inner.consume(amount);
    }
}
