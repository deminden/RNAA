use anyhow::{Context, Result, bail};
use std::{
    fs::File,
    io::{BufReader, Read, Seek, SeekFrom},
    path::Path,
};

/// Fixed-buffer reader with nested section bounds. Serialized counts never size
/// an allocation; a count must first fit its enclosing physical byte range.
pub(super) struct Reader<'a> {
    input: BufReader<File>,
    pub pos: u64,
    pub end: u64,
    check: &'a dyn Fn() -> Result<()>,
    calls: usize,
    checked_at: u64,
}
impl<'a> Reader<'a> {
    pub fn open(path: &Path, check: &'a dyn Fn() -> Result<()>) -> Result<Self> {
        check()?;
        let file = File::open(path)?;
        let metadata = file.metadata()?;
        if !metadata.is_file() {
            bail!("index preflight requires a regular file");
        }
        Ok(Self {
            input: BufReader::with_capacity(64 * 1024, file),
            pos: 0,
            end: metadata.len(),
            check,
            calls: 0,
            checked_at: 0,
        })
    }
    fn poll(&mut self) -> Result<()> {
        self.calls += 1;
        if self.calls.is_multiple_of(4096)
            || self.pos.saturating_sub(self.checked_at) >= 8 * 1024 * 1024
        {
            (self.check)()?;
            self.checked_at = self.pos;
        }
        Ok(())
    }
    pub fn fits(&self, count: u64, width: u64) -> Result<u64> {
        let bytes = count
            .checked_mul(width)
            .context("index preflight count overflows")?;
        if bytes > self.end.saturating_sub(self.pos) {
            bail!(
                "index preflight count exceeds its section at byte {}",
                self.pos
            );
        }
        Ok(bytes)
    }
    pub fn skip(&mut self, bytes: u64) -> Result<()> {
        self.fits(bytes, 1)?;
        self.pos += bytes;
        self.input.seek(SeekFrom::Start(self.pos))?;
        self.poll()
    }
    pub fn seek(&mut self, position: u64) -> Result<()> {
        if position > self.end {
            bail!("index preflight seek exceeds its section");
        }
        self.input.seek(SeekFrom::Start(position))?;
        self.pos = position;
        self.poll()
    }
    pub fn bytes<const N: usize>(&mut self) -> Result<[u8; N]> {
        self.fits(N as u64, 1)?;
        let mut bytes = [0; N];
        self.input.read_exact(&mut bytes)?;
        self.pos += N as u64;
        self.poll()?;
        Ok(bytes)
    }
    pub fn u8(&mut self) -> Result<u8> {
        Ok(self.bytes::<1>()?[0])
    }
    pub fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.bytes()?))
    }
    pub fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.bytes()?))
    }
    pub fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.bytes()?))
    }
    pub fn contains(&mut self, bytes: u64, needle: &[u8]) -> Result<bool> {
        self.fits(bytes, 1)?;
        let mut matched = 0;
        for read in 0..bytes {
            let byte = self.u8()?;
            matched = if byte == needle[matched] {
                matched + 1
            } else {
                usize::from(byte == needle[0])
            };
            if matched == needle.len() {
                self.skip(bytes - read - 1)?;
                return Ok(true);
            }
        }
        Ok(false)
    }
    pub fn section<T>(
        &mut self,
        bytes: u64,
        work: impl FnOnce(&mut Self) -> Result<T>,
    ) -> Result<T> {
        self.fits(bytes, 1)?;
        let parent = self.end;
        self.end = self.pos + bytes;
        let result = work(self);
        let end = self.end;
        self.end = parent;
        let result = result?;
        self.seek(end)?;
        Ok(result)
    }
}

pub(super) fn add(total: &mut u64, value: u64) -> Result<()> {
    *total = total
        .checked_add(value)
        .context("index preflight decoded count overflows")?;
    Ok(())
}
