//! Exact clean-kmer cardinality using bounded sorting and external merge runs.
//! Only the cardinality escapes; neither a graph nor an all-kmer set is retained.
use anyhow::{Context, Result};
use rnaa_core::resources::{
    DiskRequest, ExecutionResources, ResourceAdmissionError, ResourceRequest, ResourceReservation,
};
use std::{
    cmp::Reverse,
    collections::BinaryHeap,
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Write},
    path::PathBuf,
};

const WORDS: usize = 1 << 19;
const FAN_IN: usize = 32;
const WINDOW: u64 = 16 * 1024 * 1024;
const METADATA: u64 = 4 * 1024 * 1024;
type Heads = BinaryHeap<Reverse<(u64, usize)>>;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy)]
struct Run {
    id: u64,
    words: u64,
}
struct Cursor {
    reader: BufReader<File>,
    remaining: u64,
}
impl Cursor {
    fn next(&mut self) -> Result<Option<u64>> {
        if self.remaining == 0 {
            return Ok(None);
        }
        let mut bytes = [0; 8];
        self.reader.read_exact(&mut bytes)?;
        self.remaining -= 1;
        Ok(Some(u64::from_le_bytes(bytes)))
    }
}
pub(super) struct UniqueWords<'a> {
    directory: tempfile::TempDir,
    words: Vec<u64>,
    levels: Vec<Vec<Run>>,
    sequence: u64,
    writes: u64,
    promised_until: u64,
    check: &'a dyn Fn() -> Result<()>,
    // Last field: scratch files and all buffers end before their owner releases.
    guard: ResourceReservation,
}
impl<'a> UniqueWords<'a> {
    pub fn new(
        context: &ExecutionResources,
        memory: u64,
        workspace: &std::path::Path,
        check: &'a dyn Fn() -> Result<()>,
    ) -> Result<Self> {
        let guard = context.reserve(
            "fasta_plan",
            ResourceRequest {
                cpu_threads: 1,
                memory_bytes: memory,
                disks: vec![DiskRequest {
                    path: workspace.into(),
                    bytes: WINDOW + METADATA,
                }],
            },
            &|| check(),
        )?;
        Ok(Self {
            directory: tempfile::Builder::new()
                .prefix("rnaa-reference-plan-")
                .tempdir_in(workspace)?,
            words: Vec::with_capacity(WORDS),
            levels: Vec::new(),
            sequence: 0,
            writes: 0,
            promised_until: WINDOW,
            check,
            guard,
        })
    }
    fn path(&self, run: Run) -> PathBuf {
        self.directory.path().join(format!("{}.keys", run.id))
    }
    fn new_run(&mut self) -> Result<Run> {
        let id = self.sequence;
        self.sequence = self
            .sequence
            .checked_add(1)
            .context("reference planning run counter overflows")?;
        Ok(Run { id, words: 0 })
    }
    fn write(&mut self, output: &mut BufWriter<File>, bytes: &[u8]) -> Result<()> {
        let next = self
            .writes
            .checked_add(bytes.len() as u64)
            .context("reference planning writes overflow")?;
        if next > self.promised_until {
            (self.check)()?;
            output.flush()?;
            let claim = self.guard.snapshot().requested.clone();
            if self
                .guard
                .try_resize(claim)
                .map_err(ResourceAdmissionError::new)?
                .is_some()
            {
                return Err(ResourceAdmissionError::new(anyhow::anyhow!(
                    "reference planner scratch cannot grow within shared disk capacity"
                ))
                .into());
            }
            self.promised_until = self
                .writes
                .checked_add(WINDOW)
                .context("reference planning window overflows")?;
        }
        output.write_all(bytes)?;
        self.writes = next;
        Ok(())
    }
    pub fn push(&mut self, word: u64) -> Result<()> {
        self.words.push(word);
        if self.words.len() == WORDS {
            self.flush()?;
        }
        Ok(())
    }
    fn flush(&mut self) -> Result<()> {
        if self.words.is_empty() {
            return Ok(());
        }
        (self.check)()?;
        let mut words = std::mem::take(&mut self.words);
        words.sort_unstable();
        words.dedup();
        let mut run = self.new_run()?;
        let mut output = BufWriter::with_capacity(65536, File::create(self.path(run))?);
        let mut bytes = [0_u8; 65536];
        for chunk in words.chunks(bytes.len() / 8) {
            for (word, slot) in chunk.iter().zip(bytes.as_chunks_mut::<8>().0) {
                slot.copy_from_slice(&word.to_le_bytes());
            }
            self.write(&mut output, &bytes[..chunk.len() * 8])?;
        }
        output.flush()?;
        run.words = words.len() as u64;
        words.clear();
        self.words = words;
        self.add(0, run)
    }
    fn add(&mut self, level: usize, run: Run) -> Result<()> {
        if level >= self.levels.len() {
            self.levels.push(Vec::new());
        }
        self.levels[level].push(run);
        if self.levels[level].len() == FAN_IN {
            let inputs = std::mem::take(&mut self.levels[level]);
            let merged = self.merge(&inputs)?;
            self.add(level + 1, merged)?;
        }
        Ok(())
    }
    fn cursors(&self, runs: &[Run]) -> Result<(Vec<Cursor>, Heads)> {
        let mut inputs = Vec::with_capacity(runs.len());
        let mut heap = BinaryHeap::new();
        for (i, &run) in runs.iter().enumerate() {
            let mut input = Cursor {
                reader: BufReader::new(File::open(self.path(run))?),
                remaining: run.words,
            };
            if let Some(word) = input.next()? {
                heap.push(Reverse((word, i)));
            }
            inputs.push(input);
        }
        Ok((inputs, heap))
    }
    fn merge(&mut self, runs: &[Run]) -> Result<Run> {
        (self.check)()?;
        let (mut inputs, mut heap) = self.cursors(runs)?;
        let mut output_run = self.new_run()?;
        let mut output = BufWriter::with_capacity(65536, File::create(self.path(output_run))?);
        let mut last = None;
        let mut observed = 0_u64;
        while let Some(Reverse((word, i))) = heap.pop() {
            if observed.is_multiple_of(1 << 20) {
                (self.check)()?;
            }
            observed += 1;
            if last != Some(word) {
                self.write(&mut output, &word.to_le_bytes())?;
                output_run.words += 1;
                last = Some(word);
            }
            if let Some(word) = inputs[i].next()? {
                heap.push(Reverse((word, i)));
            }
        }
        output.flush()?;
        drop(output);
        drop(inputs);
        for &run in runs {
            fs::remove_file(self.path(run))?;
        }
        Ok(output_run)
    }
    pub fn finish(mut self) -> Result<u64> {
        self.flush()?;
        let mut runs = self.levels.iter().flatten().copied().collect::<Vec<_>>();
        while runs.len() > FAN_IN {
            let mut next = Vec::new();
            for inputs in runs.chunks(FAN_IN) {
                next.push(if inputs.len() == 1 {
                    inputs[0]
                } else {
                    self.merge(inputs)?
                });
            }
            runs = next;
        }
        let (mut inputs, mut heap) = self.cursors(&runs)?;
        let (mut last, mut count, mut observed) = (None, 0_u64, 0_u64);
        while let Some(Reverse((word, i))) = heap.pop() {
            if observed.is_multiple_of(1 << 20) {
                (self.check)()?;
            }
            observed += 1;
            if last != Some(word) {
                count += 1;
                last = Some(word);
            }
            if let Some(word) = inputs[i].next()? {
                heap.push(Reverse((word, i)));
            }
        }
        Ok(count)
    }
}
