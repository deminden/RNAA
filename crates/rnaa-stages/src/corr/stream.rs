//! Deterministic row-major output with bounded tiles and row stripes.
use super::prepared::Prepared;
use anyhow::Result;
use rnaa_core::model::OutputMode;
use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    io::{self, Write},
    path::Path,
};

pub(super) struct OutputBudget<W> {
    inner: W,
    remaining: u64,
}

impl<W> OutputBudget<W> {
    pub(super) fn new(inner: W, remaining: u64) -> Self {
        Self { inner, remaining }
    }
}

impl<W: Write> Write for OutputBudget<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if bytes.len() as u64 > self.remaining {
            return Err(io::Error::other(
                "correlation output exceeded its reserved byte limit",
            ));
        }
        let written = self.inner.write(bytes)?;
        self.remaining -= written as u64;
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

#[derive(Copy, Clone)]
pub(super) struct Edge {
    pub(super) target: usize,
    pub(super) value: f64,
}

pub(super) fn offer(heap: &mut BinaryHeap<Edge>, edge: Edge, k: usize) {
    if k == 0 {
        return;
    }
    if heap.len() < k {
        heap.push(edge);
    } else if edge < *heap.peek().expect("nonempty edge heap") {
        *heap.peek_mut().expect("nonempty edge heap") = edge;
    }
}
impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Edge {}
impl PartialOrd for Edge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Edge {
    fn cmp(&self, other: &Self) -> Ordering {
        // Heap root is the worst retained edge: smallest magnitude, last target.
        other
            .value
            .abs()
            .total_cmp(&self.value.abs())
            .then(self.target.cmp(&other.target))
    }
}

#[derive(Default, serde::Serialize)]
pub(super) struct Statistics {
    pub evaluated_entries: u64,
    pub undefined_entries: u64,
    pub peak_block_entries: usize,
    pub peak_retained_edges: usize,
    pub peak_stripe_entries: usize,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn write(
    path: &Path,
    genes: &[String],
    prepared: &Prepared,
    mode: &OutputMode,
    block_rows: usize,
    output_budget: u64,
    pool: &rayon::ThreadPool,
    check: &impl Fn() -> Result<()>,
) -> Result<Statistics> {
    check()?;
    let output = OutputBudget {
        inner: std::fs::File::create(path)?,
        remaining: output_budget,
    };
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(output);
    let mut statistics = Statistics::default();
    let mut tiles = 0usize;
    if matches!(mode, OutputMode::Dense) {
        writer.write_record(std::iter::once("gene_id").chain(genes.iter().map(String::as_str)))?;
    } else {
        writer.write_record(["gene_a", "gene_b", "r"])?;
    }
    let k = match mode {
        OutputMode::TopK { k } => (*k).min(genes.len() - 1),
        _ => 0,
    };
    for start in (0..genes.len()).step_by(block_rows) {
        check()?;
        let end = (start + block_rows).min(genes.len());
        let mut heaps = (start..end)
            .map(|_| BinaryHeap::<Edge>::new())
            .collect::<Vec<_>>();
        let mut stripes = vec![Vec::<Edge>::new(); end - start];
        let first_column = if matches!(mode, OutputMode::Threshold { .. }) {
            start
        } else {
            0
        };
        for j in (first_column..genes.len()).step_by(block_rows) {
            if tiles.is_multiple_of(16) {
                check()?;
            }
            let column_end = (j + block_rows).min(genes.len());
            let block = pool.install(|| prepared.block(start..end, j..column_end));
            tiles += 1;
            statistics.peak_block_entries = statistics.peak_block_entries.max(block.len());
            for ((row, column), &value) in block.indexed_iter() {
                let (gene, target) = (start + row, j + column);
                statistics.evaluated_entries += 1;
                if !value.is_finite() {
                    statistics.undefined_entries += 1;
                }
                let edge = Edge { target, value };
                match mode {
                    OutputMode::Dense => stripes[row].push(edge),
                    OutputMode::Threshold { min_abs_r }
                        if target > gene && value.is_finite() && value.abs() >= *min_abs_r =>
                    {
                        stripes[row].push(edge)
                    }
                    OutputMode::TopK { .. } if target != gene && value.is_finite() && k > 0 => {
                        let heap = &mut heaps[row];
                        offer(heap, edge, k);
                    }
                    _ => {}
                }
            }
        }
        statistics.peak_retained_edges = statistics
            .peak_retained_edges
            .max(heaps.iter().map(BinaryHeap::len).sum());
        statistics.peak_stripe_entries = statistics
            .peak_stripe_entries
            .max(stripes.iter().map(Vec::len).sum());
        for (row, (heap, mut stripe)) in heaps.into_iter().zip(stripes).enumerate() {
            let gene = &genes[start + row];
            if matches!(mode, OutputMode::Dense) {
                writer.write_field(gene)?;
                for edge in stripe {
                    writer.write_field(edge.value.to_string())?;
                }
                writer.write_record(std::iter::empty::<&str>())?;
                continue;
            }
            if matches!(mode, OutputMode::TopK { .. }) {
                stripe = heap.into_vec();
                stripe.sort();
            }
            for edge in stripe {
                writer.write_record([
                    gene.as_str(),
                    genes[edge.target].as_str(),
                    &edge.value.to_string(),
                ])?;
            }
        }
    }
    check()?;
    writer.flush()?;
    Ok(statistics)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rnaa_core::model::CorrelationMethod;

    #[test]
    fn primary_output_never_exceeds_the_reserved_bytes() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let path = temp.path().join("network.tsv");
        let genes = ["a\tquoted", "b\"quoted", "c"].map(str::to_owned);
        let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build()?;
        let prepared = Prepared::new(
            vec![
                vec![1.0, 2.0, 3.0],
                vec![3.0, 2.0, 1.0],
                vec![1.0, 3.0, 2.0],
            ],
            CorrelationMethod::Pearson,
        );
        for mode in [
            OutputMode::Dense,
            OutputMode::TopK { k: 1 },
            OutputMode::Threshold { min_abs_r: 0.0 },
        ] {
            write(&path, &genes, &prepared, &mode, 2, u64::MAX, &pool, &|| {
                Ok(())
            })?;
            let expected = std::fs::read(&path)?;
            let bytes = expected.len() as u64;
            write(&path, &genes, &prepared, &mode, 2, bytes, &pool, &|| Ok(()))?;
            assert_eq!(std::fs::read(&path)?, expected);
            let error = write(
                &path,
                &genes,
                &prepared,
                &mode,
                2,
                bytes - 1,
                &pool,
                &|| Ok(()),
            )
            .err()
            .expect("output exceeded limit");
            assert!(error.to_string().contains("reserved byte limit"));
            assert!(std::fs::metadata(&path)?.len() < bytes);
        }
        Ok(())
    }

    #[test]
    fn preflight_budget_covers_each_output_mode_with_quoted_identifiers() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let matrix_path = temp.path().join("matrix.tsv");
        let path = temp.path().join("network.tsv");
        let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build()?;
        let matrix = rnaa_formats::matrix::MatrixData {
            row_ids: vec!["a\tquoted".into(), "b\"quoted".into(), "c\nquoted".into()],
            col_ids: vec!["s1".into(), "s2".into(), "s3".into()],
            values: vec![
                vec![1.0, 2.0, 3.0],
                vec![3.0, 2.0, 1.0],
                vec![1.0, 3.0, 2.0],
            ],
        };
        rnaa_formats::matrix::write_matrix_tsv(&matrix_path, &matrix)?;
        let prepared = Prepared::new(matrix.values, CorrelationMethod::Pearson);
        for mode in [
            OutputMode::Dense,
            OutputMode::TopK { k: 2 },
            OutputMode::Threshold { min_abs_r: 0.0 },
        ] {
            let estimate = super::super::limits::CorrelationLimits::default().preflight(
                &matrix_path,
                CorrelationMethod::Pearson,
                &mode,
                1,
            )?;
            write(
                &path,
                &matrix.row_ids,
                &prepared,
                &mode,
                2,
                u64::try_from(estimate.output_bytes)?,
                &pool,
                &|| Ok(()),
            )?;
            assert!(u128::from(std::fs::metadata(&path)?.len()) <= estimate.output_bytes);
        }
        Ok(())
    }

    #[test]
    fn cancellation_is_checked_between_tiles() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let path = temp.path().join("network.tsv");
        let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build()?;
        let genes = (0..8).map(|i| format!("g{i}")).collect::<Vec<_>>();
        let prepared = Prepared::new(
            (0..8)
                .map(|i| vec![i as f64, (i + 1) as f64, (i + 2) as f64])
                .collect(),
            CorrelationMethod::Pearson,
        );
        let calls = std::cell::Cell::new(0);
        let result = write(
            &path,
            &genes,
            &prepared,
            &OutputMode::Dense,
            2,
            1024 * 1024,
            &pool,
            &|| {
                calls.set(calls.get() + 1);
                if calls.get() == 5 {
                    anyhow::bail!("cancelled between tiles");
                }
                Ok(())
            },
        );
        assert!(
            result
                .err()
                .expect("must cancel")
                .to_string()
                .contains("cancelled between tiles")
        );
        assert_eq!(calls.get(), 5);
        Ok(())
    }
}
