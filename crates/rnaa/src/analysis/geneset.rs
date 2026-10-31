//! Gene-set parsing, selection and materialization.
use anyhow::Context;
use anyhow::{Result, bail};
use regex::Regex;
use rnaa_core::model::OutputMode;
use rnaa_core::paths::ProjectPaths;
use std::cmp::Ordering;
use std::collections::{BTreeSet, BinaryHeap};
use std::fs;
use std::path::{Path, PathBuf};

pub(crate) fn parse_output_mode(spec: &str) -> Result<OutputMode> {
    let normalized = spec.trim().to_ascii_lowercase();
    if normalized == "dense" {
        return Ok(OutputMode::Dense);
    }
    if let Some(value) = normalized.strip_prefix("edges:topk=") {
        let k = value
            .parse::<usize>()
            .with_context(|| format!("invalid topk value in output spec '{spec}'"))?;
        return Ok(OutputMode::TopK { k: k.max(1) });
    }
    if let Some(value) = normalized.strip_prefix("edges:threshold=") {
        let min_abs_r = value
            .parse::<f64>()
            .with_context(|| format!("invalid threshold value in output spec '{spec}'"))?;
        return Ok(OutputMode::Threshold { min_abs_r });
    }
    bail!(
        "unsupported output spec '{}'; expected dense | edges:topk=N | edges:threshold=R",
        spec
    )
}

pub(crate) fn materialize_geneset_matrix(
    paths: &ProjectPaths,
    input_matrix: &Path,
    geneset: &str,
    check: &impl Fn() -> Result<()>,
) -> Result<PathBuf> {
    // The adjustment stage reads and validates the input. Avoid reading the
    // complete matrix twice when there is no selection to materialize.
    if geneset.trim().is_empty() || geneset.trim().eq_ignore_ascii_case("all") {
        return Ok(input_matrix.to_path_buf());
    }
    let selection = Selector::parse(geneset)?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input_matrix)?;
    let header = reader.headers()?.clone();
    let samples = header.iter().skip(1).collect::<Vec<_>>();
    if samples.is_empty()
        || samples.iter().any(|id| id.trim().is_empty())
        || samples.iter().collect::<BTreeSet<_>>().len() != samples.len()
    {
        bail!("matrix contains duplicate or empty sample IDs");
    }
    let mut seen = BTreeSet::new();
    let mut selected = Vec::new();
    let mut ranked = BinaryHeap::new();
    let mut record = csv::StringRecord::new();
    let mut rows = 0usize;
    while reader.read_record(&mut record)? {
        if rows.is_multiple_of(256) {
            check()?;
        }
        if record.is_empty() {
            continue;
        }
        let gene = record.get(0).unwrap_or_default();
        if gene.trim().is_empty() || !seen.insert(gene.to_owned()) {
            bail!("matrix contains duplicate or empty gene IDs");
        }
        let values = parse_values(&record, samples.len())?;
        let position = record
            .position()
            .context("matrix record has no position")?
            .clone();
        match &selection {
            Selector::TopVariance(limit) => {
                let mean = values.iter().copied().sum::<f64>() / values.len() as f64;
                let variance = values
                    .iter()
                    .map(|value| {
                        let d = *value - mean;
                        d * d
                    })
                    .sum::<f64>()
                    / values.len() as f64;
                let candidate = RankedRow {
                    variance,
                    index: rows,
                    position,
                };
                if ranked.len() < *limit {
                    ranked.push(candidate);
                } else if ranked.peek().is_some_and(|worst| candidate < *worst) {
                    ranked.pop();
                    ranked.push(candidate);
                }
            }
            Selector::MeanGreater(threshold) => {
                if values.iter().copied().sum::<f64>() / values.len() as f64 > *threshold {
                    selected.push(position);
                }
            }
            Selector::Genes(genes) if genes.contains(gene) => selected.push(position),
            Selector::Pattern(pattern) if pattern.is_match(gene) => selected.push(position),
            _ => {}
        }
        rows += 1;
    }
    if rows == 0 {
        bail!("matrix must have nonempty rows/columns and one value row per gene ID");
    }
    check()?;
    if rows < 2 {
        return Ok(input_matrix.to_path_buf());
    }
    if matches!(selection, Selector::TopVariance(_)) {
        // BinaryHeap orders the least desirable row first. Reverse that order
        // to match the stable descending-variance selection used previously.
        selected = ranked
            .into_sorted_vec()
            .into_iter()
            .map(|row| row.position)
            .collect();
    }
    if selected.len() < 2 {
        bail!(
            "geneset '{}' selected fewer than 2 genes ({})",
            geneset,
            selected.len()
        );
    }
    if selected.len() == rows {
        return Ok(input_matrix.to_path_buf());
    }
    let out_path = paths.corr_dir.join("selected_matrix.tsv");
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&out_path)?;
    let mut output_header = vec!["gene_id".to_string()];
    output_header.extend(samples.iter().map(|id| (*id).to_owned()));
    writer.write_record(&output_header)?;
    for (index, position) in selected.into_iter().enumerate() {
        if index.is_multiple_of(256) {
            check()?;
        }
        if reader.position().byte() != position.byte() {
            reader.seek(position)?;
        }
        if !reader.read_record(&mut record)? {
            bail!("matrix changed while selecting genes");
        }
        let values = parse_values(&record, samples.len())?;
        let mut output = Vec::with_capacity(values.len() + 1);
        output.push(record.get(0).unwrap_or_default().to_owned());
        output.extend(values.iter().map(|value| value.to_string()));
        writer.write_record(&output)?;
    }
    check()?;
    writer.flush()?;
    Ok(out_path)
}

fn parse_values(record: &csv::StringRecord, samples: usize) -> Result<Vec<f64>> {
    let values = record
        .iter()
        .skip(1)
        .map(|field| {
            field
                .parse::<f64>()
                .with_context(|| format!("failed to parse float '{field}'"))
        })
        .collect::<Result<Vec<_>>>()?;
    if values.len() != samples || values.iter().any(|value| !value.is_finite()) {
        bail!("matrix values must be finite and match the sample dimensions");
    }
    Ok(values)
}

enum Selector {
    TopVariance(usize),
    MeanGreater(f64),
    Genes(BTreeSet<String>),
    Pattern(Regex),
}

impl Selector {
    fn parse(geneset: &str) -> Result<Self> {
        let spec = geneset.trim();
        if let Some(value) = spec.strip_prefix("topvar:") {
            let n = value
                .trim()
                .parse::<usize>()
                .with_context(|| format!("invalid topvar geneset: '{geneset}'"))?;
            return Ok(Self::TopVariance(n));
        }
        if let Some(value) = spec.strip_prefix("mean_gt:") {
            let threshold = value
                .trim()
                .parse::<f64>()
                .with_context(|| format!("invalid mean_gt geneset: '{geneset}'"))?;
            return Ok(Self::MeanGreater(threshold));
        }
        if let Some(value) = spec.strip_prefix("genes:") {
            let gene_file = PathBuf::from(value.trim());
            let content = fs::read_to_string(&gene_file)
                .with_context(|| format!("failed to read geneset file {}", gene_file.display()))?;
            let gene_set = content
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .map(str::to_owned)
                .collect();
            return Ok(Self::Genes(gene_set));
        }
        if let Some(value) = spec.strip_prefix("regex:") {
            let pattern = Regex::new(value.trim())
                .with_context(|| format!("invalid regex geneset: '{}'", geneset))?;
            return Ok(Self::Pattern(pattern));
        }
        bail!("unsupported geneset spec '{}'", geneset)
    }
}

#[derive(Debug)]
struct RankedRow {
    variance: f64,
    index: usize,
    position: csv::Position,
}
impl PartialEq for RankedRow {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for RankedRow {}
impl PartialOrd for RankedRow {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for RankedRow {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .variance
            .total_cmp(&self.variance)
            .then_with(|| self.index.cmp(&other.index))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rnaa_formats::matrix::read_matrix_tsv;

    fn select(paths: &ProjectPaths, input: &Path, spec: &str) -> Result<PathBuf> {
        materialize_geneset_matrix(paths, input, spec, &|| Ok(()))
    }

    #[test]
    fn streaming_selection_preserves_order_ties_and_numeric_values() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let paths = ProjectPaths::new(temp.path());
        fs::create_dir_all(&paths.corr_dir)?;
        let input = temp.path().join("matrix.tsv");
        fs::write(
            &input,
            "gene_id\tA\tB\nlow\t1\t1\ntie_first\t0\t4\nhigh\t-3\t3\ntie_second\t2\t6\nmid\t1\t3\n",
        )?;
        let genes = temp.path().join("genes.txt");
        fs::write(&genes, "mid\ntie_first\n")?;
        let cases = [
            (
                "topvar:3".to_owned(),
                vec!["high", "tie_first", "tie_second"],
            ),
            ("topvar:2".to_owned(), vec!["high", "tie_first"]),
            (
                "mean_gt:1".to_owned(),
                vec!["tie_first", "tie_second", "mid"],
            ),
            ("regex:^tie_".to_owned(), vec!["tie_first", "tie_second"]),
            (
                format!("genes:{}", genes.display()),
                vec!["tie_first", "mid"],
            ),
        ];
        let source = read_matrix_tsv(&input)?;
        for (spec, expected) in cases {
            let selected = read_matrix_tsv(&select(&paths, &input, &spec)?)?;
            assert_eq!(selected.row_ids, expected, "{spec}");
            for (gene, values) in selected.row_ids.iter().zip(&selected.values) {
                let index = source.row_ids.iter().position(|item| item == gene).unwrap();
                assert_eq!(*values, source.values[index], "{spec}: {gene}");
            }
        }
        assert_eq!(select(&paths, &input, "topvar:20")?, input);
        Ok(())
    }

    #[test]
    fn selection_validates_unselected_rows_and_rejects_too_few() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let paths = ProjectPaths::new(temp.path());
        fs::create_dir_all(&paths.corr_dir)?;
        let input = temp.path().join("matrix.tsv");
        for invalid in [
            "gene_id\tA\tB\ng1\t1\t2\ng2\t3\t4\nbad\tNaN\t0\n",
            "gene_id\tA\tB\ng1\t1\t2\ng2\t3\t4\ng1\t0\t0\n",
            "gene_id\tA\tA\ng1\t1\t2\ng2\t3\t4\n",
            "gene_id\tA\tB\ng1\t1\t2\ng2\t3\n",
        ] {
            fs::write(&input, invalid)?;
            assert!(select(&paths, &input, "regex:^g").is_err());
        }
        fs::write(&input, "gene_id\tA\tB\ng1\t1\t2\ng2\t3\t4\n")?;
        assert!(select(&paths, &input, "regex:^g1$").is_err());
        Ok(())
    }

    #[test]
    fn cancellation_is_checked_during_scan_and_output() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let paths = ProjectPaths::new(temp.path());
        fs::create_dir_all(&paths.corr_dir)?;
        let input = temp.path().join("matrix.tsv");
        let mut content = String::from("gene_id\tA\tB\n");
        for index in 0..600 {
            content.push_str(&format!("g{index}\t{index}\t{}\n", index + 1));
        }
        fs::write(&input, content)?;
        for fail_at in [2, 5] {
            let calls = std::cell::Cell::new(0);
            let result = materialize_geneset_matrix(&paths, &input, "topvar:300", &|| {
                calls.set(calls.get() + 1);
                if calls.get() == fail_at {
                    bail!("cancelled for test");
                }
                Ok(())
            });
            assert!(
                result
                    .unwrap_err()
                    .to_string()
                    .contains("cancelled for test")
            );
            assert_eq!(calls.get(), fail_at);
        }
        Ok(())
    }

    #[test]
    fn top_variance_matches_stable_full_sort_across_limits() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let paths = ProjectPaths::new(temp.path());
        fs::create_dir_all(&paths.corr_dir)?;
        let input = temp.path().join("matrix.tsv");
        let mut content = String::from("gene_id\tA\tB\tC\n");
        for index in 0..73 {
            let value = index % 7;
            content.push_str(&format!("g{index}\t{value}\t0\t{}\n", -value));
        }
        content.push_str("extreme_a\t1e308\t-1e308\t0\nextreme_b\t-1e308\t1e308\t0\n");
        fs::write(&input, content)?;
        let original = read_matrix_tsv(&input)?;
        let mut expected = original
            .values
            .iter()
            .enumerate()
            .map(|(index, values)| {
                let mean = values.iter().sum::<f64>() / values.len() as f64;
                let variance = values
                    .iter()
                    .map(|value| (value - mean) * (value - mean))
                    .sum::<f64>()
                    / values.len() as f64;
                (index, variance)
            })
            .collect::<Vec<_>>();
        expected.sort_by(|a, b| b.1.total_cmp(&a.1));
        for limit in [2, 3, 7, 8, 17, 50, 74] {
            let output = read_matrix_tsv(&select(&paths, &input, &format!("topvar:{limit}"))?)?;
            let wanted = expected
                .iter()
                .take(limit)
                .map(|(index, _)| original.row_ids[*index].as_str())
                .collect::<Vec<_>>();
            assert_eq!(output.row_ids, wanted, "limit {limit}");
        }
        Ok(())
    }
}
