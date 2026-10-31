//! Pool fractional transcript counts before gene import, scaling or rounding.
use super::{ImportedCounts, import_kallisto_abundances};
use anyhow::{Context, Result, bail};
use serde::Deserialize;
use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};

#[derive(Debug, Clone)]
pub struct RunAbundance {
    pub run_id: String,
    pub path: PathBuf,
}
#[derive(Debug, Clone)]
pub struct QuantifiedSample {
    pub sample_id: String,
    pub runs: Vec<RunAbundance>,
}

#[derive(Deserialize)]
struct Row {
    target_id: String,
    length: f64,
    eff_length: f64,
    est_counts: f64,
    tpm: f64,
}
#[derive(Default)]
struct Pooled {
    length: f64,
    counts: f64,
    rate: f64,
    lengths: f64,
}

/// RNAA pooled-estimates policy v1:
/// C_t = sum_r C_tr; M_t = sum_r C_tr / L_tr; L_t = C_t / M_t;
/// TPM_t = 1e6 M_t / sum_t M_t. All-zero targets use mean effective length.
/// This pools estimated transcript mass, not raw reads or bootstrap replicates.
/// Single-run samples retain their original quantifier TPM/length estimates.
pub fn import_grouped_kallisto_abundances(
    samples: &[QuantifiedSample],
    mapping: &Path,
    mode: &str,
) -> Result<ImportedCounts> {
    let mut seen_samples = BTreeSet::new();
    let mut seen_runs = BTreeSet::new();
    for sample in samples {
        if sample.sample_id.trim().is_empty()
            || !seen_samples.insert(&sample.sample_id)
            || sample.runs.is_empty()
        {
            bail!("biological sample IDs must be unique/nonempty and have contributing runs");
        }
        for run in &sample.runs {
            if run.run_id.trim().is_empty() || !seen_runs.insert(&run.run_id) {
                bail!("each sequencing run must contribute to exactly one biological sample");
            }
        }
    }
    let scratch = tempfile::tempdir()?;
    let mut ordered = samples.iter().collect::<Vec<_>>();
    ordered.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));
    let mut paths = Vec::new();
    for (index, sample) in ordered.into_iter().enumerate() {
        let path = if sample.runs.len() == 1 {
            sample.runs[0].path.clone()
        } else {
            let output = scratch.path().join(format!("pooled_{index}.tsv"));
            pool(sample, &output)?;
            output
        };
        paths.push((sample.sample_id.clone(), path));
    }
    import_kallisto_abundances(&paths, mapping, mode)
}

fn pool(sample: &QuantifiedSample, output: &Path) -> Result<()> {
    let mut targets = BTreeMap::<String, Pooled>::new();
    let mut ordered = sample.runs.iter().collect::<Vec<_>>();
    ordered.sort_by(|a, b| a.run_id.cmp(&b.run_id));
    for (index, run) in ordered.into_iter().enumerate() {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&run.path)
            .with_context(|| format!("failed to read abundance for run {}", run.run_id))?;
        let header = reader.headers()?;
        let required = ["target_id", "length", "eff_length", "est_counts", "tpm"]
            .into_iter()
            .collect::<BTreeSet<_>>();
        if header.len() != 5 || header.iter().collect::<BTreeSet<_>>() != required {
            bail!("invalid or duplicate abundance columns for {}", run.run_id);
        }
        let mut seen = BTreeSet::new();
        let mut normalized = BTreeSet::new();
        for row in reader.deserialize::<Row>() {
            let row = row?;
            let normalized_id = rnaa_core::util::normalize_feature_id(&row.target_id);
            if normalized_id.is_empty()
                || !normalized.insert(normalized_id)
                || !seen.insert(row.target_id.clone())
            {
                bail!(
                    "duplicate/empty normalized transcript in run {}",
                    run.run_id
                );
            }
            if !row.length.is_finite()
                || row.length <= 0.0
                || !row.eff_length.is_finite()
                || row.eff_length < 0.0
                || (row.eff_length == 0.0 && (row.est_counts > 0.0 || row.tpm > 0.0))
                || !row.est_counts.is_finite()
                || row.est_counts < 0.0
                || !row.tpm.is_finite()
                || row.tpm < 0.0
            {
                bail!(
                    "invalid abundance for {} in run {}",
                    row.target_id,
                    run.run_id
                );
            }
            if index > 0 && !targets.contains_key(&row.target_id) {
                bail!(
                    "transcript set differs among runs of biological sample {}",
                    sample.sample_id
                );
            }
            let value = targets.entry(row.target_id).or_default();
            if index > 0 && value.length != row.length {
                bail!(
                    "transcript lengths differ among runs of biological sample {}",
                    sample.sample_id
                );
            }
            value.length = row.length;
            value.counts += row.est_counts;
            if row.est_counts > 0.0 {
                value.rate += row.est_counts / row.eff_length;
            }
            value.lengths += row.eff_length;
        }
        if seen.len() != targets.len() || seen.is_empty() {
            bail!("incomplete transcript set in run {}", run.run_id);
        }
    }
    let rate_sum = targets.values().map(|value| value.rate).sum::<f64>();
    if !rate_sum.is_finite() {
        bail!("pooled transcript mass overflowed");
    }
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output)?;
    writer.write_record(["target_id", "length", "eff_length", "est_counts", "tpm"])?;
    for (id, value) in targets {
        let effective_length = if value.rate > 0.0 {
            value.counts / value.rate
        } else if value.counts == 0.0 {
            value.lengths / sample.runs.len() as f64
        } else {
            bail!("pooled transcript mass underflowed for {id}");
        };
        let tpm = if rate_sum > 0.0 {
            value.rate / rate_sum * 1e6
        } else {
            0.0
        };
        if !value.counts.is_finite()
            || !effective_length.is_finite()
            || effective_length < 0.0
            || (effective_length == 0.0 && value.counts > 0.0)
            || !tpm.is_finite()
        {
            bail!("invalid pooled abundance for {id}");
        }
        writer.write_record([
            id,
            value.length.to_string(),
            effective_length.to_string(),
            value.counts.to_string(),
            tpm.to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests;
