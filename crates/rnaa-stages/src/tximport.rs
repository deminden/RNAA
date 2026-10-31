//! Native Rust equivalent of the kallisto-to-gene import performed by tximport.
//!
//! Counts are assembled in deterministic transcript and sample order. The
//! abundance table emitted by kallistors is retained as the exchange contract;
//! HDF5 remains the canonical per-run quantification artifact.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use csv::ReaderBuilder;
use rnaa_core::util::normalize_feature_id;
use serde::{Deserialize, Serialize};
mod grouped;
pub use grouped::{QuantifiedSample, RunAbundance, import_grouped_kallisto_abundances};

/// Bump whenever import or offset semantics change.
pub const IMPORT_ALGORITHM_VERSION: u32 = 4;

#[derive(Debug, Clone)]
pub struct ImportedCounts {
    pub gene_ids: Vec<String>,
    pub sample_ids: Vec<String>,
    /// Gene-major, sample-minor integer matrix matching rsdeseq2 CountMatrix.
    pub counts: Vec<u32>,
    /// Optional tximport-style gene/sample length normalization factors.
    pub normalization_factors: Option<Vec<f64>>,
    pub diagnostics: Vec<ImportDiagnostics>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ImportDiagnostics {
    pub sample_id: String,
    pub mapped_transcripts: usize,
    pub unmapped_transcripts: usize,
    pub mapped_count_fraction: f64,
}

#[derive(Debug, Deserialize)]
struct AbundanceRow {
    target_id: String,
    eff_length: f64,
    est_counts: f64,
    tpm: f64,
}

/// Aggregate kallisto-compatible abundance tables to genes.
///
/// Supports gene-level `no`, `scaledTPM`, and `lengthScaledTPM`.
/// `dtuScaledTPM` requires transcript-level output and is rejected here.
pub fn import_kallisto_abundances(
    samples: &[(String, PathBuf)],
    tx2gene_path: &Path,
    counts_from_abundance: &str,
) -> Result<ImportedCounts> {
    if samples.len() < 2 {
        bail!("need at least two quantified samples for differential expression");
    }
    let sample_names = samples.iter().map(|(id, _)| id).collect::<BTreeSet<_>>();
    if sample_names.len() != samples.len() || sample_names.iter().any(|id| id.trim().is_empty()) {
        bail!("sample identifiers must be unique and nonempty");
    }
    let mode = counts_from_abundance.to_ascii_lowercase();
    if mode == "dtuscaledtpm" {
        bail!(
            "dtuScaledTPM requires transcript-level output; use no, scaledTPM, or lengthScaledTPM for gene-level DE"
        );
    }
    if !matches!(mode.as_str(), "no" | "scaledtpm" | "lengthscaledtpm") {
        bail!("unsupported counts_from_abundance mode '{counts_from_abundance}'");
    }

    let mut tx2gene_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(tx2gene_path)
        .with_context(|| {
            format!(
                "failed to open transcript-to-gene table {}",
                tx2gene_path.display()
            )
        })?;
    let headers = tx2gene_reader.headers()?.clone();
    if headers.len() < 2 {
        bail!("transcript-to-gene table must have at least two columns");
    }
    let mut tx2gene = HashMap::new();
    for row in tx2gene_reader.records() {
        let row = row?;
        let tx = normalize_feature_id(row.get(0).unwrap_or_default());
        let gene = normalize_feature_id(row.get(1).unwrap_or_default());
        if tx.is_empty() || gene.is_empty() {
            bail!("transcript-to-gene table contains an empty identifier");
        }
        if let Some(previous) = tx2gene.insert(tx.clone(), gene.clone())
            && previous != gene
        {
            bail!("conflicting transcript-to-gene mapping for '{tx}': '{previous}' and '{gene}'");
        }
    }
    if tx2gene.is_empty() {
        bail!("transcript-to-gene table contains no usable mappings");
    }

    let mut genes = BTreeSet::new();
    let mut sample_values = Vec::with_capacity(samples.len());
    let mut transcript_set = None;
    let mut diagnostics = Vec::new();
    let mut fallback_lengths: BTreeMap<String, (f64, usize)> = BTreeMap::new();
    for (sample_id, abundance_path) in samples {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(abundance_path)
            .with_context(|| {
                format!(
                    "failed to open abundance table {}",
                    abundance_path.display()
                )
            })?;
        let mut gene_counts: BTreeMap<String, f64> = BTreeMap::new();
        let mut gene_tpm: BTreeMap<String, f64> = BTreeMap::new();
        let mut gene_weighted_length: BTreeMap<String, f64> = BTreeMap::new();
        let mut library_size = 0.0;
        let mut all_counts = 0.0;
        let mut seen = BTreeSet::new();
        let mut mapped = 0;
        let mut unmapped = 0;
        for row in reader.deserialize::<AbundanceRow>() {
            let row = row.with_context(|| {
                format!(
                    "invalid kallisto abundance row in {}",
                    abundance_path.display()
                )
            })?;
            let tx = normalize_feature_id(&row.target_id);
            if tx.is_empty() || !seen.insert(tx.clone()) {
                bail!("empty or duplicate normalized transcript '{tx}' in sample '{sample_id}'");
            }
            if !row.est_counts.is_finite()
                || row.est_counts < 0.0
                || !row.tpm.is_finite()
                || row.tpm < 0.0
                || !row.eff_length.is_finite()
                || row.eff_length < 0.0
                || (row.eff_length == 0.0 && (row.est_counts > 0.0 || row.tpm > 0.0))
            {
                bail!(
                    "invalid abundance values for transcript '{}' in sample '{sample_id}'",
                    row.target_id
                );
            }
            all_counts += row.est_counts;
            let Some(gene) = tx2gene.get(&tx) else {
                unmapped += 1;
                continue;
            };
            mapped += 1;
            let fallback = fallback_lengths.entry(gene.clone()).or_default();
            fallback.0 += row.eff_length;
            fallback.1 += 1;
            *gene_counts.entry(gene.clone()).or_default() += row.est_counts;
            *gene_tpm.entry(gene.clone()).or_default() += row.tpm;
            *gene_weighted_length.entry(gene.clone()).or_default() += row.tpm * row.eff_length;
            library_size += row.est_counts;
            genes.insert(gene.clone());
        }
        if let Some(expected) = &transcript_set {
            if expected != &seen {
                bail!(
                    "transcript set differs between abundance tables for sample '{sample_id}'; check reference identity"
                );
            }
        } else {
            transcript_set = Some(seen);
        }
        if gene_counts.is_empty() {
            bail!("no abundance transcripts mapped to genes for sample '{sample_id}'");
        }
        if !all_counts.is_finite() || !library_size.is_finite() {
            bail!("abundance totals overflowed for sample '{sample_id}'");
        }
        diagnostics.push(ImportDiagnostics {
            sample_id: sample_id.clone(),
            mapped_transcripts: mapped,
            unmapped_transcripts: unmapped,
            mapped_count_fraction: if all_counts > 0.0 {
                library_size / all_counts
            } else {
                1.0
            },
        });
        let mut gene_lengths = BTreeMap::new();
        for gene in &genes {
            let tpm = gene_tpm.get(gene).copied().unwrap_or(0.0);
            let length = if tpm > 0.0 {
                gene_weighted_length.get(gene).copied().unwrap_or_default() / tpm
            } else {
                f64::NAN
            };
            gene_lengths.insert(gene.clone(), length);
        }
        sample_values.push((gene_counts, gene_tpm, gene_lengths, library_size));
    }

    let gene_ids = genes.into_iter().collect::<Vec<_>>();
    let sample_ids = samples
        .iter()
        .map(|(sample, _)| sample.clone())
        .collect::<Vec<_>>();
    let n_genes = gene_ids.len();
    let n_samples = sample_ids.len();
    // Impute absent weighted lengths as tximport does: geometric mean of
    // expressed samples, or the unweighted transcript/sample mean if all zero.
    let mut raw_lengths = vec![0.0; n_genes * n_samples];
    let mut average_lengths = vec![0.0; n_genes];
    for (gene_idx, gene) in gene_ids.iter().enumerate() {
        let observed = sample_values
            .iter()
            .map(|(_, _, lengths, _)| lengths[gene])
            .filter(|v| v.is_finite() && *v > 0.0)
            .collect::<Vec<_>>();
        let replacement = if observed.is_empty() {
            let (sum, count) = fallback_lengths[gene];
            sum / count as f64
        } else {
            (observed.iter().map(|v| v.ln()).sum::<f64>() / observed.len() as f64).exp()
        };
        for (sample_idx, (_, _, lengths, _)) in sample_values.iter().enumerate() {
            let length = lengths[gene];
            let length = if length.is_finite() {
                length
            } else {
                replacement
            };
            if !length.is_finite() || length <= 0.0 {
                bail!("invalid aggregated effective length for '{gene}'");
            }
            raw_lengths[gene_idx * n_samples + sample_idx] = length;
            average_lengths[gene_idx] += length / n_samples as f64;
        }
    }
    let mut counts = vec![0_u32; n_genes * n_samples];
    for (sample_idx, (raw, tpm, _, library_size)) in sample_values.iter().enumerate() {
        let scaled_sum = gene_ids
            .iter()
            .enumerate()
            .map(|(gene_idx, gene)| {
                tpm[gene]
                    * if mode == "lengthscaledtpm" {
                        average_lengths[gene_idx]
                    } else {
                        1.0
                    }
            })
            .sum::<f64>();
        if !scaled_sum.is_finite() {
            bail!("scaled abundance totals overflowed");
        }
        for (gene_idx, gene) in gene_ids.iter().enumerate() {
            let value = if mode == "no" {
                raw[gene]
            } else if scaled_sum > 0.0 {
                tpm[gene]
                    * if mode == "lengthscaledtpm" {
                        average_lengths[gene_idx]
                    } else {
                        1.0
                    }
                    * library_size
                    / scaled_sum
            } else {
                0.0
            };
            if !value.is_finite() || value < 0.0 || value > u32::MAX as f64 {
                bail!("gene count is outside rsdeseq2's supported u32 range for '{gene}'");
            }
            // R's round() uses ties-to-even, including after transcript aggregation.
            counts[gene_idx * n_samples + sample_idx] = value.round_ties_even() as u32;
        }
    }
    let normalization_factors = if mode == "no" {
        Some(length_normalization_factors(
            &raw_lengths,
            n_genes,
            n_samples,
        ))
    } else {
        None
    };
    Ok(ImportedCounts {
        gene_ids,
        sample_ids,
        counts,
        normalization_factors,
        diagnostics,
    })
}

fn length_normalization_factors(lengths: &[f64], n_genes: usize, n_samples: usize) -> Vec<f64> {
    let mut factors = vec![1.0; lengths.len()];
    for gene in 0..n_genes {
        let row = &lengths[gene * n_samples..(gene + 1) * n_samples];
        let geometric_mean = (row
            .iter()
            .filter(|v| **v > 0.0)
            .map(|v| v.ln())
            .sum::<f64>()
            / row.iter().filter(|v| **v > 0.0).count().max(1) as f64)
            .exp();
        for sample in 0..n_samples {
            factors[gene * n_samples + sample] = row[sample] / geometric_mean;
        }
    }
    factors
}

/// Complete DESeq2-style offsets: median ratios on length-corrected counts,
/// followed by row geometric centering of the product. Never round corrected counts.
pub fn complete_length_offsets(
    counts: &[u32],
    factors: &mut [f64],
    n_samples: usize,
) -> Result<Vec<f64>> {
    if n_samples == 0
        || counts.len() != factors.len()
        || !counts.len().is_multiple_of(n_samples)
        || factors.iter().any(|v| !v.is_finite() || *v <= 0.0)
    {
        bail!("invalid count/length normalization matrix");
    }
    let mut ratios = vec![Vec::new(); n_samples];
    for (row, lengths) in counts.chunks(n_samples).zip(factors.chunks(n_samples)) {
        if row.contains(&0) {
            continue;
        }
        let logs = row
            .iter()
            .zip(lengths)
            .map(|(count, length)| (*count as f64 / length).ln())
            .collect::<Vec<_>>();
        let mean = logs.iter().sum::<f64>() / n_samples as f64;
        for (sample, log) in logs.iter().enumerate() {
            ratios[sample].push(log - mean);
        }
    }
    let mut size_factors = Vec::new();
    for mut values in ratios {
        if values.is_empty() {
            bail!(
                "cannot estimate ratio size factors: every gene contains a zero; supply a less sparse cohort"
            );
        }
        values.sort_by(f64::total_cmp);
        let middle = values.len() / 2;
        let median = if values.len().is_multiple_of(2) {
            (values[middle - 1] + values[middle]) / 2.0
        } else {
            values[middle]
        };
        size_factors.push(median.exp());
    }
    let center = (size_factors.iter().map(|v| v.ln()).sum::<f64>() / n_samples as f64).exp();
    for value in &mut size_factors {
        *value /= center;
    }
    for row in factors.chunks_mut(n_samples) {
        for (value, sf) in row.iter_mut().zip(&size_factors) {
            *value *= sf;
        }
        let center = (row.iter().map(|v| v.ln()).sum::<f64>() / n_samples as f64).exp();
        for value in row {
            *value /= center;
        }
    }
    Ok(size_factors)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn import(rows: &[&str], mapping: &str, mode: &str) -> Result<ImportedCounts> {
        let temp = TempDir::new()?;
        let map = temp.path().join("tx2gene.tsv");
        std::fs::write(&map, format!("transcript_id\tgene_id\n{mapping}"))?;
        let samples = rows
            .iter()
            .enumerate()
            .map(|(idx, rows)| {
                let path = temp.path().join(format!("{idx}.tsv"));
                std::fs::write(
                    &path,
                    format!("target_id\teff_length\test_counts\ttpm\n{rows}"),
                )?;
                Ok((format!("S{idx}"), path))
            })
            .collect::<Result<Vec<_>>>()?;
        import_kallisto_abundances(&samples, &map, mode)
    }

    #[test]
    fn length_scaled_tpm_uses_across_sample_average_lengths() {
        let a = "t1\t100\t100\t500000\nt2\t100\t100\t500000\n";
        let b = "t2\t100\t50\t500000\nt1\t300\t150\t500000\n";
        let data = import(&[a, b], "t1\tG1\nt2\tG2\n", "lengthScaledTPM").unwrap();
        // Arithmetic oracle: mean lengths [200,100], TPM weights 2:1,
        // each library has 200 counts. R round(c(400/3,200/3)) = [133,67].
        assert_eq!(data.counts, [133, 133, 67, 67]);
        assert!(data.normalization_factors.is_none());
        let reverse = import(&[b, a], "t1\tG1\nt2\tG2\n", "lengthScaledTPM").unwrap();
        assert_eq!(data.counts, reverse.counts);
    }

    #[test]
    fn zero_abundance_lengths_are_imputed_from_expressed_samples() {
        let a = "t1\t100\t0\t0\nt2\t500\t0\t0\nt3\t100\t100\t1000000\n";
        let b = "t1\t100\t0\t0\nt2\t500\t100\t500000\nt3\t100\t100\t500000\n";
        let data = import(&[a, b], "t1\tG1\nt2\tG1\nt3\tG2\n", "no").unwrap();
        for value in data.normalization_factors.unwrap() {
            assert!((value - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn all_zero_gene_uses_unweighted_transcript_sample_average() {
        let a = "t1\t100\t0\t0\nt2\t500\t0\t0\nt3\t100\t100\t1000000\n";
        let b = "t1\t300\t0\t0\nt2\t700\t0\t0\nt3\t100\t100\t1000000\n";
        let data = import(&[a, b], "t1\tG1\nt2\tG1\nt3\tG2\n", "no").unwrap();
        assert_eq!(data.counts, [0, 0, 100, 100]);
        for value in data.normalization_factors.unwrap() {
            assert!((value - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn gencode_headers_and_par_y_identifiers_map_consistently() {
        let a = "ENST1.2|ENSG1.3|label|\t100\t2.5\t500000\nENST1.2_PAR_Y|ENSG1.3_PAR_Y|label|\t100\t3.5\t500000\n";
        let data = import(&[a, a], "ENST1\tENSG1\nENST1_PAR_Y\tENSG1_PAR_Y\n", "no").unwrap();
        assert_eq!(data.gene_ids, ["ENSG1", "ENSG1_PAR_Y"]);
        assert_eq!(data.counts, [2, 2, 4, 4]);
        assert_eq!(data.diagnostics[0].mapped_transcripts, 2);
    }

    #[test]
    fn rejects_ambiguous_or_inconsistent_inputs() {
        let a = "t1\t100\t100\t1000000\n";
        for (rows, map, mode) in [
            (vec![a, a], "t1\tG1\nt1\tG2\n", "no"),
            (vec![a, "t2\t100\t100\t1000000\n"], "t1\tG1\nt2\tG2\n", "no"),
            (vec![a, a], "t1\tG1\n", "dtuScaledTPM"),
            (
                vec!["t1.1\t100\t1\t1\nt1.2\t100\t1\t1\n"; 2],
                "t1\tG1\n",
                "no",
            ),
        ] {
            assert!(import(&rows, map, mode).is_err());
        }
    }

    #[test]
    fn offsets_estimate_library_sizes_after_length_correction() {
        let mut factors = vec![0.5, 2.0, 0.5, 2.0];
        let sizes = complete_length_offsets(&[100, 400, 200, 800], &mut factors, 2).unwrap();
        for value in sizes {
            assert!((value - 1.0).abs() < 1e-12);
        }
        for (actual, expected) in factors.iter().zip([0.5, 2.0, 0.5, 2.0]) {
            assert!((actual - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn offset_estimation_rejects_sparse_or_invalid_matrices() {
        assert!(complete_length_offsets(&[0, 1, 1, 0], &mut [1.0; 4], 2).is_err());
        assert!(complete_length_offsets(&[1, 1], &mut [0.0, 1.0], 2).is_err());
        assert!(complete_length_offsets(&[1], &mut [1.0], 2).is_err());
    }
}
