use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail};
use mincorr::{bicor, hellcor, kendall, pearson, spearman};
use nalgebra::{DMatrix, DVector};
use ndarray::Array2;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{AdjustedMatrix, CorrelationMethod, OutputMode};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::{Correlator, MatrixAdjuster};
use rnaa_core::util::{file_size, now_rfc3339, write_json_pretty};
use rnaa_formats::matrix::{read_matrix_tsv, write_matrix_tsv};
use serde_json::json;

#[derive(Debug, Clone, Default)]
pub struct Residualizer;

impl MatrixAdjuster for Residualizer {
    fn adjust(
        &self,
        matrix_path: &Path,
        samplesheet_path: &Path,
        model: &str,
        paths: &ProjectPaths,
    ) -> Result<AdjustedMatrix> {
        let matrix = read_matrix_tsv(matrix_path)?;
        if matrix.col_ids.is_empty() {
            bail!("input matrix has no samples: {}", matrix_path.display());
        }
        let sample_metadata = load_samplesheet(samplesheet_path)?;
        let design_columns = parse_formula_columns(model);
        let design = build_design_matrix(&matrix.col_ids, &sample_metadata, &design_columns)?;

        let adjusted_values = matrix
            .values
            .iter()
            .map(|row| residualize_row(row, &design))
            .collect::<Vec<_>>();
        let adjusted = rnaa_formats::matrix::MatrixData {
            row_ids: matrix.row_ids.clone(),
            col_ids: matrix.col_ids.clone(),
            values: adjusted_values,
        };
        let out_path = paths.corr_dir.join("adjusted_matrix.tsv");
        write_matrix_tsv(&out_path, &adjusted)?;

        Ok(AdjustedMatrix {
            matrix_path: out_path,
            genes: adjusted.row_ids,
            samples: adjusted.col_ids,
        })
    }
}

#[derive(Debug, Clone, Default)]
pub struct MinCorrCorrelator;

impl Correlator for MinCorrCorrelator {
    fn correlate(
        &self,
        matrix: &AdjustedMatrix,
        method: CorrelationMethod,
        output_mode: &OutputMode,
        paths: &ProjectPaths,
        project_id: &str,
    ) -> Result<Vec<std::path::PathBuf>> {
        let started_at = now_rfc3339();
        let matrix_data = read_matrix_tsv(&matrix.matrix_path)?;
        if matrix_data.row_ids.len() < 2 {
            bail!(
                "need at least 2 genes for correlation, got {}",
                matrix_data.row_ids.len()
            );
        }

        let n_rows = matrix_data.row_ids.len();
        let n_cols = matrix_data.col_ids.len();
        let flat = matrix_data
            .values
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let data = Array2::<f64>::from_shape_vec((n_rows, n_cols), flat)
            .context("matrix shape mismatch while preparing correlation input")?;

        let corr = match method {
            CorrelationMethod::Pearson => pearson::correlation_matrix(&data),
            CorrelationMethod::Spearman => spearman::correlation_matrix(&data),
            CorrelationMethod::Kendall => kendall::correlation_matrix(&data),
            CorrelationMethod::Bicor => bicor::correlation_matrix(&data),
            CorrelationMethod::Hellcor => hellcor::correlation_matrix(&data),
        };

        let corr_dir = paths.corr_project_dir(project_id);
        fs::create_dir_all(&corr_dir)
            .with_context(|| format!("failed to create {}", corr_dir.display()))?;

        let mut outputs = Vec::new();
        match output_mode {
            OutputMode::TopK { k } => {
                let out_path = corr_dir.join("edges_topk.tsv");
                write_topk_edges(&out_path, &matrix_data.row_ids, &corr, *k)?;
                outputs.push(out_path);
            }
            OutputMode::Threshold { min_abs_r } => {
                let out_path = corr_dir.join("edges_threshold.tsv");
                write_threshold_edges(&out_path, &matrix_data.row_ids, &corr, *min_abs_r)?;
                outputs.push(out_path);
            }
            OutputMode::Dense => {
                let out_path = corr_dir.join("corr_dense.tsv");
                write_dense(&out_path, &matrix_data.row_ids, &corr)?;
                outputs.push(out_path);
            }
        }

        let manifest_path = paths.manifest_path("corr", project_id);
        let output_artifacts = outputs
            .iter()
            .map(|path| ManifestArtifact {
                kind: "CORR_OUTPUT".to_string(),
                path: path.display().to_string(),
                checksum_type: "none".to_string(),
                checksum: String::new(),
                bytes: file_size(path).unwrap_or_default(),
            })
            .collect::<Vec<_>>();
        let manifest = StageManifest {
            schema_version: 1,
            stage: "corr".to_string(),
            id: project_id.to_string(),
            started_at,
            finished_at: now_rfc3339(),
            parameters: json!({
                "method": method.to_string(),
                "output_mode": output_mode
            }),
            tool_versions: BTreeMap::new(),
            input_artifacts: vec![ManifestArtifact {
                kind: "ADJUSTED_MATRIX".to_string(),
                path: matrix.matrix_path.display().to_string(),
                checksum_type: "none".to_string(),
                checksum: String::new(),
                bytes: file_size(&matrix.matrix_path).unwrap_or_default(),
            }],
            output_artifacts,
            notes: vec!["Correlation matrix computed via mincorr crate integration.".to_string()],
        };
        write_json_pretty(&manifest_path, &manifest)?;
        outputs.push(manifest_path);
        Ok(outputs)
    }
}

fn load_samplesheet(path: &Path) -> Result<BTreeMap<String, BTreeMap<String, String>>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read samplesheet {}", path.display()))?;
    let headers = reader
        .headers()
        .with_context(|| format!("failed to parse samplesheet header {}", path.display()))?
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let run_idx = headers
        .iter()
        .position(|name| name == "run_accession")
        .ok_or_else(|| anyhow::anyhow!("samplesheet missing run_accession column"))?;
    let mut records = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let run = row
            .get(run_idx)
            .ok_or_else(|| anyhow::anyhow!("invalid samplesheet row: missing run_accession"))?;
        let mut map = BTreeMap::new();
        for (idx, name) in headers.iter().enumerate() {
            if let Some(value) = row.get(idx) {
                map.insert(name.clone(), value.to_string());
            }
        }
        records.insert(run.to_string(), map);
    }
    Ok(records)
}

fn parse_formula_columns(model: &str) -> Vec<String> {
    let mut columns = model
        .split(|ch: char| !ch.is_ascii_alphanumeric() && ch != '_')
        .filter(|token| !token.is_empty())
        .filter(|token| !token.chars().all(|ch| ch.is_ascii_digit()))
        .filter(|token| *token != "Intercept")
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    columns.sort();
    columns.dedup();
    columns
}

fn build_design_matrix(
    sample_ids: &[String],
    sample_metadata: &BTreeMap<String, BTreeMap<String, String>>,
    design_columns: &[String],
) -> Result<DMatrix<f64>> {
    let mut columns = Vec::<Vec<f64>>::new();
    columns.push(vec![1.0; sample_ids.len()]);

    for design_col in design_columns {
        let mut values = Vec::with_capacity(sample_ids.len());
        for sample_id in sample_ids {
            let metadata = sample_metadata.get(sample_id).ok_or_else(|| {
                anyhow::anyhow!("sample '{}' not found in samplesheet metadata", sample_id)
            })?;
            let value = metadata
                .get(design_col)
                .ok_or_else(|| anyhow::anyhow!("missing '{}' in samplesheet", design_col))?
                .to_string();
            values.push(value);
        }

        let numeric = values
            .iter()
            .map(|value| value.parse::<f64>())
            .collect::<std::result::Result<Vec<_>, _>>();
        if let Ok(numeric) = numeric {
            columns.push(numeric);
            continue;
        }

        let levels = values.iter().cloned().collect::<BTreeSet<_>>();
        let level_values = levels.into_iter().collect::<Vec<_>>();
        if level_values.len() <= 1 {
            continue;
        }
        for level in level_values.iter().skip(1) {
            let one_hot = values
                .iter()
                .map(|value| if value == level { 1.0 } else { 0.0 })
                .collect::<Vec<_>>();
            columns.push(one_hot);
        }
    }

    let n = sample_ids.len();
    let p = columns.len();
    let mut data = Vec::with_capacity(n * p);
    for row_idx in 0..n {
        for column in &columns {
            data.push(column[row_idx]);
        }
    }
    Ok(DMatrix::<f64>::from_row_slice(n, p, &data))
}

fn residualize_row(row: &[f64], design: &DMatrix<f64>) -> Vec<f64> {
    if row.is_empty() || design.nrows() != row.len() {
        return row.to_vec();
    }
    let y = DVector::<f64>::from_column_slice(row);
    if let Some(beta) = design.clone().qr().solve(&y) {
        let fitted = design * beta;
        let residual = y - fitted;
        residual.iter().copied().collect()
    } else {
        let mean = row.iter().sum::<f64>() / row.len() as f64;
        row.iter().map(|value| value - mean).collect()
    }
}

fn write_topk_edges(path: &Path, genes: &[String], corr: &Array2<f64>, k: usize) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    writer.write_record(["gene_a", "gene_b", "r"])?;
    for i in 0..genes.len() {
        let mut row = (0..genes.len())
            .filter(|&j| j != i)
            .map(|j| (j, corr[[i, j]]))
            .filter(|(_, r)| r.is_finite())
            .collect::<Vec<_>>();
        row.sort_by(|left, right| right.1.abs().total_cmp(&left.1.abs()));
        for (j, r) in row.into_iter().take(k) {
            writer.write_record([genes[i].as_str(), genes[j].as_str(), &r.to_string()])?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_threshold_edges(
    path: &Path,
    genes: &[String],
    corr: &Array2<f64>,
    min_abs_r: f64,
) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    writer.write_record(["gene_a", "gene_b", "r"])?;
    for i in 0..genes.len() {
        for j in (i + 1)..genes.len() {
            let r = corr[[i, j]];
            if r.is_finite() && r.abs() >= min_abs_r {
                writer.write_record([genes[i].as_str(), genes[j].as_str(), &r.to_string()])?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_dense(path: &Path, genes: &[String], corr: &Array2<f64>) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut header = vec!["gene_id".to_string()];
    header.extend(genes.iter().cloned());
    writer.write_record(header)?;
    for i in 0..genes.len() {
        let mut record = vec![genes[i].clone()];
        record.extend((0..genes.len()).map(|j| corr[[i, j]].to_string()));
        writer.write_record(record)?;
    }
    writer.flush()?;
    Ok(())
}
