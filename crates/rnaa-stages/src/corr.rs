use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

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
use rsfgsea::prelude::{RankedList, ScoreType, read_gmt, run_gsea};
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

#[derive(Debug, Clone)]
pub struct ModuleGseaParams {
    pub gmt_path: PathBuf,
    pub top_genes: usize,
    pub min_module_size: usize,
    pub permutations: usize,
    pub seed: u64,
}

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
    let svd = design.clone().svd(true, true);
    if let Ok(beta) = svd.solve(&y, 1e-12) {
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

pub fn run_module_gsea(
    matrix: &AdjustedMatrix,
    paths: &ProjectPaths,
    project_id: &str,
    params: &ModuleGseaParams,
) -> Result<Vec<PathBuf>> {
    let started_at = now_rfc3339();
    let matrix_data = read_matrix_tsv(&matrix.matrix_path)?;
    if matrix_data.row_ids.len() < 2 {
        bail!("module-gsea requires at least 2 genes");
    }
    if matrix_data.col_ids.len() < 3 {
        bail!("module-gsea requires at least 3 samples for t-statistics");
    }

    let n_rows = matrix_data.row_ids.len();
    let n_cols = matrix_data.col_ids.len();
    let flat = matrix_data
        .values
        .iter()
        .flat_map(|row| row.iter().copied())
        .collect::<Vec<_>>();
    let data = Array2::<f64>::from_shape_vec((n_rows, n_cols), flat)
        .context("matrix shape mismatch while preparing module-gsea")?;
    let corr = spearman::correlation_matrix(&data);

    let modules = detect_modules(
        &corr,
        params.min_module_size.max(2),
        params.top_genes.max(2).min(n_rows),
    );
    if modules.is_empty() {
        bail!(
            "no modules detected with min size {}; try reducing --module-min-size",
            params.min_module_size
        );
    }

    let pathways = read_gmt(params.gmt_path.to_string_lossy().as_ref()).with_context(|| {
        format!(
            "failed to parse module gmt at {}",
            params.gmt_path.display()
        )
    })?;
    if pathways.pathways.is_empty() {
        bail!("module gmt has no pathways: {}", params.gmt_path.display());
    }

    let corr_dir = paths.corr_project_dir(project_id);
    fs::create_dir_all(&corr_dir)
        .with_context(|| format!("failed to create {}", corr_dir.display()))?;
    let module_table_path = corr_dir.join("modules_top_genes.tsv");
    let t_values_path = corr_dir.join("module_t_values.tsv");
    let fgsea_path = corr_dir.join("module_fgsea.tsv");
    let manifest_path = corr_dir.join("module_fgsea_manifest.json");

    let mut module_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&module_table_path)
        .with_context(|| format!("failed to create {}", module_table_path.display()))?;
    module_writer.write_record([
        "module_id",
        "module_size",
        "gene",
        "rank_in_module",
        "connectivity",
    ])?;

    let mut t_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&t_values_path)
        .with_context(|| format!("failed to create {}", t_values_path.display()))?;
    t_writer.write_record(["module_id", "gene", "t_value"])?;

    let mut fgsea_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&fgsea_path)
        .with_context(|| format!("failed to create {}", fgsea_path.display()))?;
    fgsea_writer.write_record([
        "module_id",
        "module_size",
        "pathway",
        "pathway_size",
        "es",
        "nes",
        "p_value",
        "padj",
        "log2err",
        "leading_edge_size",
    ])?;

    for (module_idx, module) in modules.iter().enumerate() {
        let module_id = format!("M{:04}", module_idx + 1);
        let eigengene = module_eigengene(&matrix_data.values, module);
        let mut ranked = Vec::with_capacity(matrix_data.row_ids.len());
        for (gene_idx, gene) in matrix_data.row_ids.iter().enumerate() {
            let r = spearman_corr(&matrix_data.values[gene_idx], &eigengene);
            let t = corr_to_t_stat(r, matrix_data.col_ids.len());
            ranked.push((gene.clone(), t));
            t_writer.write_record([module_id.as_str(), gene.as_str(), &t.to_string()])?;
        }

        let genes = ranked
            .iter()
            .map(|(gene, _)| gene.clone())
            .collect::<Vec<_>>();
        let scores = ranked.iter().map(|(_, score)| *score).collect::<Vec<_>>();
        let ranked_list = RankedList::new(genes, scores);
        let max_size = ranked_list.len().saturating_sub(1).max(2);
        let min_size = 10.min(max_size);
        let gsea_results = run_gsea(
            &ranked_list,
            &pathways.pathways,
            params.permutations.max(100),
            params.seed,
            min_size,
            max_size,
            1e-50,
            ScoreType::Std,
            1.0,
        );

        for (rank, &(gene_idx, connectivity)) in module.iter().enumerate() {
            module_writer.write_record([
                module_id.as_str(),
                &module.len().to_string(),
                matrix_data.row_ids[gene_idx].as_str(),
                &(rank + 1).to_string(),
                &connectivity.to_string(),
            ])?;
        }
        for item in gsea_results {
            fgsea_writer.write_record([
                module_id.as_str(),
                &module.len().to_string(),
                item.pathway_name.as_str(),
                &item.size.to_string(),
                &item.es.to_string(),
                &item.nes.unwrap_or(f64::NAN).to_string(),
                &item.p_value.to_string(),
                &item.padj.unwrap_or(f64::NAN).to_string(),
                &item.log2err.unwrap_or(f64::NAN).to_string(),
                &item.leading_edge.len().to_string(),
            ])?;
        }
    }
    module_writer.flush()?;
    t_writer.flush()?;
    fgsea_writer.flush()?;

    let outputs = vec![
        module_table_path.clone(),
        t_values_path.clone(),
        fgsea_path.clone(),
    ];
    let manifest = StageManifest {
        schema_version: 1,
        stage: "module_gsea".to_string(),
        id: project_id.to_string(),
        started_at,
        finished_at: now_rfc3339(),
        parameters: json!({
            "method": "spearman",
            "top_genes_per_module": params.top_genes,
            "min_module_size": params.min_module_size,
            "permutations": params.permutations,
            "seed": params.seed,
            "gmt": params.gmt_path.display().to_string(),
        }),
        tool_versions: BTreeMap::new(),
        input_artifacts: vec![ManifestArtifact {
            kind: "ADJUSTED_MATRIX".to_string(),
            path: matrix.matrix_path.display().to_string(),
            checksum_type: "none".to_string(),
            checksum: String::new(),
            bytes: file_size(&matrix.matrix_path).unwrap_or_default(),
        }],
        output_artifacts: outputs
            .iter()
            .map(|path| ManifestArtifact {
                kind: "MODULE_GSEA_OUTPUT".to_string(),
                path: path.display().to_string(),
                checksum_type: "none".to_string(),
                checksum: String::new(),
                bytes: file_size(path).unwrap_or_default(),
            })
            .collect(),
        notes: vec![
            "Modules are connected components from mutual Spearman kNN graph.".to_string(),
            "Per-module ranking scores are Spearman-correlation-derived t-values.".to_string(),
        ],
    };
    write_json_pretty(&manifest_path, &manifest)?;

    let mut all_outputs = outputs;
    all_outputs.push(manifest_path);
    Ok(all_outputs)
}

fn detect_modules(corr: &Array2<f64>, min_size: usize, top_genes: usize) -> Vec<Vec<(usize, f64)>> {
    let n = corr.shape()[0];
    let neighbor_k = 20.min(n.saturating_sub(1)).max(2);
    let mut top_neighbors = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, corr[[i, j]]))
            .filter(|(_, r)| r.is_finite() && *r > 0.0)
            .collect::<Vec<_>>();
        row.sort_by(|left, right| right.1.total_cmp(&left.1));
        top_neighbors.push(row.into_iter().take(neighbor_k).collect::<Vec<_>>());
    }

    let mut adjacency = vec![Vec::<usize>::new(); n];
    let neighbor_sets = top_neighbors
        .iter()
        .map(|items| items.iter().map(|(idx, _)| *idx).collect::<HashSet<_>>())
        .collect::<Vec<_>>();
    for i in 0..n {
        for (j, _) in &top_neighbors[i] {
            if neighbor_sets[*j].contains(&i) {
                adjacency[i].push(*j);
                adjacency[*j].push(i);
            }
        }
    }

    let mut visited = vec![false; n];
    let mut modules = Vec::new();
    for i in 0..n {
        if visited[i] {
            continue;
        }
        let mut stack = vec![i];
        visited[i] = true;
        let mut component = Vec::new();
        while let Some(node) = stack.pop() {
            component.push(node);
            for &next in &adjacency[node] {
                if !visited[next] {
                    visited[next] = true;
                    stack.push(next);
                }
            }
        }
        if component.len() < min_size {
            continue;
        }
        let mut ranked = component
            .iter()
            .map(|&gene_idx| {
                let connectivity = component
                    .iter()
                    .filter(|&&other| other != gene_idx)
                    .map(|&other| corr[[gene_idx, other]].max(0.0))
                    .sum::<f64>();
                (gene_idx, connectivity)
            })
            .collect::<Vec<_>>();
        ranked.sort_by(|left, right| right.1.total_cmp(&left.1));
        ranked.truncate(top_genes.min(ranked.len()));
        modules.push(ranked);
    }
    modules.sort_by_key(|module| std::cmp::Reverse(module.len()));
    modules
}

fn module_eigengene(values: &[Vec<f64>], module: &[(usize, f64)]) -> Vec<f64> {
    if module.is_empty() {
        return Vec::new();
    }
    let sample_count = values[module[0].0].len();
    let mut eigengene = vec![0.0_f64; sample_count];
    for &(gene_idx, _) in module {
        for (sample_idx, value) in values[gene_idx].iter().enumerate() {
            eigengene[sample_idx] += *value;
        }
    }
    let denom = module.len() as f64;
    for item in &mut eigengene {
        *item /= denom;
    }
    eigengene
}

fn rank_with_ties(values: &[f64]) -> Vec<f64> {
    let mut indexed = values
        .iter()
        .enumerate()
        .map(|(idx, value)| (idx, *value))
        .collect::<Vec<_>>();
    indexed.sort_by(|left, right| left.1.total_cmp(&right.1));
    let mut ranks = vec![0.0_f64; values.len()];
    let mut i = 0_usize;
    while i < indexed.len() {
        let mut j = i + 1;
        while j < indexed.len() && indexed[j].1.total_cmp(&indexed[i].1) == Ordering::Equal {
            j += 1;
        }
        let avg_rank = (i + j - 1) as f64 / 2.0 + 1.0;
        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }
        i = j;
    }
    ranks
}

fn pearson_corr(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return 0.0;
    }
    let mean_x = x.iter().sum::<f64>() / x.len() as f64;
    let mean_y = y.iter().sum::<f64>() / y.len() as f64;
    let mut num = 0.0;
    let mut den_x = 0.0;
    let mut den_y = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        num += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
    }
    let den = (den_x * den_y).sqrt();
    if den == 0.0 { 0.0 } else { num / den }
}

fn spearman_corr(x: &[f64], y: &[f64]) -> f64 {
    let rx = rank_with_ties(x);
    let ry = rank_with_ties(y);
    pearson_corr(&rx, &ry)
}

fn corr_to_t_stat(r: f64, n: usize) -> f64 {
    if n <= 2 || !r.is_finite() {
        return 0.0;
    }
    let capped = r.clamp(-0.999_999_9, 0.999_999_9);
    capped * (((n - 2) as f64) / (1.0 - capped * capped)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residualize_row_handles_non_square_design() {
        let design = DMatrix::<f64>::from_row_slice(
            3,
            4,
            &[
                1.0, 0.0, 0.0, 1.0, //
                1.0, 1.0, 0.0, 0.0, //
                1.0, 0.0, 1.0, 0.0,
            ],
        );
        let row = vec![1.0, 2.0, 3.0];
        let residual = residualize_row(&row, &design);
        assert_eq!(residual.len(), row.len());
        assert!(residual.iter().all(|value| value.is_finite()));
    }
}
