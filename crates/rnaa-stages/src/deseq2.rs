use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use csv::{ReaderBuilder, WriterBuilder};
use rnaa_core::config::ProjectConfig;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{ContrastSpec, ReferenceBundle};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::DifferentialExpression;
use rnaa_core::util::{compute_sha256, file_size, now_rfc3339, write_json_pretty};
use rsdeseq2::prelude::{CountMatrix, DeseqBuilder, RowMajorMatrix};
use serde_json::json;

use crate::tximport::{
    IMPORT_ALGORITHM_VERSION, QuantifiedSample, RunAbundance, complete_length_offsets,
    import_grouped_kallisto_abundances,
};
mod fitting;
mod results;
/// Fitting policy changes invalidate normalization/DE requests as well as results.
pub const DE_ALGORITHM_VERSION: u32 = 7;
pub const DISPERSION_FIT_POLICY: &str = "parametric_with_deseq2_local_fallback_v1";

#[derive(Debug, Clone, Default)]
pub struct RustDeseq2Runner;

impl DifferentialExpression for RustDeseq2Runner {
    fn deseq2(
        &self,
        project_id: &str,
        samplesheet_path: &Path,
        reference: &ReferenceBundle,
        design: &str,
        contrasts: &[ContrastSpec],
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<Vec<PathBuf>> {
        if !samplesheet_path.exists() {
            bail!(
                "samplesheet missing at {}; run `rnaa resolve` first",
                samplesheet_path.display()
            );
        }
        let started_at = now_rfc3339();
        let de_dir = paths.de_project_dir(project_id);
        fs::create_dir_all(&de_dir)
            .with_context(|| format!("failed to create {}", de_dir.display()))?;
        let manifest_path = de_dir.join("deseq2_manifest.json");
        let sample_plan = rnaa_core::sample_plan::SamplePlan::from_sheet_with_metadata(
            samplesheet_path,
            design,
            contrasts,
            &config.metadata,
        )?;
        let (headers, rows) = sample_plan.model_table()?;
        let prepared = crate::design::PreparedDesign::new(
            &headers,
            &rows,
            &sample_plan.sample_ids(),
            design,
            &config.metadata,
            contrasts,
            crate::design::Purpose::DifferentialExpression,
        )?;
        let design_path = de_dir.join("model_design.json");
        write_json_pretty(&design_path, &prepared.contract)?;
        let groups = sample_plan.samples.iter().map(|sample| {
            let runs = sample.runs.iter().map(|run| {
                let path = paths.quant_run_dir(run, &config.quant.engine).join("abundance.tsv");
                if !path.is_file() || file_size(&path).unwrap_or_default() == 0 { bail!("kallisto abundance table is missing or empty for run '{run}' at {}", path.display()); }
                Ok(RunAbundance { run_id: run.clone(), path })
            }).collect::<Result<_>>()?;
            Ok(QuantifiedSample { sample_id: sample.sample_id.clone(), runs })
        }).collect::<Result<Vec<_>>>()?;
        let abundance_paths = groups
            .iter()
            .flat_map(|sample| &sample.runs)
            .map(|run| (run.run_id.clone(), run.path.clone()))
            .collect::<Vec<_>>();
        let imported = import_grouped_kallisto_abundances(
            &groups,
            &reference.tx2gene_path,
            &config.deseq2.counts_from_abundance,
        )?;
        let sample_metadata_path = de_dir.join("sample_metadata.tsv");
        sample_plan.write_sheet(&sample_metadata_path)?;
        let counts = CountMatrix::from_row_major_u32_with_names(
            imported.gene_ids.len(),
            imported.sample_ids.len(),
            imported.counts.clone(),
            Some(imported.gene_ids.clone()),
            Some(imported.sample_ids.clone()),
        )?;
        let mut builder = DeseqBuilder::new().threads(1);
        let mut imported_size_factors = None;
        if let Some(mut values) = imported.normalization_factors {
            // DESeq2 combines tximport's gene/sample length correction with
            // sample size factors; rsdeseq2 correctly treats a supplied
            // normalization matrix as the complete offset matrix.
            let size_factors =
                complete_length_offsets(&imported.counts, &mut values, imported.sample_ids.len())?;
            let factors = RowMajorMatrix::from_row_major(
                imported.gene_ids.len(),
                imported.sample_ids.len(),
                values,
            )?;
            builder = builder.normalization_factors(factors);
            imported_size_factors = Some(size_factors);
        }
        let builder = builder.try_model_frame(prepared.frame.clone())?;
        let (fit, dispersion_fit) = fitting::fit_formula(&builder, &counts, design)
            .context("rsdeseq2 failed to fit the Rust DESeq2 workflow")?;
        let normalized = fit.normalized_counts(&counts)?;
        let transformed = match config.deseq2.transform.as_str() {
            "vst" => fit.vst(&counts)?,
            "rlog" => fit.rlog(&counts)?.transformed,
            other => bail!("unsupported DE transform '{other}'; expected 'vst' or 'rlog'"),
        };
        let raw_matrix = matrix_f64(
            &imported.gene_ids,
            &imported.sample_ids,
            &imported
                .counts
                .iter()
                .map(|value| *value as f64)
                .collect::<Vec<_>>(),
        );
        let norm_matrix = matrix_f64(
            &imported.gene_ids,
            &imported.sample_ids,
            normalized.as_slice(),
        );
        let vst_matrix = matrix_f64(
            &imported.gene_ids,
            &imported.sample_ids,
            transformed.as_slice(),
        );
        let raw_path = de_dir.join("gene_counts.tsv");
        let norm_path = de_dir.join("gene_norm_counts.tsv");
        let vst_path = de_dir.join("vst.tsv");
        write_matrix(&raw_path, &raw_matrix)?;
        write_matrix(&norm_path, &norm_matrix)?;
        write_matrix(&vst_path, &vst_matrix)?;
        let size_factors_path = de_dir.join("size_factors.tsv");
        write_size_factors(
            &size_factors_path,
            &imported.sample_ids,
            imported_size_factors
                .as_deref()
                .unwrap_or(&fit.size_factors),
        )?;

        let annotation = reference
            .gene_annotation_path
            .as_ref()
            .map(|path| read_annotation(path))
            .transpose()?
            .unwrap_or_default();
        let annotated_raw = de_dir.join("gene_counts_annotated.tsv");
        let annotated_norm = de_dir.join("gene_norm_counts_annotated.tsv");
        let annotated_vst = de_dir.join("vst_annotated.tsv");
        write_annotated_matrix(&annotated_raw, &raw_matrix, &annotation)?;
        write_annotated_matrix(&annotated_norm, &norm_matrix, &annotation)?;
        write_annotated_matrix(&annotated_vst, &vst_matrix, &annotation)?;

        let mut outputs = vec![
            sample_metadata_path,
            design_path,
            raw_path,
            annotated_raw,
            norm_path,
            annotated_norm,
            vst_path,
            annotated_vst,
            size_factors_path,
        ];
        let mut result_policies = BTreeMap::new();
        for contrast in &prepared.contract.contrasts {
            let name = &contrast.name;
            let mut results = fit
                .wald_results_factor_level_contrast(
                    &counts,
                    rsdeseq2::contrasts::FactorLevelContrast {
                        factor: &contrast.factor,
                        numerator: &contrast.numerator,
                        denominator: &contrast.denominator,
                        // Infer from actual coding, retaining full indicators and
                        // releveling. Explicit levels preserve character-contrast
                        // all-zero handling instead of a generic numeric shortcut.
                        reference: None,
                        sample_levels: &contrast.sample_levels,
                    },
                    Some(&imported.gene_ids),
                )
                .with_context(|| format!("failed to assemble contrast {name}"))?;
            result_policies.insert(
                name,
                results::finalize(&mut results, &fit, &counts, &prepared.contract)?,
            );
            let result_path = de_dir.join(format!("de_{name}.tsv"));
            write_results(&result_path, &results.rows)?;
            let annotated_result_path = de_dir.join(format!("de_{name}_annotated.tsv"));
            write_annotated_results(&annotated_result_path, &results.rows, &annotation)?;
            outputs.push(result_path);
            outputs.push(annotated_result_path);
        }
        outputs.push(manifest_path.clone());
        let mut input_paths = vec![
            ("SAMPLESHEET", samplesheet_path),
            ("TX2GENE", reference.tx2gene_path.as_path()),
        ];
        if let Some(annotation) = &reference.gene_annotation_path {
            input_paths.push(("GENE_ANNOTATION", annotation));
        }
        input_paths.extend(
            abundance_paths
                .iter()
                .map(|(_, path)| ("ABUNDANCE_TSV", path.as_path())),
        );
        let input_artifacts = input_paths
            .into_iter()
            .map(|(kind, path)| {
                Ok(ManifestArtifact {
                    kind: kind.to_string(),
                    path: path.display().to_string(),
                    checksum_type: "sha256".into(),
                    checksum: compute_sha256(path)?,
                    bytes: file_size(path)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let manifest = StageManifest {
            schema_version: 1,
            stage: "deseq2".to_string(),
            id: project_id.to_string(),
            started_at,
            finished_at: now_rfc3339(),
            parameters: json!({
                "design": design,
                "transform": config.deseq2.transform,
                "counts_from_abundance": config.deseq2.counts_from_abundance,
                "sample_plan": sample_plan,
                "pooling_policy": "fractional_transcript_mass_v1",
                "quant_engine": "kallistors",
                "de_engine": "rsdeseq2",
                "de_engine_version": crate::engines::RSDESEQ2,
                "import_algorithm_version": IMPORT_ALGORITHM_VERSION,
                "de_algorithm_version":DE_ALGORITHM_VERSION,
                "dispersion_fit_policy":DISPERSION_FIT_POLICY,
                "dispersion_fit":dispersion_fit,
                "result_policy":"deseq2_cooks_then_independent_filtering_v1",
                "result_policies":result_policies,
                "import_diagnostics": imported.diagnostics,
                "samples": imported.sample_ids,
                "genes": imported.gene_ids.len(),
                "contrasts": contrasts.iter().map(|contrast| format!("{}:{}-{}", contrast.factor, contrast.level_a, contrast.level_b)).collect::<Vec<_>>(),
                "input_checksums": {
                    "samplesheet_sha256": compute_sha256(samplesheet_path)?,
                    "tx2gene_sha256": compute_sha256(&reference.tx2gene_path)?,
                    "abundance_sha256": abundance_paths.iter().map(|(run, path)| Ok((run.clone(), compute_sha256(path)?))).collect::<Result<BTreeMap<_, _>>>()?
                }
            }),
            tool_versions: BTreeMap::from([
                ("rsdeseq2".to_string(), crate::engines::RSDESEQ2.to_string()),
                ("kallistors".to_string(), crate::engines::KALLISTORS.to_string()),
            ]),
            input_artifacts,
            output_artifacts: outputs.iter().filter(|path| path.as_path() != manifest_path.as_path()).map(|path| Ok(ManifestArtifact {
                kind: "DE_OUTPUT".to_string(), path: path.display().to_string(),
                checksum_type: "sha256".to_string(), checksum: compute_sha256(path)?,
                bytes: file_size(path).unwrap_or_default(),
            })).collect::<Result<Vec<_>>>()?,
            notes: vec!["Native Rust transcript-to-gene import and rsdeseq2 fitting; no R runtime or Bioconductor packages required.".to_string()],
        };
        write_json_pretty(&manifest_path, &manifest)?;
        Ok(outputs)
    }
}

struct Table {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn read_table(path: &Path) -> Result<Table> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to open {}", path.display()))?;
    let headers = reader.headers()?.iter().map(str::to_string).collect();
    let rows = reader
        .records()
        .map(|row| Ok(row?.iter().map(str::to_string).collect()))
        .collect::<Result<Vec<Vec<String>>>>()?;
    Ok(Table { headers, rows })
}

fn matrix_f64(
    rows: &[String],
    cols: &[String],
    values: &[f64],
) -> rnaa_formats::matrix::MatrixData {
    rnaa_formats::matrix::MatrixData {
        row_ids: rows.to_vec(),
        col_ids: cols.to_vec(),
        values: values.chunks(cols.len()).map(|row| row.to_vec()).collect(),
    }
}
fn write_matrix(path: &Path, matrix: &rnaa_formats::matrix::MatrixData) -> Result<()> {
    rnaa_formats::matrix::write_matrix_tsv(path, matrix)
}
fn write_size_factors(path: &Path, samples: &[String], values: &[f64]) -> Result<()> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record(["sample", "size_factor"])?;
    for (sample, value) in samples.iter().zip(values) {
        writer.write_record([sample, &value.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}
fn read_annotation(path: &Path) -> Result<HashMap<String, (String, String)>> {
    let table = read_table(path)?;
    let gene_id = table
        .headers
        .iter()
        .position(|header| header == "gene_id")
        .ok_or_else(|| anyhow::anyhow!("gene annotation is missing gene_id"))?;
    let gene_name = table
        .headers
        .iter()
        .position(|header| header == "gene_name");
    let biotype = table
        .headers
        .iter()
        .position(|header| header == "gene_biotype");
    Ok(table
        .rows
        .iter()
        .map(|row| {
            (
                row.get(gene_id).cloned().unwrap_or_default(),
                (
                    gene_name
                        .and_then(|idx| row.get(idx))
                        .cloned()
                        .unwrap_or_default(),
                    biotype
                        .and_then(|idx| row.get(idx))
                        .cloned()
                        .unwrap_or_default(),
                ),
            )
        })
        .collect())
}
fn write_annotated_matrix(
    path: &Path,
    matrix: &rnaa_formats::matrix::MatrixData,
    annotation: &HashMap<String, (String, String)>,
) -> Result<()> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    let mut headers = vec![
        "gene_id".to_string(),
        "gene_symbol".to_string(),
        "gene_biotype".to_string(),
    ];
    headers.extend(matrix.col_ids.clone());
    writer.write_record(headers)?;
    for (gene, values) in matrix.row_ids.iter().zip(&matrix.values) {
        let (symbol, biotype) = annotation.get(gene).cloned().unwrap_or_default();
        let mut row = vec![gene.clone(), symbol, biotype];
        row.extend(values.iter().map(ToString::to_string));
        writer.write_record(row)?;
    }
    writer.flush()?;
    Ok(())
}
fn write_results(path: &Path, rows: &[rsdeseq2::prelude::DeseqResultRow]) -> Result<()> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record([
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    ])?;
    for row in rows {
        writer.write_record([
            row.gene.as_deref().unwrap_or(""),
            &row.base_mean.to_string(),
            &fmt_optional(row.log2_fold_change),
            &fmt_optional(row.lfc_se),
            &fmt_optional(row.stat),
            &fmt_optional(row.pvalue),
            &fmt_optional(row.padj),
        ])?;
    }
    writer.flush()?;
    Ok(())
}
fn write_annotated_results(
    path: &Path,
    rows: &[rsdeseq2::prelude::DeseqResultRow],
    annotation: &HashMap<String, (String, String)>,
) -> Result<()> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record([
        "gene_id",
        "gene_symbol",
        "gene_biotype",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    ])?;
    for row in rows {
        let gene = row.gene.as_deref().unwrap_or("");
        let (symbol, biotype) = annotation.get(gene).cloned().unwrap_or_default();
        writer.write_record([
            gene,
            &symbol,
            &biotype,
            &row.base_mean.to_string(),
            &fmt_optional(row.log2_fold_change),
            &fmt_optional(row.lfc_se),
            &fmt_optional(row.stat),
            &fmt_optional(row.pvalue),
            &fmt_optional(row.padj),
        ])?;
    }
    writer.flush()?;
    Ok(())
}
fn fmt_optional(value: Option<f64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}
