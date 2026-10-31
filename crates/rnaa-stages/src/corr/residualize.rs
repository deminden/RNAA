//! Design-aware matrix adjustment and residualization.
use anyhow::{Context, Result, bail};
use nalgebra::DMatrix;
use rnaa_core::model::AdjustedMatrix;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::MatrixAdjuster;
use rnaa_formats::matrix::{read_matrix_tsv, write_matrix_tsv};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

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
        self.adjust_typed(
            matrix_path,
            samplesheet_path,
            model,
            paths,
            &rnaa_core::metadata::MetadataConfig::default(),
        )
    }
}

impl Residualizer {
    pub fn adjust_typed(
        &self,
        matrix_path: &Path,
        samplesheet_path: &Path,
        model: &str,
        paths: &ProjectPaths,
        schema: &rnaa_core::metadata::MetadataConfig,
    ) -> Result<AdjustedMatrix> {
        let matrix = read_matrix_tsv(matrix_path)?;
        if matrix.col_ids.is_empty() {
            bail!("input matrix has no samples: {}", matrix_path.display());
        }
        let sample_metadata = load_samplesheet(samplesheet_path)?;
        let prepared = build_design(&matrix.col_ids, &sample_metadata, model, schema)?;
        let design = DMatrix::from_row_slice(
            matrix.col_ids.len(),
            prepared.contract.coefficient_names.len(),
            &prepared.contract.matrix,
        );
        rnaa_core::util::write_json_pretty(
            &paths.corr_dir.join("model_design.json"),
            &prepared.contract,
        )?;
        let adjusted_values = residualize_matrix(&matrix.values, &design)?;
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
        .position(|name| name == "sample_id")
        .or_else(|| headers.iter().position(|name| name == "run_accession"))
        .ok_or_else(|| anyhow::anyhow!("samplesheet missing sample_id or run_accession column"))?;
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
        if run.trim().is_empty() || records.contains_key(run) {
            bail!("duplicate or empty samplesheet run identifier '{run}'");
        }
        records.insert(run.to_string(), map);
    }
    Ok(records)
}

pub fn prepare_residualization(
    matrix: &Path,
    sheet: &Path,
    formula: &str,
    schema: &rnaa_core::metadata::MetadataConfig,
) -> Result<crate::design::PreparedDesign> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(matrix)?;
    let ids = reader
        .headers()?
        .iter()
        .skip(1)
        .map(str::to_string)
        .collect::<Vec<_>>();
    build_design(&ids, &load_samplesheet(sheet)?, formula, schema)
}

/// Prepare the same design while keeping matrix-header reads within an admitted size.
pub fn prepare_residualization_exact(
    matrix: &Path,
    admitted_bytes: u64,
    sheet: &Path,
    formula: &str,
    schema: &rnaa_core::metadata::MetadataConfig,
) -> Result<crate::design::PreparedDesign> {
    let ids = rnaa_formats::matrix::read_matrix_columns_exact(matrix, admitted_bytes)?;
    build_design(&ids, &load_samplesheet(sheet)?, formula, schema)
}

pub(super) fn prepare_samples(
    ids: &[String],
    sheet: &Path,
    formula: &str,
    schema: &rnaa_core::metadata::MetadataConfig,
) -> Result<crate::design::PreparedDesign> {
    build_design(ids, &load_samplesheet(sheet)?, formula, schema)
}

fn build_design(
    sample_ids: &[String],
    sample_metadata: &BTreeMap<String, BTreeMap<String, String>>,
    formula: &str,
    schema: &rnaa_core::metadata::MetadataConfig,
) -> Result<crate::design::PreparedDesign> {
    let headers = sample_metadata
        .values()
        .flat_map(|row| row.keys().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let rows = sample_ids
        .iter()
        .map(|sample| {
            let metadata = sample_metadata
                .get(sample)
                .ok_or_else(|| anyhow::anyhow!("sample '{sample}' missing from metadata"))?;
            Ok(headers
                .iter()
                .map(|header| metadata.get(header).cloned().unwrap_or_default())
                .collect::<Vec<_>>())
        })
        .collect::<Result<Vec<_>>>()?;
    crate::design::PreparedDesign::new(
        &headers,
        &rows,
        sample_ids,
        formula,
        schema,
        &[],
        crate::design::Purpose::Correlation,
    )
}

fn residualize_matrix(rows: &[Vec<f64>], design: &DMatrix<f64>) -> Result<Vec<Vec<f64>>> {
    if rows
        .iter()
        .any(|row| row.len() != design.nrows() || row.iter().any(|value| !value.is_finite()))
    {
        bail!("residualization requires finite values and matching sample dimensions");
    }
    let basis = residual_basis(design)?;
    Ok(rows
        .iter()
        .map(|row| residualize_row(row, basis.as_ref()))
        .collect())
}

pub(super) fn residual_basis(design: &DMatrix<f64>) -> Result<Option<(DMatrix<f64>, usize)>> {
    if design.ncols() == 0 {
        return Ok(None);
    }
    // Decompose once for the whole cohort; never silently replace the requested fit.
    let svd = design.clone().svd(true, true);
    let tolerance = svd.singular_values.iter().copied().fold(0.0_f64, f64::max)
        * (design.nrows().max(design.ncols()) as f64)
        * f64::EPSILON
        * 100.0;
    let rank = svd
        .singular_values
        .iter()
        .filter(|value| **value > tolerance)
        .count();
    if rank < design.ncols() {
        bail!("adjustment design is rank deficient; remove confounded or constant covariates");
    }
    if design.nrows().saturating_sub(rank) < 2 {
        bail!("adjustment design leaves fewer than two residual degrees of freedom");
    }
    let u = svd
        .u
        .ok_or_else(|| anyhow::anyhow!("SVD did not produce a basis"))?;
    Ok(Some((u, rank)))
}

pub(super) fn residualize_row(row: &[f64], basis: Option<&(DMatrix<f64>, usize)>) -> Vec<f64> {
    let Some((u, rank)) = basis else {
        return row.to_vec();
    };
    let basis = u.columns(0, *rank);
    let y = nalgebra::DVector::from_column_slice(row);
    let residual = &y - basis * (basis.transpose() * &y);
    residual.iter().copied().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formula_parser_respects_interactions_and_intercept_removal() {
        let ids = (0..8).map(|i| format!("S{i}")).collect::<Vec<_>>();
        let metadata = ids
            .iter()
            .enumerate()
            .map(|(i, id)| {
                (
                    id.clone(),
                    BTreeMap::from([
                        (
                            "batch".into(),
                            if i % 4 < 2 { "A".into() } else { "B".into() },
                        ),
                        (
                            "condition".into(),
                            if i % 2 == 0 { "C".into() } else { "D".into() },
                        ),
                        ("x".into(), i.to_string()),
                    ]),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let mut schema = rnaa_core::metadata::MetadataConfig::default();
        schema.columns.insert(
            "x".into(),
            rnaa_core::metadata::ColumnSpec::new(rnaa_core::metadata::ColumnType::Numeric),
        );
        let interaction = build_design(&ids, &metadata, "~ batch * condition", &schema).unwrap();
        assert_eq!(interaction.contract.coefficient_names.len(), 4);
        let no_intercept = build_design(&ids, &metadata, "~ 0 + x", &schema).unwrap();
        assert_eq!(no_intercept.contract.coefficient_names.len(), 1);
        for i in 0..8 {
            assert_eq!(no_intercept.contract.matrix[i], i as f64);
        }
        assert!(build_design(&ids, &metadata, "~ missing", &schema).is_err());
        assert!(build_design(&ids, &metadata, "~ offset(x)", &schema).is_err());
    }

    #[test]
    fn rejects_saturated_design_and_removes_requested_effects() {
        let design = DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0]);
        let rows = vec![vec![1.0, 3.0, 10.0, 14.0]];
        let adjusted = residualize_matrix(&rows, &design).unwrap();
        for (actual, expected) in adjusted[0].iter().zip([-1.0, 1.0, -2.0, 2.0]) {
            assert!((actual - expected).abs() < 1e-10);
        }
        assert!(residualize_matrix(&[vec![1.0, 2.0]], &DMatrix::identity(2, 2)).is_err());
    }

    #[test]
    fn rowwise_projection_matches_the_original_view_calculation_exactly() -> Result<()> {
        let design = DMatrix::from_row_slice(
            6,
            3,
            &[
                1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 0.0, 1.0,
                2.0, 1.0,
            ],
        );
        let svd = design.clone().svd(true, true);
        let u = svd.u.unwrap();
        let original_basis = u.columns(0, 3);
        let new_basis = residual_basis(&design)?;
        for i in 0..80 {
            let row = (0..6)
                .map(|j| ((i * 37 + j * j * 13) % 107) as f64 / 7.0 - 4.0)
                .collect::<Vec<_>>();
            let y = nalgebra::DVector::from_column_slice(&row);
            let original = &y - original_basis * (original_basis.transpose() * &y);
            let actual = residualize_row(&row, new_basis.as_ref());
            assert_eq!(
                actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                original.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
        }
        Ok(())
    }
}
