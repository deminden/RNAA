//! Correlation computation and edge/dense output serialization.
use anyhow::{Context, Result, bail};
use rayon::ThreadPoolBuilder;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{AdjustedMatrix, CorrelationMethod, OutputMode};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::Correlator;
use rnaa_core::util::{file_size, now_rfc3339, write_json_pretty};
use rnaa_formats::matrix::read_matrix_tsv;
use serde_json::json;
use std::collections::BTreeMap;
use std::fs;

#[derive(Debug, Clone, Default)]
pub struct MinCorrCorrelator;
pub const CORRELATION_ALGORITHM_VERSION: u32 = 7;
use super::{
    limits::{CorrelationEstimate, CorrelationLimits},
    prepared::Prepared,
    stream,
};

pub struct BudgetedCorrelator {
    limits: CorrelationLimits,
    admission: Option<CorrelationEstimate>,
}
impl MinCorrCorrelator {
    pub fn with_limits(self, limits: CorrelationLimits) -> BudgetedCorrelator {
        BudgetedCorrelator {
            limits,
            admission: None,
        }
    }
}

impl BudgetedCorrelator {
    /// Require the execution preflight to match the shared-pool reservation.
    pub fn with_admission(mut self, estimate: CorrelationEstimate) -> Self {
        self.admission = Some(estimate);
        self
    }
}

impl BudgetedCorrelator {
    /// Check the caller's ownership fence between input preparation and tiles.
    #[allow(clippy::too_many_arguments)]
    pub fn correlate_checked(
        &self,
        matrix: &AdjustedMatrix,
        method: CorrelationMethod,
        output_mode: &OutputMode,
        threads: usize,
        paths: &ProjectPaths,
        project_id: &str,
        check: &impl Fn() -> Result<()>,
    ) -> Result<Vec<std::path::PathBuf>> {
        check()?;
        let started_at = now_rfc3339();
        let corr_threads = threads.max(1);
        let estimate =
            self.limits
                .preflight(&matrix.matrix_path, method, output_mode, corr_threads)?;
        if self.admission.is_some_and(|admitted| admitted != estimate) {
            bail!("correlation input or resource estimate changed after shared admission");
        }
        check()?;
        let matrix_data = read_matrix_tsv(&matrix.matrix_path)?;
        if matrix_data.row_ids.len() != estimate.genes
            || matrix_data.col_ids.len() != estimate.samples
        {
            bail!("correlation input changed after resource preflight");
        }
        check()?;
        let pool = ThreadPoolBuilder::new()
            .num_threads(corr_threads)
            .build()
            .context("failed to build correlation thread pool")?;
        let prepared = Prepared::new(matrix_data.values, method);
        check()?;
        let corr_dir = paths.corr_project_dir(project_id);
        fs::create_dir_all(&corr_dir)
            .with_context(|| format!("failed to create {}", corr_dir.display()))?;

        let filename = match output_mode {
            OutputMode::TopK { .. } => "edges_topk.tsv",
            OutputMode::Threshold { .. } => "edges_threshold.tsv",
            OutputMode::Dense => "corr_dense.tsv",
        };
        let out_path = corr_dir.join(filename);
        let output_budget = u64::try_from(estimate.output_bytes)
            .context("correlation output estimate exceeds byte accounting")?;
        let statistics = stream::write(
            &out_path,
            &matrix_data.row_ids,
            &prepared,
            output_mode,
            self.limits.block_rows,
            output_budget,
            &pool,
            check,
        )?;
        let mut outputs = vec![out_path];

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
                "output_mode": output_mode,
                "threads": corr_threads,
                "algorithm_version":CORRELATION_ALGORITHM_VERSION,
                "constant_genes":matrix_data.row_ids.iter().zip(&prepared.constant).filter(|(_,constant)|**constant).map(|(gene,_)|gene).collect::<Vec<_>>(),
                "constant_row_policy":"undefined_correlations_v1",
                "endpoint_policy":"snap_within_8_machine_epsilons_v1",
                "computation":"bounded_tiles_and_row_stripes_v1",
                "limits":self.limits,
                "resource_estimate":estimate,
                "statistics":statistics,
                "topk_tie_policy":"absolute_coefficient_then_input_target_order",
                "threshold_policy":"inclusive_absolute_coefficient_upper_triangle"
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

impl Correlator for BudgetedCorrelator {
    fn correlate(
        &self,
        matrix: &AdjustedMatrix,
        method: CorrelationMethod,
        output_mode: &OutputMode,
        threads: usize,
        paths: &ProjectPaths,
        project_id: &str,
    ) -> Result<Vec<std::path::PathBuf>> {
        self.correlate_checked(
            matrix,
            method,
            output_mode,
            threads,
            paths,
            project_id,
            &|| Ok(()),
        )
    }
}

impl Correlator for MinCorrCorrelator {
    fn correlate(
        &self,
        matrix: &AdjustedMatrix,
        method: CorrelationMethod,
        output_mode: &OutputMode,
        threads: usize,
        paths: &ProjectPaths,
        project_id: &str,
    ) -> Result<Vec<std::path::PathBuf>> {
        self.clone()
            .with_limits(CorrelationLimits::default())
            .correlate(matrix, method, output_mode, threads, paths, project_id)
    }
}

#[cfg(test)]
mod admission_tests {
    use super::*;
    use anyhow::Result;

    #[test]
    fn refuses_changed_input_before_using_an_old_reservation() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let matrix_path = temp.path().join("matrix.tsv");
        fs::write(
            &matrix_path,
            "gene_id\tA\tB\ng1\t1\t2\ng2\t3\t4\ng3\t5\t6\n",
        )?;
        let limits = CorrelationLimits::default();
        let mode = OutputMode::TopK { k: 1 };
        let estimate = limits.preflight(&matrix_path, CorrelationMethod::Pearson, &mode, 1)?;
        fs::write(
            &matrix_path,
            "gene_id\tA\tB\nlonger_gene\t1\t2\ng2\t3\t4\ng3\t5\t6\n",
        )?;
        let matrix = AdjustedMatrix {
            matrix_path,
            genes: vec![],
            samples: vec![],
        };
        let error = MinCorrCorrelator
            .with_limits(limits)
            .with_admission(estimate)
            .correlate(
                &matrix,
                CorrelationMethod::Pearson,
                &mode,
                1,
                &ProjectPaths::new(temp.path()),
                "test",
            )
            .unwrap_err();
        assert!(error.to_string().contains("changed after shared admission"));
        Ok(())
    }
}
