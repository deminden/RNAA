//! Strict module admission, bounded graph/PCA work and owned local thread budget.
use super::super::{limits::CorrelationLimits, prepared::Prepared};
use super::{
    MODULE_ALGORITHM_VERSION, ModuleGseaParams, admission, graph, output::Tables,
    pca::standardized_module_score, ranking::t_score,
};
use anyhow::{Context, Result, bail};
use rnaa_core::{
    ManifestArtifact, ProjectPaths, StageManifest,
    model::{AdjustedMatrix, CorrelationMethod},
    util::{now_rfc3339, write_json_pretty},
};
use rsfgsea::prelude::{RankedList, ScoreType, fgsea};
use serde_json::json;
use std::{collections::BTreeMap, fs, path::PathBuf};

pub fn validate_module_gsea(
    matrix: &AdjustedMatrix,
    params: &ModuleGseaParams,
    limits: CorrelationLimits,
    threads: usize,
) -> Result<()> {
    admission::load(&matrix.matrix_path, params, limits, threads).map(|_| ())
}

pub fn run_module_gsea(
    matrix: &AdjustedMatrix,
    paths: &ProjectPaths,
    project_id: &str,
    params: &ModuleGseaParams,
) -> Result<Vec<PathBuf>> {
    run_module_gsea_with_limits(
        matrix,
        paths,
        project_id,
        params,
        CorrelationLimits::default(),
        1,
    )
}

pub fn run_module_gsea_with_limits(
    matrix: &AdjustedMatrix,
    paths: &ProjectPaths,
    project_id: &str,
    params: &ModuleGseaParams,
    limits: CorrelationLimits,
    threads: usize,
) -> Result<Vec<PathBuf>> {
    let started_at = now_rfc3339();
    params.validate()?;
    let inputs = vec![
        ManifestArtifact::snapshot("ADJUSTED_MATRIX", &matrix.matrix_path)?,
        ManifestArtifact::snapshot("PATHWAYS", &params.gmt_path)?,
    ];
    let (data, pathways, admission) =
        admission::load(&matrix.matrix_path, params, limits, threads)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("failed to build module thread pool")?;
    let execute = Execution {
        data,
        pathways,
        admission,
        inputs,
        started_at,
        paths,
        project_id,
        params,
        limits,
        threads,
    };
    pool.install(|| execute.run())
}

struct Execution<'a> {
    data: rnaa_formats::matrix::MatrixData,
    pathways: Vec<rsfgsea::prelude::Pathway>,
    admission: admission::Admission,
    inputs: Vec<ManifestArtifact>,
    started_at: String,
    paths: &'a ProjectPaths,
    project_id: &'a str,
    params: &'a ModuleGseaParams,
    limits: CorrelationLimits,
    threads: usize,
}
impl Execution<'_> {
    fn run(self) -> Result<Vec<PathBuf>> {
        let Self {
            data,
            pathways,
            admission,
            inputs,
            started_at,
            paths,
            project_id,
            params,
            limits,
            threads,
        } = self;

        let prepared = Prepared::new(data.values.clone(), CorrelationMethod::Spearman);
        let (modules, graph_statistics) =
            graph::detect(&prepared, params.min_module_size, limits.block_rows);
        if modules.is_empty() {
            bail!(
                "no modules detected with min size {}; reduce module-min-size",
                params.min_module_size
            );
        }
        let directory = paths.corr_project_dir(project_id);
        fs::create_dir_all(&directory)?;
        let mut tables = Tables::open(&directory, &data.col_ids)?;
        let mut diagnostics = Vec::new();
        for (index, module) in modules.iter().enumerate() {
            let id = format!("M{:04}", index + 1);
            let rows = module
                .members
                .iter()
                .map(|&gene| data.values[gene].clone())
                .collect::<Vec<_>>();
            let score = standardized_module_score(&rows)
                .with_context(|| format!("PCA failed for module {id}"))?;
            tables.module(&id, module, &data, &score, params.top_genes)?;
            let correlations = prepared.correlations_with(&score.values);
            let mut genes = Vec::new();
            let mut scores = Vec::new();
            for (gene, &r) in data.row_ids.iter().zip(&correlations) {
                if !r.is_finite() {
                    continue;
                }
                let t = t_score(r, data.col_ids.len());
                tables.rank(&id, gene, r, t)?;
                genes.push(gene.clone());
                scores.push(t);
            }
            if genes.len() != admission.eligible_genes {
                bail!("module score produced an undefined correlation for a nonconstant gene");
            }
            let results = fgsea(
                &RankedList::new(genes, scores),
                &pathways,
                Some(params.permutations),
                0,
                params.seed,
                admission.pathway_min_size,
                admission.pathway_max_size,
                1e-50,
                ScoreType::Std,
                1.0,
            );
            if results.len() != pathways.len() {
                bail!("enrichment omitted an admitted pathway for module {id}");
            }
            tables.enrichment(&id, module.members.len(), &results)?;
            diagnostics.push(json!({"module_id":id,"module_size":module.members.len(),"reported_genes":params.top_genes.min(module.members.len()),"variance_explained":score.variance_explained,"leading_multiplicity":score.leading_multiplicity,"orientation":score.orientation}));
        }
        for input in &inputs {
            if !input.matches_file()? {
                bail!("module input changed during computation");
            }
        }
        let mut outputs = tables.finish()?;
        let manifest_path = directory.join("module_fgsea_manifest.json");
        let manifest = StageManifest {
            schema_version:1,stage:"module_gsea".into(),id:project_id.into(),started_at,finished_at:now_rfc3339(),
            parameters:json!({"algorithm_version":MODULE_ALGORITHM_VERSION,"method":"spearman","neighbor_k":20,
                "score_policy":"standardized_full_component_pca_v1","score_variance":"unit_sample_variance",
                "tied_eigenspace_policy":"first_nonzero_canonical_sample_projection_v1",
                "rank_policy":"omit_undefined_and_cap_abs_r_at_1_minus_1e_7_v1",
                "connectivity_rank_policy":"round_to_1e_minus_10_then_input_gene_order_v1",
                "top_genes_per_module":params.top_genes,"min_module_size":params.min_module_size,
                "permutations":params.permutations,"seed":params.seed,"threads":threads.max(1),"limits":limits,
                "admission":admission,"graph_statistics":graph_statistics,"modules":diagnostics}),
            tool_versions:BTreeMap::from([("mincorr".into(),crate::engines::MINCORR_REVISION.into()),("rsfgsea".into(),crate::engines::RSFGSEA.into())]),
            input_artifacts:inputs, output_artifacts:outputs.iter().map(|path| ManifestArtifact::snapshot("MODULE_GSEA_OUTPUT",path)).collect::<Result<_>>()?,
            notes:vec!["Full modules are connected components of mutual positive Spearman top-20 neighbours. Reporting limits do not affect membership or PCA scores.".into()],
        };
        write_json_pretty(&manifest_path, &manifest)?;
        outputs.push(manifest_path);
        Ok(outputs)
    }
}
