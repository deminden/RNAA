//! Module workflows keep graph membership, PCA scores and reporting independent.
mod admission;
mod graph;
mod output;
mod pca;
mod ranking;
mod workflow;
pub use pca::{ModuleScore, standardized_module_score};
use std::path::PathBuf;
pub use workflow::{run_module_gsea, run_module_gsea_with_limits, validate_module_gsea};

pub const MODULE_ALGORITHM_VERSION: u32 = 2;
#[derive(Debug, Clone)]
pub struct ModuleGseaParams {
    pub gmt_path: PathBuf,
    pub top_genes: usize,
    pub min_module_size: usize,
    pub permutations: usize,
    pub seed: u64,
    pub pathway_min_size: usize,
    pub pathway_max_size: Option<usize>,
}

impl From<&rnaa_core::config::CorrModuleConfig> for ModuleGseaParams {
    fn from(config: &rnaa_core::config::CorrModuleConfig) -> Self {
        Self {
            gmt_path: config.gmt_path.clone(),
            top_genes: config.top_genes,
            min_module_size: config.min_module_size,
            permutations: config.permutations,
            seed: config.seed,
            pathway_min_size: config.pathway_min_size,
            pathway_max_size: config.pathway_max_size,
        }
    }
}

impl From<&ModuleGseaParams> for rnaa_core::config::CorrModuleConfig {
    fn from(params: &ModuleGseaParams) -> Self {
        Self {
            gmt_path: params.gmt_path.clone(),
            top_genes: params.top_genes,
            min_module_size: params.min_module_size,
            permutations: params.permutations,
            seed: params.seed,
            pathway_min_size: params.pathway_min_size,
            pathway_max_size: params.pathway_max_size,
        }
    }
}
