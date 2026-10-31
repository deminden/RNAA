//! Matrix adjustment, correlation output and module enrichment have separate contracts.
mod adjust_stream;
mod limits;
mod modules;
mod network;
mod prepared;
mod residualize;
mod stream;

pub use adjust_stream::adjust_matrix_streaming;
pub use limits::{AdjustmentEstimate, CorrelationEstimate, CorrelationLimits, estimate_adjustment};
pub use modules::{
    MODULE_ALGORITHM_VERSION, ModuleGseaParams, ModuleScore, run_module_gsea,
    run_module_gsea_with_limits, standardized_module_score, validate_module_gsea,
};
pub use network::{BudgetedCorrelator, CORRELATION_ALGORITHM_VERSION, MinCorrCorrelator};
pub use residualize::{Residualizer, prepare_residualization, prepare_residualization_exact};
