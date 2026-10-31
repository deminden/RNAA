//! Apply DESeq2 results() policies after constructing each explicit contrast.
use anyhow::{Context, Result};
use rnaa_core::metadata::FrozenColumn;
use rsdeseq2::prelude::{
    CountMatrix, DeseqFit, DeseqResults, IndependentFilteringOptions, apply_cooks_cutoff,
    apply_cooks_cutoff_with_low_count_heuristic, apply_independent_filtering, default_cooks_cutoff,
};
use serde::Serialize;

#[derive(Serialize)]
pub(super) struct ResultPolicy {
    cooks_cutoff: Option<f64>,
    low_count_heuristic: bool,
    cooks_outliers: usize,
    independent_filter_alpha: f64,
    filter_theta: Option<f64>,
    filter_threshold: Option<f64>,
}

pub(super) fn finalize(
    results: &mut DeseqResults,
    fit: &DeseqFit,
    counts: &CountMatrix,
    contract: &crate::design::DesignContract,
) -> Result<ResultPolicy> {
    let cutoff = default_cooks_cutoff(counts.n_samples(), contract.rank)?;
    // DESeq2's special rescue applies only to a formula with one two-level factor.
    let low_count_heuristic = contract.metadata.columns.len() == 1
        && contract.metadata.columns.values().all(
            |column| matches!(column, FrozenColumn::Factor { levels, .. } if levels.len() == 2),
        );
    if low_count_heuristic {
        apply_cooks_cutoff_with_low_count_heuristic(
            results,
            cutoff,
            counts,
            fit.cooks
                .as_ref()
                .context("fit is missing Cook's distances")?,
        )?;
    } else {
        apply_cooks_cutoff(results, cutoff)?;
    }
    let options = IndependentFilteringOptions::default();
    let filtering = apply_independent_filtering(results, &options)?;
    Ok(ResultPolicy {
        cooks_cutoff: cutoff,
        low_count_heuristic,
        cooks_outliers: results
            .rows
            .iter()
            .filter(|row| row.cooks_outlier == Some(true))
            .count(),
        independent_filter_alpha: options.alpha,
        filter_theta: filtering.filter_theta,
        filter_threshold: filtering.filter_threshold,
    })
}
