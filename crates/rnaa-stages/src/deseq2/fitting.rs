//! DESeq2's documented parametric-trend fallback, without masking other failures.
use anyhow::{Context, Result};
use rsdeseq2::prelude::{CountMatrix, DeseqBuilder, DeseqError, DeseqFit, FitType};
use serde::Serialize;

#[derive(Serialize)]
pub(super) struct DispersionFit {
    requested: &'static str,
    selected: &'static str,
    fallback_reason: Option<String>,
}
fn parametric_failure(error: &DeseqError) -> bool {
    match error {
        DeseqError::InvalidDispersion { reason } => {
            reason.contains("parametric dispersion") || reason == "dispersion fit did not converge"
        }
        DeseqError::NonFiniteValue { context, .. } => context.starts_with("parametric dispersion"),
        _ => false,
    }
}
pub(super) fn fit_formula(
    builder: &DeseqBuilder,
    counts: &CountMatrix,
    formula: &str,
) -> Result<(DeseqFit, DispersionFit)> {
    match builder.fit_formula_with_results(counts, formula) {
        Ok((fit, _)) => Ok((
            fit,
            DispersionFit {
                requested: "parametric",
                selected: "parametric",
                fallback_reason: None,
            },
        )),
        Err(error) if parametric_failure(&error) => {
            let reason = error.to_string();
            let (fit,_)=builder.clone().fit_type(FitType::Local).fit_formula_with_results(counts,formula)
                .with_context(||format!("parametric trend failed ({reason}); DESeq2 local-trend fallback also failed"))?;
            Ok((
                fit,
                DispersionFit {
                    requested: "parametric",
                    selected: "local",
                    fallback_reason: Some(reason),
                },
            ))
        }
        Err(error) => Err(error.into()),
    }
}
