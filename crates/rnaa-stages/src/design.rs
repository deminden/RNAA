//! Shared typed model construction and statistical preflight for DE/correlation.
mod contrasts;
mod scientific;
use anyhow::{Context, Result, bail};
use rnaa_core::{
    ContrastSpec,
    metadata::{FrozenColumn, FrozenMetadata, MetadataConfig},
};
use rsdeseq2::prelude::{FormulaFactorColumn, FormulaModelFrame, FormulaNumericColumn};
pub use scientific::ScientificDesign;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Purpose {
    DifferentialExpression,
    Correlation,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedContrast {
    pub name: String,
    pub factor: String,
    pub source_factor: String,
    pub numerator: String,
    pub denominator: String,
    pub reference: Option<String>,
    pub sample_levels: Vec<String>,
    pub weights: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignContract {
    pub schema_version: u32,
    pub formula: String,
    pub purpose: Purpose,
    pub samples: Vec<String>,
    pub metadata: FrozenMetadata,
    pub coefficient_names: Vec<String>,
    pub matrix: Vec<f64>,
    pub offsets: Vec<f64>,
    pub factor_references: BTreeMap<String, String>,
    pub rank: usize,
    pub residual_degrees_of_freedom: usize,
    pub contrasts: Vec<ResolvedContrast>,
}
#[derive(Debug)]
pub struct PreparedDesign {
    pub contract: DesignContract,
    pub(crate) frame: FormulaModelFrame,
}

impl PreparedDesign {
    pub fn new(
        headers: &[String],
        rows: &[Vec<String>],
        samples: &[String],
        formula: &str,
        schema: &MetadataConfig,
        contrasts: &[ContrastSpec],
        purpose: Purpose,
    ) -> Result<Self> {
        if samples.len() != rows.len()
            || samples.iter().collect::<BTreeSet<_>>().len() != samples.len()
            || samples.iter().any(|id| id.trim().is_empty())
        {
            bail!("design sample identities must be unique and match metadata rows");
        }
        if purpose == Purpose::Correlation && rsdeseq2::design::formula_has_offset_terms(formula)? {
            bail!("offset terms are not supported for linear residualization");
        }
        let metadata = FrozenMetadata::resolve(headers, rows, formula, schema)?;
        let mut frame = FormulaModelFrame::default();
        for (name, column) in &metadata.columns {
            match column {
                FrozenColumn::Numeric { values } => {
                    frame.numeric_covariates.push(FormulaNumericColumn {
                        name: name.clone(),
                        values: values.clone(),
                    })
                }
                FrozenColumn::Factor {
                    values,
                    levels,
                    reference,
                } => frame.factors.push(FormulaFactorColumn {
                    name: name.clone(),
                    sample_levels: values.clone(),
                    levels: Some(levels.clone()),
                    reference: Some(reference.clone()),
                }),
            }
        }
        if frame.factors.is_empty() && frame.numeric_covariates.is_empty() {
            // The engine needs the sample dimension even for an intercept-only formula.
            frame.numeric_covariates.push(FormulaNumericColumn {
                name: "__rnaa_sample_count".into(),
                values: vec![0.0; samples.len()],
            });
        }
        let expanded = rsdeseq2::design::expanded_formula_design_with_offsets_from_model_frame(
            formula, &frame,
        )
        .context("invalid typed design formula")?;
        let design = &expanded.design.standard_design;
        let rank = design.rank()?;
        if rank != design.n_coefficients() {
            bail!(
                "design matrix is rank deficient: rank {rank} for {} coefficients; check confounded covariates",
                design.n_coefficients()
            );
        }
        if design.n_samples() <= rank {
            bail!(
                "design is saturated: {} biological samples for {rank} coefficients; residual degrees of freedom must be positive",
                design.n_samples()
            );
        }
        if purpose == Purpose::DifferentialExpression && rank == 0 {
            bail!(
                "DE requires at least one model coefficient; use ~ 1 for normalization without covariates"
            );
        }
        let resolved =
            contrasts::resolve(design, &expanded.model_frame, headers, schema, contrasts)?;
        let mut factor_references = expanded
            .model_frame
            .resolved_factor_references()?
            .into_iter()
            .map(|item| (item.factor.to_string(), item.reference.to_string()))
            .collect::<BTreeMap<_, _>>();
        for contrast in &resolved {
            if let Some(reference) = &contrast.reference {
                factor_references.insert(contrast.source_factor.clone(), reference.clone());
                factor_references.insert(contrast.factor.clone(), reference.clone());
            } else {
                factor_references.remove(&contrast.source_factor);
                factor_references.remove(&contrast.factor);
            }
        }
        let contract = DesignContract {
            schema_version: 1,
            formula: formula.into(),
            purpose,
            samples: samples.to_vec(),
            metadata,
            coefficient_names: design
                .coefficient_names()
                .context("typed model coefficient names are missing")?
                .to_vec(),
            matrix: design.matrix().as_slice().to_vec(),
            offsets: expanded.offsets,
            factor_references,
            rank,
            residual_degrees_of_freedom: samples.len() - rank,
            contrasts: resolved,
        };
        Ok(Self { contract, frame })
    }
}

#[cfg(test)]
mod tests;
