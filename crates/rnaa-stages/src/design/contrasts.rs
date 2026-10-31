//! Resolve categorical views against actual model coefficients, including
//! releveling, relabeling and full-indicator designs without an intercept.
use super::ResolvedContrast;
use anyhow::{Context, Result, bail};
use rnaa_core::{
    ContrastSpec,
    metadata::{MetadataConfig, variables},
};
use rsdeseq2::{
    design::{DesignMatrix, FormulaModelFrame},
    prelude::ResultsContrast,
};
use std::collections::BTreeSet;

pub(super) fn resolve(
    design: &DesignMatrix,
    frame: &FormulaModelFrame,
    headers: &[String],
    schema: &MetadataConfig,
    definitions: &[ContrastSpec],
) -> Result<Vec<ResolvedContrast>> {
    let mut resolved = Vec::new();
    let mut names = BTreeSet::new();
    for definition in definitions {
        if definition.name.trim().is_empty()
            || definition
                .name
                .chars()
                .any(|c| c.is_control() || c == '/' || c == '\\')
            || !names.insert(&definition.name)
        {
            bail!("contrast names must be unique, nonempty safe filename components");
        }
        let origin = BTreeSet::from([definition.factor.clone()]);
        let mut candidates = Vec::new();
        let mut reason = "factor is absent or is not categorical in this model".to_string();
        for column in &frame.factors {
            if column.name != definition.factor
                && !variables(&format!("~ {}", column.name), headers, schema)
                    .is_ok_and(|names| names == origin)
            {
                continue;
            }
            let request = ResultsContrast::Character {
                factor: column.name.clone(),
                numerator: definition.level_a.clone(),
                denominator: definition.level_b.clone(),
                reference: None,
            };
            let levels = match rsdeseq2::contrasts::factor_level_contrast_from_model_frame(
                &request, frame,
            ) {
                Ok(Some(levels)) => levels,
                Ok(None) => continue,
                Err(error) => {
                    reason = error.to_string();
                    continue;
                }
            };
            let actual = ResultsContrast::Character {
                factor: levels.factor.into(),
                numerator: levels.numerator.into(),
                denominator: levels.denominator.into(),
                reference: None,
            };
            let numeric = match rsdeseq2::contrasts::resolve_results_contrast(design, &actual) {
                Ok(value) => value.numeric,
                Err(error) => {
                    reason = error.to_string();
                    continue;
                }
            };
            if numeric.iter().all(|value| *value == 0.0) {
                continue;
            }
            candidates.push(ResolvedContrast {
                name: definition.name.clone(),
                factor: levels.factor.into(),
                source_factor: definition.factor.clone(),
                numerator: levels.numerator.into(),
                denominator: levels.denominator.into(),
                reference: unique_reference(design, levels.factor, levels.sample_levels),
                sample_levels: levels.sample_levels.to_vec(),
                weights: numeric,
            });
        }
        if candidates.len() > 1 {
            bail!(
                "contrast '{}' has multiple estimable views of '{}'; use one unambiguous factor coding",
                definition.name,
                definition.factor
            );
        }
        let actual = candidates.pop().with_context(|| {
            format!(
                "contrast '{}' ({} {} vs {}) is not estimable: {reason}",
                definition.name, definition.factor, definition.level_a, definition.level_b
            )
        })?;
        for level in [&actual.numerator, &actual.denominator] {
            let count = actual
                .sample_levels
                .iter()
                .filter(|value| *value == level)
                .count();
            if count < schema.minimum_replicates {
                bail!(
                    "contrast '{}' level '{level}' has {count} biological samples; minimum_replicates is {}",
                    definition.name,
                    schema.minimum_replicates
                );
            }
        }
        resolved.push(actual);
    }
    Ok(resolved)
}

/// Prove a treatment baseline against actual coefficient coding. A full
/// indicator model has no baseline: assuming one drops a contrast coefficient.
fn unique_reference(
    design: &rsdeseq2::design::DesignMatrix,
    factor: &str,
    values: &[String],
) -> Option<String> {
    let levels = values.iter().cloned().collect::<BTreeSet<_>>();
    let candidates = levels
        .iter()
        .filter(|reference| {
            levels
                .iter()
                .filter(|other| *other != *reference)
                .all(|other| {
                    let auto = ResultsContrast::Character {
                        factor: factor.into(),
                        numerator: other.clone(),
                        denominator: (*reference).clone(),
                        reference: None,
                    };
                    let explicit = ResultsContrast::Character {
                        factor: factor.into(),
                        numerator: other.clone(),
                        denominator: (*reference).clone(),
                        reference: Some((*reference).clone()),
                    };
                    match (
                        rsdeseq2::contrasts::resolve_results_contrast(design, &auto),
                        rsdeseq2::contrasts::resolve_results_contrast(design, &explicit),
                    ) {
                        (Ok(a), Ok(b)) => a.numeric == b.numeric,
                        _ => false,
                    }
                })
        })
        .collect::<Vec<_>>();
    if candidates.len() == 1 {
        Some(candidates[0].clone())
    } else {
        None
    }
}
