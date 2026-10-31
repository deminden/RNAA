//! Values that determine the fitted model, separate from source/schema witnesses.
use super::{DesignContract, Purpose, ResolvedContrast};
use rnaa_core::metadata::FrozenColumn;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScientificDesign {
    pub schema_version: u32,
    pub formula: String,
    pub purpose: Purpose,
    pub samples: Vec<String>,
    pub columns: BTreeMap<String, FrozenColumn>,
    pub coefficient_names: Vec<String>,
    pub matrix: Vec<f64>,
    pub offsets: Vec<f64>,
    pub factor_references: BTreeMap<String, String>,
    pub rank: usize,
    pub residual_degrees_of_freedom: usize,
    pub contrasts: Vec<ResolvedContrast>,
}

impl DesignContract {
    pub fn scientific(&self) -> ScientificDesign {
        ScientificDesign {
            schema_version: 1,
            formula: self.formula.clone(),
            purpose: self.purpose,
            samples: self.samples.clone(),
            columns: self.metadata.columns.clone(),
            coefficient_names: self.coefficient_names.clone(),
            matrix: self.matrix.clone(),
            offsets: self.offsets.clone(),
            factor_references: self.factor_references.clone(),
            rank: self.rank,
            residual_degrees_of_freedom: self.residual_degrees_of_freedom,
            contrasts: self.contrasts.clone(),
        }
    }
}
