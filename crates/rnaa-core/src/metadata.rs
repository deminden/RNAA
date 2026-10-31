//! Explicit statistical metadata roles and frozen factor/numeric values.
mod formula;
use anyhow::{Context, Result, bail};
pub use formula::variables;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ColumnType {
    Factor,
    Numeric,
    Identifier,
    Ignore,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ColumnSpec {
    #[serde(rename = "type")]
    pub kind: ColumnType,
    #[serde(default)]
    pub levels: Vec<String>,
    #[serde(default)]
    pub reference: Option<String>,
}
impl ColumnSpec {
    pub fn new(kind: ColumnType) -> Self {
        Self {
            kind,
            levels: Vec::new(),
            reference: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MetadataConfig {
    #[serde(default)]
    pub columns: BTreeMap<String, ColumnSpec>,
    #[serde(default)]
    pub missing_values: Vec<String>,
    #[serde(default = "minimum_replicates")]
    pub minimum_replicates: usize,
}
fn minimum_replicates() -> usize {
    2
}
impl Default for MetadataConfig {
    fn default() -> Self {
        Self {
            columns: BTreeMap::new(),
            missing_values: Vec::new(),
            minimum_replicates: minimum_replicates(),
        }
    }
}
impl MetadataConfig {
    pub fn column(&self, name: &str) -> Option<ColumnSpec> {
        self.columns.get(name).cloned().or_else(|| {
            if matches!(name, "condition" | "batch") {
                Some(ColumnSpec::new(ColumnType::Factor))
            } else if matches!(
                name,
                "project_id" | "study_accession" | "experiment_accession"
            ) || protected(name)
            {
                Some(ColumnSpec::new(ColumnType::Identifier))
            } else {
                None
            }
        })
    }
    pub fn validate(&self) -> Result<()> {
        if self.minimum_replicates == 0 {
            bail!("metadata.minimum_replicates must be positive");
        }
        if self
            .missing_values
            .iter()
            .any(|value| value != value.trim())
            || self.missing_values.iter().collect::<BTreeSet<_>>().len()
                != self.missing_values.len()
        {
            bail!("metadata missing_values must be unique and have no surrounding whitespace");
        }
        for (name, spec) in &self.columns {
            if name.is_empty()
                || name != name.trim()
                || (matches!(spec.kind, ColumnType::Factor | ColumnType::Numeric)
                    && !name.chars().enumerate().all(|(i, c)| {
                        c.is_ascii_alphabetic() || c == '_' || (i > 0 && c.is_ascii_digit())
                    }))
            {
                bail!("invalid metadata column identifier '{name}'");
            }
            if protected(name) && !matches!(spec.kind, ColumnType::Identifier | ColumnType::Ignore)
            {
                bail!("sample/run identity or path '{name}' cannot be a model covariate");
            }
            if spec.kind != ColumnType::Factor
                && (!spec.levels.is_empty() || spec.reference.is_some())
            {
                bail!("only factor metadata '{name}' may declare levels/reference");
            }
            if spec
                .levels
                .iter()
                .any(|value| value.is_empty() || value != value.trim())
                || spec.levels.iter().collect::<BTreeSet<_>>().len() != spec.levels.len()
            {
                bail!("factor '{name}' levels must be unique, nonempty and trimmed");
            }
            if let Some(reference) = &spec.reference
                && (reference.is_empty()
                    || reference != reference.trim()
                    || (!spec.levels.is_empty() && !spec.levels.contains(reference)))
            {
                bail!("invalid reference level '{reference}' for factor '{name}'");
            }
        }
        Ok(())
    }
    pub fn normalize(&self, name: &str, raw: &str) -> Result<String> {
        let spec = self.column(name).with_context(|| format!("undeclared covariate '{name}'; declare metadata.columns.{name}.type as factor or numeric"))?;
        if !matches!(spec.kind, ColumnType::Factor | ColumnType::Numeric) {
            bail!(
                "metadata column '{name}' is {:?}, so it cannot be a model covariate",
                spec.kind
            );
        }
        let value = raw.trim();
        if value.is_empty() || self.missing_values.iter().any(|missing| missing == value) {
            bail!("design variable '{name}' contains missing values");
        }
        if spec.kind == ColumnType::Numeric {
            let number = value.parse::<f64>().with_context(|| {
                format!("numeric metadata '{name}' has invalid value '{value}'")
            })?;
            if !number.is_finite() {
                bail!("numeric metadata '{name}' must contain finite values");
            }
            return Ok(if number == 0.0 {
                "0".into()
            } else {
                number.to_string()
            });
        }
        if !spec.levels.is_empty() && !spec.levels.iter().any(|level| level == value) {
            bail!("undeclared level '{value}' for factor '{name}'");
        }
        Ok(value.into())
    }
}

fn protected(name: &str) -> bool {
    matches!(
        name,
        "run_accession"
            | "sample_accession"
            | "sample_id"
            | "biological_sample_id"
            | "identity_source"
            | "run_accessions"
            | "read1_path"
            | "read2_path"
    )
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum FrozenColumn {
    Factor {
        values: Vec<String>,
        levels: Vec<String>,
        reference: String,
    },
    Numeric {
        values: Vec<f64>,
    },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FrozenMetadata {
    pub schema: MetadataConfig,
    pub columns: BTreeMap<String, FrozenColumn>,
}
impl FrozenMetadata {
    pub fn resolve(
        headers: &[String],
        rows: &[Vec<String>],
        formula: &str,
        schema: &MetadataConfig,
    ) -> Result<Self> {
        schema.validate()?;
        if rows.is_empty()
            || rows.iter().any(|row| row.len() != headers.len())
            || headers.iter().collect::<BTreeSet<_>>().len() != headers.len()
        {
            bail!("model metadata must have unique columns and rectangular nonempty rows");
        }
        let mut columns = BTreeMap::new();
        for name in variables(formula, headers, schema)? {
            let index = headers
                .iter()
                .position(|header| header == &name)
                .context("design column disappeared")?;
            let spec = schema.column(&name).context("undeclared model column")?;
            let values = rows
                .iter()
                .map(|row| schema.normalize(&name, &row[index]))
                .collect::<Result<Vec<_>>>()?;
            let column = if spec.kind == ColumnType::Numeric {
                FrozenColumn::Numeric {
                    values: values
                        .iter()
                        .map(|value| value.parse::<f64>())
                        .collect::<std::result::Result<_, _>>()?,
                }
            } else {
                let observed = values.iter().cloned().collect::<BTreeSet<_>>();
                let declared = if spec.levels.is_empty() {
                    observed.iter().cloned().collect::<Vec<_>>()
                } else {
                    spec.levels.clone()
                };
                let reference = spec.reference.unwrap_or_else(|| declared[0].clone());
                if !observed.contains(&reference) {
                    bail!(
                        "reference level '{reference}' for factor '{name}' is absent from the biological sample cohort"
                    );
                }
                let levels = declared
                    .into_iter()
                    .filter(|level| observed.contains(level))
                    .collect();
                FrozenColumn::Factor {
                    values,
                    levels,
                    reference,
                }
            };
            columns.insert(name, column);
        }
        Ok(Self {
            schema: schema.clone(),
            columns,
        })
    }
}

#[cfg(test)]
mod tests;
