//! One statistical unit per biological sample, with complete run membership.
use crate::{ContrastSpec, cohorts::IdentitySource};
use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnalysisSample {
    pub sample_id: String,
    pub identity_source: IdentitySource,
    pub runs: Vec<String>,
    pub metadata: BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SamplePlan {
    pub schema_version: u32,
    pub samples: Vec<AnalysisSample>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleMembership {
    pub sample_id: String,
    pub identity_source: IdentitySource,
    pub runs: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleMembershipPlan {
    pub schema_version: u32,
    pub samples: Vec<SampleMembership>,
}

struct SourceRun {
    run: String,
    identity_source: IdentitySource,
    metadata: BTreeMap<String, String>,
}

impl SamplePlan {
    pub fn membership(&self) -> SampleMembershipPlan {
        SampleMembershipPlan {
            schema_version: 1,
            samples: self
                .samples
                .iter()
                .map(|sample| SampleMembership {
                    sample_id: sample.sample_id.clone(),
                    identity_source: sample.identity_source.clone(),
                    runs: sample.runs.clone(),
                })
                .collect(),
        }
    }
    pub fn from_sheet(path: &Path, design: &str, contrasts: &[ContrastSpec]) -> Result<Self> {
        Self::from_sheet_with_metadata(
            path,
            design,
            contrasts,
            &crate::metadata::MetadataConfig::default(),
        )
    }

    pub fn from_sheet_with_metadata(
        path: &Path,
        design: &str,
        contrasts: &[ContrastSpec],
        schema: &crate::metadata::MetadataConfig,
    ) -> Result<Self> {
        schema.validate()?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .with_context(|| format!("failed to read sample metadata {}", path.display()))?;
        let headers = reader.headers()?.clone();
        if headers.iter().any(|name| name.trim().is_empty())
            || headers.iter().collect::<BTreeSet<_>>().len() != headers.len()
        {
            bail!("sample metadata columns must be unique and nonempty");
        }
        let run_col = headers
            .iter()
            .position(|name| name == "run_accession")
            .context("sample metadata missing run_accession")?;
        let sample_col = headers.iter().position(|name| name == "sample_accession");
        let declared_col = headers
            .iter()
            .position(|name| name == "biological_sample_id");
        let source_col = headers.iter().position(|name| name == "identity_source");
        let header_names = headers.iter().map(str::to_string).collect::<Vec<_>>();
        let mut required = crate::metadata::variables(design, &header_names, schema)?;
        for contrast in contrasts {
            if !headers.iter().any(|name| name == contrast.factor) {
                bail!(
                    "contrast factor '{}' is missing from sample metadata",
                    contrast.factor
                );
            }
            required.insert(contrast.factor.clone());
        }
        let mut groups: BTreeMap<String, Vec<SourceRun>> = BTreeMap::new();
        let mut seen = BTreeSet::new();
        for record in reader.records() {
            let record = record?;
            let value = |column: Option<usize>| {
                column
                    .and_then(|i| record.get(i))
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
            };
            let run = record.get(run_col).unwrap_or_default().trim();
            if run.is_empty() || !seen.insert(run.to_string()) {
                bail!("duplicate or empty run ID '{run}' in sample metadata");
            }
            let (id, mut source) = if let Some(id) = value(declared_col) {
                (id.to_string(), IdentitySource::Declared)
            } else if let Some(id) = value(sample_col) {
                (id.to_ascii_uppercase(), IdentitySource::SampleAccession)
            } else {
                (run.to_string(), IdentitySource::RunFallback)
            };
            if let Some(declared_source) = value(source_col) {
                source = serde_json::from_value(serde_json::Value::String(declared_source.into()))
                    .context("invalid biological identity_source")?;
            }
            let metadata = headers
                .iter()
                .enumerate()
                .filter(|(_, name)| !reserved(name))
                .map(|(i, name)| {
                    let raw = record.get(i).unwrap_or_default();
                    let value = if required.contains(name) {
                        schema.normalize(name, raw)?
                    } else {
                        raw.to_string()
                    };
                    Ok((name.to_string(), value))
                })
                .collect::<Result<BTreeMap<_, _>>>()?;
            groups.entry(id).or_default().push(SourceRun {
                run: run.into(),
                identity_source: source,
                metadata,
            });
        }
        let mut samples = Vec::new();
        for (sample_id, mut rows) in groups {
            rows.sort_by(|a, b| a.run.cmp(&b.run));
            let source = rows[0].identity_source.clone();
            if rows.iter().any(|row| row.identity_source != source) {
                bail!("biological sample '{sample_id}' has conflicting identity sources");
            }
            let mut metadata = rows[0].metadata.clone();
            for (column, value) in &mut metadata {
                let agrees = rows
                    .iter()
                    .all(|row| row.metadata.get(column).is_some_and(|other| other == value));
                if required.contains(column) && (!agrees || value.trim().is_empty()) {
                    bail!(
                        "biological sample '{sample_id}' has conflicting or missing '{column}' values across contributing runs; cannot assign one model row"
                    );
                }
                if !agrees {
                    value.clear();
                }
            }
            samples.push(AnalysisSample {
                sample_id,
                identity_source: source,
                runs: rows.into_iter().map(|row| row.run).collect(),
                metadata,
            });
        }
        if samples.len() < 2 {
            bail!(
                "need at least two biological samples for analysis; found {}",
                samples.len()
            );
        }
        Ok(Self {
            schema_version: 1,
            samples,
        })
    }

    pub fn sample_ids(&self) -> Vec<String> {
        self.samples
            .iter()
            .map(|sample| sample.sample_id.clone())
            .collect()
    }

    pub fn model_table(&self) -> Result<(Vec<String>, Vec<Vec<String>>)> {
        let columns = self
            .samples
            .iter()
            .flat_map(|sample| sample.metadata.keys().cloned())
            .collect::<BTreeSet<_>>();
        let mut header = vec![
            "sample_id".into(),
            "biological_sample_id".into(),
            "identity_source".into(),
            "run_accessions".into(),
        ];
        header.extend(columns.iter().cloned());
        let rows = self
            .samples
            .iter()
            .map(|sample| {
                let mut row = vec![
                    sample.sample_id.clone(),
                    sample.sample_id.clone(),
                    sample.identity_source.as_str().into(),
                    serde_json::to_string(&sample.runs)?,
                ];
                row.extend(
                    columns
                        .iter()
                        .map(|name| sample.metadata.get(name).cloned().unwrap_or_default()),
                );
                Ok(row)
            })
            .collect::<Result<_>>()?;
        Ok((header, rows))
    }

    pub fn write_sheet(&self, path: &Path) -> Result<()> {
        let (header, rows) = self.model_table()?;
        let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;
        writer.write_record(header)?;
        for row in rows {
            writer.write_record(row)?;
        }
        writer.flush()?;
        Ok(())
    }
}

fn reserved(name: &str) -> bool {
    matches!(
        name,
        "run_accession"
            | "sample_id"
            | "biological_sample_id"
            | "identity_source"
            | "run_accessions"
    )
}

#[cfg(test)]
mod tests;
