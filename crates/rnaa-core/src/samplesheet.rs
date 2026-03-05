use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use csv::{ReaderBuilder, WriterBuilder};
use serde::Serialize;

use crate::model::{ArtifactKind, ArtifactRecord, RunRecord, SamplesheetRow};

const STABLE_COLUMNS: &[&str] = &[
    "project_id",
    "study_accession",
    "run_accession",
    "sample_accession",
    "experiment_accession",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "read1_path",
    "read2_path",
    "condition",
    "batch",
];

#[derive(Debug, Clone, Serialize)]
struct SamplesheetSchema {
    schema_version: u32,
    format: &'static str,
    stable_columns: Vec<&'static str>,
    notes: Vec<&'static str>,
}

#[derive(Debug, Clone)]
pub struct ColumnMap {
    pub source_to_canonical: BTreeMap<String, String>,
}

impl ColumnMap {
    pub fn canonical_for(&self, source: &str) -> Option<&str> {
        self.source_to_canonical.get(source).map(String::as_str)
    }
}

pub fn load_or_init_column_map(path: &Path, metadata_columns: &[String]) -> Result<ColumnMap> {
    if !path.exists() {
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .with_context(|| format!("failed to create {}", path.display()))?;
        writer.write_record(["source_column", "canonical_field", "suggestion", "notes"])?;
        let mut sorted = metadata_columns.to_vec();
        sorted.sort();
        for column in sorted {
            let suggestion = suggest_canonical_field(&column).unwrap_or_default();
            let note = if suggestion.is_empty() {
                "edit canonical_field to map this column into RNAA"
            } else {
                "suggested from column name"
            };
            writer.write_record([column.as_str(), "", suggestion, note])?;
        }
        writer.flush()?;
    }

    let mut source_to_canonical = BTreeMap::new();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read {}", path.display()))?;
    for row in reader.records() {
        let row = row?;
        let source = row.get(0).unwrap_or("").trim();
        let canonical = row.get(1).unwrap_or("").trim();
        if !source.is_empty() && !canonical.is_empty() {
            source_to_canonical.insert(source.to_string(), canonical.to_string());
        }
    }
    Ok(ColumnMap {
        source_to_canonical,
    })
}

pub fn write_samplesheet_schema(path: &Path) -> Result<()> {
    let schema = SamplesheetSchema {
        schema_version: 1,
        format: "tsv",
        stable_columns: STABLE_COLUMNS.to_vec(),
        notes: vec![
            "Rows are keyed by run_accession.",
            "Additional metadata columns are preserved verbatim after the stable columns.",
            "condition and batch may be supplied directly or derived through metadata/column_map.tsv.",
        ],
    };
    let json = serde_json::to_string_pretty(&schema)?;
    fs::write(path, json).with_context(|| format!("failed to write {}", path.display()))
}

pub fn write_samplesheet(
    path: &Path,
    project_id: &str,
    runs: &[RunRecord],
    artifacts: &[ArtifactRecord],
    column_map: &ColumnMap,
) -> Result<Vec<SamplesheetRow>> {
    let artifacts_by_run = group_artifacts_by_run(artifacts);
    let rows = build_samplesheet_rows(project_id, runs, &artifacts_by_run, column_map);

    let mut extra_columns = BTreeSet::new();
    for row in &rows {
        extra_columns.extend(row.metadata.keys().cloned());
    }

    let mut header: Vec<String> = STABLE_COLUMNS.iter().map(ToString::to_string).collect();
    header.extend(extra_columns.iter().cloned());

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    writer.write_record(&header)?;

    for row in &rows {
        let mut record = vec![
            row.project_id.clone(),
            row.study_accession.clone(),
            row.run_accession.clone(),
            row.sample_accession.clone().unwrap_or_default(),
            row.experiment_accession.clone().unwrap_or_default(),
            row.library_layout.clone(),
            row.instrument_platform.clone().unwrap_or_default(),
            row.instrument_model.clone().unwrap_or_default(),
            row.read1_path.clone().unwrap_or_default(),
            row.read2_path.clone().unwrap_or_default(),
            row.condition.clone().unwrap_or_default(),
            row.batch.clone().unwrap_or_default(),
        ];
        record.extend(
            extra_columns
                .iter()
                .map(|column| row.metadata.get(column).cloned().unwrap_or_default()),
        );
        writer.write_record(record)?;
    }
    writer.flush()?;

    Ok(rows)
}

pub fn build_samplesheet_rows(
    project_id: &str,
    runs: &[RunRecord],
    artifacts_by_run: &HashMap<String, Vec<ArtifactRecord>>,
    column_map: &ColumnMap,
) -> Vec<SamplesheetRow> {
    let mut rows = Vec::with_capacity(runs.len());

    for run in runs {
        let mut metadata = run.metadata.clone();
        let mut condition = metadata.remove("condition");
        let mut batch = metadata.remove("batch");

        for (source, canonical) in &column_map.source_to_canonical {
            if let Some(value) = metadata.get(source).cloned() {
                match canonical.as_str() {
                    "condition" if condition.is_none() => condition = Some(value),
                    "batch" if batch.is_none() => batch = Some(value),
                    _ => {}
                }
            }
        }

        let read1_path = artifacts_by_run
            .get(&run.run_accession)
            .and_then(|artifacts| {
                artifacts
                    .iter()
                    .find(|artifact| {
                        artifact.kind == ArtifactKind::FastqR1
                            || artifact.kind == ArtifactKind::FastqSingle
                    })
                    .map(|artifact| artifact.path.clone())
            });

        let read2_path = artifacts_by_run
            .get(&run.run_accession)
            .and_then(|artifacts| {
                artifacts
                    .iter()
                    .find(|artifact| artifact.kind == ArtifactKind::FastqR2)
                    .map(|artifact| artifact.path.clone())
            });

        rows.push(SamplesheetRow {
            project_id: project_id.to_string(),
            study_accession: run.study_accession.clone(),
            run_accession: run.run_accession.clone(),
            sample_accession: run.sample_accession.clone(),
            experiment_accession: run.experiment_accession.clone(),
            library_layout: run.library_layout.as_str().to_string(),
            instrument_platform: run.instrument_platform.clone(),
            instrument_model: run.instrument_model.clone(),
            read1_path,
            read2_path,
            condition,
            batch,
            metadata,
        });
    }

    rows.sort_by(|left, right| left.run_accession.cmp(&right.run_accession));
    rows
}

fn group_artifacts_by_run(artifacts: &[ArtifactRecord]) -> HashMap<String, Vec<ArtifactRecord>> {
    let mut grouped: HashMap<String, Vec<ArtifactRecord>> = HashMap::new();
    for artifact in artifacts {
        if let Some(run_accession) = &artifact.run_accession {
            grouped
                .entry(run_accession.clone())
                .or_default()
                .push(artifact.clone());
        }
    }
    grouped
}

fn suggest_canonical_field(column: &str) -> Option<&'static str> {
    let normalized = column.to_ascii_lowercase();
    if normalized == "condition"
        || normalized.contains("condition")
        || normalized.contains("group")
        || normalized.contains("treatment")
    {
        return Some("condition");
    }
    if normalized == "batch"
        || normalized.contains("batch")
        || normalized.contains("lane")
        || normalized.contains("center")
    {
        return Some("batch");
    }
    None
}
