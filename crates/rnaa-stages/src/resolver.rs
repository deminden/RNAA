use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;
use std::str::FromStr;

use anyhow::{Context, Result, anyhow, bail};
use csv::ReaderBuilder;
use reqwest::blocking::Client;
use reqwest::header::{ACCEPT, USER_AGENT};
use rnaa_core::config::{MetadataMergeStrategy, ProjectConfig};
use rnaa_core::db::Database;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{
    InputRecord, InputType, LibraryLayout, MetadataOverrideRecord, RemoteFile, RemoteFileKind,
    ResolvedProject, RunRecord,
};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::samplesheet::{
    load_or_init_column_map, write_samplesheet, write_samplesheet_schema,
};
use rnaa_core::traits::Resolver;
use rnaa_core::util::{normalize_url, now_rfc3339, sanitize_basename, write_json_pretty};
use serde_json::{Map, Value, json};

const ENA_FILE_REPORT_URL: &str = "https://www.ebi.ac.uk/ena/portal/api/filereport";

#[derive(Debug, Clone)]
pub struct EnaResolver {
    client: Client,
    fields: Vec<String>,
}

impl EnaResolver {
    pub fn new(fields: Vec<String>) -> Result<Self> {
        let client = Client::builder()
            .timeout(std::time::Duration::from_secs(90))
            .build()
            .context("failed to construct HTTP client")?;
        Ok(Self { client, fields })
    }

    pub fn resolve_project(
        &self,
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<ResolvedProject> {
        let inputs = db.list_inputs()?;
        if inputs.is_empty() {
            bail!(
                "no inputs found; add accessions/samples with `rnaa add --id ... --sample ...` first"
            );
        }
        let selectors = split_inputs(&inputs);
        if selectors.query_inputs.is_empty() {
            bail!(
                "no resolvable accession inputs found; add at least one run/study/project/sample accession"
            );
        }
        let include_filters = parse_filter_predicates(&selectors.include_filters_raw)?;
        let exclude_filters = parse_filter_predicates(&selectors.exclude_filters_raw)?;

        let overrides = load_override_tables(&db.list_metadata_overrides()?)?;
        let mut runs = Vec::new();

        for input in selectors.query_inputs {
            let mut resolved = self.resolve_input(input)?;
            for run in &mut resolved {
                apply_overrides(run, &overrides, config.resolver.merge_strategy);
            }
            runs.extend(resolved);
        }
        runs = filter_runs_by_selectors(
            runs,
            &selectors.include_samples,
            &selectors.exclude_samples,
            &include_filters,
            &exclude_filters,
        );
        if runs.is_empty() {
            bail!(
                "no runs matched active selectors (include_samples={}, exclude_samples={}, include_filters={}, exclude_filters={})",
                selectors.include_samples.len(),
                selectors.exclude_samples.len(),
                include_filters.len(),
                exclude_filters.len()
            );
        }
        let metadata_columns = runs
            .iter()
            .flat_map(|run| run.metadata.keys().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();

        runs.sort_by(|left, right| left.run_accession.cmp(&right.run_accession));
        db.upsert_runs(&runs)?;
        db.append_event(
            "resolve",
            None,
            "resolved inputs to run metadata",
            json!({
                "inputs": inputs.iter().map(|input| input.input_id.clone()).collect::<Vec<_>>(),
                "run_count": runs.len(),
                "include_samples": selectors.include_samples.iter().cloned().collect::<Vec<_>>(),
                "exclude_samples": selectors.exclude_samples.iter().cloned().collect::<Vec<_>>(),
                "include_filters": selectors.include_filters_raw,
                "exclude_filters": selectors.exclude_filters_raw
            }),
        )?;
        let column_map = load_or_init_column_map(&paths.column_map_path(), &metadata_columns)?;
        let artifacts = db.list_artifacts()?;
        let samplesheet_rows = write_samplesheet(
            &paths.samplesheet_path(),
            db.project_id(),
            &runs,
            &artifacts,
            &column_map,
        )?;
        write_samplesheet_schema(&paths.samplesheet_schema_path())?;

        let manifest_path = paths.manifest_path("resolve", db.project_id());
        let manifest = StageManifest {
            schema_version: 1,
            stage: "resolve".to_string(),
            id: db.project_id().to_string(),
            started_at: now_rfc3339(),
            finished_at: now_rfc3339(),
            parameters: json!({
                "provider": "ena",
                "fields": self.fields,
                "merge_strategy": config.resolver.merge_strategy,
                "samplesheet_rows": samplesheet_rows.len(),
                "include_samples": selectors.include_samples.iter().cloned().collect::<Vec<_>>(),
                "exclude_samples": selectors.exclude_samples.iter().cloned().collect::<Vec<_>>(),
                "include_filters": selectors.include_filters_raw,
                "exclude_filters": selectors.exclude_filters_raw
            }),
            tool_versions: BTreeMap::new(),
            input_artifacts: Vec::new(),
            output_artifacts: vec![ManifestArtifact {
                kind: "SAMPLESHEET".to_string(),
                path: paths.samplesheet_path().display().to_string(),
                checksum_type: "none".to_string(),
                checksum: String::new(),
                bytes: std::fs::metadata(paths.samplesheet_path())?.len(),
            }],
            notes: vec![
                "Resolver response fields are preserved in run.metadata with an `ena_` prefix."
                    .to_string(),
                "Canonical study_accession prefers secondary_study_accession when ENA provides one."
                    .to_string(),
            ],
        };
        write_json_pretty(&manifest_path, &manifest)?;

        Ok(ResolvedProject {
            runs,
            metadata_columns,
        })
    }

    pub fn resolve_input(&self, input: &InputRecord) -> Result<Vec<RunRecord>> {
        let response = self.fetch_input(&input.input_id)?;
        let rows = Self::parse_report_rows(&response, &input.project_id)?;
        if rows.is_empty() {
            bail!("ENA returned no runs for {}", input.input_id);
        }
        Ok(rows)
    }

    pub fn fetch_input(&self, accession: &str) -> Result<String> {
        let mut url = reqwest::Url::parse(ENA_FILE_REPORT_URL)?;
        url.query_pairs_mut()
            .append_pair("accession", accession)
            .append_pair("result", "read_run")
            .append_pair("format", "json")
            .append_pair("fields", &self.fields.join(","));
        self.client
            .get(url)
            .header(USER_AGENT, "RNAA/0.1.0")
            .header(ACCEPT, "application/json")
            .send()
            .and_then(reqwest::blocking::Response::error_for_status)
            .context("ENA filereport request failed")?
            .text()
            .context("failed to read ENA response body")
    }

    pub fn parse_report_rows(body: &str, project_id: &str) -> Result<Vec<RunRecord>> {
        let rows: Vec<Map<String, Value>> =
            serde_json::from_str(body).context("failed to parse ENA JSON response")?;
        let mut parsed = Vec::new();
        for row in rows {
            parsed.push(parse_row_map(project_id, row)?);
        }
        Ok(parsed)
    }
}

impl Resolver for EnaResolver {
    fn resolve(&self, inputs: &[InputRecord]) -> Result<ResolvedProject> {
        let selectors = split_inputs(inputs);
        let include_filters = parse_filter_predicates(&selectors.include_filters_raw)?;
        let exclude_filters = parse_filter_predicates(&selectors.exclude_filters_raw)?;
        let mut runs = Vec::new();
        for input in selectors.query_inputs {
            let resolved = self.resolve_input(input)?;
            runs.extend(resolved);
        }
        runs = filter_runs_by_selectors(
            runs,
            &selectors.include_samples,
            &selectors.exclude_samples,
            &include_filters,
            &exclude_filters,
        );
        let metadata_columns = runs
            .iter()
            .flat_map(|run| run.metadata.keys().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        Ok(ResolvedProject {
            runs,
            metadata_columns,
        })
    }
}

#[derive(Debug, Default)]
struct InputSelectors<'a> {
    query_inputs: Vec<&'a InputRecord>,
    include_samples: BTreeSet<String>,
    exclude_samples: BTreeSet<String>,
    include_filters_raw: Vec<String>,
    exclude_filters_raw: Vec<String>,
}

fn split_inputs(inputs: &[InputRecord]) -> InputSelectors<'_> {
    let mut query_inputs = Vec::new();
    let mut include_samples = BTreeSet::new();
    let mut exclude_samples = BTreeSet::new();
    let mut include_filters = Vec::new();
    let mut exclude_filters = Vec::new();
    for input in inputs {
        match input.input_type {
            InputType::SampleInclude => {
                include_samples.insert(normalize_sample_selector(&input.input_id));
            }
            InputType::SampleExclude => {
                exclude_samples.insert(normalize_sample_selector(&input.input_id));
            }
            InputType::FilterInclude => include_filters.push(input.input_id.clone()),
            InputType::FilterExclude => exclude_filters.push(input.input_id.clone()),
            _ => query_inputs.push(input),
        }
    }
    InputSelectors {
        query_inputs,
        include_samples,
        exclude_samples,
        include_filters_raw: include_filters,
        exclude_filters_raw: exclude_filters,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PredicateOp {
    Eq,
    Neq,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FieldPredicate {
    field: String,
    op: PredicateOp,
    value: String,
}

fn parse_filter_predicates(raw: &[String]) -> Result<Vec<FieldPredicate>> {
    raw.iter()
        .map(|item| parse_filter_predicate(item))
        .collect()
}

fn parse_filter_predicate(spec: &str) -> Result<FieldPredicate> {
    let trimmed = spec.trim();
    if let Some((field, value)) = trimmed.split_once("!=") {
        return build_filter_predicate(field, PredicateOp::Neq, value, spec);
    }
    if let Some((field, value)) = trimmed.split_once('=') {
        return build_filter_predicate(field, PredicateOp::Eq, value, spec);
    }
    bail!(
        "invalid filter predicate '{}'; expected 'field=value' or 'field!=value'",
        spec
    )
}

fn build_filter_predicate(
    field: &str,
    op: PredicateOp,
    value: &str,
    original: &str,
) -> Result<FieldPredicate> {
    let field = field.trim().to_ascii_lowercase();
    let value = value.trim().to_ascii_lowercase();
    if field.is_empty() || value.is_empty() {
        bail!(
            "invalid filter predicate '{}'; field and value must be non-empty",
            original
        );
    }
    Ok(FieldPredicate { field, op, value })
}

fn filter_runs_by_selectors(
    runs: Vec<RunRecord>,
    include_samples: &BTreeSet<String>,
    exclude_samples: &BTreeSet<String>,
    include_filters: &[FieldPredicate],
    exclude_filters: &[FieldPredicate],
) -> Vec<RunRecord> {
    runs.into_iter()
        .filter(|run| {
            let sample = run
                .sample_accession
                .as_deref()
                .map(normalize_sample_selector)
                .unwrap_or_default();
            (include_samples.is_empty() || include_samples.contains(&sample))
                && (exclude_samples.is_empty() || !exclude_samples.contains(&sample))
                && include_filters
                    .iter()
                    .all(|predicate| run_matches_predicate(run, predicate))
                && !exclude_filters
                    .iter()
                    .any(|predicate| run_matches_predicate(run, predicate))
        })
        .collect()
}

fn normalize_sample_selector(value: &str) -> String {
    value.trim().to_ascii_uppercase()
}

fn run_matches_predicate(run: &RunRecord, predicate: &FieldPredicate) -> bool {
    let Some(field_value) = run_field_value(run, &predicate.field) else {
        return false;
    };
    let normalized = field_value.trim().to_ascii_lowercase();
    match predicate.op {
        PredicateOp::Eq => normalized == predicate.value,
        PredicateOp::Neq => normalized != predicate.value,
    }
}

fn run_field_value<'a>(run: &'a RunRecord, field: &str) -> Option<&'a str> {
    match field {
        "run_accession" => Some(run.run_accession.as_str()),
        "sample_accession" => run.sample_accession.as_deref(),
        "study_accession" => Some(run.study_accession.as_str()),
        "source_study_accession" => run.source_study_accession.as_deref(),
        "experiment_accession" => run.experiment_accession.as_deref(),
        "library_layout" => Some(run.library_layout.as_str()),
        "instrument_platform" => run.instrument_platform.as_deref(),
        "instrument_model" => run.instrument_model.as_deref(),
        _ => run.metadata.get(field).map(String::as_str),
    }
}

#[derive(Debug, Clone)]
struct OverrideRow {
    run_accession: Option<String>,
    sample_accession: Option<String>,
    values: BTreeMap<String, String>,
    merge_strategy: MetadataMergeStrategy,
}

fn parse_row_map(project_id: &str, row: Map<String, Value>) -> Result<RunRecord> {
    let run_accession = required_string(&row, "run_accession")?;
    let primary_study = opt_string(&row, "study_accession");
    let secondary_study = opt_string(&row, "secondary_study_accession");
    let study_accession = secondary_study
        .clone()
        .or(primary_study.clone())
        .ok_or_else(|| anyhow!("ENA row missing study accession for {}", run_accession))?;
    let source_study_accession = match (&primary_study, &secondary_study) {
        (Some(primary), Some(_)) => Some(primary.clone()),
        _ => None,
    };
    let library_layout = opt_string(&row, "library_layout")
        .as_deref()
        .map(LibraryLayout::from_str)
        .transpose()?
        .unwrap_or_default();
    let fastq_files = expand_remote_files(
        &opt_string(&row, "fastq_ftp").unwrap_or_default(),
        &opt_string(&row, "fastq_md5").unwrap_or_default(),
        &opt_string(&row, "fastq_bytes").unwrap_or_default(),
        library_layout,
        false,
    );
    let sra_files = expand_remote_files(
        &opt_string(&row, "submitted_ftp").unwrap_or_default(),
        &opt_string(&row, "submitted_md5").unwrap_or_default(),
        &opt_string(&row, "submitted_bytes").unwrap_or_default(),
        library_layout,
        true,
    );
    let remote_files = if fastq_files.is_empty() {
        sra_files
    } else {
        fastq_files.into_iter().chain(sra_files).collect()
    };

    let metadata = row
        .iter()
        .filter_map(|(key, value)| {
            if matches!(
                key.as_str(),
                "run_accession"
                    | "study_accession"
                    | "secondary_study_accession"
                    | "sample_accession"
                    | "experiment_accession"
                    | "library_layout"
                    | "instrument_platform"
                    | "instrument_model"
                    | "fastq_ftp"
                    | "fastq_md5"
                    | "fastq_bytes"
                    | "submitted_ftp"
                    | "submitted_md5"
                    | "submitted_bytes"
            ) {
                return None;
            }
            let stringified = json_value_to_string(value);
            (!stringified.is_empty()).then(|| (format!("ena_{key}"), stringified))
        })
        .chain(
            primary_study
                .iter()
                .map(|value| ("ena_primary_study_accession".to_string(), value.to_string())),
        )
        .chain(secondary_study.iter().map(|value| {
            (
                "ena_secondary_study_accession".to_string(),
                value.to_string(),
            )
        }))
        .collect::<BTreeMap<_, _>>();

    Ok(RunRecord {
        project_id: project_id.to_string(),
        study_accession,
        source_study_accession,
        run_accession,
        sample_accession: opt_string(&row, "sample_accession")
            .or_else(|| opt_string(&row, "secondary_sample_accession")),
        experiment_accession: opt_string(&row, "experiment_accession"),
        library_layout,
        instrument_platform: opt_string(&row, "instrument_platform"),
        instrument_model: opt_string(&row, "instrument_model"),
        remote_files,
        metadata,
        state: rnaa_core::state::RunState::Resolved,
        last_error: None,
        updated_at: now_rfc3339(),
    })
}

fn expand_remote_files(
    urls_text: &str,
    md5_text: &str,
    bytes_text: &str,
    layout: LibraryLayout,
    submitted: bool,
) -> Vec<RemoteFile> {
    let urls = split_field(urls_text);
    let md5s = split_field(md5_text);
    let bytes = split_field(bytes_text);

    urls.into_iter()
        .enumerate()
        .filter_map(|(index, url)| {
            if url.is_empty() {
                return None;
            }
            let normalized = normalize_url(&url);
            let basename = normalized
                .rsplit('/')
                .next()
                .map(sanitize_basename)
                .unwrap_or_else(|| format!("file_{index}"));
            let kind = if submitted && basename.ends_with(".sra") {
                RemoteFileKind::Sra
            } else if layout == LibraryLayout::Paired {
                if basename.contains("_1.") || basename.contains("R1") {
                    RemoteFileKind::FastqR1
                } else if basename.contains("_2.") || basename.contains("R2") {
                    RemoteFileKind::FastqR2
                } else if submitted {
                    RemoteFileKind::Submitted
                } else {
                    RemoteFileKind::Fastq
                }
            } else if submitted && basename.ends_with(".sra") {
                RemoteFileKind::Sra
            } else {
                RemoteFileKind::Fastq
            };

            Some(RemoteFile {
                url: normalized,
                md5: md5s.get(index).cloned().filter(|value| !value.is_empty()),
                bytes: bytes.get(index).and_then(|value| value.parse::<u64>().ok()),
                kind,
                basename,
            })
        })
        .collect()
}

fn split_field(value: &str) -> Vec<String> {
    value
        .split(';')
        .map(str::trim)
        .filter(|field| !field.is_empty())
        .map(ToString::to_string)
        .collect()
}

fn required_string(row: &Map<String, Value>, key: &str) -> Result<String> {
    opt_string(row, key).ok_or_else(|| anyhow!("ENA row missing required field {key}"))
}

fn opt_string(row: &Map<String, Value>, key: &str) -> Option<String> {
    row.get(key)
        .map(json_value_to_string)
        .filter(|value| !value.is_empty())
}

fn json_value_to_string(value: &Value) -> String {
    match value {
        Value::Null => String::new(),
        Value::String(text) => text.trim().to_string(),
        Value::Number(number) => number.to_string(),
        Value::Bool(boolean) => boolean.to_string(),
        other => other.to_string(),
    }
}

fn load_override_tables(overrides: &[MetadataOverrideRecord]) -> Result<Vec<OverrideRow>> {
    let mut rows = Vec::new();
    for override_record in overrides {
        let path = Path::new(&override_record.stored_path);
        rows.extend(load_override_table(path, override_record.merge_strategy)?);
    }
    Ok(rows)
}

fn load_override_table(
    path: &Path,
    merge_strategy: MetadataMergeStrategy,
) -> Result<Vec<OverrideRow>> {
    let content =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
    let delimiter = if path.extension().and_then(|ext| ext.to_str()) == Some("csv") {
        b','
    } else if content
        .lines()
        .next()
        .is_some_and(|line| line.contains('\t'))
    {
        b'\t'
    } else {
        b','
    };

    let mut reader = ReaderBuilder::new()
        .delimiter(delimiter)
        .from_reader(content.as_bytes());
    let headers = reader.headers()?.clone();
    let run_idx = headers.iter().position(|name| name == "run_accession");
    let sample_idx = headers.iter().position(|name| name == "sample_accession");
    if run_idx.is_none() && sample_idx.is_none() {
        bail!(
            "metadata override {} must contain run_accession or sample_accession",
            path.display()
        );
    }

    let mut overrides = Vec::new();
    for record in reader.records() {
        let record = record?;
        let mut values = BTreeMap::new();
        for (index, field) in record.iter().enumerate() {
            let header = headers.get(index).unwrap_or_default();
            if field.trim().is_empty() || header == "run_accession" || header == "sample_accession"
            {
                continue;
            }
            values.insert(header.to_string(), field.to_string());
        }
        overrides.push(OverrideRow {
            run_accession: run_idx
                .and_then(|index| record.get(index))
                .map(ToString::to_string),
            sample_accession: sample_idx
                .and_then(|index| record.get(index))
                .map(ToString::to_string),
            values,
            merge_strategy,
        });
    }
    Ok(overrides)
}

fn apply_overrides(
    run: &mut RunRecord,
    overrides: &[OverrideRow],
    _default_strategy: MetadataMergeStrategy,
) {
    for override_row in overrides.iter().filter(|override_row| {
        override_row.run_accession.as_deref() == Some(run.run_accession.as_str())
            || override_row.sample_accession.as_deref() == run.sample_accession.as_deref()
    }) {
        merge_metadata(
            &mut run.metadata,
            &override_row.values,
            override_row.merge_strategy,
        );
    }
}

pub fn merge_metadata(
    target: &mut BTreeMap<String, String>,
    incoming: &BTreeMap<String, String>,
    strategy: MetadataMergeStrategy,
) {
    for (key, value) in incoming {
        match strategy {
            MetadataMergeStrategy::Override | MetadataMergeStrategy::PreferLocal => {
                target.insert(key.clone(), value.clone());
            }
            MetadataMergeStrategy::PreferRemote => {
                target.entry(key.clone()).or_insert_with(|| value.clone());
            }
            MetadataMergeStrategy::UnionWithPrefix => {
                if let Some(existing) = target.get(key).cloned() {
                    if existing != *value {
                        target.insert(format!("remote_{key}"), existing);
                        target.insert(format!("local_{key}"), value.clone());
                    }
                } else {
                    target.insert(key.clone(), value.clone());
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rnaa_core::config::MetadataMergeStrategy;
    use rnaa_core::model::InputType;

    #[test]
    fn parses_ena_json_fixture() {
        let fixture = include_str!("../../../fixtures/ena_srp114962.json");
        let runs = EnaResolver::parse_report_rows(fixture, "project-1").unwrap();
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0].study_accession, "SRP114962");
        assert_eq!(
            runs[0].source_study_accession.as_deref(),
            Some("PRJNA396019")
        );
        assert!(
            runs[0]
                .remote_files
                .iter()
                .any(|file| file.kind == RemoteFileKind::FastqR1
                    || file.kind == RemoteFileKind::Fastq)
        );
    }

    #[test]
    fn merges_metadata_prefer_local() {
        let mut remote = BTreeMap::from([
            ("condition".to_string(), "remote".to_string()),
            ("center".to_string(), "ebi".to_string()),
        ]);
        let incoming = BTreeMap::from([
            ("condition".to_string(), "local".to_string()),
            ("batch".to_string(), "b1".to_string()),
        ]);
        merge_metadata(&mut remote, &incoming, MetadataMergeStrategy::PreferLocal);
        assert_eq!(remote.get("condition").unwrap(), "local");
        assert_eq!(remote.get("batch").unwrap(), "b1");
    }

    #[test]
    fn merges_metadata_union_with_prefix() {
        let mut remote = BTreeMap::from([("condition".to_string(), "remote".to_string())]);
        let incoming = BTreeMap::from([("condition".to_string(), "local".to_string())]);
        merge_metadata(
            &mut remote,
            &incoming,
            MetadataMergeStrategy::UnionWithPrefix,
        );
        assert_eq!(remote.get("remote_condition").unwrap(), "remote");
        assert_eq!(remote.get("local_condition").unwrap(), "local");
    }

    #[test]
    fn filters_runs_by_sample_selectors() {
        let fixture = include_str!("../../../fixtures/ena_srp114962.json");
        let runs = EnaResolver::parse_report_rows(fixture, "project-1").unwrap();
        let include_samples = BTreeSet::from_iter([String::from("SAMN09909193")]);
        let exclude_samples = BTreeSet::new();
        let include_filters = Vec::new();
        let exclude_filters = Vec::new();
        let included = filter_runs_by_selectors(
            runs.clone(),
            &include_samples,
            &exclude_samples,
            &include_filters,
            &exclude_filters,
        );
        assert_eq!(included.len(), 1);
        assert_eq!(
            included[0].sample_accession.as_deref(),
            Some("SAMN09909193")
        );

        let include_samples = BTreeSet::new();
        let exclude_samples = BTreeSet::from_iter([String::from("SAMN09909193")]);
        let excluded = filter_runs_by_selectors(
            runs,
            &include_samples,
            &exclude_samples,
            &include_filters,
            &exclude_filters,
        );
        assert_eq!(excluded.len(), 1);
        assert_eq!(
            excluded[0].sample_accession.as_deref(),
            Some("SAMN09909196")
        );
    }

    #[test]
    fn splits_selector_inputs_from_query_inputs() {
        let inputs = vec![
            InputRecord {
                project_id: "p".to_string(),
                input_id: "SRP1".to_string(),
                input_type: InputType::Study,
                added_at: "now".to_string(),
            },
            InputRecord {
                project_id: "p".to_string(),
                input_id: "samn1".to_string(),
                input_type: InputType::SampleInclude,
                added_at: "now".to_string(),
            },
            InputRecord {
                project_id: "p".to_string(),
                input_id: "samn2".to_string(),
                input_type: InputType::SampleExclude,
                added_at: "now".to_string(),
            },
            InputRecord {
                project_id: "p".to_string(),
                input_id: "instrument_platform=illumina".to_string(),
                input_type: InputType::FilterInclude,
                added_at: "now".to_string(),
            },
            InputRecord {
                project_id: "p".to_string(),
                input_id: "sample_accession!=SAMN0".to_string(),
                input_type: InputType::FilterExclude,
                added_at: "now".to_string(),
            },
        ];
        let selectors = split_inputs(&inputs);
        assert_eq!(selectors.query_inputs.len(), 1);
        assert_eq!(selectors.query_inputs[0].input_id, "SRP1");
        assert!(selectors.include_samples.contains("SAMN1"));
        assert!(selectors.exclude_samples.contains("SAMN2"));
        assert_eq!(
            selectors.include_filters_raw,
            vec!["instrument_platform=illumina"]
        );
        assert_eq!(
            selectors.exclude_filters_raw,
            vec!["sample_accession!=SAMN0"]
        );
    }

    #[test]
    fn parses_and_applies_where_predicates() {
        let fixture = include_str!("../../../fixtures/ena_srp114962.json");
        let runs = EnaResolver::parse_report_rows(fixture, "project-1").unwrap();
        let include_filters =
            parse_filter_predicates(&[String::from("sample_accession=SAMN09909193")]).unwrap();
        let exclude_filters = Vec::new();
        let selected = filter_runs_by_selectors(
            runs.clone(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            &include_filters,
            &exclude_filters,
        );
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].run_accession, "SRR7750979");

        let include_filters = Vec::new();
        let exclude_filters =
            parse_filter_predicates(&[String::from("sample_accession=SAMN09909193")]).unwrap();
        let selected = filter_runs_by_selectors(
            runs,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &include_filters,
            &exclude_filters,
        );
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].run_accession, "SRR7750980");
    }
}
