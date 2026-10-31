mod overrides;
mod rows;
mod selectors;

pub use overrides::merge_metadata;
use overrides::{apply_overrides, load_override_tables};
use rows::parse_row_map;
use selectors::{filter_runs_by_selectors, parse_filter_predicates, split_inputs};

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
    load_or_init_column_map,
    publication::{ConflictPolicy, PreparedPublication, recover_view},
    write_samplesheet_schema,
};
use rnaa_core::traits::Resolver;
use rnaa_core::util::{normalize_url, now_rfc3339, sanitize_basename, write_json_pretty};
use serde_json::{Map, Value, json};

const ENA_FILE_REPORT_URL: &str = "https://www.ebi.ac.uk/ena/portal/api/filereport";

#[derive(Debug, Clone)]
pub struct EnaResolver {
    client: Client,
    fields: Vec<String>,
    metadata_conflicts: ConflictPolicy,
}

impl EnaResolver {
    pub fn new(fields: Vec<String>) -> Result<Self> {
        let client = Client::builder()
            .timeout(std::time::Duration::from_secs(90))
            .build()
            .context("failed to construct HTTP client")?;
        Ok(Self {
            client,
            fields,
            metadata_conflicts: ConflictPolicy::Fail,
        })
    }

    pub fn with_metadata_conflicts(mut self, policy: ConflictPolicy) -> Self {
        self.metadata_conflicts = policy;
        self
    }

    pub fn resolve_project(
        &self,
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<ResolvedProject> {
        let mut lease = rnaa_core::jobs::JobLease::try_acquire(db, "resolve", None, &[])?
            .context("another analysis or resolver owns this project; retry after it finishes")?;
        let result = self.resolve_project_owned(&lease.database(), paths, config, |input| {
            self.resolve_input(input)
        });
        match result {
            Ok(resolved) => {
                lease.finish(rnaa_core::jobs::JobOutcome::Succeeded, None)?;
                Ok(resolved)
            }
            Err(error) => {
                lease.finish(
                    rnaa_core::jobs::JobOutcome::Failed,
                    Some(&format!("{error:#}")),
                )?;
                Err(error)
            }
        }
    }

    fn resolve_project_owned(
        &self,
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        resolve_input: impl Fn(&InputRecord) -> Result<Vec<RunRecord>>,
    ) -> Result<ResolvedProject> {
        let expected_head = db.active_cohort()?.map(|cohort| cohort.cohort_id);
        let selector_contract = db.selector_contract()?;
        selector_contract.verify_config(config)?;
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

        let override_records = db.list_metadata_overrides()?;
        let override_inputs = override_records
            .iter()
            .map(|record| {
                ManifestArtifact::snapshot("METADATA_OVERRIDE", Path::new(&record.stored_path))
            })
            .collect::<Result<Vec<_>>>()?;
        let overrides = load_override_tables(&override_records)?;
        let mut runs = Vec::new();

        for input in selectors.query_inputs {
            let mut resolved = resolve_input(input)?;
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
        // Overlapping accessions commonly return the same run. Ambiguous
        // metadata must not be resolved by incidental provider response order.
        let mut unique: BTreeMap<String, RunRecord> = BTreeMap::new();
        for run in runs {
            if let Some(previous) = unique.get(&run.run_accession)
                && serde_json::to_value(previous)? != serde_json::to_value(&run)?
            {
                // Provider timestamps are local observation times, not metadata.
                let mut previous: RunRecord = previous.clone();
                previous.updated_at.clone_from(&run.updated_at);
                if serde_json::to_value(previous)? != serde_json::to_value(&run)? {
                    bail!(
                        "conflicting metadata for overlapping run {}",
                        run.run_accession
                    );
                }
            }
            unique.insert(run.run_accession.clone(), run);
        }
        let mut runs = unique.into_values().collect::<Vec<_>>();
        let metadata_columns = runs
            .iter()
            .flat_map(|run| run.metadata.keys().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();

        runs.sort_by(|left, right| left.run_accession.cmp(&right.run_accession));
        for input in &override_inputs {
            if !input.matches_file()? {
                bail!(
                    "metadata override changed during resolution: {}",
                    input.path
                );
            }
        }
        let disk = ProjectConfig::load(&paths.config_path)?;
        selector_contract.verify_config(&disk)?;
        let column_map = load_or_init_column_map(&paths.column_map_path(), &metadata_columns)?;
        let artifacts = db.list_artifacts()?;
        let metadata = PreparedPublication::prepare_with_sources(
            db,
            paths,
            &runs,
            &artifacts,
            &column_map,
            self.metadata_conflicts,
            &override_inputs,
        )?;
        let cohort = db.publish_resolved_cohort_with(
            &runs,
            &selector_contract,
            expected_head.as_deref(),
            |db| metadata.commit_catalog_preverified(db),
        )?;
        recover_view(db, paths, self.metadata_conflicts)?;
        db.append_event(
            "resolve",
            None,
            "resolved inputs to run metadata",
            json!({
                "inputs": inputs.iter().map(|input| input.input_id.clone()).collect::<Vec<_>>(),
                "run_count": runs.len(),
                "cohort_id": cohort.cohort_id,
                "include_samples": selectors.include_samples.iter().cloned().collect::<Vec<_>>(),
                "exclude_samples": selectors.exclude_samples.iter().cloned().collect::<Vec<_>>(),
                "include_filters": selectors.include_filters_raw,
                "exclude_filters": selectors.exclude_filters_raw
            }),
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
                "samplesheet_rows": runs.len(),
                "cohort_id": cohort.cohort_id,
                "include_samples": selectors.include_samples.iter().cloned().collect::<Vec<_>>(),
                "exclude_samples": selectors.exclude_samples.iter().cloned().collect::<Vec<_>>(),
                "include_filters": selectors.include_filters_raw,
                "exclude_filters": selectors.exclude_filters_raw
            }),
            tool_versions: BTreeMap::new(),
            input_artifacts: metadata.publication.sources.clone(),
            output_artifacts: vec![
                metadata.publication.generated.clone(),
                metadata.publication.merged.clone(),
                metadata.publication.column_map.clone(),
            ],
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

#[cfg(test)]
mod tests;
