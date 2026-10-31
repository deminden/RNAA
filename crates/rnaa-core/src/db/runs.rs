//! Run ingestion, normalized metadata and legacy metadata backfill.
use super::Database;
use super::from_json_text;
use crate::config::MetadataMergeStrategy;
use crate::model::{InputRecord, InputType, LibraryLayout, MetadataOverrideRecord, RunRecord};
use crate::state::RunState;
use crate::util::now_rfc3339;
use anyhow::{Context, Result};
use rusqlite::{Connection, OptionalExtension, Row, params};
use std::path::Path;
use std::str::FromStr;

impl Database {
    pub fn add_input(&self, input_id: &str, input_type: InputType) -> Result<()> {
        self.write(|conn| {
            let changed = conn.execute(
                "INSERT OR IGNORE INTO inputs (project_id, input_id, input_type, added_at) VALUES (?1, ?2, ?3, ?4)",
                params![self.project_id, input_id, input_type.as_str(), now_rfc3339()],
            )
            .with_context(|| format!("failed to add input {input_id}"))?;
            if changed > 0 {
                crate::cohorts::invalidate_analysis(conn, self.project_id(), "cohort selectors changed; resolve again")?;
            }
            Ok(())
        })
    }

    pub fn list_inputs(&self) -> Result<Vec<InputRecord>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT project_id, input_id, input_type, added_at FROM inputs ORDER BY added_at, input_id",
            )?;
            let rows = stmt.query_map([], |row| {
                Ok((
                    row.get::<_, String>(0)?,
                    row.get::<_, String>(1)?,
                    row.get::<_, String>(2)?,
                    row.get::<_, String>(3)?,
                ))
            })?;

            let mut inputs = Vec::new();
            for row in rows {
                let (project_id, input_id, input_type, added_at) = row?;
                inputs.push(InputRecord {
                    project_id,
                    input_id,
                    input_type: InputType::from_str(&input_type)?,
                    added_at,
                });
            }
            Ok(inputs)
        })
    }

    pub fn add_metadata_override(
        &self,
        source_path: &Path,
        stored_path: &Path,
        merge_strategy: MetadataMergeStrategy,
    ) -> Result<()> {
        self.write(|conn| {
            conn.execute(
                "INSERT INTO metadata_overrides (project_id, source_path, stored_path, merge_strategy, added_at) VALUES (?1, ?2, ?3, ?4, ?5)",
                params![
                    self.project_id,
                    source_path.display().to_string(),
                    stored_path.display().to_string(),
                    serde_json::to_string(&merge_strategy)?,
                    now_rfc3339()
                ],
            )
            .with_context(|| {
                format!(
                    "failed to register metadata override {}",
                    stored_path.display()
                )
            })?;
            crate::cohorts::invalidate_analysis(conn, self.project_id(), "metadata overrides changed; resolve again")?;
            Ok(())
        })
    }

    pub fn list_metadata_overrides(&self) -> Result<Vec<MetadataOverrideRecord>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT id, project_id, source_path, stored_path, merge_strategy, added_at FROM metadata_overrides ORDER BY id",
            )?;
            let rows = stmt.query_map([], |row| {
                Ok((
                    row.get::<_, i64>(0)?,
                    row.get::<_, String>(1)?,
                    row.get::<_, String>(2)?,
                    row.get::<_, String>(3)?,
                    row.get::<_, String>(4)?,
                    row.get::<_, String>(5)?,
                ))
            })?;

            let mut overrides = Vec::new();
            for row in rows {
                let (id, project_id, source_path, stored_path, merge_strategy, added_at) = row?;
                overrides.push(MetadataOverrideRecord {
                    id,
                    project_id,
                    source_path,
                    stored_path,
                    merge_strategy: serde_json::from_str(&merge_strategy)
                        .context("failed to parse stored merge strategy")?,
                    added_at,
                });
            }
            Ok(overrides)
        })
    }

    pub fn upsert_runs(&self, runs: &[RunRecord]) -> Result<()> {
        self.write(|conn| {
            for run in runs {
                Self::upsert_run_tx(conn, run)?;
            }
            Ok(())
        })
    }

    pub fn upsert_run(&self, run: &RunRecord) -> Result<()> {
        self.write(|conn| Self::upsert_run_tx(conn, run))
    }

    pub fn get_run(&self, run_accession: &str) -> Result<Option<RunRecord>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT project_id, study_accession, source_study_accession, run_accession, sample_accession, experiment_accession, library_layout, instrument_platform, instrument_model, remote_files_json, metadata_json, state, last_error, updated_at FROM runs WHERE run_accession = ?1",
            )?;
            let raw = stmt
                .query_row(params![run_accession], parse_raw_run)
                .optional()?;
            let Some(raw) = raw else {
                return Ok(None);
            };
            let mut run: RunRecord = raw.try_into()?;
            let metadata = Self::load_run_metadata_map(conn, &run.run_accession)?;
            if !metadata.is_empty() {
                run.metadata = metadata;
            }
            Ok(Some(run))
        })
    }

    pub fn list_runs(&self) -> Result<Vec<RunRecord>> {
        self.list_runs_filtered(None)
    }

    pub fn list_runs_in_states(&self, states: &[RunState]) -> Result<Vec<RunRecord>> {
        if states.is_empty() {
            return Ok(Vec::new());
        }
        let placeholders = std::iter::repeat_n("?", states.len())
            .collect::<Vec<_>>()
            .join(",");
        let sql = format!(
            "SELECT project_id, study_accession, source_study_accession, run_accession, sample_accession, experiment_accession, library_layout, instrument_platform, instrument_model, remote_files_json, metadata_json, state, last_error, updated_at FROM runs WHERE state IN ({placeholders}) ORDER BY run_accession"
        );
        self.read(|conn| {
            let mut stmt = conn.prepare(&sql)?;
            let state_strings = states
                .iter()
                .map(|state| state.as_str().to_string())
                .collect::<Vec<_>>();
            let rows = stmt.query_map(
                rusqlite::params_from_iter(state_strings.iter()),
                parse_raw_run,
            )?;
            let metadata_by_run = Self::load_all_run_metadata_maps(conn)?;
            let mut runs = Vec::new();
            for row in rows {
                let mut run: RunRecord = row?.try_into()?;
                if let Some(metadata) = metadata_by_run.get(&run.run_accession) {
                    run.metadata = metadata.clone();
                }
                runs.push(run);
            }
            Ok(runs)
        })
    }

    pub fn list_metadata_columns(&self) -> Result<Vec<String>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT DISTINCT field_name
                 FROM run_metadata
                 WHERE project_id = ?1
                 ORDER BY field_name",
            )?;
            let rows = stmt.query_map(params![self.project_id], |row| row.get::<_, String>(0))?;
            let mut columns = Vec::new();
            for row in rows {
                columns.push(row?);
            }
            Ok(columns)
        })
    }

    pub fn set_run_state(
        &self,
        run_accession: &str,
        state: RunState,
        last_error: Option<&str>,
    ) -> Result<()> {
        self.write(|conn| {
            conn.execute(
                "UPDATE runs SET state = ?1, last_error = ?2, updated_at = ?3 WHERE run_accession = ?4",
                params![state.as_str(), last_error, now_rfc3339(), run_accession],
            )
            .with_context(|| format!("failed to update state for {run_accession}"))?;
            Ok(())
        })
    }

    fn list_runs_filtered(&self, _filter: Option<&str>) -> Result<Vec<RunRecord>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT project_id, study_accession, source_study_accession, run_accession, sample_accession, experiment_accession, library_layout, instrument_platform, instrument_model, remote_files_json, metadata_json, state, last_error, updated_at FROM runs ORDER BY run_accession",
            )?;
            let rows = stmt.query_map([], parse_raw_run)?;
            let metadata_by_run = Self::load_all_run_metadata_maps(conn)?;
            let mut runs = Vec::new();
            for row in rows {
                let mut run: RunRecord = row?.try_into()?;
                if let Some(metadata) = metadata_by_run.get(&run.run_accession) {
                    run.metadata = metadata.clone();
                }
                runs.push(run);
            }
            Ok(runs)
        })
    }

    fn upsert_run_tx(conn: &Connection, run: &RunRecord) -> Result<()> {
        let urls = run
            .remote_files
            .iter()
            .map(|file| file.url.clone())
            .collect::<Vec<_>>();
        let md5s = run
            .remote_files
            .iter()
            .map(|file| file.md5.clone())
            .collect::<Vec<_>>();
        let bytes = run
            .remote_files
            .iter()
            .map(|file| file.bytes)
            .collect::<Vec<_>>();
        conn.execute(
            "INSERT INTO runs (project_id, study_accession, source_study_accession, run_accession, sample_accession, experiment_accession, library_layout, instrument_platform, instrument_model, urls_json, md5_json, bytes_json, remote_files_json, metadata_json, state, last_error, updated_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16, ?17)
             ON CONFLICT(run_accession) DO UPDATE SET
               study_accession = excluded.study_accession,
               source_study_accession = excluded.source_study_accession,
               sample_accession = excluded.sample_accession,
               experiment_accession = excluded.experiment_accession,
               library_layout = excluded.library_layout,
               instrument_platform = excluded.instrument_platform,
               instrument_model = excluded.instrument_model,
               urls_json = excluded.urls_json,
               md5_json = excluded.md5_json,
               bytes_json = excluded.bytes_json,
               remote_files_json = excluded.remote_files_json,
               metadata_json = excluded.metadata_json,
               state = runs.state,
               last_error = runs.last_error,
               updated_at = excluded.updated_at",
            params![
                run.project_id,
                run.study_accession,
                run.source_study_accession,
                run.run_accession,
                run.sample_accession,
                run.experiment_accession,
                run.library_layout.as_str(),
                run.instrument_platform,
                run.instrument_model,
                serde_json::to_string(&urls)?,
                serde_json::to_string(&md5s)?,
                serde_json::to_string(&bytes)?,
                serde_json::to_string(&run.remote_files)?,
                serde_json::to_string(&run.metadata)?,
                run.state.as_str(),
                run.last_error,
                run.updated_at,
            ],
        )
        .with_context(|| format!("failed to upsert run {}", run.run_accession))?;
        Self::replace_run_metadata_tx(conn, &run.project_id, &run.run_accession, &run.metadata)?;
        Ok(())
    }

    fn replace_run_metadata_tx(
        conn: &Connection,
        project_id: &str,
        run_accession: &str,
        metadata: &std::collections::BTreeMap<String, String>,
    ) -> Result<()> {
        conn.execute(
            "DELETE FROM run_metadata WHERE project_id = ?1 AND run_accession = ?2",
            params![project_id, run_accession],
        )
        .with_context(|| format!("failed to clear normalized metadata for {run_accession}"))?;
        if metadata.is_empty() {
            return Ok(());
        }
        let mut stmt = conn.prepare(
            "INSERT INTO run_metadata (project_id, run_accession, field_name, field_value, updated_at)
             VALUES (?1, ?2, ?3, ?4, ?5)",
        )?;
        let updated_at = now_rfc3339();
        for (field_name, field_value) in metadata {
            stmt.execute(params![
                project_id,
                run_accession,
                field_name,
                field_value,
                updated_at,
            ])
            .with_context(|| {
                format!("failed to write normalized metadata {field_name} for {run_accession}")
            })?;
        }
        Ok(())
    }

    fn load_run_metadata_map(
        conn: &Connection,
        run_accession: &str,
    ) -> Result<std::collections::BTreeMap<String, String>> {
        let mut stmt = conn.prepare(
            "SELECT field_name, field_value
             FROM run_metadata
             WHERE run_accession = ?1
             ORDER BY field_name",
        )?;
        let rows = stmt.query_map(params![run_accession], |row| {
            Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
        })?;
        let mut metadata = std::collections::BTreeMap::new();
        for row in rows {
            let (field_name, field_value) = row?;
            metadata.insert(field_name, field_value);
        }
        Ok(metadata)
    }

    fn load_all_run_metadata_maps(
        conn: &Connection,
    ) -> Result<std::collections::BTreeMap<String, std::collections::BTreeMap<String, String>>>
    {
        let mut stmt = conn.prepare(
            "SELECT run_accession, field_name, field_value
             FROM run_metadata
             ORDER BY run_accession, field_name",
        )?;
        let rows = stmt.query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, String>(1)?,
                row.get::<_, String>(2)?,
            ))
        })?;
        let mut metadata_by_run = std::collections::BTreeMap::new();
        for row in rows {
            let (run_accession, field_name, field_value) = row?;
            metadata_by_run
                .entry(run_accession)
                .or_insert_with(std::collections::BTreeMap::new)
                .insert(field_name, field_value);
        }
        Ok(metadata_by_run)
    }

    pub(super) fn backfill_run_metadata(conn: &Connection) -> Result<()> {
        // Acquire write ownership before reading. A live SELECT cursor followed
        // by an autocommit write otherwise races other processes in WAL mode and
        // can fail immediately with SQLITE_BUSY_SNAPSHOT despite busy_timeout.
        let transaction =
            rusqlite::Transaction::new_unchecked(conn, rusqlite::TransactionBehavior::Immediate)?;
        let mut stmt = transaction.prepare(
            "SELECT project_id, run_accession, metadata_json
             FROM runs
             WHERE metadata_json <> '{}' AND NOT EXISTS (
               SELECT 1
               FROM run_metadata
               WHERE run_metadata.project_id = runs.project_id
                 AND run_metadata.run_accession = runs.run_accession
             )
             ORDER BY run_accession",
        )?;
        let rows = stmt.query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, String>(1)?,
                row.get::<_, String>(2)?,
            ))
        })?;
        let records = rows.collect::<rusqlite::Result<Vec<_>>>()?;
        drop(stmt);
        for (project_id, run_accession, metadata_json) in records {
            let metadata =
                from_json_text::<std::collections::BTreeMap<String, String>>(&metadata_json)
                    .with_context(|| {
                        format!(
                            "failed to backfill normalized metadata for legacy run {run_accession}"
                        )
                    })?;
            if !metadata.is_empty() {
                Self::replace_run_metadata_tx(
                    &transaction,
                    &project_id,
                    &run_accession,
                    &metadata,
                )?;
            }
        }
        transaction.commit()?;
        Ok(())
    }
}

#[derive(Debug)]
struct RawRun {
    project_id: String,
    study_accession: String,
    source_study_accession: Option<String>,
    run_accession: String,
    sample_accession: Option<String>,
    experiment_accession: Option<String>,
    library_layout: String,
    instrument_platform: Option<String>,
    instrument_model: Option<String>,
    remote_files_json: String,
    metadata_json: String,
    state: String,
    last_error: Option<String>,
    updated_at: String,
}

impl TryFrom<RawRun> for RunRecord {
    type Error = anyhow::Error;

    fn try_from(value: RawRun) -> Result<Self> {
        Ok(Self {
            project_id: value.project_id,
            study_accession: value.study_accession,
            source_study_accession: value.source_study_accession,
            run_accession: value.run_accession,
            sample_accession: value.sample_accession,
            experiment_accession: value.experiment_accession,
            library_layout: LibraryLayout::from_str(&value.library_layout)?,
            instrument_platform: value.instrument_platform,
            instrument_model: value.instrument_model,
            remote_files: from_json_text(&value.remote_files_json)?,
            metadata: from_json_text(&value.metadata_json)?,
            state: RunState::from_str(&value.state)?,
            last_error: value.last_error,
            updated_at: value.updated_at,
        })
    }
}

fn parse_raw_run(row: &Row<'_>) -> rusqlite::Result<RawRun> {
    Ok(RawRun {
        project_id: row.get(0)?,
        study_accession: row.get(1)?,
        source_study_accession: row.get(2)?,
        run_accession: row.get(3)?,
        sample_accession: row.get(4)?,
        experiment_accession: row.get(5)?,
        library_layout: row.get(6)?,
        instrument_platform: row.get(7)?,
        instrument_model: row.get(8)?,
        remote_files_json: row.get(9)?,
        metadata_json: row.get(10)?,
        state: row.get(11)?,
        last_error: row.get(12)?,
        updated_at: row.get(13)?,
    })
}
