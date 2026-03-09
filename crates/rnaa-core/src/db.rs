use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::Duration;

use anyhow::{Context, Result, anyhow, bail};
use rusqlite::{Connection, OptionalExtension, Row, params};
use serde::de::DeserializeOwned;

use crate::config::{MetadataMergeStrategy, ProjectConfig};
use crate::model::{
    ArtifactRecord, ContrastSpec, EventRecord, InputRecord, InputType, LibraryLayout,
    MetadataOverrideRecord, PreprocessQcRecord, ProjectRecord, ProjectStageRecord, RunRecord,
    SharedBlobRecord,
};
use crate::state::{ProjectStageState, RunState};
use crate::util::now_rfc3339;

const MIGRATIONS: &[(&str, &str)] = &[
    ("0001_init", include_str!("../migrations/0001_init.sql")),
    (
        "0002_shared_store",
        include_str!("../migrations/0002_shared_store.sql"),
    ),
    (
        "0003_project_stages",
        include_str!("../migrations/0003_project_stages.sql"),
    ),
    (
        "0004_run_metadata",
        include_str!("../migrations/0004_run_metadata.sql"),
    ),
    (
        "0005_preprocess_qc",
        include_str!("../migrations/0005_preprocess_qc.sql"),
    ),
];

#[derive(Debug, Clone)]
pub struct Database {
    path: PathBuf,
    project_id: String,
}

impl Database {
    pub fn create(root: &Path, config: &ProjectConfig) -> Result<Self> {
        let db_path = root.join("state.sqlite");
        let conn = Self::connect_path(&db_path)?;
        Self::migrate(&conn)?;
        Self::backfill_run_metadata(&conn)?;

        let project_id = uuid::Uuid::new_v4().to_string();
        let project = ProjectRecord {
            project_id: project_id.clone(),
            root_dir: root.display().to_string(),
            created_at: now_rfc3339(),
            config_json: serde_json::to_string_pretty(config)?,
        };

        conn.execute(
            "INSERT INTO projects (project_id, root_dir, created_at, config_json) VALUES (?1, ?2, ?3, ?4)",
            params![
                project.project_id,
                project.root_dir,
                project.created_at,
                project.config_json
            ],
        )
        .context("failed to insert project record")?;

        Ok(Self {
            path: db_path,
            project_id,
        })
    }

    pub fn open(root: &Path) -> Result<Self> {
        let db_path = root.join("state.sqlite");
        if !db_path.exists() {
            bail!("RNAA project database not found at {}", db_path.display());
        }
        let conn = Self::connect_path(&db_path)?;
        Self::migrate(&conn)?;
        Self::backfill_run_metadata(&conn)?;
        let project_id: String = conn
            .query_row("SELECT project_id FROM projects LIMIT 1", [], |row| {
                row.get(0)
            })
            .context("failed to load project_id from database")?;
        Ok(Self {
            path: db_path,
            project_id,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn project_id(&self) -> &str {
        &self.project_id
    }

    pub fn project(&self) -> Result<ProjectRecord> {
        let conn = self.connect()?;
        conn.query_row(
            "SELECT project_id, root_dir, created_at, config_json FROM projects LIMIT 1",
            [],
            |row| {
                Ok(ProjectRecord {
                    project_id: row.get(0)?,
                    root_dir: row.get(1)?,
                    created_at: row.get(2)?,
                    config_json: row.get(3)?,
                })
            },
        )
        .context("failed to load project record")
    }

    pub fn set_project_config(&self, config: &ProjectConfig) -> Result<()> {
        let config_json = serde_json::to_string_pretty(config)?;
        let conn = self.connect()?;
        conn.execute(
            "UPDATE projects SET config_json = ?1 WHERE project_id = ?2",
            params![config_json, self.project_id],
        )
        .context("failed to update project config")?;
        Ok(())
    }

    pub fn add_input(&self, input_id: &str, input_type: InputType) -> Result<()> {
        let conn = self.connect()?;
        conn.execute(
            "INSERT OR IGNORE INTO inputs (project_id, input_id, input_type, added_at) VALUES (?1, ?2, ?3, ?4)",
            params![self.project_id, input_id, input_type.as_str(), now_rfc3339()],
        )
        .with_context(|| format!("failed to add input {input_id}"))?;
        Ok(())
    }

    pub fn list_inputs(&self) -> Result<Vec<InputRecord>> {
        let conn = self.connect()?;
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
    }

    pub fn add_metadata_override(
        &self,
        source_path: &Path,
        stored_path: &Path,
        merge_strategy: MetadataMergeStrategy,
    ) -> Result<()> {
        let conn = self.connect()?;
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
        Ok(())
    }

    pub fn list_metadata_overrides(&self) -> Result<Vec<MetadataOverrideRecord>> {
        let conn = self.connect()?;
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
    }

    pub fn upsert_runs(&self, runs: &[RunRecord]) -> Result<()> {
        let mut conn = self.connect()?;
        let tx = conn.transaction()?;
        for run in runs {
            Self::upsert_run_tx(&tx, run)?;
        }
        tx.commit()?;
        Ok(())
    }

    pub fn upsert_run(&self, run: &RunRecord) -> Result<()> {
        let conn = self.connect()?;
        Self::upsert_run_tx(&conn, run)
    }

    pub fn get_run(&self, run_accession: &str) -> Result<Option<RunRecord>> {
        let conn = self.connect()?;
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
        let metadata = Self::load_run_metadata_map(&conn, &run.run_accession)?;
        if !metadata.is_empty() {
            run.metadata = metadata;
        }
        Ok(Some(run))
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
        let conn = self.connect()?;
        let mut stmt = conn.prepare(&sql)?;
        let state_strings = states
            .iter()
            .map(|state| state.as_str().to_string())
            .collect::<Vec<_>>();
        let rows = stmt.query_map(
            rusqlite::params_from_iter(state_strings.iter()),
            parse_raw_run,
        )?;
        let metadata_by_run = Self::load_all_run_metadata_maps(&conn)?;
        let mut runs = Vec::new();
        for row in rows {
            let mut run: RunRecord = row?.try_into()?;
            if let Some(metadata) = metadata_by_run.get(&run.run_accession) {
                run.metadata = metadata.clone();
            }
            runs.push(run);
        }
        Ok(runs)
    }

    pub fn list_metadata_columns(&self) -> Result<Vec<String>> {
        let conn = self.connect()?;
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
    }

    pub fn set_run_state(
        &self,
        run_accession: &str,
        state: RunState,
        last_error: Option<&str>,
    ) -> Result<()> {
        let conn = self.connect()?;
        conn.execute(
            "UPDATE runs SET state = ?1, last_error = ?2, updated_at = ?3 WHERE run_accession = ?4",
            params![state.as_str(), last_error, now_rfc3339(), run_accession],
        )
        .with_context(|| format!("failed to update state for {run_accession}"))?;
        Ok(())
    }

    pub fn record_artifact(&self, artifact: &ArtifactRecord) -> Result<()> {
        let conn = self.connect()?;
        conn.execute(
            "INSERT OR REPLACE INTO artifacts (project_id, run_accession, kind, path, blob_id, shared_path, checksum_type, checksum, bytes, created_at) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)",
            params![
                artifact.project_id,
                artifact.run_accession,
                artifact.kind.as_str(),
                artifact.path,
                artifact.blob_id,
                artifact.shared_path,
                artifact.checksum_type,
                artifact.checksum,
                artifact.bytes as i64,
                artifact.created_at
            ],
        )
        .with_context(|| format!("failed to record artifact {}", artifact.path))?;
        if let Some(blob_id) = artifact.blob_id.as_deref() {
            conn.execute(
                "INSERT OR REPLACE INTO shared_blob_refs (blob_id, project_id, artifact_path, referenced_at) VALUES (?1, ?2, ?3, ?4)",
                params![blob_id, artifact.project_id, artifact.path, now_rfc3339()],
            )
            .with_context(|| {
                format!(
                    "failed to link artifact {} to shared blob {blob_id}",
                    artifact.path
                )
            })?;
        }
        Ok(())
    }

    pub fn record_shared_blob(&self, blob: &SharedBlobRecord) -> Result<()> {
        let conn = self.connect()?;
        conn.execute(
            "INSERT OR REPLACE INTO shared_blobs (blob_id, storage_path, checksum_type, checksum, bytes, created_at) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![
                blob.blob_id,
                blob.storage_path,
                blob.checksum_type,
                blob.checksum,
                blob.bytes as i64,
                blob.created_at
            ],
        )
        .with_context(|| format!("failed to record shared blob {}", blob.blob_id))?;
        Ok(())
    }

    pub fn get_shared_blob(&self, blob_id: &str) -> Result<Option<SharedBlobRecord>> {
        let conn = self.connect()?;
        let mut stmt = conn.prepare(
            "SELECT blob_id, storage_path, checksum_type, checksum, bytes, created_at FROM shared_blobs WHERE blob_id = ?1",
        )?;
        stmt.query_row(params![blob_id], parse_shared_blob)
            .optional()
            .context("failed to load shared blob")
    }

    pub fn list_artifacts(&self) -> Result<Vec<ArtifactRecord>> {
        self.list_artifacts_for_run(None)
    }

    pub fn list_artifacts_for_run(
        &self,
        run_accession: Option<&str>,
    ) -> Result<Vec<ArtifactRecord>> {
        let conn = self.connect()?;
        let sql = if run_accession.is_some() {
            "SELECT project_id, run_accession, kind, path, blob_id, shared_path, checksum_type, checksum, bytes, created_at FROM artifacts WHERE run_accession = ?1 ORDER BY kind, path"
        } else {
            "SELECT project_id, run_accession, kind, path, blob_id, shared_path, checksum_type, checksum, bytes, created_at FROM artifacts ORDER BY kind, path"
        };
        let mut stmt = conn.prepare(sql)?;
        let rows = if let Some(run_accession) = run_accession {
            stmt.query_map(params![run_accession], parse_artifact)?
        } else {
            stmt.query_map([], parse_artifact)?
        };
        let mut artifacts = Vec::new();
        for row in rows {
            artifacts.push(row?);
        }
        Ok(artifacts)
    }

    pub fn delete_artifact(&self, path: &str) -> Result<()> {
        let conn = self.connect()?;
        conn.execute("DELETE FROM artifacts WHERE path = ?1", params![path])
            .with_context(|| format!("failed to delete artifact record {path}"))?;
        Ok(())
    }

    pub fn record_preprocess_qc(&self, qc: &PreprocessQcRecord) -> Result<()> {
        let conn = self.connect()?;
        conn.execute(
            "INSERT INTO preprocess_qc (
               project_id, run_accession, report_path, mode, threads, total_reads_before,
               passed_reads, failed_reads, pass_rate, failed_rate, low_quality_reads,
               low_complexity_reads, too_many_n_reads, too_short_reads, duplicated_reads,
               duplication_rate, adapter_trimmed_reads, adapter_trimmed_bases, gate_passed,
               gate_reason, updated_at
             ) VALUES (
               ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16,
               ?17, ?18, ?19, ?20, ?21
             )
             ON CONFLICT(project_id, run_accession) DO UPDATE SET
               report_path = excluded.report_path,
               mode = excluded.mode,
               threads = excluded.threads,
               total_reads_before = excluded.total_reads_before,
               passed_reads = excluded.passed_reads,
               failed_reads = excluded.failed_reads,
               pass_rate = excluded.pass_rate,
               failed_rate = excluded.failed_rate,
               low_quality_reads = excluded.low_quality_reads,
               low_complexity_reads = excluded.low_complexity_reads,
               too_many_n_reads = excluded.too_many_n_reads,
               too_short_reads = excluded.too_short_reads,
               duplicated_reads = excluded.duplicated_reads,
               duplication_rate = excluded.duplication_rate,
               adapter_trimmed_reads = excluded.adapter_trimmed_reads,
               adapter_trimmed_bases = excluded.adapter_trimmed_bases,
               gate_passed = excluded.gate_passed,
               gate_reason = excluded.gate_reason,
               updated_at = excluded.updated_at",
            params![
                qc.project_id,
                qc.run_accession,
                qc.report_path,
                qc.mode,
                qc.threads as i64,
                qc.total_reads_before as i64,
                qc.passed_reads as i64,
                qc.failed_reads as i64,
                qc.pass_rate,
                qc.failed_rate,
                qc.low_quality_reads as i64,
                qc.low_complexity_reads as i64,
                qc.too_many_n_reads as i64,
                qc.too_short_reads as i64,
                qc.duplicated_reads as i64,
                qc.duplication_rate,
                qc.adapter_trimmed_reads as i64,
                qc.adapter_trimmed_bases as i64,
                if qc.gate_passed { 1_i64 } else { 0_i64 },
                qc.gate_reason,
                qc.updated_at,
            ],
        )
        .with_context(|| format!("failed to record preprocess QC for {}", qc.run_accession))?;
        Ok(())
    }

    pub fn get_preprocess_qc(&self, run_accession: &str) -> Result<Option<PreprocessQcRecord>> {
        let conn = self.connect()?;
        conn.query_row(
            "SELECT project_id, run_accession, report_path, mode, threads, total_reads_before,
                    passed_reads, failed_reads, pass_rate, failed_rate, low_quality_reads,
                    low_complexity_reads, too_many_n_reads, too_short_reads, duplicated_reads,
                    duplication_rate, adapter_trimmed_reads, adapter_trimmed_bases, gate_passed,
                    gate_reason, updated_at
             FROM preprocess_qc
             WHERE project_id = ?1 AND run_accession = ?2",
            params![self.project_id, run_accession],
            parse_preprocess_qc,
        )
        .optional()
        .context("failed to load preprocess qc")
    }

    pub fn list_preprocess_qc(&self) -> Result<Vec<PreprocessQcRecord>> {
        let conn = self.connect()?;
        let mut stmt = conn.prepare(
            "SELECT project_id, run_accession, report_path, mode, threads, total_reads_before,
                    passed_reads, failed_reads, pass_rate, failed_rate, low_quality_reads,
                    low_complexity_reads, too_many_n_reads, too_short_reads, duplicated_reads,
                    duplication_rate, adapter_trimmed_reads, adapter_trimmed_bases, gate_passed,
                    gate_reason, updated_at
             FROM preprocess_qc
             WHERE project_id = ?1
             ORDER BY run_accession",
        )?;
        let rows = stmt.query_map(params![self.project_id], parse_preprocess_qc)?;
        let mut records = Vec::new();
        for row in rows {
            records.push(row?);
        }
        Ok(records)
    }

    pub fn replace_contrasts(&self, contrasts: &[ContrastSpec]) -> Result<()> {
        let mut conn = self.connect()?;
        let tx = conn.transaction()?;
        tx.execute(
            "DELETE FROM contrasts WHERE project_id = ?1",
            params![self.project_id],
        )?;
        for contrast in contrasts {
            tx.execute(
                "INSERT INTO contrasts (project_id, name, factor, level_a, level_b) VALUES (?1, ?2, ?3, ?4, ?5)",
                params![
                    self.project_id,
                    contrast.name,
                    contrast.factor,
                    contrast.level_a,
                    contrast.level_b
                ],
            )?;
        }
        tx.commit()?;
        Ok(())
    }

    pub fn list_contrasts(&self) -> Result<Vec<ContrastSpec>> {
        let conn = self.connect()?;
        let mut stmt = conn.prepare(
            "SELECT name, factor, level_a, level_b FROM contrasts WHERE project_id = ?1 ORDER BY name",
        )?;
        let rows = stmt.query_map(params![self.project_id], |row| {
            Ok(ContrastSpec {
                name: row.get(0)?,
                factor: row.get(1)?,
                level_a: row.get(2)?,
                level_b: row.get(3)?,
            })
        })?;
        let mut contrasts = Vec::new();
        for row in rows {
            contrasts.push(row?);
        }
        Ok(contrasts)
    }

    pub fn append_event(
        &self,
        stage: &str,
        run_accession: Option<&str>,
        message: &str,
        context: serde_json::Value,
    ) -> Result<()> {
        let conn = self.connect()?;
        conn.execute(
            "INSERT INTO events (project_id, ts, run_accession, stage, message, context_json) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            params![
                self.project_id,
                now_rfc3339(),
                run_accession,
                stage,
                message,
                serde_json::to_string(&context)?
            ],
        )
        .with_context(|| format!("failed to append event for stage {stage}"))?;
        Ok(())
    }

    pub fn set_project_stage_state(
        &self,
        stage: &str,
        state: ProjectStageState,
        last_error: Option<&str>,
    ) -> Result<()> {
        let conn = self.connect()?;
        conn.execute(
            "INSERT INTO project_stages (project_id, stage, state, last_error, updated_at)
             VALUES (?1, ?2, ?3, ?4, ?5)
             ON CONFLICT(project_id, stage) DO UPDATE SET
               state = excluded.state,
               last_error = excluded.last_error,
               updated_at = excluded.updated_at",
            params![
                self.project_id,
                stage,
                state.as_str(),
                last_error,
                now_rfc3339()
            ],
        )
        .with_context(|| format!("failed to set project stage state for {stage}"))?;
        Ok(())
    }

    pub fn get_project_stage_state(&self, stage: &str) -> Result<Option<ProjectStageRecord>> {
        let conn = self.connect()?;
        conn.query_row(
            "SELECT project_id, stage, state, last_error, updated_at FROM project_stages WHERE project_id = ?1 AND stage = ?2",
            params![self.project_id, stage],
            parse_project_stage,
        )
        .optional()
        .context("failed to load project stage")
    }

    pub fn list_project_stage_states(&self) -> Result<Vec<ProjectStageRecord>> {
        let conn = self.connect()?;
        let mut stmt = conn.prepare(
            "SELECT project_id, stage, state, last_error, updated_at FROM project_stages WHERE project_id = ?1 ORDER BY stage",
        )?;
        let rows = stmt.query_map(params![self.project_id], parse_project_stage)?;
        let mut stages = Vec::new();
        for row in rows {
            stages.push(row?);
        }
        Ok(stages)
    }

    pub fn list_events(&self, limit: usize) -> Result<Vec<EventRecord>> {
        let conn = self.connect()?;
        let mut stmt = conn.prepare(
            "SELECT id, project_id, ts, run_accession, stage, message, context_json
             FROM events
             ORDER BY id DESC
             LIMIT ?1",
        )?;
        let rows = stmt.query_map(params![limit as i64], parse_event)?;
        let mut events = Vec::new();
        for row in rows {
            events.push(row?);
        }
        Ok(events)
    }

    pub fn state_counts(&self) -> Result<Vec<(RunState, u64)>> {
        let conn = self.connect()?;
        let mut stmt =
            conn.prepare("SELECT state, COUNT(*) FROM runs GROUP BY state ORDER BY state")?;
        let rows = stmt.query_map([], |row| {
            Ok((row.get::<_, String>(0)?, row.get::<_, i64>(1)?))
        })?;
        let mut counts = Vec::new();
        for row in rows {
            let (state, count) = row?;
            counts.push((RunState::from_str(&state)?, count as u64));
        }
        Ok(counts)
    }

    pub fn last_errors(&self, limit: usize) -> Result<Vec<(String, String, String)>> {
        let conn = self.connect()?;
        let mut stmt = conn.prepare(
            "SELECT run_accession, state, last_error FROM runs WHERE last_error IS NOT NULL ORDER BY updated_at DESC LIMIT ?1",
        )?;
        let rows = stmt.query_map(params![limit as i64], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, String>(1)?,
                row.get::<_, String>(2)?,
            ))
        })?;
        let mut errors = Vec::new();
        for row in rows {
            errors.push(row?);
        }
        Ok(errors)
    }

    pub fn disk_usage_bytes(&self) -> Result<u64> {
        let conn = self.connect()?;
        let total: Option<i64> =
            conn.query_row("SELECT SUM(bytes) FROM artifacts", [], |row| row.get(0))?;
        Ok(total.unwrap_or_default() as u64)
    }

    fn list_runs_filtered(&self, _filter: Option<&str>) -> Result<Vec<RunRecord>> {
        let conn = self.connect()?;
        let mut stmt = conn.prepare(
            "SELECT project_id, study_accession, source_study_accession, run_accession, sample_accession, experiment_accession, library_layout, instrument_platform, instrument_model, remote_files_json, metadata_json, state, last_error, updated_at FROM runs ORDER BY run_accession",
        )?;
        let rows = stmt.query_map([], parse_raw_run)?;
        let metadata_by_run = Self::load_all_run_metadata_maps(&conn)?;
        let mut runs = Vec::new();
        for row in rows {
            let mut run: RunRecord = row?.try_into()?;
            if let Some(metadata) = metadata_by_run.get(&run.run_accession) {
                run.metadata = metadata.clone();
            }
            runs.push(run);
        }
        Ok(runs)
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

    fn connect(&self) -> Result<Connection> {
        Self::connect_path(&self.path)
    }

    fn connect_path(path: &Path) -> Result<Connection> {
        let conn =
            Connection::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        conn.busy_timeout(Duration::from_secs(30))?;
        conn.execute_batch(
            "PRAGMA journal_mode = WAL;
             PRAGMA foreign_keys = ON;
             PRAGMA synchronous = NORMAL;
             PRAGMA temp_store = MEMORY;",
        )?;
        Ok(conn)
    }

    fn migrate(conn: &Connection) -> Result<()> {
        conn.execute(
            "CREATE TABLE IF NOT EXISTS schema_migrations (version TEXT PRIMARY KEY, applied_at TEXT NOT NULL)",
            [],
        )?;
        for (version, sql) in MIGRATIONS {
            let already_applied: Option<String> = conn
                .query_row(
                    "SELECT version FROM schema_migrations WHERE version = ?1",
                    params![version],
                    |row| row.get(0),
                )
                .optional()?;
            if already_applied.is_none() {
                conn.execute_batch(sql)
                    .with_context(|| format!("failed to apply migration {version}"))?;
                conn.execute(
                    "INSERT INTO schema_migrations (version, applied_at) VALUES (?1, ?2)",
                    params![version, now_rfc3339()],
                )?;
            }
        }
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

    fn backfill_run_metadata(conn: &Connection) -> Result<()> {
        let mut stmt = conn.prepare(
            "SELECT project_id, run_accession, metadata_json
             FROM runs
             WHERE NOT EXISTS (
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
        for row in rows {
            let (project_id, run_accession, metadata_json) = row?;
            let metadata =
                from_json_text::<std::collections::BTreeMap<String, String>>(&metadata_json)
                    .with_context(|| {
                        format!(
                            "failed to backfill normalized metadata for legacy run {run_accession}"
                        )
                    })?;
            Self::replace_run_metadata_tx(conn, &project_id, &run_accession, &metadata)?;
        }
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

fn parse_artifact(row: &Row<'_>) -> rusqlite::Result<ArtifactRecord> {
    let kind_text: String = row.get(2)?;
    let bytes: i64 = row.get(8)?;
    let kind = crate::model::ArtifactKind::from_str(&kind_text).map_err(|err| {
        rusqlite::Error::FromSqlConversionFailure(
            2,
            rusqlite::types::Type::Text,
            Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                err.to_string(),
            )),
        )
    })?;
    Ok(ArtifactRecord {
        project_id: row.get(0)?,
        run_accession: row.get(1)?,
        kind,
        path: row.get(3)?,
        blob_id: row.get(4)?,
        shared_path: row.get(5)?,
        checksum_type: row.get(6)?,
        checksum: row.get(7)?,
        bytes: bytes as u64,
        created_at: row.get(9)?,
    })
}

fn parse_shared_blob(row: &Row<'_>) -> rusqlite::Result<SharedBlobRecord> {
    let bytes: i64 = row.get(4)?;
    Ok(SharedBlobRecord {
        blob_id: row.get(0)?,
        storage_path: row.get(1)?,
        checksum_type: row.get(2)?,
        checksum: row.get(3)?,
        bytes: bytes as u64,
        created_at: row.get(5)?,
    })
}

fn parse_event(row: &Row<'_>) -> rusqlite::Result<EventRecord> {
    let context_text: String = row.get(6)?;
    let context = serde_json::from_str(&context_text).unwrap_or(serde_json::Value::Null);
    Ok(EventRecord {
        id: row.get(0)?,
        project_id: row.get(1)?,
        ts: row.get(2)?,
        run_accession: row.get(3)?,
        stage: row.get(4)?,
        message: row.get(5)?,
        context,
    })
}

fn parse_project_stage(row: &Row<'_>) -> rusqlite::Result<ProjectStageRecord> {
    let state_text: String = row.get(2)?;
    let state = ProjectStageState::from_str(&state_text).map_err(|err| {
        rusqlite::Error::FromSqlConversionFailure(
            2,
            rusqlite::types::Type::Text,
            Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                err.to_string(),
            )),
        )
    })?;
    Ok(ProjectStageRecord {
        project_id: row.get(0)?,
        stage: row.get(1)?,
        state,
        last_error: row.get(3)?,
        updated_at: row.get(4)?,
    })
}

fn parse_preprocess_qc(row: &Row<'_>) -> rusqlite::Result<PreprocessQcRecord> {
    Ok(PreprocessQcRecord {
        project_id: row.get(0)?,
        run_accession: row.get(1)?,
        report_path: row.get(2)?,
        mode: row.get(3)?,
        threads: row.get::<_, i64>(4)? as u64,
        total_reads_before: row.get::<_, i64>(5)? as u64,
        passed_reads: row.get::<_, i64>(6)? as u64,
        failed_reads: row.get::<_, i64>(7)? as u64,
        pass_rate: row.get(8)?,
        failed_rate: row.get(9)?,
        low_quality_reads: row.get::<_, i64>(10)? as u64,
        low_complexity_reads: row.get::<_, i64>(11)? as u64,
        too_many_n_reads: row.get::<_, i64>(12)? as u64,
        too_short_reads: row.get::<_, i64>(13)? as u64,
        duplicated_reads: row.get::<_, i64>(14)? as u64,
        duplication_rate: row.get(15)?,
        adapter_trimmed_reads: row.get::<_, i64>(16)? as u64,
        adapter_trimmed_bases: row.get::<_, i64>(17)? as u64,
        gate_passed: row.get::<_, i64>(18)? != 0,
        gate_reason: row.get(19)?,
        updated_at: row.get(20)?,
    })
}

fn from_json_text<T: DeserializeOwned>(text: &str) -> Result<T> {
    serde_json::from_str(text).map_err(|err| anyhow!("failed to parse JSON '{text}': {err}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn run_state_transitions_are_persisted() {
        let temp = TempDir::new().unwrap();
        let config = ProjectConfig::default();
        let db = Database::create(temp.path(), &config).unwrap();
        let run = RunRecord {
            project_id: db.project_id().to_string(),
            study_accession: "SRP1".to_string(),
            source_study_accession: Some("PRJNA1".to_string()),
            run_accession: "SRR1".to_string(),
            sample_accession: Some("SRS1".to_string()),
            experiment_accession: Some("SRX1".to_string()),
            library_layout: LibraryLayout::Paired,
            instrument_platform: Some("ILLUMINA".to_string()),
            instrument_model: Some("HiSeq".to_string()),
            remote_files: vec![],
            metadata: BTreeMap::new(),
            state: RunState::Resolved,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        db.upsert_run(&run).unwrap();

        db.set_run_state("SRR1", RunState::Downloading, None)
            .unwrap();
        db.set_run_state("SRR1", RunState::Verified, None).unwrap();

        let persisted = db.get_run("SRR1").unwrap().unwrap();
        assert_eq!(persisted.state, RunState::Verified);
    }

    #[test]
    fn upsert_run_preserves_existing_state_for_resolve_refresh() {
        let temp = TempDir::new().unwrap();
        let config = ProjectConfig::default();
        let db = Database::create(temp.path(), &config).unwrap();
        let mut run = RunRecord {
            project_id: db.project_id().to_string(),
            study_accession: "SRP1".to_string(),
            source_study_accession: Some("PRJNA1".to_string()),
            run_accession: "SRR1".to_string(),
            sample_accession: Some("SRS1".to_string()),
            experiment_accession: Some("SRX1".to_string()),
            library_layout: LibraryLayout::Paired,
            instrument_platform: Some("ILLUMINA".to_string()),
            instrument_model: Some("HiSeq".to_string()),
            remote_files: vec![],
            metadata: BTreeMap::from([("condition".to_string(), "A".to_string())]),
            state: RunState::Resolved,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        db.upsert_run(&run).unwrap();
        db.set_run_state("SRR1", RunState::QuantDone, None).unwrap();

        run.metadata
            .insert("condition".to_string(), "B".to_string());
        run.state = RunState::Resolved;
        run.updated_at = now_rfc3339();
        db.upsert_run(&run).unwrap();

        let persisted = db.get_run("SRR1").unwrap().unwrap();
        assert_eq!(persisted.state, RunState::QuantDone);
        assert_eq!(
            persisted.metadata.get("condition").map(String::as_str),
            Some("B")
        );
    }

    #[test]
    fn normalized_run_metadata_is_persisted_and_queryable() {
        let temp = TempDir::new().unwrap();
        let config = ProjectConfig::default();
        let db = Database::create(temp.path(), &config).unwrap();
        let run = RunRecord {
            project_id: db.project_id().to_string(),
            study_accession: "SRP1".to_string(),
            source_study_accession: Some("PRJNA1".to_string()),
            run_accession: "SRR1".to_string(),
            sample_accession: Some("SRS1".to_string()),
            experiment_accession: Some("SRX1".to_string()),
            library_layout: LibraryLayout::Paired,
            instrument_platform: Some("ILLUMINA".to_string()),
            instrument_model: Some("HiSeq".to_string()),
            remote_files: vec![],
            metadata: BTreeMap::from([
                ("condition".to_string(), "A".to_string()),
                ("batch".to_string(), "B1".to_string()),
                (
                    "ena_library_source".to_string(),
                    "TRANSCRIPTOMIC".to_string(),
                ),
            ]),
            state: RunState::Resolved,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        db.upsert_run(&run).unwrap();

        let persisted = db.get_run("SRR1").unwrap().unwrap();
        assert_eq!(
            persisted.metadata.get("condition").map(String::as_str),
            Some("A")
        );
        assert_eq!(
            persisted.metadata.get("batch").map(String::as_str),
            Some("B1")
        );
        assert_eq!(
            db.list_metadata_columns().unwrap(),
            vec![
                "batch".to_string(),
                "condition".to_string(),
                "ena_library_source".to_string()
            ]
        );
    }

    #[test]
    fn legacy_metadata_json_is_backfilled_into_normalized_table_on_open() {
        let temp = TempDir::new().unwrap();
        let config = ProjectConfig::default();
        let db = Database::create(temp.path(), &config).unwrap();
        let run = RunRecord {
            project_id: db.project_id().to_string(),
            study_accession: "SRP1".to_string(),
            source_study_accession: Some("PRJNA1".to_string()),
            run_accession: "SRR1".to_string(),
            sample_accession: Some("SRS1".to_string()),
            experiment_accession: Some("SRX1".to_string()),
            library_layout: LibraryLayout::Paired,
            instrument_platform: Some("ILLUMINA".to_string()),
            instrument_model: Some("HiSeq".to_string()),
            remote_files: vec![],
            metadata: BTreeMap::from([
                ("condition".to_string(), "A".to_string()),
                ("donor".to_string(), "D1".to_string()),
            ]),
            state: RunState::Resolved,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        db.upsert_run(&run).unwrap();

        let conn = Connection::open(db.path()).unwrap();
        conn.execute(
            "DELETE FROM run_metadata WHERE project_id = ?1 AND run_accession = ?2",
            params![db.project_id(), "SRR1"],
        )
        .unwrap();
        let deleted_count: i64 = conn
            .query_row(
                "SELECT COUNT(*) FROM run_metadata WHERE run_accession = ?1",
                params!["SRR1"],
                |row| row.get(0),
            )
            .unwrap();
        assert_eq!(deleted_count, 0);

        let reopened = Database::open(temp.path()).unwrap();
        let persisted = reopened.get_run("SRR1").unwrap().unwrap();
        assert_eq!(
            persisted.metadata.get("condition").map(String::as_str),
            Some("A")
        );
        assert_eq!(
            persisted.metadata.get("donor").map(String::as_str),
            Some("D1")
        );
        assert_eq!(
            reopened.list_metadata_columns().unwrap(),
            vec!["condition".to_string(), "donor".to_string()]
        );
    }

    #[test]
    fn artifact_blob_linkage_is_persisted() {
        let temp = TempDir::new().unwrap();
        let config = ProjectConfig::default();
        let db = Database::create(temp.path(), &config).unwrap();
        let artifact_path = temp.path().join("quant").join("SRR1").join("abundance.h5");
        fs::create_dir_all(artifact_path.parent().unwrap()).unwrap();
        fs::write(&artifact_path, b"test").unwrap();

        let blob = SharedBlobRecord {
            blob_id: "sha256:abc123".to_string(),
            storage_path: "/tmp/shared/blobs/sha256/ab/c1/abc123".to_string(),
            checksum_type: "sha256".to_string(),
            checksum: "abc123".to_string(),
            bytes: 4,
            created_at: now_rfc3339(),
        };
        db.record_shared_blob(&blob).unwrap();
        db.record_artifact(&ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: Some("SRR1".to_string()),
            kind: crate::model::ArtifactKind::QuantAbundanceH5,
            path: artifact_path.display().to_string(),
            blob_id: Some(blob.blob_id.clone()),
            shared_path: Some(blob.storage_path.clone()),
            checksum_type: "sha256".to_string(),
            checksum: "abc123".to_string(),
            bytes: 4,
            created_at: now_rfc3339(),
        })
        .unwrap();

        let artifacts = db.list_artifacts_for_run(Some("SRR1")).unwrap();
        assert_eq!(artifacts.len(), 1);
        assert_eq!(artifacts[0].blob_id.as_deref(), Some("sha256:abc123"));
        assert_eq!(
            artifacts[0].shared_path.as_deref(),
            Some("/tmp/shared/blobs/sha256/ab/c1/abc123")
        );

        let loaded_blob = db.get_shared_blob("sha256:abc123").unwrap().unwrap();
        assert_eq!(loaded_blob.checksum, "abc123");
    }

    #[test]
    fn project_stage_states_are_persisted() {
        let temp = TempDir::new().unwrap();
        let config = ProjectConfig::default();
        let db = Database::create(temp.path(), &config).unwrap();

        db.set_project_stage_state("normalize", ProjectStageState::Running, None)
            .unwrap();
        db.set_project_stage_state("normalize", ProjectStageState::Failed, Some("interrupted"))
            .unwrap();

        let stage = db.get_project_stage_state("normalize").unwrap().unwrap();
        assert_eq!(stage.state, ProjectStageState::Failed);
        assert_eq!(stage.last_error.as_deref(), Some("interrupted"));
    }
}
