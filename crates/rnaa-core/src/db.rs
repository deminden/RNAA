//! Database lifecycle, transactions and project-wide records; domain queries are delegated.
mod analysis;
mod artifacts;
#[cfg(test)]
mod migration_atomicity_tests;
mod qc;
mod run_progress;
pub use run_progress::LegacyRunAnalysisState;
mod runs;
#[cfg(test)]
mod tests;
use crate::config::ProjectConfig;
use crate::model::{EventRecord, ProjectRecord};
use crate::state::RunState;
use crate::util::now_rfc3339;
use anyhow::{Context, Result, anyhow, bail};
use rusqlite::{Connection, OptionalExtension, Row, params};
use serde::de::DeserializeOwned;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::{Arc, Mutex};
use std::time::Duration;

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
    (
        "0006_retained_inputs",
        include_str!("../migrations/0006_retained_inputs.sql"),
    ),
    ("0007_jobs", include_str!("../migrations/0007_jobs.sql")),
    (
        "0008_generations",
        include_str!("../migrations/0008_generations.sql"),
    ),
    (
        "0009_generation_pruning",
        include_str!("../migrations/0009_generation_pruning.sql"),
    ),
    (
        "0010_eligibility",
        include_str!("../migrations/0010_eligibility.sql"),
    ),
    (
        "0011_reference_selection",
        include_str!("../migrations/0011_reference_selection.sql"),
    ),
    (
        "0012_cohorts",
        include_str!("../migrations/0012_cohorts.sql"),
    ),
    (
        "0013_metadata_publications",
        include_str!("../migrations/0013_metadata_publications.sql"),
    ),
    (
        "0014_generation_lineage",
        include_str!("../migrations/0014_generation_lineage.sql"),
    ),
    (
        "0015_generation_request_lookup",
        include_str!("../migrations/0015_generation_request_lookup.sql"),
    ),
    (
        "0016_scientific_requests",
        include_str!("../migrations/0016_scientific_requests.sql"),
    ),
    (
        "0017_generation_retirement",
        include_str!("../migrations/0017_generation_retirement.sql"),
    ),
    (
        "0018_file_trash",
        include_str!("../migrations/0018_file_trash.sql"),
    ),
    (
        "0019_lineage_import",
        include_str!("../migrations/0019_lineage_import.sql"),
    ),
    (
        "0020_unpublished_trash",
        include_str!("../migrations/0020_unpublished_trash.sql"),
    ),
    (
        "0021_run_progress",
        include_str!("../migrations/0021_run_progress.sql"),
    ),
    (
        "0022_export_recovery",
        include_str!("../migrations/0022_export_recovery.sql"),
    ),
];

#[derive(Debug, Clone)]
pub struct Database {
    path: PathBuf,
    project_id: String,
    pub(crate) job: Option<crate::jobs::JobIdentity>,
    transaction: Option<Arc<Mutex<AtomicConnection>>>,
    pub(crate) generation_maintenance: bool,
}

#[derive(Debug)]
struct AtomicConnection {
    connection: Connection,
    active: bool,
}

struct AtomicScope(Arc<Mutex<AtomicConnection>>);
impl Drop for AtomicScope {
    fn drop(&mut self) {
        let mut session = self
            .0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if session.active {
            let _ = session.connection.execute_batch("ROLLBACK");
            session.active = false;
        }
    }
}

impl Database {
    pub fn create(root: &Path, config: &ProjectConfig) -> Result<Self> {
        let db_path = root.join("state.sqlite");
        let _maintenance_guard = crate::generations::pins::read_lock(&db_path)?;
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
            job: None,
            transaction: None,
            generation_maintenance: false,
        })
    }

    pub fn open(root: &Path) -> Result<Self> {
        let db_path = root.join("state.sqlite");
        if !db_path.exists() {
            bail!("RNAA project database not found at {}", db_path.display());
        }
        // Opening can migrate/backfill the catalog, so it participates in the
        // same maintenance exclusion as ordinary database mutations.
        let _maintenance_guard = crate::generations::pins::read_lock(&db_path)?;
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
            job: None,
            transaction: None,
            generation_maintenance: false,
        })
    }

    /// Reference work has an independent resource and heartbeat. Never detach an
    /// atomic publication context or bypass a cancelled parent's ownership.
    pub(crate) fn for_reference_job(&self) -> Result<Self> {
        if self.transaction.is_some() {
            anyhow::bail!("cannot start reference work inside a publication transaction");
        }
        self.check_job()?;
        let mut independent = self.clone();
        independent.job = None;
        Ok(independent)
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn project_id(&self) -> &str {
        &self.project_id
    }
    pub(crate) fn supports_migration(version: &str) -> bool {
        MIGRATIONS.iter().any(|(known, _)| *known == version)
    }

    pub fn project(&self) -> Result<ProjectRecord> {
        self.read(|conn| {
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
        })
    }

    pub fn set_project_config(&self, config: &ProjectConfig) -> Result<()> {
        self.write(|conn| {
            let config_json = serde_json::to_string_pretty(config)?;
            conn.execute(
                "UPDATE projects SET config_json = ?1 WHERE project_id = ?2",
                params![config_json, self.project_id],
            )
            .context("failed to update project config")?;
            Ok(())
        })
    }

    pub fn append_event(
        &self,
        stage: &str,
        run_accession: Option<&str>,
        message: &str,
        context: serde_json::Value,
    ) -> Result<()> {
        self.write(|conn| {
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
        })
    }

    pub fn list_events(&self, limit: usize) -> Result<Vec<EventRecord>> {
        self.read(|conn| {
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
        })
    }

    pub fn state_counts(&self) -> Result<Vec<(RunState, u64)>> {
        self.read(|conn| {
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
        })
    }

    pub fn last_errors(&self, limit: usize) -> Result<Vec<(String, String, String)>> {
        self.read(|conn| {
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
        })
    }

    pub(crate) fn connect(&self) -> Result<Connection> {
        Self::connect_path(&self.path)
    }

    pub(crate) fn read<T>(&self, operation: impl FnOnce(&Connection) -> Result<T>) -> Result<T> {
        if let Some(transaction) = &self.transaction {
            let session = transaction
                .lock()
                .map_err(|_| anyhow!("atomic database session poisoned"))?;
            if !session.active {
                bail!("atomic database session is closed");
            }
            return operation(&session.connection);
        }
        operation(&self.connect()?)
    }

    /// Group existing database operations in one transaction. Context handles
    /// retained after this callback are closed, and panics roll the transaction back.
    /// Compose a catalog operation with an existing publication transaction.
    pub(crate) fn atomic_or_current<T>(
        &self,
        operation: impl FnOnce(&Self) -> Result<T>,
    ) -> Result<T> {
        if self.transaction.is_some() {
            operation(self)
        } else {
            self.atomic(operation)
        }
    }

    pub fn atomic<T>(&self, operation: impl FnOnce(&Self) -> Result<T>) -> Result<T> {
        if self.transaction.is_some() {
            bail!("cannot nest an atomic database session");
        }
        let _generation_guard = if self.generation_maintenance {
            None
        } else {
            Some(crate::generations::pins::read_lock(&self.path)?)
        };
        let connection = self.connect()?;
        connection.execute_batch("BEGIN IMMEDIATE")?;
        let shared = Arc::new(Mutex::new(AtomicConnection {
            connection,
            active: true,
        }));
        let _scope = AtomicScope(Arc::clone(&shared));
        let mut context = self.clone();
        context.transaction = Some(Arc::clone(&shared));
        let result = operation(&context)?;
        let mut session = shared
            .lock()
            .map_err(|_| anyhow!("atomic database session poisoned"))?;
        session.connection.execute_batch("COMMIT")?;
        session.active = false;
        Ok(result)
    }

    // Validate ownership and mutate under the same write transaction. A stale
    // worker cannot pass a check and then overwrite a newer attempt's records.
    pub(crate) fn write<T>(&self, operation: impl FnOnce(&Connection) -> Result<T>) -> Result<T> {
        let _generation_guard = if self.generation_maintenance {
            None
        } else {
            Some(crate::generations::pins::read_lock(&self.path)?)
        };
        if self.transaction.is_some() {
            return self.read(|connection| {
                if let Some(job) = &self.job {
                    job.check(connection)?;
                }
                operation(connection)
            });
        }
        let mut connection = self.connect()?;
        let transaction =
            connection.transaction_with_behavior(rusqlite::TransactionBehavior::Immediate)?;
        if let Some(job) = &self.job {
            job.check(&transaction)?;
        }
        let result = operation(&transaction)?;
        transaction.commit()?;
        Ok(result)
    }

    fn connect_path(path: &Path) -> Result<Connection> {
        let conn =
            Connection::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        conn.busy_timeout(Duration::from_secs(30))?;
        conn.execute_batch(
            "PRAGMA journal_mode = WAL;
             PRAGMA foreign_keys = ON;
             PRAGMA synchronous = FULL;
             PRAGMA temp_store = MEMORY;",
        )?;
        Ok(conn)
    }

    fn migrate(conn: &Connection) -> Result<()> {
        // An upgrade must not rewrite old analysis labels while an old worker
        // still owns the cohort. Poll outside SQLite's writer transaction so a
        // competing initializer can finish and ordinary new-version opens proceed.
        let _analysis_upgrade = Self::run_progress_upgrade_guard(conn)?;
        let transaction =
            rusqlite::Transaction::new_unchecked(conn, rusqlite::TransactionBehavior::Immediate)?;
        transaction.execute(
            "CREATE TABLE IF NOT EXISTS schema_migrations (version TEXT PRIMARY KEY, applied_at TEXT NOT NULL)",
            [],
        )?;
        for (version, sql) in MIGRATIONS {
            let already_applied: Option<String> = transaction
                .query_row(
                    "SELECT version FROM schema_migrations WHERE version = ?1",
                    params![version],
                    |row| row.get(0),
                )
                .optional()?;
            if already_applied.is_none() {
                transaction
                    .execute_batch(sql)
                    .with_context(|| format!("failed to apply migration {version}"))?;
                transaction.execute(
                    "INSERT INTO schema_migrations (version, applied_at) VALUES (?1, ?2)",
                    params![version, now_rfc3339()],
                )?;
            }
        }
        transaction.commit()?;
        Ok(())
    }
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

fn from_json_text<T: DeserializeOwned>(text: &str) -> Result<T> {
    serde_json::from_str(text).map_err(|err| anyhow!("failed to parse JSON '{text}': {err}"))
}
