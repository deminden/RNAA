//! Durable local worker ownership. Database fencing and file writer exclusion are
//! complementary: expiry cannot permit a second writer while a paused process is
//! still inside a scientific engine that does not expose cancellation callbacks.
mod database;
mod locking;

pub(crate) use locking::analysis_upgrade_lock;
use locking::{resource_key, writer_lock};

use crate::{Database, RunState, util::now_rfc3339};
use anyhow::{Context, Result, bail};
use rusqlite::{Connection, OptionalExtension, params};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::fs::{File, OpenOptions, TryLockError};
use std::sync::mpsc::{self, RecvTimeoutError, Sender};
use std::thread::{self, JoinHandle};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

const LEASE: Duration = Duration::from_secs(45);
const HEARTBEAT: Duration = Duration::from_secs(5);

#[derive(Debug, Clone)]
pub(crate) struct JobIdentity {
    pub attempt_id: String,
    owner: String,
    token: i64,
}

impl JobIdentity {
    pub(crate) fn check(&self, connection: &Connection) -> Result<()> {
        let valid: bool = connection.query_row(
            "SELECT EXISTS(SELECT 1 FROM job_attempts a JOIN job_resources r
             ON r.project_id = a.project_id AND r.resource = a.resource
             WHERE a.attempt_id = ?1 AND a.owner = ?2 AND a.fencing_token = ?3
             AND r.fencing_token = a.fencing_token AND a.state = 'RUNNING'
             AND a.cancel_requested = 0 AND a.lease_expires_ms > ?4)",
            params![self.attempt_id, self.owner, self.token, now_ms()?],
            |row| row.get(0),
        )?;
        if !valid {
            bail!(
                "job {} lost ownership, expired, or was cancelled",
                self.attempt_id
            );
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct JobAttempt {
    pub attempt_id: String,
    pub resource: String,
    pub stage: String,
    pub run_accession: Option<String>,
    pub request_sha256: Option<String>,
    pub owner: String,
    pub fencing_token: i64,
    pub state: String,
    pub cancel_requested: bool,
    pub started_at: String,
    pub heartbeat_ms: i64,
    pub lease_expires_ms: i64,
    pub finished_at: Option<String>,
    pub last_error: Option<String>,
}

#[derive(Debug, Clone, Copy)]
pub enum JobOutcome {
    Succeeded,
    Failed,
    Interrupted,
}

pub struct JobLease {
    _generation_pin: std::sync::Arc<crate::generations::GenerationStorePin>,
    db: Database,
    identity: JobIdentity,
    // Never remove lock files: unlinking would allow two independent lock inodes.
    _writer_lock: File,
    stop: Option<Sender<()>>,
    heartbeat: Option<JoinHandle<()>>,
    finished: bool,
}

impl JobLease {
    /// None means another worker owns the resource, or the queued run is no longer
    /// in a ready state. Readiness, the claim and the RUNNING state are atomic.
    pub fn try_acquire(
        db: &Database,
        stage: &str,
        run: Option<&str>,
        ready: &[RunState],
    ) -> Result<Option<Self>> {
        Self::acquire_with_timing(db, stage, run, ready, LEASE, HEARTBEAT)
    }

    pub fn try_acquire_reference(db: &Database, reference: &str) -> Result<Option<Self>> {
        if std::path::Path::new(reference).components().count() != 1
            || !matches!(
                std::path::Path::new(reference).components().next(),
                Some(std::path::Component::Normal(_))
            )
        {
            bail!("reference subject must be a single path component");
        }
        let independent = db.for_reference_job()?;
        Self::acquire_resource(
            &independent,
            "refs",
            None,
            &[],
            LEASE,
            HEARTBEAT,
            &format!("reference:{reference}"),
        )
    }

    fn acquire_with_timing(
        db: &Database,
        stage: &str,
        run: Option<&str>,
        ready: &[RunState],
        ttl: Duration,
        interval: Duration,
    ) -> Result<Option<Self>> {
        Self::acquire_resource(db, stage, run, ready, ttl, interval, &resource_key(run))
    }

    fn acquire_resource(
        db: &Database,
        stage: &str,
        run: Option<&str>,
        ready: &[RunState],
        ttl: Duration,
        interval: Duration,
        resource: &str,
    ) -> Result<Option<Self>> {
        if db.job.is_some() {
            bail!("cannot nest a worker lease");
        }
        let generation_pin = crate::generations::GenerationStorePin::acquire(db.path(), stage)?;
        let running = running_state(stage)?;
        if matches!(
            stage,
            "normalize" | "deseq2" | "corr" | "refs" | "resolve" | "metadata" | "export" | "survey"
        ) != run.is_none()
        {
            bail!("analysis jobs require a project subject; ingest jobs require a run");
        }
        if run.is_some() && ready.is_empty() {
            bail!("run jobs require explicit ready states");
        }
        if stage == "refs" && !resource.starts_with("reference:") {
            bail!("reference jobs require a reference resource");
        }
        let Some(lock) = writer_lock(db, resource)? else {
            return Ok(None);
        };
        let identity = db.write(|connection| {
            if !expire_resource(connection, db, resource)? {
                return Ok(None);
            }
            if let Some(run) = run {
                let state: Option<String> = connection.query_row(
                    "SELECT state FROM runs WHERE project_id = ?1 AND run_accession = ?2",
                    params![db.project_id(), run], |row| row.get(0),
                ).optional()?;
                if !ready.iter().any(|ready| Some(ready.as_str()) == state.as_deref()) {
                    return Ok(None);
                }
            }
            connection.execute(
                "INSERT INTO job_resources(project_id, resource, fencing_token) VALUES (?1, ?2, 1)
                 ON CONFLICT(project_id, resource) DO UPDATE SET fencing_token = fencing_token + 1",
                params![db.project_id(), resource],
            )?;
            let token = connection.query_row(
                "SELECT fencing_token FROM job_resources WHERE project_id = ?1 AND resource = ?2",
                params![db.project_id(), resource], |row| row.get(0),
            )?;
            let identity = JobIdentity {
                attempt_id: uuid::Uuid::new_v4().to_string(),
                owner: format!("{}:{}", std::process::id(), uuid::Uuid::new_v4()),
                token,
            };
            let now = now_ms()?;
            connection.execute(
                "INSERT INTO job_attempts(attempt_id, project_id, resource, stage, run_accession,
                 owner, fencing_token, state, started_at, heartbeat_ms, lease_expires_ms)
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, 'RUNNING', ?8, ?9, ?10)",
                params![identity.attempt_id, db.project_id(), resource, stage, run,
                    identity.owner, identity.token, now_rfc3339(), now, now + duration_ms(ttl)?],
            )?;
            if let (Some(run), Some(running)) = (run, running) {
                connection.execute(
                    "UPDATE runs SET state = ?1, last_error = NULL, updated_at = ?2 WHERE run_accession = ?3",
                    params![running.as_str(), now_rfc3339(), run],
                )?;
            }
            job_event(connection, db, stage, run, "job claimed", &identity.attempt_id)?;
            Ok(Some(identity))
        })?;
        let Some(identity) = identity else {
            return Ok(None);
        };
        let (sender, receiver) = mpsc::channel();
        let renewal_db = db.clone();
        let renewal_identity = identity.clone();
        let heartbeat = thread::spawn(move || {
            while matches!(
                receiver.recv_timeout(interval),
                Err(RecvTimeoutError::Timeout)
            ) {
                let renewed = renewal_db.write(|connection| {
                    renewal_identity.check(connection)?;
                    let now = now_ms()?;
                    connection.execute(
                        "UPDATE job_attempts SET heartbeat_ms = ?1, lease_expires_ms = ?2 WHERE attempt_id = ?3",
                        params![now, now + duration_ms(ttl)?, renewal_identity.attempt_id],
                    )?;
                    Ok(())
                });
                if let Err(error) = renewed {
                    let terminal =
                        renewal_db
                            .read(|connection| {
                                Ok(connection.query_row(
                        "SELECT state <> 'RUNNING' FROM job_attempts WHERE attempt_id = ?1",
                        params![renewal_identity.attempt_id], |row| row.get::<_, bool>(0))?)
                            })
                            .unwrap_or(false);
                    if !terminal {
                        tracing::warn!(attempt = %renewal_identity.attempt_id, "job heartbeat stopped: {error:#}");
                    }
                    break;
                }
            }
        });
        Ok(Some(Self {
            _generation_pin: generation_pin,
            db: db.clone(),
            identity,
            _writer_lock: lock,
            stop: Some(sender),
            heartbeat: Some(heartbeat),
            finished: false,
        }))
    }

    /// Mutations through this handle check the fencing token inside their write
    /// transaction. Copies retained after completion lose their write authority.
    pub fn database(&self) -> Database {
        let mut db = self.db.clone();
        db.job = Some(self.identity.clone());
        db
    }

    pub fn attempt_id(&self) -> &str {
        &self.identity.attempt_id
    }

    pub fn finish(&mut self, outcome: JobOutcome, message: Option<&str>) -> Result<()> {
        self.stop_heartbeat();
        let already_committed = self.db.read(|connection| {
            Ok(connection.query_row("SELECT EXISTS(SELECT 1 FROM job_attempts a JOIN job_resources r
                ON r.project_id = a.project_id AND r.resource = a.resource WHERE a.attempt_id = ?1
                AND a.owner = ?2 AND a.fencing_token = ?3 AND r.fencing_token = a.fencing_token AND a.state = ?4)",
                params![self.identity.attempt_id, self.identity.owner, self.identity.token, match outcome {
                    JobOutcome::Succeeded => "SUCCEEDED", JobOutcome::Failed => "FAILED", JobOutcome::Interrupted => "INTERRUPTED",
                }], |row| row.get::<_, bool>(0))?)
        })?;
        if !already_committed {
            self.terminate(outcome, message)?;
        }
        self.finished = true;
        Ok(())
    }

    /// End ownership and release the physical writer lock before starting a
    /// different operation on this resource. Completed fenced handles stay inert.
    pub fn finish_and_release(
        mut self,
        outcome: JobOutcome,
        message: Option<&str>,
    ) -> Result<Database> {
        self.finish(outcome, message)?;
        Ok(self.db.clone())
    }

    fn stop_heartbeat(&mut self) {
        if let Some(sender) = self.stop.take() {
            let _ = sender.send(());
        }
        if let Some(thread) = self.heartbeat.take() {
            let _ = thread.join();
        }
    }

    fn terminate(&self, outcome: JobOutcome, message: Option<&str>) -> Result<()> {
        self.db
            .write(|connection| self.identity.finish(&self.db, connection, outcome, message))
    }
}

impl JobIdentity {
    fn finish(
        &self,
        db: &Database,
        connection: &Connection,
        outcome: JobOutcome,
        message: Option<&str>,
    ) -> Result<()> {
        if matches!(outcome, JobOutcome::Succeeded) {
            self.check(connection)?;
        }
        let job = connection
            .query_row(
                "SELECT a.stage, a.run_accession, a.cancel_requested FROM job_attempts a
                 JOIN job_resources r ON r.project_id = a.project_id AND r.resource = a.resource
                 WHERE a.attempt_id = ?1 AND a.owner = ?2 AND a.fencing_token = ?3
                 AND r.fencing_token = a.fencing_token AND a.state = 'RUNNING'",
                params![self.attempt_id, self.owner, self.token],
                |row| {
                    Ok((
                        row.get::<_, String>(0)?,
                        row.get::<_, Option<String>>(1)?,
                        row.get::<_, bool>(2)?,
                    ))
                },
            )
            .optional()?
            .context("attempt no longer owns the job")?;
        let status = if job.2 {
            "CANCELLED"
        } else {
            match outcome {
                JobOutcome::Succeeded => "SUCCEEDED",
                JobOutcome::Failed => "FAILED",
                JobOutcome::Interrupted => "INTERRUPTED",
            }
        };
        let detail = if job.2 {
            Some("cancellation requested".to_string())
        } else if let Some(message) = message {
            Some(message.to_string())
        } else if status == "FAILED" {
            let error: Option<String> = if let Some(run) = &job.1 {
                connection
                    .query_row(
                        "SELECT last_error FROM runs WHERE run_accession = ?1",
                        params![run],
                        |row| row.get(0),
                    )
                    .optional()?
                    .flatten()
            } else {
                connection.query_row("SELECT last_error FROM project_stages WHERE project_id = ?1 AND stage = ?2", params![db.project_id(), job.0], |row| row.get(0)).optional()?.flatten()
            };
            Some(error.unwrap_or_else(|| "stage reported failure".into()))
        } else {
            None
        };
        connection.execute(
                "UPDATE job_attempts SET state = ?1, finished_at = ?2, last_error = ?3 WHERE attempt_id = ?4",
                params![status, now_rfc3339(), detail, self.attempt_id],
            )?;
        if status != "SUCCEEDED" {
            reset_running(
                connection,
                &job.0,
                job.1.as_deref(),
                detail.as_deref().unwrap_or(status),
            )?;
        }
        job_event(
            connection,
            db,
            &job.0,
            job.1.as_deref(),
            status,
            &self.attempt_id,
        )
    }
}

impl Drop for JobLease {
    fn drop(&mut self) {
        self.stop_heartbeat();
        if !self.finished
            && let Err(error) = self.terminate(
                JobOutcome::Interrupted,
                Some("worker exited before completing the attempt"),
            )
        {
            tracing::warn!(attempt = %self.identity.attempt_id, "could not close interrupted job: {error:#}");
        }
    }
}

/// False means an unexpired lease still exists. The caller already holds the
/// physical writer lock, so expired attempts can safely be fenced and recovered.
fn expire_resource(connection: &Connection, db: &Database, resource: &str) -> Result<bool> {
    let job = connection.query_row(
        "SELECT attempt_id, stage, run_accession, lease_expires_ms, cancel_requested FROM job_attempts
         WHERE project_id = ?1 AND resource = ?2 AND state = 'RUNNING'",
        params![db.project_id(), resource],
        |row| Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?, row.get::<_, Option<String>>(2)?, row.get::<_, i64>(3)?, row.get::<_, bool>(4)?)),
    ).optional()?;
    if let Some((attempt, stage, run, expiry, cancelled)) = job {
        if expiry > now_ms()? {
            return Ok(false);
        }
        connection.execute(
            "UPDATE job_attempts SET state = ?1, finished_at = ?2, last_error = 'worker lease expired' WHERE attempt_id = ?3",
            params![if cancelled { "CANCELLED" } else { "INTERRUPTED" }, now_rfc3339(), attempt],
        )?;
        reset_running(connection, &stage, run.as_deref(), "worker lease expired")?;
        job_event(
            connection,
            db,
            &stage,
            run.as_deref(),
            "expired worker recovered",
            &attempt,
        )?;
    }
    Ok(true)
}

fn running_state(stage: &str) -> Result<Option<RunState>> {
    Ok(Some(match stage {
        "retention" | "refs" | "resolve" | "metadata" | "normalize" | "deseq2" | "corr"
        | "export" | "survey" => {
            return Ok(None);
        }
        "download" => RunState::Downloading,
        "preprocess" => RunState::PreprocessRunning,
        "quant" => RunState::QuantRunning,
        _ => bail!("unsupported job stage {stage}"),
    }))
}

fn reset_running(
    connection: &Connection,
    stage: &str,
    run: Option<&str>,
    message: &str,
) -> Result<()> {
    if run.is_none() {
        connection.execute(
            "UPDATE project_stages SET state = 'FAILED', last_error = ?1, updated_at = ?2 WHERE state = 'RUNNING' AND stage = ?3",
            params![message, now_rfc3339(), stage],
        )?;
        return Ok(());
    }
    let Some(running) = running_state(stage)? else {
        return Ok(());
    };
    let retry = match stage {
        "download" => RunState::DownloadFailed,
        "preprocess" => RunState::PreprocessFailed,
        "quant" => RunState::QuantFailed,
        _ => unreachable!(),
    };
    connection.execute(
        "UPDATE runs SET state = ?1, last_error = ?2, updated_at = ?3 WHERE state = ?4
         AND (?5 IS NULL OR run_accession = ?5)",
        params![
            retry.as_str(),
            message,
            now_rfc3339(),
            running.as_str(),
            run
        ],
    )?;
    Ok(())
}

fn job_event(
    connection: &Connection,
    db: &Database,
    stage: &str,
    run: Option<&str>,
    message: &str,
    attempt: &str,
) -> Result<()> {
    connection.execute(
        "INSERT INTO events(project_id, ts, stage, run_accession, message, context_json) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
        params![db.project_id(), now_rfc3339(), stage, run, message, serde_json::json!({"attempt_id":attempt}).to_string()],
    )?;
    Ok(())
}

fn now_ms() -> Result<i64> {
    Ok(i64::try_from(
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis(),
    )?)
}
fn duration_ms(duration: Duration) -> Result<i64> {
    Ok(i64::try_from(duration.as_millis())?)
}

#[cfg(test)]
mod tests;
