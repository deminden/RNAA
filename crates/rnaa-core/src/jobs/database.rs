use super::*;

impl Database {
    /// Finish the attempt as the last mutation in an atomic publication session.
    pub fn complete_owned_job(&self, outcome: JobOutcome, message: Option<&str>) -> Result<()> {
        let identity = self
            .job
            .as_ref()
            .context("publication requires an owned job")?;
        self.write(|connection| identity.finish(self, connection, outcome, message))
    }

    /// Reconciliation shares the worker lock; it cannot promote, demote or prune
    /// artifacts while a worker is still writing the same run.
    pub fn with_idle_run<T>(
        &self,
        run: &str,
        operation: impl FnOnce(&Self) -> Result<T>,
    ) -> Result<Option<T>> {
        self.with_idle_resource(&resource_key(Some(run)), operation)
    }

    pub fn with_idle_analysis<T>(
        &self,
        operation: impl FnOnce(&Self) -> Result<T>,
    ) -> Result<Option<T>> {
        self.with_idle_resource(&resource_key(None), operation)
    }

    fn with_idle_resource<T>(
        &self,
        resource: &str,
        operation: impl FnOnce(&Self) -> Result<T>,
    ) -> Result<Option<T>> {
        let _generation_pin =
            crate::generations::GenerationStorePin::acquire(self.path(), resource)?;
        let Some(_lock) = writer_lock(self, resource)? else {
            return Ok(None);
        };
        if !self.write(|connection| expire_resource(connection, self, resource))? {
            return Ok(None);
        }
        operation(self).map(Some)
    }

    pub fn check_job(&self) -> Result<()> {
        if let Some(identity) = &self.job {
            self.read(|connection| identity.check(connection))?;
        }
        Ok(())
    }

    pub fn set_job_request(&self, digest: &str) -> Result<()> {
        let Some(identity) = &self.job else {
            return Ok(());
        };
        if digest.len() != 64 || !digest.bytes().all(|byte| byte.is_ascii_hexdigit()) {
            bail!("job request identity must be a SHA-256 digest");
        }
        self.write(|connection| {
            let changed = connection.execute(
                "UPDATE job_attempts SET request_sha256 = ?1 WHERE attempt_id = ?2
                 AND (request_sha256 IS NULL OR request_sha256 = ?1)",
                params![digest, identity.attempt_id],
            )?;
            if changed != 1 {
                bail!("an attempt cannot change its scientific request identity");
            }
            Ok(())
        })
    }

    pub fn cancel_job(&self, attempt: &str) -> Result<bool> {
        self.write(|connection| Ok(connection.execute(
            "UPDATE job_attempts SET cancel_requested = 1 WHERE project_id = ?1 AND attempt_id = ?2 AND state = 'RUNNING'",
            params![self.project_id(), attempt],
        )? == 1))
    }

    pub fn list_job_attempts(&self, limit: usize) -> Result<Vec<JobAttempt>> {
        self.read(|connection| {
        let mut statement = connection.prepare(
            "SELECT attempt_id, resource, stage, run_accession, request_sha256, owner, fencing_token,
             state, cancel_requested, started_at, heartbeat_ms, lease_expires_ms, finished_at, last_error
             FROM job_attempts WHERE project_id = ?1 ORDER BY rowid DESC LIMIT ?2",
        )?;
        let rows = statement.query_map(
            params![self.project_id(), i64::try_from(limit).unwrap_or(i64::MAX)],
            |row| {
                Ok(JobAttempt {
                    attempt_id: row.get(0)?,
                    resource: row.get(1)?,
                    stage: row.get(2)?,
                    run_accession: row.get(3)?,
                    request_sha256: row.get(4)?,
                    owner: row.get(5)?,
                    fencing_token: row.get(6)?,
                    state: row.get(7)?,
                    cancel_requested: row.get(8)?,
                    started_at: row.get(9)?,
                    heartbeat_ms: row.get(10)?,
                    lease_expires_ms: row.get(11)?,
                    finished_at: row.get(12)?,
                    last_error: row.get(13)?,
                })
            },
        )?;
        rows.map(|row| Ok(row?)).collect()
        })
    }

    /// Recovery and claims use the same writer lock and immediate transaction.
    /// Live processes are left alone, even if a long pause expired their lease.
    pub fn recover_jobs(&self) -> Result<usize> {
        // Recovery visits only active resources, not the full cohort's metadata
        // or the unbounded history of completed attempts.
        let mut resources = std::collections::BTreeSet::new();
        {
            let connection = self.connect()?;
            let mut statement = connection.prepare("SELECT run_accession, state FROM runs
                WHERE state IN ('DOWNLOADING','PREPROCESS_RUNNING','QUANT_RUNNING','DE_RUNNING','CORR_RUNNING')")?;
            let rows = statement.query_map([], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
            })?;
            for row in rows {
                let (run, state) = row?;
                resources.insert(if matches!(state.as_str(), "DE_RUNNING" | "CORR_RUNNING") {
                    resource_key(None)
                } else {
                    resource_key(Some(&run))
                });
            }
            let active_analysis: bool = connection.query_row("SELECT EXISTS(SELECT 1 FROM project_stages WHERE project_id = ?1 AND state = 'RUNNING')",
                params![self.project_id()], |row| row.get(0))?;
            if active_analysis {
                resources.insert(resource_key(None));
            }
            let mut statement = connection.prepare(
                "SELECT resource FROM job_attempts WHERE project_id = ?1 AND state = 'RUNNING'",
            )?;
            let rows =
                statement.query_map(params![self.project_id()], |row| row.get::<_, String>(0))?;
            for row in rows {
                resources.insert(row?);
            }
        }
        let mut recovered = 0;
        for resource in resources {
            let Some(_lock) = writer_lock(self, &resource)? else {
                continue;
            };
            recovered += self.write(|connection| {
                let before = connection.total_changes();
                if !expire_resource(connection, self, &resource)? {
                    return Ok(0);
                }
                // Legacy RUNNING states have no attempt row. Migrate them only
                // while holding the same resource lock used by all new workers.
                if let Some(run) = resource.strip_prefix("run:") {
                    for stage in ["download", "preprocess", "quant"] {
                        reset_running(
                            connection,
                            stage,
                            Some(run),
                            "interrupted worker recovered",
                        )?;
                    }
                } else {
                    for stage in ["normalize", "deseq2", "corr"] {
                        reset_running(connection, stage, None, "interrupted worker recovered")?;
                    }
                }
                let changed = connection.total_changes() > before;
                if changed {
                    connection.execute(
                        "INSERT INTO events(project_id, ts, stage, message, context_json)
                        VALUES (?1, ?2, 'recovery', 'interrupted work recovered', ?3)",
                        params![
                            self.project_id(),
                            now_rfc3339(),
                            serde_json::json!({"resource":resource}).to_string()
                        ],
                    )?;
                }
                Ok(usize::from(changed))
            })?;
        }
        Ok(recovered)
    }
}
