use super::*;
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Debug, Clone, Copy)]
pub struct GcOptions {
    pub apply: bool,
    pub purge_expired: bool,
    pub trash_days: u64,
}
#[derive(Debug, Serialize)]
pub struct GcEntry {
    pub sha256: String,
    pub bytes: u64,
    pub state: String,
    pub action: String,
    pub owners: u64,
    pub legacy_hold: bool,
    pub diagnostics: Vec<String>,
}
#[derive(Debug, Serialize)]
pub struct GcReport {
    pub applied: bool,
    pub entries: Vec<GcEntry>,
    pub errors: Vec<String>,
}
fn now_ms() -> Result<i64> {
    Ok(i64::try_from(
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis(),
    )?)
}

impl SharedStore {
    fn state(&self, sha256: &str) -> Result<(u64, String, bool, Option<i64>)> {
        Ok(self.connection()?.query_row(
            "SELECT bytes,state,legacy_hold,purge_after_ms FROM blobs WHERE sha256=?1",
            [sha256],
            |row| Ok((read_u64(row, 0)?, row.get(1)?, row.get(2)?, row.get(3)?)),
        )?)
    }
    /// The blob lock remains held throughout each filesystem transition. Durable
    /// intent precedes rename/unlink; retries complete the same operation.
    pub(super) fn recover_blob(&self, sha256: &str) -> Result<()> {
        let state: Option<(u64, String)> = self
            .connection()?
            .query_row(
                "SELECT bytes,state FROM blobs WHERE sha256=?1",
                [sha256],
                |row| Ok((read_u64(row, 0)?, row.get(1)?)),
            )
            .optional()?;
        let Some((bytes, state)) = state else {
            return Ok(());
        };
        let live = self.blob_path(sha256)?;
        let trash = self.trash_path(sha256)?;
        match state.as_str() {
            "TRASHING" | "RESTORING" => {
                if state == "TRASHING" {
                    self.verify_legacy_import(sha256)?;
                    self.require_unowned(sha256)?;
                    if live.exists() {
                        self.require_unlinked_blob(&live)?;
                    } else if trash.exists() {
                        self.require_unlinked_blob(&trash)?;
                    }
                }
                let (from, to, completed) = if state == "TRASHING" {
                    (&live, &trash, "TRASHED")
                } else {
                    (&trash, &live, "LIVE")
                };
                match (from.exists(), to.exists()) {
                    (true, false) => {
                        self.verify_file(from, sha256, bytes)?;
                        self.move_file(from, to)?;
                    }
                    (false, true) => {
                        self.verify_file(to, sha256, bytes)?;
                    }
                    _ => bail!("ambiguous shared store move for {sha256}; retain for diagnosis"),
                }
                self.connection()?.execute(
                    "UPDATE blobs SET state=?1 WHERE sha256=?2 AND state=?3",
                    params![completed, sha256, state],
                )?;
            }
            "DELETING" => {
                self.verify_legacy_import(sha256)?;
                self.require_unowned(sha256)?;
                if live.exists() {
                    bail!("refusing shared purge with an unexpected live copy");
                }
                if trash.exists() {
                    self.verify_file(&trash, sha256, bytes)?;
                    self.require_unlinked_blob(&trash)?;
                    fs::remove_file(&trash)?;
                    File::open(trash.parent().unwrap())?.sync_all()?;
                }
                self.finish_blob_deletion(sha256, "DELETING")?;
            }
            "LIVE" | "TRASHED" | "DELETED" => {}
            _ => bail!("invalid shared blob state"),
        }
        Ok(())
    }
    pub fn restore(&self, sha256: &str) -> Result<()> {
        let _lock = self
            .lock(sha256, false)?
            .context("shared blob is busy; retry restore")?;
        self.recover_blob(sha256)?;
        self.restore_locked(sha256)
    }
    pub(super) fn restore_locked(&self, sha256: &str) -> Result<()> {
        let (bytes, state, _, _) = self.state(sha256)?;
        if state == "LIVE" {
            return self.verify_file(&self.blob_path(sha256)?, sha256, bytes);
        }
        if state != "TRASHED" {
            bail!("only a retained trash blob can be restored");
        }
        self.verify_file(&self.trash_path(sha256)?, sha256, bytes)?;
        self.connection()?.execute(
            "UPDATE blobs SET state='RESTORING' WHERE sha256=?1 AND state='TRASHED'",
            [sha256],
        )?;
        self.recover_blob(sha256)
    }
    /// Dry runs do not reconcile owners/pins or complete pending filesystem work.
    /// Unindexed files and legacy holds are never collection candidates.
    pub fn gc(&self, options: GcOptions) -> Result<GcReport> {
        let retention_ms = i64::try_from(
            options
                .trash_days
                .checked_mul(86_400_000)
                .context("trash retention is too large")?,
        )?;
        let mut report = GcReport {
            applied: options.apply,
            entries: Vec::new(),
            errors: Vec::new(),
        };
        for status in self.status()? {
            let mut entry = GcEntry {
                sha256: status.sha256.clone(),
                bytes: status.bytes,
                state: status.state.clone(),
                action: String::new(),
                owners: status.owners,
                legacy_hold: status.legacy_hold,
                diagnostics: Vec::new(),
            };
            if let Err(error) = self.collect_one(&mut entry, options, retention_ms) {
                entry.action = "held_error".into();
                entry.diagnostics.push(format!("{error:#}"));
                report.errors.push(format!("{}: {error:#}", entry.sha256));
            }
            report.entries.push(entry);
        }
        Ok(report)
    }
    fn collect_one(
        &self,
        entry: &mut GcEntry,
        options: GcOptions,
        retention_ms: i64,
    ) -> Result<()> {
        let sha256 = &entry.sha256;
        let Some(_lock) = self.lock(sha256, false)? else {
            entry.action = "held_busy".into();
            return Ok(());
        };
        self.verify_legacy_import(sha256)?;
        if options.apply {
            self.recover_blob(sha256)?;
        }
        let (bytes, state, legacy, purge_after) = self.state(sha256)?;
        entry.state = state.clone();
        let (owners, errors) = self.reconcile_owners(sha256, options.apply)?;
        entry.owners = owners;
        entry.diagnostics = errors;
        if options.apply {
            // Exclusive acquisition proves all durable reader rows are stale.
            self.connection()?
                .execute("DELETE FROM pins WHERE sha256=?1", [sha256])?;
        }
        if owners > 0 {
            entry.action = "held_owned".into();
            return Ok(());
        }
        if legacy {
            entry.action = "held_legacy".into();
            return Ok(());
        }
        let path = if state == "TRASHED" {
            self.trash_path(sha256)?
        } else {
            self.blob_path(sha256)?
        };
        if path.exists() && self.has_additional_hardlinks(&path)? {
            entry.action = "held_external_links".into();
            return Ok(());
        }
        if state == "LIVE" {
            let path = self.blob_path(sha256)?;
            if !path.exists() {
                entry.action = "retire_empty".into();
                if options.apply {
                    self.finish_blob_deletion(sha256, "LIVE")?;
                }
            } else {
                self.verify_file(&path, sha256, bytes)?;
                entry.action = "trash".into();
                if options.apply {
                    // All owner reconciliation is complete before this short mark.
                    let mut connection = self.connection()?;
                    let transaction =
                        connection.transaction_with_behavior(TransactionBehavior::Immediate)?;
                    let protected:bool=transaction.query_row("SELECT EXISTS(SELECT 1 FROM owners WHERE sha256=?1) OR EXISTS(SELECT 1 FROM pins WHERE sha256=?1)",[sha256],|row|row.get(0))?;
                    if protected {
                        bail!("shared ownership changed before collection");
                    }
                    let now = now_ms()?;
                    let expiry = now
                        .checked_add(retention_ms)
                        .context("trash expiry overflow")?;
                    transaction.execute("UPDATE blobs SET state='TRASHING',trashed_at_ms=?1,purge_after_ms=?2 WHERE sha256=?3 AND state='LIVE'",params![now,expiry,sha256])?;
                    transaction.commit()?;
                    self.recover_blob(sha256)?;
                }
            }
        } else if state == "TRASHED" {
            self.verify_file(&self.trash_path(sha256)?, sha256, bytes)?;
            let expired = purge_after.context("trash blob has no retention deadline")? <= now_ms()?;
            if options.purge_expired && expired {
                entry.action = "purge".into();
                if options.apply {
                    self.connection()?.execute(
                        "UPDATE blobs SET state='DELETING' WHERE sha256=?1 AND state='TRASHED'",
                        [sha256],
                    )?;
                    self.recover_blob(sha256)?;
                }
            } else {
                entry.action = "held_trash_window".into();
            }
        } else if state == "DELETED" {
            entry.action = "already_deleted".into();
        } else {
            entry.action = "resume_pending".into();
        }
        if options.apply {
            entry.state = self.state(sha256)?.1;
        }
        Ok(())
    }
    fn require_unlinked_blob(&self, path: &Path) -> Result<()> {
        if self.has_additional_hardlinks(path)? {
            bail!("shared blob has additional hardlinks; preserve unknown or legacy views");
        }
        Ok(())
    }
    fn require_unowned(&self, sha: &str) -> Result<()> {
        let (owners, errors) = self.reconcile_owners(sha, false)?;
        if owners > 0 || !errors.is_empty() {
            bail!("shared ownership changed before pending collection; keep the blob protected");
        }
        Ok(())
    }
    fn finish_blob_deletion(&self, sha: &str, previous: &str) -> Result<()> {
        let mut conn = self.connection()?;
        let tx = conn.transaction_with_behavior(TransactionBehavior::Immediate)?;
        if tx.execute(
            "UPDATE blobs SET state='DELETED',legacy_import_id=NULL WHERE sha256=?1 AND state=?2",
            params![sha, previous],
        )? != 1
        {
            bail!("shared deletion state changed");
        }
        tx.execute("DELETE FROM legacy_blob_imports WHERE sha256=?1", [sha])?;
        tx.commit()?;
        Ok(())
    }
    fn has_additional_hardlinks(&self, path: &Path) -> Result<bool> {
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            Ok(path.symlink_metadata()?.nlink() > 1)
        }
        #[cfg(not(unix))]
        {
            let _ = path;
            Ok(false)
        }
    }
}
