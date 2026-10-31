//! Store-wide ownership with a pinned handoff from project publication.
//! Blob locks are always taken before short catalog transactions. Project SQLite
//! reads never run inside a store transaction; unresolved projects retain owners.
mod files;
mod gc;
mod legacy;
mod ownership;
pub use legacy::{LegacyImportEntry, LegacyImportOptions, LegacyImportReport};
#[cfg(test)]
mod tests;
use crate::{
    Database, SharedBlobRecord,
    util::{compute_sha256, file_size, now_rfc3339, publish_file},
};
use anyhow::{Context, Result, bail};
pub use gc::{GcEntry, GcOptions, GcReport};
use rusqlite::{Connection, OptionalExtension, TransactionBehavior, params};
use serde::Serialize;
use std::{
    fs::{self, File},
    path::{Path, PathBuf},
    time::Duration,
};

#[derive(Clone, Debug)]
pub struct SharedStore {
    root: PathBuf,
}

fn read_u64(row: &rusqlite::Row<'_>, index: usize) -> rusqlite::Result<u64> {
    let value: i64 = row.get(index)?;
    u64::try_from(value).map_err(|error| {
        rusqlite::Error::FromSqlConversionFailure(
            index,
            rusqlite::types::Type::Integer,
            Box::new(error),
        )
    })
}

/// Holds writer exclusion until the project transaction and store handoff finish.
/// A crash leaves a durable PREPARED owner which collection must reconcile.
pub struct BlobPublication {
    _project_pin: std::sync::Arc<crate::generations::GenerationStorePin>,
    store: SharedStore,
    owner: ownership::Owner,
    record: SharedBlobRecord,
    _lock: File,
}
impl BlobPublication {
    pub fn record(&self) -> &SharedBlobRecord {
        &self.record
    }
    pub fn finish(&self) -> Result<()> {
        let ownership::Proof::Present(generation) = ownership::prove(&self.store, &self.owner)?
        else {
            bail!("shared blob publication is not committed in its project");
        };
        let changed = self.store.connection()?.execute(
            "UPDATE owners SET state='COMMITTED',generation_id=?1
             WHERE sha256=?2 AND project_db=?3 AND project_id=?4 AND artifact_path=?5",
            params![
                generation,
                self.owner.sha256,
                self.owner.project_db,
                self.owner.project_id,
                self.owner.artifact_path
            ],
        )?;
        if changed != 1 {
            bail!("shared publication reservation changed before finalization");
        }
        Ok(())
    }
}

/// A reader keeps both an inspectable durable pin and a shared OS lock. Expiry or
/// a stopped heartbeat can never authorize deletion underneath a live reader.
pub struct BlobReader {
    store: SharedStore,
    pin: String,
    record: SharedBlobRecord,
    _lock: File,
}
impl BlobReader {
    pub fn record(&self) -> &SharedBlobRecord {
        &self.record
    }
}
impl Drop for BlobReader {
    fn drop(&mut self) {
        if let Ok(connection) = self.store.connection() {
            let _ = connection.execute("DELETE FROM pins WHERE pin_id=?1", [&self.pin]);
        }
    }
}

#[derive(Debug, Serialize)]
pub struct BlobStatus {
    pub sha256: String,
    pub bytes: u64,
    pub state: String,
    pub legacy_hold: bool,
    pub owners: u64,
    pub prepared_owners: u64,
    pub reader_pins: u64,
    pub trashed_at_ms: Option<i64>,
    pub purge_after_ms: Option<i64>,
    pub legacy_import_id: Option<String>,
}

impl SharedStore {
    pub fn open(root: &Path) -> Result<Self> {
        fs::create_dir_all(root)?;
        let store = Self {
            root: fs::canonicalize(root)?,
        };
        for name in ["blobs/sha256", "trash", ".locks", ".imports"] {
            store.checked_path(Path::new(name), true)?;
        }
        // WAL mode changes precede the migration transaction and need their own
        // exclusion. Release this lock before any blob publication/GC locking.
        let _initialization_lock = store.initialization_lock()?;
        let path = store.root.join("catalog.sqlite3");
        if path
            .symlink_metadata()
            .is_ok_and(|metadata| !metadata.is_file())
        {
            bail!("shared catalog must be a regular file");
        }
        let mut connection = store.connection()?;
        connection
            .pragma_update(None, "journal_mode", "WAL")
            .context("failed to select shared catalog WAL mode")?;
        let transaction = connection.transaction_with_behavior(TransactionBehavior::Immediate)?;
        // Read the version after taking SQLite's writer: simultaneous initializers
        // must observe the upgrade completed by the previous lock holder.
        let version: i64 =
            transaction.pragma_query_value(None, "user_version", |row| row.get(0))?;
        if version > 2 {
            bail!("unsupported shared store schema {version}");
        }
        if version < 1 {
            transaction.execute_batch(include_str!("schema.sql"))?;
        }
        if version < 2 {
            transaction.execute_batch(include_str!("legacy_schema.sql"))?;
            transaction.execute_batch("ALTER TABLE blobs ADD COLUMN legacy_import_id TEXT REFERENCES legacy_imports(import_id);")?;
        }
        transaction.commit()?;
        Ok(store)
    }
    pub fn root(&self) -> &Path {
        &self.root
    }
    fn connection(&self) -> Result<Connection> {
        let connection = Connection::open(self.root.join("catalog.sqlite3"))?;
        connection.busy_timeout(Duration::from_secs(15))?;
        connection.pragma_update(None, "foreign_keys", true)?;
        connection.pragma_update(None, "synchronous", "FULL")?;
        Ok(connection)
    }
    pub fn status(&self) -> Result<Vec<BlobStatus>> {
        let connection = self.connection()?;
        let mut query = connection.prepare(
            "SELECT b.sha256,b.bytes,b.state,b.legacy_hold,b.trashed_at_ms,
            (SELECT count(*) FROM owners o WHERE o.sha256=b.sha256),
            (SELECT count(*) FROM owners o WHERE o.sha256=b.sha256 AND o.state='PREPARED'),
            (SELECT count(*) FROM pins p WHERE p.sha256=b.sha256),b.purge_after_ms,b.legacy_import_id
            FROM blobs b ORDER BY b.sha256",
        )?;
        Ok(query
            .query_map([], |row| {
                Ok(BlobStatus {
                    sha256: row.get(0)?,
                    bytes: read_u64(row, 1)?,
                    state: row.get(2)?,
                    legacy_hold: row.get(3)?,
                    trashed_at_ms: row.get(4)?,
                    owners: read_u64(row, 5)?,
                    prepared_owners: read_u64(row, 6)?,
                    reader_pins: read_u64(row, 7)?,
                    purge_after_ms: row.get(8)?,
                    legacy_import_id: row.get(9)?,
                })
            })?
            .collect::<std::result::Result<Vec<_>, _>>()?)
    }
    pub fn prepare(
        &self,
        db: &Database,
        source: &Path,
        sha256: &str,
        bytes: u64,
    ) -> Result<BlobPublication> {
        let project_pin =
            crate::generations::GenerationStorePin::acquire(db.path(), "shared blob publication")?;
        db.check_job()?;
        files::validate_hash(sha256)?;
        let started = std::time::Instant::now();
        let lock = loop {
            db.check_job()?;
            if let Some(lock) = self.lock(sha256, false)? {
                break lock;
            }
            if started.elapsed() >= Duration::from_secs(15) {
                bail!("shared blob is busy; retry publication");
            }
            std::thread::sleep(Duration::from_millis(25));
        };
        self.recover_blob(sha256)?;
        let path = self.blob_path(sha256)?;
        let project_db = fs::canonicalize(db.path())?.display().to_string();
        let owner = ownership::Owner {
            sha256: sha256.into(),
            project_db,
            project_id: db.project_id().into(),
            artifact_path: source.display().to_string(),
            state: "PREPARED".into(),
            generation_id: ownership::source_generation(db, source)?,
        };
        let mut connection = self.connection()?;
        let transaction = connection.transaction_with_behavior(TransactionBehavior::Immediate)?;
        let existing: Option<(u64, String)> = transaction
            .query_row(
                "SELECT bytes,state FROM blobs WHERE sha256=?1",
                [sha256],
                |row| Ok((read_u64(row, 0)?, row.get(1)?)),
            )
            .optional()?;
        if let Some((length, state)) = &existing {
            if *length != bytes {
                bail!("shared catalog size differs from requested blob");
            }
            if state == "TRASHED" {
                drop(transaction);
                drop(connection);
                self.restore_locked(sha256)?;
            } else if state != "LIVE" && state != "DELETED" {
                bail!("shared blob has an unresolved operation");
            } else {
                transaction.commit()?;
            }
        } else {
            transaction.execute("INSERT INTO blobs(sha256,bytes,state,legacy_hold,created_at) VALUES (?1,?2,'LIVE',?3,?4)",
                params![sha256,i64::try_from(bytes)?,path.exists(),now_rfc3339()])?;
            transaction.commit()?;
        }
        // This durable reservation precedes the first shared payload byte.
        let mut connection = self.connection()?;
        let transaction = connection.transaction_with_behavior(TransactionBehavior::Immediate)?;
        transaction.execute(
            "INSERT INTO owners(sha256,project_db,project_id,artifact_path,state,created_at,generation_id)
            VALUES (?1,?2,?3,?4,'PREPARED',?5,?6) ON CONFLICT DO NOTHING",
            params![
                sha256,
                owner.project_db,
                owner.project_id,
                owner.artifact_path,
                now_rfc3339()
                ,owner.generation_id
            ],
        )?;
        transaction.execute(
            "UPDATE blobs SET state='LIVE',trashed_at_ms=NULL,purge_after_ms=NULL WHERE sha256=?1 AND state='DELETED'",
            [sha256],
        )?;
        transaction.commit()?;
        if !path.exists() {
            let parent = path.parent().context("shared blob has no parent")?;
            self.checked_path(parent.strip_prefix(&self.root)?, true)?;
            publish_file(&path, |temporary| {
                fs::copy(source, temporary)?;
                Ok(())
            })?;
        }
        self.verify_file(&path, sha256, bytes)?;
        db.check_job()?;
        Ok(BlobPublication {
            _project_pin: project_pin,
            store: self.clone(),
            owner,
            record: self.record(sha256, bytes, &path),
            _lock: lock,
        })
    }
    pub fn pin(&self, sha256: &str, owner: &str) -> Result<BlobReader> {
        files::validate_hash(sha256)?;
        let lock = self
            .lock(sha256, true)?
            .context("shared blob is being published or collected; retry")?;
        let connection = self.connection()?;
        let (bytes, state): (u64, String) = connection.query_row(
            "SELECT bytes,state FROM blobs WHERE sha256=?1",
            [sha256],
            |row| Ok((read_u64(row, 0)?, row.get(1)?)),
        )?;
        if state != "LIVE" {
            bail!("shared blob is not live; restore it before reading");
        }
        let path = self.blob_path(sha256)?;
        self.verify_file(&path, sha256, bytes)?;
        let pin = uuid::Uuid::new_v4().to_string();
        connection.execute(
            "INSERT INTO pins(pin_id,sha256,owner,created_at) VALUES (?1,?2,?3,?4)",
            params![pin, sha256, owner, now_rfc3339()],
        )?;
        Ok(BlobReader {
            store: self.clone(),
            pin,
            record: self.record(sha256, bytes, &path),
            _lock: lock,
        })
    }
    fn verify_file(&self, path: &Path, sha256: &str, bytes: u64) -> Result<()> {
        if !path.symlink_metadata()?.is_file()
            || file_size(path)? != bytes
            || compute_sha256(path)? != sha256
        {
            bail!(
                "shared blob failed size/hash verification: {}",
                path.display()
            );
        }
        Ok(())
    }
    fn record(&self, sha256: &str, bytes: u64, path: &Path) -> SharedBlobRecord {
        SharedBlobRecord {
            blob_id: format!("sha256:{sha256}"),
            storage_path: path.display().to_string(),
            checksum_type: "sha256".into(),
            checksum: sha256.into(),
            bytes,
            created_at: now_rfc3339(),
        }
    }
}
