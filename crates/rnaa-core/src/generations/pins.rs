//! Project maintenance exclusion. Locks precede SQLite transactions everywhere.
use crate::{Database, util::now_rfc3339};
use anyhow::{Context, Result, bail};
use rusqlite::{Connection, OpenFlags, params};
use std::{
    fs::{self, File, OpenOptions, TryLockError},
    path::{Path, PathBuf},
    sync::Arc,
    time::Duration,
};

fn lock_file(db_path: &Path) -> Result<File> {
    let root = fs::canonicalize(
        db_path
            .parent()
            .context("generation database has no root")?,
    )?;
    let directory = root.join(".rnaa-locks");
    fs::create_dir_all(&directory)?;
    if directory.symlink_metadata()?.file_type().is_symlink() {
        bail!("generation lock directory must not be a symlink");
    }
    let path = directory.join("generation-store.lock");
    if path.symlink_metadata().is_ok_and(|meta| !meta.is_file()) {
        bail!("generation lock must be a regular file");
    }
    Ok(OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(path)?)
}
pub(crate) fn read_lock(db_path: &Path) -> Result<File> {
    let lock = lock_file(db_path)?;
    let started = std::time::Instant::now();
    loop {
        match lock.try_lock_shared() {
            Ok(()) => return Ok(lock),
            Err(TryLockError::WouldBlock) if started.elapsed() < Duration::from_secs(15) => {
                std::thread::sleep(Duration::from_millis(25))
            }
            Err(TryLockError::WouldBlock) => bail!("generation maintenance is busy; retry"),
            Err(TryLockError::Error(error)) => return Err(error.into()),
        }
    }
}
pub(super) fn maintenance_lock(db: &Database) -> Result<Option<File>> {
    if db.job.is_some() {
        bail!("generation maintenance cannot run inside an owned worker");
    }
    maintenance_lock_path(db.path())
}
/// Legacy catalog inspection uses the same exclusion without migrating or writing
/// the inspected SQLite database. Callers must first verify the catalog exists.
pub(crate) fn maintenance_lock_path(path: &Path) -> Result<Option<File>> {
    let lock = lock_file(path)?;
    match lock.try_lock() {
        Ok(()) => Ok(Some(lock)),
        Err(TryLockError::WouldBlock) => Ok(None),
        Err(TryLockError::Error(error)) => Err(error.into()),
    }
}
#[derive(Debug)]
pub struct GenerationStorePin {
    path: PathBuf,
    pin: String,
    _lock: File,
}
impl GenerationStorePin {
    pub fn acquire(db_path: &Path, owner: &str) -> Result<Arc<Self>> {
        let lock = read_lock(db_path)?;
        let pin = uuid::Uuid::new_v4().to_string();
        let connection = Connection::open_with_flags(db_path, OpenFlags::SQLITE_OPEN_READ_WRITE)?;
        connection.busy_timeout(Duration::ZERO)?;
        // A containing publication transaction may hold the SQLite writer. The
        // OS pin still protects its reader; durable row insertion is best-effort
        // in that case, while normal snapshots record an inspectable pin.
        match connection.execute(
            "INSERT INTO generation_store_pins(pin_id,owner,created_at) VALUES (?1,?2,?3)",
            params![pin, owner, now_rfc3339()],
        ) {
            Ok(_) => {}
            Err(rusqlite::Error::SqliteFailure(error, _))
                if matches!(
                    error.code,
                    rusqlite::ErrorCode::DatabaseBusy | rusqlite::ErrorCode::DatabaseLocked
                ) => {}
            Err(error) => return Err(error.into()),
        }
        Ok(Arc::new(Self {
            path: db_path.to_path_buf(),
            pin,
            _lock: lock,
        }))
    }
}
impl Drop for GenerationStorePin {
    fn drop(&mut self) {
        if let Ok(connection) =
            Connection::open_with_flags(&self.path, OpenFlags::SQLITE_OPEN_READ_WRITE)
        {
            let _ = connection.busy_timeout(Duration::ZERO);
            let _ = connection.execute(
                "DELETE FROM generation_store_pins WHERE pin_id=?1",
                [&self.pin],
            );
        }
    }
}
