use super::*;

pub(super) fn resource_key(run: Option<&str>) -> String {
    run.map(|run| format!("run:{run}"))
        .unwrap_or_else(|| "analysis".into())
}

pub(super) fn writer_lock(db: &Database, resource: &str) -> Result<Option<File>> {
    writer_lock_path(db.path(), resource)
}

pub(crate) fn analysis_upgrade_lock(path: &std::path::Path) -> Result<Option<File>> {
    writer_lock_path(path, &resource_key(None))
}

fn writer_lock_path(path: &std::path::Path, resource: &str) -> Result<Option<File>> {
    let root = path
        .parent()
        .context("database has no parent directory")?
        .join(".rnaa-locks");
    std::fs::create_dir_all(&root)?;
    let digest = hex::encode(Sha256::digest(resource.as_bytes()));
    let file = OpenOptions::new()
        .create(true)
        .truncate(false)
        .read(true)
        .write(true)
        .open(root.join(format!("{digest}.lock")))?;
    match file.try_lock() {
        Ok(()) => Ok(Some(file)),
        Err(TryLockError::WouldBlock) => Ok(None),
        Err(TryLockError::Error(error)) => Err(error).context("failed to lock job writer"),
    }
}
