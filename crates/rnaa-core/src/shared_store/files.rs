use super::*;
use std::fs::{OpenOptions, TryLockError};

pub(super) fn validate_hash(value: &str) -> Result<()> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        bail!("shared blob identity must be a lowercase SHA-256");
    }
    Ok(())
}
impl SharedStore {
    /// Serialize journal-mode selection as well as migrations. SQLite's busy
    /// handler need not wait on a conflicting journal-mode upgrade during open.
    pub(super) fn initialization_lock(&self) -> Result<File> {
        let path = self.checked_path(Path::new(".locks/catalog-init.lock"), false)?;
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(path)?;
        let started = std::time::Instant::now();
        loop {
            match file.try_lock() {
                Ok(()) => return Ok(file),
                Err(TryLockError::WouldBlock) if started.elapsed() < Duration::from_secs(15) => {
                    std::thread::sleep(Duration::from_millis(25));
                }
                Err(TryLockError::WouldBlock) => {
                    bail!("shared catalog initialization is busy; retry opening the store")
                }
                Err(TryLockError::Error(error)) => return Err(error.into()),
            }
        }
    }
    pub(super) fn checked_path(&self, relative: &Path, create: bool) -> Result<PathBuf> {
        let mut path = self.root.clone();
        for component in relative.components() {
            let std::path::Component::Normal(name) = component else {
                bail!("invalid shared store relative path");
            };
            path.push(name);
            match path.symlink_metadata() {
                Ok(metadata) if metadata.file_type().is_symlink() => {
                    bail!("shared store paths must not traverse symlinks")
                }
                Ok(_) => {}
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                    if create {
                        match fs::create_dir(&path) {
                            Ok(()) => {}
                            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
                            Err(error) => return Err(error.into()),
                        }
                        let metadata = path.symlink_metadata()?;
                        if !metadata.is_dir() || metadata.file_type().is_symlink() {
                            bail!("shared directory was replaced while creating it");
                        }
                    }
                }
                Err(error) => return Err(error.into()),
            }
        }
        Ok(path)
    }
    pub(super) fn blob_path(&self, sha256: &str) -> Result<PathBuf> {
        validate_hash(sha256)?;
        self.checked_path(
            Path::new(&format!(
                "blobs/sha256/{}/{}/{sha256}",
                &sha256[..2],
                &sha256[2..4]
            )),
            false,
        )
    }
    pub(super) fn trash_path(&self, sha256: &str) -> Result<PathBuf> {
        validate_hash(sha256)?;
        self.checked_path(Path::new(&format!("trash/{sha256}.blob")), false)
    }
    pub(super) fn lock(&self, sha256: &str, shared: bool) -> Result<Option<File>> {
        validate_hash(sha256)?;
        let path = self.checked_path(Path::new(&format!(".locks/{sha256}.lock")), false)?;
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(path)?;
        let result = if shared {
            file.try_lock_shared()
        } else {
            file.try_lock()
        };
        match result {
            Ok(()) => Ok(Some(file)),
            Err(TryLockError::WouldBlock) => Ok(None),
            Err(TryLockError::Error(error)) => Err(error.into()),
        }
    }
    pub(super) fn move_file(&self, from: &Path, to: &Path) -> Result<()> {
        if to.exists() {
            bail!("shared move destination already exists");
        }
        self.checked_path(
            to.parent()
                .context("shared path has no parent")?
                .strip_prefix(&self.root)?,
            true,
        )?;
        fs::rename(from, to)?;
        File::open(from.parent().unwrap())?.sync_all()?;
        File::open(to.parent().unwrap())?.sync_all()?;
        Ok(())
    }
}
