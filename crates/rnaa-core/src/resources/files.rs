use super::OwnerFileIdentity;
use anyhow::{Context, Result, bail};
use std::{
    fs::{self, File, OpenOptions, TryLockError},
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

pub(super) fn existing_path(path: &Path) -> Result<PathBuf> {
    let mut current = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };
    loop {
        match fs::canonicalize(&current) {
            Ok(path) => return Ok(path),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                if !current.pop() {
                    return Err(error).context("disk reservation has no existing ancestor");
                }
            }
            Err(error) => return Err(error.into()),
        }
    }
}

pub(super) fn regular(path: &Path, create: bool) -> Result<File> {
    if let Ok(metadata) = path.symlink_metadata()
        && (!metadata.is_file() || metadata.file_type().is_symlink())
    {
        bail!(
            "resource coordination files must be regular files: {}",
            path.display()
        );
    }
    let mut options = OpenOptions::new();
    options
        .read(true)
        .write(true)
        .create(create)
        .truncate(false);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.mode(0o600);
    }
    let file = options
        .open(path)
        .with_context(|| format!("cannot open resource coordination file {}", path.display()))?;
    let named = path.symlink_metadata()?;
    if !file.metadata()?.is_file()
        || named.file_type().is_symlink()
        || identity(&file.metadata()?)? != identity(&named)?
    {
        bail!("resource coordination file changed while opening it");
    }
    Ok(file)
}

#[cfg(unix)]
pub(super) fn identity(metadata: &fs::Metadata) -> Result<OwnerFileIdentity> {
    use std::os::unix::fs::MetadataExt;
    if !metadata.is_file() || metadata.nlink() != 1 {
        bail!("resource owner proof requires an unaliased regular file");
    }
    Ok(OwnerFileIdentity {
        device: metadata.dev(),
        inode: metadata.ino(),
    })
}
#[cfg(not(unix))]
pub(super) fn identity(_metadata: &fs::Metadata) -> Result<OwnerFileIdentity> {
    bail!("resource owner identity is supported on Unix hosts")
}

pub(super) fn remove_owner(path: &Path, expected: OwnerFileIdentity) -> Result<()> {
    let named = match path.symlink_metadata() {
        Ok(named) => named,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error.into()),
    };
    if named.file_type().is_symlink() || named.len() != 0 || identity(&named)? != expected {
        bail!("resource owner lock file identity changed; its pathname is held");
    }
    fs::remove_file(path)?;
    Ok(())
}

pub(super) fn exclusive(path: &Path) -> Result<File> {
    let file = regular(path, true)?;
    let deadline = Instant::now() + Duration::from_secs(15);
    loop {
        match file.try_lock() {
            Ok(()) => return Ok(file),
            Err(TryLockError::WouldBlock) if Instant::now() < deadline => {
                std::thread::sleep(Duration::from_millis(10))
            }
            Err(TryLockError::WouldBlock) => bail!("resource catalog is busy; retry"),
            Err(TryLockError::Error(error)) => return Err(error.into()),
        }
    }
}
pub(super) fn try_owner(path: &Path, expected: Option<OwnerFileIdentity>) -> Result<Option<File>> {
    let file = regular(path, false)?;
    let metadata = file.metadata()?;
    let actual = identity(&metadata)?;
    if metadata.len() != 0 || expected.is_some_and(|expected| actual != expected) {
        bail!("resource owner lock file identity changed; its record is held");
    }
    match file.try_lock() {
        Ok(()) => Ok(Some(file)),
        Err(TryLockError::WouldBlock) => Ok(None),
        Err(TryLockError::Error(error)) => Err(error.into()),
    }
}
