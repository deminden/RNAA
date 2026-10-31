use super::*;
use std::fs::File;

impl Database {
    fn trash_directory(&self, relative: &Path, create: bool) -> Result<PathBuf> {
        let mut path = fs::canonicalize(self.path().parent().context("project has no root")?)?;
        for part in relative.components() {
            let std::path::Component::Normal(name) = part else {
                bail!("invalid file trash path");
            };
            path.push(name);
            match path.symlink_metadata() {
                Ok(metadata) if !metadata.is_dir() || metadata.file_type().is_symlink() => {
                    bail!("file trash paths require real directories")
                }
                Ok(_) => {}
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                    if create {
                        fs::create_dir(&path)?;
                    }
                }
                Err(error) => return Err(error.into()),
            }
        }
        Ok(path)
    }
    pub(super) fn file_receipt_path(&self, id: &str, create: bool) -> Result<PathBuf> {
        uuid::Uuid::parse_str(id)?;
        let path = self
            .trash_directory(Path::new(".rnaa_support/file_trash_receipts"), create)?
            .join(format!("{id}.json"));
        if path.symlink_metadata().is_ok_and(|meta| !meta.is_file()) {
            bail!("file trash receipt must be a regular file");
        }
        Ok(path)
    }
    pub(super) fn file_trash_path(&self, id: &str, create: bool) -> Result<PathBuf> {
        uuid::Uuid::parse_str(id)?;
        Ok(self
            .trash_directory(Path::new(&format!("trash/files/{id}")), create)?
            .join("data"))
    }
    pub(super) fn original_file_path(&self, value: &Descriptor, create: bool) -> Result<PathBuf> {
        let root = fs::canonicalize(self.path().parent().context("project has no root")?)?;
        let path = Path::new(&value.canonical_source);
        let parent = path.parent().context("cleanup source has no parent")?;
        self.trash_directory(parent.strip_prefix(&root)?, create)?;
        Ok(path.to_path_buf())
    }
    pub(super) fn copy_staging_path(&self, value: &Descriptor, restore: bool) -> Result<PathBuf> {
        let destination = if restore {
            self.original_file_path(value, false)?
        } else {
            self.file_trash_path(&value.trash_id, false)?
        };
        let path = destination
            .parent()
            .context("cleanup destination has no parent")?
            .join(format!(".rnaa-file-copy-{}", value.trash_id));
        // Reserved by the durable intent and checked absent before marking.
        if let Ok(metadata) = path.symlink_metadata() {
            if !metadata.is_file() {
                bail!("cleanup copy staging requires a regular file");
            }
            #[cfg(unix)]
            {
                use std::os::unix::fs::MetadataExt;
                if metadata.nlink() != 1
                    && !(metadata.nlink() == 2
                        && destination.symlink_metadata().is_ok_and(|target| {
                            target.is_file()
                                && target.dev() == metadata.dev()
                                && target.ino() == metadata.ino()
                        }))
                {
                    bail!("cleanup staging has unexpected hardlink aliases");
                }
            }
        }
        Ok(path)
    }
    pub(super) fn verify_trash_file(&self, path: &Path, value: &Descriptor) -> Result<()> {
        if !path.symlink_metadata()?.is_file()
            || fs::metadata(path)?.len() != value.source.bytes
            || compute_sha256(path)? != value.source.checksum
        {
            bail!("file cleanup bytes changed: {}", path.display());
        }
        Ok(())
    }
    pub(super) fn verify_trash_directory(&self, id: &str) -> Result<()> {
        let path = self.file_trash_path(id, false)?;
        if let Some(parent) = path.parent()
            && parent.exists()
        {
            let pending: bool = self.read(|conn| Ok(conn.query_row("SELECT EXISTS(SELECT 1 FROM file_trash WHERE trash_id=?1 AND state IN ('MOVING','RESTORING'))", [id], |row| row.get(0))?))?;
            for entry in fs::read_dir(parent)? {
                let entry = entry?;
                if (entry.file_name() != "data"
                    && !(pending && entry.file_name() == format!(".rnaa-file-copy-{id}").as_str()))
                    || !entry.file_type()?.is_file()
                {
                    bail!("file trash directory contains unexpected data");
                }
            }
        }
        Ok(())
    }
    pub(super) fn move_cleanup_file(&self, value: &Descriptor, restore: bool) -> Result<()> {
        if let Some(producer) = &value.producer
            && self.generation_retired(producer)?
        {
            bail!("restore the file's producer generation first");
        }
        let original = self.original_file_path(value, restore)?;
        let trash = self.file_trash_path(&value.trash_id, !restore)?;
        self.verify_trash_directory(&value.trash_id)?;
        let (from, to) = if restore {
            (&trash, &original)
        } else {
            (&original, &trash)
        };
        match (from.exists(), to.exists()) {
            (true, false) => {
                self.verify_trash_file(from, value)?;
                // Creating a second name refuses an existing destination, unlike
                // rename. The durable intent also recovers a crash between names.
                match fs::hard_link(from, to) {
                    Ok(()) => {}
                    Err(error) if error.kind() == std::io::ErrorKind::CrossesDevices => {
                        let staging = self.copy_staging_path(value, restore)?;
                        let mut output = fs::OpenOptions::new()
                            .write(true)
                            .create(true)
                            .truncate(true)
                            .open(&staging)?;
                        std::io::copy(&mut File::open(from)?, &mut output)?;
                        output.sync_all()?;
                        self.verify_trash_file(&staging, value)?;
                        fs::hard_link(&staging, to)?;
                        fs::remove_file(staging)?;
                    }
                    Err(error) => return Err(error.into()),
                }
                File::open(to.parent().unwrap())?.sync_all()?;
                self.verify_trash_file(to, value)?;
                self.verify_trash_file(from, value)?;
                fs::remove_file(from)?;
                self.verify_trash_file(to, value)?;
            }
            (true, true) => {
                // A cross-device copy may have synced its destination before a
                // crash. Exact matching bytes and the durable move intent permit
                // finishing its source unlink without replacing either copy.
                self.verify_trash_file(from, value)?;
                self.verify_trash_file(to, value)?;
                fs::remove_file(from)?;
            }
            (false, true) => self.verify_trash_file(to, value)?,
            _ => bail!("cleanup input is missing from both original and trash locations"),
        }
        let staging = self.copy_staging_path(value, restore)?;
        if staging.exists() {
            fs::remove_file(staging)?;
        }
        File::open(from.parent().unwrap())?.sync_all()?;
        File::open(to.parent().unwrap())?.sync_all()?;
        Ok(())
    }
}
