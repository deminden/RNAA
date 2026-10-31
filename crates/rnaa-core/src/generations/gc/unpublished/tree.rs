//! Whole-container observations include unsealed files and empty directories.
use super::*;
use std::fs::File;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct TreeEntry {
    pub path: PathBuf,
    pub directory: bool,
    pub bytes: u64,
    pub sha256: Option<String>,
}
pub(super) fn observe(root: &Path, flush: bool) -> Result<Vec<TreeEntry>> {
    let mut entries = Vec::new();
    let mut directories = Vec::new();
    for entry in walkdir::WalkDir::new(root).follow_links(false) {
        let entry = entry?;
        let metadata = entry.path().symlink_metadata()?;
        if metadata.file_type().is_symlink() || (!metadata.is_dir() && !metadata.is_file()) {
            bail!("unpublished cleanup requires real directories and regular files");
        }
        if metadata.is_dir() {
            directories.push(entry.path().to_path_buf());
        } else if flush {
            File::open(entry.path())?.sync_all()?;
        }
        if entry.depth() > 0 {
            entries.push(TreeEntry {
                path: entry.path().strip_prefix(root)?.into(),
                directory: metadata.is_dir(),
                bytes: if metadata.is_dir() { 0 } else { metadata.len() },
                sha256: if metadata.is_file() {
                    Some(compute_sha256(entry.path())?)
                } else {
                    None
                },
            });
        }
    }
    if flush {
        for directory in directories.into_iter().rev() {
            File::open(directory)?.sync_all()?;
        }
    }
    entries.sort_by(|a, b| a.path.cmp(&b.path));
    Ok(entries)
}
pub(super) fn verify(root: &Path, expected: &[TreeEntry], partial: bool) -> Result<()> {
    let actual = observe(root, false)?;
    if if partial {
        let expected = expected
            .iter()
            .map(|entry| (&entry.path, entry))
            .collect::<BTreeMap<_, _>>();
        actual
            .iter()
            .any(|entry| expected.get(&entry.path).copied() != Some(entry))
    } else {
        actual != expected
    } {
        bail!("unpublished container changed from its trash receipt");
    }
    Ok(())
}

/// Check each parent before creating a child; never create through a symlink.
pub(super) fn checked_directory(db: &Database, names: &[&str], create: bool) -> Result<PathBuf> {
    let mut path = fs::canonicalize(db.path().parent().context("project root missing")?)?;
    for name in names {
        path.push(name);
        match path.symlink_metadata() {
            Ok(meta) if !meta.is_dir() || meta.file_type().is_symlink() => {
                bail!("unpublished cleanup paths must be real directories")
            }
            Ok(_) => {}
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                if create {
                    fs::create_dir(&path)?;
                    File::open(path.parent().unwrap())?.sync_all()?;
                }
            }
            Err(error) => return Err(error.into()),
        }
    }
    Ok(path)
}
pub(super) fn container(db: &Database, id: &str, trash: bool, create: bool) -> Result<PathBuf> {
    uuid::Uuid::parse_str(id)?;
    let parent = checked_directory(
        db,
        if trash {
            &["trash", "unpublished"]
        } else {
            &[".rnaa_generations"]
        },
        create,
    )?;
    let path = parent.join(id);
    match path.symlink_metadata() {
        Ok(meta) if !meta.is_dir() || meta.file_type().is_symlink() => {
            bail!("unpublished container must be a real directory")
        }
        Ok(_) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error.into()),
    }
    Ok(path)
}
pub(super) fn sync_move(source: &Path, target: &Path) -> Result<()> {
    fs::rename(source, target)?;
    File::open(source.parent().context("source parent missing")?)?.sync_all()?;
    File::open(target.parent().context("target parent missing")?)?.sync_all()?;
    Ok(())
}
