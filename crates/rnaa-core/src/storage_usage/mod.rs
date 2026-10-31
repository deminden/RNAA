//! Observed filesystem usage; logical categories are distinct from allocated inodes.
mod catalog;
#[cfg(test)]
mod tests;
use crate::{Database, generations::GenerationStorePin, util::now_rfc3339};
use anyhow::{Context, Result};
use catalog::Catalog;
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

#[derive(Default, Serialize)]
pub struct CategoryUsage {
    pub file_paths: u64,
    pub logical_bytes: u64,
}
#[derive(Serialize)]
pub struct PathUsage {
    pub path: PathBuf,
    pub category: String,
    pub kind: String,
    pub logical_bytes: u64,
    pub allocated_bytes: Option<u64>,
    pub hardlinks: Option<u64>,
}
#[derive(Serialize)]
pub struct StorageUsageReport {
    pub schema_version: u32,
    pub started_at: String,
    pub finished_at: String,
    pub roots: Vec<PathBuf>,
    pub complete: bool,
    pub catalog_changed: bool,
    pub shared_catalog_changed: Option<bool>,
    pub active_attempts_observed: bool,
    pub regular_file_paths: u64,
    pub directories: u64,
    pub symlinks: u64,
    pub directory_aliases_skipped: u64,
    pub other_entries: u64,
    pub logical_bytes: u64,
    pub distinct_regular_files: Option<u64>,
    pub inode_allocated_bytes: Option<u64>,
    pub observed_hardlink_aliases: Option<u64>,
    pub unobserved_hardlinks: Option<u64>,
    pub categories: BTreeMap<String, CategoryUsage>,
    pub paths: Option<Vec<PathUsage>>,
    pub errors: Vec<String>,
    pub measurement: String,
}
#[cfg(unix)]
type FileKey = (u64, u64);
#[cfg(not(unix))]
type FileKey = PathBuf;
#[cfg(unix)]
fn file_key(metadata: &fs::Metadata, _path: &Path) -> FileKey {
    use std::os::unix::fs::MetadataExt;
    (metadata.dev(), metadata.ino())
}
#[cfg(not(unix))]
fn file_key(_metadata: &fs::Metadata, path: &Path) -> FileKey {
    path.to_path_buf()
}
#[cfg(unix)]
fn allocation(metadata: &fs::Metadata) -> Option<u64> {
    use std::os::unix::fs::MetadataExt;
    metadata.blocks().checked_mul(512)
}
#[cfg(not(unix))]
fn allocation(_metadata: &fs::Metadata) -> Option<u64> {
    None
}
#[cfg(unix)]
fn hardlinks(metadata: &fs::Metadata) -> Option<u64> {
    use std::os::unix::fs::MetadataExt;
    Some(metadata.nlink())
}
#[cfg(not(unix))]
fn hardlinks(_metadata: &fs::Metadata) -> Option<u64> {
    None
}
struct ObservedInode {
    allocated: Option<u64>,
    bytes: u64,
    regular: bool,
    #[cfg(unix)]
    links: Option<u64>,
    paths: u64,
}

impl Database {
    /// Pins project retirement, but does not stop writers or create a filesystem
    /// snapshot. Concurrent activity and incomplete traversal are explicit output.
    pub fn storage_usage(
        &self,
        shared_root: Option<&Path>,
        details: bool,
    ) -> Result<StorageUsageReport> {
        let _pin = GenerationStorePin::acquire(self.path(), "storage accounting")?;
        let project = fs::canonicalize(
            self.path()
                .parent()
                .context("storage accounting project root missing")?,
        )?;
        let shared = shared_root.map(fs::canonicalize).transpose()?;
        let mut catalog = Catalog::capture(self, &project, shared.as_deref())?;
        let mut roots = vec![project.clone()];
        if let Some(shared) = &shared
            && !shared.starts_with(&project)
        {
            if project.starts_with(shared) {
                roots = vec![shared.clone()];
            } else {
                roots.push(shared.clone());
            }
        }
        roots.sort();
        let mut report = StorageUsageReport {
            schema_version: 1,
            started_at: now_rfc3339(), finished_at: String::new(), roots: roots.clone(),
            complete: true, catalog_changed: false, shared_catalog_changed: None,
            active_attempts_observed: catalog.active,
            regular_file_paths: 0, directories: 0, symlinks: 0, directory_aliases_skipped: 0,
            other_entries: 0, logical_bytes: 0, distinct_regular_files: None,
            inode_allocated_bytes: None, observed_hardlink_aliases: None, unobserved_hardlinks: None,
            categories: BTreeMap::new(), paths: details.then(Vec::new), errors: Vec::new(),
            measurement: concat!(
                "Observed regular-file logical bytes; Unix allocation counts distinct (device,inode) ",
                "entries once, including directory/link metadata. Hardlink aliases are deduplicated; ",
                "compression, reflinks and shared extents are not free-space estimates. This is an ",
                "observation, not an admission reservation or atomic filesystem snapshot."
            ).into(),
        };
        let mut inodes = BTreeMap::<FileKey, ObservedInode>::new();
        let mut visited = BTreeSet::<PathBuf>::new();
        let mut pending = roots;
        while let Some(path) = pending.pop() {
            if !visited.insert(path.clone()) {
                continue;
            }
            let metadata = match path.symlink_metadata() {
                Ok(value) => value,
                Err(error) => {
                    report.errors.push(format!("{}: {error}", path.display()));
                    continue;
                }
            };
            let link = metadata.file_type().is_symlink();
            let regular = metadata.is_file() && !link;
            let directory = metadata.is_dir() && !link;
            let kind = if regular {
                "file"
            } else if directory {
                "directory"
            } else if link {
                "symlink"
            } else {
                "other"
            };
            let allocated = allocation(&metadata);
            let links = hardlinks(&metadata);
            let bytes = if regular { metadata.len() } else { 0 };
            let mut directory_alias = false;
            match inodes.entry(file_key(&metadata, &path)) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(ObservedInode {
                        allocated,
                        bytes,
                        regular,
                        #[cfg(unix)]
                        links,
                        paths: 1,
                    });
                }
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    directory_alias = directory;
                    let previous = entry.get_mut();
                    previous.paths += 1;
                    if previous.bytes != bytes
                        || previous.allocated != allocated
                        || previous.regular != regular
                    {
                        report.errors.push(format!(
                            "{}: inode changed during accounting",
                            path.display()
                        ));
                    }
                }
            }
            let category = catalog.classify(&path);
            if regular {
                report.regular_file_paths += 1;
                report.logical_bytes = report
                    .logical_bytes
                    .checked_add(bytes)
                    .context("storage logical-byte total overflow")?;
                let value = report.categories.entry(category.clone()).or_default();
                value.file_paths += 1;
                value.logical_bytes = value
                    .logical_bytes
                    .checked_add(bytes)
                    .context("storage category total overflow")?;
            } else if directory {
                report.directories += 1;
            } else if link {
                report.symlinks += 1;
            } else {
                report.other_entries += 1;
            }
            if let Some(paths) = &mut report.paths {
                paths.push(PathUsage {
                    path: path.clone(),
                    category,
                    kind: kind.into(),
                    logical_bytes: bytes,
                    allocated_bytes: allocated,
                    hardlinks: links,
                });
            }
            if directory_alias {
                report.directory_aliases_skipped += 1;
            }
            if directory && !directory_alias {
                match fs::read_dir(&path) {
                    Ok(entries) => {
                        let mut children = Vec::new();
                        for entry in entries {
                            match entry {
                                Ok(entry) => children.push(entry.path()),
                                Err(error) => {
                                    report.errors.push(format!("{}: {error}", path.display()))
                                }
                            }
                        }
                        children.sort();
                        pending.extend(children.into_iter().rev());
                    }
                    Err(error) => report.errors.push(format!("{}: {error}", path.display())),
                }
            }
        }
        #[cfg(unix)]
        {
            let regular = inodes
                .values()
                .filter(|value| value.regular)
                .collect::<Vec<_>>();
            report.distinct_regular_files = Some(regular.len() as u64);
            report.inode_allocated_bytes = inodes
                .values()
                .try_fold(0u64, |total, value| total.checked_add(value.allocated?));
            report.observed_hardlink_aliases = Some(
                regular
                    .iter()
                    .map(|value| value.paths.saturating_sub(1))
                    .sum(),
            );
            report.unobserved_hardlinks = Some(
                regular
                    .iter()
                    .map(|value| value.links.unwrap_or_default().saturating_sub(value.paths))
                    .sum(),
            );
        }
        #[cfg(not(unix))]
        {
            let _ = inodes;
        }
        let (changed, shared_changed, active) = catalog.finish(self)?;
        report.catalog_changed = changed;
        report.shared_catalog_changed = shared_changed;
        report.active_attempts_observed |= active;
        report.complete = report.errors.is_empty();
        report.finished_at = now_rfc3339();
        if let Some(paths) = &mut report.paths {
            paths.sort_by(|a, b| a.path.cmp(&b.path));
        }
        Ok(report)
    }
}
