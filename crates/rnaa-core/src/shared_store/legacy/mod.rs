//! Import store ownership from an explicitly declared inventory of idle projects.
mod catalog;
mod projects;
#[cfg(test)]
mod tests;
use super::*;
use crate::{ArtifactRecord, generations::canonical_missing_path};
use projects::ProjectEvidence;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
fn digest(bytes: impl AsRef<[u8]>) -> String {
    hex::encode(Sha256::digest(bytes.as_ref()))
}
fn canonical(path: &Path) -> Result<PathBuf> {
    canonical_missing_path(path)
}

#[derive(Clone, Copy)]
pub struct LegacyImportOptions {
    pub apply: bool,
    /// Operator assertion that all legacy projects using this store are supplied.
    pub complete_project_inventory: bool,
}
#[derive(Serialize)]
pub struct LegacyImportEntry {
    pub sha256: String,
    pub bytes: Option<u64>,
    pub action: String,
    pub owners: u64,
    pub releases_legacy_hold: bool,
    pub observed_hardlinks: Option<u64>,
    pub accounted_hardlinks: u64,
    pub reasons: Vec<String>,
}
#[derive(Serialize)]
pub struct LegacyImportReport {
    pub schema_version: u32,
    pub import_id: Option<String>,
    pub applied: bool,
    pub complete: bool,
    pub entries: Vec<LegacyImportEntry>,
    pub busy_projects: Vec<PathBuf>,
    pub unindexed_paths: Vec<PathBuf>,
    pub errors: Vec<String>,
}
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Descriptor {
    schema_version: u32,
    import_id: String,
    store_root: PathBuf,
    complete_project_inventory: bool,
    projects: Vec<ProjectEvidence>,
    blobs: Vec<BlobEvidence>,
    created_at: String,
}
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct BlobEvidence {
    sha256: String,
    bytes: u64,
    owners: Vec<ImportedOwner>,
    observed_hardlinks: Option<u64>,
    accounted_hardlinks: u64,
    was_legacy_hold: bool,
    release_hold: bool,
}
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ImportedOwner {
    project_db: PathBuf,
    project_id: String,
    artifact_path: String,
    generation_id: Option<String>,
}
impl ImportedOwner {
    fn owner(&self, sha: &str) -> ownership::Owner {
        ownership::Owner {
            sha256: sha.into(),
            project_db: self.project_db.display().to_string(),
            project_id: self.project_id.clone(),
            artifact_path: self.artifact_path.clone(),
            state: "COMMITTED".into(),
            generation_id: self.generation_id.clone(),
        }
    }
}
#[cfg(unix)]
fn inode(path: &Path) -> Result<(u64, u64, u64)> {
    use std::os::unix::fs::MetadataExt;
    let meta = path.symlink_metadata()?;
    if !meta.is_file() {
        bail!("legacy shared inode must be a regular file");
    }
    Ok((meta.dev(), meta.ino(), meta.nlink()))
}
fn hardlink_census(path: &Path, sources: &[PathBuf]) -> Result<(Option<u64>, u64)> {
    #[cfg(unix)]
    {
        let (device, number, links) = inode(path)?;
        let mut paths = BTreeSet::from([fs::canonicalize(path)?]);
        for source in sources {
            match inode(source) {
                Ok((dev, ino, _)) if dev == device && ino == number => {
                    paths.insert(fs::canonicalize(source)?);
                }
                Ok(_) => {}
                Err(_) if !source.exists() => {}
                Err(error) => return Err(error),
            }
        }
        Ok((Some(links), u64::try_from(paths.len())?))
    }
    #[cfg(not(unix))]
    {
        let _ = (path, sources);
        Ok((None, 0))
    }
}
impl SharedStore {
    pub fn import_legacy(
        &self,
        projects: &[PathBuf],
        options: LegacyImportOptions,
    ) -> Result<LegacyImportReport> {
        if projects.is_empty() {
            bail!("provide the explicit legacy project inventory with --project");
        }
        let mut report = LegacyImportReport {
            schema_version: 1,
            import_id: None,
            applied: options.apply,
            complete: false,
            entries: Vec::new(),
            busy_projects: Vec::new(),
            unindexed_paths: Vec::new(),
            errors: Vec::new(),
        };
        let mut paths = BTreeSet::new();
        for project in projects {
            let path = fs::canonicalize(project)?.join("state.sqlite");
            if !path.symlink_metadata()?.is_file() {
                bail!("declared legacy project has no regular state.sqlite");
            }
            if !paths.insert(path) {
                bail!("legacy project inventory contains duplicate paths");
            }
        }
        let mut project_locks = Vec::new();
        for path in &paths {
            let Some(lock) = crate::generations::pins::maintenance_lock_path(path)? else {
                report.busy_projects.push(path.clone());
                continue;
            };
            project_locks.push(lock);
        }
        if !report.busy_projects.is_empty() {
            return Ok(report);
        }
        let mut evidence = Vec::new();
        let mut ids = BTreeSet::new();
        for path in &paths {
            match self.inspect_legacy_project(path) {
                Ok(project) => {
                    if !ids.insert(project.project_id.clone()) {
                        report.errors.push("declared project identities are duplicated; resolve copied catalogs before import".into());
                    }
                    evidence.push(project);
                }
                Err(error) => report.errors.push(format!("{}: {error:#}", path.display())),
            }
        }
        let inventory_valid = report.errors.is_empty();
        let omissions = catalog::known_projects(&self.connection()?)?
            .into_iter()
            .filter(|(path, id)| {
                !evidence
                    .iter()
                    .any(|project| project.database == *path && project.project_id == *id)
            })
            .map(|(path, _)| format!("undeclared_known_project:{}", path.display()))
            .collect::<Vec<_>>();
        let mut candidates = self
            .status()?
            .into_iter()
            .filter(|status| status.legacy_hold)
            .map(|status| status.sha256)
            .collect::<BTreeSet<_>>();
        for project in &evidence {
            for source in &project.sources {
                candidates.insert(source.artifact.checksum.clone());
            }
        }
        for entry in walkdir::WalkDir::new(self.root.join("blobs/sha256")).follow_links(false) {
            let entry = entry?;
            if entry.file_type().is_dir() {
                continue;
            }
            let name = entry.file_name().to_string_lossy().into_owned();
            if entry.file_type().is_file()
                && files::validate_hash(&name).is_ok()
                && self.blob_path(&name).is_ok_and(|path| path == entry.path())
            {
                candidates.insert(name);
            } else {
                report.unindexed_paths.push(entry.path().into());
            }
        }
        let mut locks = Vec::new();
        let mut blobs = Vec::new();
        for sha in candidates {
            let mut entry = LegacyImportEntry {
                sha256: sha.clone(),
                bytes: None,
                action: String::new(),
                owners: 0,
                releases_legacy_hold: false,
                observed_hardlinks: None,
                accounted_hardlinks: 0,
                reasons: Vec::new(),
            };
            let result = (|| -> Result<()> {
                let Some(lock) = self.lock(&sha, false)? else {
                    entry.action = "held_busy".into();
                    return Ok(());
                };
                locks.push(lock);
                self.verify_legacy_import(&sha)?;
                let existing: Option<(u64, String, bool)> = self
                    .connection()?
                    .query_row(
                        "SELECT bytes,state,legacy_hold FROM blobs WHERE sha256=?1",
                        [&sha],
                        |row| Ok((read_u64(row, 0)?, row.get(1)?, row.get(2)?)),
                    )
                    .optional()?;
                if existing
                    .as_ref()
                    .is_some_and(|(_, state, _)| state != "LIVE")
                {
                    bail!(
                        "restore or finish pending shared maintenance before importing this blob"
                    );
                }
                let blob = self.blob_path(&sha)?;
                let bytes = file_size(&blob)?;
                self.verify_file(&blob, &sha, bytes)?;
                if existing
                    .as_ref()
                    .is_some_and(|(expected, _, _)| *expected != bytes)
                {
                    bail!("legacy shared catalog size differs from payload");
                }
                let was_legacy_hold = existing.as_ref().is_none_or(|(_, _, legacy)| *legacy);
                entry.bytes = Some(bytes);
                let mut owners = Vec::new();
                let mut physical = Vec::new();
                for project in &evidence {
                    for source in project
                        .sources
                        .iter()
                        .filter(|source| source.artifact.checksum == sha)
                    {
                        if source.artifact.bytes != bytes {
                            bail!("legacy project size differs from shared blob");
                        }
                        owners.push(ImportedOwner {
                            project_db: project.database.clone(),
                            project_id: project.project_id.clone(),
                            artifact_path: source.artifact.path.clone(),
                            generation_id: source.generation_id.clone(),
                        });
                        physical.push(source.physical_path.clone());
                    }
                }
                for owner in self.owners(&sha)? {
                    if !evidence.iter().any(|project| {
                        project.database == Path::new(&owner.project_db)
                            && project.project_id == owner.project_id
                    }) {
                        entry
                            .reasons
                            .push(format!("undeclared_known_project:{}", owner.project_db));
                    }
                }
                let (links, accounted) = hardlink_census(&blob, &physical)?;
                entry.reasons.extend(omissions.iter().cloned());
                entry.observed_hardlinks = links;
                entry.accounted_hardlinks = accounted;
                if !options.complete_project_inventory {
                    entry
                        .reasons
                        .push("complete_project_inventory_not_asserted".into());
                }
                if !inventory_valid {
                    entry
                        .reasons
                        .push("incomplete_or_invalid_project_inventory".into());
                }
                if links.is_none_or(|links| links != accounted) {
                    entry
                        .reasons
                        .push("unexplained_or_unavailable_hardlinks".into());
                }
                entry.owners = u64::try_from(owners.len())?;
                entry.releases_legacy_hold = entry.reasons.is_empty();
                let already_imported: bool = self.connection()?.query_row(
                    "SELECT EXISTS(SELECT 1 FROM legacy_blob_imports WHERE sha256=?1)",
                    [&sha],
                    |row| row.get(0),
                )?;
                if already_imported && !was_legacy_hold && !entry.releases_legacy_hold {
                    entry.action = "held_existing_import".into();
                    return Ok(());
                }
                entry.action = "import".into();
                blobs.push(BlobEvidence {
                    sha256: sha.clone(),
                    bytes,
                    owners,
                    observed_hardlinks: links,
                    accounted_hardlinks: accounted,
                    was_legacy_hold,
                    release_hold: entry.releases_legacy_hold || !was_legacy_hold,
                });
                Ok(())
            })();
            if let Err(error) = result {
                entry.action = "held_error".into();
                entry.reasons.push(format!("{error:#}"));
                report.errors.push(format!("{sha}: {error:#}"));
            }
            report.entries.push(entry);
        }
        if options.apply && !blobs.is_empty() {
            for project in &evidence {
                if self
                    .inspect_legacy_project(&project.database)?
                    .catalog_sha256
                    != project.catalog_sha256
                {
                    bail!("legacy project changed before import publication");
                }
            }
            let descriptor = Descriptor {
                schema_version: 1,
                import_id: uuid::Uuid::new_v4().to_string(),
                store_root: self.root.clone(),
                complete_project_inventory: options.complete_project_inventory,
                projects: evidence,
                blobs,
                created_at: now_rfc3339(),
            };
            self.publish_legacy_import(&descriptor)?;
            report.import_id = Some(descriptor.import_id.clone());
            for entry in &mut report.entries {
                if entry.action == "import" {
                    entry.action = "imported".into();
                }
            }
        }
        report.complete = report.errors.is_empty()
            && report.unindexed_paths.is_empty()
            && report
                .entries
                .iter()
                .all(|entry| entry.releases_legacy_hold);
        Ok(report)
    }
}
