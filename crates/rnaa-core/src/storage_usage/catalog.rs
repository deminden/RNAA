use super::*;
use anyhow::bail;
use rusqlite::{Connection, OpenFlags};
pub(super) struct Catalog {
    project: PathBuf,
    shared: Option<PathBuf>,
    generations: BTreeMap<String, String>,
    retired: BTreeSet<String>,
    file_trash: BTreeSet<String>,
    unpublished_trash: BTreeSet<String>,
    artifacts: BTreeSet<PathBuf>,
    blobs: BTreeMap<PathBuf, String>,
    project_epoch: i64,
    pub active: bool,
    shared_connection: Option<Connection>,
    shared_version: Option<i64>,
}
fn epoch(db: &Database) -> Result<(i64, bool)> {
    db.read(|conn|Ok(conn.query_row("SELECT (SELECT coalesce(max(id),0) FROM events WHERE project_id=?1),EXISTS(SELECT 1 FROM job_attempts WHERE project_id=?1 AND state='RUNNING')",[db.project_id()],|row|Ok((row.get(0)?,row.get(1)?)))?))
}
impl Catalog {
    pub fn capture(db: &Database, project: &Path, shared: Option<&Path>) -> Result<Self> {
        let (project_epoch, active) = epoch(db)?;
        let generations = db.read(|conn| {
            let mut query = conn.prepare(
                "SELECT generation_id,state FROM artifact_generations WHERE project_id=?1",
            )?;
            Ok(query
                .query_map([db.project_id()], |row| {
                    Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
                })?
                .collect::<std::result::Result<BTreeMap<_, _>, _>>()?)
        })?;
        let retired = db.read(|conn| {
            let mut query = conn
                .prepare("SELECT generation_id FROM generation_retirements WHERE project_id=?1")?;
            Ok(query
                .query_map([db.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<BTreeSet<_>, _>>()?)
        })?;
        let file_trash = db.read(|conn| {
            let mut query = conn.prepare("SELECT trash_id FROM file_trash WHERE project_id=?1")?;
            Ok(query
                .query_map([db.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<BTreeSet<_>, _>>()?)
        })?;
        let unpublished_trash = db.read(|conn| {
            let mut query = conn.prepare(
                "SELECT generation_id FROM unpublished_generation_trash WHERE project_id=?1",
            )?;
            Ok(query
                .query_map([db.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<BTreeSet<_>, _>>()?)
        })?;
        let artifacts = db
            .list_artifacts()?
            .into_iter()
            .map(|record| crate::generations::canonical_missing_path(Path::new(&record.path)))
            .collect::<Result<BTreeSet<_>>>()?;
        let mut value = Self {
            project: project.into(),
            shared: shared.map(PathBuf::from),
            generations,
            retired,
            file_trash,
            unpublished_trash,
            artifacts,
            blobs: BTreeMap::new(),
            project_epoch,
            active,
            shared_connection: None,
            shared_version: None,
        };
        if let Some(shared) = shared {
            let path = shared.join("catalog.sqlite3");
            if path.exists() {
                if !path.symlink_metadata()?.is_file() {
                    bail!("shared accounting catalog must be a regular file");
                }
                let conn = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
                conn.busy_timeout(std::time::Duration::from_secs(5))?;
                let schema: i64 =
                    conn.pragma_query_value(None, "user_version", |row| row.get(0))?;
                if !matches!(schema, 1 | 2) {
                    bail!("unsupported shared accounting catalog schema {schema}");
                }
                value.shared_version =
                    Some(conn.pragma_query_value(None, "data_version", |row| row.get(0))?);
                {
                    let mut query = conn.prepare("SELECT sha256 FROM blobs")?;
                    let hashes = query
                        .query_map([], |row| row.get::<_, String>(0))?
                        .collect::<std::result::Result<Vec<_>, _>>()?;
                    for hash in hashes {
                        if hash.len() != 64
                            || !hash
                                .bytes()
                                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
                        {
                            bail!("invalid shared accounting blob identity");
                        }
                        value.blobs.insert(
                            shared.join(format!(
                                "blobs/sha256/{}/{}/{hash}",
                                &hash[..2],
                                &hash[2..4]
                            )),
                            "shared_blob".into(),
                        );
                        value.blobs.insert(
                            shared.join(format!("trash/{hash}.blob")),
                            "shared_blob_trash".into(),
                        );
                    }
                }
                value.shared_connection = Some(conn);
            }
        }
        Ok(value)
    }
    pub fn classify(&self, path: &Path) -> String {
        if let Some(category) = self.blobs.get(path) {
            return category.clone();
        }
        if let Some(shared) = &self.shared
            && let Ok(relative) = path.strip_prefix(shared)
            && (relative.starts_with("blobs")
                || relative.starts_with("trash")
                    && relative.extension().is_some_and(|ext| ext == "blob"))
        {
            return "unindexed_shared".into();
        }
        if let Ok(relative) = path.strip_prefix(&self.project) {
            let parts = relative
                .components()
                .map(|part| part.as_os_str().to_string_lossy().into_owned())
                .collect::<Vec<_>>();
            if parts
                .first()
                .is_some_and(|part| part == ".rnaa_generations")
            {
                return match parts
                    .get(1)
                    .and_then(|id| self.generations.get(id))
                    .map(String::as_str)
                {
                    Some("COMMITTED") => "committed_generation",
                    Some("ABANDONED") => "abandoned_generation",
                    Some("BUILDING" | "SEALED") => "unpublished_generation",
                    _ => "unindexed_generation",
                }
                .into();
            }
            if parts.first().is_some_and(|part| part == "trash") {
                return if parts.get(1).is_some_and(|part| part == "generations")
                    && parts.get(2).is_some_and(|id| self.retired.contains(id))
                {
                    "generation_trash"
                } else if parts.get(1).is_some_and(|part| part == "files")
                    && parts.get(2).is_some_and(|id| self.file_trash.contains(id))
                {
                    "input_trash"
                } else if parts.get(1).is_some_and(|part| part == "unpublished")
                    && parts
                        .get(2)
                        .is_some_and(|id| self.unpublished_trash.contains(id))
                {
                    "unpublished_trash"
                } else {
                    "unindexed_trash"
                }
                .into();
            }
            if self.artifacts.contains(path) {
                return "legacy_registered".into();
            }
            if parts.first().is_some_and(|part| part == ".rnaa_transfers") {
                return "transfer_cache".into();
            }
            if parts.first().is_some_and(|part| part == ".rnaa_support") {
                return "support_and_temporary".into();
            }
            if parts
                .first()
                .is_some_and(|part| matches!(part.as_str(), "metadata" | "logs" | ".rnaa-locks"))
                || parts.len() == 1 && parts[0].starts_with("state.sqlite")
                || relative == Path::new("rnaa.toml")
            {
                return "project_metadata".into();
            }
            return "legacy_or_unindexed".into();
        }
        "shared_metadata_or_unindexed".into()
    }
    pub fn finish(&mut self, db: &Database) -> Result<(bool, Option<bool>, bool)> {
        let (project, active) = epoch(db)?;
        let shared = self
            .shared_connection
            .as_ref()
            .map(|conn| -> Result<bool> {
                let version: i64 =
                    conn.pragma_query_value(None, "data_version", |row| row.get(0))?;
                Ok(Some(version) != self.shared_version)
            })
            .transpose()?;
        Ok((project != self.project_epoch, shared, active))
    }
}
