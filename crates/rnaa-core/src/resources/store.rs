use super::{
    OwnerFileIdentity, PoolSnapshot, ReservationSnapshot, ResourceRecovery, ResourceRequest, files,
};
use anyhow::{Context, Result, bail};
use rusqlite::{Connection, OptionalExtension};
use std::{
    fs::{self, File},
    path::{Path, PathBuf},
    time::Duration,
};

#[derive(Debug, Clone)]
pub struct ResourcePool {
    pub(super) root: PathBuf,
}
pub(super) struct Catalog {
    pub connection: Connection,
    _guard: File,
}
pub(super) struct Reclaimed {
    path: PathBuf,
    file: File,
    identity: OwnerFileIdentity,
}
impl ResourcePool {
    pub fn open(root: &Path) -> Result<Self> {
        let mut builder = fs::DirBuilder::new();
        builder.recursive(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::DirBuilderExt;
            builder.mode(0o700);
        }
        builder.create(root)?;
        let pool = Self {
            root: fs::canonicalize(root)?,
        };
        let owners = pool.root.join("owners");
        if let Ok(metadata) = owners.symlink_metadata()
            && (!metadata.is_dir() || metadata.file_type().is_symlink())
        {
            bail!("resource owner directory must not be a symlink");
        }
        builder.create(&owners)?;
        pool.catalog()?;
        Ok(pool)
    }
    pub fn root(&self) -> &Path {
        &self.root
    }

    pub(super) fn catalog(&self) -> Result<Catalog> {
        let guard = files::exclusive(&self.root.join("catalog.lock"))?;
        let path = self.root.join("reservations.sqlite");
        if let Ok(metadata) = path.symlink_metadata()
            && (!metadata.is_file() || metadata.file_type().is_symlink())
        {
            bail!("resource catalog must be a regular file");
        }
        let mut connection = Connection::open(path)?;
        connection.busy_timeout(Duration::from_secs(15))?;
        connection.execute_batch(
            "PRAGMA journal_mode=WAL; PRAGMA synchronous=FULL; PRAGMA foreign_keys=ON;",
        )?;
        let transaction = connection.transaction()?;
        transaction.execute_batch("CREATE TABLE IF NOT EXISTS pool_schema(singleton INTEGER PRIMARY KEY CHECK(singleton=1), version INTEGER NOT NULL);")?;
        let version = transaction
            .query_row(
                "SELECT version FROM pool_schema WHERE singleton=1",
                [],
                |row| row.get::<_, i64>(0),
            )
            .optional()?;
        match version {
            None => {
                transaction
                    .execute("INSERT INTO pool_schema(singleton,version) VALUES(1,1)", [])?;
            }
            Some(1) => {}
            Some(version) => bail!("unsupported resource catalog schema {version}"),
        }
        transaction.execute_batch("CREATE TABLE IF NOT EXISTS reservations(reservation_id TEXT PRIMARY KEY, descriptor_json TEXT NOT NULL);")?;
        transaction.commit()?;
        Ok(Catalog {
            connection,
            _guard: guard,
        })
    }

    pub(super) fn owner_path(&self, id: &str) -> Result<PathBuf> {
        let parsed = uuid::Uuid::parse_str(id)
            .context("resource catalog has an invalid reservation identity")?;
        if parsed.to_string() != id {
            bail!("resource reservation identity is not canonical");
        }
        let owners = self.root.join("owners");
        let metadata = owners.symlink_metadata()?;
        if !metadata.is_dir() || metadata.file_type().is_symlink() {
            bail!("resource owner directory changed; ownership records are held");
        }
        Ok(owners.join(format!("{id}.lock")))
    }
    pub(super) fn normalize_request(
        &self,
        mut request: ResourceRequest,
    ) -> Result<ResourceRequest> {
        for disk in &mut request.disks {
            disk.path = files::existing_path(&disk.path)?;
        }
        Ok(request)
    }

    /// The caller holds the catalog lock and a SQLite transaction throughout.
    /// Never infer an owner's death from a clock or a process identifier.
    pub(super) fn live_rows(
        &self,
        connection: &Connection,
    ) -> Result<(Vec<ReservationSnapshot>, Vec<Reclaimed>)> {
        let mut statement = connection.prepare(
            "SELECT reservation_id,descriptor_json FROM reservations ORDER BY reservation_id",
        )?;
        let rows = statement
            .query_map([], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
            })?
            .collect::<rusqlite::Result<Vec<_>>>()?;
        drop(statement);
        let mut live = Vec::new();
        let mut reclaimed = Vec::new();
        for (id, json) in rows {
            let path = self.owner_path(&id)?;
            let row: ReservationSnapshot = serde_json::from_str(&json)
                .context("resource reservation descriptor is corrupt")?;
            if row.reservation_id != id {
                bail!("resource reservation descriptor identity differs from its catalog key");
            }
            match files::try_owner(&path, Some(row.owner_file_identity)).with_context(|| {
                format!(
                    "cannot establish ownership of resource reservation {id}; its record is held"
                )
            })? {
                Some(file) => {
                    connection
                        .execute("DELETE FROM reservations WHERE reservation_id=?1", [&id])?;
                    reclaimed.push(Reclaimed {
                        path,
                        file,
                        identity: row.owner_file_identity,
                    });
                }
                None => {
                    live.push(row);
                }
            }
        }
        Ok((live, reclaimed))
    }
    pub(super) fn clean_reclaimed(&self, reclaimed: Vec<Reclaimed>) {
        for row in reclaimed {
            drop(row.file);
            if let Err(error) = files::remove_owner(&row.path, row.identity) {
                tracing::warn!(path=%row.path.display(),%error,"retired resource lock file remains for maintenance");
            }
        }
    }

    pub fn snapshot(&self) -> Result<PoolSnapshot> {
        let mut catalog = self.catalog()?;
        let transaction = catalog.connection.transaction()?;
        let (live, reclaimed) = self.live_rows(&transaction)?;
        let snapshot = PoolSnapshot::from_rows(live)?;
        transaction.commit()?;
        self.clean_reclaimed(reclaimed);
        Ok(snapshot)
    }

    /// Reclaim dead owners and empty protocol lock files left before a catalog
    /// commit. Unknown names, nonempty files and still-locked orphans are held.
    pub fn recover(&self) -> Result<ResourceRecovery> {
        let mut catalog = self.catalog()?;
        let transaction = catalog.connection.transaction()?;
        let (live, reclaimed) = self.live_rows(&transaction)?;
        let _ = PoolSnapshot::from_rows(live.clone())?;
        let mut result = ResourceRecovery {
            abandoned_reservations: reclaimed.len(),
            orphan_lock_files: 0,
        };
        transaction.commit()?;
        self.clean_reclaimed(reclaimed);
        let ids = live
            .into_iter()
            .map(|row| row.reservation_id)
            .collect::<std::collections::BTreeSet<_>>();
        for entry in fs::read_dir(self.root.join("owners"))? {
            let entry = entry?;
            let name = entry.file_name();
            let Some(name) = name.to_str() else { continue };
            let Some(id) = name.strip_suffix(".lock") else {
                continue;
            };
            if !uuid::Uuid::parse_str(id).is_ok_and(|parsed| parsed.to_string() == id)
                || ids.contains(id)
            {
                continue;
            }
            let path = self.owner_path(id)?;
            let metadata = path.symlink_metadata()?;
            if !metadata.is_file() || metadata.file_type().is_symlink() || metadata.len() != 0 {
                continue;
            }
            if let Some(file) = files::try_owner(&path, None)? {
                let identity = files::identity(&file.metadata()?)?;
                drop(file);
                files::remove_owner(&path, identity)?;
                result.orphan_lock_files += 1;
            }
        }
        Ok(result)
    }
}
