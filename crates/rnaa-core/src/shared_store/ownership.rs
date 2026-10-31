use super::*;
use crate::ArtifactRecord;
use rusqlite::OpenFlags;
// 0022 adds only export publication bookkeeping and no shared-blob references.
const LAST_SUPPORTED_PROJECT_MIGRATION: u32 = 22;
fn table_exists(conn: &Connection, name: &str) -> Result<bool> {
    Ok(conn.query_row(
        "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name=?1)",
        [name],
        |row| row.get(0),
    )?)
}
/// Ownership inspection must opt in when a future project schema changes its
/// reference contracts. Missing/unknown migration history is not absence proof.
pub(super) fn project_schema(conn: &Connection) -> Result<Vec<String>> {
    let mut query = conn.prepare("SELECT version FROM schema_migrations ORDER BY version")?;
    let versions = query
        .query_map([], |row| row.get::<_, String>(0))?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    if !versions.iter().any(|version| version == "0001_init")
        || !versions
            .iter()
            .any(|version| version == "0002_shared_store")
        || versions.iter().any(|version| {
            !Database::supports_migration(version)
                || version
                    .split('_')
                    .next()
                    .and_then(|value| value.parse::<u32>().ok())
                    .is_none_or(|number| number > LAST_SUPPORTED_PROJECT_MIGRATION)
        })
    {
        bail!("shared ownership cannot inspect this project schema; keep its ownership protected");
    }
    Ok(versions)
}
fn shared_path_matches(store: &SharedStore, path: Option<&str>, sha: &str) -> Result<bool> {
    let Some(path) = path else {
        return Ok(false);
    };
    Ok(Path::new(path).is_absolute()
        && crate::generations::canonical_missing_path(Path::new(path))? == store.blob_path(sha)?)
}

#[derive(Clone)]
pub(super) struct Owner {
    pub sha256: String,
    pub project_db: String,
    pub project_id: String,
    pub artifact_path: String,
    pub state: String,
    pub generation_id: Option<String>,
}
pub(super) enum Proof {
    Present(Option<String>),
    Absent,
}

pub(super) fn source_generation(db: &Database, source: &Path) -> Result<Option<String>> {
    let source = fs::canonicalize(source)?;
    let root = fs::canonicalize(
        db.path()
            .parent()
            .context("shared owner has no project root")?,
    )?;
    let managed = root.join(".rnaa_generations");
    let Ok(relative) = source.strip_prefix(managed) else {
        return Ok(None);
    };
    let generation = relative
        .components()
        .next()
        .context("managed shared source has no generation")?
        .as_os_str()
        .to_str()
        .context("invalid managed generation identity")?
        .to_string();
    let payload:String=db.read(|connection| {
        Ok(connection.query_row("SELECT payload_path FROM artifact_generations WHERE project_id=?1 AND generation_id=?2",
            params![db.project_id(),generation],|row|row.get(0))?)
    }).context("managed shared source has no project generation")?;
    if !source.starts_with(payload) {
        bail!("shared source differs from its project generation");
    }
    Ok(Some(generation))
}

/// Read a consistent project snapshot without creating/migrating/recovering it.
/// Missing, substituted, busy or corrupt project catalogs are errors, never proof
/// that their last known ownership can be discarded.
pub(super) fn prove(store: &SharedStore, owner: &Owner) -> Result<Proof> {
    let expected_bytes: u64 = store.connection()?.query_row(
        "SELECT bytes FROM blobs WHERE sha256=?1",
        [&owner.sha256],
        |row| read_u64(row, 0),
    )?;
    let mut connection =
        Connection::open_with_flags(&owner.project_db, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    connection.busy_timeout(Duration::from_secs(5))?;
    let transaction = connection.transaction()?;
    project_schema(&transaction)?;
    let project: String =
        transaction.query_row("SELECT project_id FROM projects LIMIT 1", [], |row| {
            row.get(0)
        })?;
    if project != owner.project_id {
        bail!("shared owner project identity changed");
    }
    let matches = |record: &ArtifactRecord| -> Result<bool> {
        Ok(record.project_id == owner.project_id
            && record.path == owner.artifact_path
            && record.blob_id.as_deref() == Some(format!("sha256:{}", owner.sha256).as_str())
            && record.checksum_type == "sha256"
            && record.checksum == owner.sha256
            && record.bytes == expected_bytes
            && shared_path_matches(store, record.shared_path.as_deref(), &owner.sha256)?)
    };
    if let Some(generation) = &owner.generation_id {
        let (state,artifacts):(String,Option<String>)=transaction.query_row(
            "SELECT state,artifacts_json FROM artifact_generations WHERE project_id=?1 AND generation_id=?2",
            params![owner.project_id,generation],|row|Ok((row.get(0)?,row.get(1)?)))
            .context("shared generation ownership witness is missing")?;
        let retirement: Option<(String, String, String)> = if table_exists(
            &transaction,
            "generation_retirements",
        )? {
            transaction.query_row(
            "SELECT state,inventory_sha256,artifacts_sha256 FROM generation_retirements WHERE generation_id=?1 AND project_id=?2",
            params![generation,owner.project_id],|row|Ok((row.get(0)?,row.get(1)?,row.get(2)?))).optional()?
        } else {
            None
        };
        if let Some((retired, inventory_hash, artifacts_hash)) = retirement
            && retired == "PURGED"
        {
            let active:bool=transaction.query_row("SELECT EXISTS(SELECT 1 FROM generation_heads WHERE generation_id=?1) OR EXISTS(SELECT 1 FROM reference_selection WHERE generation_id=?1) OR EXISTS(SELECT 1 FROM generation_keeps WHERE generation_id=?1) OR EXISTS(SELECT 1 FROM artifacts WHERE project_id=?2 AND path=?3)",params![generation,owner.project_id,owner.artifact_path],|row|row.get(0))?;
            if active {
                bail!("purged generation still has active project references");
            }
            use sha2::{Digest, Sha256};
            let inventory: String = transaction.query_row(
                "SELECT inventory_json FROM artifact_generations WHERE generation_id=?1",
                [generation],
                |row| row.get(0),
            )?;
            if hex::encode(Sha256::digest(inventory.as_bytes())) != inventory_hash
                || artifacts.as_ref().is_none_or(|json| {
                    hex::encode(Sha256::digest(json.as_bytes())) != artifacts_hash
                })
            {
                bail!("shared ownership retirement receipt differs from its generation");
            }
            return Ok(Proof::Absent);
        }
        if state == "COMMITTED" {
            for record in serde_json::from_str::<Vec<ArtifactRecord>>(
                &artifacts.context("committed shared generation has no artifact registry")?,
            )? {
                if matches(&record)? {
                    return Ok(Proof::Present(Some(generation.clone())));
                }
            }
        }
    } else if table_exists(&transaction, "artifact_generations")? {
        let mut query=transaction.prepare("SELECT generation_id,artifacts_json FROM artifact_generations WHERE project_id=?1 AND state='COMMITTED'")?;
        let rows = query.query_map([&owner.project_id], |row| {
            Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
        })?;
        for row in rows {
            let (generation, artifacts) = row?;
            for record in serde_json::from_str::<Vec<ArtifactRecord>>(&artifacts)? {
                if matches(&record)? {
                    return Ok(Proof::Present(Some(generation)));
                }
            }
        }
    }
    struct ActiveBlob {
        id: Option<String>,
        kind: String,
        hash: String,
        path: Option<String>,
        bytes: u64,
    }
    let active:Option<ActiveBlob>=transaction.query_row(
        "SELECT blob_id,checksum_type,checksum,shared_path,bytes FROM artifacts WHERE project_id=?1 AND path=?2",
        params![owner.project_id,owner.artifact_path],|row|Ok(ActiveBlob{id:row.get(0)?,kind:row.get(1)?,hash:row.get(2)?,path:row.get(3)?,bytes:read_u64(row,4)?})).optional()?;
    if let Some(active) = active
        && active.id.as_deref() == Some(format!("sha256:{}", owner.sha256).as_str())
    {
        if active.kind != "sha256"
            || active.hash != owner.sha256
            || active.bytes != expected_bytes
            || !shared_path_matches(store, active.path.as_deref(), &owner.sha256)?
        {
            bail!("active shared ownership metadata differs from its blob");
        }
        return Ok(Proof::Present(owner.generation_id.clone()));
    }
    if owner.state == "COMMITTED" && owner.generation_id.is_some() {
        bail!("committed shared generation ownership changed without a retirement receipt");
    }
    Ok(Proof::Absent)
}

impl SharedStore {
    pub(super) fn owners(&self, sha256: &str) -> Result<Vec<Owner>> {
        let connection = self.connection()?;
        let mut query=connection.prepare("SELECT project_db,project_id,artifact_path,state,generation_id FROM owners WHERE sha256=?1 ORDER BY project_db,artifact_path")?;
        Ok(query
            .query_map([sha256], |row| {
                Ok(Owner {
                    sha256: sha256.into(),
                    project_db: row.get(0)?,
                    project_id: row.get(1)?,
                    artifact_path: row.get(2)?,
                    state: row.get(3)?,
                    generation_id: row.get(4)?,
                })
            })?
            .collect::<std::result::Result<Vec<_>, _>>()?)
    }
    pub(super) fn reconcile_owners(&self, sha256: &str, apply: bool) -> Result<(u64, Vec<String>)> {
        let mut present = 0;
        let mut errors = Vec::new();
        for owner in self.owners(sha256)? {
            let proof = prove(self, &owner);
            match proof {
                Ok(Proof::Present(generation)) => {
                    present += 1;
                    if apply {
                        self.connection()?.execute("UPDATE owners SET state='COMMITTED',generation_id=?1 WHERE sha256=?2 AND project_db=?3 AND project_id=?4 AND artifact_path=?5",
                        params![generation,sha256,owner.project_db,owner.project_id,owner.artifact_path])?;
                    }
                }
                Ok(Proof::Absent) => {
                    if apply {
                        self.connection()?.execute("DELETE FROM owners WHERE sha256=?1 AND project_db=?2 AND project_id=?3 AND artifact_path=?4",params![sha256,owner.project_db,owner.project_id,owner.artifact_path])?;
                    }
                }
                Err(error) => {
                    present += 1;
                    errors.push(format!("{}: {error:#}", owner.project_db));
                }
            }
        }
        Ok((present, errors))
    }
}
