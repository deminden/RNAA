//! Read-only snapshots support pre-generation project catalogs without migrations.
use super::*;
use crate::{ArtifactKind, SharedBlobRecord};
use rusqlite::OpenFlags;
use std::str::FromStr;

#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct ProjectEvidence {
    pub database: PathBuf,
    pub project_id: String,
    pub catalog_sha256: String,
    pub sources: Vec<Source>,
}
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Source {
    pub artifact: ArtifactRecord,
    pub generation_id: Option<String>,
    pub physical_path: PathBuf,
}
pub(super) fn table(conn: &Connection, name: &str) -> Result<bool> {
    Ok(conn.query_row(
        "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name=?1)",
        [name],
        |row| row.get(0),
    )?)
}
fn artifacts(conn: &Connection) -> Result<Vec<ArtifactRecord>> {
    let mut query=conn.prepare("SELECT project_id,run_accession,kind,path,blob_id,shared_path,checksum_type,checksum,bytes,created_at FROM artifacts ORDER BY path")?;
    Ok(query
        .query_map([], |row| {
            let kind: String = row.get(2)?;
            let kind = ArtifactKind::from_str(&kind).map_err(|error| {
                rusqlite::Error::FromSqlConversionFailure(
                    2,
                    rusqlite::types::Type::Text,
                    error.into(),
                )
            })?;
            Ok(ArtifactRecord {
                project_id: row.get(0)?,
                run_accession: row.get(1)?,
                kind,
                path: row.get(3)?,
                blob_id: row.get(4)?,
                shared_path: row.get(5)?,
                checksum_type: row.get(6)?,
                checksum: row.get(7)?,
                bytes: read_u64(row, 8)?,
                created_at: row.get(9)?,
            })
        })?
        .collect::<std::result::Result<_, _>>()?)
}
#[derive(Serialize)]
struct GenerationEvidence {
    id: String,
    inventory: String,
    registry: String,
    retirement: Option<(String, String, String)>,
}
impl SharedStore {
    pub(super) fn inspect_legacy_project(&self, path: &Path) -> Result<ProjectEvidence> {
        if !path.symlink_metadata()?.is_file() {
            bail!("legacy project database must be a regular file");
        }
        let mut conn = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
        conn.busy_timeout(Duration::from_secs(5))?;
        let conn = conn.transaction()?;
        let schema = ownership::project_schema(&conn)?;
        let ids = {
            let mut query = conn.prepare("SELECT project_id FROM projects")?;
            query
                .query_map([], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<Vec<_>, _>>()?
        };
        if ids.len() != 1 {
            bail!("legacy project catalog must contain exactly one project");
        }
        let id = &ids[0];
        let root = path.parent().context("legacy project has no root")?;
        let active = artifacts(&conn)?;
        let blobs = {
            let mut query=conn.prepare("SELECT blob_id,storage_path,checksum_type,checksum,bytes,created_at FROM shared_blobs ORDER BY blob_id")?;
            query
                .query_map([], |row| {
                    Ok(SharedBlobRecord {
                        blob_id: row.get(0)?,
                        storage_path: row.get(1)?,
                        checksum_type: row.get(2)?,
                        checksum: row.get(3)?,
                        bytes: read_u64(row, 4)?,
                        created_at: row.get(5)?,
                    })
                })?
                .collect::<std::result::Result<Vec<_>, _>>()?
        };
        let refs = {
            let mut query=conn.prepare("SELECT blob_id,project_id,artifact_path FROM shared_blob_refs ORDER BY blob_id,project_id,artifact_path")?;
            query
                .query_map([], |row| {
                    Ok((
                        row.get::<_, String>(0)?,
                        row.get::<_, String>(1)?,
                        row.get::<_, String>(2)?,
                    ))
                })?
                .collect::<std::result::Result<Vec<_>, _>>()?
        };
        let mut generations = Vec::new();
        if table(&conn, "artifact_generations")? {
            let rows = {
                let mut query=conn.prepare("SELECT generation_id,inventory_json,artifacts_json FROM artifact_generations WHERE project_id=?1 AND state='COMMITTED' ORDER BY generation_id")?;
                query
                    .query_map([id], |row| {
                        Ok((
                            row.get::<_, String>(0)?,
                            row.get::<_, String>(1)?,
                            row.get::<_, String>(2)?,
                        ))
                    })?
                    .collect::<std::result::Result<Vec<_>, _>>()?
            };
            for (generation, inventory, registry) in rows {
                let retirement = if table(&conn, "generation_retirements")? {
                    conn.query_row("SELECT state,inventory_sha256,artifacts_sha256 FROM generation_retirements WHERE generation_id=?1",[&generation],|row|Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,String>(2)?))).optional()?
                } else {
                    None
                };
                generations.push(GenerationEvidence {
                    id: generation,
                    inventory,
                    registry,
                    retirement,
                });
            }
        }
        let fingerprint = digest(serde_json::to_vec(&(
            &schema,
            &ids,
            &active,
            &blobs,
            &refs,
            &generations,
        ))?);
        let mut sources = BTreeMap::new();
        for generation in &generations {
            uuid::Uuid::parse_str(&generation.id)?;
            let seal: serde_json::Value = serde_json::from_str(&generation.inventory)?;
            let payload = root
                .join(".rnaa_generations")
                .join(&generation.id)
                .join("payload");
            if seal["schema_version"] != 1
                || seal["head"]["generation_id"] != generation.id
                || seal["head"]["payload_path"] != payload.display().to_string()
            {
                bail!("legacy shared generation has an unsupported or inconsistent seal");
            }
            let registry: Vec<ArtifactRecord> = serde_json::from_str(&generation.registry)?;
            let retired = generation.retirement.as_ref().map(|row| row.0.as_str());
            if let Some((_, inventory, registry)) = &generation.retirement
                && (digest(&generation.inventory) != *inventory
                    || digest(&generation.registry) != *registry)
            {
                bail!("legacy shared retirement metadata differs from its original generation");
            }
            if !matches!(retired, None | Some("TRASHED" | "PURGED")) {
                bail!("finish pending project generation maintenance before shared import");
            }
            for record in registry {
                if !self.legacy_record_belongs(&record)? {
                    continue;
                }
                let file = Path::new(&record.path);
                let relative = file
                    .strip_prefix(&payload)
                    .context("shared generation artifact escapes its payload")?;
                let matches = seal["files"]
                    .as_array()
                    .context("legacy shared seal has no inventory")?
                    .iter()
                    .any(|entry| {
                        entry["path"] == record.path
                            && entry["checksum_type"] == "sha256"
                            && entry["checksum"] == record.checksum
                            && entry["bytes"].as_u64() == Some(record.bytes)
                    });
                if record.project_id != *id || !matches {
                    bail!("shared artifact differs from its sealed generation");
                }
                if retired == Some("PURGED") {
                    continue;
                }
                let physical = if retired == Some("TRASHED") {
                    root.join("trash/generations")
                        .join(&generation.id)
                        .join("payload")
                        .join(relative)
                } else {
                    file.into()
                };
                let witness = if retired == Some("TRASHED") {
                    root.join("trash/generations")
                        .join(&generation.id)
                        .join("seal.json")
                } else {
                    payload.parent().unwrap().join("seal.json")
                };
                if !witness.symlink_metadata()?.is_file()
                    || serde_json::from_slice::<serde_json::Value>(&fs::read(witness)?)? != seal
                {
                    bail!("shared generation seal witness differs from its catalog");
                }
                sources.insert(
                    (record.blob_id.clone(), record.path.clone()),
                    Source {
                        artifact: record,
                        generation_id: Some(generation.id.clone()),
                        physical_path: physical,
                    },
                );
            }
        }
        for record in &active {
            if !self.legacy_record_belongs(record)? {
                continue;
            }
            if record.project_id != *id {
                bail!("shared artifact project identity differs");
            }
            let key = (record.blob_id.clone(), record.path.clone());
            if let Some(source) = sources.get(&key) {
                if serde_json::to_value(&source.artifact)? != serde_json::to_value(record)? {
                    bail!("active shared view differs from its historical registration");
                }
            } else {
                if Path::new(&record.path).starts_with(root.join(".rnaa_generations")) {
                    bail!("active shared managed view has no committed owner");
                }
                sources.insert(
                    key,
                    Source {
                        artifact: record.clone(),
                        generation_id: None,
                        physical_path: PathBuf::from(&record.path),
                    },
                );
            }
        }
        for (blob, project, artifact) in &refs {
            if project != id {
                bail!("shared project reference has a foreign project identity");
            }
            if sources.contains_key(&(Some(blob.clone()), artifact.clone())) {
                continue;
            }
            if blobs.iter().any(|item| {
                item.blob_id == *blob && Path::new(&item.storage_path).starts_with(self.root())
            }) {
                bail!("legacy shared reference has no matching active or retained artifact");
            }
        }
        for source in sources.values() {
            let record = &source.artifact;
            let blob = blobs
                .iter()
                .find(|blob| Some(&blob.blob_id) == record.blob_id.as_ref())
                .context("legacy shared artifact has no local blob witness")?;
            if blob.checksum_type != "sha256"
                || blob.checksum != record.checksum
                || blob.bytes != record.bytes
                || canonical(Path::new(&blob.storage_path))? != self.blob_path(&record.checksum)?
            {
                bail!("legacy shared blob descriptor differs from its artifact");
            }
            if !Path::new(&record.path).is_absolute() {
                bail!("legacy shared artifact paths must be absolute");
            }
            if source.physical_path.exists() {
                self.verify_file(&source.physical_path, &record.checksum, record.bytes)?;
            } else if source.generation_id.is_some() {
                bail!("retained shared generation artifact is missing");
            }
        }
        Ok(ProjectEvidence {
            database: path.into(),
            project_id: id.clone(),
            catalog_sha256: fingerprint,
            sources: sources.into_values().collect(),
        })
    }
    fn legacy_record_belongs(&self, record: &ArtifactRecord) -> Result<bool> {
        let Some(path) = record.shared_path.as_ref() else {
            if record.blob_id.is_some() {
                bail!("legacy shared artifact has no shared path");
            }
            return Ok(false);
        };
        if !Path::new(path).is_absolute() {
            bail!("legacy shared paths must be absolute");
        }
        let path = canonical(Path::new(path))?;
        if !path.starts_with(self.root()) {
            return Ok(false);
        }
        files::validate_hash(&record.checksum)?;
        if record.checksum_type != "sha256"
            || record.blob_id.as_deref() != Some(format!("sha256:{}", record.checksum).as_str())
            || path != self.blob_path(&record.checksum)?
        {
            bail!("legacy shared artifact identity/path is inconsistent");
        }
        Ok(true)
    }
}
