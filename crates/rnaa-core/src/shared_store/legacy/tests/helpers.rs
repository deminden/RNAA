use super::*;
pub(super) struct Project {
    pub root: PathBuf,
    pub database: PathBuf,
    pub source: PathBuf,
    pub sha: String,
}
pub(super) fn legacy_project(
    root: &Path,
    store: &SharedStore,
    payload: &[u8],
    hardlink: bool,
) -> Result<Project> {
    fs::create_dir_all(root)?;
    let root = fs::canonicalize(root)?;
    let database = root.join("state.sqlite");
    let conn = Connection::open(&database)?;
    conn.execute_batch(include_str!("../../../../migrations/0001_init.sql"))?;
    conn.execute_batch(include_str!("../../../../migrations/0002_shared_store.sql"))?;
    conn.execute_batch(
        "INSERT INTO schema_migrations VALUES ('0001_init','old'),('0002_shared_store','old');",
    )?;
    let id = uuid::Uuid::new_v4().to_string();
    conn.execute("INSERT INTO projects(project_id,root_dir,created_at,config_json) VALUES (?1,?2,'old','{}')",params![id,root.display().to_string()])?;
    let source = root.join("abundance.h5");
    let sha = digest(payload);
    let blob = store.blob_path(&sha)?;
    store.checked_path(blob.parent().unwrap().strip_prefix(store.root())?, true)?;
    if !blob.exists() {
        fs::write(&blob, payload)?;
    }
    if hardlink {
        fs::hard_link(&blob, &source)?;
    } else {
        fs::copy(&blob, &source)?;
    }
    let blob_id = format!("sha256:{sha}");
    conn.execute("INSERT INTO shared_blobs(blob_id,storage_path,checksum_type,checksum,bytes,created_at) VALUES (?1,?2,'sha256',?3,?4,'old')",params![blob_id,blob.display().to_string(),sha,i64::try_from(payload.len())?])?;
    conn.execute("INSERT INTO artifacts(project_id,run_accession,kind,path,checksum_type,checksum,bytes,created_at,blob_id,shared_path) VALUES (?1,NULL,?2,?3,'sha256',?4,?5,'old',?6,?7)",params![id,crate::ArtifactKind::QuantAbundanceH5.as_str(),source.display().to_string(),sha,i64::try_from(payload.len())?,blob_id,blob.display().to_string()])?;
    conn.execute("INSERT INTO shared_blob_refs(blob_id,project_id,artifact_path,referenced_at) VALUES (?1,?2,?3,'old')",params![blob_id,id,source.display().to_string()])?;
    drop(conn);
    Ok(Project {
        root,
        database,
        source,
        sha,
    })
}
pub(super) fn remove_registration(project: &Project) -> Result<()> {
    let conn = Connection::open(&project.database)?;
    conn.execute("DELETE FROM shared_blob_refs", [])?;
    conn.execute("DELETE FROM artifacts", [])?;
    Ok(())
}
pub(super) fn opts(apply: bool, complete: bool) -> LegacyImportOptions {
    LegacyImportOptions {
        apply,
        complete_project_inventory: complete,
    }
}
pub(super) fn gc(apply: bool) -> GcOptions {
    GcOptions {
        apply,
        purge_expired: false,
        trash_days: 7,
    }
}
