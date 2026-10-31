use super::*;
pub(super) fn known_projects(conn: &Connection) -> Result<BTreeSet<(PathBuf, String)>> {
    let mut query=conn.prepare("SELECT project_db,project_id FROM owners UNION SELECT project_db,project_id FROM legacy_projects")?;
    Ok(query
        .query_map([], |row| {
            Ok((
                PathBuf::from(row.get::<_, String>(0)?),
                row.get::<_, String>(1)?,
            ))
        })?
        .collect::<std::result::Result<_, _>>()?)
}
impl SharedStore {
    fn import_witness(&self, id: &str) -> Result<PathBuf> {
        uuid::Uuid::parse_str(id)?;
        let path = self.checked_path(Path::new(&format!(".imports/{id}.json")), false)?;
        if path.symlink_metadata().is_ok_and(|meta| !meta.is_file()) {
            bail!("shared import witness must be a regular file");
        }
        Ok(path)
    }
    pub(super) fn publish_legacy_import(&self, descriptor: &Descriptor) -> Result<()> {
        for blob in &descriptor.blobs {
            let path = self.blob_path(&blob.sha256)?;
            self.verify_file(&path, &blob.sha256, blob.bytes)?;
            let physical = descriptor
                .projects
                .iter()
                .flat_map(|project| &project.sources)
                .filter(|source| source.artifact.checksum == blob.sha256)
                .map(|source| source.physical_path.clone())
                .collect::<Vec<_>>();
            if hardlink_census(&path, &physical)?
                != (blob.observed_hardlinks, blob.accounted_hardlinks)
            {
                bail!("legacy shared links changed before import publication");
            }
            File::open(&path)?.sync_all()?;
        }
        let json = serde_json::to_string(descriptor)?;
        crate::util::publish_file(&self.import_witness(&descriptor.import_id)?, |temporary| {
            fs::write(temporary, &json)?;
            Ok(())
        })?;
        let mut conn = self.connection()?;
        let tx = conn.transaction_with_behavior(TransactionBehavior::Immediate)?;
        if descriptor
            .blobs
            .iter()
            .any(|blob| blob.was_legacy_hold && blob.release_hold)
            && known_projects(&tx)?.iter().any(|(path, id)| {
                !descriptor
                    .projects
                    .iter()
                    .any(|project| project.database == *path && project.project_id == *id)
            })
        {
            bail!("known shared project inventory changed before hold release");
        }
        tx.execute("INSERT INTO legacy_imports(import_id,descriptor_json,descriptor_sha256,created_at) VALUES (?1,?2,?3,?4)",params![descriptor.import_id,json,digest(&json),descriptor.created_at])?;
        for project in &descriptor.projects {
            tx.execute("INSERT INTO legacy_projects(project_db,project_id,first_import_id) VALUES (?1,?2,?3) ON CONFLICT DO NOTHING",params![project.database.display().to_string(),project.project_id,descriptor.import_id])?;
        }
        for blob in &descriptor.blobs {
            tx.execute("INSERT INTO blobs(sha256,bytes,state,legacy_hold,created_at) VALUES (?1,?2,'LIVE',1,?3) ON CONFLICT DO NOTHING",params![blob.sha256,i64::try_from(blob.bytes)?,descriptor.created_at])?;
            let (bytes, state, legacy): (u64, String, bool) = tx.query_row(
                "SELECT bytes,state,legacy_hold FROM blobs WHERE sha256=?1",
                [&blob.sha256],
                |row| Ok((read_u64(row, 0)?, row.get(1)?, row.get(2)?)),
            )?;
            if bytes != blob.bytes || state != "LIVE" || legacy != blob.was_legacy_hold {
                bail!("shared blob catalog changed before legacy import publication");
            }
            for owner in &blob.owners {
                let previous:Option<Option<String>>=tx.query_row("SELECT generation_id FROM owners WHERE sha256=?1 AND project_db=?2 AND project_id=?3 AND artifact_path=?4",params![blob.sha256,owner.project_db.display().to_string(),owner.project_id,owner.artifact_path],|row|row.get(0)).optional()?;
                if previous
                    .flatten()
                    .is_some_and(|generation| Some(&generation) != owner.generation_id.as_ref())
                {
                    bail!("shared owner generation differs from imported evidence");
                }
                tx.execute("INSERT INTO owners(sha256,project_db,project_id,artifact_path,state,generation_id,created_at) VALUES (?1,?2,?3,?4,'COMMITTED',?5,?6) ON CONFLICT(sha256,project_db,project_id,artifact_path) DO UPDATE SET state='COMMITTED',generation_id=excluded.generation_id",params![blob.sha256,owner.project_db.display().to_string(),owner.project_id,owner.artifact_path,owner.generation_id,descriptor.created_at])?;
            }
            tx.execute("INSERT INTO legacy_blob_imports(sha256,import_id,released_hold) VALUES (?1,?2,?3) ON CONFLICT(sha256) DO UPDATE SET import_id=excluded.import_id,released_hold=excluded.released_hold",params![blob.sha256,descriptor.import_id,blob.release_hold])?;
            tx.execute(
                "UPDATE blobs SET legacy_import_id=?1 WHERE sha256=?2",
                params![descriptor.import_id, blob.sha256],
            )?;
            if blob.release_hold {
                tx.execute(
                    "UPDATE blobs SET legacy_hold=0 WHERE sha256=?1",
                    [&blob.sha256],
                )?;
            }
        }
        tx.commit()?;
        Ok(())
    }
    /// Immutable receipt integrity and still-present historical owners both gate
    /// collection. A missing owner row cannot erase a migrated project's claim.
    pub(in crate::shared_store) fn verify_legacy_import(&self, sha: &str) -> Result<()> {
        let stored:Option<(String,String,String,bool)>=self.connection()?.query_row("SELECT i.import_id,i.descriptor_json,i.descriptor_sha256,b.released_hold FROM legacy_blob_imports b JOIN legacy_imports i ON i.import_id=b.import_id WHERE b.sha256=?1",[sha],|row|Ok((row.get(0)?,row.get(1)?,row.get(2)?,row.get(3)?))).optional()?;
        let marker: Option<String> = self
            .connection()?
            .query_row(
                "SELECT legacy_import_id FROM blobs WHERE sha256=?1",
                [sha],
                |row| row.get::<_, Option<String>>(0),
            )
            .optional()?
            .flatten();
        let Some((id, json, hash, released)) = stored else {
            if marker.is_some() {
                bail!("legacy shared import catalog index is missing");
            }
            return Ok(());
        };
        if marker.as_deref() != Some(id.as_str()) {
            bail!("legacy shared import marker differs from its catalog index");
        }
        let descriptor: Descriptor = serde_json::from_str(&json)?;
        if descriptor.schema_version != 1
            || descriptor.import_id != id
            || descriptor.store_root != self.root
            || digest(&json) != hash
            || digest(fs::read(self.import_witness(&id)?)?) != hash
        {
            bail!("shared legacy import receipt or filesystem witness changed");
        }
        let matches = descriptor
            .blobs
            .iter()
            .filter(|blob| blob.sha256 == sha)
            .collect::<Vec<_>>();
        if matches.len() != 1 {
            bail!("shared legacy import omits or duplicates the blob");
        }
        let blob = matches[0];
        let (bytes, legacy): (u64, bool) = self.connection()?.query_row(
            "SELECT bytes,legacy_hold FROM blobs WHERE sha256=?1",
            [sha],
            |row| Ok((read_u64(row, 0)?, row.get(1)?)),
        )?;
        if blob.bytes != bytes
            || blob.release_hold != released
            || (!legacy && !released)
            || (blob.was_legacy_hold
                && released
                && (!descriptor.complete_project_inventory
                    || blob.observed_hardlinks != Some(blob.accounted_hardlinks)
                    || blob.accounted_hardlinks == 0))
        {
            bail!("shared legacy hold release differs from its migration proof");
        }
        let mut identities = BTreeSet::new();
        for owner in &blob.owners {
            if !identities.insert((&owner.project_db, &owner.project_id, &owner.artifact_path)) {
                bail!("duplicate imported shared owner");
            }
            let project = descriptor
                .projects
                .iter()
                .find(|project| {
                    project.database == owner.project_db && project.project_id == owner.project_id
                })
                .context("imported owner has no project inventory witness")?;
            if !project.sources.iter().any(|source| {
                source.artifact.checksum == sha
                    && source.artifact.bytes == bytes
                    && source.artifact.path == owner.artifact_path
                    && source.generation_id == owner.generation_id
            }) {
                bail!("imported owner differs from its project snapshot");
            }
            if matches!(
                ownership::prove(self, &owner.owner(sha))?,
                ownership::Proof::Present(_)
            ) {
                let recorded:bool=self.connection()?.query_row("SELECT EXISTS(SELECT 1 FROM owners WHERE sha256=?1 AND project_db=?2 AND project_id=?3 AND artifact_path=?4)",params![sha,owner.project_db.display().to_string(),owner.project_id,owner.artifact_path],|row|row.get(0))?;
                if !recorded {
                    bail!(
                        "migrated shared ownership row is missing while its project still references the blob"
                    );
                }
            }
        }
        Ok(())
    }
}
