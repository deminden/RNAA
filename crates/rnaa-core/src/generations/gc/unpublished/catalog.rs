use super::*;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct GenerationRow {
    pub id: String,
    pub project: String,
    pub group: String,
    pub stage: String,
    pub subject: String,
    pub attempt: String,
    pub payload: PathBuf,
    pub request: Option<String>,
    pub state: String,
    pub inventory: Option<String>,
    pub artifacts: Option<String>,
    pub created_at: String,
    pub committed_at: Option<String>,
}
impl GenerationRow {
    pub fn unpublished(db: &Database) -> Result<BTreeMap<String, Self>> {
        db.read(|conn| {
            let mut query = conn.prepare("SELECT generation_id,project_id,output_group,stage,subject,attempt_id,payload_path,request_sha256,state,inventory_json,artifacts_json,created_at,committed_at FROM artifact_generations WHERE project_id=?1 AND state<>'COMMITTED' ORDER BY generation_id")?;
            Ok(query.query_map([db.project_id()], |row| Ok(Self {
                id:row.get(0)?,project:row.get(1)?,group:row.get(2)?,stage:row.get(3)?,subject:row.get(4)?,attempt:row.get(5)?,payload:PathBuf::from(row.get::<_,String>(6)?),request:row.get(7)?,state:row.get(8)?,inventory:row.get(9)?,artifacts:row.get(10)?,created_at:row.get(11)?,committed_at:row.get(12)?,
            }))?.map(|row|row.map(|item|(item.id.clone(),item))).collect::<std::result::Result<_,_>>()?)
        })
    }
    pub fn verify_unpublished(&self, db: &Database) -> Result<()> {
        if self.project != db.project_id()
            || !matches!(self.state.as_str(), "BUILDING" | "SEALED" | "ABANDONED")
            || self.committed_at.is_some()
            || self.payload != tree::container(db, &self.id, false, false)?.join("payload")
        {
            bail!("generation is not a tracked unpublished container");
        }
        db.read(|conn| {
            let (project,stage,state):(String,String,String)=conn.query_row("SELECT project_id,stage,state FROM job_attempts WHERE attempt_id=?1",[&self.attempt],|row|Ok((row.get(0)?,row.get(1)?,row.get(2)?)))?;
            if project != self.project || stage != self.stage || state == "RUNNING" {
                bail!("owning attempt is active or does not match the unpublished generation");
            }
            let referenced:bool=conn.query_row(
                "SELECT EXISTS(SELECT 1 FROM generation_heads WHERE generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_keeps WHERE generation_id=?1)
                 OR EXISTS(SELECT 1 FROM reference_selection WHERE generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_retirements WHERE generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_files WHERE generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_lineage_imports WHERE generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_input_bindings WHERE generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_binding_pruning WHERE producer_generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_pruned_files WHERE consumer_generation_id=?1 OR producer_generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_dependencies WHERE consumer_generation_id=?1 OR producer_generation_id=?1)
                 OR EXISTS(SELECT 1 FROM generation_import_dependencies WHERE consumer_generation_id=?1 OR producer_generation_id=?1)",
                [&self.id],|row|row.get(0))?;
            if referenced {bail!("unpublished generation has published catalog references");}
            Ok(())
        })
    }
}

#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Descriptor {
    pub schema_version: u32,
    pub operation_id: String,
    pub generation: GenerationRow,
    pub original_state: String,
    pub tree: Vec<tree::TreeEntry>,
    pub trashed_at_ms: i64,
    pub purge_after_ms: i64,
}
#[derive(Clone)]
pub(super) struct Receipt {
    pub descriptor: Descriptor,
    pub state: String,
}
pub(super) fn witness_path(db: &Database, operation: &str, create: bool) -> Result<PathBuf> {
    uuid::Uuid::parse_str(operation)?;
    let path = tree::checked_directory(db, &[".rnaa_support", "unpublished_trash"], create)?
        .join(format!("{operation}.json"));
    if path.symlink_metadata().is_ok_and(|meta| !meta.is_file()) {
        bail!("unpublished trash witness must be a regular file");
    }
    Ok(path)
}
impl Receipt {
    pub fn all(db: &Database) -> Result<BTreeMap<String, Self>> {
        let stored=db.read(|conn| {
            let mut query=conn.prepare("SELECT generation_id,project_id,state,descriptor_json,descriptor_sha256,trashed_at_ms,purge_after_ms FROM unpublished_generation_trash ORDER BY generation_id")?;
            Ok(query.query_map([],|row|Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,String>(2)?,row.get::<_,String>(3)?,row.get::<_,String>(4)?,row.get::<_,i64>(5)?,row.get::<_,i64>(6)?)))?.collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        let generations = GenerationRow::unpublished(db)?;
        let mut result = BTreeMap::new();
        for (id, project, state, json, sha, trashed, deadline) in stored {
            let descriptor: Descriptor = serde_json::from_str(&json)?;
            if descriptor.schema_version != 1
                || project != db.project_id()
                || descriptor.generation.id != id
                || descriptor.generation.state != "ABANDONED"
                || !matches!(
                    descriptor.original_state.as_str(),
                    "BUILDING" | "SEALED" | "ABANDONED"
                )
                || descriptor.trashed_at_ms != trashed
                || descriptor.purge_after_ms != deadline
                || deadline < trashed
                || digest(&json) != sha
                || digest(fs::read(witness_path(
                    db,
                    &descriptor.operation_id,
                    false,
                )?)?)
                    != sha
                || generations.get(&id) != Some(&descriptor.generation)
            {
                bail!("unpublished trash receipt differs from its catalog or witness: {id}");
            }
            descriptor.generation.verify_unpublished(db)?;
            let mut paths = BTreeSet::new();
            for entry in &descriptor.tree {
                if entry.path.as_os_str().is_empty()
                    || entry
                        .path
                        .components()
                        .any(|part| !matches!(part, std::path::Component::Normal(_)))
                    || !paths.insert(&entry.path)
                    || if entry.directory {
                        entry.bytes != 0 || entry.sha256.is_some()
                    } else {
                        entry
                            .sha256
                            .as_ref()
                            .is_none_or(|sha| !super::super::super::valid_sha256(sha))
                    }
                {
                    bail!("unpublished trash contains an invalid relative inventory");
                }
            }
            result.insert(id, Self { descriptor, state });
        }
        Ok(result)
    }
    pub fn transition(&self, db: &Database, next: &str) -> Result<()> {
        db.atomic(|db| {
            db.write(|conn| {
                if conn.execute("UPDATE unpublished_generation_trash SET state=?1 WHERE generation_id=?2 AND state=?3",params![next,self.descriptor.generation.id,self.state])?!=1 {bail!("unpublished trash state changed");}
                Ok(())
            })?;
            db.append_event("cleanup",None,"unpublished generation cleanup transition",serde_json::json!({"generation_id":self.descriptor.generation.id,"from":self.state,"to":next}))
        })
    }
}
