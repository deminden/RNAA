//! Reversible generation retirement. A project-wide exclusive lock excludes live
//! workers, path snapshots and database mutations throughout planning and moves.
mod files;
mod graph;
mod unpublished;
pub use unpublished::{UnpublishedCleanupEntry, UnpublishedCleanupReport};
#[cfg(test)]
mod tests;
use super::{Seal, inventory};
use crate::{
    ArtifactRecord, Database, ManifestArtifact,
    util::{compute_sha256, now_rfc3339},
};
use anyhow::{Context, Result, bail};
use graph::Graph;
use rusqlite::{OptionalExtension, params};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

fn digest(bytes: impl AsRef<[u8]>) -> String {
    hex::encode(Sha256::digest(bytes.as_ref()))
}
fn now_ms() -> Result<i64> {
    Ok(i64::try_from(
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_millis(),
    )?)
}
#[derive(Clone, Copy)]
pub struct GenerationGcOptions {
    pub apply: bool,
    pub purge_expired: bool,
    pub trash_days: u64,
}
#[derive(Serialize)]
pub struct GenerationGcEntry {
    pub generation_id: String,
    pub group: String,
    pub subject: String,
    pub state: String,
    pub action: String,
    pub reasons: Vec<String>,
    pub payload_bytes: Option<u64>,
}
#[derive(Serialize)]
pub struct GenerationGcReport {
    pub applied: bool,
    pub busy: bool,
    pub entries: Vec<GenerationGcEntry>,
    pub errors: Vec<String>,
}

#[derive(Clone, Serialize, Deserialize)]
pub(super) struct Retirement {
    id: String,
    state: String,
    inventory_sha256: String,
    artifacts_sha256: String,
    payload_inventory_json: String,
    seal_sha256: String,
    operation_id: String,
    trashed_at_ms: i64,
    purge_after_ms: i64,
}
impl Retirement {
    fn all(db: &Database) -> Result<BTreeMap<String, Self>> {
        db.read(|conn| {
            let mut query=conn.prepare("SELECT generation_id,state,inventory_sha256,artifacts_sha256,payload_inventory_json,seal_sha256,operation_id,trashed_at_ms,purge_after_ms FROM generation_retirements WHERE project_id=?1")?;
            Ok(query.query_map([db.project_id()],|row|Ok(Self {id:row.get(0)?,state:row.get(1)?,inventory_sha256:row.get(2)?,artifacts_sha256:row.get(3)?,payload_inventory_json:row.get(4)?,seal_sha256:row.get(5)?,operation_id:row.get(6)?,trashed_at_ms:row.get(7)?,purge_after_ms:row.get(8)?}))?
                .map(|row|row.map(|retired|(retired.id.clone(),retired))).collect::<std::result::Result<BTreeMap<_,_>,_>>()?)
        })
    }
    fn verify_catalog(&self, db: &Database, graph: &Graph) -> Result<()> {
        let seal = graph
            .nodes
            .get(&self.id)
            .context("retirement has no committed generation")?;
        if self.inventory_sha256 != digest(serde_json::to_vec(seal)?)
            || self.artifacts_sha256
                != digest(
                    graph
                        .records
                        .get(&self.id)
                        .context("retirement registry is missing")?,
                )
            || self.purge_after_ms < self.trashed_at_ms
            || uuid::Uuid::parse_str(&self.operation_id).is_err()
        {
            bail!(
                "generation retirement receipt differs from its catalog: {}",
                self.id
            );
        }
        let observed: Vec<ManifestArtifact> = serde_json::from_str(&self.payload_inventory_json)?;
        if observed.iter().any(|file| !seal.files.contains(file)) {
            bail!("retirement contains unsealed payloads");
        }
        let _ = db;
        Ok(())
    }
}

impl Database {
    pub fn generation_storage_status(&self) -> Result<serde_json::Value> {
        let (pins,keeps)=self.read(|conn| {
            let mut query=conn.prepare("SELECT pin_id,owner,created_at FROM generation_store_pins ORDER BY created_at,pin_id")?;
            let pins=query.query_map([],|row|Ok(serde_json::json!({"pin_id":row.get::<_,String>(0)?,"owner":row.get::<_,String>(1)?,"created_at":row.get::<_,String>(2)?})))?.collect::<std::result::Result<Vec<_>,_>>()?;
            let mut query=conn.prepare("SELECT generation_id,label,created_at FROM generation_keeps ORDER BY generation_id")?;
            let keeps=query.query_map([],|row|Ok(serde_json::json!({"generation_id":row.get::<_,String>(0)?,"label":row.get::<_,String>(1)?,"created_at":row.get::<_,String>(2)?})))?.collect::<std::result::Result<Vec<_>,_>>()?;
            Ok((pins,keeps))
        })?;
        let retired=Retirement::all(self)?.values().map(|row|serde_json::json!({"generation_id":row.id,"state":row.state,"trashed_at_ms":row.trashed_at_ms,"purge_after_ms":row.purge_after_ms,"operation_id":row.operation_id})).collect::<Vec<_>>();
        Ok(
            serde_json::json!({"schema_version":1,"pin_observations":pins,"keeps":keeps,"retirements":retired}),
        )
    }
    pub fn generation_keep(&self, id: &str, label: Option<&str>) -> Result<()> {
        self.atomic(|db| {
            if db.committed_seal(id)?.is_none() {bail!("generation is not committed");}
            if let Some(label)=label {
                if db.generation_retired(id)? {bail!("restore the generation before keeping it");}
                db.write(|conn| {conn.execute("INSERT INTO generation_keeps(generation_id,label,created_at) VALUES (?1,?2,?3) ON CONFLICT(generation_id) DO UPDATE SET label=excluded.label",params![id,label,now_rfc3339()])?;Ok(())})
            } else {db.write(|conn| {conn.execute("DELETE FROM generation_keeps WHERE generation_id=?1",[id])?;Ok(())})}
        })
    }
    pub(crate) fn generation_retired(&self, id: &str) -> Result<bool> {
        self.read(|conn| {
            Ok(conn.query_row(
                "SELECT EXISTS(SELECT 1 FROM generation_retirements WHERE generation_id=?1)",
                [id],
                |row| row.get(0),
            )?)
        })
    }
    pub fn generation_gc(&self, options: GenerationGcOptions) -> Result<GenerationGcReport> {
        let mut report = GenerationGcReport {
            applied: options.apply,
            busy: false,
            entries: Vec::new(),
            errors: Vec::new(),
        };
        let Some(_lock) = super::pins::maintenance_lock(self)? else {
            report.busy = true;
            return Ok(report);
        };
        let mut db = self.clone();
        db.generation_maintenance = true;
        if options.apply {
            db.recover_retirements()?;
            db.write(|conn| {
                conn.execute("DELETE FROM generation_store_pins", [])?;
                Ok(())
            })?;
        }
        let graph = Graph::capture(&db)?;
        let mut reachable = graph.closure(graph.roots.keys().cloned());
        let mut observed = BTreeMap::new();
        let mut invalid = BTreeSet::new();
        // A damaged result stays in place and also keeps its proof dependencies.
        for (id, seal) in &graph.nodes {
            if !reachable.contains(id) && !graph.retired.contains_key(id) {
                match db.observe_live_generation(seal) {
                    Ok(files) => {
                        observed.insert(id.clone(), files);
                    }
                    Err(error) => {
                        invalid.insert(id.clone());
                        report.errors.push(format!("{id}: {error:#}"));
                        reachable.extend(graph.closure([id.clone()]));
                    }
                }
            }
        }
        let now = now_ms()?;
        let purge_candidates = graph
            .retired
            .iter()
            .filter(|(_, row)| row.state == "TRASHED" && row.purge_after_ms <= now)
            .map(|(id, _)| id.clone())
            .collect::<BTreeSet<_>>();
        let purge_roots = graph
            .nodes
            .keys()
            .filter(|id| {
                !purge_candidates.contains(*id)
                    && graph
                        .retired
                        .get(*id)
                        .is_none_or(|row| row.state != "PURGED")
            })
            .cloned()
            .chain(graph.roots.keys().cloned());
        let purge_protected = graph.closure(purge_roots);
        let mut trash = Vec::new();
        let mut purge = Vec::new();
        for (id, seal) in &graph.nodes {
            let retired = graph.retired.get(id);
            let mut entry = GenerationGcEntry {
                generation_id: id.clone(),
                group: seal.head.output_group.as_str().into(),
                subject: seal.head.subject.clone(),
                state: retired.map_or("LIVE", |row| row.state.as_str()).into(),
                action: String::new(),
                reasons: graph.roots.get(id).into_iter().flatten().cloned().collect(),
                payload_bytes: observed
                    .get(id)
                    .map(|files| files.iter().map(|file| file.bytes).sum()),
            };
            if let Some(row) = retired {
                entry.action = if row.state == "PURGED" {
                    "already_purged"
                } else if row.state != "TRASHED" {
                    "resume_pending"
                } else if options.purge_expired
                    && purge_candidates.contains(id)
                    && !purge_protected.contains(id)
                {
                    purge.push(id.clone());
                    "purge"
                } else {
                    "held_trash_window_or_dependency"
                }
                .into();
            } else if invalid.contains(id) {
                entry.action = "held_invalid".into();
            } else if reachable.contains(id) {
                entry.action = "held_reachable".into();
                if entry.reasons.is_empty() {
                    entry.reasons.push("dependency_or_retention_proof".into());
                }
            } else if observed.contains_key(id) {
                entry.action = "trash".into();
                trash.push(id.clone());
            } else {
                entry.action = "held_invalid".into();
            }
            report.entries.push(entry);
        }
        if options.apply {
            db.mark_generation_trash(&graph, &trash, &observed, options.trash_days)?;
            db.recover_retirements()?;
            // Only trash present at the beginning of this operation can be purged.
            db.mark_generation_purge(&graph, &purge)?;
            db.recover_retirements()?;
            let states = Retirement::all(&db)?;
            for entry in &mut report.entries {
                entry.state = states
                    .get(&entry.generation_id)
                    .map_or("LIVE", |row| row.state.as_str())
                    .into();
            }
        }
        Ok(report)
    }
    pub fn generation_restore(&self, id: &str) -> Result<()> {
        let _lock = super::pins::maintenance_lock(self)?
            .context("generation store is busy; retry restore")?;
        let mut db = self.clone();
        db.generation_maintenance = true;
        db.recover_retirements()?;
        let graph = Graph::capture(&db)?;
        if !graph.nodes.contains_key(id) {
            bail!("unknown committed generation");
        }
        let needed = graph.closure([id.into()]);
        let mut restoring = Vec::new();
        for id in needed {
            if let Some(row) = graph.retired.get(&id) {
                if row.state != "TRASHED" {
                    bail!("generation or required dependency has already been purged: {id}");
                }
                db.verify_retired_generation(row, &graph, true)?;
                restoring.push(id);
            }
        }
        let operation = uuid::Uuid::new_v4().to_string();
        db.atomic(|db|db.write(|conn| {for id in restoring {conn.execute("UPDATE generation_retirements SET state='RESTORING',operation_id=?1 WHERE generation_id=?2 AND state='TRASHED'",params![operation,id])?;}Ok(())}))?;
        db.recover_retirements()
    }
}
