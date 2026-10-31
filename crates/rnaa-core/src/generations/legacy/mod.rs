//! Additive ancestry receipts; original generation inventories stay immutable.
mod catalog;
mod evidence;
#[cfg(test)]
mod tests;
use super::{GenerationDependency, OutputGroup, Seal};
use crate::{
    ArtifactKind, ArtifactRecord, Database, ManifestArtifact, StageManifest,
    util::{compute_sha256, now_rfc3339, write_json_pretty},
};
use anyhow::{Context, Result, bail};
use rusqlite::{OptionalExtension, params};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fs,
    path::{Path, PathBuf},
};
fn digest(bytes: impl AsRef<[u8]>) -> String {
    hex::encode(Sha256::digest(bytes.as_ref()))
}
fn identity(a: &ManifestArtifact, b: &ManifestArtifact) -> bool {
    a.path == b.path
        && a.checksum_type == b.checksum_type
        && a.checksum == b.checksum
        && a.bytes == b.bytes
}
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Evidence {
    file: ManifestArtifact,
    json: String,
}
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct ImportedLineage {
    schema_version: u32,
    import_id: String,
    project_id: String,
    generation_id: String,
    inventory_sha256: String,
    registry_sha256: String,
    seal_sha256: String,
    evidence: Vec<Evidence>,
    pub(super) inputs: Vec<ManifestArtifact>,
    pub(super) ancestry_inputs: Vec<ManifestArtifact>,
    pub(super) dependencies: Vec<GenerationDependency>,
    created_at: String,
}
#[derive(Serialize)]
pub struct LineageImportEntry {
    pub generation_id: String,
    pub group: String,
    pub subject: String,
    pub action: String,
    pub dependency_count: usize,
    pub reason: Option<String>,
}
#[derive(Serialize)]
pub struct LineageImportReport {
    pub applied: bool,
    pub busy: bool,
    pub complete: bool,
    pub entries: Vec<LineageImportEntry>,
    pub errors: Vec<String>,
}

impl Database {
    pub fn import_generation_lineage(&self, apply: bool) -> Result<LineageImportReport> {
        let mut report = LineageImportReport {
            applied: apply,
            busy: false,
            complete: false,
            entries: Vec::new(),
            errors: Vec::new(),
        };
        let Some(_lock) = super::pins::maintenance_lock(self)? else {
            report.busy = true;
            return Ok(report);
        };
        let mut db = self.clone();
        db.generation_maintenance = true;
        let ids=db.read(|conn| {let mut query=conn.prepare("SELECT generation_id FROM artifact_generations WHERE project_id=?1 AND state='COMMITTED' ORDER BY generation_id")?;Ok(query.query_map([db.project_id()],|row|row.get::<_,String>(0))?.collect::<std::result::Result<Vec<_>,_>>()?)})?;
        let mut plans = Vec::new();
        let mut graph = BTreeMap::<String, Vec<String>>::new();
        for id in ids {
            let seal = db
                .committed_seal(&id)?
                .context("legacy generation disappeared")?;
            let result =
                (|| -> Result<(String, Vec<GenerationDependency>, Option<ImportedLineage>)> {
                    if !db.lineage_matches_catalog(&seal)? {
                        bail!("generation catalog integrity failed");
                    }
                    if seal.lineage_version == 1 {
                        return Ok(("native_recorded".into(), seal.dependencies.clone(), None));
                    }
                    if let Some(imported) = db.imported_lineage(&seal)? {
                        return Ok(("already_imported".into(), imported.dependencies, None));
                    }
                    if db.generation_retired(&id)? {
                        bail!("restore this retired legacy generation before ancestry import");
                    }
                    if !db.verify_seal(&seal)? {
                        bail!(
                            "legacy generation payload or pruning proof failed verification; restore or reconstruct its inputs first"
                        );
                    }
                    let evidence = db.legacy_manifest_evidence(&seal)?;
                    let inputs = db.legacy_inputs(&seal, &evidence, None)?;
                    let managed = fs::canonicalize(
                        db.path().parent().context("legacy project root missing")?,
                    )?
                    .join(".rnaa_generations");
                    for input in &inputs.ancestry {
                        let canonical = super::canonical_missing_path(Path::new(&input.path))?;
                        if canonical.starts_with(&managed) && canonical != Path::new(&input.path) {
                            bail!(
                                "legacy managed source alias lacks an immutable canonical producer witness: {}",
                                input.path
                            );
                        }
                    }
                    let upstream = inputs
                        .ancestry
                        .iter()
                        .filter(|input| {
                            !Path::new(&input.path).starts_with(&seal.head.payload_path)
                        })
                        .cloned()
                        .collect::<Vec<_>>();
                    let dependencies =
                        db.resolve_generation_dependencies(&upstream, &seal.head.payload_path)?;
                    let (inventory, registry) = db.legacy_catalog_documents(&id)?;
                    let imported = ImportedLineage {
                        schema_version: 1,
                        import_id: uuid::Uuid::new_v4().to_string(),
                        project_id: db.project_id().into(),
                        generation_id: id.clone(),
                        inventory_sha256: digest(inventory),
                        registry_sha256: digest(registry),
                        seal_sha256: digest(serde_json::to_vec(&seal)?),
                        evidence,
                        inputs: inputs.inputs,
                        ancestry_inputs: inputs.ancestry,
                        dependencies: dependencies.clone(),
                        created_at: now_rfc3339(),
                    };
                    db.verify_lineage_import(&seal, &imported, false)?;
                    Ok(("import".into(), dependencies, Some(imported)))
                })();
            match result {
                Ok((action, dependencies, imported)) => {
                    report.entries.push(LineageImportEntry {
                        generation_id: id.clone(),
                        group: seal.head.output_group.as_str().into(),
                        subject: seal.head.subject.clone(),
                        action,
                        dependency_count: dependencies.len(),
                        reason: None,
                    });
                    graph.insert(
                        id,
                        dependencies
                            .into_iter()
                            .map(|edge| edge.producer_generation_id)
                            .collect(),
                    );
                    if let Some(imported) = imported {
                        plans.push(imported);
                    }
                }
                Err(error) => {
                    let reason = format!("{error:#}");
                    report.errors.push(format!("{id}: {reason}"));
                    report.entries.push(LineageImportEntry {
                        generation_id: id.clone(),
                        group: seal.head.output_group.as_str().into(),
                        subject: seal.head.subject,
                        action: "held".into(),
                        dependency_count: 0,
                        reason: Some(reason),
                    });
                    graph.insert(id, Vec::new());
                }
            }
        }
        verify_acyclic(&graph)?;
        if apply && !plans.is_empty() {
            for plan in &plans {
                let path = db.lineage_import_path(&plan.import_id, true)?;
                if path.exists() {
                    bail!("lineage import witness already exists");
                }
                write_json_pretty(&path, plan)?;
            }
            db.atomic(|db| {
                for plan in &plans {
                    let seal=db.committed_seal(&plan.generation_id)?.context("legacy generation disappeared before import")?;
                    db.verify_lineage_import(&seal,plan,false)?;
                    let json=serde_json::to_string(plan)?;
                    db.write(|conn| {
                        conn.execute("INSERT INTO generation_lineage_imports(generation_id,project_id,import_id,descriptor_json,descriptor_sha256,created_at) VALUES (?1,?2,?3,?4,?5,?6)",params![plan.generation_id,db.project_id(),plan.import_id,json,digest(json.as_bytes()),plan.created_at])?;
                        for edge in &plan.dependencies {conn.execute("INSERT INTO generation_import_dependencies(consumer_generation_id,producer_generation_id,producer_seal_sha256,role,path,artifact_json) VALUES (?1,?2,?3,?4,?5,?6)",params![plan.generation_id,edge.producer_generation_id,edge.producer_seal_sha256,edge.input.kind,edge.input.path,serde_json::to_string(&edge.input)?])?;}
                        Ok(())
                    })?;
                    db.append_event("lineage_import",None,"legacy ancestry imported without changing its seal",serde_json::json!({"generation_id":plan.generation_id,"import_id":plan.import_id,"dependencies":plan.dependencies.len()}))?;
                }
                Ok(())
            })?;
            for entry in &mut report.entries {
                if entry.action == "import" {
                    entry.action = "imported".into();
                }
            }
        }
        report.complete = report.errors.is_empty();
        Ok(report)
    }
}
fn verify_acyclic(graph: &BTreeMap<String, Vec<String>>) -> Result<()> {
    let mut incoming = graph
        .keys()
        .map(|id| (id.clone(), 0usize))
        .collect::<BTreeMap<_, _>>();
    for targets in graph.values() {
        for target in targets {
            *incoming
                .get_mut(target)
                .context("legacy ancestry producer is not committed")? += 1;
        }
    }
    let mut queue = incoming
        .iter()
        .filter(|(_, n)| **n == 0)
        .map(|(id, _)| id.clone())
        .collect::<VecDeque<_>>();
    let mut count = 0;
    while let Some(id) = queue.pop_front() {
        count += 1;
        for target in &graph[&id] {
            let n = incoming.get_mut(target).unwrap();
            *n -= 1;
            if *n == 0 {
                queue.push_back(target.clone());
            }
        }
    }
    if count != graph.len() {
        bail!("legacy ancestry import would create a dependency cycle");
    }
    Ok(())
}
