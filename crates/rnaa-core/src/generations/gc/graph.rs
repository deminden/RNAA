use super::*;
use crate::eligibility::Eligibility;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

pub(super) struct Graph {
    pub nodes: BTreeMap<String, Seal>,
    pub records: BTreeMap<String, String>,
    pub edges: BTreeMap<String, BTreeSet<String>>,
    pub roots: BTreeMap<String, BTreeSet<String>>,
    pub retired: BTreeMap<String, Retirement>,
    pub legacy: bool,
}
impl Graph {
    pub fn capture(db: &Database) -> Result<Self> {
        let rows=db.read(|conn| {
            let mut query=conn.prepare("SELECT generation_id,artifacts_json FROM artifact_generations WHERE project_id=?1 AND state='COMMITTED'")?;
            Ok(query.query_map([db.project_id()],|row|Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?)))?.collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        let mut graph = Self {
            nodes: BTreeMap::new(),
            records: BTreeMap::new(),
            edges: BTreeMap::new(),
            roots: BTreeMap::new(),
            retired: Retirement::all(db)?,
            legacy: false,
        };
        for (id, records) in rows {
            let seal = db
                .committed_seal(&id)?
                .context("generation disappeared from cleanup snapshot")?;
            if !db.lineage_matches_catalog(&seal)? {
                bail!("generation lineage catalog is invalid: {id}");
            }
            let artifacts: Vec<ArtifactRecord> = serde_json::from_str(&records)?;
            if artifacts.is_empty()
                || artifacts.iter().any(|record| {
                    record.project_id != db.project_id()
                        || !seal.files.iter().any(|file| {
                            file.path == record.path
                                && file.checksum == record.checksum
                                && file.bytes == record.bytes
                        })
                })
            {
                bail!("generation registry is invalid: {id}");
            }
            let (recorded, dependencies) = db.effective_dependencies(&seal)?;
            graph.legacy |= !recorded;
            graph.edges.insert(
                id.clone(),
                dependencies
                    .iter()
                    .map(|edge| edge.producer_generation_id.clone())
                    .collect(),
            );
            graph.records.insert(id.clone(), records);
            graph.nodes.insert(id, seal);
        }
        for retired in graph.retired.values() {
            retired.verify_catalog(db, &graph)?;
        }
        let heads = db.generation_heads()?;
        for head in heads {
            graph.root(&head.generation_id, "current_head")?;
        }
        let selection: Option<String> = db.read(|conn| {
            Ok(conn
                .query_row(
                    "SELECT generation_id FROM reference_selection WHERE project_id=?1",
                    [db.project_id()],
                    |row| row.get(0),
                )
                .optional()?)
        })?;
        if let Some(id) = selection {
            db.selected_reference()?
                .context("selected reference disappeared")?;
            graph.root(&id, "selected_reference")?;
        }
        let kept = db.read(|conn| {
            let mut query = conn.prepare("SELECT generation_id,label FROM generation_keeps")?;
            Ok(query
                .query_map([], |row| {
                    Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
                })?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        for (id, label) in kept {
            graph.root(&id, &format!("keep:{label}"))?;
        }
        let active = db.read(|conn| {
            let mut query = conn.prepare("SELECT path FROM artifacts WHERE project_id=?1")?;
            Ok(query
                .query_map([db.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        let owners = graph
            .nodes
            .iter()
            .flat_map(|(id, seal)| {
                seal.files
                    .iter()
                    .map(move |file| (file.path.clone(), id.clone()))
            })
            .collect::<BTreeMap<_, _>>();
        let managed = fs::canonicalize(db.path().parent().context("project has no root")?)?
            .join(".rnaa_generations");
        // A recorded external source can be retargeted to a managed payload after
        // publication. Protect its present alias target without inventing a new
        // historical producer edge or changing the consumer's immutable lineage.
        let mut aliases = Vec::new();
        for (id, seal) in &graph.nodes {
            for input in db.ancestry_inputs(seal)? {
                let canonical = super::super::canonical_missing_path(Path::new(&input.path))?;
                if canonical != Path::new(&input.path) && canonical.starts_with(&managed) {
                    let owner = owners.get(canonical.to_string_lossy().as_ref()).context("recorded source alias points to a managed path without committed ownership")?;
                    aliases.push((id.clone(), owner.clone()));
                }
            }
        }
        for (consumer, producer) in aliases {
            graph.edge(&consumer, &producer)?;
        }
        for path in active {
            let canonical = super::super::canonical_missing_path(Path::new(&path))?;
            if let Some(id) = owners.get(canonical.to_string_lossy().as_ref()) {
                graph.root(id, "active_artifact")?;
            } else if canonical.starts_with(&managed) {
                bail!("active managed artifact has no generation owner: {path}");
            }
        }
        let decisions=db.read(|conn| {
            let mut query=conn.prepare("SELECT d.descriptor_json,d.fingerprint FROM eligibility_heads h JOIN eligibility_decisions d ON d.decision_id=h.decision_id WHERE h.project_id=?1")?;
            Ok(query.query_map([db.project_id()],|row|Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?)))?.collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        for (json, fingerprint) in decisions {
            let eligibility: Eligibility = serde_json::from_str(&json)?;
            if eligibility.fingerprint()? != fingerprint {
                bail!("QC decision differs from its recorded fingerprint");
            }
            if let Some(id) = eligibility
                .evidence
                .and_then(|evidence| evidence.producer_generation)
            {
                graph.root(&id, "qc_evidence")?;
            }
        }
        let binding_ids = db.read(|conn| {
            let mut query = conn
                .prepare("SELECT binding_id FROM generation_input_bindings WHERE project_id=?1")?;
            Ok(query
                .query_map([db.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        let mut bindings = BTreeMap::new();
        for id in binding_ids {
            let binding = db
                .read_input_binding(&id)?
                .context("cleanup input binding disappeared")?;
            db.verify_input_binding_catalog(&binding)?;
            for edge in &binding.dependencies {
                graph.edge(&binding.generation_id, &edge.producer_generation_id)?;
            }
            bindings.insert(id, binding.generation_id);
        }
        let pruning=db.read(|conn| {
            let mut query=conn.prepare("SELECT producer_generation_id,consumer_generation_id FROM generation_pruned_files WHERE project_id=?1")?;
            Ok(query.query_map([db.project_id()],|row|Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?)))?.collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        for (producer, consumer) in pruning {
            graph.edge(&producer, &consumer)?;
        }
        let binding_pruning = db.read(|conn| {
            let mut query = conn.prepare(
                "SELECT producer_generation_id,binding_id FROM generation_binding_pruning",
            )?;
            Ok(query
                .query_map([], |row| {
                    Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
                })?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        for (producer, binding) in binding_pruning {
            graph.edge(
                &producer,
                bindings
                    .get(&binding)
                    .context("pruning binding is missing")?,
            )?;
        }
        for link in db.file_trash_generation_links()? {
            if link.restorable {
                graph.root(&link.consumer, "restorable_file_trash")?;
                if let Some(producer) = &link.producer {
                    graph.root(producer, "restorable_file_trash")?;
                }
            }
            if link.absence {
                graph.root(&link.consumer, "file_trash_input_proof")?;
                if let Some(producer) = &link.producer {
                    graph.edge(producer, &link.consumer)?;
                }
            }
        }
        let retained = db.read(|conn| {
            let mut query = conn.prepare(
                "SELECT DISTINCT generation_sha256 FROM retained_inputs WHERE project_id=?1",
            )?;
            Ok(query
                .query_map([db.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        for digest in retained {
            let witnesses = graph
                .nodes
                .iter()
                .filter(|(_, seal)| {
                    super::super::quant_manifest_digest(seal) == Some(digest.as_str())
                })
                .map(|(id, _)| id.clone())
                .collect::<Vec<_>>();
            if witnesses.is_empty() {
                graph.legacy = true;
            }
            for id in witnesses {
                graph.root(&id, "retained_input_proof")?;
            }
        }
        if graph.legacy {
            for id in graph.nodes.keys().cloned().collect::<Vec<_>>() {
                if graph
                    .retired
                    .get(&id)
                    .is_none_or(|row| row.state != "PURGED")
                {
                    graph.root(&id, "unproven_legacy_ancestry")?;
                }
            }
        }
        Ok(graph)
    }
    fn root(&mut self, id: &str, reason: &str) -> Result<()> {
        if !self.nodes.contains_key(id) {
            bail!("cleanup root has no committed generation: {id}");
        }
        self.roots
            .entry(id.into())
            .or_default()
            .insert(reason.into());
        Ok(())
    }
    fn edge(&mut self, from: &str, to: &str) -> Result<()> {
        if !self.nodes.contains_key(from) || !self.nodes.contains_key(to) {
            bail!("cleanup dependency has no committed generation");
        }
        self.edges.entry(from.into()).or_default().insert(to.into());
        Ok(())
    }
    pub fn closure(&self, roots: impl IntoIterator<Item = String>) -> BTreeSet<String> {
        let mut visited = BTreeSet::new();
        let mut pending = VecDeque::from_iter(roots);
        while let Some(id) = pending.pop_front() {
            if visited.insert(id.clone()) {
                pending.extend(self.edges.get(&id).into_iter().flatten().cloned());
            }
        }
        visited
    }
}
