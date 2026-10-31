//! Protection is independent of trash eligibility and rebuilt for each maintenance pass.
use super::*;

pub(super) struct Protection {
    paths: BTreeSet<PathBuf>,
    legacy: bool,
}
impl Protection {
    pub fn capture(db: &Database) -> Result<Self> {
        let graph = Graph::capture(db)?;
        let mut value = Self {
            paths: BTreeSet::new(),
            legacy: graph.legacy,
        };
        let project = fs::canonicalize(db.path().parent().context("project root missing")?)?;
        let owned = graph
            .nodes
            .values()
            .flat_map(|seal| seal.files.iter())
            .map(|file| Path::new(&file.path))
            .collect::<BTreeSet<_>>();
        for record in db.list_artifacts()? {
            value.add(&project, &record.path)?;
            if let Some(path) = record.shared_path {
                value.add(&project, &path)?;
            }
            // An unsealed legacy result has no complete input contract. Until its
            // own adoption proves ancestry, it cannot authorize unpublished GC.
            let path = super::super::super::canonical_missing_path(Path::new(&record.path))?;
            if !owned.contains(path.as_path()) {
                value.legacy = true;
            }
        }
        for seal in graph.nodes.values() {
            for input in db.ancestry_inputs(seal)? {
                value.add(&project, &input.path)?;
            }
        }
        let bindings = db.read(|conn| {
            let mut query = conn
                .prepare("SELECT binding_id FROM generation_input_bindings WHERE project_id=?1")?;
            Ok(query
                .query_map([db.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        for id in bindings {
            let binding = db
                .read_input_binding(&id)?
                .context("unpublished cleanup binding disappeared")?;
            for input in binding.inputs {
                value.add(&project, &input.path)?;
            }
        }
        let retained = db.read(|conn| {
            let mut query = conn.prepare("SELECT path FROM retained_inputs WHERE project_id=?1")?;
            Ok(query
                .query_map([db.project_id()], |row| row.get::<_, String>(0))?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        for path in retained {
            value.add(&project, &path)?;
        }
        // Check both persisted settings and the current on-disk configuration.
        // Generic string traversal also covers newly added path settings.
        value.json(&project, &serde_json::from_str(&db.project()?.config_json)?)?;
        let config = project.join("rnaa.toml");
        if config.exists() {
            let parsed: toml::Value = toml::from_str(&fs::read_to_string(config)?)?;
            value.json(&project, &serde_json::to_value(parsed)?)?;
        }
        Ok(value)
    }
    fn json(&mut self, project: &Path, value: &serde_json::Value) -> Result<()> {
        match value {
            serde_json::Value::String(path) => self.add(project, path)?,
            serde_json::Value::Array(values) => {
                for value in values {
                    self.json(project, value)?;
                }
            }
            serde_json::Value::Object(values) => {
                for value in values.values() {
                    self.json(project, value)?;
                }
            }
            _ => {}
        }
        Ok(())
    }
    fn add(&mut self, project: &Path, raw: &str) -> Result<()> {
        if raw.is_empty() || (!Path::new(raw).is_absolute() && raw.contains("://")) {
            return Ok(());
        }
        let path = Path::new(raw);
        let path = if path.is_absolute() {
            path.to_path_buf()
        } else {
            project.join(path)
        };
        self.paths
            .insert(super::super::super::canonical_missing_path(&path)?);
        Ok(())
    }
    pub fn reasons(&self, db: &Database, id: &str) -> Result<Vec<String>> {
        let mut reasons = Vec::new();
        if self.legacy {
            reasons.push("unproven_legacy_ancestry".into());
        }
        for container in [
            tree::container(db, id, false, false)?,
            tree::container(db, id, true, false)?,
        ] {
            for path in &self.paths {
                if path.starts_with(&container) || container.starts_with(path) {
                    reasons.push(format!("referenced_path:{}", path.display()));
                }
            }
        }
        reasons.sort();
        reasons.dedup();
        Ok(reasons)
    }
}
