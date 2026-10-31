use super::*;
use std::fs::File;

impl Database {
    pub(super) fn generation_container(
        &self,
        id: &str,
        trash: bool,
        create_parent: bool,
    ) -> Result<PathBuf> {
        if uuid::Uuid::parse_str(id).is_err() {
            bail!("invalid generation retirement identity");
        }
        let root = fs::canonicalize(self.path().parent().context("project has no root")?)?;
        let mut parent = root;
        for name in if trash {
            vec!["trash", "generations"]
        } else {
            vec![".rnaa_generations"]
        } {
            parent.push(name);
            if create_parent {
                fs::create_dir_all(&parent)?;
            }
            if parent.symlink_metadata()?.file_type().is_symlink() {
                bail!("generation cleanup must not traverse symlinks");
            }
        }
        let path = parent.join(id);
        if path
            .symlink_metadata()
            .is_ok_and(|metadata| !metadata.is_dir() || metadata.file_type().is_symlink())
        {
            bail!("generation container must be an independent directory");
        }
        Ok(path)
    }
    fn validate_container(&self, path: &Path, partial: bool) -> Result<()> {
        for entry in fs::read_dir(path)? {
            let entry = entry?;
            if entry.file_type()?.is_symlink()
                || !matches!(entry.file_name().to_str(), Some("payload" | "seal.json"))
            {
                bail!("generation container has unexpected entries");
            }
        }
        if !partial && (!path.join("payload").is_dir() || !path.join("seal.json").is_file()) {
            bail!("generation container is incomplete");
        }
        Ok(())
    }
    pub(super) fn observe_live_generation(&self, seal: &Seal) -> Result<Vec<ManifestArtifact>> {
        if !self.verify_seal(seal)? {
            bail!("generation payload or pruning proof failed verification");
        }
        let root = self.generation_container(&seal.head.generation_id, false, false)?;
        self.validate_container(&root, false)?;
        let witness: Seal = serde_json::from_slice(&fs::read(root.join("seal.json"))?)?;
        if witness != *seal {
            bail!("generation seal witness changed");
        }
        inventory(&seal.head.payload_path, false)
    }
    pub(super) fn verify_retired_generation(
        &self,
        row: &Retirement,
        graph: &Graph,
        trash: bool,
    ) -> Result<()> {
        self.verify_retired_files(row, graph, trash, false)
    }
    fn verify_retired_files(
        &self,
        row: &Retirement,
        graph: &Graph,
        trash: bool,
        partial: bool,
    ) -> Result<()> {
        row.verify_catalog(self, graph)?;
        let root = self.generation_container(&row.id, trash, false)?;
        self.validate_container(&root, partial)?;
        let seal_path = root.join("seal.json");
        if (!partial || seal_path.exists()) && compute_sha256(&seal_path)? != row.seal_sha256 {
            bail!("retired generation seal witness changed");
        }
        let expected: Vec<ManifestArtifact> = serde_json::from_str(&row.payload_inventory_json)?;
        let payload = root.join("payload");
        let mut actual = if payload.exists() {
            inventory(&payload, false)?
        } else {
            Vec::new()
        };
        let original = &graph
            .nodes
            .get(&row.id)
            .context("retired generation disappeared")?
            .head
            .payload_path;
        for file in &mut actual {
            let relative = Path::new(&file.path).strip_prefix(&payload)?;
            file.path = original.join(relative).display().to_string();
        }
        if partial {
            if actual.iter().any(|file| !expected.contains(file)) {
                bail!("purging generation has unexpected or changed files");
            }
        } else if actual != expected {
            bail!("retired generation inventory changed");
        }
        Ok(())
    }
    pub(super) fn mark_generation_trash(
        &self,
        graph: &Graph,
        ids: &[String],
        observed: &BTreeMap<String, Vec<ManifestArtifact>>,
        days: u64,
    ) -> Result<()> {
        let now = now_ms()?;
        let deadline = now
            .checked_add(i64::try_from(
                days.checked_mul(86_400_000)
                    .context("generation trash retention is too large")?,
            )?)
            .context("generation trash deadline overflow")?;
        let operation = uuid::Uuid::new_v4().to_string();
        let mut rows = Vec::new();
        for id in ids {
            let seal = &graph.nodes[id];
            let root = self.generation_container(id, false, false)?;
            rows.push(Retirement {
                id: id.clone(),
                state: "TRASHING".into(),
                inventory_sha256: digest(serde_json::to_vec(seal)?),
                artifacts_sha256: digest(&graph.records[id]),
                payload_inventory_json: serde_json::to_string(&observed[id])?,
                seal_sha256: compute_sha256(&root.join("seal.json"))?,
                operation_id: operation.clone(),
                trashed_at_ms: now,
                purge_after_ms: deadline,
            });
        }
        self.atomic(|db| {
            for row in rows {
                db.write(|conn| {
                    conn.execute("INSERT INTO generation_retirements(generation_id,project_id,state,inventory_sha256,artifacts_sha256,payload_inventory_json,seal_sha256,operation_id,trashed_at_ms,purge_after_ms) VALUES (?1,?2,'TRASHING',?3,?4,?5,?6,?7,?8,?9)",
                        params![row.id,db.project_id(),row.inventory_sha256,row.artifacts_sha256,row.payload_inventory_json,row.seal_sha256,row.operation_id,row.trashed_at_ms,row.purge_after_ms])?;
                    Ok(())
                })?;
                db.append_event("cleanup",None,"generation retirement marked",serde_json::json!({"generation_id":row.id,"operation_id":operation,"purge_after_ms":deadline}))?;
            }
            Ok(())
        })
    }
    pub(super) fn mark_generation_purge(&self, graph: &Graph, ids: &[String]) -> Result<()> {
        for id in ids {
            self.verify_retired_generation(&graph.retired[id], graph, true)?;
        }
        let operation = uuid::Uuid::new_v4().to_string();
        self.atomic(|db|db.write(|conn| {
            for id in ids {conn.execute("UPDATE generation_retirements SET state='PURGING',operation_id=?1 WHERE generation_id=?2 AND state='TRASHED'",params![operation,id])?;}
            Ok(())
        }))
    }
    pub(super) fn move_generation(
        &self,
        row: &Retirement,
        graph: &Graph,
        restore: bool,
    ) -> Result<()> {
        let source = self.generation_container(&row.id, restore, false)?;
        let target = self.generation_container(&row.id, !restore, true)?;
        match (source.exists(), target.exists()) {
            (true, false) => {
                self.verify_retired_generation(row, graph, restore)?;
                fs::rename(&source, &target)?;
                File::open(source.parent().unwrap())?.sync_all()?;
                File::open(target.parent().unwrap())?.sync_all()?;
            }
            (false, true) => self.verify_retired_generation(row, graph, !restore)?,
            _ => bail!("ambiguous generation retirement/restore paths"),
        }
        Ok(())
    }
    pub(super) fn recover_retirements(&self) -> Result<()> {
        let graph = Graph::capture(self)?;
        let mut restores: BTreeMap<String, Vec<&Retirement>> = BTreeMap::new();
        for row in graph.retired.values() {
            match row.state.as_str() {
                "TRASHING" => {
                    self.move_generation(row, &graph, false)?;
                    self.write(|conn|{conn.execute("UPDATE generation_retirements SET state='TRASHED' WHERE generation_id=?1 AND state='TRASHING'",[&row.id])?;Ok(())})?;
                }
                "RESTORING" => {
                    restores
                        .entry(row.operation_id.clone())
                        .or_default()
                        .push(row);
                }
                "PURGING" => {
                    let root = self.generation_container(&row.id, true, false)?;
                    if root.exists() {
                        self.verify_retired_files(row, &graph, true, true)?;
                        fs::remove_dir_all(&root)?;
                        File::open(root.parent().unwrap())?.sync_all()?;
                    }
                    self.write(|conn|{conn.execute("UPDATE generation_retirements SET state='PURGED' WHERE generation_id=?1 AND state='PURGING'",[&row.id])?;Ok(())})?;
                }
                "TRASHED" | "PURGED" => {}
                _ => bail!("unsupported generation retirement state"),
            }
        }
        for (operation, rows) in restores {
            for row in &rows {
                self.move_generation(row, &graph, true)?;
            }
            for row in &rows {
                if !self.verify_seal(&graph.nodes[&row.id])? {
                    bail!("restored generation or retention evidence failed verification");
                }
            }
            self.atomic(|db|db.write(|conn|{conn.execute("DELETE FROM generation_retirements WHERE state='RESTORING' AND operation_id=?1",[operation])?;Ok(())}))?;
        }
        Ok(())
    }
}
