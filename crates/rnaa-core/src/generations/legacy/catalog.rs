use super::*;
impl Database {
    pub(super) fn legacy_catalog_documents(&self, id: &str) -> Result<(String, String)> {
        self.read(|conn|Ok(conn.query_row("SELECT inventory_json,artifacts_json FROM artifact_generations WHERE project_id=?1 AND generation_id=?2 AND state='COMMITTED'",params![self.project_id(),id],|row|Ok((row.get(0)?,row.get(1)?)))?))
    }
    pub(super) fn lineage_import_path(&self, id: &str, create: bool) -> Result<PathBuf> {
        uuid::Uuid::parse_str(id)?;
        let mut path = fs::canonicalize(
            self.path()
                .parent()
                .context("lineage import project root missing")?,
        )?;
        for name in [".rnaa_support", "lineage_imports"] {
            path.push(name);
            match path.symlink_metadata() {
                Ok(meta) if !meta.is_dir() || meta.file_type().is_symlink() => {
                    bail!("lineage import witness paths must use real directories")
                }
                Ok(_) => {}
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                    if create {
                        fs::create_dir(&path)?;
                    }
                }
                Err(error) => return Err(error.into()),
            }
        }
        path.push(format!("{id}.json"));
        if path.symlink_metadata().is_ok_and(|meta| !meta.is_file()) {
            bail!("lineage import witness must be a regular file");
        }
        Ok(path)
    }
    pub(super) fn verify_lineage_import(
        &self,
        seal: &Seal,
        imported: &ImportedLineage,
        stored: bool,
    ) -> Result<()> {
        let (inventory, registry) = self.legacy_catalog_documents(&seal.head.generation_id)?;
        let derived =
            self.legacy_inputs(seal, &imported.evidence, Some(&imported.ancestry_inputs))?;
        if imported.schema_version != 1
            || uuid::Uuid::parse_str(&imported.import_id).is_err()
            || seal.lineage_version != 0
            || !seal.dependencies.is_empty()
            || imported.project_id != self.project_id()
            || imported.generation_id != seal.head.generation_id
            || imported.inventory_sha256 != digest(&inventory)
            || imported.registry_sha256 != digest(registry)
            || imported.seal_sha256 != digest(serde_json::to_vec(seal)?)
            || imported.inputs != derived.inputs
            || imported.ancestry_inputs != derived.ancestry
        {
            bail!("lineage import differs from its original sealed evidence");
        }
        let registry_records: Vec<ArtifactRecord> =
            serde_json::from_str(&self.legacy_catalog_documents(&seal.head.generation_id)?.1)?;
        let expected = registry_records
            .iter()
            .filter(|record| record.kind == ArtifactKind::Manifest)
            .map(|record| record.path.clone())
            .collect::<BTreeSet<_>>();
        let observed = imported
            .evidence
            .iter()
            .map(|item| item.file.path.clone())
            .collect::<BTreeSet<_>>();
        if expected != observed || observed.len() != imported.evidence.len() {
            bail!("lineage import omits or duplicates sealed manifest evidence");
        }
        let native: i64 = self.read(|conn| {
            Ok(conn.query_row(
                "SELECT count(*) FROM generation_dependencies WHERE consumer_generation_id=?1",
                [&seal.head.generation_id],
                |row| row.get(0),
            )?)
        })?;
        if native != 0 {
            bail!("legacy generation has unexpected native dependency rows");
        }
        let mut unique = BTreeSet::new();
        for edge in &imported.dependencies {
            if edge.producer_generation_id == seal.head.generation_id
                || !imported.ancestry_inputs.contains(&edge.input)
                || !unique.insert((&edge.input.kind, &edge.input.path))
            {
                bail!("invalid imported producer edge");
            }
        }
        self.verify_dependency_catalog(&imported.dependencies)?;
        // Resolve exact indexed ownership without touching present-day input bytes.
        // A missing managed producer can never become an external source silently.
        let mut expected = Vec::new();
        let root = fs::canonicalize(self.path().parent().context("lineage root missing")?)?
            .join(".rnaa_generations");
        for input in &imported.ancestry_inputs {
            if Path::new(&input.path).starts_with(&seal.head.payload_path) {
                if !seal.files.iter().any(|file| identity(file, input)) {
                    bail!("consumer-owned legacy input is not sealed");
                }
                continue;
            }
            let owners=self.read(|conn| {let mut query=conn.prepare("SELECT g.inventory_json FROM generation_files f JOIN artifact_generations g ON g.generation_id=f.generation_id WHERE g.project_id=?1 AND g.state='COMMITTED' AND f.path=?2 AND f.sha256=?3 AND f.bytes=?4")?;Ok(query.query_map(params![self.project_id(),input.path,input.checksum,i64::try_from(input.bytes)?],|row|row.get::<_,String>(0))?.collect::<std::result::Result<Vec<_>,_>>()?)})?;
            if owners.len() > 1 {
                bail!("legacy input has ambiguous producer ownership");
            }
            if let Some(owner) = owners.first() {
                let producer: Seal = serde_json::from_str(owner)?;
                if !producer.files.iter().any(|file| identity(file, input)) {
                    bail!("legacy producer differs from its ownership index");
                }
                expected.push(GenerationDependency {
                    producer_generation_id: producer.head.generation_id.clone(),
                    producer_seal_sha256: digest(serde_json::to_vec(&producer)?),
                    input: input.clone(),
                });
            } else if Path::new(&input.path).starts_with(&root) {
                bail!("legacy managed input has no exact indexed producer");
            }
        }
        expected.sort_by(|a, b| {
            (&a.producer_generation_id, &a.input.kind, &a.input.path).cmp(&(
                &b.producer_generation_id,
                &b.input.kind,
                &b.input.path,
            ))
        });
        if expected != imported.dependencies {
            bail!("imported ancestry omits or changes a producer");
        }
        if stored {
            let rows=self.read(|conn| {let mut query=conn.prepare("SELECT producer_generation_id,producer_seal_sha256,role,path,artifact_json FROM generation_import_dependencies WHERE consumer_generation_id=?1 ORDER BY producer_generation_id,role,path")?;Ok(query.query_map([&seal.head.generation_id],|row|Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,String>(2)?,row.get::<_,String>(3)?,row.get::<_,String>(4)?)))?.collect::<std::result::Result<Vec<_>,_>>()?)})?;
            let recorded = rows
                .into_iter()
                .map(
                    |(id, sha, role, path, json)| -> Result<GenerationDependency> {
                        let input: ManifestArtifact = serde_json::from_str(&json)?;
                        if input.kind != role || input.path != path {
                            bail!("imported edge index differs from its input");
                        }
                        Ok(GenerationDependency {
                            producer_generation_id: id,
                            producer_seal_sha256: sha,
                            input,
                        })
                    },
                )
                .collect::<Result<Vec<_>>>()?;
            if recorded != expected {
                bail!("imported lineage catalog differs from its receipt");
            }
        }
        Ok(())
    }
    pub(in crate::generations) fn imported_lineage(
        &self,
        seal: &Seal,
    ) -> Result<Option<ImportedLineage>> {
        let row:Option<(String,String,String,String)>=self.read(|conn|Ok(conn.query_row("SELECT import_id,descriptor_json,descriptor_sha256,created_at FROM generation_lineage_imports WHERE generation_id=?1 AND project_id=?2",params![seal.head.generation_id,self.project_id()],|row|Ok((row.get(0)?,row.get(1)?,row.get(2)?,row.get(3)?))).optional()?))?;
        row.map(|(id, json, sha, created)| {
            let imported: ImportedLineage = serde_json::from_str(&json)?;
            if serde_json::to_string(&imported)? != json
                || digest(json.as_bytes()) != sha
                || imported.import_id != id
                || imported.created_at != created
            {
                bail!("lineage import receipt failed integrity verification");
            }
            self.verify_lineage_import(seal, &imported, true)?;
            let witness: serde_json::Value =
                serde_json::from_slice(&fs::read(self.lineage_import_path(&id, false)?)?)?;
            if witness != serde_json::to_value(&imported)? {
                bail!("lineage import physical witness changed");
            }
            Ok(imported)
        })
        .transpose()
    }
    pub(in crate::generations) fn effective_inputs(
        &self,
        seal: &Seal,
    ) -> Result<Vec<ManifestArtifact>> {
        Ok(self
            .imported_lineage(seal)?
            .map(|value| value.inputs)
            .unwrap_or_else(|| seal.inputs.clone()))
    }
    pub(in crate::generations) fn ancestry_inputs(
        &self,
        seal: &Seal,
    ) -> Result<Vec<ManifestArtifact>> {
        Ok(self
            .imported_lineage(seal)?
            .map(|value| value.ancestry_inputs)
            .unwrap_or_else(|| seal.inputs.clone()))
    }
    pub(in crate::generations) fn effective_dependencies(
        &self,
        seal: &Seal,
    ) -> Result<(bool, Vec<GenerationDependency>)> {
        if let Some(imported) = self.imported_lineage(seal)? {
            return Ok((true, imported.dependencies));
        }
        Ok((seal.lineage_version == 1, seal.dependencies.clone()))
    }
}
