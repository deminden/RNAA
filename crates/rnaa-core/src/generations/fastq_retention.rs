use super::*;

impl Database {
    pub fn record_fastq_retention(
        &self,
        run: &str,
        manifest: &str,
        artifact: &ManifestArtifact,
    ) -> Result<()> {
        if !is_trimmed_fastq(&artifact.kind)
            && !matches!(
                artifact.kind.as_str(),
                "FASTQ_R1" | "FASTQ_R2" | "FASTQ_SINGLE"
            )
        {
            bail!("only FASTQ artifacts can be pruned");
        }
        let root = fs::canonicalize(self.path().parent().context("project has no root")?)?;
        let store = root.join(".rnaa_generations");
        let source = Path::new(&artifact.path);
        let present = artifact.matches_file()?;
        if source.exists() && !present {
            bail!("retention input changed: {}", artifact.path);
        }
        if !present
            && !self
                .quant_retention_evidence(run)?
                .iter()
                .any(|item| item.generation_sha256 == manifest && item.artifact == *artifact)
        {
            bail!(
                "retention input changed or disappeared without evidence: {}",
                artifact.path
            );
        }
        let Ok(relative) = source.strip_prefix(&store) else {
            return self.record_retained_input(run, manifest, artifact);
        };
        let id = relative
            .components()
            .next()
            .and_then(|part| part.as_os_str().to_str())
            .context("generation input has no ID")?;
        let producer = self
            .committed_seal(id)?
            .context("retention producer is not committed")?;
        if !matches!(
            producer.head.output_group,
            OutputGroup::Download | OutputGroup::Preprocess
        ) || producer.head.subject != run
            || !prunable_fastq(producer.head.output_group, &artifact.kind)
            || !producer.files.iter().any(|file| same_file(file, artifact))
            || self.producer_fastq(&producer, artifact)?.as_ref() != Some(artifact)
            || !self.verify_seal(&producer)?
        {
            bail!("retention does not match a verified FASTQ producer");
        }
        let consumer = self
            .current_seal(OutputGroup::Quant, run)?
            .context("generation pruning requires a committed quantification consumer")?;
        let binding = self.generation_input_binding(OutputGroup::Quant, run)?;
        if !self.verify_seal(&consumer)?
            || (!self.effective_inputs(&consumer)?.contains(artifact)
                && !binding
                    .as_ref()
                    .is_some_and(|binding| binding.inputs.contains(artifact)))
            || quant_manifest_digest(&consumer) != Some(manifest)
        {
            bail!("retention does not match the consuming quantification generation");
        }
        self.atomic(|db| {
            if db.generation_head(OutputGroup::Quant, run)?.as_ref() != Some(&consumer.head) {
                bail!("quantification changed before retention authorization");
            }
            if db.generation_input_binding(OutputGroup::Quant, run)? != binding {
                bail!("quantification input binding changed before retention authorization");
            }
            db.record_retained_input(run, manifest, artifact)?;
            if !db.effective_inputs(&consumer)?.contains(artifact) {
                return db.record_binding_pruning(&producer.head, binding.as_ref().context("retention lacks its input binding")?, manifest, artifact);
            }
            db.write(|connection| {
                connection.execute("INSERT INTO generation_pruned_files(project_id, run_accession, producer_generation_id,
                    consumer_generation_id, quant_manifest_sha256, path, artifact_json, recorded_at)
                    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8) ON CONFLICT DO NOTHING",
                    params![db.project_id(), run, producer.head.generation_id, consumer.head.generation_id, manifest,
                        artifact.path, serde_json::to_string(artifact)?, now_rfc3339()])?;
                Ok(())
            })
        })
    }
}
