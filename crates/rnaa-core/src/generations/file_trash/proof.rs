use super::*;

impl Database {
    pub(super) fn incomplete_cleanup_groups(
        &self,
        values: &[Descriptor],
    ) -> Result<BTreeMap<String, String>> {
        let family = |kind: &str| match kind {
            "FASTQ_R1" | "FASTQ_R2" | "FASTQ_SINGLE" => Some(0),
            "FASTQ_TRIMMED_R1" | "FASTQ_TRIMMED_R2" | "FASTQ_TRIMMED_SINGLE" => Some(1),
            _ => None,
        };
        let mut groups = BTreeMap::<(String, u8), Vec<&Descriptor>>::new();
        for value in values {
            if let Some(family) = family(&value.source.kind) {
                groups
                    .entry((value.run.clone(), family))
                    .or_default()
                    .push(value);
            }
        }
        let mut held = BTreeMap::new();
        for ((run, group), members) in groups {
            let inputs =
                if let Some(binding) = self.generation_input_binding(OutputGroup::Quant, &run)? {
                    binding.inputs
                } else {
                    self.effective_inputs(
                        &self
                            .committed_seal(&members[0].consumer)?
                            .context("cleanup count consumer disappeared")?,
                    )?
                };
            let mut expected = inputs
                .into_iter()
                .filter(|input| family(&input.kind) == Some(group))
                .collect::<Vec<_>>();
            let mut actual = members
                .iter()
                .map(|value| value.source.clone())
                .collect::<Vec<_>>();
            expected.sort_by(|a, b| (&a.kind, &a.path).cmp(&(&b.kind, &b.path)));
            actual.sort_by(|a, b| (&a.kind, &a.path).cmp(&(&b.kind, &b.path)));
            if expected != actual {
                held.insert(run, "registered FASTQs are not the complete original/trimmed input group consumed by the selected counts".into());
            }
        }
        Ok(held)
    }
    pub(super) fn verify_trash_descriptor(&self, value: &Descriptor) -> Result<()> {
        let root = fs::canonicalize(self.path().parent().context("project has no root")?)?;
        let source = Path::new(&value.canonical_source);
        if value.schema_version != 1
            || value.project_id != self.project_id()
            || uuid::Uuid::parse_str(&value.trash_id).is_err()
            || uuid::Uuid::parse_str(&value.operation_id).is_err()
            || value.record.project_id != self.project_id()
            || value.record.run_accession.as_deref() != Some(&value.run)
            || !selected(value.mode, value.record.kind)
            || value.source.kind != value.record.kind.as_str()
            || value.source.path != value.record.path
            || value.source.checksum_type != "sha256"
            || value.source.checksum.len() != 64
            || !value
                .source
                .checksum
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit())
            || value.purge_after_ms < value.created_at_ms
            || value.source.bytes != value.record.bytes
            || !source.starts_with(&root)
            || source
                .components()
                .any(|part| matches!(part, std::path::Component::ParentDir))
            || ["trash", ".rnaa_support", ".rnaa_transfers"]
                .iter()
                .any(|name| source.starts_with(root.join(name)))
        {
            bail!("invalid file cleanup ownership descriptor");
        }
        let consumer = self
            .committed_seal(&value.consumer)?
            .context("file cleanup count consumer is missing")?;
        if consumer.head.output_group != OutputGroup::Quant
            || consumer.head.subject != value.run
            || digest(serde_json::to_vec(&consumer)?) != value.consumer_seal
            || quant_manifest_digest(&consumer) != Some(value.consumer_manifest.as_str())
            || !self.lineage_matches_catalog(&consumer)?
        {
            bail!("file cleanup count consumer proof changed");
        }
        if fastq(value.record.kind) {
            if let Some(binding) = &value.binding {
                let binding = self
                    .read_input_binding(binding)?
                    .context("file cleanup input binding is missing")?;
                if binding.generation_id != value.consumer
                    || Some(binding.digest()?) != value.binding_sha256
                    || !binding.inputs.contains(&value.source)
                {
                    bail!("file cleanup binding differs from consumed inputs");
                }
            } else if value.binding_sha256.is_some()
                || !self.effective_inputs(&consumer)?.contains(&value.source)
            {
                bail!("file cleanup lacks exact count input evidence");
            }
        } else {
            if value.record.kind != ArtifactKind::Sra
                || value.binding.is_some()
                || value.binding_sha256.is_some()
            {
                bail!("archive cleanup requires an unused archive beside verified counts");
            }
            for input in &self.effective_inputs(&consumer)? {
                if canonical_missing_path(Path::new(&input.path))? == source {
                    bail!("archive cleanup cannot remove a consumed count input");
                }
            }
        }
        match (&value.producer, &value.producer_seal) {
            (Some(id), Some(hash)) => {
                let producer = self
                    .committed_seal(id)?
                    .context("file cleanup producer is missing")?;
                if producer.head.subject != value.run
                    || digest(serde_json::to_vec(&producer)?) != *hash
                    || !self.lineage_matches_catalog(&producer)?
                    || !producer.files.iter().any(|file| {
                        file.path == value.canonical_source
                            && file.checksum == value.source.checksum
                            && file.bytes == value.source.bytes
                    })
                    || !(if value.record.kind == ArtifactKind::Sra {
                        producer.head.output_group == OutputGroup::Download
                    } else {
                        super::super::prunable_fastq(
                            producer.head.output_group,
                            value.source.kind.as_str(),
                        )
                    })
                {
                    bail!("file cleanup producer differs from sealed input ownership");
                }
                let records: String = self.read(|conn| {
                    Ok(conn.query_row(
                        "SELECT artifacts_json FROM artifact_generations WHERE generation_id=?1",
                        [id],
                        |row| row.get(0),
                    )?)
                })?;
                if !serde_json::from_str::<Vec<ArtifactRecord>>(&records)?
                    .iter()
                    .any(|record| {
                        record.path == value.canonical_source
                            && record.kind == value.record.kind
                            && record.checksum == value.source.checksum
                            && record.bytes == value.source.bytes
                    })
                {
                    bail!("file cleanup producer kind or registration changed");
                }
            }
            (None, None) => {
                if source.starts_with(root.join(".rnaa_generations")) {
                    bail!("managed cleanup file has no sealed producer");
                }
            }
            _ => bail!("incomplete file cleanup producer proof"),
        }
        Ok(())
    }
    pub(super) fn prepare_trash_descriptor(
        &self,
        record: &ArtifactRecord,
        operation: &str,
        mode: CleanupMode,
        days: u64,
        multiplicity: &BTreeMap<String, usize>,
    ) -> Result<Descriptor> {
        let run = record
            .run_accession
            .as_deref()
            .context("cleanup input has no run")?;
        let current = self.get_run(run)?.context("cleanup run is missing")?;
        let canonical = canonical_missing_path(Path::new(&record.path))?;
        if multiplicity.get(canonical.to_string_lossy().as_ref()) != Some(&1) {
            bail!("source has multiple active registrations; ownership is not exclusive");
        }
        let metadata = fs::symlink_metadata(&canonical)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            if metadata.nlink() != 1 {
                bail!("cleanup requires an exclusively owned file with no hardlink aliases");
            }
        }
        if !metadata.is_file() {
            bail!("cleanup source must be a regular file");
        }
        let source = ManifestArtifact::snapshot(record.kind.as_str(), Path::new(&record.path))?;
        if source.bytes != record.bytes
            || match record.checksum_type.as_str() {
                "sha256" => source.checksum != record.checksum,
                "md5" => !crate::util::compute_md5(Path::new(&record.path))?
                    .eq_ignore_ascii_case(&record.checksum),
                _ => true,
            }
        {
            bail!("cleanup input changed since registration");
        }
        let consumer = self
            .current_seal(OutputGroup::Quant, run)?
            .context("cleanup requires a committed native count result")?;
        if self.generation_retired(&consumer.head.generation_id)? || !self.verify_seal(&consumer)? {
            bail!("count generation failed cleanup verification");
        }
        let manifest = quant_manifest_digest(&consumer)
            .context("count result lacks a sealed manifest")?
            .to_string();
        let (binding, binding_sha256) =
            if fastq(record.kind) && !self.effective_inputs(&consumer)?.contains(&source) {
                let binding = self
                    .generation_input_binding(OutputGroup::Quant, run)?
                    .context("cleanup input is absent from count provenance")?;
                self.verify_input_binding(&binding)?;
                if !binding.inputs.contains(&source) {
                    bail!("cleanup input is absent from the selected count binding");
                }
                (Some(binding.binding_id.clone()), Some(binding.digest()?))
            } else {
                (None, None)
            };
        let owners=self.read(|conn| {
            let mut query=conn.prepare("SELECT f.generation_id FROM generation_files f JOIN artifact_generations g ON g.generation_id=f.generation_id WHERE g.project_id=?1 AND g.state='COMMITTED' AND f.path=?2 AND f.sha256=?3 AND f.bytes=?4")?;
            Ok(query.query_map(params![self.project_id(),canonical.to_string_lossy(),source.checksum,i64::try_from(source.bytes)?],|row|row.get::<_,String>(0))?.collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        if owners.len() > 1 {
            bail!("cleanup input has ambiguous producer ownership");
        }
        let producer = owners.first().cloned();
        let producer_seal = if let Some(id) = &producer {
            let producer = self
                .committed_seal(id)?
                .context("cleanup producer disappeared")?;
            if self.generation_retired(id)? || !self.verify_seal(&producer)? {
                bail!("cleanup producer failed verification");
            }
            Some(digest(serde_json::to_vec(&producer)?))
        } else {
            None
        };
        let created = now_ms()?;
        let expiry = created
            .checked_add(i64::try_from(
                days.checked_mul(86_400_000)
                    .context("file trash retention is too large")?,
            )?)
            .context("file trash expiry overflow")?;
        let value = Descriptor {
            schema_version: 1,
            trash_id: uuid::Uuid::new_v4().to_string(),
            operation_id: operation.into(),
            project_id: self.project_id().into(),
            run: run.into(),
            layout: current.library_layout,
            source_selection: digest(serde_json::to_vec(&current.remote_files)?),
            record: record.clone(),
            source,
            canonical_source: canonical.display().to_string(),
            producer,
            producer_seal,
            consumer: consumer.head.generation_id.clone(),
            consumer_seal: digest(serde_json::to_vec(&consumer)?),
            consumer_manifest: manifest,
            binding,
            binding_sha256,
            mode,
            created_at_ms: created,
            purge_after_ms: expiry,
        };
        self.verify_trash_descriptor(&value)?;
        Ok(value)
    }
    pub(in crate::generations) fn file_trash_pruning_authorized(
        &self,
        producer: &Seal,
        file: &ManifestArtifact,
    ) -> Result<bool> {
        for row in self.file_trash_rows_filtered(None, Some(&file.path))? {
            let value = &row.descriptor;
            if row.active_absence()
                && value.producer.as_deref() == Some(&producer.head.generation_id)
                && value.canonical_source == file.path
                && value.source.checksum == file.checksum
                && value.source.bytes == file.bytes
            {
                return Ok(true);
            }
        }
        Ok(false)
    }
    pub(crate) fn file_trash_retained_inputs(
        &self,
        run: &str,
    ) -> Result<Vec<crate::fingerprint::RetainedInput>> {
        Ok(self
            .file_trash_rows_filtered(Some(run), None)?
            .into_iter()
            .filter(|row| {
                row.active_absence()
                    && row.descriptor.run == run
                    && fastq(row.descriptor.record.kind)
            })
            .map(|row| crate::fingerprint::RetainedInput {
                generation_sha256: row.descriptor.consumer_manifest,
                artifact: row.descriptor.source,
            })
            .collect())
    }
    pub(in crate::generations) fn file_trash_binding_retention_verified(
        &self,
        run: &str,
        manifest: &str,
        artifact: &ManifestArtifact,
    ) -> Result<bool> {
        Ok(self
            .file_trash_rows_filtered(Some(run), None)?
            .iter()
            .any(|row| {
                row.active_absence()
                    && row.descriptor.run == run
                    && row.descriptor.consumer_manifest == manifest
                    && row.descriptor.binding.is_some()
                    && row.descriptor.source == *artifact
            }))
    }
    pub(in crate::generations) fn file_trash_generation_links(
        &self,
    ) -> Result<Vec<GenerationLink>> {
        Ok(self
            .file_trash_rows()?
            .into_iter()
            .map(|row| {
                let pin = matches!(
                    row.state.as_str(),
                    "MOVING" | "TRASHED" | "RESTORING" | "PURGING"
                );
                let absence = row.active_absence();
                GenerationLink {
                    producer: row.descriptor.producer,
                    consumer: row.descriptor.consumer,
                    restorable: pin,
                    absence,
                }
            })
            .collect())
    }
}
