use super::*;

impl Database {
    pub fn generation_inputs(
        &self,
        group: OutputGroup,
        subject: &str,
    ) -> Result<Option<Vec<ManifestArtifact>>> {
        if let Some(binding) = self.generation_input_binding(group, subject)? {
            return Ok(Some(binding.inputs));
        }
        self.current_seal(group, subject)?
            .map(|seal| self.effective_inputs(&seal))
            .transpose()
    }

    pub fn verified_retained_inputs(&self, run: &str) -> Result<Vec<ManifestArtifact>> {
        let retained = self.retained_inputs(run)?;
        if retained.is_empty() {
            return Ok(vec![]);
        }
        let ids = self.read(|connection| {
            let mut query = connection.prepare(
                "SELECT generation_id FROM artifact_generations WHERE project_id=?1
                AND output_group='quant' AND subject=?2 AND state='COMMITTED'",
            )?;
            Ok(query
                .query_map(params![self.project_id(), run], |row| {
                    row.get::<_, String>(0)
                })?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        let mut verified = Vec::new();
        for id in ids {
            let Some(seal) = self.committed_seal(&id)? else {
                continue;
            };
            for item in &retained {
                if quant_manifest_digest(&seal) == Some(item.generation_sha256.as_str())
                    && self.effective_inputs(&seal)?.contains(&item.artifact)
                    && !verified.contains(&item.artifact)
                {
                    verified.push(item.artifact.clone());
                }
            }
        }
        for item in &retained {
            if !verified.contains(&item.artifact)
                && self.binding_retention_verified(run, &item.generation_sha256, &item.artifact)?
            {
                verified.push(item.artifact.clone());
            }
        }
        Ok(verified)
    }

    pub fn generation_heads(&self) -> Result<Vec<GenerationHead>> {
        self.read(|connection| {
            let mut query = connection.prepare("SELECT h.generation_id, h.output_group, h.subject,
                CASE WHEN g.state = 'COMMITTED' AND g.project_id = h.project_id AND g.output_group = h.output_group
                    AND g.subject = h.subject THEN g.inventory_json END
                FROM generation_heads h LEFT JOIN artifact_generations g ON h.generation_id = g.generation_id
                WHERE h.project_id = ?1 ORDER BY h.output_group, h.subject")?;
            query.query_map(params![self.project_id()], |row| Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?,
                row.get::<_, String>(2)?, row.get::<_, Option<String>>(3)?)))?
                .map(|row| {
                    let (id, group, subject, json) = row?;
                    let group: OutputGroup = serde_json::from_value(serde_json::Value::String(group))?;
                    Ok(self.validate_seal(&id, group, &subject, json.as_deref())?.head)
                }).collect()
        })
    }

    pub fn generation_head(
        &self,
        group: OutputGroup,
        subject: &str,
    ) -> Result<Option<GenerationHead>> {
        Ok(self.current_seal(group, subject)?.map(|seal| seal.head))
    }

    pub(super) fn current_seal(&self, group: OutputGroup, subject: &str) -> Result<Option<Seal>> {
        self.read(|connection| {
            let stored: Option<(String, Option<String>)> = connection.query_row("SELECT h.generation_id,
                CASE WHEN g.state = 'COMMITTED' AND g.project_id = h.project_id AND g.output_group = h.output_group AND g.subject = h.subject THEN g.inventory_json END
                FROM generation_heads h LEFT JOIN artifact_generations g ON h.generation_id = g.generation_id
                WHERE h.project_id = ?1 AND h.output_group = ?2 AND h.subject = ?3",
                params![self.project_id(), group.as_str(), subject], |row| Ok((row.get(0)?, row.get(1)?))).optional()?;
            let Some((id, json)) = stored else { return Ok(None) };
            Ok(Some(self.validate_seal(&id, group, subject, json.as_deref())?))
        })
    }

    fn validate_seal(
        &self,
        id: &str,
        group: OutputGroup,
        subject: &str,
        json: Option<&str>,
    ) -> Result<Seal> {
        let json =
            json.context("current generation is incomplete or belongs to another output group")?;
        let seal: Seal = serde_json::from_str(json)?;
        let expected = fs::canonicalize(self.path().parent().context("project has no root")?)?
            .join(".rnaa_generations")
            .join(id)
            .join("payload");
        let mut unique = std::collections::BTreeSet::new();
        if uuid::Uuid::parse_str(id).is_err()
            || seal.schema_version != 1
            || seal.head.generation_id != id
            || seal.head.output_group != group
            || seal.head.subject != subject
            || seal.head.payload_path != expected
            || seal.files.is_empty()
            || seal.files.iter().any(|file| {
                !Path::new(&file.path).starts_with(&expected)
                    || !unique.insert(&file.path)
                    || file.checksum_type != "sha256"
                    || !valid_sha256(&file.checksum)
                    || Path::new(&file.path)
                        .components()
                        .any(|part| matches!(part, std::path::Component::ParentDir))
            })
            || !valid_sha256(&seal.head.request_sha256)
            || seal
                .inputs
                .iter()
                .any(|input| input.checksum_type != "sha256" || !valid_sha256(&input.checksum))
        {
            bail!("invalid committed generation inventory");
        }
        Ok(seal)
    }

    /// A published generation requires its durable seal witness and the exact
    /// complete file inventory to verify. Legacy files alone are not a generation.
    pub fn verify_generation(&self, group: OutputGroup, subject: &str) -> Result<bool> {
        let _pin = if self.generation_maintenance {
            None
        } else {
            Some(pins::read_lock(self.path())?)
        };
        let Some(seal) = self.current_seal(group, subject)? else {
            return Ok(false);
        };
        if self.generation_retired(&seal.head.generation_id)? {
            return Ok(false);
        }
        self.verify_seal(&seal)
    }

    pub fn verify_generation_head(&self, head: &GenerationHead) -> Result<bool> {
        let _pin = if self.generation_maintenance {
            None
        } else {
            Some(pins::read_lock(self.path())?)
        };
        if self.generation_retired(&head.generation_id)? {
            return Ok(false);
        }
        let Some(seal) = self.committed_seal(&head.generation_id)? else {
            return Ok(false);
        };
        if seal.head != *head {
            return Ok(false);
        }
        self.verify_seal(&seal)
    }

    pub(super) fn verify_seal(&self, seal: &Seal) -> Result<bool> {
        if !self.lineage_matches_catalog(seal)? {
            return Ok(false);
        }
        if !seal.head.payload_path.is_dir() {
            return Ok(false);
        }
        let witness = fs::read(
            seal.head
                .payload_path
                .parent()
                .context("generation has no root")?
                .join("seal.json"),
        )
        .ok()
        .and_then(|bytes| serde_json::from_slice::<Seal>(&bytes).ok());
        if witness.as_ref() != Some(seal) {
            return Ok(false);
        }
        let actual = inventory(&seal.head.payload_path, false)?;
        let expected: std::collections::BTreeMap<_, _> =
            seal.files.iter().map(|file| (&file.path, file)).collect();
        let available: std::collections::BTreeMap<_, _> =
            actual.iter().map(|file| (&file.path, file)).collect();
        // Pruning can explain absence only. It must never hide an altered file,
        // unexpected payload, or missing report/manifest.
        for file in &actual {
            if expected.get(&file.path).copied() != Some(file) {
                return Ok(false);
            }
        }
        for file in &seal.files {
            if !available.contains_key(&file.path) && !self.pruning_authorized(seal, file)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(super) fn committed_seal(&self, id: &str) -> Result<Option<Seal>> {
        self.read(|connection| {
            let stored: Option<(String, String, String)> = connection
                .query_row(
                    "SELECT output_group, subject, inventory_json FROM artifact_generations
                 WHERE project_id = ?1 AND generation_id = ?2 AND state = 'COMMITTED'",
                    params![self.project_id(), id],
                    |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
                )
                .optional()?;
            stored
                .map(|(group, subject, json)| {
                    let group = serde_json::from_value(serde_json::Value::String(group))?;
                    self.validate_seal(id, group, &subject, Some(&json))
                })
                .transpose()
        })
    }

    fn pruning_authorized(&self, producer: &Seal, file: &ManifestArtifact) -> Result<bool> {
        if self.file_trash_pruning_authorized(producer, file)? {
            return Ok(true);
        }
        if !matches!(
            producer.head.output_group,
            OutputGroup::Download | OutputGroup::Preprocess
        ) {
            return Ok(false);
        }
        let Some(expected) = self.producer_fastq(producer, file)? else {
            return Ok(false);
        };
        let receipts = self.read(|connection| {
            let mut query = connection.prepare(
                "SELECT consumer_generation_id, quant_manifest_sha256, artifact_json
                FROM generation_pruned_files WHERE project_id = ?1 AND run_accession = ?2
                AND producer_generation_id = ?3 AND path = ?4",
            )?;
            Ok(query
                .query_map(
                    params![
                        self.project_id(),
                        producer.head.subject,
                        producer.head.generation_id,
                        file.path
                    ],
                    |row| {
                        Ok((
                            row.get::<_, String>(0)?,
                            row.get::<_, String>(1)?,
                            row.get::<_, String>(2)?,
                        ))
                    },
                )?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        for (consumer, manifest, json) in receipts {
            let artifact: ManifestArtifact = serde_json::from_str(&json)?;
            if artifact != expected {
                continue;
            }
            let Some(consumer) = self.committed_seal(&consumer)? else {
                continue;
            };
            if consumer.head.output_group == OutputGroup::Quant
                && consumer.head.subject == producer.head.subject
                && self.effective_inputs(&consumer)?.contains(&artifact)
                && quant_manifest_digest(&consumer) == Some(manifest.as_str())
            {
                return Ok(true);
            }
        }
        self.binding_pruning_authorized(&producer.head, &expected, None)
    }

    pub(super) fn producer_fastq(
        &self,
        producer: &Seal,
        file: &ManifestArtifact,
    ) -> Result<Option<ManifestArtifact>> {
        let records = self.read(|connection| {
            let json: String = connection.query_row(
                "SELECT artifacts_json FROM artifact_generations
                WHERE project_id = ?1 AND generation_id = ?2 AND state = 'COMMITTED'",
                params![self.project_id(), producer.head.generation_id],
                |row| row.get(0),
            )?;
            Ok(serde_json::from_str::<Vec<ArtifactRecord>>(&json)?)
        })?;
        Ok(records
            .into_iter()
            .filter(|record| {
                record.project_id == self.project_id()
                    && record.run_accession.as_deref() == Some(&producer.head.subject)
                    && prunable_fastq(producer.head.output_group, record.kind.as_str())
            })
            .map(|record| ManifestArtifact {
                kind: record.kind.as_str().into(),
                path: record.path,
                checksum_type: record.checksum_type,
                checksum: record.checksum,
                bytes: record.bytes,
            })
            .find(|record| same_file(record, file)))
    }

    /// Persist removal evidence before unlinking. A file inside a generation also
    /// requires a committed consuming quant generation with the exact input digest.
    /// Reports and manifests can never receive a pruning authorization.
    pub(super) fn generation_artifacts(
        &self,
        group: OutputGroup,
        subject: &str,
    ) -> Result<Option<Vec<ArtifactRecord>>> {
        self.read(|connection| {
            let records: Option<String> = connection.query_row("SELECT g.artifacts_json FROM artifact_generations g JOIN generation_heads h
                ON h.generation_id = g.generation_id WHERE h.project_id = ?1 AND h.output_group = ?2 AND h.subject = ?3 AND g.state = 'COMMITTED'",
                params![self.project_id(), group.as_str(), subject], |row| row.get(0)).optional()?;
            records.map(|records| Ok(serde_json::from_str(&records)?)).transpose()
        })
    }

    pub fn recover_generations(&self) -> Result<usize> {
        self.write(|connection| {
            Ok(connection.execute(
                "UPDATE artifact_generations SET state = 'ABANDONED'
            WHERE state IN ('BUILDING','SEALED') AND EXISTS(SELECT 1 FROM job_attempts j
                WHERE j.attempt_id = artifact_generations.attempt_id AND j.state <> 'RUNNING')",
                [],
            )?)
        })
    }
}
