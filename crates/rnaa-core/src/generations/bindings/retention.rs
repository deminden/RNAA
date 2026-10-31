//! Retention evidence for inputs adopted by a later reuse attempt.
use super::*;
use crate::generations::quant_manifest_digest;
impl Database {
    /// Preserve strict original-manifest receipts while allowing a verified
    /// current binding to refer to inputs pruned by another committed count result.
    pub fn quant_retention_evidence(
        &self,
        run: &str,
    ) -> Result<Vec<crate::fingerprint::RetainedInput>> {
        let mut evidence = self.retained_inputs(run)?;
        if let Some(binding) = self.generation_input_binding(OutputGroup::Quant, run)? {
            let consumer = self
                .committed_seal(&binding.generation_id)?
                .context("binding consumer disappeared")?;
            let manifest = quant_manifest_digest(&consumer)
                .context("binding consumer lacks a quant manifest")?;
            for input in self.verified_retained_inputs(run)? {
                if binding.inputs.contains(&input)
                    && !evidence
                        .iter()
                        .any(|item| item.generation_sha256 == manifest && item.artifact == input)
                {
                    evidence.push(crate::fingerprint::RetainedInput {
                        generation_sha256: manifest.into(),
                        artifact: input,
                    });
                }
            }
        }
        Ok(evidence)
    }

    pub(in crate::generations) fn binding_pruning_authorized(
        &self,
        producer: &GenerationHead,
        artifact: &ManifestArtifact,
        required_manifest: Option<&str>,
    ) -> Result<bool> {
        let receipts = self.read(|conn| {
            let mut query = conn.prepare("SELECT binding_id,binding_sha256,quant_manifest_sha256,artifact_json FROM generation_binding_pruning
                WHERE producer_generation_id=?1 AND path=?2")?;
            Ok(query.query_map(params![producer.generation_id,artifact.path], |row| Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,String>(2)?,row.get::<_,String>(3)?)))?
                .collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        for (id, sha, manifest, json) in receipts {
            if required_manifest.is_some_and(|required| required != manifest) {
                continue;
            }
            if serde_json::from_str::<ManifestArtifact>(&json)? != *artifact {
                continue;
            }
            let Some(binding) = self.read_input_binding(&id)? else {
                continue;
            };
            let Some(consumer) = self.committed_seal(&binding.generation_id)? else {
                continue;
            };
            if binding.digest()? == sha
                && consumer.head.output_group == OutputGroup::Quant
                && consumer.head.subject == producer.subject
                && quant_manifest_digest(&consumer) == Some(manifest.as_str())
                && binding.inputs.contains(artifact)
                && binding.dependencies.iter().any(|dependency| {
                    dependency.producer_generation_id == producer.generation_id
                        && dependency.input == *artifact
                })
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub(in crate::generations) fn binding_retention_verified(
        &self,
        run: &str,
        manifest: &str,
        artifact: &ManifestArtifact,
    ) -> Result<bool> {
        if self.file_trash_binding_retention_verified(run, manifest, artifact)? {
            return Ok(true);
        }
        let owners = self.read(|conn| {
            let mut query=conn.prepare("SELECT g.generation_id FROM generation_files f JOIN artifact_generations g ON g.generation_id=f.generation_id
                WHERE g.project_id=?1 AND g.subject=?2 AND g.state='COMMITTED' AND f.path=?3 AND f.sha256=?4 AND f.bytes=?5")?;
            Ok(query.query_map(params![self.project_id(),run,artifact.path,artifact.checksum,i64::try_from(artifact.bytes)?],|row| row.get::<_,String>(0))?
                .collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        for id in owners {
            let producer = self
                .committed_generation_head(&id)?
                .context("retained binding producer disappeared")?;
            if self.binding_pruning_authorized(&producer, artifact, Some(manifest))? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub(in crate::generations) fn record_binding_pruning(
        &self,
        producer: &GenerationHead,
        binding: &GenerationInputBinding,
        manifest: &str,
        artifact: &ManifestArtifact,
    ) -> Result<()> {
        self.verify_input_binding_catalog(binding)?;
        self.write(|conn| {
            conn.execute("INSERT INTO generation_binding_pruning(producer_generation_id,path,binding_id,binding_sha256,quant_manifest_sha256,artifact_json)
                VALUES (?1,?2,?3,?4,?5,?6) ON CONFLICT DO NOTHING",params![producer.generation_id,artifact.path,binding.binding_id,binding.digest()?,manifest,serde_json::to_string(artifact)?])?;
            Ok(())
        })
    }
}
