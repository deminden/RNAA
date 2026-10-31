//! Immutable use-time bindings for scientifically identical files at new paths.
mod retention;
use super::{GenerationDependency, GenerationHead, OutputGroup, VerifiedGenerationReuse};
use crate::{Database, ManifestArtifact, fingerprint::InputIdentity};
use anyhow::{Context, Result, bail};
use rusqlite::{OptionalExtension, params};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenerationInputBinding {
    pub schema_version: u32,
    pub binding_id: String,
    pub generation_id: String,
    pub request_sha256: String,
    pub attempt_id: String,
    pub inputs: Vec<ManifestArtifact>,
    pub dependencies: Vec<GenerationDependency>,
}
impl GenerationInputBinding {
    pub(in crate::generations) fn digest(&self) -> Result<String> {
        Ok(hex::encode(Sha256::digest(serde_json::to_vec(self)?)))
    }
}

impl Database {
    /// Bind current physical witnesses to an old request without editing its
    /// original manifest. Exact semantic identities and verified owners are required.
    pub fn bind_generation_reuse(
        &self,
        mut reuse: VerifiedGenerationReuse,
        inputs: &[ManifestArtifact],
    ) -> Result<VerifiedGenerationReuse> {
        self.check_job()?;
        let head = reuse.head();
        let binding = GenerationInputBinding {
            schema_version: 1,
            binding_id: uuid::Uuid::new_v4().to_string(),
            generation_id: head.generation_id.clone(),
            request_sha256: head.request_sha256.clone(),
            attempt_id: self
                .job
                .as_ref()
                .context("input binding requires an owned attempt")?
                .attempt_id
                .clone(),
            inputs: inputs.to_vec(),
            dependencies: self.resolve_generation_dependencies(inputs, &head.payload_path)?,
        };
        self.verify_input_binding(&binding)?;
        reuse.input_binding = Some(binding);
        Ok(reuse)
    }

    /// Rehash bound producer payloads and validate every missing-input receipt.
    /// Call outside writer transactions; catalog-only validation fences selection.
    pub fn verify_input_binding(&self, binding: &GenerationInputBinding) -> Result<()> {
        let _pin = if self.generation_maintenance {
            None
        } else {
            Some(super::pins::read_lock(self.path())?)
        };
        let head = self
            .committed_generation_head(&binding.generation_id)?
            .context("binding consumer disappeared")?;
        self.verify_input_binding_catalog(binding)?;
        for dependency in &binding.dependencies {
            let producer = self
                .committed_generation_head(&dependency.producer_generation_id)?
                .context("input producer disappeared")?;
            if !self.verify_generation_head(&producer)? {
                bail!("reuse input producer failed verification");
            }
        }
        for input in &binding.inputs {
            if input.matches_file()? {
                continue;
            }
            if Path::new(&input.path).exists() {
                bail!("reuse input changed: {}", input.path);
            }
            if binding
                .dependencies
                .iter()
                .any(|dependency| dependency.input == *input)
            {
                continue;
            }
            let retained = self.verified_retained_inputs(&head.subject)?;
            if !retained.contains(input) {
                bail!(
                    "reuse input disappeared without producer/retention evidence: {}",
                    input.path
                );
            }
        }
        Ok(())
    }

    pub(super) fn verify_input_binding_catalog(
        &self,
        binding: &GenerationInputBinding,
    ) -> Result<()> {
        let consumer = self
            .committed_seal(&binding.generation_id)?
            .context("input binding consumer is not committed")?;
        let request = self
            .scientific_request(&binding.request_sha256)?
            .context("input binding lacks its full scientific request")?;
        let mut identities = binding
            .inputs
            .iter()
            .map(InputIdentity::from)
            .collect::<Vec<_>>();
        identities.sort();
        let mut expected = request.inputs.clone();
        expected.sort();
        let mut unique = std::collections::BTreeSet::new();
        if binding.schema_version != 1
            || uuid::Uuid::parse_str(&binding.binding_id).is_err()
            || binding.request_sha256 != consumer.head.request_sha256
            || request.stage != consumer.head.stage
            || !self.binding_request_inputs_match(
                binding,
                &identities,
                &expected,
                &request.stage,
            )?
            || binding.inputs.is_empty()
            || binding.inputs.iter().any(|input| {
                input.checksum_type != "sha256" || !unique.insert((&input.kind, &input.path))
            })
            || binding.dependencies.iter().any(|dependency| {
                !binding.inputs.contains(&dependency.input)
                    || dependency.producer_generation_id == binding.generation_id
            })
            || binding.inputs.iter().any(|input| {
                Path::new(&input.path)
                    .components()
                    .any(|part| part.as_os_str() == ".rnaa_generations")
                    && binding
                        .dependencies
                        .iter()
                        .filter(|dependency| dependency.input == *input)
                        .count()
                        != 1
            })
        {
            bail!("input binding differs from its committed scientific request");
        }
        self.verify_dependency_catalog(&binding.dependencies)
    }

    fn binding_request_inputs_match(
        &self,
        binding: &GenerationInputBinding,
        actual: &[InputIdentity],
        expected: &[InputIdentity],
        stage: &str,
    ) -> Result<bool> {
        if actual == expected {
            return Ok(true);
        }
        let raw = |role: &str| matches!(role, "FASTQ_SINGLE" | "FASTQ_R1" | "FASTQ_R2");
        let trimmed = |role: &str| {
            matches!(
                role,
                "FASTQ_TRIMMED_SINGLE" | "FASTQ_TRIMMED_R1" | "FASTQ_TRIMMED_R2"
            )
        };
        // A fit made from retained trimmed reads can omit originals from its
        // scientific key. Current originals remain retention/source observations
        // only when their entire identity equals that trimmed producer's inputs.
        if stage != "quant"
            || expected.iter().any(|input| raw(&input.role))
            || !expected.iter().any(|input| trimmed(&input.role))
        {
            return Ok(false);
        }
        let projected = actual
            .iter()
            .filter(|input| !raw(&input.role))
            .cloned()
            .collect::<Vec<_>>();
        if projected != expected {
            return Ok(false);
        }
        let owners = binding
            .dependencies
            .iter()
            .filter(|dependency| trimmed(&dependency.input.kind))
            .map(|dependency| &dependency.producer_generation_id)
            .collect::<std::collections::BTreeSet<_>>();
        if owners.len() != 1 {
            return Ok(false);
        }
        let producer = self
            .committed_seal(owners.first().context("trimmed input lacks its producer")?)?
            .context("trimmed input producer disappeared")?;
        if producer.head.output_group != OutputGroup::Preprocess {
            return Ok(false);
        }
        let mut originals = producer
            .inputs
            .iter()
            .filter(|input| raw(&input.kind))
            .map(InputIdentity::from)
            .collect::<Vec<_>>();
        originals.sort();
        let observed = actual
            .iter()
            .filter(|input| raw(&input.role))
            .cloned()
            .collect::<Vec<_>>();
        Ok(!originals.is_empty() && originals == observed)
    }

    pub(super) fn select_input_binding(
        &self,
        head: &GenerationHead,
        binding: Option<&GenerationInputBinding>,
    ) -> Result<()> {
        self.write(|conn| {
            conn.execute("DELETE FROM generation_binding_heads WHERE project_id=?1 AND output_group=?2 AND subject=?3",
                params![self.project_id(),head.output_group.as_str(),head.subject])?;
            if let Some(binding) = binding {
                if self.job.as_ref().map(|job| job.attempt_id.as_str()) != Some(binding.attempt_id.as_str())
                    || binding.generation_id != head.generation_id { bail!("input binding does not belong to its reuse attempt"); }
                conn.execute("INSERT INTO generation_input_bindings(binding_id,project_id,generation_id,attempt_id,binding_sha256,binding_json)
                    VALUES (?1,?2,?3,?4,?5,?6)",params![binding.binding_id,self.project_id(),binding.generation_id,binding.attempt_id,binding.digest()?,serde_json::to_string(binding)?])?;
                conn.execute("INSERT INTO generation_binding_heads(project_id,output_group,subject,binding_id) VALUES (?1,?2,?3,?4)",
                    params![self.project_id(),head.output_group.as_str(),head.subject,binding.binding_id])?;
            }
            Ok(())
        })
    }

    pub fn generation_input_binding(
        &self,
        group: OutputGroup,
        subject: &str,
    ) -> Result<Option<GenerationInputBinding>> {
        let id: Option<String> = self.read(|conn| {
            Ok(conn
                .query_row(
                    "SELECT binding_id FROM generation_binding_heads
            WHERE project_id=?1 AND output_group=?2 AND subject=?3",
                    params![self.project_id(), group.as_str(), subject],
                    |row| row.get(0),
                )
                .optional()?)
        })?;
        id.map(|id| {
            let binding = self
                .read_input_binding(&id)?
                .context("selected input binding disappeared")?;
            let head = self
                .generation_head(group, subject)?
                .context("input binding has no selected generation")?;
            if binding.generation_id != head.generation_id {
                bail!("selected input binding belongs to another generation");
            }
            Ok(binding)
        })
        .transpose()
    }

    pub(super) fn read_input_binding(&self, id: &str) -> Result<Option<GenerationInputBinding>> {
        let stored: Option<(String, String, String, String)> = self.read(|conn| {
            Ok(conn
                .query_row(
                    "SELECT binding_json,binding_sha256,generation_id,attempt_id
            FROM generation_input_bindings WHERE project_id=?1 AND binding_id=?2",
                    params![self.project_id(), id],
                    |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?)),
                )
                .optional()?)
        })?;
        stored
            .map(|(json, sha, generation, attempt)| {
                let binding: GenerationInputBinding = serde_json::from_str(&json)?;
                if binding.binding_id != id
                    || binding.generation_id != generation
                    || binding.attempt_id != attempt
                    || binding.digest()? != sha
                {
                    bail!("input binding failed its immutable receipt verification");
                }
                self.verify_input_binding_catalog(&binding)?;
                Ok(binding)
            })
            .transpose()
    }
}
