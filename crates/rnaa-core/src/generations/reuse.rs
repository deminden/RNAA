//! Historical lookup and fenced selection of an already committed payload.
use super::{GenerationHead, OutputGroup, Seal};
use crate::{ArtifactRecord, Database, jobs::JobOutcome};
use anyhow::{Context, Result, bail};
use rusqlite::{OptionalExtension, params};

/// Prepared outside the writer transaction. Callers must additionally validate the
/// stage's current scientific request and admission rules, just as for a fresh fit.
/// Private fields prevent constructing a token from unverified catalog metadata.
pub struct VerifiedGenerationReuse {
    _pin: std::sync::Arc<super::GenerationStorePin>,
    project_id: String,
    seal: Seal,
    records: Vec<ArtifactRecord>,
    records_json: String,
    expected_head: Option<GenerationHead>,
    pub(super) input_binding: Option<super::GenerationInputBinding>,
}

impl VerifiedGenerationReuse {
    pub fn head(&self) -> &GenerationHead {
        &self.seal.head
    }

    pub fn artifact_count(&self) -> usize {
        self.records.len()
    }
}

impl Database {
    /// All committed request variants for one output subject, including results
    /// that are no longer selected. Payload validity requires explicit verification.
    pub fn committed_generation_history(
        &self,
        group: OutputGroup,
        subject: &str,
    ) -> Result<Vec<GenerationHead>> {
        let ids = self.read(|conn| {
            let mut query = conn.prepare(
                "SELECT generation_id FROM artifact_generations
                WHERE project_id=?1 AND output_group=?2 AND subject=?3 AND state='COMMITTED'
                AND NOT EXISTS(SELECT 1 FROM generation_retirements r WHERE r.generation_id=artifact_generations.generation_id)
                ORDER BY committed_at DESC,generation_id DESC",
            )?;
            Ok(query
                .query_map(params![self.project_id(), group.as_str(), subject], |row| {
                    row.get::<_, String>(0)
                })?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        ids.iter()
            .map(|id| {
                Ok(self
                    .committed_seal(id)?
                    .context("committed history disappeared")?
                    .head)
            })
            .collect()
    }

    /// Exact scoped request lookup, preferring the selected result and then the
    /// most recently committed candidate. This is discovery, not validity proof.
    pub fn generation_candidates(
        &self,
        group: OutputGroup,
        stage: &str,
        subject: &str,
        request: &str,
    ) -> Result<Vec<GenerationHead>> {
        let ids = self.read(|conn| {
            let mut query = conn.prepare(
                "SELECT g.generation_id FROM artifact_generations g
                 LEFT JOIN generation_heads h ON h.project_id=g.project_id
                 AND h.output_group=g.output_group AND h.subject=g.subject
                 AND h.generation_id=g.generation_id
                 WHERE g.project_id=?1 AND g.output_group=?2 AND g.stage=?3
                 AND g.subject=?4 AND g.request_sha256=?5 AND g.state='COMMITTED'
                 AND NOT EXISTS(SELECT 1 FROM generation_retirements r WHERE r.generation_id=g.generation_id)
                 ORDER BY (h.generation_id IS NOT NULL) DESC, g.committed_at DESC,
                 g.generation_id DESC",
            )?;
            Ok(query
                .query_map(
                    params![self.project_id(), group.as_str(), stage, subject, request],
                    |row| row.get::<_, String>(0),
                )?
                .collect::<std::result::Result<Vec<_>, _>>()?)
        })?;
        ids.iter()
            .map(|id| {
                let head = self
                    .committed_seal(id)?
                    .context("reuse candidate disappeared")?
                    .head;
                if head.stage != stage || head.request_sha256 != request {
                    bail!("indexed request differs from its committed seal");
                }
                Ok(head)
            })
            .collect()
    }

    /// Hash payloads before the short selection transaction. Pruned registered
    /// FASTQs are restored to the registry only when present. Their sealed pruning
    /// authorization remains mandatory; reports and count outputs cannot be pruned.
    pub fn prepare_generation_reuse(
        &self,
        head: &GenerationHead,
    ) -> Result<Option<VerifiedGenerationReuse>> {
        let pin = super::GenerationStorePin::acquire(self.path(), "prepared generation reuse")?;
        if self.generation_retired(&head.generation_id)? {
            return Ok(None);
        }
        let expected_head = self.generation_head(head.output_group, &head.subject)?;
        let Some(seal) = self.committed_seal(&head.generation_id)? else {
            return Ok(None);
        };
        if seal.head != *head || !self.verify_seal(&seal)? {
            return Ok(None);
        }
        let records_json: String = self.read(|conn| {
            Ok(conn.query_row("SELECT artifacts_json FROM artifact_generations WHERE project_id=?1 AND generation_id=?2 AND state='COMMITTED'",
                params![self.project_id(), head.generation_id], |row| row.get(0))?)
        })?;
        let mut records: Vec<ArtifactRecord> = serde_json::from_str(&records_json)?;
        let mut paths = std::collections::BTreeSet::new();
        if records.is_empty()
            || records.iter().any(|record| {
                record.project_id != self.project_id()
                    || record.run_accession.as_deref()
                        != match head.output_group {
                            OutputGroup::Download
                            | OutputGroup::Preprocess
                            | OutputGroup::Quant => Some(head.subject.as_str()),
                            OutputGroup::Refs | OutputGroup::De | OutputGroup::Corr => None,
                        }
                    || !paths.insert(&record.path)
                    || !seal.files.iter().any(|file| {
                        file.path == record.path
                            && file.checksum_type == record.checksum_type
                            && file.checksum == record.checksum
                            && file.bytes == record.bytes
                    })
            })
        {
            return Ok(None);
        }
        records.retain(|record| std::path::Path::new(&record.path).is_file());
        if records.is_empty() {
            return Ok(None);
        }
        Ok(Some(VerifiedGenerationReuse {
            _pin: pin,
            project_id: self.project_id().into(),
            seal,
            records,
            records_json,
            expected_head,
            input_binding: None,
        }))
    }

    /// Select a verified historical generation without modifying any payload or
    /// producer metadata. Registry, head, callback and new attempt success commit
    /// together. The callback must perform only catalog/state/event operations.
    pub fn reattach_generation<T>(
        &self,
        reuse: VerifiedGenerationReuse,
        complete: impl FnOnce(&Database) -> Result<T>,
    ) -> Result<T> {
        if reuse.project_id != self.project_id() {
            bail!("verified generation belongs to another project");
        }
        let identity = self
            .job
            .as_ref()
            .context("generation reuse requires an owned job")?;
        // Other validation may have happened since preparation. Recheck outside
        // the writer transaction, as for fresh generation publication.
        if !self.verify_seal(&reuse.seal)? {
            bail!("historical payload changed after verification");
        }
        if let Some(binding) = &reuse.input_binding {
            self.verify_input_binding(binding)?;
        }
        let head = &reuse.seal.head;
        self.atomic(|db| {
            db.read(|conn| {
                let (stage, run, resource, request): (String, Option<String>, String, Option<String>) = conn.query_row(
                    "SELECT stage,run_accession,resource,request_sha256 FROM job_attempts WHERE attempt_id=?1",
                    [&identity.attempt_id], |row| Ok((row.get(0)?,row.get(1)?,row.get(2)?,row.get(3)?)))?;
                let scope_matches = match head.output_group {
                    OutputGroup::Download | OutputGroup::Preprocess | OutputGroup::Quant => run.as_deref() == Some(&head.subject),
                    OutputGroup::Refs => run.is_none() && resource.strip_prefix("reference:") == Some(&head.subject),
                    OutputGroup::De | OutputGroup::Corr => run.is_none() && head.subject == db.project_id(),
                };
                if stage != head.stage || !scope_matches || request.as_deref() != Some(&head.request_sha256) {
                    bail!("historical generation does not match the owning attempt's scope and request");
                }
                let stored: Option<(String,String)> = conn.query_row(
                    "SELECT inventory_json,artifacts_json FROM artifact_generations WHERE project_id=?1 AND generation_id=?2
                     AND output_group=?3 AND stage=?4 AND subject=?5 AND request_sha256=?6 AND state='COMMITTED'",
                    params![db.project_id(),head.generation_id,head.output_group.as_str(),head.stage,head.subject,head.request_sha256],
                    |row| Ok((row.get(0)?,row.get(1)?))).optional()?;
                let Some((seal_json, records_json)) = stored else { bail!("reuse candidate is no longer committed"); };
                if serde_json::from_str::<Seal>(&seal_json)? != reuse.seal || records_json != reuse.records_json {
                    bail!("reuse candidate changed after verification");
                }
                Ok(())
            })?;
            if db.generation_head(head.output_group, &head.subject)? != reuse.expected_head {
                bail!("selected generation changed after reuse preparation");
            }
            if !db.lineage_matches_catalog(&reuse.seal)? {
                bail!("historical generation lineage changed after verification");
            }
            if let Some(binding) = &reuse.input_binding {
                db.verify_input_binding_catalog(binding)?;
            }
            db.retire_generation_artifacts(head.output_group, &head.subject)?;
            for record in &reuse.records { db.record_artifact(record)?; }
            db.write(|conn| {
                conn.execute("INSERT INTO generation_heads(project_id,output_group,subject,generation_id) VALUES (?1,?2,?3,?4)
                    ON CONFLICT(project_id,output_group,subject) DO UPDATE SET generation_id=excluded.generation_id",
                    params![db.project_id(),head.output_group.as_str(),head.subject,head.generation_id])?;
                Ok(())
            })?;
            db.select_input_binding(head, reuse.input_binding.as_ref())?;
            let result = complete(db)?;
            db.append_event(&head.stage, match head.output_group {
                OutputGroup::Download | OutputGroup::Preprocess | OutputGroup::Quant => Some(&head.subject),
                _ => None,
            }, "selected verified committed generation", serde_json::json!({
                "generation_id": head.generation_id,
                "producer_attempt_id": head.attempt_id,
                "reuse_attempt_id": identity.attempt_id,
                "request_sha256": head.request_sha256,
                "previous_generation_id": reuse.expected_head.as_ref().map(|h| &h.generation_id),
            }))?;
            db.complete_owned_job(JobOutcome::Succeeded, None)?;
            Ok(result)
        })
    }
}
