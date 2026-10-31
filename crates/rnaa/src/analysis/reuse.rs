//! Shared current-request admission and immutable historical result selection.
use crate::quality;
use anyhow::Result;
use rnaa_core::{
    Database, ProjectConfig, ProjectPaths, RunRecord,
    generations::OutputGroup,
    references::{ReferenceScope, ReferenceSelection},
    state::ProjectStageState,
    util::file_size,
};
use rnaa_stages::analysis_request::AnalysisRequest;
use serde::Serialize;
use std::path::PathBuf;

pub(super) struct AnalysisReuse<'a> {
    pub paths: &'a ProjectPaths,
    pub db: &'a Database,
    pub config: &'a ProjectConfig,
    pub runs: &'a [RunRecord],
    pub eligibility: &'a quality::Snapshot,
    pub reference: &'a ReferenceSelection,
    pub group: OutputGroup,
    pub stage: &'a str,
}

pub(super) struct OutputContract {
    pub descriptor: PathBuf,
    pub manifest: PathBuf,
    pub required: Vec<PathBuf>,
}

impl AnalysisReuse<'_> {
    pub fn reconcile<P: Serialize>(
        &self,
        request: &AnalysisRequest<P>,
        outputs: impl Fn(&ProjectPaths) -> OutputContract,
    ) -> Result<Option<(u64, u64)>> {
        let Self {
            db,
            config,
            eligibility,
            reference,
            ..
        } = self;
        reference.verify_current(db, config, ReferenceScope::Analysis)?;
        eligibility.verify_current(db, config)?;
        request.verify_current_settings(db)?;
        let digest = request.request.digest()?;
        db.set_job_request(&digest)?;
        for head in db.generation_candidates(self.group, self.stage, db.project_id(), &digest)? {
            let Some(verified) = db.prepare_generation_reuse(&head)? else {
                continue;
            };
            let contract = outputs(&head.paths(self.paths)?);
            let manifest_stage = if self.group == OutputGroup::De {
                "deseq2"
            } else {
                "corr"
            };
            if !request.matches_descriptor(&contract.descriptor)?
                || !request.matches(
                    &contract.manifest,
                    manifest_stage,
                    db.project_id(),
                    &contract.required,
                )?
                || contract
                    .required
                    .iter()
                    .any(|path| file_size(path).unwrap_or_default() == 0)
            {
                continue;
            }
            reference.verify_current(db, config, ReferenceScope::Analysis)?;
            eligibility.verify_current(db, config)?;
            request.verify_current_settings(db)?;
            request.verify_inputs()?;
            let metadata_witnesses = request.retain_metadata_witnesses(
                &self.paths.metadata_dir.join("analysis_observations"),
            )?;
            let observation = request.configuration_observation()?;
            let recorded = verified.artifact_count() as u64;
            db.reattach_generation(verified, |db| {
                request.verify_settings_catalog(db)?;
                reference.verify_catalog(db, config, ReferenceScope::Analysis)?;
                eligibility.verify_catalog(db, config)?;
                db.set_project_stage_state(self.stage, ProjectStageState::Done, None)?;
                db.append_event(
                    self.stage,
                    None,
                    "stage reconciled from existing outputs",
                    serde_json::json!({
                        "generation_id": head.generation_id,
                        "request_manifest": contract.descriptor,
                        "metadata_witnesses": metadata_witnesses,
                        "configuration_observation": observation,
                        "outputs_recorded": recorded,
                        "contributing_runs": self.runs.len(),
                    }),
                )?;
                Ok(())
            })?;
            return Ok(Some((contract.required.len() as u64 + 1, recorded)));
        }
        Ok(None)
    }
}
