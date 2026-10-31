//! Plan missing count inputs or QC report reconstruction independently.
use super::inputs::{load_raw_fastq_artifacts, load_trimmed_fastq_artifacts};
use crate::quality;
use anyhow::Result;
use rnaa_core::{
    Database, ProjectConfig, ProjectPaths, ReferenceBundle, RunRecord, RunState,
    eligibility::EligibilityState,
    generations::{GenerationHead, OutputGroup},
};
use rnaa_stages::{download::ReqwestDownloader, fastq::validate_fastqs};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Purpose {
    Counts,
    QcReport,
}
impl Purpose {
    fn as_str(self) -> &'static str {
        match self {
            Self::Counts => "counts",
            Self::QcReport => "qc_report",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Action {
    QuantifyRaw,
    PreprocessRaw,
    QuantifyTrimmed,
    RestoreTrash,
    Download,
    Blocked,
}
impl Action {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::QuantifyRaw => "quantify_raw",
            Self::PreprocessRaw => "preprocess_raw",
            Self::QuantifyTrimmed => "quantify_trimmed",
            Self::RestoreTrash => "restore_trash",
            Self::Download => "download",
            Self::Blocked => "blocked",
        }
    }
    fn next(self) -> RunState {
        match self {
            Self::QuantifyRaw | Self::PreprocessRaw | Self::RestoreTrash => RunState::Verified,
            Self::QuantifyTrimmed => RunState::PreprocessDone,
            Self::Download => RunState::Resolved,
            Self::Blocked => RunState::QuantFailed,
        }
    }
}
pub(crate) struct Plan {
    pub(crate) run: String,
    pub(crate) action: Action,
    target: String,
    source_request: String,
    previous: Option<GenerationHead>,
    reason: String,
    available_inputs: Vec<rnaa_core::ManifestArtifact>,
    fingerprint: String,
    purpose: Purpose,
}
impl Plan {
    pub(crate) fn capture(
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        reference: &ReferenceBundle,
        run: &RunRecord,
        purpose: Purpose,
    ) -> Result<Self> {
        let target = rnaa_stages::quant::quant_policy_fingerprint(run, reference, config)?;
        let raw = load_raw_fastq_artifacts(db, &run.run_accession)?;
        let trimmed = load_trimmed_fastq_artifacts(db, &run.run_accession)?;
        let valid = |group| -> Result<bool> {
            Ok(db.generation_head(group, &run.run_accession)?.is_none()
                || db.verify_generation(group, &run.run_accession)?)
        };
        let qc_rejected = config.quant.preprocess
            && db
                .eligibility(&run.run_accession)?
                .is_some_and(|value| value.descriptor.state == EligibilityState::Rejected);
        let source_current = rnaa_stages::quant::source_binding_current(db, run, config)?;
        let trimmed_ready = source_current
            && config.quant.preprocess
            && quality::allowed(db, config, &run.run_accession)?
            && valid(OutputGroup::Preprocess)?
            && rnaa_stages::preprocess::retained_processing_compatible(run, paths, config)?
            && validate_fastqs(run.library_layout, &trimmed, &|| db.check_job()).is_ok();
        let raw_ready = source_current
            && !trimmed_ready
            && valid(OutputGroup::Download)?
            && validate_fastqs(run.library_layout, &raw, &|| db.check_job()).is_ok();
        let downloader = ReqwestDownloader;
        let (action, reason) = if qc_rejected {
            (
                Action::Blocked,
                "current QC policy rejects the verified preprocessing evidence",
            )
        } else if trimmed_ready {
            (
                Action::QuantifyTrimmed,
                "compatible retained trimmed reads are available",
            )
        } else if raw_ready {
            (
                if config.quant.preprocess {
                    Action::PreprocessRaw
                } else {
                    Action::QuantifyRaw
                },
                "a complete verified original FASTQ group is available",
            )
        } else if source_current && db.restorable_file_trash(&run.run_accession)? {
            (
                Action::RestoreTrash,
                "restore a complete verified read group from project trash",
            )
        } else if downloader.fastq_recovery_available(run, config.download.prefer) {
            (
                Action::Download,
                "required reads are missing, incomplete or invalid; restage the complete selected FASTQ source",
            )
        } else {
            (
                Action::Blocked,
                "required reads are unavailable and no complete supported FASTQ source is selected; restore/register original reads or resolve a downloadable FASTQ source",
            )
        };
        let available_inputs = match action {
            Action::QuantifyTrimmed => rnaa_stages::quant::snapshot_quant_inputs(&trimmed)?,
            Action::PreprocessRaw | Action::QuantifyRaw => {
                rnaa_stages::quant::snapshot_quant_inputs(&raw)?
            }
            _ => Vec::new(),
        };
        let source_request = downloader.source_request_digest(run, config.download.prefer)?;
        let previous = db.generation_head(OutputGroup::Quant, &run.run_accession)?;
        let fingerprint = rnaa_core::fingerprint::StageRequest {
            schema_version:1, stage:"dependency_recovery".into(),algorithm_version:3,engines:Default::default(),
            parameters:serde_json::json!({"purpose":purpose.as_str(),"action":action.as_str(),"target_policy":target,"source_request":source_request,"previous":previous,"input_observations":available_inputs}),
            inputs:available_inputs.iter().map(rnaa_core::fingerprint::InputIdentity::from).collect(),dependencies:Default::default(),
        }.digest()?;
        Ok(Self {
            run: run.run_accession.clone(),
            action,
            target,
            source_request,
            previous,
            available_inputs,
            fingerprint,
            reason: reason.into(),
            purpose,
        })
    }

    pub(crate) fn publish(&self, db: &Database, run: &RunRecord) -> Result<()> {
        let reason = format!(
            "{} recovery {} [{}]: {}",
            self.purpose.as_str(),
            self.action.as_str(),
            &self.fingerprint[..12],
            self.reason
        );
        db.atomic(|db| {
            if db.generation_head(OutputGroup::Quant, &self.run)? != self.previous {
                anyhow::bail!("count generation changed while planning recovery; retry");
            }
            let next = if self.purpose == Purpose::QcReport && self.action == Action::Blocked {
                // Missing QC proof does not make the verified count payload fail.
                run.state
            } else {
                self.action.next()
            };
            if run.state != next || run.last_error.as_deref() != Some(&reason) {
                db.set_run_state(&self.run, next, Some(&reason))?;
                for stage in ["normalize", "deseq2", "corr"] {
                    if db.get_project_stage_state(stage)?.is_some() {
                        db.set_project_stage_state(
                            stage,
                            rnaa_core::ProjectStageState::Pending,
                            Some(&reason),
                        )?;
                    }
                }
                let message = match self.purpose {
                    Purpose::Counts => "planned count dependency recovery",
                    Purpose::QcReport => "planned QC report recovery",
                };
                db.append_event(
                    "recovery",
                    Some(&self.run),
                    message,
                    serde_json::json!({
                        "schema_version":2,
                        "purpose":self.purpose.as_str(),
                        "action":self.action.as_str(),
                        "reason":self.reason,
                        "target_policy_sha256":self.target,
                        "download_source_sha256":self.source_request,
                        "previous_quant_generation":self.previous,
                        "plan_sha256":self.fingerprint,
                        "available_inputs":self.available_inputs,
                    }),
                )?;
            }
            Ok(())
        })
    }
}

/// Only recovery plans restage inputs for `quant`; unrelated newly resolved runs
/// stay under the explicit downloader/orchestrator workflow.
pub(crate) fn restage(
    plans: &[Plan],
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
) -> Result<std::collections::HashSet<String>> {
    let mut failed = std::collections::HashSet::new();
    for plan in plans {
        if plan.action == Action::Blocked {
            failed.insert(plan.run.clone());
        }
        if plan.action == Action::RestoreTrash {
            if let Err(error) = restore_trash(db, paths, config, plan) {
                db.append_event(
                    "recovery",
                    Some(&plan.run),
                    "input trash restoration failed",
                    serde_json::json!({"error":format!("{error:#}")}),
                )?;
                failed.insert(plan.run.clone());
            }
            continue;
        }
        if plan.action != Action::Download {
            continue;
        }
        if let Err(error) = ReqwestDownloader.recover_inputs(db, paths, config, &plan.run) {
            db.append_event(
                "recovery",
                Some(&plan.run),
                "input restaging failed",
                serde_json::json!({"error":format!("{error:#}")}),
            )?;
            failed.insert(plan.run.clone());
        }
    }
    Ok(failed)
}

pub(crate) fn restore_trash(
    db: &Database,
    _paths: &ProjectPaths,
    config: &ProjectConfig,
    plan: &Plan,
) -> Result<()> {
    let run = db
        .get_run(&plan.run)?
        .ok_or_else(|| anyhow::anyhow!("recovery run disappeared"))?;
    if !rnaa_stages::quant::source_binding_current(db, &run, config)? {
        anyhow::bail!("source selection changed before trash recovery; reconcile again");
    }
    if !db.restore_run_file_trash(&plan.run)? {
        anyhow::bail!("no complete current read group remains restorable");
    }
    Ok(())
}
