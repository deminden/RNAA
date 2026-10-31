//! Historical ingest selection with current byte identity and owned publication.
use super::source::{Source, identities, ordered, raw, trimmed};
use crate::{quality, reference};
use anyhow::{Result, bail};
use rnaa_core::{
    Database, ManifestArtifact, ProjectConfig, ProjectPaths, ReferenceBundle, RunRecord, RunState,
    StageManifest,
    generations::{GenerationHead, OutputGroup, VerifiedGenerationReuse},
    jobs::{JobLease, JobOutcome},
    references::ReferenceScope,
};
use std::{fs, path::Path};

/// Completed counts must describe the selected source/processing values even
/// when every physical input was removed under an older retention receipt.
/// Report integrity remains a separate QC admission decision.
pub(super) fn current_count_sources_match(
    db: &Database,
    config: &ProjectConfig,
    run: &RunRecord,
    inputs: &[ManifestArtifact],
) -> Result<bool> {
    let Some(source) = Source::capture(db, run, config)? else {
        return Ok(false);
    };
    let originals = inputs
        .iter()
        .filter(|input| raw(input))
        .cloned()
        .collect::<Vec<_>>();
    if !originals.is_empty() && identities(&originals) != identities(&source.raw) {
        return Ok(false);
    }
    let consumed_trimmed = inputs
        .iter()
        .filter(|input| trimmed(input))
        .cloned()
        .collect::<Vec<_>>();
    if !config.quant.preprocess {
        return Ok(consumed_trimmed.is_empty() && identities(&originals) == identities(&source.raw));
    }
    if let Some(head) = db.generation_head(OutputGroup::Preprocess, &run.run_accession)? {
        if identities(
            &db.generation_inputs(OutputGroup::Preprocess, &run.run_accession)?
                .unwrap_or_default(),
        ) != identities(&source.raw)
        {
            return Ok(false);
        }
        let current = db
            .generation_output_artifacts(&head)?
            .into_iter()
            .filter(trimmed)
            .collect::<Vec<_>>();
        if identities(&current) != identities(&consumed_trimmed) {
            return Ok(false);
        }
    }
    Ok(true)
}

fn index_history(
    db: &Database,
    paths: &ProjectPaths,
    run: &RunRecord,
    config: &ProjectConfig,
    group: OutputGroup,
) -> Result<()> {
    for head in db.committed_generation_history(group, &run.run_accession)? {
        if db.scientific_request(&head.request_sha256)?.is_some() {
            continue;
        }
        let pinned = head.paths(paths)?;
        let path = match group {
            OutputGroup::Quant => pinned
                .quant_run_dir(&run.run_accession, &config.quant.engine)
                .join("rnaa_quant_manifest.json"),
            OutputGroup::Preprocess => pinned
                .preprocess_run_dir(&run.run_accession, &config.quant.engine)
                .join("rnaa_preprocess_manifest.json"),
            _ => unreachable!(),
        };
        let Some(manifest) = fs::read(path)
            .ok()
            .and_then(|bytes| serde_json::from_slice::<StageManifest>(&bytes).ok())
        else {
            continue;
        };
        let request = if group == OutputGroup::Quant {
            manifest.parameters["request"].clone()
        } else {
            serde_json::to_value(rnaa_stages::preprocess::request_from_inputs(
                run,
                &manifest.input_artifacts,
                config,
            )?)?
        };
        let Ok(parsed) = serde_json::from_value::<
            rnaa_core::fingerprint::StageRequest<serde_json::Value>,
        >(request.clone()) else {
            continue;
        };
        if parsed.digest()? == head.request_sha256 && parsed.stage == head.stage {
            db.index_generation_request(&head, &request)?;
        }
    }
    Ok(())
}

/// Rebind a matching preprocessing result to identical newer original reads.
pub(super) fn preprocessing(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    run: &RunRecord,
) -> Result<bool> {
    preprocessing_impl(db, paths, config, run, false)
}

pub(super) fn preprocessing_owned(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    run: &RunRecord,
) -> Result<bool> {
    preprocessing_impl(db, paths, config, run, true)
}

fn preprocessing_impl(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    run: &RunRecord,
    already_owned: bool,
) -> Result<bool> {
    if !config.quant.preprocess {
        return Ok(false);
    }
    let Some(source) = Source::capture(db, run, config)? else {
        return Ok(false);
    };
    if ordered(
        db.generation_inputs(OutputGroup::Preprocess, &run.run_accession)?
            .unwrap_or_default(),
    ) == source.raw
        && quality::allowed(db, config, &run.run_accession)?
    {
        return Ok(false);
    }
    index_history(db, paths, run, config, OutputGroup::Preprocess)?;
    let request = rnaa_stages::preprocess::request_from_inputs(run, &source.raw, config)?;
    let digest = request.digest()?;
    for head in db.generation_candidates(
        OutputGroup::Preprocess,
        "preprocess",
        &run.run_accession,
        &digest,
    )? {
        // An execution worker needs actual trimmed reads for a new count fit.
        // Evidence-only selection is useful before historical count lookup, but
        // cannot replace reconstruction when the requested count result is absent.
        if already_owned
            && db
                .generation_output_artifacts(&head)?
                .iter()
                .filter(|input| trimmed(input))
                .any(|input| !Path::new(&input.path).is_file())
        {
            continue;
        }
        let Some(verified) = db.prepare_generation_reuse(&head)? else {
            continue;
        };
        let pinned = head.paths(paths)?;
        let report_path = pinned
            .preprocess_run_dir(&run.run_accession, &config.quant.engine)
            .join("fasterp.json");
        let qc = super::preprocess::preprocess_qc_for_policy(
            db.project_id(),
            &run.run_accession,
            &report_path,
            config,
        )?;
        let report = ManifestArtifact::snapshot("PREPROCESS_REPORT", &report_path)?;
        let assessment = rnaa_core::eligibility::Assessment::from_qc(
            config,
            rnaa_core::eligibility::Evidence {
                report,
                producer_generation: Some(head.generation_id.clone()),
            },
            qc.clone(),
        );
        if !assessment.descriptor.eligible() {
            continue;
        }
        let mut lease = if already_owned {
            None
        } else {
            let Some(lease) =
                JobLease::try_acquire(db, "preprocess", Some(&run.run_accession), &[run.state])?
            else {
                return Ok(false);
            };
            Some(lease)
        };
        let owned = lease
            .as_ref()
            .map(JobLease::database)
            .unwrap_or_else(|| db.clone());
        let result = (|| {
            owned.set_job_request(&digest)?;
            let verified = owned.bind_generation_reuse(verified, &source.raw)?;
            source.verify_current(&owned, run, config)?;
            super::preprocess::verify_settings(&owned, config)?;
            owned.reattach_generation(verified, |db| {
                source.verify_catalog(db, run, config)?;
                super::preprocess::verify_catalog(db, config)?;
                let changed = db.record_eligibility(&run.run_accession, &assessment)?;
                if changed.changed {
                    db.record_preprocess_qc(&qc)?;
                }
                db.set_run_state(&run.run_accession, RunState::PreprocessDone, None)
            })
        })();
        if let Some(lease) = &mut lease {
            match &result {
                Ok(()) => lease.finish(JobOutcome::Succeeded, None)?,
                Err(error) => lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))?,
            }
        }
        result?;
        return Ok(true);
    }
    Ok(false)
}

pub(super) struct CountCandidate {
    verified: VerifiedGenerationReuse,
    source: Source,
    inputs: Vec<ManifestArtifact>,
    preprocessing: Option<GenerationHead>,
    preprocessing_binding: Option<rnaa_core::generations::GenerationInputBinding>,
    eligibility: Option<String>,
}

pub(super) fn count_candidate(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &ReferenceBundle,
    run: &RunRecord,
) -> Result<Option<CountCandidate>> {
    if !quality::allowed(db, config, &run.run_accession)? {
        return Ok(None);
    }
    let Some(source) = Source::capture(db, run, config)? else {
        return Ok(None);
    };
    let preprocessing = db.generation_head(OutputGroup::Preprocess, &run.run_accession)?;
    let preprocessing_binding =
        db.generation_input_binding(OutputGroup::Preprocess, &run.run_accession)?;
    let mut inputs = source.raw.clone();
    if config.quant.preprocess {
        let Some(head) = &preprocessing else {
            return Ok(None);
        };
        if !db.verify_generation_head(head)? {
            return Ok(None);
        }
        let pinned = head.paths(paths)?;
        let manifest: StageManifest = serde_json::from_slice(&fs::read(
            pinned
                .preprocess_run_dir(&run.run_accession, &config.quant.engine)
                .join("rnaa_preprocess_manifest.json"),
        )?)?;
        let bypass = manifest.parameters["bypass"] == true && !config.quant.preprocess_strict;
        if (!rnaa_stages::preprocess::processing_identity_matches(run, &manifest) && !bypass)
            || identities(
                &db.generation_inputs(OutputGroup::Preprocess, &run.run_accession)?
                    .unwrap_or_default(),
            ) != identities(&source.raw)
        {
            return Ok(None);
        }
        let trimmed = db
            .generation_output_artifacts(head)?
            .into_iter()
            .filter(trimmed)
            .collect::<Vec<_>>();
        if trimmed.is_empty() && !bypass {
            return Ok(None);
        }
        inputs.extend(trimmed);
    }
    index_history(db, paths, run, config, OutputGroup::Quant)?;
    let policy = rnaa_stages::quant::quant_policy_fingerprint(run, reference, config)?;
    for head in
        db.generation_policy_candidates(OutputGroup::Quant, "quant", &run.run_accession, &policy)?
    {
        let Some(verified) = db.prepare_generation_reuse(&head)? else {
            continue;
        };
        if rnaa_stages::quant::historical_quant_artifacts(
            run,
            reference,
            &head.paths(paths)?,
            config,
            &inputs,
        )?
        .is_none()
        {
            continue;
        }
        return Ok(Some(CountCandidate {
            verified,
            source,
            inputs,
            preprocessing,
            preprocessing_binding,
            eligibility: db
                .eligibility(&run.run_accession)?
                .map(|value| value.fingerprint),
        }));
    }
    Ok(None)
}

pub(super) fn select_counts(
    db: &Database,
    config: &ProjectConfig,
    reference: &ReferenceBundle,
    run: &RunRecord,
    candidate: CountCandidate,
) -> Result<()> {
    let guard = reference::pin(db, reference, config, ReferenceScope::Index)?;
    db.set_job_request(&candidate.verified.head().request_sha256)?;
    let verified = db.bind_generation_reuse(candidate.verified, &candidate.inputs)?;
    candidate.source.verify_current(db, run, config)?;
    guard.verify_current(db, config, ReferenceScope::Index)?;
    if let Some(head) = &candidate.preprocessing
        && config.quant.preprocess
        && !db.verify_generation_head(head)?
    {
        bail!("preprocessing changed before count reuse");
    }
    db.reattach_generation(verified, |db| {
        candidate.source.verify_catalog(db, run, config)?;
        guard.verify_catalog(db, config, ReferenceScope::Index)?;
        if db.generation_head(OutputGroup::Preprocess, &run.run_accession)?
            != candidate.preprocessing
            || db.generation_input_binding(OutputGroup::Preprocess, &run.run_accession)?
                != candidate.preprocessing_binding
            || db
                .eligibility(&run.run_accession)?
                .map(|value| value.fingerprint)
                != candidate.eligibility
        {
            bail!("preprocessing or QC admission changed before count reuse");
        }
        db.set_run_state(&run.run_accession, RunState::QuantDone, None)?;
        for stage in ["normalize", "deseq2", "corr"] {
            db.set_project_stage_state(
                stage,
                rnaa_core::ProjectStageState::Pending,
                Some("selected count request changed"),
            )?;
        }
        db.append_event(
            "quant",
            Some(&run.run_accession),
            "quantification reused from verified history",
            serde_json::json!({"input_witnesses":candidate.inputs,"source_observations":candidate.source.raw}),
        )
    })
}

pub(super) fn counts(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &ReferenceBundle,
    run: &RunRecord,
) -> Result<bool> {
    let Some(candidate) = count_candidate(db, paths, config, reference, run)? else {
        return Ok(false);
    };
    let Some(mut lease) =
        JobLease::try_acquire(db, "quant", Some(&run.run_accession), &[run.state])?
    else {
        return Ok(false);
    };
    match select_counts(&lease.database(), config, reference, run, candidate) {
        Ok(()) => {
            lease.finish(JobOutcome::Succeeded, None)?;
            Ok(true)
        }
        Err(error) => {
            lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))?;
            Err(error)
        }
    }
}
