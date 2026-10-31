//! Quantification execution, immutable publication and cache reconciliation.
use crate::ingest::finish_boolean_job;
use crate::ingest::inputs::{
    load_quant_input_fastqs, load_raw_fastq_artifacts, load_trimmed_fastq_artifacts,
};
use crate::ingest::retention::run_fastq_retention;
use crate::run_state::is_quant_done_state;
use crate::storage::prepare_shared_blob;
use crate::{quality, reference};
use anyhow::Context;
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::generations::{Generation, OutputGroup};
use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::model::{ArtifactKind, ArtifactRecord, QuantArtifacts, RunRecord, VerifiedFile};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::references::{ReferenceScope, ReferenceSelection};
use rnaa_core::shared_store::BlobPublication;
use rnaa_core::state::{ProjectStageState, RunState};
use rnaa_core::util::{compute_sha256, file_size, now_rfc3339};
use rnaa_stages::quant::{QuantIndexGroup, attach_quant_request, snapshot_quant_inputs};
use std::fs;
use std::path::PathBuf;

pub(crate) fn process_quant_run(
    mut lease: JobLease,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::ReferenceBundle,
    index_group: &QuantIndexGroup,
    run: RunRecord,
) -> Result<bool> {
    let run_accession = run.run_accession.clone();
    let result = process_quant_run_owned(
        &lease.database(),
        paths,
        config,
        reference,
        index_group,
        run,
        lease.attempt_id(),
    );
    let success = finish_boolean_job(&mut lease, result)?;
    if success {
        let db = lease.finish_and_release(JobOutcome::Succeeded, None)?;
        run_fastq_retention(&db, paths, config, &run_accession)
            .context("quantification committed, but FASTQ retention failed")?;
    }
    Ok(success)
}

fn process_quant_run_owned(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::model::ReferenceBundle,
    index_group: &QuantIndexGroup,
    run: RunRecord,
    attempt_id: &str,
) -> Result<bool> {
    let eligibility = quality::refresh(db, config, &run.run_accession)?;
    if !eligibility.current.descriptor.eligible() {
        bail!(
            "run {} is ineligible: {}",
            run.run_accession,
            eligibility
                .current
                .descriptor
                .reason
                .as_deref()
                .unwrap_or("QC evidence unavailable")
        );
    }
    let reference_guard = reference::pin(db, reference, config, ReferenceScope::Index)?;
    rnaa_stages::quant::verify_quant_current(db, &run, config)?;
    if !rnaa_stages::quant::source_binding_current(db, &run, config)? {
        bail!(
            "download source selection is unverified; recover the selected FASTQ source before quantification"
        );
    }
    if let Some(candidate) = super::reuse::count_candidate(db, paths, config, reference, &run)? {
        super::reuse::select_counts(db, config, reference, &run, candidate)?;
        return Ok(true);
    }
    let current_paths = paths.for_quant_run(db, &run.run_accession)?;
    let paths = &current_paths;
    let generation_valid =
        quant_generations_valid(db, &run.run_accession, config.quant.preprocess)?;

    if let Some(artifacts) = if generation_valid {
        reconcile_current_quant_artifacts(db, &run, reference, paths, config)?
    } else {
        None
    } {
        let manifest: rnaa_core::StageManifest =
            serde_json::from_slice(&fs::read(&artifacts.manifest_path)?)?;
        let digest = manifest.parameters["request_sha256"]
            .as_str()
            .ok_or_else(|| {
                anyhow::anyhow!("reconciled quantification lacks its request identity")
            })?;
        db.set_job_request(digest)?;
        let prepared = prepare_quant_success(db, config, &run.run_accession, artifacts)?;
        reference_guard.verify_current(db, config, ReferenceScope::Index)?;
        rnaa_stages::quant::verify_quant_current(db, &run, config)?;
        db.atomic(|db| {
            rnaa_stages::quant::verify_quant_catalog(db, &run, config)?;
            reference_guard.verify_catalog(db, config, ReferenceScope::Index)?;
            register_quant_success(
                db,
                &run.run_accession,
                &prepared,
                config.quant.preprocess,
                "quantification reused from existing outputs",
            )?;
            db.complete_owned_job(JobOutcome::Succeeded, None)
        })?;
        finish_shared_publication(&prepared)?;
        return Ok(true);
    }

    db.set_run_state(&run.run_accession, RunState::QuantRunning, None)?;
    db.append_event(
        "quant",
        Some(&run.run_accession),
        "quantification started",
        serde_json::json!({
            "engine": config.quant.engine,
            "threads": config.quant.threads,
            "workers": config.quant.workers,
            "preprocess": config.quant.preprocess,
            "preprocess_strict": config.quant.preprocess_strict
        }),
    )?;
    let quant_fastqs = if config.quant.preprocess {
        load_quant_input_fastqs(db, &run.run_accession, true)?
    } else {
        load_raw_fastq_artifacts(db, &run.run_accession)?
    };
    if quant_fastqs.is_empty() {
        let message = format!("no FASTQ artifacts found for {}", run.run_accession);
        db.set_run_state(&run.run_accession, RunState::QuantFailed, Some(&message))?;
        db.append_event(
            "quant",
            Some(&run.run_accession),
            "quantification failed",
            serde_json::json!({ "error": message }),
        )?;
        return Ok(false);
    }
    if config.quant.preprocess
        && db
            .generation_head(OutputGroup::Preprocess, &run.run_accession)?
            .is_some()
        && !db.verify_generation(OutputGroup::Preprocess, &run.run_accession)?
    {
        bail!(
            "preprocessing generation failed verification for {}",
            run.run_accession
        );
    }

    if db
        .generation_head(OutputGroup::Download, &run.run_accession)?
        .is_some()
        && !db.verify_generation(OutputGroup::Download, &run.run_accession)?
    {
        bail!("download generation failed verification; rerun download before quantification");
    }
    let mut generation = Generation::begin(db, OutputGroup::Quant, "quant")?;
    let output_paths = generation.paths(paths)?;
    let result: Result<QuantArtifacts> = (|| {
        let inputs =
            snapshot_quant_inputs(&quant_registry_inputs(db, &run.run_accession, config)?)?;
        let digest = rnaa_stages::quant::quant_request_digest(&run, reference, config, &inputs)?;
        db.set_job_request(&digest)?;
        db.check_job()?;
        let admission = rnaa_core::resources::ExecutionResources::open(
            &config.resources,
            rnaa_core::resources::ReservationOwner {
                project_root: fs::canonicalize(&paths.root)?,
                attempt_id: attempt_id.into(),
                stage: "quant".into(),
                run_accession: Some(run.run_accession.clone()),
            },
        )?;
        let artifacts = index_group.quantify_with_resources(
            &run,
            &quant_fastqs,
            reference,
            &output_paths,
            config,
            Some(&admission),
            &|| db.check_job(),
        )?;
        db.check_job()?;
        attach_quant_request(&artifacts, &run, reference, config, inputs)?;
        let mut manifest: rnaa_core::StageManifest =
            serde_json::from_slice(&fs::read(&artifacts.manifest_path)?)?;
        manifest.parameters["reference_binding"] = serde_json::to_value(&reference_guard)?;
        rnaa_core::util::write_json_pretty(&artifacts.manifest_path, &manifest)?;
        Ok(artifacts)
    })();
    match result {
        Ok(artifacts) => {
            let prepared = prepare_quant_success(db, config, &run.run_accession, artifacts)?;
            let required = prepared
                .records
                .iter()
                .map(|record| PathBuf::from(&record.path))
                .collect::<Vec<_>>();
            let manifest: rnaa_core::StageManifest =
                serde_json::from_slice(&fs::read(&prepared.artifacts.manifest_path)?)?;
            reference_guard.verify_current(db, config, ReferenceScope::Index)?;
            rnaa_stages::quant::verify_quant_current(db, &run, config)?;
            let mut lineage_inputs = manifest.input_artifacts.clone();
            lineage_inputs.push(rnaa_core::ManifestArtifact::snapshot(
                "REFERENCE_INDEX",
                &reference.kallisto_index_path,
            )?);
            generation.seal_request(
                &required,
                &rnaa_stages::quant::quant_request_value(
                    &run,
                    reference,
                    config,
                    &manifest.input_artifacts,
                )?,
                &lineage_inputs,
            )?;
            generation.commit_prepared(
                &prepared.records,
                |db| register_quant_blob(db, &prepared),
                |db| {
                    rnaa_stages::quant::verify_quant_catalog(db, &run, config)?;
                    reference_guard.verify_catalog(db, config, ReferenceScope::Index)?;
                    complete_quant_success(
                        db,
                        &run.run_accession,
                        &prepared,
                        config.quant.preprocess,
                        "quantification completed",
                    )
                },
            )?;
            finish_shared_publication(&prepared)?;
            Ok(true)
        }
        Err(err) => {
            let message = format!("{err:#}");
            db.set_run_state(&run.run_accession, RunState::QuantFailed, Some(&message))?;
            db.append_event(
                "quant",
                Some(&run.run_accession),
                "quantification failed",
                serde_json::json!({ "error": message }),
            )?;
            Ok(false)
        }
    }
}

pub(crate) fn quant_registry_inputs(
    db: &Database,
    run: &str,
    config: &ProjectConfig,
) -> Result<Vec<VerifiedFile>> {
    let mut inputs = load_raw_fastq_artifacts(db, run)?;
    if config.quant.preprocess {
        inputs.extend(load_trimmed_fastq_artifacts(db, run)?);
    }
    Ok(inputs)
}

pub(crate) fn reconcile_current_quant_artifacts(
    db: &Database,
    run: &RunRecord,
    reference: &rnaa_core::ReferenceBundle,
    paths: &ProjectPaths,
    config: &ProjectConfig,
) -> Result<Option<QuantArtifacts>> {
    let binding = db.generation_input_binding(OutputGroup::Quant, &run.run_accession)?;
    let manifest_path = paths
        .quant_run_dir(&run.run_accession, &config.quant.engine)
        .join("rnaa_quant_manifest.json");
    let Some(manifest) = fs::read(&manifest_path)
        .ok()
        .and_then(|bytes| serde_json::from_slice::<rnaa_core::StageManifest>(&bytes).ok())
    else {
        return Ok(None);
    };
    let inputs = binding
        .as_ref()
        .map(|binding| binding.inputs.as_slice())
        .unwrap_or(&manifest.input_artifacts);
    if !super::reuse::current_count_sources_match(db, config, run, inputs)? {
        return Ok(None);
    }
    rnaa_stages::quant::reconcile_quant_artifacts_with_inputs(
        run,
        reference,
        paths,
        config,
        &quant_registry_inputs(db, &run.run_accession, config)?,
        &db.quant_retention_evidence(&run.run_accession)?,
        binding.as_ref().map(|binding| binding.inputs.as_slice()),
    )
}

pub(crate) fn quant_generations_valid(db: &Database, run: &str, preprocess: bool) -> Result<bool> {
    for group in [
        OutputGroup::Quant,
        OutputGroup::Preprocess,
        OutputGroup::Download,
    ] {
        if group == OutputGroup::Preprocess && !preprocess {
            continue;
        }
        if db.generation_head(group, run)?.is_some() && !db.verify_generation(group, run)? {
            return Ok(false);
        }
    }
    Ok(true)
}

pub(crate) fn reconcile_quant_state(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::ReferenceBundle,
) -> Result<Vec<crate::ingest::recovery::Plan>> {
    Ok(db
        .with_idle_analysis(|db| reconcile_quant_state_idle(db, paths, config, reference))?
        .unwrap_or_default())
}

fn reconcile_quant_state_idle(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::model::ReferenceBundle,
) -> Result<Vec<crate::ingest::recovery::Plan>> {
    let runs = db.list_active_runs()?;
    let mut recovery = Vec::new();
    for queued in runs {
        super::reuse::preprocessing(db, paths, config, &queued)?;
        let (valid, plan) = db
            .with_idle_run(&queued.run_accession, |db| {
                reconcile_quant_run_idle(db, paths, config, reference, &queued.run_accession)
            })?
            .unwrap_or((false, None));
        if plan.is_some() {
            let run = db
                .get_run(&queued.run_accession)?
                .context("recovery run disappeared")?;
            if super::reuse::counts(db, paths, config, reference, &run)? {
                run_fastq_retention(db, paths, config, &queued.run_accession)?;
                continue;
            }
            // History verification releases the run lock. Another producer may
            // have restored reads or completed counts in that interval. Replan
            // under exclusion instead of publishing the earlier recovery action.
            if let Some((fresh_valid, fresh_plan)) =
                db.with_idle_run(&queued.run_accession, |db| {
                    let (valid, plan) = reconcile_quant_run_idle(
                        db,
                        paths,
                        config,
                        reference,
                        &queued.run_accession,
                    )?;
                    if let Some(plan) = &plan {
                        let run = db
                            .get_run(&queued.run_accession)?
                            .context("recovery run disappeared")?;
                        plan.publish(db, &run)?;
                    }
                    Ok((valid, plan))
                })?
            {
                if let Some(plan) = fresh_plan {
                    recovery.push(plan);
                }
                if fresh_valid {
                    run_fastq_retention(db, paths, config, &queued.run_accession)?;
                }
            }
        }
        if valid {
            run_fastq_retention(db, paths, config, &queued.run_accession)?;
        }
    }
    Ok(recovery)
}

fn reconcile_quant_run_idle(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::ReferenceBundle,
    run_accession: &str,
) -> Result<(bool, Option<crate::ingest::recovery::Plan>)> {
    let Some(run) = db.get_run(run_accession)? else {
        return Ok((false, None));
    };
    // Another command may have replaced the report or source since the outer
    // refresh. Recovery decisions must inspect admission under the same run
    // exclusion that protects planning, rather than reuse that older decision.
    quality::refresh(db, config, run_accession)?;
    let reference_guard = reference::pin(db, reference, config, ReferenceScope::Index)?;
    let paths = paths.for_quant_run(db, run_accession)?;
    let complete = is_quant_done_state(run.state);
    // QC/report admission is independent of the validity of already-produced counts.
    let generation_valid = quant_generations_valid(db, run_accession, false)?;
    let reconciled =
        if generation_valid && rnaa_stages::quant::source_binding_current(db, &run, config)? {
            reconcile_current_quant_artifacts(db, &run, reference, &paths, config)?
        } else {
            None
        };
    if let Some(artifacts) = reconciled {
        if !quality::allowed(db, config, run_accession)? {
            if db
                .eligibility(run_accession)?
                .is_some_and(|value| value.descriptor.needs_preprocessing())
            {
                let plan = crate::ingest::recovery::Plan::capture(
                    db,
                    &paths,
                    config,
                    reference,
                    &run,
                    crate::ingest::recovery::Purpose::QcReport,
                )?;
                return Ok((false, Some(plan)));
            }
            return Ok((false, None));
        }
        if !complete {
            finalize_quant_success(
                db,
                config,
                &run,
                artifacts,
                &reference_guard,
                false,
                "quantification reconciled from existing outputs",
            )?;
        }
        return Ok((true, None));
    }
    let has_previous = complete
        || db
            .generation_head(OutputGroup::Quant, run_accession)?
            .is_some();
    if has_previous {
        let plan = crate::ingest::recovery::Plan::capture(
            db,
            &paths,
            config,
            reference,
            &run,
            crate::ingest::recovery::Purpose::Counts,
        )?;
        return Ok((false, Some(plan)));
    }
    Ok((false, None))
}

struct PreparedQuant {
    artifacts: QuantArtifacts,
    records: Vec<ArtifactRecord>,
    shared_blob: Option<BlobPublication>,
}

fn prepare_quant_success(
    db: &Database,
    config: &ProjectConfig,
    run_accession: &str,
    artifacts: QuantArtifacts,
) -> Result<PreparedQuant> {
    let abundance_checksum = compute_sha256(&artifacts.abundance_h5)?;
    let abundance_bytes = file_size(&artifacts.abundance_h5)?;
    let shared_blob = prepare_shared_blob(
        db,
        config,
        &artifacts.abundance_h5,
        &abundance_checksum,
        abundance_bytes,
    )?;
    let mut records = Vec::new();
    for (kind, path) in [
        (
            ArtifactKind::QuantAbundanceH5,
            artifacts.abundance_h5.clone(),
        ),
        (
            ArtifactKind::QuantAbundanceTsv,
            artifacts.out_dir.join("abundance.tsv"),
        ),
        (ArtifactKind::QuantRunInfo, artifacts.run_info_json.clone()),
        (ArtifactKind::Manifest, artifacts.manifest_path.clone()),
    ] {
        let snapshot = rnaa_core::ManifestArtifact::snapshot(kind.as_str(), &path)?;
        let shared = shared_blob
            .as_ref()
            .filter(|_| kind == ArtifactKind::QuantAbundanceH5)
            .map(BlobPublication::record);
        records.push(ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some(run_accession.into()),
            kind,
            path: snapshot.path,
            blob_id: shared.map(|blob| blob.blob_id.clone()),
            shared_path: shared.map(|blob| blob.storage_path.clone()),
            checksum_type: snapshot.checksum_type,
            checksum: snapshot.checksum,
            bytes: snapshot.bytes,
            created_at: now_rfc3339(),
        });
    }
    Ok(PreparedQuant {
        artifacts,
        records,
        shared_blob,
    })
}

fn register_quant_blob(db: &Database, prepared: &PreparedQuant) -> Result<()> {
    if let Some(blob) = &prepared.shared_blob {
        db.record_shared_blob(blob.record())?;
    }
    Ok(())
}

fn finish_shared_publication(prepared: &PreparedQuant) -> Result<()> {
    if let Some(blob) = &prepared.shared_blob {
        blob.finish()?;
    }
    Ok(())
}

fn complete_quant_success(
    db: &Database,
    run_accession: &str,
    prepared: &PreparedQuant,
    preprocessed: bool,
    event_message: &str,
) -> Result<()> {
    db.set_run_state(run_accession, RunState::QuantDone, None)?;
    for stage in ["normalize", "deseq2", "corr"] {
        db.set_project_stage_state(
            stage,
            ProjectStageState::Pending,
            Some("quantified cohort changed"),
        )?;
    }
    db.append_event(
        "quant",
        Some(run_accession),
        event_message,
        serde_json::json!({
            "out_dir": prepared.artifacts.out_dir.display().to_string(),
            "abundance_h5": prepared.artifacts.abundance_h5.display().to_string(),
            "run_info_json": prepared.artifacts.run_info_json.display().to_string(),
            "shared_blob_id": prepared.shared_blob.as_ref().map(|blob| &blob.record().blob_id),
            "shared_blob_path": prepared.shared_blob.as_ref().map(|blob| &blob.record().storage_path),
            "preprocessed": preprocessed,
        }),
    )
}

fn register_quant_success(
    db: &Database,
    run_accession: &str,
    prepared: &PreparedQuant,
    preprocessed: bool,
    event_message: &str,
) -> Result<()> {
    register_quant_blob(db, prepared)?;
    for record in &prepared.records {
        db.record_artifact(record)?;
    }
    complete_quant_success(db, run_accession, prepared, preprocessed, event_message)
}

fn finalize_quant_success(
    db: &Database,
    config: &ProjectConfig,
    run: &RunRecord,
    artifacts: QuantArtifacts,
    reference_guard: &ReferenceSelection,
    preprocessed: bool,
    event_message: &str,
) -> Result<()> {
    let run_accession = &run.run_accession;
    let prepared = prepare_quant_success(db, config, run_accession, artifacts)?;
    reference_guard.verify_current(db, config, ReferenceScope::Index)?;
    rnaa_stages::quant::verify_quant_current(db, run, config)?;
    db.atomic(|db| {
        rnaa_stages::quant::verify_quant_catalog(db, run, config)?;
        reference_guard.verify_catalog(db, config, ReferenceScope::Index)?;
        register_quant_success(db, run_accession, &prepared, preprocessed, event_message)
    })?;
    finish_shared_publication(&prepared)
}
