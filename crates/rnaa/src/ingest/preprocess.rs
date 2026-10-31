//! Preprocessing execution and atomic publication of QC evidence.
use crate::ingest::finish_boolean_job;
use crate::ingest::inputs::load_raw_fastq_artifacts;
use crate::ingest::source::{Source, ordered};
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::generations::{Generation, OutputGroup};
use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::model::{ArtifactKind, ArtifactRecord, RunRecord};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::RunState;
use rnaa_core::util::now_rfc3339;
use rnaa_stages::preprocess::FasterpInProcessPreprocessor;
use rnaa_stages::quant::snapshot_quant_inputs;
use serde_json::json;
use std::path::{Path, PathBuf};

pub(crate) fn process_preprocess_run(
    mut lease: JobLease,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    run: RunRecord,
) -> Result<bool> {
    let result =
        process_preprocess_run_owned(&lease.database(), paths, config, run, lease.attempt_id());
    finish_boolean_job(&mut lease, result)
}

fn process_preprocess_run_owned(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    run: RunRecord,
    attempt_id: &str,
) -> Result<bool> {
    if super::reuse::preprocessing_owned(db, paths, config, &run)? {
        return Ok(true);
    }
    if db
        .generation_head(OutputGroup::Download, &run.run_accession)?
        .is_some()
        && !db.verify_generation(OutputGroup::Download, &run.run_accession)?
    {
        return fail_preprocess(
            db,
            &run.run_accession,
            "download generation failed verification; rerun download before preprocessing",
        );
    }
    let raw_fastqs = load_raw_fastq_artifacts(db, &run.run_accession)?;
    if raw_fastqs.is_empty() {
        return fail_preprocess(
            db,
            &run.run_accession,
            "no raw FASTQ artifacts available for preprocessing",
        );
    }
    let frozen_inputs = snapshot_quant_inputs(&raw_fastqs)?;
    let source = Source::capture(db, &run, config)?
        .ok_or_else(|| anyhow::anyhow!("preprocessing source is unavailable or unverified"))?;
    if source.raw != ordered(frozen_inputs.clone()) {
        bail!("preprocessing requires the complete current original FASTQ group");
    }
    let request =
        rnaa_stages::preprocess::request_digest_from_inputs(&run, &frozen_inputs, config)?;
    db.set_job_request(&request)?;
    let paths = paths.for_quant_run(db, &run.run_accession)?;
    let generation_valid = db
        .generation_head(OutputGroup::Preprocess, &run.run_accession)?
        .is_none()
        || db.verify_generation(OutputGroup::Preprocess, &run.run_accession)?;
    if generation_valid
        && let Some(artifacts) = rnaa_stages::preprocess::reconcile_preprocess_artifacts(
            &run,
            &raw_fastqs,
            &paths,
            config,
        )?
    {
        for input in &frozen_inputs {
            if !input.matches_file()? {
                bail!("preprocessing input changed during cache validation");
            }
        }
        let prepared = prepare_preprocess_success(
            db,
            config,
            &run.run_accession,
            artifacts,
            db.generation_head(OutputGroup::Preprocess, &run.run_accession)?
                .map(|head| head.generation_id),
        )?;
        let outcome = if prepared.qc.gate_passed {
            JobOutcome::Succeeded
        } else {
            JobOutcome::Failed
        };
        source.verify_current(db, &run, config)?;
        verify_settings(db, config)?;
        db.atomic(|db| {
            source.verify_catalog(db, &run, config)?;
            verify_catalog(db, config)?;
            for record in &prepared.records {
                db.record_artifact(record)?;
            }
            complete_preprocess(db, &run.run_accession, &prepared)?;
            db.complete_owned_job(outcome, prepared.qc.gate_reason.as_deref())
        })?;
        return Ok(prepared.qc.gate_passed);
    }
    db.append_event("preprocess", Some(&run.run_accession), "preprocessing started",
        json!({"strict":config.quant.preprocess_strict, "max_input_mb":config.quant.preprocess_max_input_mb,
            "threads":config.quant.preprocess_threads}))?;
    let mut generation = Generation::begin(db, OutputGroup::Preprocess, "preprocess")?;
    let output_paths = generation.paths(&paths)?;
    db.check_job()?;
    let result = (|| {
        let admission = rnaa_core::resources::ExecutionResources::open(
            &config.resources,
            rnaa_core::resources::ReservationOwner {
                project_root: std::fs::canonicalize(&paths.root)?,
                attempt_id: attempt_id.into(),
                stage: "preprocess".into(),
                run_accession: Some(run.run_accession.clone()),
            },
        )?;
        FasterpInProcessPreprocessor.preprocess_with_resources(
            &run,
            &raw_fastqs,
            &output_paths,
            config,
            Some(&admission),
            &|| db.check_job(),
        )
    })();
    let result = match result {
        Err(error)
            if !config.quant.preprocess_strict
                && error
                    .downcast_ref::<rnaa_core::resources::ResourceAdmissionError>()
                    .is_none() =>
        {
            rnaa_stages::preprocess::bypass(
                &run,
                &output_paths,
                config,
                &frozen_inputs,
                &format!("{error:#}"),
            )
        }
        other => other,
    };
    let artifacts = match result {
        Ok(artifacts) => artifacts,
        Err(error) => {
            drop(generation);
            return fail_preprocess(db, &run.run_accession, &format!("{error:#}"));
        }
    };
    let prepared = prepare_preprocess_success(
        db,
        config,
        &run.run_accession,
        artifacts,
        Some(generation.head().generation_id.clone()),
    )?;
    let required = prepared
        .records
        .iter()
        .filter(|record| record.bytes > 0)
        .map(|record| PathBuf::from(&record.path))
        .collect::<Vec<_>>();
    generation.seal_request(
        &required,
        &rnaa_stages::preprocess::request_from_inputs(&run, &frozen_inputs, config)?,
        &frozen_inputs,
    )?;
    let outcome = if prepared.qc.gate_passed {
        JobOutcome::Succeeded
    } else {
        JobOutcome::Failed
    };
    source.verify_current(db, &run, config)?;
    verify_settings(db, config)?;
    generation.commit_outcome(
        &prepared.records,
        |_| Ok(()),
        |db| {
            source.verify_catalog(db, &run, config)?;
            verify_catalog(db, config)?;
            complete_preprocess(db, &run.run_accession, &prepared)
        },
        outcome,
        prepared.qc.gate_reason.as_deref(),
    )?;
    Ok(prepared.qc.gate_passed)
}

fn fail_preprocess(db: &Database, run: &str, message: &str) -> Result<bool> {
    db.atomic(|db| {
        db.set_run_state(run, RunState::PreprocessFailed, Some(message))?;
        db.append_event(
            "preprocess",
            Some(run),
            "preprocessing failed",
            json!({"error":message}),
        )?;
        db.complete_owned_job(JobOutcome::Failed, Some(message))
    })?;
    Ok(false)
}

/// Source observations fence scientific source/layout/QC settings. These checks
/// additionally bind the size limit for both new execution and historical reuse.
pub(super) fn verify_catalog(db: &Database, config: &ProjectConfig) -> Result<()> {
    let latest: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
    if latest.quant.preprocess_max_input_mb != config.quant.preprocess_max_input_mb {
        bail!("preprocessing size policy changed during the attempt");
    }
    Ok(())
}

pub(super) fn verify_settings(db: &Database, config: &ProjectConfig) -> Result<()> {
    verify_catalog(db, config)?;
    let latest =
        ProjectConfig::load(&ProjectPaths::new(Path::new(&db.project()?.root_dir)).config_path)?;
    if latest.quant.preprocess_max_input_mb != config.quant.preprocess_max_input_mb {
        bail!("preprocessing size policy changed on disk during the attempt");
    }
    Ok(())
}

struct PreparedPreprocess {
    artifacts: rnaa_core::PreprocessArtifacts,
    records: Vec<ArtifactRecord>,
    qc: rnaa_core::PreprocessQcRecord,
    assessment: rnaa_core::eligibility::Assessment,
}

fn prepare_preprocess_success(
    db: &Database,
    config: &ProjectConfig,
    run: &str,
    artifacts: rnaa_core::PreprocessArtifacts,
    producer: Option<String>,
) -> Result<PreparedPreprocess> {
    let mut records = Vec::new();
    let outputs = artifacts
        .fastqs
        .iter()
        .map(|file| (file.artifact_kind, file.path.clone()))
        .chain([
            (
                ArtifactKind::PreprocessReport,
                artifacts.report_json.clone(),
            ),
            (
                ArtifactKind::Manifest,
                artifacts.out_dir.join("rnaa_preprocess_manifest.json"),
            ),
        ]);
    for (kind, path) in outputs {
        let snapshot = rnaa_core::ManifestArtifact::snapshot(kind.as_str(), &path)?;
        records.push(ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some(run.into()),
            kind,
            path: snapshot.path,
            blob_id: None,
            shared_path: None,
            checksum_type: snapshot.checksum_type,
            checksum: snapshot.checksum,
            bytes: snapshot.bytes,
            created_at: now_rfc3339(),
        });
    }
    let qc = preprocess_qc_for_policy(db.project_id(), run, &artifacts.report_json, config)?;
    let report = records
        .iter()
        .find(|record| record.kind == ArtifactKind::PreprocessReport)
        .expect("report is prepared");
    let evidence = rnaa_core::eligibility::Evidence {
        report: rnaa_core::ManifestArtifact {
            kind: report.kind.as_str().into(),
            path: report.path.clone(),
            checksum_type: report.checksum_type.clone(),
            checksum: report.checksum.clone(),
            bytes: report.bytes,
        },
        producer_generation: producer,
    };
    let assessment = rnaa_core::eligibility::Assessment::from_qc(config, evidence, qc.clone());
    Ok(PreparedPreprocess {
        artifacts,
        records,
        qc,
        assessment,
    })
}

fn complete_preprocess(db: &Database, run: &str, prepared: &PreparedPreprocess) -> Result<()> {
    let qc = &prepared.qc;
    let change = db.record_eligibility(run, &prepared.assessment)?;
    if change.changed {
        db.record_preprocess_qc(qc)?;
    }
    db.set_run_state(run, RunState::PreprocessDone, None)?;
    let message = if qc.mode == "bypass" && qc.gate_passed {
        "preprocessing bypassed; using raw FASTQ"
    } else if qc.gate_passed {
        "preprocessing completed"
    } else {
        "preprocessing rejected by QC gate"
    };
    db.append_event("preprocess", Some(run), message, json!({
        "tool":prepared.artifacts.tool_name, "version":prepared.artifacts.tool_version, "reused":prepared.artifacts.reused,
        "passed_reads":qc.passed_reads, "failed_reads":qc.failed_reads, "report_json":prepared.artifacts.report_json,
        "qc_gate_passed":qc.gate_passed, "qc_gate_reason":qc.gate_reason,
        "outputs":prepared.artifacts.fastqs.iter().map(|file| file.path.display().to_string()).collect::<Vec<_>>()
    }))
}

pub(super) fn preprocess_qc_for_policy(
    project: &str,
    run: &str,
    report: &Path,
    config: &ProjectConfig,
) -> Result<rnaa_core::PreprocessQcRecord> {
    let mut qc = parse_preprocess_qc_report(project, run, report, &config.quant.qc_gate)?;
    if qc.mode == "bypass" {
        qc.gate_passed = !config.quant.preprocess_strict;
        qc.gate_reason = config
            .quant
            .preprocess_strict
            .then(|| "preprocessing bypass is disabled by current policy".into());
    }
    Ok(qc)
}

pub(crate) fn parse_preprocess_qc_report(
    project: &str,
    run: &str,
    report: &Path,
    gate: &rnaa_core::PreprocessQcGate,
) -> Result<rnaa_core::PreprocessQcRecord> {
    rnaa_core::eligibility::parse_report(project, run, report, gate)
}
