//! Normalization and DE request snapshots, execution, reuse and publication.
use crate::analysis::inputs::{
    bind_metadata, read_samplesheet_header, validate_deseq2_quant_inputs, validate_design_formula,
    write_filtered_deseq2_samplesheet,
};
use crate::analysis::prepare_analysis_artifacts;
use crate::cli::{Deseq2Args, NormalizeArgs, parse_cli_contrasts};
use crate::ingest::inputs::recover_interrupted_work;
use crate::{quality, reference};
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::generations::{Generation, OutputGroup};
use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::model::{ContrastSpec, RunRecord};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::references::{ReferenceScope, ReferenceSelection};
use rnaa_core::state::{ProjectStageState, RunState};
use rnaa_core::traits::{DifferentialExpression, ReferenceManager};
use rnaa_stages::analysis_request::deseq_request;
use rnaa_stages::deseq2::RustDeseq2Runner;
use rnaa_stages::refs::EnsemblReferenceManager;
use std::fs;
use std::path::Path;

pub(crate) fn cmd_deseq2(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: Deseq2Args,
) -> Result<()> {
    recover_interrupted_work(db)?;
    let contrasts = if args.contrast.is_empty() {
        db.list_contrasts()?
    } else {
        parse_cli_contrasts(&args.contrast)?
    };
    let (_outputs, _recorded) =
        run_deseq2_stage(paths, db, config, args.design, contrasts, "deseq2", true)?;
    Ok(())
}

pub(crate) fn cmd_normalize(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: NormalizeArgs,
) -> Result<()> {
    recover_interrupted_work(db)?;
    let (outputs, recorded) = run_deseq2_stage(
        paths,
        db,
        config,
        args.design,
        Vec::new(),
        "normalize",
        false,
    )?;
    println!("normalize_outputs\t{}", outputs);
    println!("normalize_recorded\t{}", recorded);
    Ok(())
}

fn run_deseq2_stage(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    design: Option<String>,
    contrasts: Vec<ContrastSpec>,
    stage: &str,
    store_contrasts: bool,
) -> Result<(u64, u64)> {
    let mut lease = JobLease::try_acquire(db, stage, None, &[])?.ok_or_else(|| {
        anyhow::anyhow!("another analysis worker owns this project; inspect `rnaa jobs list`")
    })?;
    let result = run_deseq2_stage_owned(
        paths,
        &lease.database(),
        config,
        design,
        contrasts,
        stage,
        store_contrasts,
    );
    match result {
        Ok(value) => {
            lease.finish(JobOutcome::Succeeded, None)?;
            Ok(value)
        }
        Err(error) => {
            let _ = lease.database().set_project_stage_state(
                stage,
                ProjectStageState::Failed,
                Some(&format!("{error:#}")),
            );
            lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))?;
            Err(error)
        }
    }
}

fn run_deseq2_stage_owned(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    design_override: Option<String>,
    contrasts: Vec<ContrastSpec>,
    stage_name: &str,
    store_contrasts: bool,
) -> Result<(u64, u64)> {
    let current_paths = paths.with_generation_heads(db)?;
    let paths = &current_paths;
    rnaa_core::samplesheet::publication::recover_view(db, paths, Default::default())?;
    db.set_project_stage_state(stage_name, ProjectStageState::Pending, None)?;
    let runs = db.list_active_runs_in_states(&[
        RunState::QuantDone,
        RunState::DeRunning,
        RunState::DeFailed,
        RunState::DeDone,
        RunState::CorrRunning,
        RunState::CorrFailed,
        RunState::CorrDone,
        RunState::Cleaned,
    ])?;
    if runs.is_empty() {
        println!("{stage_name}\tno_quantified_runs");
        return Ok((0, 0));
    }

    let canonical_samplesheet_path = paths.samplesheet_path();
    if !canonical_samplesheet_path.exists() {
        bail!(
            "samplesheet missing at {}; run `rnaa resolve` first",
            canonical_samplesheet_path.display()
        );
    }
    let header = read_samplesheet_header(&canonical_samplesheet_path)?;
    let design = design_override.unwrap_or_else(|| config.deseq2.design.clone());
    validate_design_formula(&header, &design, &config.metadata)?;
    for contrast in &contrasts {
        if !header.contains(&contrast.factor) {
            bail!(
                "contrast factor '{}' not found in samplesheet columns",
                contrast.factor
            );
        }
    }

    config.deseq2.design = design.clone();
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    if store_contrasts {
        db.replace_contrasts(&contrasts)?;
    }

    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, db, &config)?;
    let reference_guard = reference::pin(db, &reference, &config, ReferenceScope::Analysis)?;
    let validated_runs = validate_deseq2_quant_inputs(db, paths, &config, &reference, &runs)?;
    let eligibility = quality::Snapshot::capture(db, &config)?;
    let valid_ids = validated_runs
        .into_iter()
        .map(|run| run.run_accession)
        .collect::<std::collections::BTreeSet<_>>();
    let validated_runs = db
        .list_active_runs()?
        .into_iter()
        .filter(|run| {
            valid_ids.contains(&run.run_accession)
                && eligibility.eligible.contains(&run.run_accession)
        })
        .collect::<Vec<_>>();
    eligibility.verify_catalog(db, &config)?;
    if validated_runs.is_empty() {
        bail!(
            "no valid quantified runs match current reference, scientific settings and QC policy for {stage_name}; inspect run errors"
        );
    }
    if validated_runs.len() < 2 {
        bail!(
            "need at least 2 valid quantified runs for {stage_name}; found {}",
            validated_runs.len()
        );
    }
    let request_inputs = tempfile::Builder::new()
        .prefix("rnaa-de-inputs-")
        .tempdir_in(&paths.support_dir)?;
    let filtered_samplesheet_path = request_inputs.path().join("samplesheet_quantified.tsv");
    let metadata_source =
        write_filtered_deseq2_samplesheet(paths, &validated_runs, &filtered_samplesheet_path)?;
    if let Some((outputs, recorded)) = reconcile_deseq2_stage(
        paths,
        db,
        &validated_runs,
        &reference,
        &design,
        &contrasts,
        &config,
        store_contrasts,
        &filtered_samplesheet_path,
        &eligibility,
        &reference_guard,
        &metadata_source,
    )? {
        return Ok((outputs, recorded));
    }

    let mut analysis_request = deseq_request(
        paths,
        &filtered_samplesheet_path,
        &reference,
        &design,
        &contrasts,
        &config,
        stage_name,
    )?;
    eligibility.bind(&mut analysis_request)?;
    bind_metadata(&mut analysis_request, &metadata_source)?;
    reference::bind(&mut analysis_request, &reference_guard)?;
    reference_guard.verify_current(db, &config, ReferenceScope::Analysis)?;
    eligibility.verify_current(db, &config)?;
    analysis_request.verify_current_settings(db)?;
    db.set_job_request(&analysis_request.request.digest()?)?;

    let mut generation = Generation::begin(db, OutputGroup::De, stage_name)?;
    let output_paths = generation.paths(paths)?;
    let paths = &output_paths;
    let generation_sheet = paths
        .de_project_dir(db.project_id())
        .join("samplesheet_quantified.tsv");
    fs::create_dir_all(generation_sheet.parent().unwrap())?;
    fs::copy(&filtered_samplesheet_path, &generation_sheet)?;
    let filtered_samplesheet_path = generation_sheet;
    analysis_request.relocate_input("SAMPLESHEET", &filtered_samplesheet_path)?;
    crate::analysis::inputs::freeze_metadata_source(
        &mut analysis_request,
        &metadata_source,
        &generation.head().payload_path,
    )?;
    let runner = RustDeseq2Runner;
    let target_runs = validated_runs
        .iter()
        .map(|run| run.run_accession.clone())
        .collect::<Vec<_>>();
    db.set_project_stage_state(stage_name, ProjectStageState::Running, None)?;
    db.append_event(
        stage_name,
        None,
        "stage started",
        serde_json::json!({
            "design": design,
            "contrasts": contrasts,
            "runs": target_runs.len(),
        }),
    )?;

    let outputs = match runner
        .deseq2(
            db.project_id(),
            &filtered_samplesheet_path,
            &reference,
            &design,
            &contrasts,
            paths,
            &config,
        )
        .and_then(|outputs| {
            analysis_request.seal(
                &paths
                    .de_project_dir(db.project_id())
                    .join("deseq2_manifest.json"),
                &outputs,
            )?;
            Ok(outputs)
        }) {
        Ok(outputs) => outputs,
        Err(err) => {
            let message = format!("{err:#}");
            db.set_project_stage_state(stage_name, ProjectStageState::Failed, Some(&message))?;
            db.append_event(
                stage_name,
                None,
                "stage failed",
                serde_json::json!({ "error": message }),
            )?;
            return Err(err);
        }
    };

    let de_request_manifest = paths
        .de_project_dir(db.project_id())
        .join("rnaa_deseq2_request.json");
    analysis_request.write_descriptor(&de_request_manifest)?;
    let mut registered_outputs = outputs.clone();
    registered_outputs.push(de_request_manifest.clone());
    let records = prepare_analysis_artifacts(db, &registered_outputs, OutputGroup::De)?;
    let recorded = records.len() as u64;
    reference_guard.verify_current(db, &config, ReferenceScope::Analysis)?;
    eligibility.verify_current(db, &config)?;
    analysis_request.verify_current_settings(db)?;
    analysis_request.verify_inputs()?;
    generation.seal_request(
        &registered_outputs,
        &analysis_request.request,
        &analysis_request.inputs,
    )?;
    generation.commit(&records, |db| {
        analysis_request.verify_settings_catalog(db)?;
        reference_guard.verify_catalog(db, &config, ReferenceScope::Analysis)?;
        eligibility.verify_catalog(db, &config)?;
        db.set_project_stage_state(stage_name, ProjectStageState::Done, None)?;
        db.append_event(
            stage_name,
            None,
            "stage completed",
            serde_json::json!({
                "outputs_recorded": recorded,
                "contributing_runs": target_runs.len(),
            }),
        )?;

        Ok(())
    })?;

    Ok((outputs.len() as u64, recorded))
}

#[allow(clippy::too_many_arguments)]
fn reconcile_deseq2_stage(
    paths: &ProjectPaths,
    db: &Database,
    runs: &[RunRecord],
    reference: &rnaa_core::model::ReferenceBundle,
    design: &str,
    contrasts: &[ContrastSpec],
    config: &ProjectConfig,
    store_contrasts: bool,
    samplesheet_path: &Path,
    eligibility: &quality::Snapshot,
    reference_guard: &ReferenceSelection,
    metadata_source: &rnaa_core::ManifestArtifact,
) -> Result<Option<(u64, u64)>> {
    let stage = if store_contrasts {
        "deseq2"
    } else {
        "normalize"
    };
    let mut request = deseq_request(
        paths,
        samplesheet_path,
        reference,
        design,
        contrasts,
        config,
        stage,
    )?;
    eligibility.bind(&mut request)?;
    bind_metadata(&mut request, metadata_source)?;
    reference::bind(&mut request, reference_guard)?;
    super::reuse::AnalysisReuse {
        paths,
        db,
        config,
        runs,
        eligibility,
        reference: reference_guard,
        group: OutputGroup::De,
        stage,
    }
    .reconcile(&request, |paths| {
        let scientific_manifest_path = paths
            .de_project_dir(db.project_id())
            .join("deseq2_manifest.json");
        let mut required = vec![
            paths
                .de_project_dir(db.project_id())
                .join("model_design.json"),
            paths
                .de_project_dir(db.project_id())
                .join("sample_metadata.tsv"),
            paths
                .de_project_dir(db.project_id())
                .join("gene_counts.tsv"),
            paths
                .de_project_dir(db.project_id())
                .join("gene_counts_annotated.tsv"),
            paths
                .de_project_dir(db.project_id())
                .join("gene_norm_counts.tsv"),
            paths
                .de_project_dir(db.project_id())
                .join("gene_norm_counts_annotated.tsv"),
            paths.de_project_dir(db.project_id()).join("vst.tsv"),
            paths
                .de_project_dir(db.project_id())
                .join("vst_annotated.tsv"),
            paths
                .de_project_dir(db.project_id())
                .join("deseq2_manifest.json"),
        ];
        if paths
            .de_project_dir(db.project_id())
            .join("vst.rds")
            .exists()
        {
            required.push(paths.de_project_dir(db.project_id()).join("vst.rds"));
        }
        for contrast in contrasts {
            required.push(
                paths
                    .de_project_dir(db.project_id())
                    .join(format!("de_{}.tsv", contrast.name)),
            );
            required.push(
                paths
                    .de_project_dir(db.project_id())
                    .join(format!("de_{}_annotated.tsv", contrast.name)),
            );
        }
        super::reuse::OutputContract {
            descriptor: paths
                .de_project_dir(db.project_id())
                .join("rnaa_deseq2_request.json"),
            manifest: scientific_manifest_path,
            required,
        }
    })
}
