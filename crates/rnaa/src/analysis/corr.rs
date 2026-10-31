//! Correlation and enrichment request snapshots, execution, reuse and publication.
use crate::analysis::annotations::annotate_corr_outputs;
use crate::analysis::geneset::{materialize_geneset_matrix, parse_output_mode};
use crate::analysis::inputs::bind_metadata;
use crate::analysis::prepare_analysis_artifacts;
use crate::cli::CorrArgs;
use crate::ingest::inputs::recover_interrupted_work;
use crate::resources::resolve_scheduler_budget;
use crate::{quality, reference};
use anyhow::Context;
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::generations::{Generation, OutputGroup};
use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::model::{CorrelationMethod, OutputMode, RunRecord};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::references::{ReferenceScope, ReferenceSelection};
use rnaa_core::resources::{DiskRequest, ExecutionResources, ReservationOwner, ResourceRequest};
use rnaa_core::state::{ProjectStageState, RunState};
use rnaa_core::traits::ReferenceManager;
use rnaa_stages::analysis_request::{
    AnalysisRequest, CorrOptions, correlation_request, correlation_request_exact,
};
use rnaa_stages::corr::{
    CorrelationLimits, MinCorrCorrelator, ModuleGseaParams, adjust_matrix_streaming,
    estimate_adjustment, run_module_gsea_with_limits, validate_module_gsea,
};
use rnaa_stages::refs::EnsemblReferenceManager;
use serde_json::json;
use std::fs;
use std::path::{Path, PathBuf};

fn bind_matrix_header(
    request: &mut AnalysisRequest<CorrOptions>,
    witness: Option<&rnaa_core::ManifestArtifact>,
) -> Result<()> {
    if let Some(witness) = witness {
        let mut matrices = request.inputs.iter().filter(|input| input.kind == "MATRIX");
        if matrices.next() != Some(witness)
            || matrices.next().is_some()
            || !witness.matches_file()?
        {
            bail!("custom matrix changed while preparing its sample metadata");
        }
        request.observe_input(witness.clone());
    }
    Ok(())
}

fn capture_corr_request(
    matrix: &Path,
    sheet: &Path,
    annotation: Option<&Path>,
    config: &ProjectConfig,
    modules: Option<&ModuleGseaParams>,
    matrix_witness: Option<&rnaa_core::ManifestArtifact>,
) -> Result<AnalysisRequest<CorrOptions>> {
    let mut request = match matrix_witness {
        Some(witness) => {
            correlation_request_exact(matrix, sheet, annotation, config, modules, witness.bytes)?
        }
        None => correlation_request(matrix, sheet, annotation, config, modules)?,
    };
    bind_matrix_header(&mut request, matrix_witness)?;
    Ok(request)
}

fn corr_module_params(args: &CorrArgs) -> Option<ModuleGseaParams> {
    args.module_gmt.as_ref().map(|path| ModuleGseaParams {
        gmt_path: path.clone(),
        top_genes: args.module_top,
        min_module_size: args.module_min_size,
        permutations: args.fgsea_permutations,
        seed: args.fgsea_seed,
        pathway_min_size: args.fgsea_min_size,
        pathway_max_size: args.fgsea_max_size,
    })
}

#[allow(clippy::too_many_arguments)]
fn reconcile_corr_stage(
    paths: &ProjectPaths,
    db: &Database,
    runs: &[RunRecord],
    reference: &rnaa_core::ReferenceBundle,
    matrix_path: &Path,
    config: &ProjectConfig,
    args: &CorrArgs,
    eligibility: &quality::Snapshot,
    reference_guard: &ReferenceSelection,
    metadata_sheet: &Path,
    metadata_source: &rnaa_core::ManifestArtifact,
    matrix_witness: Option<&rnaa_core::ManifestArtifact>,
) -> Result<Option<(u64, u64)>> {
    let module_params = corr_module_params(args);
    let mut request = capture_corr_request(
        matrix_path,
        metadata_sheet,
        reference.gene_annotation_path.as_deref(),
        config,
        module_params.as_ref(),
        matrix_witness,
    )?;
    eligibility.bind(&mut request)?;
    bind_metadata(&mut request, metadata_source)?;
    reference::bind(&mut request, reference_guard)?;
    let output_mode = parse_output_mode(&config.corr.output)?;
    super::reuse::AnalysisReuse {
        paths,
        db,
        config,
        runs,
        eligibility,
        reference: reference_guard,
        group: OutputGroup::Corr,
        stage: "corr",
    }
    .reconcile(&request, |paths| {
        let corr_dir = paths.corr_project_dir(db.project_id());
        let mut outputs = vec![
            paths.corr_dir.join("adjusted_matrix.tsv"),
            paths.manifest_path("corr", db.project_id()),
        ];
        match &output_mode {
            OutputMode::TopK { .. } => {
                outputs.push(corr_dir.join("edges_topk.tsv"));
                outputs.push(corr_dir.join("edges_topk_annotated.tsv"));
            }
            OutputMode::Threshold { .. } => {
                outputs.push(corr_dir.join("edges_threshold.tsv"));
                outputs.push(corr_dir.join("edges_threshold_annotated.tsv"));
            }
            OutputMode::Dense => outputs.push(corr_dir.join("corr_dense.tsv")),
        }
        if args.module_gmt.is_some() {
            outputs.push(corr_dir.join("modules_top_genes.tsv"));
            outputs.push(corr_dir.join("modules_top_genes_annotated.tsv"));
            outputs.push(corr_dir.join("module_t_values.tsv"));
            outputs.push(corr_dir.join("module_t_values_annotated.tsv"));
            outputs.push(corr_dir.join("module_fgsea.tsv"));
            outputs.push(corr_dir.join("module_fgsea_manifest.json"));
            outputs.push(corr_dir.join("module_membership.tsv"));
            outputs.push(corr_dir.join("module_eigengenes.tsv"));
        }
        if config.corr.adjust == "residualize" {
            outputs.push(paths.corr_dir.join("model_design.json"));
        }
        super::reuse::OutputContract {
            descriptor: corr_dir.join("rnaa_corr_request.json"),
            manifest: paths.manifest_path("corr", db.project_id()),
            required: outputs,
        }
    })
}

pub(crate) fn cmd_corr(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: CorrArgs,
) -> Result<()> {
    let mut lease = JobLease::try_acquire(db, "corr", None, &[])?.ok_or_else(|| {
        anyhow::anyhow!("another analysis worker owns this project; inspect `rnaa jobs list`")
    })?;
    let result = cmd_corr_owned(paths, &lease.database(), config, args, lease.attempt_id());
    match result {
        Ok(()) => lease.finish(JobOutcome::Succeeded, None),
        Err(error) => {
            let _ = lease.database().set_project_stage_state(
                "corr",
                ProjectStageState::Failed,
                Some(&format!("{error:#}")),
            );
            lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))?;
            Err(error)
        }
    }
}

fn cmd_corr_owned(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    mut args: CorrArgs,
    attempt: &str,
) -> Result<()> {
    let current_paths = paths.with_generation_heads(db)?;
    let paths = &current_paths;
    rnaa_core::samplesheet::publication::recover_view(db, paths, Default::default())?;
    recover_interrupted_work(db)?;
    db.set_project_stage_state("corr", ProjectStageState::Pending, None)?;
    if let Some(adjust) = &args.adjust {
        config.corr.adjust = adjust.clone();
    }
    if let Some(model) = &args.model {
        config.corr.model = model.clone();
    }
    if let Some(method) = &args.method {
        config.corr.method = method.clone();
    }
    if let Some(geneset) = &args.geneset {
        config.corr.geneset = geneset.clone();
    }
    if let Some(out_spec) = &args.out_spec {
        config.corr.output = out_spec.clone();
    }
    config.corr.matrix = if let Some(matrix) = &mut args.matrix {
        let selected = PathBuf::from(&*matrix)
            .canonicalize()
            .context("failed to locate correlation matrix")?;
        *matrix = selected
            .to_str()
            .context("correlation matrix path is not valid UTF-8")?
            .to_owned();
        let selected_de = paths
            .de_project_dir(db.project_id())
            .join("vst.tsv")
            .canonicalize()
            .ok();
        (selected_de.as_ref() != Some(&selected)).then_some(selected)
    } else {
        None
    };
    if let Some(path) = &mut args.module_gmt {
        *path = path
            .canonicalize()
            .context("failed to locate module pathway file")?;
    }
    let modules = corr_module_params(&args);
    if let Some(params) = &modules {
        params.validate()?;
    }
    config.corr.modules = modules.as_ref().map(Into::into);
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, db, &config)?;

    let reference_guard = reference::pin(db, &reference, &config, ReferenceScope::Analysis)?;
    let runs = db.list_active_runs_in_states(&[
        RunState::QuantDone,
        RunState::Cleaned,
        RunState::DeDone,
        RunState::CorrFailed,
        RunState::CorrRunning,
        RunState::CorrDone,
    ])?;
    if runs.is_empty() {
        bail!("no quantified runs found; quantify the cohort before correlation");
    }

    let eligibility = quality::Snapshot::capture(db, &config)?;
    let completed = runs
        .into_iter()
        .map(|run| run.run_accession)
        .collect::<std::collections::BTreeSet<_>>();
    let mut runs = db
        .list_active_runs()?
        .into_iter()
        .filter(|run| {
            completed.contains(&run.run_accession)
                && eligibility.eligible.contains(&run.run_accession)
        })
        .collect::<Vec<_>>();
    eligibility.verify_catalog(db, &config)?;
    if runs.is_empty() {
        bail!("no quantified runs are currently eligible; reevaluate QC before correlation");
    }
    let metadata_staging = tempfile::Builder::new()
        .prefix("rnaa-corr-metadata-")
        .tempdir_in(&paths.support_dir)?;
    let matrix_path = args
        .matrix
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or_else(|| paths.de_project_dir(db.project_id()).join("vst.tsv"));
    if !matrix_path.exists() {
        bail!(
            "missing matrix at {}; run `rnaa deseq2` first or pass --matrix",
            matrix_path.display()
        );
    }
    let matrix_path = matrix_path
        .canonicalize()
        .context("failed to pin the matrix's canonical input path")?;

    // An explicit path to the selected DE matrix has the same producer proof
    // as the default. Its biological columns need not include inactive history.
    let current_de_matrix = paths.de_project_dir(db.project_id()).join("vst.tsv");
    let uses_de_matrix = args.matrix.is_none()
        || matrix_path
            .canonicalize()
            .ok()
            .zip(current_de_matrix.canonicalize().ok())
            .is_some_and(|(given, selected)| given == selected);

    if uses_de_matrix && !db.verify_generation(OutputGroup::De, db.project_id())? {
        bail!(
            "current DE generation failed integrity verification; rerun deseq2 or normalize before correlation"
        );
    }

    let mut matrix_witness = None;
    let (metadata_sheet, metadata_source) = if uses_de_matrix {
        let descriptor = crate::analysis::producer::load_descriptor(paths, db, &reference_guard)?;
        let admitted = crate::analysis::inputs::validate_deseq2_quant_inputs(
            db, paths, &config, &reference, &runs,
        )?;
        if admitted.len() != runs.len() {
            bail!(
                "DE count dependencies are stale; recover/recompute quantification and DE before correlation"
            );
        }
        crate::analysis::producer::current_metadata(
            paths,
            db,
            &config,
            &eligibility,
            descriptor,
            &runs,
            metadata_staging.path(),
        )?
    } else {
        let witness = rnaa_core::ManifestArtifact::snapshot("MATRIX", &matrix_path)?;
        let columns = rnaa_formats::matrix::read_matrix_columns_exact(&matrix_path, witness.bytes)?;
        if !witness.matches_file()? {
            bail!("custom matrix changed while preparing its sample metadata");
        }
        let prepared = crate::analysis::custom_matrix::prepare(
            db,
            paths,
            &config,
            &runs,
            &columns,
            metadata_staging.path(),
        )?;
        runs = prepared.runs;
        matrix_witness = Some(witness);
        (prepared.sheet, prepared.source)
    };
    if let Some((_outputs, _recorded)) = reconcile_corr_stage(
        paths,
        db,
        &runs,
        &reference,
        &matrix_path,
        &config,
        &args,
        &eligibility,
        &reference_guard,
        &metadata_sheet,
        &metadata_source,
        matrix_witness.as_ref(),
    )? {
        return Ok(());
    }

    let module_params = corr_module_params(&args);
    let mut analysis_request = capture_corr_request(
        &matrix_path,
        &metadata_sheet,
        reference.gene_annotation_path.as_deref(),
        &config,
        module_params.as_ref(),
        matrix_witness.as_ref(),
    )?;
    eligibility.bind(&mut analysis_request)?;
    bind_metadata(&mut analysis_request, &metadata_source)?;
    reference::bind(&mut analysis_request, &reference_guard)?;
    reference_guard.verify_current(db, &config, ReferenceScope::Analysis)?;
    eligibility.verify_current(db, &config)?;
    analysis_request.verify_current_settings(db)?;
    db.set_job_request(&analysis_request.request.digest()?)?;

    // Selection currently has cooperative CPU admission only.
    let resources = ExecutionResources::open(
        &config.resources,
        ReservationOwner {
            project_root: fs::canonicalize(&paths.root)?,
            attempt_id: attempt.into(),
            stage: "corr".into(),
            run_accession: None,
        },
    )?;
    let threads = resolve_scheduler_budget(&config)
        .corr_threads
        .min(usize::try_from(resources.limits().cpu_threads).unwrap_or(usize::MAX));
    let mut cpu = Some(resources.reserve(
        "compute",
        ResourceRequest {
            cpu_threads: threads as u64,
            memory_bytes: 0,
            disks: vec![],
        },
        &|| db.check_job(),
    )?);
    db.check_job()?;

    let mut generation = Generation::begin(db, OutputGroup::Corr, "corr")?;
    let mut output_paths = generation.paths(paths)?;
    fs::create_dir_all(&output_paths.corr_dir)?;
    output_paths.metadata_dir = generation.head().payload_path.join("metadata");
    fs::create_dir_all(&output_paths.metadata_dir)?;
    let frozen_metadata = output_paths.metadata_dir.join(
        metadata_sheet
            .file_name()
            .context("metadata sheet has no file name")?,
    );
    fs::copy(&metadata_sheet, &frozen_metadata)?;
    let paths = &output_paths;
    let corr_request_manifest = paths
        .corr_project_dir(db.project_id())
        .join("rnaa_corr_request.json");
    analysis_request.relocate_input("METADATA", &frozen_metadata)?;
    crate::analysis::inputs::freeze_metadata_source(
        &mut analysis_request,
        &metadata_source,
        &generation.head().payload_path,
    )?;
    db.set_project_stage_state("corr", ProjectStageState::Running, None)?;
    db.append_event(
        "corr",
        None,
        "corr stage started",
        serde_json::json!({
            "matrix": matrix_path.display().to_string(),
            "model": config.corr.model,
            "method": config.corr.method,
            "output": config.corr.output,
            "geneset": config.corr.geneset,
            "gene_annotation": reference
                .gene_annotation_path
                .as_ref()
                .map(|path| path.display().to_string()),
            "module_gmt": args.module_gmt.as_ref().map(|path| path.display().to_string()),
            "module_top": args.module_top,
            "module_min_size": args.module_min_size,
            "fgsea_permutations": args.fgsea_permutations,
        }),
    )?;

    let computation: Result<(Vec<PathBuf>, rnaa_core::AdjustedMatrix)> = (|| {
        let selected_matrix =
            materialize_geneset_matrix(paths, &matrix_path, &config.corr.geneset, &|| {
                db.check_job()
            })?;
        let residualize = match config.corr.adjust.as_str() {
            "residualize" => true,
            "none" => false,
            other => bail!("unsupported adjustment '{other}'; expected none or residualize"),
        };
        let adjustment_estimate =
            estimate_adjustment(&selected_matrix, &frozen_metadata, residualize)?;
        drop(cpu.take());
        let adjustment = resources.reserve(
            "adjust",
            ResourceRequest {
                cpu_threads: threads as u64,
                memory_bytes: adjustment_estimate.memory_bytes,
                disks: vec![DiskRequest {
                    path: paths.corr_dir.clone(),
                    bytes: adjustment_estimate.output_bytes,
                }],
            },
            &|| db.check_job(),
        )?;
        if estimate_adjustment(&selected_matrix, &frozen_metadata, residualize)?
            != adjustment_estimate
        {
            bail!("adjustment input or resource estimate changed after shared admission");
        }
        let mut adjusted = adjust_matrix_streaming(
            &selected_matrix,
            &frozen_metadata,
            &config.corr.adjust,
            &config.corr.model,
            paths,
            &config.metadata,
            adjustment_estimate,
            &|| db.check_job(),
        )?;
        // Later stages consume the matrix by path; do not retain its identifiers
        // while waiting for the next phase's independent resource claim.
        adjusted.genes = Vec::new();
        adjusted.samples = Vec::new();
        drop(adjustment);
        cpu = Some(resources.reserve(
            "prepare-network",
            ResourceRequest {
                cpu_threads: threads as u64,
                memory_bytes: 0,
                disks: vec![],
            },
            &|| db.check_job(),
        )?);

        let module_params = corr_module_params(&args);
        let module_limits = CorrelationLimits::from(&config.corr);
        let module_threads = threads;
        if let Some(params) = &module_params {
            validate_module_gsea(&adjusted, params, module_limits, module_threads)?;
        }
        let method = config.corr.method.parse::<CorrelationMethod>()?;
        let output_mode = parse_output_mode(&config.corr.output)?;
        let estimate =
            module_limits.preflight(&adjusted.matrix_path, method, &output_mode, threads)?;
        let memory_bytes = u64::try_from(
            estimate
                .working_bytes
                .checked_add(64 * 1024 * 1024)
                .context("correlation memory estimate overflow")?,
        )
        .context("correlation memory estimate exceeds resource accounting")?;
        let output_bytes = u64::try_from(estimate.output_bytes)
            .context("correlation output estimate exceeds resource accounting")?;
        // Never wait for the larger claim while holding the preliminary CPU claim.
        drop(cpu.take());
        let network = resources.reserve(
            "network",
            ResourceRequest {
                cpu_threads: threads as u64,
                memory_bytes,
                disks: vec![DiskRequest {
                    path: paths.corr_dir.clone(),
                    bytes: output_bytes,
                }],
            },
            &|| db.check_job(),
        )?;
        let correlator = MinCorrCorrelator
            .with_limits(CorrelationLimits::from(&config.corr))
            .with_admission(estimate);
        let mut outputs = correlator.correlate_checked(
            &adjusted,
            method,
            &output_mode,
            threads,
            paths,
            db.project_id(),
            &|| db.check_job(),
        )?;
        drop(network);
        cpu = Some(resources.reserve(
            "finish",
            ResourceRequest {
                cpu_threads: threads as u64,
                memory_bytes: 0,
                disks: vec![],
            },
            &|| db.check_job(),
        )?);
        if let Some(params) = &module_params {
            let module_outputs = run_module_gsea_with_limits(
                &adjusted,
                paths,
                db.project_id(),
                params,
                module_limits,
                module_threads,
            )
            .with_context(|| {
                format!(
                    "module enrichment failed for gmt {}",
                    params.gmt_path.display()
                )
            })?;
            outputs.extend(module_outputs);
        }
        if let Some(annotation_path) = &reference.gene_annotation_path {
            let annotated = annotate_corr_outputs(&outputs, annotation_path)?;
            outputs.extend(annotated);
        }

        outputs.push(adjusted.matrix_path.clone());
        if config.corr.adjust == "residualize" {
            outputs.push(paths.corr_dir.join("model_design.json"));
        }
        analysis_request.seal(&paths.manifest_path("corr", db.project_id()), &outputs)?;
        analysis_request.write_descriptor(&corr_request_manifest)?;

        Ok((outputs, adjusted))
    })();
    let (outputs, _adjusted) = match computation {
        Ok(result) => result,
        Err(error) => {
            let message = format!("{error:#}");
            db.set_project_stage_state("corr", ProjectStageState::Failed, Some(&message))?;
            db.append_event("corr", None, "corr stage failed", json!({"error":message}))?;
            return Err(error);
        }
    };
    let mut registered_outputs = outputs.clone();
    registered_outputs.push(corr_request_manifest.clone());
    let records = prepare_analysis_artifacts(db, &registered_outputs, OutputGroup::Corr)?;
    let recorded = records.len() as u64;
    db.check_job()?;
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
        db.set_project_stage_state("corr", ProjectStageState::Done, None)?;
        db.append_event(
            "corr",
            None,
            "corr stage completed",
            serde_json::json!({
                "outputs_recorded": recorded,
                "contributing_runs": runs.len()
            }),
        )?;
        Ok(())
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn custom_header_witness_rejects_same_size_matrix_replacement() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let matrix = dir.path().join("matrix.tsv");
        let sheet = dir.path().join("sheet.tsv");
        fs::write(&matrix, "gene_id\tA\tB\ng1\t1\t2\ng2\t3\t4\n")?;
        fs::write(&sheet, "sample_id\nA\nB\n")?;
        let witness = rnaa_core::ManifestArtifact::snapshot("MATRIX", &matrix)?;
        fs::write(&matrix, "gene_id\tA\tC\ng1\t1\t2\ng2\t3\t4\n")?;
        let mut config = ProjectConfig::default();
        config.corr.adjust = "none".into();
        let mut request = correlation_request(&matrix, &sheet, None, &config, None)?;
        assert!(bind_matrix_header(&mut request, Some(&witness)).is_err());
        let current = rnaa_core::ManifestArtifact::snapshot("MATRIX", &matrix)?;
        bind_matrix_header(&mut request, Some(&current))?;
        fs::write(&matrix, "gene_id\tA\tD\ng1\t1\t2\ng2\t3\t4\n")?;
        assert!(request.verify_inputs().is_err());
        Ok(())
    }
}
