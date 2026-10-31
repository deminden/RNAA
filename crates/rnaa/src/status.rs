//! Human-readable project status and environment diagnostics.
mod json;
use crate::ingest::inputs::{reconcile_verified_download_state, recover_interrupted_work};
use crate::progress::{compute_progress_snapshot, format_eta, ratio_pct};
use crate::storage::{configured_storage_limit_bytes, project_active_storage_bytes};
use crate::{quality, reference};
use anyhow::Result;
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::util::dir_size;

pub(crate) fn cmd_status(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    json: bool,
) -> Result<()> {
    recover_interrupted_work(db)?;
    if db
        .verify_cohort_config(config)
        .and_then(|()| db.list_active_runs())
        .is_ok()
    {
        reconcile_verified_download_state(db)?;
        let _ = db.with_idle_analysis(|db| quality::refresh_all(db, config))?;
    }
    if json {
        return json::print_status(paths, db, config);
    }
    for head in db.generation_heads()? {
        println!(
            "generation\t{}\t{}\t{}\t{}",
            head.output_group.as_str(),
            head.subject,
            head.generation_id,
            head.payload_path.display()
        );
    }
    let reference = reference::refresh_status(db, config)?;
    println!(
        "reference_ready\t{}",
        if reference.reason.is_none() {
            "yes"
        } else {
            "no"
        }
    );
    if let Some(selected) = reference.selection {
        println!(
            "selected_reference\t{}\t{}",
            selected.reference_id, selected.generation.generation_id
        );
        println!("reference_index_sha256\t{}", selected.index.checksum);
    }
    if let Some(reason) = reference.reason {
        println!("reference_reason\t{reason}");
    }
    let project = db.project()?;
    let inputs = db.list_inputs()?;
    let overrides = db.list_metadata_overrides()?;
    let disk_usage = db.disk_usage_bytes()?;
    let active_storage_bytes = project_active_storage_bytes(paths);
    let storage_limit_bytes = configured_storage_limit_bytes(config)?;
    let catalog_runs = db.list_runs()?;
    let (runs, cohort_error) = match db
        .verify_cohort_config(config)
        .and_then(|()| db.list_active_runs())
    {
        Ok(runs) => (runs, None),
        Err(error) => (Vec::new(), Some(format!("{error:#}"))),
    };
    println!("runs_catalog_total\t{}", catalog_runs.len());
    println!(
        "cohort_ready\t{}",
        if cohort_error.is_none() { "yes" } else { "no" }
    );
    if let Some(cohort) = db.active_cohort()? {
        println!("active_cohort\t{}", cohort.cohort_id);
        for member in cohort.members {
            println!(
                "cohort_member\t{}\t{}\t{:?}",
                member.run, member.sample.id, member.sample.source
            );
        }
    }
    if let Some(error) = &cohort_error {
        println!("cohort_reason\t{error}");
    }
    let artifacts = db.list_artifacts()?;
    let events = db.list_events(10_000)?;
    let stages = db.list_project_stage_states()?;
    for stage in &stages {
        println!(
            "cohort_stage\t{}\t{}\t{}",
            stage.stage, stage.state, stage.updated_at
        );
        if let Some(error) = &stage.last_error {
            println!(
                "cohort_stage_error\t{}\t{}",
                stage.stage,
                error.replace(['\n', '\r', '\t'], " ")
            );
        }
    }
    println!("run_progress_model\tingest_and_quant_only");
    println!(
        "quant_match_cache_mb_per_run\t{}",
        config.quant.match_cache_mb
    );
    let legacy_states = db.legacy_run_analysis_states()?;
    println!("legacy_run_analysis_states\t{}", legacy_states.len());
    for receipt in legacy_states {
        println!(
            "legacy_run_analysis_state\t{}\t{}\t{}\t{}",
            receipt.run_accession,
            receipt.analysis_state,
            receipt.source_updated_at,
            receipt.projected_at
        );
    }
    let mut counts = std::collections::BTreeMap::new();
    for run in &runs {
        *counts.entry(run.state.as_str()).or_insert(0_usize) += 1;
    }
    let progress = compute_progress_snapshot(&runs, &artifacts, &events, &stages, config);
    let qc_failed = if cohort_error.is_none() {
        quality::blocked_count(db, config)?
    } else {
        0
    };
    let qc_skipped = qc_failed;

    println!("project_id\t{}", project.project_id);
    println!("root\t{}", project.root_dir);
    println!("inputs\t{}", inputs.len());
    println!("runs_total\t{}", progress.total_runs);
    println!("metadata_overrides\t{}", overrides.len());
    println!("samplesheet\t{}", paths.samplesheet_path().display());
    println!(
        "disk_usage_bytes\t{}",
        disk_usage.max(dir_size(&paths.raw_dir))
    );
    println!("active_storage_bytes\t{}", active_storage_bytes);
    if let Some(limit) = storage_limit_bytes {
        println!("storage_limit_bytes\t{}", limit);
        println!(
            "storage_limit_utilization_pct\t{:.1}",
            if limit == 0 {
                0.0
            } else {
                (active_storage_bytes as f64 / limit as f64) * 100.0
            }
        );
    } else {
        println!("storage_limit_bytes\tunlimited");
    }
    println!(
        "progress_download\t{}/{} ({:.1}%)",
        progress.download_done,
        progress.total_runs,
        ratio_pct(progress.download_done, progress.total_runs)
    );
    if progress.preprocess_enabled {
        println!(
            "progress_preprocess\t{}/{} ({:.1}%)",
            progress.preprocess_done,
            progress.total_runs,
            ratio_pct(progress.preprocess_done, progress.total_runs)
        );
    } else {
        println!("progress_preprocess\toff");
    }
    println!("preprocess_qc_records\t{}", db.list_preprocess_qc()?.len());
    println!("preprocess_qc_failed\t{}", qc_failed);
    println!("runs_skipped_qc\t{}", qc_skipped);
    for run in &runs {
        if let Some(decision) = db.eligibility(&run.run_accession)? {
            println!(
                "qc_eligibility\t{}\t{:?}\t{}\t{}",
                run.run_accession,
                decision.descriptor.state,
                decision.fingerprint,
                decision.descriptor.reason.as_deref().unwrap_or("")
            );
        }
    }

    println!(
        "progress_quant\t{}/{} ({:.1}%)",
        progress.quant_done,
        progress.total_runs,
        ratio_pct(progress.quant_done, progress.total_runs)
    );
    println!(
        "progress_deseq2\t{}/1 ({:.1}%)",
        progress.normalize_done as u8,
        if progress.normalize_done { 100.0 } else { 0.0 }
    );
    println!(
        "progress_corr\t{}/1 ({:.1}%)",
        progress.corr_stage_done as u8,
        if progress.corr_stage_done { 100.0 } else { 0.0 }
    );
    println!("downloaded_bytes\t{}", progress.downloaded_bytes);
    println!(
        "expected_total_download_bytes\t{}",
        progress.expected_total_bytes
    );
    if progress.download_rate_bps > 0.0 {
        println!(
            "download_speed_mbps\t{:.2}",
            progress.download_rate_bps / (1024.0 * 1024.0)
        );
    } else {
        println!("download_speed_mbps\tunknown");
    }
    if progress.quant_rate_runs_per_sec > 0.0 {
        println!(
            "quant_rate_runs_per_hour\t{:.2}",
            progress.quant_rate_runs_per_sec * 3600.0
        );
    } else {
        println!("quant_rate_runs_per_hour\tunknown");
    }
    println!(
        "eta_download\t{}",
        format_eta(
            progress.download_eta,
            progress.download_done >= progress.total_runs
        )
    );
    println!(
        "eta_processing\t{}",
        format_eta(progress.processing_eta, progress.active_stage == "complete")
    );
    println!(
        "eta_total\t{}",
        format_eta(progress.total_eta, progress.active_stage == "complete")
    );
    println!("limiting_stage\t{}", progress.limiting_stage);
    println!("active\t{}", progress.active_detail);

    if counts.is_empty() {
        println!("states\t0");
    } else {
        println!("state_counts");
        for (state, count) in counts {
            println!("{}\t{}", state, count);
        }
    }

    let errors = db.last_errors(10)?;
    if !errors.is_empty() {
        println!("last_errors");
        for (run, state, error) in errors {
            println!("{run}\t{state}\t{error}");
        }
    }
    Ok(())
}
