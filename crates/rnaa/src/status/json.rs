//! Structured status assembled from the same durable projections as text status.
use crate::progress::compute_progress_snapshot;
use crate::storage::{configured_storage_limit_bytes, project_active_storage_bytes};
use crate::{quality, reference};
use anyhow::Result;
use rnaa_core::{Database, config::ProjectConfig, paths::ProjectPaths, util::dir_size};
use serde_json::{Value, json};
use std::{collections::BTreeMap, time::Duration};

fn seconds(value: Option<Duration>) -> Option<f64> {
    value.map(|duration| duration.as_secs_f64())
}

pub(super) fn print_status(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
) -> Result<()> {
    let reference = reference::refresh_status(db, config)?;
    let selected = reference.selection.as_ref().map(|value| {
        json!({
            "reference_id": value.reference_id,
            "generation_id": value.generation.generation_id,
            "index_sha256": value.index.checksum,
        })
    });
    let project = db.project()?;
    let catalog_runs = db.list_runs()?;
    let cohort = db.active_cohort()?;
    let runs = db
        .verify_cohort_config(config)
        .and_then(|()| db.list_active_runs());
    let (cohort_error, active_runs, progress, qc_failed, state_counts) = match runs {
        Ok(runs) => {
            let artifacts = db.list_artifacts()?;
            let events = db.list_events(10_000)?;
            let stages = db.list_project_stage_states()?;
            let progress = compute_progress_snapshot(&runs, &artifacts, &events, &stages, config);
            let mut counts = BTreeMap::new();
            let mut active = Vec::with_capacity(runs.len());
            for run in &runs {
                *counts.entry(run.state.as_str()).or_insert(0_u64) += 1;
                let eligibility = db.eligibility(&run.run_accession)?;
                active.push(json!({"run_accession":run.run_accession, "sample_accession":run.sample_accession,
                    "state":run.state.as_str(), "last_error":run.last_error, "updated_at":run.updated_at,
                    "eligibility":eligibility}));
            }
            let result = json!({
                "total_runs":progress.total_runs, "download_done":progress.download_done,
                "preprocess_enabled":progress.preprocess_enabled, "preprocess_done":progress.preprocess_done,
                "quant_done":progress.quant_done, "normalize_done":progress.normalize_done,
                "correlation_done":progress.corr_stage_done,
                "downloaded_bytes":progress.downloaded_bytes, "expected_total_download_bytes":progress.expected_total_bytes,
                "download_rate_bytes_per_second":(progress.download_rate_bps > 0.0).then_some(progress.download_rate_bps),
                "quant_rate_runs_per_second":(progress.quant_rate_runs_per_sec > 0.0).then_some(progress.quant_rate_runs_per_sec),
                "download_eta_seconds":seconds(progress.download_eta), "processing_eta_seconds":seconds(progress.processing_eta),
                "total_eta_seconds":seconds(progress.total_eta), "active_stage":progress.active_stage,
                "active_detail":progress.active_detail, "limiting_stage":progress.limiting_stage,
            });
            (
                None,
                Some(active),
                Some(result),
                Some(quality::blocked_count(db, config)?),
                Some(counts),
            )
        }
        Err(error) => (Some(format!("{error:#}")), None, None, None, None),
    };
    let stages: Vec<Value> = db
        .list_project_stage_states()?
        .into_iter()
        .map(|stage| {
            json!({
                "stage":stage.stage, "state":stage.state.as_str(), "last_error":stage.last_error,
                "updated_at":stage.updated_at
            })
        })
        .collect();
    let errors: Vec<Value> = db
        .last_errors(10)?
        .into_iter()
        .map(|(run, state, error)| json!({"run_accession":run,"state":state,"error":error}))
        .collect();
    let report = json!({
        "schema_version":1,
        "project":{"id":project.project_id,"root":project.root_dir,"inputs":db.list_inputs()?.len(),
            "metadata_overrides":db.list_metadata_overrides()?.len(),"samplesheet":paths.samplesheet_path()},
        "reference":{"ready":reference.reason.is_none(),"selected":selected,"reason":reference.reason},
        "cohort":{"ready":cohort_error.is_none(),"reason":cohort_error,"active":cohort},
        "runs_catalog_total":catalog_runs.len(),"active_runs":active_runs,"state_counts":state_counts,
        "progress":progress,"qc_failed":qc_failed,"qc_records":db.list_preprocess_qc()?.len(),
        "stages":stages,"generation_heads":db.generation_heads()?,"recent_errors":errors,
        "storage":{"disk_usage_bytes":db.disk_usage_bytes()?.max(dir_size(&paths.raw_dir)),
            "active_bytes":project_active_storage_bytes(paths),"limit_bytes":configured_storage_limit_bytes(config)?},
        "quant_match_cache_mb_per_run":config.quant.match_cache_mb,
        "consistency":"observational; reconciliation may update state and concurrent work may change it"
    });
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(())
}
