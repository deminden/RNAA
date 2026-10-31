//! Progress, throughput and ETA projections from durable project records.
use crate::quality;
use crate::resources::SchedulerBudget;
use crate::run_state::{is_download_complete_state, is_preprocess_done_state, is_quant_done_state};
use anyhow::Result;
use chrono::{DateTime, Utc};
use rnaa_core::Database;
use rnaa_core::config::{DownloadPreference, ProjectConfig};
use rnaa_core::model::{
    ArtifactKind, ArtifactRecord, EventRecord, ProjectStageRecord, RemoteFile, RemoteFileKind,
    RunRecord,
};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::{ProjectStageState, RunState};
use std::time::Duration;

pub(crate) fn print_run_progress(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    scheduler: &SchedulerBudget,
    elapsed: Duration,
    last_snapshot: &mut Option<String>,
) -> Result<()> {
    let runs = db.list_active_runs()?;
    let artifacts = db.list_artifacts()?;
    let events = db.list_events(5_000)?;
    let stages = db.list_project_stage_states()?;
    let skipped_qc = quality::blocked_count(db, config)?;
    let progress = compute_progress_snapshot(&runs, &artifacts, &events, &stages, config);
    let eta_total_text = format_eta(progress.total_eta, progress.active_stage == "complete");
    let eta_download_text = format_eta(
        progress.download_eta,
        progress.download_done >= progress.total_runs,
    );
    let eta_processing_text =
        format_eta(progress.processing_eta, progress.active_stage == "complete");
    let mut lines = vec![
        format!("Project\t{}", db.project_id()),
        format!("Elapsed\t{}", format_duration(elapsed)),
        format!(
            "Stage\t{} | download {}/{} | preprocess {} | quant {}/{} | normalize {} | corr {}",
            progress.active_stage,
            progress.download_done,
            progress.total_runs,
            progress.preprocess_display(),
            progress.quant_done,
            progress.total_runs,
            if progress.normalize_done {
                "1/1"
            } else {
                "0/1"
            },
            if progress.corr_stage_done {
                "1/1"
            } else {
                "0/1"
            }
        ),
        format!(
            "ETA\t{} total | {} download | {} processing | limiting {}",
            eta_total_text, eta_download_text, eta_processing_text, progress.limiting_stage
        ),
        format!("Active\t{}", progress.active_detail),
        format!(
            "Workers\tbudget {} | download {} | preprocess up_to {}x{} | quant up_to {}x{} | deseq2 {} | corr {}{}",
            scheduler.cpu_budget,
            scheduler.download_workers,
            scheduler.preprocess_workers,
            scheduler.preprocess_threads,
            scheduler.quant_workers,
            scheduler.quant_threads,
            scheduler.deseq2_threads,
            scheduler.corr_threads,
            if scheduler.requested_worker_budget > 0 {
                String::new()
            } else {
                format!(" | reserve {}", scheduler.reserve_workers)
            }
        ),
        format!("Skipped\t{} run(s) excluded by QC gate", skipped_qc),
    ];
    if progress.download_rate_bps > 0.0 {
        lines.push(format!(
            "Speed\t{:.2} MiB/s download",
            progress.download_rate_bps / (1024.0 * 1024.0)
        ));
    }
    lines.push(format!(
        "Paths\t{} | {}",
        paths.samplesheet_path().display(),
        paths.logs_dir.display()
    ));
    let snapshot = lines.join("\n");
    if last_snapshot.as_deref() != Some(snapshot.as_str()) {
        println!("{snapshot}");
        *last_snapshot = Some(snapshot);
    }
    Ok(())
}

pub(crate) fn format_eta(value: Option<Duration>, completed: bool) -> String {
    match value {
        Some(duration) if completed => format_duration(duration),
        Some(duration) if duration > Duration::from_secs(0) => format_duration(duration),
        _ if completed => "00:00:00".to_string(),
        _ => "estimating".to_string(),
    }
}

pub(crate) struct ProgressSnapshot {
    pub(crate) total_runs: u64,
    pub(crate) download_done: u64,
    pub(crate) preprocess_enabled: bool,
    pub(crate) preprocess_done: u64,
    pub(crate) quant_done: u64,
    pub(crate) normalize_done: bool,
    pub(crate) corr_stage_done: bool,
    pub(crate) downloaded_bytes: u64,
    pub(crate) expected_total_bytes: u64,
    pub(crate) download_rate_bps: f64,
    pub(crate) quant_rate_runs_per_sec: f64,
    pub(crate) download_eta: Option<Duration>,
    pub(crate) processing_eta: Option<Duration>,
    pub(crate) total_eta: Option<Duration>,
    pub(crate) limiting_stage: &'static str,
    pub(crate) active_stage: &'static str,
    pub(crate) active_detail: String,
}

pub(crate) fn compute_progress_snapshot(
    runs: &[RunRecord],
    artifacts: &[ArtifactRecord],
    events: &[EventRecord],
    project_stages: &[ProjectStageRecord],
    config: &ProjectConfig,
) -> ProgressSnapshot {
    // Historical, excluded runs retain their records but cannot inflate the
    // active cohort's progress, byte totals or ingest throughput observations.
    let active = runs
        .iter()
        .map(|run| run.run_accession.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    let artifacts = artifacts
        .iter()
        .filter(|artifact| {
            artifact
                .run_accession
                .as_deref()
                .is_some_and(|run| active.contains(run))
        })
        .collect::<Vec<_>>();
    let events = events
        .iter()
        .filter(|event| {
            event.run_accession.as_deref().map_or(
                !matches!(event.stage.as_str(), "download" | "preprocess" | "quant"),
                |run| active.contains(run),
            )
        })
        .collect::<Vec<_>>();
    let artifacts = artifacts.as_slice();
    let events = events.as_slice();
    let total_runs = runs.len() as u64;
    let download_active = runs
        .iter()
        .filter(|run| matches!(run.state, RunState::Downloading))
        .count() as u64;
    let preprocess_enabled = config.quant.preprocess;
    let preprocess_running = runs
        .iter()
        .filter(|run| matches!(run.state, RunState::PreprocessRunning))
        .count() as u64;
    let quant_running = runs
        .iter()
        .filter(|run| matches!(run.state, RunState::QuantRunning))
        .count() as u64;
    let preprocess_done = if preprocess_enabled {
        runs.iter()
            .filter(|run| is_preprocess_done_state(run.state))
            .count() as u64
    } else {
        total_runs
    };
    let quant_done = runs
        .iter()
        .filter(|run| is_quant_done_state(run.state))
        .count() as u64;
    let normalize_done =
        total_runs > 0 && is_project_stage_done(project_stages, &["deseq2", "normalize"]);
    let corr_stage_done = total_runs > 0 && is_project_stage_done(project_stages, &["corr"]);
    let raw_downloaded_runs = artifacts
        .iter()
        .filter(|artifact| {
            matches!(
                artifact.kind,
                ArtifactKind::FastqR1
                    | ArtifactKind::FastqR2
                    | ArtifactKind::FastqSingle
                    | ArtifactKind::Sra
            )
        })
        .filter_map(|artifact| artifact.run_accession.clone())
        .collect::<std::collections::BTreeSet<_>>();
    let state_downloaded_runs = runs
        .iter()
        .filter(|run| is_download_complete_state(run.state))
        .map(|run| run.run_accession.clone())
        .collect::<std::collections::BTreeSet<_>>();
    let download_done = raw_downloaded_runs.union(&state_downloaded_runs).count() as u64;
    let downloaded_bytes = artifacts
        .iter()
        .filter(|artifact| {
            matches!(
                artifact.kind,
                ArtifactKind::FastqR1
                    | ArtifactKind::FastqR2
                    | ArtifactKind::FastqSingle
                    | ArtifactKind::Sra
            )
        })
        .map(|artifact| artifact.bytes)
        .sum::<u64>();
    let expected_total_bytes =
        estimate_total_download_bytes(runs, config.download.prefer).max(downloaded_bytes);
    let remaining_download_bytes = expected_total_bytes.saturating_sub(downloaded_bytes);
    let download_rate_bps = estimate_download_speed_bps(events);
    let preprocess_rate_runs_per_sec = estimate_stage_runs_per_sec(
        events,
        "preprocess",
        "preprocessing started",
        "preprocessing completed",
    );
    let quant_rate_runs_per_sec = estimate_stage_runs_per_sec(
        events,
        "quant",
        "quantification started",
        "quantification completed",
    );
    let de_avg_secs = estimate_stage_avg_duration_secs(
        events,
        "deseq2",
        "deseq2 stage started",
        "deseq2 stage completed",
    );
    let corr_avg_secs = estimate_stage_avg_duration_secs(
        events,
        "corr",
        "corr stage started",
        "corr stage completed",
    );
    let download_eta = duration_from_rate(remaining_download_bytes as f64, download_rate_bps);
    let preprocess_remaining_runs = total_runs.saturating_sub(preprocess_done);
    let preprocess_eta = if preprocess_enabled {
        duration_from_rate(
            preprocess_remaining_runs as f64,
            preprocess_rate_runs_per_sec,
        )
    } else {
        Some(Duration::from_secs(0))
    };
    let quant_remaining_runs = total_runs.saturating_sub(quant_done);
    let quant_eta = duration_from_rate(quant_remaining_runs as f64, quant_rate_runs_per_sec);
    let de_eta = if normalize_done {
        Some(Duration::from_secs(0))
    } else {
        de_avg_secs.map(Duration::from_secs_f64)
    };
    let corr_eta = if corr_stage_done {
        Some(Duration::from_secs(0))
    } else {
        corr_avg_secs.map(Duration::from_secs_f64)
    };
    let processing_eta = match (preprocess_eta, quant_eta, de_eta, corr_eta) {
        (Some(p), Some(q), Some(d), Some(c)) => Some(p + q + d + c),
        (Some(p), Some(q), _, _) => Some(p + q),
        (Some(p), _, Some(d), Some(c)) => Some(p + d + c),
        (Some(p), _, Some(d), _) => Some(p + d),
        (Some(p), _, _, Some(c)) => Some(p + c),
        (Some(p), _, _, _) => Some(p),
        (None, Some(q), Some(d), Some(c)) => Some(q + d + c),
        (None, Some(q), _, _) => Some(q),
        (None, _, Some(d), Some(c)) => Some(d + c),
        (None, _, Some(d), _) => Some(d),
        (None, _, _, Some(c)) => Some(c),
        _ => None,
    };
    let total_eta = match (download_eta, processing_eta) {
        (Some(d), Some(p)) => Some(std::cmp::max(d, p)),
        (Some(d), None) => Some(d),
        (None, Some(p)) => Some(p),
        (None, None) => None,
    };
    let limiting_stage = match (download_eta, processing_eta) {
        (Some(d), Some(p)) if d >= p => "download",
        (Some(_), Some(_)) => "processing",
        (Some(_), None) => "download",
        (None, Some(_)) => "processing",
        (None, None) => "unknown",
    };
    let running_analysis = project_stages.iter().find_map(|stage| {
        if stage.state != ProjectStageState::Running {
            return None;
        }
        match stage.stage.as_str() {
            "normalize" => Some("normalization"),
            "deseq2" => Some("differential expression"),
            "corr" => Some("correlation"),
            _ => None,
        }
    });
    let active_stage = if quant_running > 0 {
        "quantification"
    } else if preprocess_running > 0 {
        "preprocessing"
    } else if download_active > 0 {
        "download"
    } else if total_runs == 0 {
        "waiting"
    } else if let Some(stage) = running_analysis {
        stage
    } else if corr_stage_done {
        "complete"
    } else {
        "waiting"
    };
    let mut active_detail = describe_active_work(
        runs,
        total_runs,
        download_active,
        preprocess_enabled,
        preprocess_running,
        preprocess_done,
        quant_running,
        quant_done,
        active_stage,
    );
    if active_stage == "waiting" && total_runs > 0 {
        let failed = project_stages
            .iter()
            .filter(|stage| {
                stage.state == ProjectStageState::Failed
                    && matches!(stage.stage.as_str(), "normalize" | "deseq2" | "corr")
            })
            .map(|stage| stage.stage.as_str())
            .collect::<Vec<_>>();
        if !failed.is_empty() {
            active_detail = format!("waiting with failed cohort analysis: {}", failed.join(", "));
        } else if quant_done == total_runs {
            active_detail = if normalize_done && !corr_stage_done {
                "quantification complete; waiting for correlation".into()
            } else {
                "quantification complete; waiting for cohort analysis".into()
            };
        }
    }

    ProgressSnapshot {
        total_runs,
        download_done,
        preprocess_enabled,
        preprocess_done,
        quant_done,
        normalize_done,
        corr_stage_done,
        downloaded_bytes,
        expected_total_bytes,
        download_rate_bps,
        quant_rate_runs_per_sec,
        download_eta,
        processing_eta,
        total_eta,
        limiting_stage,
        active_stage,
        active_detail,
    }
}

#[allow(clippy::too_many_arguments)]
fn describe_active_work(
    runs: &[RunRecord],
    total_runs: u64,
    download_active: u64,
    preprocess_enabled: bool,
    preprocess_running: u64,
    preprocess_done: u64,
    quant_running: u64,
    quant_done: u64,
    active_stage: &str,
) -> String {
    match active_stage {
        "download" => format!(
            "downloading {download_active} run(s); {}/{} completed",
            runs.iter()
                .filter(|run| is_download_complete_state(run.state))
                .count(),
            total_runs
        ),
        "preprocessing" => {
            if preprocess_enabled {
                format!(
                    "preprocessing {preprocess_running} run(s); {preprocess_done}/{total_runs} complete"
                )
            } else {
                "preprocessing disabled".to_string()
            }
        }
        "quantification" => {
            format!("quantifying {quant_running} run(s); {quant_done}/{total_runs} complete")
        }
        "normalization" | "differential expression" | "correlation" => {
            format!("running {active_stage} for a cohort of {total_runs} runs")
        }
        "complete" => format!("all stages complete for {} runs", total_runs),
        _ => {
            let failed = runs
                .iter()
                .filter(|run| {
                    matches!(
                        run.state,
                        RunState::DownloadFailed
                            | RunState::VerifyFailed
                            | RunState::PreprocessFailed
                            | RunState::QuantFailed
                            | RunState::DeFailed
                            | RunState::CorrFailed
                    )
                })
                .count();
            if failed > 0 {
                format!("waiting with {failed} failed run(s) needing retry")
            } else {
                "waiting for next eligible work item".to_string()
            }
        }
    }
}

fn is_project_stage_done(stages: &[ProjectStageRecord], names: &[&str]) -> bool {
    names.iter().any(|target| {
        stages
            .iter()
            .find(|stage| stage.stage == *target)
            .is_some_and(|stage| stage.state == ProjectStageState::Done)
    })
}

pub(crate) fn ratio_pct(done: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        (done as f64) * 100.0 / (total as f64)
    }
}

fn estimate_total_download_bytes(runs: &[RunRecord], prefer: DownloadPreference) -> u64 {
    runs.iter()
        .map(|run| estimate_run_download_bytes(run, prefer))
        .sum()
}

fn estimate_run_download_bytes(run: &RunRecord, prefer: DownloadPreference) -> u64 {
    let fastq_bytes = run
        .remote_files
        .iter()
        .filter(|file| {
            matches!(
                file.kind,
                RemoteFileKind::Fastq | RemoteFileKind::FastqR1 | RemoteFileKind::FastqR2
            )
        })
        .map(|file| file.bytes.unwrap_or_default())
        .sum::<u64>();
    let sra_bytes = run
        .remote_files
        .iter()
        .filter(|file| file.kind == RemoteFileKind::Sra)
        .map(|file| file.bytes.unwrap_or_default())
        .sum::<u64>();
    let preferred = match prefer {
        DownloadPreference::Fastq if fastq_bytes > 0 => fastq_bytes,
        DownloadPreference::Sra if sra_bytes > 0 => sra_bytes,
        DownloadPreference::Fastq if sra_bytes > 0 => sra_bytes,
        DownloadPreference::Sra if fastq_bytes > 0 => fastq_bytes,
        _ => 0,
    };
    if preferred > 0 {
        return preferred;
    }
    // If remote sizes are missing, fall back to a conservative estimate from any listed files.
    run.remote_files
        .iter()
        .filter_map(|file: &RemoteFile| file.bytes)
        .sum()
}

fn estimate_download_speed_bps(events: &[&EventRecord]) -> f64 {
    let mut total_bytes = 0_f64;
    let mut total_seconds = 0_f64;
    for run in distinct_event_runs(events, "download", "download started", "download verified") {
        if let (Some(start), Some(done)) = (
            find_event_ts(events, "download", "download started", Some(&run)),
            find_event_ts(events, "download", "download verified", Some(&run)),
        ) && done > start
        {
            let secs = (done - start).num_milliseconds() as f64 / 1000.0;
            if secs > 0.0 {
                let run_bytes = events
                    .iter()
                    .find(|event| {
                        event.stage == "download"
                            && event.message == "download verified"
                            && event.run_accession.as_deref() == Some(run.as_str())
                    })
                    .and_then(|event| event.context.get("bytes"))
                    .and_then(|value| value.as_f64())
                    .unwrap_or(0.0);
                total_seconds += secs;
                total_bytes += run_bytes;
            }
        }
    }

    if total_bytes > 0.0 && total_seconds > 0.0 {
        return total_bytes / total_seconds;
    }

    // Fallback: coarse stage-level estimate from first start to last completion when per-run bytes are unavailable.
    let first_start = earliest_event_ts(events, "download", "download started");
    let last_done = latest_event_ts(events, "download", "download verified");
    let completed = events
        .iter()
        .filter(|event| event.stage == "download" && event.message == "download verified")
        .count();
    if completed == 0 {
        return 0.0;
    }
    let avg_bytes_per_run = events
        .iter()
        .filter(|event| event.stage == "download" && event.message == "download verified")
        .filter_map(|event| event.context.get("bytes").and_then(|value| value.as_f64()))
        .sum::<f64>()
        / (completed as f64);
    match (first_start, last_done) {
        (Some(start), Some(done)) if done > start => {
            let secs = (done - start).num_milliseconds() as f64 / 1000.0;
            if secs > 0.0 {
                (avg_bytes_per_run * completed as f64) / secs
            } else {
                0.0
            }
        }
        _ => 0.0,
    }
}

fn estimate_stage_runs_per_sec(
    events: &[&EventRecord],
    stage: &str,
    started_message: &str,
    completed_message: &str,
) -> f64 {
    let starts = events
        .iter()
        .filter(|event| event.stage == stage && event.message == started_message)
        .count() as f64;
    let completed = events
        .iter()
        .filter(|event| event.stage == stage && event.message == completed_message)
        .count() as f64;
    if starts == 0.0 || completed == 0.0 {
        return 0.0;
    }
    let first_start = earliest_event_ts(events, stage, started_message);
    let last_done = latest_event_ts(events, stage, completed_message);
    match (first_start, last_done) {
        (Some(start), Some(done)) if done > start => {
            let secs = (done - start).num_milliseconds() as f64 / 1000.0;
            if secs > 0.0 { completed / secs } else { 0.0 }
        }
        _ => 0.0,
    }
}

fn estimate_stage_avg_duration_secs(
    events: &[&EventRecord],
    stage: &str,
    started_message: &str,
    completed_message: &str,
) -> Option<f64> {
    let start = latest_event_ts(events, stage, started_message)?;
    let done = latest_event_ts(events, stage, completed_message)?;
    if done <= start {
        return None;
    }
    Some((done - start).num_milliseconds() as f64 / 1000.0)
}

fn duration_from_rate(work_left: f64, rate_per_sec: f64) -> Option<std::time::Duration> {
    if work_left <= 0.0 {
        return Some(std::time::Duration::from_secs(0));
    }
    if rate_per_sec <= 0.0 {
        return None;
    }
    Some(std::time::Duration::from_secs_f64(work_left / rate_per_sec))
}

fn format_duration(duration: std::time::Duration) -> String {
    let secs = duration.as_secs();
    let hours = secs / 3600;
    let minutes = (secs % 3600) / 60;
    let seconds = secs % 60;
    format!("{hours:02}:{minutes:02}:{seconds:02}")
}

fn distinct_event_runs(
    events: &[&EventRecord],
    stage: &str,
    started_message: &str,
    completed_message: &str,
) -> std::collections::BTreeSet<String> {
    events
        .iter()
        .filter(|event| {
            event.stage == stage
                && (event.message == started_message || event.message == completed_message)
        })
        .filter_map(|event| event.run_accession.clone())
        .collect()
}

fn parse_event_ts(ts: &str) -> Option<DateTime<Utc>> {
    DateTime::parse_from_rfc3339(ts)
        .ok()
        .map(|value| value.with_timezone(&Utc))
}

fn find_event_ts(
    events: &[&EventRecord],
    stage: &str,
    message: &str,
    run_accession: Option<&str>,
) -> Option<DateTime<Utc>> {
    events
        .iter()
        .find(|event| {
            event.stage == stage
                && event.message == message
                && event.run_accession.as_deref() == run_accession
        })
        .and_then(|event| parse_event_ts(&event.ts))
}

fn earliest_event_ts(events: &[&EventRecord], stage: &str, message: &str) -> Option<DateTime<Utc>> {
    events
        .iter()
        .filter(|event| event.stage == stage && event.message == message)
        .filter_map(|event| parse_event_ts(&event.ts))
        .min()
}

fn latest_event_ts(events: &[&EventRecord], stage: &str, message: &str) -> Option<DateTime<Utc>> {
    events
        .iter()
        .filter(|event| event.stage == stage && event.message == message)
        .filter_map(|event| parse_event_ts(&event.ts))
        .max()
}

impl ProgressSnapshot {
    pub(crate) fn preprocess_display(&self) -> String {
        if self.preprocess_enabled {
            format!("{}/{}", self.preprocess_done, self.total_runs)
        } else {
            "off".to_string()
        }
    }
}

#[cfg(test)]
mod tests;
