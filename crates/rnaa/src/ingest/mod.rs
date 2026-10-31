//! Ingest command orchestration; workers, scheduling and retention stay separate.
use rnaa_core::traits::ReferenceManager;
pub(crate) mod inputs;
pub(crate) mod preprocess;
pub(crate) mod quant;
pub(crate) mod recovery;
pub(crate) mod retention;
mod reuse;
pub(crate) mod scheduler;
mod source;
use crate::analysis::corr::cmd_corr;
use crate::analysis::de::cmd_deseq2;
use crate::cli::{
    CorrArgs, Deseq2Args, DownloadArgs, QuantArgs, ResolveArgs, RunArgs, parse_cleanup_mode,
    parse_cleanup_when, parse_download_preference, parse_fastq_retention,
};
use crate::ingest::inputs::{
    has_pending_download_runs, reconcile_verified_download_state, recover_interrupted_work,
};
use crate::ingest::quant::reconcile_quant_state;
use crate::ingest::scheduler::drive_quant_phase;
use crate::metadata::refresh_samplesheet;
use crate::project::cmd_resolve;
use crate::quality;
use crate::resources::resolve_scheduler_budget;
use crate::run_state::is_quant_done_state;
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::RunState;
use rnaa_core::util::parse_size_bytes;
use rnaa_stages::download::ReqwestDownloader;
use rnaa_stages::refs::EnsemblReferenceManager;
use std::thread;

pub(crate) fn cmd_download(
    paths: &ProjectPaths,
    db: Database,
    mut config: ProjectConfig,
    args: DownloadArgs,
) -> Result<()> {
    recover_interrupted_work(&db)?;
    db.set_project_config(&config)?;
    if let Some(storage_limit) = args.storage_limit {
        parse_size_bytes(&storage_limit)?;
        config.storage.max_active_size = storage_limit;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    let downloader = ReqwestDownloader;
    let prefer = args
        .prefer
        .as_deref()
        .map(parse_download_preference)
        .transpose()?;
    downloader.run_worker(
        db.clone(),
        paths.clone(),
        config.clone(),
        args.concurrency.unwrap_or(config.download.concurrency),
        args.forever,
        prefer,
    )?;
    if !args.forever {
        refresh_samplesheet(paths, &db)?;
    }
    println!("download\tcomplete");
    Ok(())
}

pub(crate) fn cmd_quant(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: QuantArgs,
) -> Result<()> {
    let result = quant_command(paths, db, config, args);
    finish_cleanup(paths, db, result)
}

fn finish_cleanup(paths: &ProjectPaths, db: &Database, result: Result<()>) -> Result<()> {
    let cleanup = ProjectConfig::load(&paths.config_path).and_then(|config| {
        rnaa_stages::cleanup::TrashCleaner.after_command(db, &config, result.is_ok())
    });
    match (result, cleanup) {
        (Err(primary), Err(cleanup)) => {
            tracing::warn!(error=%cleanup, "cleanup after failed command was deferred or held");
            Err(primary.context(format!("input cleanup also reported: {cleanup:#}")))
        }
        (Err(primary), _) => Err(primary),
        (Ok(()), cleanup) => cleanup,
    }
}

fn quant_command(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    args: QuantArgs,
) -> Result<()> {
    recover_interrupted_work(db)?;
    reconcile_verified_download_state(db)?;
    if let Some(strand) = args.strand {
        config.quant.strand = strand;
    }
    if let Some(length) = args.single_fragment_length {
        config.quant.single_fragment_length = length;
    }
    if let Some(sd) = args.single_fragment_sd {
        config.quant.single_fragment_sd = sd;
    }
    config.quant.validate()?;
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    if let Some(engine) = args.engine {
        config.quant.engine = engine;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(worker_budget) = args.worker_budget {
        if worker_budget == 0 {
            bail!("--worker-budget must be >= 1");
        }
        config.quant.worker_budget = worker_budget;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(threads) = args.threads {
        if threads == 0 {
            bail!("--threads must be >= 1");
        }
        config.quant.threads = threads;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(cache_mb) = args.match_cache_mb {
        config.quant.match_cache_mb = cache_mb;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(workers) = args.workers {
        config.quant.workers = workers;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if args.preprocess {
        config.quant.preprocess = true;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(fastq_retention) = args.fastq_retention {
        config.quant.fastq_retention = parse_fastq_retention(&fastq_retention)?;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if args.no_preprocess {
        config.quant.preprocess = false;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if args.preprocess_strict {
        config.quant.preprocess_strict = true;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if args.preprocess_bypass {
        config.quant.preprocess_strict = false;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(preprocess_max_mb) = args.preprocess_max_mb {
        config.quant.preprocess_max_input_mb = preprocess_max_mb;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(preprocess_threads) = args.preprocess_threads {
        if preprocess_threads == 0 {
            bail!("--preprocess-threads must be >= 1");
        }
        config.quant.preprocess_threads = preprocess_threads;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(cleanup) = args.cleanup {
        config.quant.cleanup = parse_cleanup_mode(&cleanup)?;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(when) = args.cleanup_on {
        config.quant.cleanup_on = parse_cleanup_when(&when)?;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(trash_days) = args.trash_days {
        config.quant.trash_days = trash_days;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }

    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, db, &config)?;
    quality::refresh_all(db, &config)?;
    let recovery = reconcile_quant_state(db, paths, &config, &reference)?;
    let restoring_trash = recovery
        .iter()
        .any(|plan| plan.action == recovery::Action::RestoreTrash);
    let mut recovery_failed = recovery::restage(&recovery, db, paths, &config)?;
    // Restoration changes physical availability while preserving scientific input
    // identity. Refresh its admission before selecting the next processing action.
    let recovery = if restoring_trash {
        quality::refresh_all(db, &config)?;
        reconcile_quant_state(db, paths, &config, &reference)?
    } else {
        recovery
    };
    let blocked_runs = recovery
        .iter()
        .filter(|plan| plan.action == recovery::Action::Blocked)
        .map(|plan| plan.run.clone())
        .collect();
    let scheduler = resolve_scheduler_budget(&config);
    let summary = drive_quant_phase(
        paths,
        db,
        &config,
        &reference,
        &scheduler,
        &blocked_runs,
        || Ok((false, false)),
    )?;

    for run in db.list_active_runs()? {
        if !is_quant_done_state(run.state)
            && db
                .generation_head(
                    rnaa_core::generations::OutputGroup::Quant,
                    &run.run_accession,
                )?
                .is_some()
        {
            recovery_failed.insert(run.run_accession);
        }
    }
    let recovery_failed = recovery_failed.len();
    if !summary.did_work {
        println!("quant\tno_ready_runs");
        let blocked = quality::blocked_count(db, &config)?;
        if blocked > 0 || recovery_failed > 0 {
            bail!(
                "{blocked} requested run(s) are ineligible under current QC policy; {recovery_failed} dependency recoveries are blocked or failed; prior counts are preserved"
            );
        }
        return Ok(());
    }

    let tracked = summary.attempted_runs;
    let runs = db.list_active_runs()?;
    let completed = runs
        .iter()
        .filter(|run| tracked.contains(&run.run_accession) && is_quant_done_state(run.state))
        .count() as u64;
    let failed = runs
        .iter()
        .filter(|run| {
            tracked.contains(&run.run_accession)
                && matches!(
                    run.state,
                    RunState::PreprocessFailed | RunState::QuantFailed
                )
        })
        .count() as u64;
    println!("quant_workers\t{}", scheduler.quant_workers);
    println!("quant_completed\t{}", completed);
    println!("quant_failed\t{}", failed);
    let blocked = quality::blocked_count(db, &config)?;
    println!("quant_ineligible\t{blocked}");
    if failed > 0 || blocked > 0 || recovery_failed > 0 {
        bail!(
            "quantification failed for {failed} requested run(s); {blocked} run(s) are ineligible under current QC policy; {recovery_failed} dependency recoveries are blocked or failed"
        );
    }
    Ok(())
}

pub(crate) fn cmd_run(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: RunArgs,
) -> Result<()> {
    let result = run_command(paths, db, config, args);
    finish_cleanup(paths, db, result)
}

fn run_command(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    args: RunArgs,
) -> Result<()> {
    if let Some(cache_mb) = args.match_cache_mb {
        config.quant.match_cache_mb = cache_mb;
    }
    if let Some(cleanup) = args.cleanup {
        config.quant.cleanup = parse_cleanup_mode(&cleanup)?;
    }
    if let Some(when) = args.cleanup_on {
        config.quant.cleanup_on = parse_cleanup_when(&when)?;
    }
    if let Some(days) = args.trash_days {
        config.quant.trash_days = days;
    }
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    if let Some(fastq_retention) = args.fastq_retention {
        config.quant.fastq_retention = parse_fastq_retention(&fastq_retention)?;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(worker_budget) = args.worker_budget {
        if worker_budget == 0 {
            bail!("--worker-budget must be >= 1");
        }
        config.quant.worker_budget = worker_budget;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(preprocess_threads) = args.preprocess_threads {
        if preprocess_threads == 0 {
            bail!("--preprocess-threads must be >= 1");
        }
        config.quant.preprocess_threads = preprocess_threads;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    recover_interrupted_work(db)?;
    reconcile_verified_download_state(db)?;
    cmd_resolve(
        paths,
        db,
        &config,
        ResolveArgs {
            force: false,
            metadata_conflicts: "fail".into(),
        },
    )?;
    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, db, &config)?;
    quality::refresh_all(db, &config)?;
    let recovery = reconcile_quant_state(db, paths, &config, &reference)?;
    let restoring_trash = recovery
        .iter()
        .any(|plan| plan.action == recovery::Action::RestoreTrash);
    for plan in recovery
        .iter()
        .filter(|plan| plan.action == recovery::Action::RestoreTrash)
    {
        recovery::restore_trash(db, paths, &config, plan)?;
    }
    let recovery = if restoring_trash {
        quality::refresh_all(db, &config)?;
        reconcile_quant_state(db, paths, &config, &reference)?
    } else {
        recovery
    };
    let blocked_runs = recovery
        .iter()
        .filter(|plan| plan.action == recovery::Action::Blocked)
        .map(|plan| plan.run.clone())
        .collect();
    let scheduler = resolve_scheduler_budget(&config);

    let mut download_handle = if args.no_download {
        let missing = recovery
            .iter()
            .filter(|plan| plan.action == recovery::Action::Download)
            .count();
        if missing > 0 {
            bail!(
                "{missing} stale count dependencies require original reads, but --no-download disables input recovery; restore the reads or allow downloads"
            );
        }
        None
    } else {
        let downloader = ReqwestDownloader;
        let db_for_download = db.clone();
        let paths_for_download = paths.clone();
        let config_for_download = config.clone();
        Some(thread::spawn(move || {
            downloader.run_worker(
                db_for_download,
                paths_for_download,
                config_for_download.clone(),
                scheduler.download_workers,
                false,
                None,
            )
        }))
    };
    drive_quant_phase(
        paths,
        db,
        &config,
        &reference,
        &scheduler,
        &blocked_runs,
        || {
            let mut did_work = false;
            if let Some(handle) = &download_handle
                && handle.is_finished()
            {
                let join_result = download_handle.take().expect("download handle present");
                let worker_result = join_result
                    .join()
                    .map_err(|_| anyhow::anyhow!("download worker thread panicked"))?;
                worker_result?;
                refresh_samplesheet(paths, db)?;
                did_work = true;
            }
            let pending_download = if args.no_download {
                false
            } else {
                has_pending_download_runs(db)?
            };
            Ok((did_work, pending_download))
        },
    )?;

    let failed = db
        .list_active_runs_in_states(&[
            RunState::DownloadFailed,
            RunState::VerifyFailed,
            RunState::PreprocessFailed,
            RunState::QuantFailed,
        ])?
        .len();
    let blocked = quality::blocked_count(db, &config)?;
    if failed > 0 || blocked > 0 {
        bail!(
            "pipeline failed for {failed} requested run(s); {blocked} run(s) are ineligible under current QC policy"
        );
    }

    cmd_deseq2(
        paths,
        db,
        config.clone(),
        Deseq2Args {
            design: None,
            contrast: Vec::new(),
        },
    )?;
    if !args.no_corr {
        cmd_corr(
            paths,
            db,
            config.clone(),
            CorrArgs {
                matrix: None,
                adjust: None,
                model: None,
                method: None,
                geneset: None,
                out_spec: None,
                module_gmt: None,
                module_top: 500,
                module_min_size: 25,
                fgsea_permutations: 1000,
                fgsea_seed: 42,
                fgsea_min_size: 10,
                fgsea_max_size: None,
            },
        )?;
    }
    println!("run\tcompleted");
    Ok(())
}

pub(crate) fn finish_boolean_job(lease: &mut JobLease, result: Result<bool>) -> Result<bool> {
    match result {
        Ok(success) => {
            lease.finish(
                if success {
                    JobOutcome::Succeeded
                } else {
                    JobOutcome::Failed
                },
                None,
            )?;
            Ok(success)
        }
        Err(error) => {
            if let Err(finish_error) = lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))
            {
                return Err(
                    error.context(format!("could not finalize failed job: {finish_error:#}"))
                );
            }
            Err(error)
        }
    }
}

#[cfg(test)]
mod finish_tests {
    use super::*;

    #[test]
    fn finalization_failure_preserves_the_original_error() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let db = Database::create(temp.path(), &ProjectConfig::default())?;
        let mut lease = JobLease::try_acquire(&db, "export", None, &[])?.expect("lease");
        lease
            .database()
            .complete_owned_job(JobOutcome::Succeeded, None)?;
        let result = finish_boolean_job(&mut lease, Err(anyhow::anyhow!("original failure")));
        let message = format!("{:#}", result.unwrap_err());
        assert!(message.contains("original failure"), "{message}");
        assert!(
            message.contains("could not finalize failed job"),
            "{message}"
        );
        Ok(())
    }
}
