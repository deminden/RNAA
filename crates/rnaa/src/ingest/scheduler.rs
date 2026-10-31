//! Worker scheduling and joins; callers receive an outcome, not mutable worker state.
use crate::ingest::preprocess::process_preprocess_run;
use crate::ingest::quant::process_quant_run;
use crate::progress::print_run_progress;
use crate::quality;
use crate::resources::SchedulerBudget;
use crate::run_state::is_quant_done_state;
use anyhow::Result;
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::jobs::JobLease;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::RunState;
use rnaa_stages::quant::QuantIndexGroup;
use std::collections::HashSet;
use std::sync::Arc;
use std::thread;
use std::time::Duration;
#[cfg(test)]
mod tests;

struct RunningTask {
    handle: thread::JoinHandle<Result<()>>,
    threads_used: usize,
    db: Database,
    attempt_id: String,
}

#[derive(Default)]
struct QuantSchedulerState {
    index_group: Option<Arc<QuantIndexGroup>>,
    preprocess_tasks: Vec<RunningTask>,
    quant_tasks: Vec<RunningTask>,
    attempted_preprocess: HashSet<String>,
    attempted_quant: HashSet<String>,
    cleanup_done: HashSet<(String, String)>,
}

impl Drop for QuantSchedulerState {
    fn drop(&mut self) {
        // An error in the driver must not detach remaining workers. Request
        // cancellation first for every sibling, then retain their ownership and
        // resource lifetimes until all actual work has stopped.
        for task in self.preprocess_tasks.iter().chain(&self.quant_tasks) {
            if let Err(error) = task.db.cancel_job(&task.attempt_id) {
                tracing::warn!(attempt=%task.attempt_id,%error,"cannot request worker cancellation during scheduler shutdown");
            }
        }
        for task in self
            .preprocess_tasks
            .drain(..)
            .chain(self.quant_tasks.drain(..))
        {
            match task.handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(error)) => {
                    tracing::warn!(attempt=%task.attempt_id,%error,"worker ended during scheduler shutdown")
                }
                Err(_) => {
                    tracing::warn!(attempt=%task.attempt_id,"worker panicked during scheduler shutdown")
                }
            }
        }
    }
}

pub(crate) struct QuantPhaseSummary {
    pub(crate) did_work: bool,
    pub(crate) attempted_runs: HashSet<String>,
}

pub(crate) fn drive_quant_phase<F>(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    reference: &rnaa_core::ReferenceBundle,
    scheduler: &SchedulerBudget,
    blocked_runs: &HashSet<String>,
    mut external_step: F,
) -> Result<QuantPhaseSummary>
where
    F: FnMut() -> Result<(bool, bool)>,
{
    let mut scheduler_state = QuantSchedulerState {
        index_group: None,
        preprocess_tasks: Vec::new(),
        quant_tasks: Vec::new(),
        attempted_preprocess: blocked_runs.clone(),
        attempted_quant: blocked_runs.clone(),
        cleanup_done: HashSet::new(),
    };
    let started_at = std::time::Instant::now();
    let mut last_progress = started_at;
    let mut last_progress_snapshot = None;
    print_run_progress(
        paths,
        db,
        config,
        scheduler,
        started_at.elapsed(),
        &mut last_progress_snapshot,
    )?;

    let mut did_work = false;
    loop {
        if last_progress.elapsed() >= Duration::from_secs(10) {
            print_run_progress(
                paths,
                db,
                config,
                scheduler,
                started_at.elapsed(),
                &mut last_progress_snapshot,
            )?;
            last_progress = std::time::Instant::now();
        }
        let scheduled = schedule_quant_iteration(
            paths,
            db,
            config,
            reference,
            scheduler,
            &mut scheduler_state,
        )?;
        did_work |= scheduled;
        let (external_did_work, pending_external_work) = external_step()?;
        did_work |= external_did_work;
        if external_did_work {
            last_progress = std::time::Instant::now();
            print_run_progress(
                paths,
                db,
                config,
                scheduler,
                started_at.elapsed(),
                &mut last_progress_snapshot,
            )?;
        }
        if quant_scheduler_is_idle(db, config, &scheduler_state)? && !pending_external_work {
            break;
        }
        if !scheduled && !external_did_work {
            thread::sleep(Duration::from_secs(2));
        }
    }

    print_run_progress(
        paths,
        db,
        config,
        scheduler,
        started_at.elapsed(),
        &mut last_progress_snapshot,
    )?;

    Ok(QuantPhaseSummary {
        did_work,
        attempted_runs: scheduler_state
            .attempted_preprocess
            .union(&scheduler_state.attempted_quant)
            .cloned()
            .collect(),
    })
}

fn schedule_quant_iteration(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    reference: &rnaa_core::ReferenceBundle,
    scheduler: &SchedulerBudget,
    state: &mut QuantSchedulerState,
) -> Result<bool> {
    reap_finished_tasks(&mut state.preprocess_tasks)?;
    reap_finished_tasks(&mut state.quant_tasks)?;
    if state.preprocess_tasks.is_empty() && state.quant_tasks.is_empty() {
        db.recover_jobs()?;
    }

    let runs = db.list_active_runs()?;
    if state.preprocess_tasks.is_empty()
        && state.quant_tasks.is_empty()
        && config.quant.cleanup != rnaa_core::CleanupMode::None
        && config.quant.cleanup_on != rnaa_core::CleanupWhen::Never
    {
        for run in &runs {
            if !is_quant_done_state(run.state) || !quality::allowed(db, config, &run.run_accession)?
            {
                continue;
            }
            let Some(head) = db.generation_head(
                rnaa_core::generations::OutputGroup::Quant,
                &run.run_accession,
            )?
            else {
                continue;
            };
            let key = (run.run_accession.clone(), head.generation_id);
            if state.cleanup_done.contains(&key) {
                continue;
            }
            match rnaa_stages::cleanup::TrashCleaner.try_completed_run(
                db,
                config,
                &run.run_accession,
            ) {
                Ok(true) => {
                    state.cleanup_done.insert(key);
                }
                Ok(false) => {}
                Err(error) => {
                    tracing::warn!(run=%run.run_accession, error=%error, "completed count inputs held during cleanup");
                    state.cleanup_done.insert(key);
                }
            }
        }
    }
    let mut preprocess_running = runs
        .iter()
        .filter(|run| matches!(run.state, RunState::PreprocessRunning))
        .count();
    let mut quant_running = runs
        .iter()
        .filter(|run| matches!(run.state, RunState::QuantRunning))
        .count();
    let preprocess_cpu_in_use = state
        .preprocess_tasks
        .iter()
        .map(|task| task.threads_used)
        .sum::<usize>();
    let quant_cpu_in_use = state
        .quant_tasks
        .iter()
        .map(|task| task.threads_used)
        .sum::<usize>();
    let mut available_cpu = scheduler
        .cpu_budget
        .saturating_sub(preprocess_cpu_in_use + quant_cpu_in_use);
    let mut did_work = false;

    if config.quant.preprocess {
        let mut preprocess_ready = Vec::new();
        for run in &runs {
            if state.attempted_preprocess.contains(&run.run_accession) {
                continue;
            }
            if matches!(run.state, RunState::Verified | RunState::PreprocessFailed)
                || (matches!(run.state, RunState::PreprocessDone) || is_quant_done_state(run.state))
                    && quality::needs_repair(db, config, run)?
            {
                preprocess_ready.push(run.clone());
            }
        }
        while available_cpu >= scheduler.preprocess_threads.max(1)
            && preprocess_running < scheduler.preprocess_workers
            && !preprocess_ready.is_empty()
        {
            let run = preprocess_ready.remove(0);
            let Some(lease) = JobLease::try_acquire(
                db,
                "preprocess",
                Some(&run.run_accession),
                &[
                    RunState::Verified,
                    RunState::PreprocessFailed,
                    RunState::PreprocessDone,
                    RunState::QuantDone,
                    RunState::DeDone,
                    RunState::DeFailed,
                    RunState::CorrDone,
                    RunState::CorrFailed,
                    RunState::Cleaned,
                ],
            )?
            else {
                continue;
            };
            state.attempted_preprocess.insert(run.run_accession.clone());
            let paths = paths.clone();
            let mut config = config.clone();
            let threads_used = scheduler.preprocess_threads.max(1);
            config.quant.preprocess_threads = threads_used;
            state.preprocess_tasks.push(RunningTask {
                db: db.clone(),
                attempt_id: lease.attempt_id().to_string(),
                handle: thread::spawn(move || {
                    process_preprocess_run(lease, &paths, &config, run).map(|_| ())
                }),
                threads_used,
            });
            preprocess_running += 1;
            available_cpu = available_cpu.saturating_sub(threads_used);
            did_work = true;
        }
    }

    let quant_states = if config.quant.preprocess {
        vec![RunState::PreprocessDone, RunState::QuantFailed]
    } else {
        vec![RunState::Verified, RunState::QuantFailed]
    };
    let quant_ready = db.list_active_runs_in_states(&quant_states)?;
    let mut allowed_quant = Vec::new();
    for run in quant_ready {
        if !state.attempted_quant.contains(&run.run_accession)
            && quality::allowed(db, config, &run.run_accession)?
        {
            allowed_quant.push(run);
        }
    }
    let mut quant_ready = allowed_quant;
    if quant_ready.is_empty()
        && state.quant_tasks.is_empty()
        && let Some(group) = &state.index_group
    {
        // No ready consumer needs an idle graph. Release it while upstream
        // producers supply another wave instead of pinning their memory budget.
        group.release_idle()?;
    }
    if !quant_ready.is_empty() && state.index_group.is_none() {
        let selected = crate::reference::pin(
            db,
            reference,
            config,
            rnaa_core::references::ReferenceScope::Index,
        )?;
        state.index_group = Some(Arc::new(QuantIndexGroup::new(selected.index)?));
    }
    while available_cpu >= scheduler.quant_threads.min(scheduler.quant_thread_cap)
        && quant_running < scheduler.quant_workers
        && !quant_ready.is_empty()
    {
        let run = quant_ready.remove(0);
        let Some(lease) =
            JobLease::try_acquire(db, "quant", Some(&run.run_accession), &quant_states)?
        else {
            continue;
        };
        state.attempted_quant.insert(run.run_accession.clone());
        let paths = paths.clone();
        let mut config = config.clone();
        let reference = reference.clone();
        let index_group = Arc::clone(state.index_group.as_ref().expect("initialized index group"));
        let threads_used = available_cpu.min(scheduler.quant_thread_cap).max(1);
        config.quant.threads = threads_used;
        state.quant_tasks.push(RunningTask {
            db: db.clone(),
            attempt_id: lease.attempt_id().to_string(),
            handle: thread::spawn(move || {
                process_quant_run(lease, &paths, &config, &reference, &index_group, run).map(|_| ())
            }),
            threads_used,
        });
        quant_running += 1;
        available_cpu = available_cpu.saturating_sub(threads_used);
        did_work = true;
    }

    Ok(did_work)
}

fn quant_scheduler_is_idle(
    db: &Database,
    config: &ProjectConfig,
    state: &QuantSchedulerState,
) -> Result<bool> {
    if !state.preprocess_tasks.is_empty() || !state.quant_tasks.is_empty() {
        return Ok(false);
    }
    for run in db.list_active_runs()? {
        if run.state == RunState::Downloading
            && db
                .generation_head(
                    rnaa_core::generations::OutputGroup::Quant,
                    &run.run_accession,
                )?
                .is_some()
        {
            return Ok(false);
        }
        if matches!(
            run.state,
            RunState::PreprocessRunning | RunState::QuantRunning
        ) {
            return Ok(false);
        }
        if state.attempted_preprocess.contains(&run.run_accession)
            && state.attempted_quant.contains(&run.run_accession)
        {
            continue;
        }
        if config.quant.preprocess {
            if !state.attempted_preprocess.contains(&run.run_accession)
                && (run.state == RunState::Verified
                    || (run.state == RunState::PreprocessDone || is_quant_done_state(run.state))
                        && quality::needs_repair(db, config, &run)?)
            {
                return Ok(false);
            }
            if run.state == RunState::PreprocessDone
                && !state.attempted_quant.contains(&run.run_accession)
                && quality::allowed(db, config, &run.run_accession)?
            {
                return Ok(false);
            }
        } else if run.state == RunState::Verified
            && !state.attempted_quant.contains(&run.run_accession)
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn reap_finished_tasks(tasks: &mut Vec<RunningTask>) -> Result<()> {
    let mut index = 0;
    while index < tasks.len() {
        if tasks[index].handle.is_finished() {
            let task = tasks.swap_remove(index);
            task.handle
                .join()
                .map_err(|_| anyhow::anyhow!("worker thread panicked"))??;
        } else {
            index += 1;
        }
    }
    Ok(())
}
