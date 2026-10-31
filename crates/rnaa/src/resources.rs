//! Shared execution budget policy for ingest and analysis.
use rnaa_core::ProjectConfig;

#[derive(Debug, Clone, Copy)]
pub(crate) struct SchedulerBudget {
    pub(crate) cpu_budget: usize,
    pub(crate) reserve_workers: usize,
    pub(crate) requested_worker_budget: usize,
    pub(crate) preprocess_workers: usize,
    pub(crate) preprocess_threads: usize,
    pub(crate) quant_workers: usize,
    pub(crate) quant_threads: usize,
    pub(crate) quant_thread_cap: usize,
    pub(crate) deseq2_threads: usize,
    pub(crate) corr_threads: usize,
    pub(crate) download_workers: usize,
}

pub(crate) fn resolve_scheduler_budget(config: &ProjectConfig) -> SchedulerBudget {
    let total_workers = std::thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1);
    let requested_worker_budget = config.quant.worker_budget.min(total_workers);
    let reserve_workers = if requested_worker_budget > 0 {
        0
    } else if total_workers > 4 {
        2
    } else {
        total_workers.saturating_sub(1)
    };
    let cpu_budget = if requested_worker_budget > 0 {
        requested_worker_budget.max(1)
    } else {
        total_workers.saturating_sub(reserve_workers).max(1)
    };
    let preprocess_threads = if config.quant.preprocess {
        config.quant.preprocess_threads.clamp(1, 16).min(cpu_budget)
    } else {
        0
    };
    let quant_thread_cap = config.quant.threads.max(1).min(cpu_budget);
    let quant_threads = quant_thread_cap;
    let quant_cpu_budget = cpu_budget;
    let quant_workers = if config.quant.workers == 0 {
        quant_cpu_budget.saturating_div(quant_thread_cap).max(1)
    } else {
        config
            .quant
            .workers
            .max(1)
            .min(quant_cpu_budget.saturating_div(quant_thread_cap).max(1))
    };
    let preprocess_workers = if config.quant.preprocess {
        cpu_budget
            .saturating_sub(quant_threads)
            .saturating_div(preprocess_threads.max(1))
            .max(1)
            .min(cpu_budget.saturating_div(preprocess_threads.max(1)).max(1))
    } else {
        0
    };
    SchedulerBudget {
        cpu_budget,
        reserve_workers,
        requested_worker_budget,
        preprocess_workers,
        preprocess_threads,
        quant_workers,
        quant_threads,
        quant_thread_cap,
        deseq2_threads: 1,
        corr_threads: if config.corr.threads == 0 {
            cpu_budget
        } else {
            config.corr.threads.max(1).min(cpu_budget)
        },
        download_workers: config.download.concurrency.max(1),
    }
}
