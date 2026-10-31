use super::{
    NativeProcessOptions, batches,
    window::{BatchWindow, PanicGuard, ScopeGuard},
};
use crate::{
    pipeline::{self, Batch, PairedBatch, PairedWorkerResult, WorkerResult},
    split::SplitConfig,
    trimming::TrimmingConfig,
};
use anyhow::{Result, anyhow};
use crossbeam_channel::bounded;
use std::{
    sync::Arc,
    thread::{self, ScopedJoinHandle},
};

fn task<T>(window: Arc<BatchWindow>, work: impl FnOnce() -> Result<T>) -> Result<T> {
    let _guard = PanicGuard(Arc::clone(&window));
    let result = work();
    if result.is_err() {
        window.cancel();
    }
    result
}
fn joined<T>(handle: ScopedJoinHandle<'_, T>, name: &str) -> Result<T> {
    handle
        .join()
        .map_err(|_| anyhow!("native preprocessing {name} thread panicked"))
}

pub(super) fn single(
    input: &str,
    output: &str,
    options: NativeProcessOptions,
    trim: TrimmingConfig,
    window: &Arc<BatchWindow>,
) -> Result<crate::processor::StreamAccumulator> {
    thread::scope(|scope| {
        let _scope_guard = ScopeGuard(Arc::clone(window));
        let (tx, rx) = bounded::<Option<Batch>>(options.max_in_flight_batches);
        let (out_tx, out_rx) = bounded::<Option<WorkerResult>>(options.max_in_flight_batches);
        let producer = thread::Builder::new().spawn_scoped(scope, {
            let window = Arc::clone(window);
            move || {
                task(Arc::clone(&window), || {
                    batches::single(input, tx, options, &window)
                })
            }
        })?;
        let mut workers = Vec::new();
        for _ in 0..options.threads {
            let rx = rx.clone();
            let out_tx = out_tx.clone();
            let trim = trim.clone();
            let window = Arc::clone(window);
            workers.push(thread::Builder::new().spawn_scoped(scope, move || {
                task(window, || {
                    pipeline::worker_thread(rx, out_tx, 15, 5, 15, 40, 0, false, 30, false, trim);
                    Ok(())
                })
            })?);
        }
        drop(rx);
        drop(out_tx);
        let merger = thread::Builder::new().spawn_scoped(scope, {
            let window = Arc::clone(window);
            move || {
                task(window, || {
                    pipeline::merger_thread(
                        out_rx,
                        output.into(),
                        options.threads,
                        Some(6),
                        SplitConfig::default(),
                        false,
                    )
                })
            }
        })?;
        let producer_result = joined(producer, "producer").and_then(|value| value);
        let mut workers_result = Ok(());
        for worker in workers {
            if let Err(error) = joined(worker, "worker").and_then(|value| value) {
                workers_result = Err(error);
            }
        }
        let merged = joined(merger, "merger")??;
        workers_result?;
        producer_result?;
        Ok(merged)
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn paired(
    input1: &str,
    input2: &str,
    output1: &str,
    output2: &str,
    options: NativeProcessOptions,
    trim: TrimmingConfig,
    trim2: TrimmingConfig,
    window: &Arc<BatchWindow>,
) -> Result<crate::processor::PairedEndAccumulator> {
    thread::scope(|scope| {
        let _scope_guard = ScopeGuard(Arc::clone(window));
        let (tx, rx) = bounded::<Option<PairedBatch>>(options.max_in_flight_batches);
        let (out_tx, out_rx) = bounded::<Option<PairedWorkerResult>>(options.max_in_flight_batches);
        let producer = thread::Builder::new().spawn_scoped(scope, {
            let window = Arc::clone(window);
            move || {
                task(Arc::clone(&window), || {
                    batches::paired(input1, input2, tx, options, &window)
                })
            }
        })?;
        let mut workers = Vec::new();
        for _ in 0..options.threads {
            let rx = rx.clone();
            let out_tx = out_tx.clone();
            let trim = trim.clone();
            let trim2 = trim2.clone();
            let window = Arc::clone(window);
            workers.push(thread::Builder::new().spawn_scoped(scope, move || {
                task(window, || {
                    pipeline::paired_worker_thread(
                        rx, out_tx, 15, 5, 15, 40, 0, false, 30, false, trim, trim2, None, None,
                        None, None,
                    );
                    Ok(())
                })
            })?);
        }
        drop(rx);
        drop(out_tx);
        let merger = thread::Builder::new().spawn_scoped(scope, {
            let window = Arc::clone(window);
            move || {
                task(window, || {
                    pipeline::paired_merger_thread(
                        out_rx,
                        output1.into(),
                        output2.into(),
                        options.threads,
                        Some(6),
                        SplitConfig::default(),
                        false,
                    )
                })
            }
        })?;
        let producer_result = joined(producer, "producer").and_then(|value| value);
        let mut workers_result = Ok(());
        for worker in workers {
            if let Err(error) = joined(worker, "worker").and_then(|value| value) {
                workers_result = Err(error);
            }
        }
        let merged = joined(merger, "merger")??;
        workers_result?;
        producer_result?;
        Ok(merged)
    })
}
