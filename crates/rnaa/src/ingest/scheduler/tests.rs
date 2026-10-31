use super::*;
use rnaa_core::{
    jobs::JobOutcome,
    resources::{ExecutionResources, ReservationOwner, ResourceConfig, ResourceRequest},
};
use std::sync::{
    atomic::{AtomicBool, Ordering},
    mpsc,
};

#[test]
fn scheduler_shutdown_cancels_and_joins_workers_before_returning() -> Result<()> {
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    let mut lease = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    let attempt_id = lease.attempt_id().to_string();
    let context = ExecutionResources::open(
        &ResourceConfig {
            pool_root: root.path().join("pool").display().to_string(),
            cpu_threads: 1,
            memory_mb: 1,
            disk_headroom_mb: 0,
            wait_seconds: 1,
        },
        ReservationOwner {
            project_root: root.path().into(),
            attempt_id: attempt_id.clone(),
            stage: "test".into(),
            run_accession: None,
        },
    )?;
    let pool = context.pool().clone();
    let ended = Arc::new(AtomicBool::new(false));
    let cancelled = Arc::new(AtomicBool::new(false));
    let observed = Arc::clone(&ended);
    let observed_cancel = Arc::clone(&cancelled);
    let (tx, rx) = mpsc::channel();
    let handle = thread::spawn(move || -> Result<()> {
        let owned = lease.database();
        let _memory = context.reserve(
            "engine",
            ResourceRequest {
                cpu_threads: 1,
                memory_bytes: 256,
                disks: Vec::new(),
            },
            &|| owned.check_job(),
        )?;
        tx.send(())?;
        let started = std::time::Instant::now();
        while owned.check_job().is_ok() && started.elapsed() < Duration::from_secs(5) {
            thread::sleep(Duration::from_millis(5));
        }
        observed_cancel.store(owned.check_job().is_err(), Ordering::SeqCst);
        lease.finish(JobOutcome::Interrupted, Some("scheduler shutdown"))?;
        observed.store(true, Ordering::SeqCst);
        Ok(())
    });
    let mut state = QuantSchedulerState::default();
    state.preprocess_tasks.push(RunningTask {
        handle,
        threads_used: 1,
        db: db.clone(),
        attempt_id,
    });
    rx.recv_timeout(Duration::from_secs(5))?;
    assert_eq!(pool.snapshot()?.reservations.len(), 1);
    drop(state);
    assert!(ended.load(Ordering::SeqCst));
    assert!(cancelled.load(Ordering::SeqCst));
    assert!(pool.snapshot()?.reservations.is_empty());
    assert!(db.list_job_attempts(1)?[0].state != "RUNNING");
    Ok(())
}
