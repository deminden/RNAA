use super::*;
use crate::{ArtifactKind, ArtifactRecord, LibraryLayout, ProjectConfig, RunRecord};
use std::collections::BTreeMap;
use std::process::{Child, Command, Stdio};
use std::sync::{Arc, Barrier};
use std::time::Instant;

fn setup() -> Result<(tempfile::TempDir, Database)> {
    let root = tempfile::tempdir()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    db.upsert_run(&RunRecord {
        project_id: db.project_id().into(),
        study_accession: "study".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: vec![],
        metadata: BTreeMap::new(),
        state: RunState::Verified,
        last_error: None,
        updated_at: now_rfc3339(),
    })?;
    Ok((root, db))
}

#[test]
fn retention_owns_the_run_without_changing_completed_progress() -> Result<()> {
    let (_root, db) = setup()?;
    db.set_run_state("RUN1", RunState::QuantDone, None)?;
    let mut lease =
        JobLease::try_acquire(&db, "retention", Some("RUN1"), &[RunState::QuantDone])?.unwrap();
    let owned = lease.database();
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantDone);
    assert!(JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::QuantDone])?.is_none());
    lease.finish(JobOutcome::Failed, Some("cleanup failed"))?;
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantDone);
    assert!(
        owned
            .append_event("retention", Some("RUN1"), "stale", serde_json::json!({}))
            .is_err()
    );
    let base = lease.finish_and_release(JobOutcome::Failed, None)?;
    let next =
        JobLease::try_acquire(&base, "retention", Some("RUN1"), &[RunState::QuantDone])?.unwrap();
    assert!(db.cancel_job(next.attempt_id())?);
    assert!(next.database().check_job().is_err());
    drop(next);
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantDone);
    assert_eq!(db.list_job_attempts(1)?[0].state, "CANCELLED");
    Ok(())
}

#[test]
fn competing_claims_execute_once_and_stale_queue_entries_do_not_restart() -> Result<()> {
    let (_root, db) = setup()?;
    let barrier = Arc::new(Barrier::new(8));
    let mut handles = Vec::new();
    for _ in 0..8 {
        let db = db.clone();
        let barrier = Arc::clone(&barrier);
        handles.push(thread::spawn(move || -> Result<bool> {
            barrier.wait();
            let Some(mut job) =
                JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::Verified])?
            else {
                return Ok(false);
            };
            job.database()
                .set_run_state("RUN1", RunState::QuantDone, None)?;
            job.finish(JobOutcome::Succeeded, None)?;
            Ok(true)
        }));
    }
    let mut claims = 0;
    for handle in handles {
        claims += usize::from(handle.join().unwrap()?);
    }
    assert_eq!(claims, 1);
    assert_eq!(db.list_job_attempts(100)?.len(), 1);
    assert!(JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::Verified])?.is_none());
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantDone);
    Ok(())
}

#[test]
fn expired_live_writer_is_not_recovered_and_stale_writes_are_fenced() -> Result<()> {
    let (_root, db) = setup()?;
    let lease = JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::Verified])?.unwrap();
    let stale = lease.database();
    stale.set_job_request(&"a".repeat(64))?;
    assert!(stale.set_job_request(&"b".repeat(64)).is_err());
    db.write(|connection| {
        connection.execute("UPDATE job_attempts SET lease_expires_ms = 0", [])?;
        Ok(())
    })?;
    assert_eq!(db.recover_jobs()?, 0);
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantRunning);
    assert!(
        stale
            .set_run_state("RUN1", RunState::QuantDone, None)
            .is_err()
    );
    assert!(
        JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::QuantRunning])?.is_none()
    );
    drop(lease);
    let mut next =
        JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::QuantFailed])?.unwrap();
    next.database()
        .set_run_state("RUN1", RunState::QuantDone, None)?;
    next.finish(JobOutcome::Succeeded, None)?;
    assert!(
        stale
            .set_run_state("RUN1", RunState::QuantFailed, Some("stale worker"))
            .is_err()
    );
    assert!(
        stale
            .append_event(
                "quant",
                Some("RUN1"),
                "stale completion",
                serde_json::json!({})
            )
            .is_err()
    );
    assert!(
        stale
            .record_artifact(&ArtifactRecord {
                project_id: db.project_id().into(),
                run_accession: Some("RUN1".into()),
                kind: ArtifactKind::Manifest,
                path: "/stale-manifest".into(),
                blob_id: None,
                shared_path: None,
                checksum_type: "sha256".into(),
                checksum: "x".into(),
                bytes: 1,
                created_at: now_rfc3339(),
            })
            .is_err()
    );
    assert!(db.list_artifacts()?.is_empty());
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantDone);
    let attempts = db.list_job_attempts(10)?;
    assert_eq!(attempts[0].fencing_token, 2);
    assert_eq!(attempts[0].state, "SUCCEEDED");
    assert_eq!(attempts[1].state, "INTERRUPTED");
    Ok(())
}

#[test]
fn heartbeat_renews_and_cancellation_rejects_completion() -> Result<()> {
    let (_root, db) = setup()?;
    let mut lease = JobLease::acquire_with_timing(
        &db,
        "quant",
        Some("RUN1"),
        &[RunState::Verified],
        Duration::from_secs(2),
        Duration::from_millis(20),
    )?
    .unwrap();
    let first = db.list_job_attempts(1)?.remove(0);
    let deadline = Instant::now() + Duration::from_secs(3);
    loop {
        let renewed = db.list_job_attempts(1)?.remove(0);
        if renewed.heartbeat_ms > first.heartbeat_ms {
            assert!(renewed.lease_expires_ms > first.lease_expires_ms);
            break;
        }
        assert!(Instant::now() < deadline, "heartbeat did not renew");
        thread::sleep(Duration::from_millis(10));
    }
    assert!(db.cancel_job(lease.attempt_id())?);
    assert!(
        lease
            .database()
            .set_run_state("RUN1", RunState::QuantDone, None)
            .is_err()
    );
    assert!(lease.finish(JobOutcome::Succeeded, None).is_err());
    drop(lease);
    let job = db.list_job_attempts(1)?.remove(0);
    assert_eq!(job.state, "CANCELLED");
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantFailed);
    assert!(!db.cancel_job(&job.attempt_id)?);
    Ok(())
}

#[test]
fn live_analysis_is_not_recovered_and_analysis_stages_exclude_each_other() -> Result<()> {
    let (_root, db) = setup()?;
    db.set_run_state("RUN1", RunState::QuantDone, None)?;
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    assert!(
        owned
            .set_run_state("RUN1", RunState::DeRunning, None)
            .is_err()
    );
    owned.set_project_stage_state("deseq2", crate::ProjectStageState::Running, None)?;
    assert_eq!(db.recover_jobs()?, 0);
    assert!(JobLease::try_acquire(&db, "corr", None, &[])?.is_none());
    assert!(JobLease::try_acquire(&db, "resolve", None, &[])?.is_none());
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantDone);
    assert!(owned.set_run_state("RUN1", RunState::DeDone, None).is_err());
    owned.set_project_stage_state("deseq2", crate::ProjectStageState::Done, None)?;
    lease.finish(JobOutcome::Succeeded, None)?;
    drop(lease);
    assert!(JobLease::try_acquire(&db, "corr", None, &[])?.is_some());
    Ok(())
}

struct KillOnDrop(Child);
impl Drop for KillOnDrop {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

#[test]
fn crash_worker_process() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_TEST_CRASH_JOB_ROOT") else {
        return Ok(());
    };
    let root = std::path::Path::new(&root);
    let db = Database::open(root)?;
    let lease = JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::Verified])?.unwrap();
    std::fs::write(root.join("worker-ready"), lease.attempt_id())?;
    loop {
        thread::park_timeout(Duration::from_secs(1));
    }
}

#[test]
fn killed_process_releases_writer_lock_and_expired_attempt_can_retry() -> Result<()> {
    let (root, db) = setup()?;
    let mut child = KillOnDrop(
        Command::new(std::env::current_exe()?)
            .args([
                "--exact",
                "jobs::tests::crash_worker_process",
                "--nocapture",
            ])
            .env("RNAA_TEST_CRASH_JOB_ROOT", root.path())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()?,
    );
    let deadline = Instant::now() + Duration::from_secs(10);
    while !root.path().join("worker-ready").exists() {
        assert!(
            child.0.try_wait()?.is_none(),
            "worker exited before claiming its job"
        );
        assert!(Instant::now() < deadline, "worker did not start");
        thread::sleep(Duration::from_millis(10));
    }
    assert_eq!(db.recover_jobs()?, 0);
    assert!(
        JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::QuantRunning])?.is_none()
    );
    child.0.kill()?;
    child.0.wait()?;
    assert_eq!(
        db.recover_jobs()?,
        0,
        "unexpired leases are not presumed abandoned"
    );
    db.write(|connection| {
        connection.execute("UPDATE job_attempts SET lease_expires_ms = 0", [])?;
        Ok(())
    })?;
    assert_eq!(db.recover_jobs()?, 1);
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantFailed);
    let mut next =
        JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::QuantFailed])?.unwrap();
    next.database()
        .set_run_state("RUN1", RunState::QuantDone, None)?;
    next.finish(JobOutcome::Succeeded, None)?;
    let jobs = db.list_job_attempts(10)?;
    assert_eq!(jobs[0].fencing_token, 2);
    assert_eq!(jobs[1].state, "INTERRUPTED");
    Ok(())
}
#[test]
fn reference_jobs_have_independent_ownership_inside_an_analysis_and_fence_stale_children()
-> Result<()> {
    let (_root, db) = setup()?;
    let mut parent = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    let owned = parent.database();
    let mut reference = JobLease::try_acquire_reference(&owned, "custom_reference")?.unwrap();
    let stale = reference.database();
    assert!(JobLease::try_acquire_reference(&owned, "custom_reference")?.is_none());
    assert!(JobLease::try_acquire(&db, "corr", None, &[])?.is_none());
    assert!(
        owned
            .atomic(|db| JobLease::try_acquire_reference(db, "another"))
            .is_err()
    );
    assert!(JobLease::try_acquire_reference(&owned, "../escape").is_err());
    reference.finish(JobOutcome::Succeeded, None)?;
    drop(reference);
    assert!(
        stale
            .append_event("refs", None, "late", serde_json::json!({}))
            .is_err()
    );
    owned.check_job()?;
    assert!(db.cancel_job(parent.attempt_id())?);
    assert!(JobLease::try_acquire_reference(&owned, "another").is_err());
    parent.finish(JobOutcome::Failed, Some("cancelled"))?;
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::Verified);
    Ok(())
}
