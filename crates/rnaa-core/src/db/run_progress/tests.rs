use super::*;
use crate::{
    LibraryLayout, ProjectConfig, RunRecord,
    jobs::{JobLease, JobOutcome},
};
use std::{collections::BTreeMap, fs};

fn legacy(root: &Path) -> Result<Database> {
    let path = root.join("state.sqlite");
    let connection = Database::connect_path(&path)?;
    for (version, sql) in MIGRATIONS
        .iter()
        .filter(|(version, _)| *version != "0021_run_progress")
    {
        connection.execute_batch(sql)?;
        connection.execute(
            "INSERT INTO schema_migrations(version,applied_at) VALUES (?1,?2)",
            params![version, now_rfc3339()],
        )?;
    }
    connection.execute("INSERT INTO projects(project_id,root_dir,created_at,config_json) VALUES('legacy',?1,?2,?3)", params![root.display().to_string(),now_rfc3339(),serde_json::to_string(&ProjectConfig::default())?])?;
    Ok(Database {
        path,
        project_id: "legacy".into(),
        job: None,
        transaction: None,
        generation_maintenance: false,
    })
}
fn run(db: &Database, id: &str, state: RunState) -> RunRecord {
    RunRecord {
        project_id: db.project_id().into(),
        study_accession: "study".into(),
        source_study_accession: None,
        run_accession: id.into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: vec![],
        metadata: BTreeMap::new(),
        state,
        last_error: Some(format!("original {state} note")),
        updated_at: "2026-10-29T10:00:00Z".into(),
    }
}

#[test]
fn legacy_analysis_labels_are_projected_with_immutable_lossless_receipts() -> Result<()> {
    let root = tempfile::tempdir()?;
    let db = legacy(root.path())?;
    for (i, state) in [
        RunState::DeRunning,
        RunState::DeDone,
        RunState::DeFailed,
        RunState::CorrRunning,
        RunState::CorrDone,
        RunState::CorrFailed,
    ]
    .into_iter()
    .enumerate()
    {
        db.upsert_run(&run(&db, &format!("RUN{i}"), state))?;
    }
    let ingest = run(&db, "DOWNLOAD", RunState::Downloading);
    db.upsert_run(&ingest)?;
    let upgraded = Database::open(root.path())?;
    let receipts = upgraded.legacy_run_analysis_states()?;
    assert_eq!(receipts.len(), 6);
    for receipt in &receipts {
        assert_eq!(receipt.source_updated_at, "2026-10-29T10:00:00Z");
        assert_eq!(
            receipt.analysis_error,
            Some(format!("original {} note", receipt.analysis_state))
        );
        let current = upgraded.get_run(&receipt.run_accession)?.unwrap();
        assert_eq!(current.state, RunState::QuantDone);
        assert!(current.last_error.is_none());
    }
    assert_eq!(
        serde_json::to_value(upgraded.get_run("DOWNLOAD")?.unwrap())?,
        serde_json::to_value(ingest)?
    );
    assert!(upgraded.generation_heads()?.is_empty());
    assert!(upgraded.list_artifacts()?.is_empty());
    assert!(upgraded.list_project_stage_states()?.is_empty());
    let connection = upgraded.connect()?;
    assert!(
        connection
            .execute(
                "UPDATE legacy_run_analysis_states SET analysis_error='edited'",
                []
            )
            .is_err()
    );
    assert!(
        connection
            .execute("DELETE FROM legacy_run_analysis_states", [])
            .is_err()
    );
    assert_eq!(
        Database::open(root.path())?
            .legacy_run_analysis_states()?
            .len(),
        6
    );
    Ok(())
}
#[test]
fn failed_projection_rolls_back_receipts_rows_guards_and_version() -> Result<()> {
    let root = tempfile::tempdir()?;
    let db = legacy(root.path())?;
    let original = run(&db, "RUN1", RunState::CorrFailed);
    db.upsert_run(&original)?;
    db.connect()?.execute_batch("CREATE TRIGGER reject_projection BEFORE UPDATE ON runs BEGIN SELECT RAISE(ABORT,'injected projection failure'); END;")?;
    assert!(Database::open(root.path()).is_err());
    assert_eq!(
        serde_json::to_value(db.get_run("RUN1")?.unwrap())?,
        serde_json::to_value(original)?
    );
    let connection = db.connect()?;
    assert_eq!(connection.query_row("SELECT count(*) FROM sqlite_master WHERE name IN ('legacy_run_analysis_states','run_progress_reject_analysis_update')",[],|row|row.get::<_,i64>(0))?,0);
    assert_eq!(
        connection.query_row(
            "SELECT count(*) FROM schema_migrations WHERE version='0021_run_progress'",
            [],
            |row| row.get::<_, i64>(0)
        )?,
        0
    );
    connection.execute_batch("DROP TRIGGER reject_projection")?;
    assert_eq!(
        Database::open(root.path())?
            .legacy_run_analysis_states()?
            .len(),
        1
    );
    Ok(())
}
#[test]
fn live_analysis_blocks_upgrade_but_migrated_readers_can_open_during_analysis() -> Result<()> {
    let root = tempfile::tempdir()?;
    let db = legacy(root.path())?;
    db.upsert_run(&run(&db, "RUN1", RunState::DeDone))?;
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    assert!(
        Database::open(root.path())
            .unwrap_err()
            .to_string()
            .contains("idle cohort analysis")
    );
    assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::DeDone);
    lease.finish(JobOutcome::Failed, Some("test release"))?;
    drop(lease);
    let upgraded = Database::open(root.path())?;
    let mut lease = JobLease::try_acquire(&upgraded, "corr", None, &[])?.unwrap();
    assert_eq!(
        Database::open(root.path())?.get_run("RUN1")?.unwrap().state,
        RunState::QuantDone
    );
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(())
}
#[test]
fn guards_reject_old_writer_labels_without_changing_ingest_state() -> Result<()> {
    let root = tempfile::tempdir()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    let original = run(&db, "RUN1", RunState::QuantDone);
    db.upsert_run(&original)?;
    for state in [
        RunState::DeRunning,
        RunState::DeDone,
        RunState::DeFailed,
        RunState::CorrRunning,
        RunState::CorrDone,
        RunState::CorrFailed,
    ] {
        assert!(
            db.set_run_state("RUN1", state, Some("analysis error"))
                .is_err()
        );
        assert!(db.upsert_run(&run(&db, "NEW", state)).is_err());
    }
    assert_eq!(
        serde_json::to_value(db.get_run("RUN1")?.unwrap())?,
        serde_json::to_value(original)?
    );
    assert!(db.get_run("NEW")?.is_none());
    assert!(db.legacy_run_analysis_states()?.is_empty());
    Ok(())
}
#[test]
fn competing_upgrades_publish_one_receipt_and_preserve_original_data() -> Result<()> {
    let root = tempfile::tempdir()?;
    let db = legacy(root.path())?;
    db.upsert_run(&run(&db, "RUN1", RunState::CorrDone))?;
    let barrier = std::sync::Arc::new(std::sync::Barrier::new(8));
    std::thread::scope(|scope| -> Result<()> {
        let workers = (0..8)
            .map(|_| {
                let barrier = barrier.clone();
                let path = root.path().to_path_buf();
                scope.spawn(move || {
                    barrier.wait();
                    Database::open(&path)
                })
            })
            .collect::<Vec<_>>();
        for worker in workers {
            assert_eq!(
                worker.join().unwrap()?.legacy_run_analysis_states()?.len(),
                1
            );
        }
        Ok(())
    })?;
    assert!(fs::metadata(root.path().join("state.sqlite"))?.is_file());
    Ok(())
}
