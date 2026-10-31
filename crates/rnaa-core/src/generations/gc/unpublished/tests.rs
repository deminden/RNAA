use super::*;
mod crashes;
use crate::{
    ProjectConfig,
    generations::{Generation, GenerationHead, OutputGroup},
    jobs::{JobLease, JobOutcome},
};

fn setup() -> Result<(tempfile::TempDir, Database)> {
    let root = tempfile::tempdir()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    Ok((root, db))
}
fn abandoned(db: &Database, sealed: bool) -> Result<GenerationHead> {
    let mut lease = JobLease::try_acquire(db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::De, "deseq2")?;
    let path = generation.head().payload_path.join("partial.tsv");
    fs::write(&path, "diagnostic data")?;
    fs::create_dir(generation.head().payload_path.join("empty"))?;
    if sealed {
        generation.seal(&[path], &"a".repeat(64))?;
    }
    let head = generation.head().clone();
    drop(generation);
    lease.finish(
        JobOutcome::Failed,
        Some("failed after producing diagnostic files"),
    )?;
    drop(lease);
    Ok(head)
}
fn options(apply: bool) -> GenerationGcOptions {
    GenerationGcOptions {
        apply,
        purge_expired: false,
        trash_days: 7,
    }
}
fn entry<'a>(
    report: &'a UnpublishedCleanupReport,
    head: &GenerationHead,
) -> &'a UnpublishedCleanupEntry {
    report
        .entries
        .iter()
        .find(|entry| entry.generation_id == head.generation_id)
        .unwrap()
}
#[test]
fn preview_trash_restore_and_explicit_expiry_preserve_unpublished_identity() -> Result<()> {
    let (root, db) = setup()?;
    let unfinished = abandoned(&db, true)?;
    let alias = root.path().join("diagnostic-hardlink");
    fs::hard_link(unfinished.payload_path.join("partial.tsv"), &alias)?;
    let seal = fs::read(unfinished.payload_path.parent().unwrap().join("seal.json"))?;
    let before = GenerationRow::unpublished(&db)?;
    let preview = db.unpublished_generation_gc(options(false))?;
    assert_eq!(entry(&preview, &unfinished).action, "trash");
    assert!(Receipt::all(&db)?.is_empty());
    assert!(unfinished.payload_path.exists());
    let applied = db.unpublished_generation_gc(options(true))?;
    assert!(applied.errors.is_empty(), "{:?}", applied.errors);
    assert_eq!(entry(&applied, &unfinished).state, "TRASHED");
    assert!(!unfinished.payload_path.exists());
    assert_eq!(GenerationRow::unpublished(&db)?, before);
    assert!(db.generation_heads()?.is_empty());
    let mut purge = options(true);
    purge.purge_expired = true;
    purge.trash_days = 0;
    assert_eq!(
        entry(&db.unpublished_generation_gc(purge)?, &unfinished).state,
        "TRASHED",
        "later options cannot shorten the recorded window"
    );
    db.restore_unpublished_generation(&unfinished.generation_id)?;
    assert!(unfinished.payload_path.join("empty").is_dir());
    assert_eq!(
        fs::read(unfinished.payload_path.parent().unwrap().join("seal.json"))?,
        seal
    );
    assert_eq!(GenerationRow::unpublished(&db)?, before);
    assert!(db.generation_heads()?.is_empty());
    assert!(
        db.committed_generation_head(&unfinished.generation_id)?
            .is_none()
    );
    let again = db.unpublished_generation_gc(purge)?;
    assert_eq!(
        entry(&again, &unfinished).state,
        "TRASHED",
        "new zero-window trash must survive its first invocation"
    );
    let usage = db.storage_usage(None, true)?;
    let json = serde_json::to_value(usage)?;
    assert!(json.to_string().contains("unpublished_trash"));
    assert_eq!(
        entry(&db.unpublished_generation_gc(purge)?, &unfinished).state,
        "PURGED"
    );
    assert!(
        db.restore_unpublished_generation(&unfinished.generation_id)
            .is_err()
    );
    assert!(
        !root
            .path()
            .join("trash/unpublished")
            .join(&unfinished.generation_id)
            .exists()
    );
    assert!(db.generation_heads()?.is_empty());
    assert_eq!(fs::read(alias)?, b"diagnostic data");
    Ok(())
}
#[test]
fn pins_live_attempts_config_sources_and_unknown_directories_are_preserved() -> Result<()> {
    let (root, db) = setup()?;
    let head = abandoned(&db, false)?;
    let pin = crate::generations::GenerationStorePin::acquire(db.path(), "test reader")?;
    assert!(db.unpublished_generation_gc(options(true))?.busy);
    drop(pin);
    let mut job = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    assert!(db.unpublished_generation_gc(options(true))?.busy);
    job.finish(JobOutcome::Succeeded, None)?;
    drop(job);
    db.write(|conn| {
        conn.execute(
            "UPDATE job_attempts SET state='RUNNING' WHERE attempt_id=?1",
            [&head.attempt_id],
        )?;
        Ok(())
    })?;
    assert_eq!(
        entry(&db.unpublished_generation_gc(options(true))?, &head).action,
        "held_invalid"
    );
    db.write(|conn| {
        conn.execute(
            "UPDATE job_attempts SET state='FAILED' WHERE attempt_id=?1",
            [&head.attempt_id],
        )?;
        Ok(())
    })?;
    let config = root.path().join("rnaa.toml");
    fs::write(
        &config,
        format!(
            "[refs]\ncustom_cdna = {:?}\n",
            head.payload_path.join("partial.tsv").display().to_string()
        ),
    )?;
    assert_eq!(
        entry(&db.unpublished_generation_gc(options(true))?, &head).action,
        "held_referenced"
    );
    // A colon is legal in a Unix filename. An absolute local path must not
    // become an ignored URL just because a later component contains "://".
    #[cfg(unix)]
    {
        fs::write(
            &config,
            format!(
                "[refs]\ncustom_cdna = {:?}\n",
                format!("{}/colon://source.fa", head.payload_path.display())
            ),
        )?;
        assert_eq!(
            entry(&db.unpublished_generation_gc(options(true))?, &head).action,
            "held_referenced"
        );
    }
    fs::remove_file(config)?;
    let orphan = root
        .path()
        .join(".rnaa_generations")
        .join(uuid::Uuid::new_v4().to_string());
    fs::create_dir(&orphan)?;
    fs::write(orphan.join("unknown"), "unowned")?;
    let report = db.unpublished_generation_gc(options(true))?;
    assert_eq!(entry(&report, &head).state, "TRASHED");
    assert!(
        report
            .entries
            .iter()
            .any(|entry| entry.path == fs::canonicalize(&orphan).unwrap()
                && entry.action == "held_unindexed")
    );
    assert!(orphan.join("unknown").exists());
    Ok(())
}
#[cfg(unix)]
#[test]
fn aliases_and_legacy_artifacts_hold_unpublished_data_and_symlinks_never_authorize_moves()
-> Result<()> {
    let (root, db) = setup()?;
    let head = abandoned(&db, false)?;
    let alias = root.path().join("source-alias");
    std::os::unix::fs::symlink(head.payload_path.join("partial.tsv"), &alias)?;
    let mut config = ProjectConfig::default();
    config.refs.custom_cdna = alias.display().to_string();
    db.set_project_config(&config)?;
    assert_eq!(
        entry(&db.unpublished_generation_gc(options(true))?, &head).action,
        "held_referenced"
    );
    db.set_project_config(&ProjectConfig::default())?;
    let other = root.path().join("other");
    fs::write(&other, "independent")?;
    db.record_artifact(&crate::ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: None,
        kind: crate::ArtifactKind::Counts,
        path: other.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".into(),
        checksum: compute_sha256(&other)?,
        bytes: 11,
        created_at: now_rfc3339(),
    })?;
    assert!(
        entry(&db.unpublished_generation_gc(options(true))?, &head)
            .reasons
            .contains(&"unproven_legacy_ancestry".into())
    );
    db.delete_artifact(other.to_str().unwrap())?;
    let inside = head.payload_path.join("escape");
    std::os::unix::fs::symlink(&other, &inside)?;
    assert_eq!(
        entry(&db.unpublished_generation_gc(options(true))?, &head).action,
        "held_invalid"
    );
    assert!(other.exists());
    assert!(head.payload_path.exists());
    assert!(Receipt::all(&db)?.is_empty());
    Ok(())
}
#[test]
fn changed_trash_or_catalog_and_transaction_failure_do_not_authorize_deletion() -> Result<()> {
    let (_root, db) = setup()?;
    let first = abandoned(&db, false)?;
    let second = abandoned(&db, true)?;
    db.write(|conn|{conn.execute_batch("CREATE TRIGGER reject_unpublished BEFORE INSERT ON unpublished_generation_trash BEGIN SELECT RAISE(ABORT,'injected'); END;")?;Ok(())})?;
    assert!(db.unpublished_generation_gc(options(true)).is_err());
    assert!(Receipt::all(&db)?.is_empty());
    assert!(first.payload_path.exists());
    assert!(second.payload_path.exists());
    db.write(|conn| {
        conn.execute_batch("DROP TRIGGER reject_unpublished;")?;
        Ok(())
    })?;
    let mut opts = options(true);
    opts.trash_days = 0;
    db.unpublished_generation_gc(opts)?;
    let path = tree::container(&db, &first.generation_id, true, false)?.join("payload/partial.tsv");
    fs::write(&path, "changed")?;
    opts.purge_expired = true;
    let report = db.unpublished_generation_gc(opts)?;
    assert_eq!(entry(&report, &first).action, "held_invalid");
    assert!(path.exists());
    assert_eq!(entry(&report, &second).state, "PURGED");
    assert!(
        db.restore_unpublished_generation(&first.generation_id)
            .is_err()
    );
    db.write(|conn| {
        conn.execute(
            "UPDATE artifact_generations SET subject='altered' WHERE generation_id=?1",
            [&first.generation_id],
        )?;
        Ok(())
    })?;
    assert!(db.unpublished_generation_gc(opts).is_err());
    assert!(path.exists());
    Ok(())
}

#[test]
fn newly_configured_sources_hold_existing_trash_and_restore_keeps_diagnostic_data() -> Result<()> {
    let (_root, db) = setup()?;
    let head = abandoned(&db, false)?;
    let mut opts = options(true);
    opts.trash_days = 0;
    db.unpublished_generation_gc(opts)?;
    let path = tree::container(&db, &head.generation_id, true, false)?.join("payload/partial.tsv");
    let mut config = ProjectConfig::default();
    config.refs.custom_cdna = path.display().to_string();
    db.set_project_config(&config)?;
    opts.purge_expired = true;
    let report = db.unpublished_generation_gc(opts)?;
    assert_eq!(
        entry(&report, &head).action,
        "held_trash_window_or_reference"
    );
    assert!(path.exists());
    db.restore_unpublished_generation(&head.generation_id)?;
    assert!(head.payload_path.join("partial.tsv").exists());
    assert!(Receipt::all(&db)?.is_empty());
    Ok(())
}
