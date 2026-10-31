//! Failure boundaries use real child processes so locks are released by death.
use super::*;
use crate::{
    RunRecord, RunState,
    generations::{Generation, GenerationHead},
    jobs::{JobLease, JobOutcome},
    util::now_rfc3339,
};

fn options(apply: bool) -> FileCleanupOptions {
    FileCleanupOptions {
        mode: CleanupMode::Fastq,
        run: None,
        apply,
        purge_expired: false,
        trash_days: 7,
        expected_policy: None,
    }
}

fn cleanup_when_idle(db: &Database, options: FileCleanupOptions) -> Result<FileCleanupReport> {
    let start = std::time::Instant::now();
    loop {
        let report = db.cleanup_files(options.clone())?;
        if !report.busy {
            return Ok(report);
        }
        assert!(
            start.elapsed() < std::time::Duration::from_secs(5),
            "cleanup remained busy"
        );
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
}
fn record(db: &Database, kind: ArtifactKind, path: &Path) -> Result<ArtifactRecord> {
    let file = ManifestArtifact::snapshot(kind.as_str(), path)?;
    Ok(ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: Some("RUN1".into()),
        kind,
        path: file.path,
        blob_id: None,
        shared_path: None,
        checksum_type: file.checksum_type,
        checksum: file.checksum,
        bytes: file.bytes,
        created_at: now_rfc3339(),
    })
}
fn publish_counts(
    db: &Database,
    inputs: &[ManifestArtifact],
    value: &str,
) -> Result<GenerationHead> {
    let mut lease = JobLease::try_acquire(
        db,
        "quant",
        Some("RUN1"),
        &[RunState::Verified, RunState::QuantDone],
    )?
    .unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::Quant, "quant")?;
    let path = generation
        .head()
        .payload_path
        .join("rnaa_quant_manifest.json");
    fs::write(&path, value)?;
    let rec = record(db, ArtifactKind::Manifest, &path)?;
    generation.seal_with_inputs(&[path], &"a".repeat(64), inputs)?;
    let head = generation.head().clone();
    generation.commit(&[rec], |db| {
        db.set_run_state("RUN1", RunState::QuantDone, None)
    })?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(head)
}
fn setup() -> Result<(tempfile::TempDir, Database, GenerationHead, GenerationHead)> {
    let root = tempfile::tempdir()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    db.upsert_run(&RunRecord {
        project_id: db.project_id().into(),
        study_accession: "study".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Paired,
        instrument_platform: None,
        instrument_model: None,
        remote_files: Vec::new(),
        metadata: Default::default(),
        state: RunState::Verified,
        last_error: None,
        updated_at: now_rfc3339(),
    })?;
    let mut lease =
        JobLease::try_acquire(&db, "download", Some("RUN1"), &[RunState::Verified])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"b".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::Download, "download")?;
    let mut records = Vec::new();
    let mut inputs = Vec::new();
    let mut required = Vec::new();
    for (kind, name) in [
        (ArtifactKind::FastqR1, "R1.fastq"),
        (ArtifactKind::FastqR2, "R2.fastq"),
    ] {
        let path = generation.head().payload_path.join(name);
        fs::write(&path, name)?;
        records.push(record(&db, kind, &path)?);
        inputs.push(ManifestArtifact::snapshot(kind.as_str(), &path)?);
        required.push(path);
    }
    generation.seal(&required, &"b".repeat(64))?;
    let producer = generation.head().clone();
    generation.commit(&records, |db| {
        db.set_run_state("RUN1", RunState::Verified, None)
    })?;
    lease.finish(JobOutcome::Succeeded, None)?;
    drop(generation);
    drop(owned);
    drop(lease);
    let counts = publish_counts(&db, &inputs, "first")?;
    Ok((root, db, producer, counts))
}
fn prepare(db: &Database) -> Result<Vec<Descriptor>> {
    let records = db.list_artifacts_for_run(Some("RUN1"))?;
    let multiplicity = records
        .iter()
        .map(|record| (record.path.clone(), 1))
        .collect();
    let operation = uuid::Uuid::new_v4().to_string();
    records
        .iter()
        .filter(|record| fastq(record.kind))
        .map(|record| {
            db.prepare_trash_descriptor(record, &operation, CleanupMode::Fastq, 0, &multiplicity)
        })
        .collect()
}

#[test]
fn pins_ownership_and_damaged_members_hold_complete_paired_groups() -> Result<()> {
    let (_root, db, producer, _) = setup()?;
    let pin = super::super::GenerationStorePin::acquire(db.path(), "test")?;
    assert!(db.cleanup_files(options(true))?.busy);
    drop(pin);
    let path = producer.payload_path.join("R2.fastq");
    fs::write(&path, "changed")?;
    let report = db.cleanup_files(options(true))?;
    assert!(!report.errors.is_empty());
    assert!(producer.payload_path.join("R1.fastq").exists());
    assert!(db.file_trash_rows()?.is_empty());
    fs::write(&path, "R2.fastq")?;
    let mate = db
        .list_artifacts_for_run(Some("RUN1"))?
        .into_iter()
        .find(|row| row.kind == ArtifactKind::FastqR2)
        .unwrap();
    db.delete_artifact(&mate.path)?;
    assert!(!db.cleanup_files(options(true))?.errors.is_empty());
    assert!(producer.payload_path.join("R1.fastq").exists());
    assert!(db.file_trash_rows()?.is_empty());
    db.record_artifact(&mate)?;
    let alias = producer
        .payload_path
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("hardlink.fastq");
    fs::hard_link(&path, &alias)?;
    assert!(!db.cleanup_files(options(true))?.errors.is_empty());
    assert!(db.file_trash_rows()?.is_empty());
    fs::remove_file(alias)?;
    assert!(db.cleanup_files(options(true))?.errors.is_empty());
    assert!(db.verify_generation_head(&producer)?);
    Ok(())
}

#[test]
fn restorable_inputs_pin_old_count_proofs_and_generation_gc_preserves_purged_absence_evidence()
-> Result<()> {
    let (_root, db, producer, counts) = setup()?;
    let mut opts = options(true);
    opts.trash_days = 0;
    db.cleanup_files(opts)?;
    publish_counts(&db, &[], "replacement independent counts")?;
    let preview = db.generation_gc(super::super::GenerationGcOptions {
        apply: false,
        purge_expired: false,
        trash_days: 7,
    })?;
    assert_eq!(
        preview
            .entries
            .iter()
            .find(|row| row.generation_id == counts.generation_id)
            .unwrap()
            .action,
        "held_reachable"
    );
    let mut opts = options(true);
    opts.mode = CleanupMode::None;
    opts.purge_expired = true;
    db.cleanup_files(opts)?;
    assert!(db.verify_generation_head(&producer)?);
    assert!(!db.restorable_file_trash("RUN1")?);
    let preview = db.generation_gc(super::super::GenerationGcOptions {
        apply: false,
        purge_expired: false,
        trash_days: 7,
    })?;
    assert_eq!(
        preview
            .entries
            .iter()
            .find(|row| row.generation_id == counts.generation_id)
            .unwrap()
            .action,
        "held_reachable"
    );
    Ok(())
}

#[test]
fn changed_receipt_or_payload_never_authorizes_purge() -> Result<()> {
    let (_root, db, _, _) = setup()?;
    let mut opts = options(true);
    opts.trash_days = 0;
    db.cleanup_files(opts)?;
    let row = db.file_trash_rows()?.remove(0);
    let path = db.file_trash_path(&row.descriptor.trash_id, false)?;
    let original = fs::read(&path)?;
    fs::write(&path, "corruption")?;
    let mut opts = options(true);
    opts.mode = CleanupMode::None;
    opts.purge_expired = true;
    assert!(db.cleanup_files(opts.clone()).is_err());
    assert_eq!(fs::read(&path)?, b"corruption");
    fs::write(&path, original)?;
    let witness = db.file_receipt_path(&row.descriptor.trash_id, false)?;
    fs::write(witness, b"{}")?;
    assert!(db.cleanup_files(opts).is_err());
    assert!(path.exists());
    Ok(())
}

#[test]
fn existing_deadlines_survive_later_options_and_pending_recovery_does_not_select_new_files()
-> Result<()> {
    let (_root, db, producer, _) = setup()?;
    let initial = cleanup_when_idle(&db, options(true))?;
    assert!(
        !initial.busy && initial.errors.is_empty(),
        "initial cleanup failed: {initial:?}"
    );
    let mut opts = options(true);
    opts.mode = CleanupMode::None;
    opts.trash_days = 0;
    opts.purge_expired = true;
    let later = cleanup_when_idle(&db, opts)?;
    assert!(
        !later.busy && later.errors.is_empty(),
        "subsequent cleanup failed: {later:?}"
    );
    let rows = db.file_trash_rows()?;
    assert_eq!(
        rows.len(),
        2,
        "unexpected trash rows after {initial:?}; {later:?}"
    );
    assert!(rows.iter().all(|row| row.state == "TRASHED"));
    db.restore_file_trash(&rows[0].descriptor.trash_id)?;
    assert!(db.recover_input_cleanup()?);
    assert!(producer.payload_path.join("R1.fastq").exists());
    assert!(
        db.file_trash_rows()?
            .iter()
            .all(|row| row.state == "RESTORED")
    );
    Ok(())
}

#[test]
fn managed_archive_absence_requires_count_proof_and_external_archives_hold_the_run_group()
-> Result<()> {
    let (_root, db, _, counts) = setup()?;
    let mut lease =
        JobLease::try_acquire(&db, "download", Some("RUN1"), &[RunState::QuantDone])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"c".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::Download, "download")?;
    let path = generation.head().payload_path.join("reads.sra");
    fs::write(&path, b"archive")?;
    let archive = record(&db, ArtifactKind::Sra, &path)?;
    generation.seal(std::slice::from_ref(&path), &"c".repeat(64))?;
    let producer = generation.head().clone();
    generation.commit(&[archive], |db| {
        db.set_run_state("RUN1", RunState::QuantDone, None)
    })?;
    lease.finish(JobOutcome::Succeeded, None)?;
    drop(generation);
    drop(owned);
    drop(lease);
    let external = tempfile::NamedTempFile::new()?;
    fs::write(external.path(), b"external archive")?;
    let external_record = record(&db, ArtifactKind::Sra, external.path())?;
    db.record_artifact(&external_record)?;
    let mut opts = options(true);
    opts.mode = CleanupMode::Sra;
    let report = db.cleanup_files(opts.clone())?;
    assert!(!report.errors.is_empty());
    assert!(path.exists());
    assert!(db.file_trash_rows()?.is_empty());
    db.delete_artifact(&external_record.path)?;
    assert!(db.cleanup_files(opts)?.errors.is_empty());
    assert!(!path.exists());
    assert!(db.verify_generation_head(&producer)?);
    assert!(db.verify_generation_head(&counts)?);
    let id = db.file_trash_rows()?[0].descriptor.trash_id.clone();
    db.restore_file_trash(&id)?;
    assert_eq!(fs::read(path)?, b"archive");
    assert_eq!(fs::read(external.path())?, b"external archive");
    Ok(())
}

#[test]
fn crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_FILE_TRASH_TEST_ROOT") else {
        return Ok(());
    };
    let db = Database::open(Path::new(&root))?;
    let _lock = super::super::pins::maintenance_lock(&db)?.context("test maintenance lock")?;
    let mut db = db;
    db.generation_maintenance = true;
    let phase = std::env::var("RNAA_FILE_TRASH_TEST_PHASE")?;
    if phase.starts_with("trash") {
        let values = prepare(&db)?;
        db.mark_file_trash(&values)?;
        if phase == "trash_moved" {
            db.move_cleanup_file(&values[0], false)?;
        }
        if phase == "trash_copied" {
            let target = db.file_trash_path(&values[0].trash_id, true)?;
            fs::copy(&values[0].canonical_source, &target)?;
            std::fs::File::open(target)?.sync_all()?;
        }
        if phase == "trash_partial_copy" {
            db.file_trash_path(&values[0].trash_id, true)?;
            fs::write(db.copy_staging_path(&values[0], false)?, "incomplete copy")?;
        }
    } else {
        let rows = db.file_trash_rows()?;
        if phase == "restore_moved" || phase == "restore_partial_copy" {
            db.write(|conn| {
                conn.execute("UPDATE file_trash SET state='RESTORING'", [])?;
                Ok(())
            })?;
            if phase == "restore_moved" {
                db.move_cleanup_file(&rows[0].descriptor, true)?;
            } else {
                fs::write(
                    db.copy_staging_path(&rows[0].descriptor, true)?,
                    "incomplete copy",
                )?;
            }
        } else if phase == "purge_unlinked" {
            db.mark_file_purge(&rows)?;
            fs::remove_file(db.file_trash_path(&rows[0].descriptor.trash_id, false)?)?;
        } else {
            bail!("unknown file cleanup crash phase");
        }
    }
    fs::write(Path::new(&root).join("ready"), b"ready")?;
    loop {
        std::thread::sleep(std::time::Duration::from_secs(1));
    }
}

#[test]
fn killed_collectors_resume_partial_paired_moves_restores_copies_and_unlinks() -> Result<()> {
    struct Child(std::process::Child);
    impl Drop for Child {
        fn drop(&mut self) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
    for phase in [
        "trash",
        "trash_moved",
        "trash_copied",
        "trash_partial_copy",
        "restore_moved",
        "restore_partial_copy",
        "purge_unlinked",
    ] {
        let (root, db, producer, counts) = setup()?;
        if matches!(
            phase,
            "restore_moved" | "restore_partial_copy" | "purge_unlinked"
        ) {
            let mut opts = options(true);
            opts.trash_days = 0;
            db.cleanup_files(opts)?;
        }
        let mut child = Child(
            std::process::Command::new(std::env::current_exe()?)
                .args([
                    "--exact",
                    "generations::file_trash::tests::crash_worker",
                    "--nocapture",
                ])
                .env("RNAA_FILE_TRASH_TEST_ROOT", root.path())
                .env("RNAA_FILE_TRASH_TEST_PHASE", phase)
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::piped())
                .spawn()?,
        );
        let start = std::time::Instant::now();
        while !root.path().join("ready").exists() {
            assert!(child.0.try_wait()?.is_none(), "worker exited at {phase}");
            assert!(
                start.elapsed() < std::time::Duration::from_secs(15),
                "worker stalled at {phase}"
            );
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        assert!(db.cleanup_files(options(false))?.busy);
        child.0.kill()?;
        child.0.wait()?;
        let mut opts = options(true);
        opts.mode = CleanupMode::None;
        db.cleanup_files(opts)?;
        let rows = db.file_trash_rows()?;
        let expected = match phase {
            "restore_moved" | "restore_partial_copy" => "RESTORED",
            "purge_unlinked" => "PURGED",
            _ => "TRASHED",
        };
        assert!(rows.iter().all(|row| row.state == expected), "{phase}");
        assert!(db.verify_generation_head(&producer)?);
        assert!(db.verify_generation_head(&counts)?);
        if phase == "restore_moved" {
            assert_eq!(
                db.list_artifacts_for_run(Some("RUN1"))?
                    .iter()
                    .filter(|row| fastq(row.kind))
                    .count(),
                2
            );
        }
    }
    Ok(())
}
