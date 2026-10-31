use super::*;
use crate::{
    ArtifactKind, ProjectConfig, ProjectPaths,
    generations::{Generation, GenerationHead, OutputGroup},
    jobs::{JobLease, JobOutcome},
    util::now_rfc3339,
};

fn setup() -> Result<(tempfile::TempDir, Database)> {
    let root = tempfile::tempdir()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    Ok((root, db))
}
fn publish(db: &Database, value: &str, inputs: &[ManifestArtifact]) -> Result<GenerationHead> {
    let mut lease = JobLease::try_acquire(db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::De, "deseq2")?;
    let path = generation.head().payload_path.join("counts.tsv");
    fs::write(&path, value)?;
    let file = ManifestArtifact::snapshot("COUNTS", &path)?;
    let record = ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: None,
        kind: ArtifactKind::Counts,
        path: file.path,
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".into(),
        checksum: file.checksum,
        bytes: file.bytes,
        created_at: now_rfc3339(),
    };
    generation.seal_with_inputs(&[path], &"a".repeat(64), inputs)?;
    let head = generation.head().clone();
    generation.commit(&[record], |_| Ok(()))?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(head)
}
fn input(head: &GenerationHead) -> Result<ManifestArtifact> {
    ManifestArtifact::snapshot("COUNTS", &head.payload_path.join("counts.tsv"))
}
fn options(apply: bool) -> GenerationGcOptions {
    GenerationGcOptions {
        apply,
        purge_expired: false,
        trash_days: 7,
    }
}
fn entry<'a>(report: &'a GenerationGcReport, head: &GenerationHead) -> &'a GenerationGcEntry {
    assert!(!report.busy, "cleanup returned busy instead of a report");
    report
        .entries
        .iter()
        .find(|row| row.generation_id == head.generation_id)
        .unwrap_or_else(|| {
            panic!(
                "cleanup report omitted generation {}: {}",
                head.generation_id,
                serde_json::to_string(&report.entries).unwrap_or_default()
            )
        })
}

fn collect_when_idle(db: &Database, opts: GenerationGcOptions) -> Result<GenerationGcReport> {
    let start = std::time::Instant::now();
    loop {
        let report = db.generation_gc(opts)?;
        if !report.busy {
            return Ok(report);
        }
        assert!(
            start.elapsed() < std::time::Duration::from_secs(5),
            "generation cleanup remained busy"
        );
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
}

#[test]
fn current_heads_and_transitive_producers_survive_while_unreachable_history_is_reversible()
-> Result<()> {
    let (root, db) = setup()?;
    let first = publish(&db, "first", &[])?;
    let second = publish(&db, "second", &[input(&first)?])?;
    let current = publish(&db, "current", &[input(&second)?])?;
    let report = db.generation_gc(options(false))?;
    assert!(
        report
            .entries
            .iter()
            .all(|row| row.action == "held_reachable")
    );
    let independent = publish(&db, "new independent current", &[])?;
    let seal = fs::read(first.payload_path.parent().unwrap().join("seal.json"))?;
    let report = db.generation_gc(options(false))?;
    assert_eq!(entry(&report, &first).action, "trash");
    assert!(first.payload_path.exists());
    assert!(Retirement::all(&db)?.is_empty());
    let report = db.generation_gc(options(true))?;
    assert!(report.errors.is_empty(), "{:?}", report.errors);
    assert_eq!(entry(&report, &independent).action, "held_reachable");
    assert_eq!(entry(&report, &first).state, "TRASHED");
    assert!(!first.payload_path.exists());
    assert!(db.verify_generation_head(&independent)?);
    assert!(first.paths(&ProjectPaths::new(root.path())).is_err());
    assert!(db.prepare_generation_reuse(&first)?.is_none());
    assert_eq!(
        db.committed_generation_history(OutputGroup::De, db.project_id())?,
        vec![independent.clone()]
    );
    db.generation_restore(&current.generation_id)?;
    for head in [&first, &second, &current] {
        assert!(db.verify_generation_head(head)?);
    }
    assert_eq!(
        fs::read(first.payload_path.parent().unwrap().join("seal.json"))?,
        seal
    );
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?,
        Some(independent)
    );
    Ok(())
}
#[test]
fn snapshots_and_workers_block_collection_and_keeps_retain_history() -> Result<()> {
    let (root, db) = setup()?;
    let first = publish(&db, "old", &[])?;
    let _current = publish(&db, "current", &[])?;
    let snapshot = first.paths(&ProjectPaths::new(root.path()))?;
    assert!(db.generation_gc(options(true))?.busy);
    let cloned = snapshot.clone();
    drop(snapshot);
    assert!(db.generation_gc(options(true))?.busy);
    drop(cloned);
    let prepared = db.prepare_generation_reuse(&first)?.unwrap();
    assert!(db.generation_gc(options(true))?.busy);
    drop(prepared);
    let mut job = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    assert!(db.generation_gc(options(true))?.busy);
    job.finish(JobOutcome::Succeeded, None)?;
    drop(job);
    db.generation_keep(&first.generation_id, Some("publication"))?;
    assert_eq!(
        entry(&collect_when_idle(&db, options(true))?, &first).action,
        "held_reachable"
    );
    db.generation_keep(&first.generation_id, None)?;
    assert_eq!(
        entry(&collect_when_idle(&db, options(true))?, &first).state,
        "TRASHED"
    );
    Ok(())
}
#[test]
fn purge_uses_recorded_deadline_and_keeps_required_trash_dependencies() -> Result<()> {
    let (_root, db) = setup()?;
    let first = publish(&db, "producer", &[])?;
    let consumer = publish(&db, "consumer", &[input(&first)?])?;
    publish(&db, "current", &[])?;
    db.generation_gc(options(true))?;
    let mut purge = options(true);
    purge.purge_expired = true;
    purge.trash_days = 0;
    assert!(
        db.generation_gc(purge)?
            .entries
            .iter()
            .all(|row| row.action != "purge")
    );
    db.write(|conn|{conn.execute("UPDATE generation_retirements SET trashed_at_ms=0,purge_after_ms=1 WHERE generation_id=?1",[&first.generation_id])?;Ok(())})?;
    assert_ne!(entry(&db.generation_gc(purge)?, &first).action, "purge");
    db.write(|conn|{conn.execute("UPDATE generation_retirements SET trashed_at_ms=0,purge_after_ms=1 WHERE generation_id=?1",[&consumer.generation_id])?;Ok(())})?;
    let report = db.generation_gc(purge)?;
    assert_eq!(entry(&report, &first).state, "PURGED");
    assert_eq!(entry(&report, &consumer).state, "PURGED");
    assert!(db.generation_restore(&consumer.generation_id).is_err());
    assert!(
        db.committed_generation_head(&first.generation_id)?
            .is_some(),
        "metadata remains an inspectable tombstone"
    );
    Ok(())
}
#[test]
fn corrupted_unreachable_payloads_keep_their_dependencies_and_are_not_deleted() -> Result<()> {
    let (_root, db) = setup()?;
    let first = publish(&db, "producer", &[])?;
    let broken = publish(&db, "broken consumer", &[input(&first)?])?;
    publish(&db, "current", &[])?;
    fs::write(broken.payload_path.join("counts.tsv"), "corrupted")?;
    let report = db.generation_gc(options(true))?;
    assert!(!report.errors.is_empty());
    assert_eq!(entry(&report, &broken).action, "held_invalid");
    assert_eq!(entry(&report, &first).action, "held_reachable");
    assert!(first.payload_path.exists());
    assert!(broken.payload_path.exists());
    Ok(())
}

#[test]
fn legacy_unknown_ancestry_and_broken_catalog_indexes_never_authorize_collection() -> Result<()> {
    let (_root, db) = setup()?;
    let old = publish(&db, "legacy", &[])?;
    publish(&db, "current", &[])?;
    let mut seal = db.committed_seal(&old.generation_id)?.unwrap();
    seal.lineage_version = 0;
    db.write(|conn| {
        conn.execute(
            "UPDATE artifact_generations SET inventory_json=?1 WHERE generation_id=?2",
            params![serde_json::to_string(&seal)?, old.generation_id],
        )?;
        Ok(())
    })?;
    crate::util::write_json_pretty(&old.payload_path.parent().unwrap().join("seal.json"), &seal)?;
    let report = db.generation_gc(options(true))?;
    assert!(
        report
            .entries
            .iter()
            .all(|entry| entry.action == "held_reachable")
    );
    assert!(
        entry(&report, &old)
            .reasons
            .contains(&"unproven_legacy_ancestry".into())
    );
    db.write(|conn| {
        conn.execute(
            "UPDATE generation_files SET sha256=?1 WHERE generation_id=?2",
            params!["0".repeat(64), old.generation_id],
        )?;
        Ok(())
    })?;
    assert!(db.generation_gc(options(true)).is_err());
    assert!(old.payload_path.exists());
    assert!(Retirement::all(&db)?.is_empty());
    Ok(())
}

#[test]
#[cfg(unix)]
fn retargeted_recorded_source_alias_protects_current_data_without_rewriting_history() -> Result<()>
{
    let (root, db) = setup()?;
    let external = tempfile::tempdir()?;
    let source = external.path().join("source");
    fs::write(&source, "original external data")?;
    let old = publish(&db, "older managed data", &[])?;
    let alias = root.path().join("source-alias");
    std::os::unix::fs::symlink(&source, &alias)?;
    let snapshot = ManifestArtifact::snapshot("SOURCE", &alias)?;
    let current = publish(&db, "current result", &[snapshot])?;
    let before = db.generation_lineage(&current)?;
    assert!(before.dependencies.is_empty());
    fs::remove_file(&alias)?;
    std::os::unix::fs::symlink(old.payload_path.join("counts.tsv"), &alias)?;
    let report = db.generation_gc(options(true))?;
    assert_eq!(entry(&report, &old).action, "held_reachable");
    assert!(old.payload_path.exists());
    assert_eq!(db.generation_lineage(&current)?, before);
    Ok(())
}

#[cfg(unix)]
#[test]
fn active_aliases_keep_their_generation_and_dangling_aliases_cannot_reactivate_trash() -> Result<()>
{
    let (root, db) = setup()?;
    let old = publish(&db, "old", &[])?;
    publish(&db, "current", &[])?;
    let alias = root.path().join("view");
    std::os::unix::fs::symlink(&old.payload_path, &alias)?;
    let file = input(&old)?;
    let record = ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: None,
        kind: ArtifactKind::Counts,
        path: alias.join("counts.tsv").display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: file.checksum_type,
        checksum: file.checksum,
        bytes: file.bytes,
        created_at: now_rfc3339(),
    };
    db.record_artifact(&record)?;
    assert_eq!(
        entry(&db.generation_gc(options(true))?, &old).action,
        "held_reachable"
    );
    db.delete_artifact(&record.path)?;
    assert_eq!(
        entry(&db.generation_gc(options(true))?, &old).state,
        "TRASHED"
    );
    assert!(db.record_artifact(&record).is_err());
    assert!(!old.payload_path.exists());
    db.generation_restore(&old.generation_id)?;
    db.record_artifact(&record)?;
    assert_eq!(
        entry(&db.generation_gc(options(true))?, &old).action,
        "held_reachable"
    );
    Ok(())
}

#[test]
fn maintenance_crash_worker() -> Result<()> {
    let Ok(root) = std::env::var("RNAA_GENERATION_GC_TEST_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let db = Database::open(root)?;
    let phase = std::env::var("RNAA_GENERATION_GC_TEST_PHASE")?;
    let id = std::env::var("RNAA_GENERATION_GC_TEST_ID")?;
    let _reader;
    let _maintenance;
    if phase == "reader" {
        _reader = Some(super::super::GenerationStorePin::acquire(
            db.path(),
            "crash reader",
        )?);
        _maintenance = None;
    } else {
        _reader = None;
        _maintenance = Some(super::super::pins::maintenance_lock(&db)?.unwrap());
        let mut db = db.clone();
        db.generation_maintenance = true;
        let graph = Graph::capture(&db)?;
        if phase.starts_with("trash") {
            let observed =
                BTreeMap::from([(id.clone(), db.observe_live_generation(&graph.nodes[&id])?)]);
            db.mark_generation_trash(&graph, std::slice::from_ref(&id), &observed, 7)?;
            if phase == "trash_moved" {
                let graph = Graph::capture(&db)?;
                db.move_generation(&graph.retired[&id], &graph, false)?;
            }
        } else if phase == "restore_moved" {
            db.write(|conn| {
                conn.execute(
                    "UPDATE generation_retirements SET state='RESTORING' WHERE generation_id=?1",
                    [&id],
                )?;
                Ok(())
            })?;
            let graph = Graph::capture(&db)?;
            db.move_generation(&graph.retired[&id], &graph, true)?;
        } else if phase == "purge_partial" {
            db.mark_generation_purge(&graph, std::slice::from_ref(&id))?;
            fs::remove_file(
                db.generation_container(&id, true, false)?
                    .join("payload/counts.tsv"),
            )?;
        } else {
            bail!("unknown cleanup crash phase");
        }
    }
    fs::write(root.join("ready"), b"ready")?;
    loop {
        std::thread::sleep(std::time::Duration::from_secs(1));
    }
}

#[test]
fn killed_readers_and_collectors_resume_durable_move_restore_and_partial_purge_intents()
-> Result<()> {
    struct Child(std::process::Child);
    impl Drop for Child {
        fn drop(&mut self) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
    for phase in [
        "reader",
        "trash",
        "trash_moved",
        "restore_moved",
        "purge_partial",
    ] {
        let (root, db) = setup()?;
        let old = publish(&db, "old", &[])?;
        let current = publish(&db, "current", &[])?;
        if matches!(phase, "restore_moved" | "purge_partial") {
            let mut opts = options(true);
            opts.trash_days = 0;
            db.generation_gc(opts)?;
        }
        let mut child = Child(
            std::process::Command::new(std::env::current_exe()?)
                .args([
                    "--exact",
                    "generations::gc::tests::maintenance_crash_worker",
                    "--nocapture",
                ])
                .env("RNAA_GENERATION_GC_TEST_ROOT", root.path())
                .env("RNAA_GENERATION_GC_TEST_PHASE", phase)
                .env("RNAA_GENERATION_GC_TEST_ID", &old.generation_id)
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::piped())
                .spawn()?,
        );
        let started = std::time::Instant::now();
        while !root.path().join("ready").exists() {
            assert!(
                child.0.try_wait()?.is_none(),
                "crash worker exited at {phase}"
            );
            assert!(started.elapsed() < std::time::Duration::from_secs(15));
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        assert!(db.generation_gc(options(false))?.busy);
        child.0.kill()?;
        child.0.wait()?;
        if phase == "restore_moved" {
            db.generation_restore(&old.generation_id)?;
            assert!(db.verify_generation_head(&old)?);
        } else {
            let report = db.generation_gc(options(true))?;
            assert!(report.errors.is_empty(), "{:?}", report.errors);
            assert_eq!(
                entry(&report, &old).state,
                if phase == "purge_partial" {
                    "PURGED"
                } else {
                    "TRASHED"
                }
            );
        }
        assert!(db.verify_generation_head(&current)?);
        assert!(
            db.generation_storage_status()?["pin_observations"]
                .as_array()
                .unwrap()
                .is_empty()
        );
    }
    Ok(())
}
