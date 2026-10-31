use super::*;
use crate::{ArtifactKind, ArtifactRecord, ProjectConfig};

fn project(root: &Path) -> Result<Database> {
    fs::create_dir_all(root)?;
    Database::create(root, &ProjectConfig::default())
}
fn source(root: &Path) -> Result<(PathBuf, String, u64)> {
    let path = root.join("abundance.h5");
    fs::write(&path, b"shared immutable test payload")?;
    Ok((path.clone(), compute_sha256(&path)?, file_size(&path)?))
}
fn register(db: &Database, publication: &BlobPublication) -> Result<()> {
    let blob = publication.record();
    db.record_shared_blob(blob)?;
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: None,
        kind: ArtifactKind::QuantAbundanceH5,
        path: publication.owner.artifact_path.clone(),
        blob_id: Some(blob.blob_id.clone()),
        shared_path: Some(blob.storage_path.clone()),
        checksum_type: "sha256".into(),
        checksum: blob.checksum.clone(),
        bytes: blob.bytes,
        created_at: now_rfc3339(),
    })?;
    publication.finish()
}
fn options(apply: bool) -> GcOptions {
    GcOptions {
        apply,
        purge_expired: false,
        trash_days: 7,
    }
}

#[test]
fn two_projects_share_one_blob_and_one_removed_view_cannot_collect_the_other() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let first = project(&root.path().join("one"))?;
    let second = project(&root.path().join("two"))?;
    let (path, hash, bytes) = source(root.path())?;
    let publication = store.prepare(&first, &path, &hash, bytes)?;
    register(&first, &publication)?;
    drop(publication);
    let publication = store.prepare(&second, &path, &hash, bytes)?;
    register(&second, &publication)?;
    drop(publication);
    first.delete_artifact(path.to_str().unwrap())?;
    let report = store.gc(options(true))?;
    assert!(report.errors.is_empty(), "{:?}", report.errors);
    assert_eq!(report.entries[0].owners, 1);
    assert_eq!(report.entries[0].action, "held_owned");
    assert!(store.blob_path(&hash)?.exists());
    Ok(())
}
#[test]
fn orphan_handoff_dry_run_trash_restore_and_windowed_purge_are_consistent() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let db = project(&root.path().join("project"))?;
    let (path, hash, bytes) = source(root.path())?;
    drop(store.prepare(&db, &path, &hash, bytes)?);
    let before = serde_json::to_value(store.status()?)?;
    assert_eq!(store.gc(options(false))?.entries[0].action, "trash");
    assert_eq!(serde_json::to_value(store.status()?)?, before);
    let report = store.gc(options(true))?;
    assert!(report.errors.is_empty(), "{:?}", report.errors);
    assert_eq!(report.entries[0].state, "TRASHED");
    assert!(!store.blob_path(&hash)?.exists());
    let mut purge = options(true);
    purge.purge_expired = true;
    assert_eq!(store.gc(purge)?.entries[0].action, "held_trash_window");
    store.restore(&hash)?;
    assert!(store.blob_path(&hash)?.exists());
    store.gc(options(true))?;
    purge.trash_days = 0;
    assert_eq!(store.gc(purge)?.entries[0].action, "held_trash_window");
    store
        .connection()?
        .execute("UPDATE blobs SET purge_after_ms=0 WHERE sha256=?1", [&hash])?;
    assert_eq!(store.gc(purge)?.entries[0].state, "DELETED");
    assert!(!store.trash_path(&hash)?.exists());
    assert!(store.restore(&hash).is_err());
    Ok(())
}
#[test]
fn live_writer_and_reader_pins_prevent_collection_and_stale_rows_reconcile() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let db = project(&root.path().join("project"))?;
    let (path, hash, bytes) = source(root.path())?;
    let publication = store.prepare(&db, &path, &hash, bytes)?;
    assert_eq!(store.gc(options(true))?.entries[0].action, "held_busy");
    assert!(store.pin(&hash, "reader").is_err());
    drop(publication);
    let reader = store.pin(&hash, "reader")?;
    assert_eq!(store.status()?[0].reader_pins, 1);
    assert_eq!(store.gc(options(true))?.entries[0].action, "held_busy");
    drop(reader);
    store.connection()?.execute(
        "INSERT INTO pins(pin_id,sha256,owner,created_at) VALUES ('stale',?1,'crashed','now')",
        [&hash],
    )?;
    assert_eq!(store.gc(options(true))?.entries[0].action, "trash");
    assert_eq!(store.status()?[0].reader_pins, 0);
    Ok(())
}
#[test]
fn interrupted_move_and_restore_intents_resume_without_overwriting_either_copy() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let db = project(&root.path().join("project"))?;
    let (path, hash, bytes) = source(root.path())?;
    drop(store.prepare(&db, &path, &hash, bytes)?);
    store.connection()?.execute(
        "UPDATE blobs SET state='TRASHING',trashed_at_ms=1,purge_after_ms=9223372036854775807 WHERE sha256=?1",
        [&hash],
    )?;
    store.move_file(&store.blob_path(&hash)?, &store.trash_path(&hash)?)?;
    assert_eq!(
        store.gc(options(false))?.entries[0].action,
        "resume_pending"
    );
    assert_eq!(store.gc(options(true))?.entries[0].state, "TRASHED");
    store.connection()?.execute(
        "UPDATE blobs SET state='RESTORING' WHERE sha256=?1",
        [&hash],
    )?;
    store.move_file(&store.trash_path(&hash)?, &store.blob_path(&hash)?)?;
    store.restore(&hash)?;
    assert_eq!(store.status()?[0].state, "LIVE");
    Ok(())
}
#[test]
fn missing_project_and_legacy_payloads_stay_protected() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let db = project(&root.path().join("project"))?;
    let (path, hash, bytes) = source(root.path())?;
    let target = store.blob_path(&hash)?;
    store.checked_path(target.parent().unwrap().strip_prefix(store.root())?, true)?;
    fs::copy(&path, &target)?;
    drop(store.prepare(&db, &path, &hash, bytes)?);
    let report = store.gc(options(true))?;
    assert_eq!(report.entries[0].action, "held_legacy");
    // A second, newly cataloged blob has an owner whose project cannot be read.
    fs::write(&path, b"another payload")?;
    let hash = compute_sha256(&path)?;
    let bytes = file_size(&path)?;
    drop(store.prepare(&db, &path, &hash, bytes)?);
    fs::rename(db.path(), root.path().join("unavailable.sqlite3"))?;
    let report = store.gc(options(true))?;
    let entry = report
        .entries
        .iter()
        .find(|entry| entry.sha256 == hash)
        .unwrap();
    assert_eq!(entry.action, "held_owned");
    assert!(!entry.diagnostics.is_empty());
    assert!(store.blob_path(&hash)?.exists());
    Ok(())
}
#[test]
fn corrupt_payloads_are_not_reused_or_silently_collected() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let db = project(&root.path().join("project"))?;
    let (path, hash, bytes) = source(root.path())?;
    drop(store.prepare(&db, &path, &hash, bytes)?);
    fs::write(store.blob_path(&hash)?, b"corrupt")?;
    assert!(store.prepare(&db, &path, &hash, bytes).is_err());
    assert!(store.pin(&hash, "reader").is_err());
    assert_eq!(store.gc(options(true))?.entries[0].action, "held_error");
    assert!(store.blob_path(&hash)?.exists());
    Ok(())
}

#[test]
fn crash_worker() -> Result<()> {
    let Ok(root) = std::env::var("RNAA_SHARED_TEST_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let store = SharedStore::open(&root.join("store"))?;
    let db = Database::open(&root.join("project"))?;
    let (path, hash, bytes) = source(root)?;
    let phase = std::env::var("RNAA_SHARED_TEST_PHASE")?;
    let _writer;
    let _reader;
    if phase == "reader" {
        _reader = Some(store.pin(&hash, "crash-test")?);
        _writer = None;
    } else {
        let publication = store.prepare(&db, &path, &hash, bytes)?;
        if phase == "committed" {
            let blob = publication.record();
            db.record_shared_blob(blob)?;
            db.record_artifact(&ArtifactRecord {
                project_id: db.project_id().into(),
                run_accession: None,
                kind: ArtifactKind::QuantAbundanceH5,
                path: publication.owner.artifact_path.clone(),
                blob_id: Some(blob.blob_id.clone()),
                shared_path: Some(blob.storage_path.clone()),
                checksum_type: "sha256".into(),
                checksum: hash,
                bytes,
                created_at: now_rfc3339(),
            })?;
        }
        _writer = Some(publication);
        _reader = None;
    }
    fs::write(root.join("ready"), b"ready")?;
    loop {
        std::thread::sleep(Duration::from_secs(1));
    }
}

#[test]
fn killed_writers_reconcile_both_sides_of_commit_and_killed_readers_release_exclusion() -> Result<()>
{
    struct Child(std::process::Child);
    impl Drop for Child {
        fn drop(&mut self) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
    for phase in ["uncommitted", "committed", "reader"] {
        let root = tempfile::tempdir()?;
        let store = SharedStore::open(&root.path().join("store"))?;
        let db = project(&root.path().join("project"))?;
        let (path, hash, bytes) = source(root.path())?;
        if phase == "reader" {
            drop(store.prepare(&db, &path, &hash, bytes)?);
        }
        let mut child = Child(
            std::process::Command::new(std::env::current_exe()?)
                .args([
                    "--exact",
                    "shared_store::tests::crash_worker",
                    "--nocapture",
                ])
                .env("RNAA_SHARED_TEST_ROOT", root.path())
                .env("RNAA_SHARED_TEST_PHASE", phase)
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::piped())
                .spawn()?,
        );
        let started = std::time::Instant::now();
        while !root.path().join("ready").exists() {
            assert!(child.0.try_wait()?.is_none(), "worker exited in {phase}");
            assert!(started.elapsed() < Duration::from_secs(15));
            std::thread::sleep(Duration::from_millis(10));
        }
        assert_eq!(store.gc(options(true))?.entries[0].action, "held_busy");
        child.0.kill()?;
        child.0.wait()?;
        let report = store.gc(options(true))?;
        assert!(report.errors.is_empty(), "{:?}", report.errors);
        if phase == "committed" {
            assert_eq!(report.entries[0].action, "held_owned");
            assert_eq!(store.status()?[0].prepared_owners, 0);
            assert!(store.blob_path(&hash)?.exists());
        } else {
            assert_eq!(report.entries[0].action, "trash");
        }
        assert_eq!(store.status()?[0].reader_pins, 0);
    }
    Ok(())
}

#[cfg(unix)]
#[test]
fn catalog_and_payload_symlinks_cannot_redirect_store_operations() -> Result<()> {
    use std::os::unix::fs::symlink;
    let root = tempfile::tempdir()?;
    let external = root.path().join("external");
    fs::write(&external, b"untouched")?;
    let catalog_root = root.path().join("catalog-link");
    fs::create_dir(&catalog_root)?;
    symlink(&external, catalog_root.join("catalog.sqlite3"))?;
    assert!(SharedStore::open(&catalog_root).is_err());
    let store = SharedStore::open(&root.path().join("store"))?;
    let db = project(&root.path().join("project"))?;
    let (path, hash, bytes) = source(root.path())?;
    let target = store.blob_path(&hash)?;
    store.checked_path(target.parent().unwrap().strip_prefix(store.root())?, true)?;
    symlink(&external, &target)?;
    assert!(store.prepare(&db, &path, &hash, bytes).is_err());
    assert_eq!(fs::read(external)?, b"untouched");
    Ok(())
}

#[test]
fn simultaneous_store_initialization_and_publications_preserve_every_owner() -> Result<()> {
    let root = tempfile::tempdir()?;
    let barrier = std::sync::Arc::new(std::sync::Barrier::new(8));
    std::thread::scope(|scope| -> Result<()> {
        let mut workers = Vec::new();
        for index in 0..8 {
            let barrier = barrier.clone();
            let path = root.path().to_path_buf();
            workers.push(scope.spawn(move || -> Result<()> {
                let project_root = path.join(format!("project-{index}"));
                let db = project(&project_root)?;
                let (source, hash, bytes) = source(&project_root)?;
                barrier.wait();
                let store = SharedStore::open(&path.join("store"))?;
                let publication = store.prepare(&db, &source, &hash, bytes)?;
                register(&db, &publication)
            }));
        }
        for worker in workers {
            worker.join().unwrap()?;
        }
        Ok(())
    })?;
    let store = SharedStore::open(&root.path().join("store"))?;
    assert_eq!(store.status()?.len(), 1);
    assert_eq!(store.status()?[0].owners, 8);
    assert_eq!(store.status()?[0].prepared_owners, 0);
    Ok(())
}

#[test]
fn catalog_initialization_exclusion_covers_connection_and_journal_mode() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore {
        root: fs::canonicalize(root.path())?,
    };
    store.checked_path(Path::new(".locks"), true)?;
    let lock = store.initialization_lock()?;
    let (started_tx, started_rx) = std::sync::mpsc::channel();
    let (finished_tx, finished_rx) = std::sync::mpsc::channel();
    let path = root.path().to_path_buf();
    let worker = std::thread::spawn(move || {
        started_tx.send(()).unwrap();
        let result = SharedStore::open(&path);
        finished_tx.send(result).unwrap();
    });
    started_rx.recv()?;
    assert!(
        finished_rx
            .recv_timeout(Duration::from_millis(150))
            .is_err(),
        "opening must wait before touching SQLite"
    );
    assert!(!root.path().join("catalog.sqlite3").exists());
    drop(lock);
    let opened = finished_rx.recv_timeout(Duration::from_secs(5))??;
    worker.join().unwrap();
    assert!(opened.status()?.is_empty());
    Ok(())
}
