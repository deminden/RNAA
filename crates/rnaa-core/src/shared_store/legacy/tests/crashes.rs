use super::*;
#[test]
fn open_crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_SHARED_OPEN_CRASH_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    fs::write(root.join("open-attempt"), b"ready")?;
    Database::open(&root.join("project"))?;
    fs::write(root.join("open-done"), b"done")?;
    loop {
        std::thread::sleep(Duration::from_secs(1));
    }
}
#[test]
fn normal_project_open_cannot_migrate_or_backfill_during_legacy_maintenance() -> Result<()> {
    struct KillOnDrop(std::process::Child);
    impl Drop for KillOnDrop {
        fn drop(&mut self) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let project = legacy_project(
        &root.path().join("project"),
        &store,
        b"unchanged legacy payload",
        false,
    )?;
    let before = compute_sha256(&project.database)?;
    let lock = crate::generations::pins::maintenance_lock_path(&project.database)?.unwrap();
    let mut child = KillOnDrop(
        std::process::Command::new(std::env::current_exe()?)
            .args([
                "--exact",
                "shared_store::legacy::tests::crashes::open_crash_worker",
                "--nocapture",
            ])
            .env("RNAA_SHARED_OPEN_CRASH_ROOT", root.path())
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::piped())
            .spawn()?,
    );
    let started = std::time::Instant::now();
    while !root.path().join("open-attempt").exists() {
        assert!(child.0.try_wait()?.is_none());
        assert!(started.elapsed() < Duration::from_secs(15));
        std::thread::sleep(Duration::from_millis(10));
    }
    std::thread::sleep(Duration::from_millis(200));
    assert!(child.0.try_wait()?.is_none());
    assert!(!root.path().join("open-done").exists());
    assert_eq!(compute_sha256(&project.database)?, before);
    drop(lock);
    while !root.path().join("open-done").exists() {
        assert!(child.0.try_wait()?.is_none());
        assert!(started.elapsed() < Duration::from_secs(15));
        std::thread::sleep(Duration::from_millis(10));
    }
    child.0.kill()?;
    child.0.wait()?;
    let conn = Connection::open_with_flags(
        &project.database,
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    )?;
    assert!(projects::table(&conn, "artifact_generations")?);
    assert_eq!(fs::read(project.source)?, b"unchanged legacy payload");
    Ok(())
}
#[test]
fn import_crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_SHARED_IMPORT_CRASH_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let store = SharedStore::open(&root.join("store"))?;
    fs::write(root.join("import-opened"), b"ready")?;
    while !root.join("import-start").exists() {
        std::thread::sleep(Duration::from_millis(10));
    }
    store.import_legacy(&[root.join("project")], opts(true, true))?;
    fs::write(root.join("import-ready"), b"done")?;
    loop {
        std::thread::sleep(Duration::from_secs(1));
    }
}
#[test]
fn real_killed_importers_leave_unpublished_witnesses_or_one_atomic_owned_catalog() -> Result<()> {
    struct KillOnDrop(std::process::Child);
    impl Drop for KillOnDrop {
        fn drop(&mut self) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
    fn wait(child: &mut KillOnDrop, predicate: impl Fn() -> bool) -> Result<()> {
        let started = std::time::Instant::now();
        while !predicate() {
            assert!(
                child.0.try_wait()?.is_none(),
                "shared import child exited before readiness"
            );
            assert!(started.elapsed() < Duration::from_secs(15));
            std::thread::sleep(Duration::from_millis(10));
        }
        Ok(())
    }
    for committed in [false, true] {
        let root = tempfile::tempdir()?;
        let store = SharedStore::open(&root.path().join("store"))?;
        let project = legacy_project(
            &root.path().join("project"),
            &store,
            b"crash-safe legacy data",
            false,
        )?;
        let mut child = KillOnDrop(
            std::process::Command::new(std::env::current_exe()?)
                .args([
                    "--exact",
                    "shared_store::legacy::tests::crashes::import_crash_worker",
                    "--nocapture",
                ])
                .env("RNAA_SHARED_IMPORT_CRASH_ROOT", root.path())
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::piped())
                .spawn()?,
        );
        wait(&mut child, || root.path().join("import-opened").exists())?;
        let writer = store.connection()?;
        if !committed {
            writer.execute_batch("BEGIN IMMEDIATE")?;
        }
        fs::write(root.path().join("import-start"), b"go")?;
        if committed {
            wait(&mut child, || root.path().join("import-ready").exists())?;
        } else {
            wait(&mut child, || {
                fs::read_dir(store.root.join(".imports"))
                    .unwrap()
                    .any(|entry| {
                        entry
                            .unwrap()
                            .path()
                            .extension()
                            .is_some_and(|extension| extension == "json")
                    })
            })?;
        }
        child.0.kill()?;
        child.0.wait()?;
        drop(writer);
        if !committed {
            assert!(store.status()?.is_empty());
            assert!(
                store
                    .import_legacy(std::slice::from_ref(&project.root), opts(true, true))?
                    .errors
                    .is_empty()
            );
        }
        assert_eq!(store.status()?[0].owners, 1);
        assert!(store.gc(gc(true))?.errors.is_empty());
        assert_eq!(fs::read(&project.source)?, b"crash-safe legacy data");
    }
    Ok(())
}
