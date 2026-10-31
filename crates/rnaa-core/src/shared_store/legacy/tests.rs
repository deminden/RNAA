use super::*;
mod crashes;
mod helpers;
use helpers::*;

#[test]
fn schema_one_upgrade_preserves_rows_and_legacy_holds() -> Result<()> {
    let root = tempfile::tempdir()?;
    let path = root.path().join("catalog.sqlite3");
    let conn = Connection::open(&path)?;
    conn.execute_batch(include_str!("../schema.sql"))?;
    let sha = "a".repeat(64);
    conn.execute(
        "INSERT INTO blobs(sha256,bytes,state,legacy_hold,created_at) VALUES (?1,3,'LIVE',1,'old')",
        [&sha],
    )?;
    drop(conn);
    let store = SharedStore::open(root.path())?;
    let schema: i64 = store
        .connection()?
        .pragma_query_value(None, "user_version", |row| row.get(0))?;
    assert_eq!(schema, 2);
    assert!(store.status()?[0].legacy_hold);
    assert_eq!(store.status()?[0].bytes, 3);
    drop(store);
    let store = SharedStore::open(root.path())?;
    assert!(store.status()?[0].legacy_hold);
    Ok(())
}
#[test]
fn read_only_old_catalog_import_preserves_bytes_and_only_explicit_inventory_releases_holds()
-> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let first = legacy_project(&root.path().join("one"), &store, b"old shared data", false)?;
    let second = legacy_project(&root.path().join("two"), &store, b"old shared data", false)?;
    let inventory = vec![first.root.clone(), second.root.clone()];
    let before = compute_sha256(&first.database)?;
    let preview = store.import_legacy(&inventory, opts(false, false))?;
    assert_eq!(preview.entries[0].owners, 2);
    assert!(!preview.complete);
    assert!(store.status()?.is_empty());
    store.import_legacy(&inventory, opts(true, false))?;
    assert!(store.status()?[0].legacy_hold);
    assert_eq!(store.status()?[0].owners, 2);
    let applied = store.import_legacy(&inventory, opts(true, true))?;
    #[cfg(unix)]
    {
        assert!(applied.complete, "{:?}", applied.errors);
        assert!(!store.status()?[0].legacy_hold);
    }
    assert!(applied.errors.is_empty(), "{:?}", applied.errors);
    assert_eq!(compute_sha256(&first.database)?, before);
    let conn =
        Connection::open_with_flags(&first.database, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    assert!(!projects::table(&conn, "artifact_generations")?);
    drop(conn);
    remove_registration(&first)?;
    assert_eq!(store.gc(gc(false))?.entries[0].action, "held_owned");
    remove_registration(&second)?;
    #[cfg(unix)]
    {
        let report = store.gc(gc(true))?;
        assert!(report.errors.is_empty(), "{:?}", report.errors);
        assert_eq!(report.entries[0].state, "TRASHED");
        store.restore(&first.sha)?;
        assert_eq!(fs::read(store.blob_path(&first.sha)?)?, b"old shared data");
    }
    assert_eq!(fs::read(first.source)?, b"old shared data");
    assert_eq!(fs::read(second.source)?, b"old shared data");
    Ok(())
}
#[cfg(unix)]
#[test]
fn unexplained_links_hold_migration_and_new_links_block_orphan_purge_after_import() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let first = legacy_project(
        &root.path().join("one"),
        &store,
        b"old hardlinked data",
        true,
    )?;
    let unknown = root.path().join("unregistered-link");
    fs::hard_link(&first.source, &unknown)?;
    let paths = vec![first.root.clone()];
    let report = store.import_legacy(&paths, opts(true, true))?;
    assert!(!report.complete);
    assert!(store.status()?[0].legacy_hold);
    assert!(
        report.entries[0]
            .reasons
            .contains(&"unexplained_or_unavailable_hardlinks".into())
    );
    fs::remove_file(&unknown)?;
    assert!(store.import_legacy(&paths, opts(true, true))?.complete);
    remove_registration(&first)?;
    assert_eq!(store.gc(gc(true))?.entries[0].action, "held_external_links");
    fs::remove_file(&first.source)?;
    let mut options = gc(true);
    options.trash_days = 0;
    store.gc(options)?;
    let trash = store.trash_path(&first.sha)?;
    fs::hard_link(&trash, &unknown)?;
    options.purge_expired = true;
    assert_eq!(store.gc(options)?.entries[0].action, "held_external_links");
    store.connection()?.execute(
        "UPDATE blobs SET state='DELETING' WHERE sha256=?1",
        [&first.sha],
    )?;
    assert_eq!(store.gc(options)?.entries[0].action, "held_error");
    assert!(trash.exists());
    fs::remove_file(unknown)?;
    assert_eq!(store.gc(options)?.entries[0].state, "DELETED");
    let marker: Option<String> = store.connection()?.query_row(
        "SELECT legacy_import_id FROM blobs WHERE sha256=?1",
        [&first.sha],
        |row| row.get(0),
    )?;
    assert!(marker.is_none());
    Ok(())
}
#[test]
fn invalid_project_payload_and_busy_pins_never_grant_new_legacy_authority() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let project = legacy_project(
        &root.path().join("project"),
        &store,
        b"original shared data",
        false,
    )?;
    let lock = crate::generations::pins::read_lock(&project.database)?;
    assert!(
        !store
            .import_legacy(std::slice::from_ref(&project.root), opts(true, true))?
            .busy_projects
            .is_empty()
    );
    drop(lock);
    fs::write(&project.source, b"changed local view")?;
    let report = store.import_legacy(std::slice::from_ref(&project.root), opts(true, true))?;
    assert!(!report.errors.is_empty());
    assert!(!report.complete);
    assert!(store.status()?[0].legacy_hold);
    assert_eq!(store.status()?[0].owners, 0);
    assert_eq!(
        fs::read(store.blob_path(&project.sha)?)?,
        b"original shared data"
    );
    fs::copy(store.blob_path(&project.sha)?, &project.source)?;
    let reader = store.pin(&project.sha, "migration test")?;
    let report = store.import_legacy(std::slice::from_ref(&project.root), opts(true, true))?;
    assert_eq!(report.entries[0].action, "held_busy");
    drop(reader);
    assert!(
        store
            .import_legacy(std::slice::from_ref(&project.root), opts(true, true))?
            .errors
            .is_empty()
    );
    Ok(())
}
#[cfg(unix)]
#[test]
fn missing_owner_index_or_witness_keeps_imported_bytes_protected() -> Result<()> {
    for fault in ["owner", "index", "witness"] {
        let root = tempfile::tempdir()?;
        let store = SharedStore::open(&root.path().join("store"))?;
        let project = legacy_project(
            &root.path().join("project"),
            &store,
            b"preserve migrated bytes",
            false,
        )?;
        assert!(
            store
                .import_legacy(std::slice::from_ref(&project.root), opts(true, true))?
                .complete
        );
        match fault {
            "owner" => {
                store.connection()?.execute("DELETE FROM owners", [])?;
            }
            "index" => {
                store
                    .connection()?
                    .execute("DELETE FROM legacy_blob_imports", [])?;
            }
            "witness" => {
                let id: String = store.connection()?.query_row(
                    "SELECT import_id FROM legacy_imports",
                    [],
                    |row| row.get(0),
                )?;
                fs::write(store.root.join(format!(".imports/{id}.json")), b"{}")?;
            }
            _ => unreachable!(),
        }
        let report = store.gc(gc(true))?;
        assert_eq!(report.entries[0].action, "held_error");
        assert!(!report.errors.is_empty());
        assert!(store.blob_path(&project.sha)?.exists());
    }
    Ok(())
}
#[test]
fn failed_catalog_transaction_publishes_no_partial_import_and_retry_never_adopts_orphan_witnesses()
-> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let project = legacy_project(
        &root.path().join("project"),
        &store,
        b"atomic imported data",
        false,
    )?;
    store.connection()?.execute_batch("CREATE TRIGGER reject_legacy_owner BEFORE INSERT ON owners BEGIN SELECT RAISE(ABORT,'injected'); END;")?;
    assert!(
        store
            .import_legacy(std::slice::from_ref(&project.root), opts(true, true))
            .is_err()
    );
    assert!(store.status()?.is_empty());
    let count: i64 =
        store
            .connection()?
            .query_row("SELECT count(*) FROM legacy_imports", [], |row| row.get(0))?;
    assert_eq!(count, 0);
    let witnesses = fs::read_dir(store.root.join(".imports"))?.count();
    assert!(witnesses > 0);
    store
        .connection()?
        .execute_batch("DROP TRIGGER reject_legacy_owner;")?;
    assert!(
        store
            .import_legacy(std::slice::from_ref(&project.root), opts(true, true))?
            .errors
            .is_empty()
    );
    assert!(fs::read_dir(store.root.join(".imports"))?.count() > witnesses);
    assert_eq!(store.status()?[0].owners, 1);
    Ok(())
}

#[cfg(unix)]
#[test]
fn partial_inventory_cannot_erase_prior_proof_or_release_unindexed_blobs_for_known_missing_projects()
-> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let first = legacy_project(&root.path().join("one"), &store, b"shared inventory", false)?;
    let second = legacy_project(&root.path().join("two"), &store, b"shared inventory", false)?;
    let paths = vec![first.root.clone(), second.root.clone()];
    assert!(store.import_legacy(&paths, opts(true, true))?.complete);
    let before: String = store.connection()?.query_row(
        "SELECT legacy_import_id FROM blobs WHERE sha256=?1",
        [&first.sha],
        |row| row.get(0),
    )?;
    let orphan = digest(b"unindexed orphan");
    let path = store.blob_path(&orphan)?;
    store.checked_path(path.parent().unwrap().strip_prefix(store.root())?, true)?;
    fs::write(&path, b"unindexed orphan")?;
    let report = store.import_legacy(std::slice::from_ref(&first.root), opts(true, true))?;
    assert!(!report.complete);
    assert_eq!(
        report
            .entries
            .iter()
            .find(|row| row.sha256 == first.sha)
            .unwrap()
            .action,
        "held_existing_import"
    );
    assert!(
        store
            .status()?
            .iter()
            .find(|blob| blob.sha256 == orphan)
            .unwrap()
            .legacy_hold
    );
    let after: String = store.connection()?.query_row(
        "SELECT legacy_import_id FROM blobs WHERE sha256=?1",
        [&first.sha],
        |row| row.get(0),
    )?;
    assert_eq!(before, after);
    assert!(store.import_legacy(&paths, opts(true, true))?.complete);
    Ok(())
}
#[cfg(unix)]
#[test]
fn unsupported_project_schema_cannot_become_owner_absence_after_import() -> Result<()> {
    let root = tempfile::tempdir()?;
    let store = SharedStore::open(&root.path().join("store"))?;
    let project = legacy_project(
        &root.path().join("project"),
        &store,
        b"future schema protection",
        false,
    )?;
    assert!(
        store
            .import_legacy(std::slice::from_ref(&project.root), opts(true, true))?
            .complete
    );
    remove_registration(&project)?;
    Connection::open(&project.database)?.execute(
        "INSERT INTO schema_migrations(version,applied_at) VALUES ('9999_unknown','future')",
        [],
    )?;
    let report = store.gc(gc(true))?;
    assert_eq!(report.entries[0].action, "held_error");
    assert!(store.blob_path(&project.sha)?.exists());
    Ok(())
}
