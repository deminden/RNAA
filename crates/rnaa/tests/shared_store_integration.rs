//! Native result publication and operator contracts for a shared store.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/recovery_fixture.rs"]
mod fixture;
use commands::{command, success};
use fixture::{head, setup};
use rnaa_core::{Database, ProjectConfig, generations::OutputGroup, shared_store::SharedStore};
use std::{fs, path::Path};

fn configure(root: &Path, shared: &Path) {
    let path = root.join("rnaa.toml");
    let mut config = ProjectConfig::load(&path).unwrap();
    config.storage.shared_root = shared.display().to_string();
    config.save(&path).unwrap();
    Database::open(root)
        .unwrap()
        .set_project_config(&config)
        .unwrap();
}
fn storage(root: &Path, args: &[&str]) -> serde_json::Value {
    let mut command_args = vec!["storage"];
    command_args.extend(args);
    let result = command(root, &command_args);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    serde_json::from_slice(&result.stdout).unwrap()
}

#[test]
fn native_publications_finalize_store_ownership_and_historical_results_survive_gc() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let shared = root.join("shared");
    configure(root, &shared);
    success(root, &["quant", "--single-fragment-length", "180"]);
    let original = head(&db, OutputGroup::Quant);
    let first = storage(root, &["status"]);
    assert!(!first["blobs"].as_array().unwrap().is_empty());
    assert!(
        first["blobs"]
            .as_array()
            .unwrap()
            .iter()
            .all(|blob| blob["prepared_owners"] == 0 && blob["owners"].as_u64().unwrap() > 0)
    );
    success(root, &["quant", "--single-fragment-length", "170"]);
    assert_ne!(head(&db, OutputGroup::Quant), original);
    let report = storage(root, &["gc", "--apply"]);
    assert!(
        report["entries"]
            .as_array()
            .unwrap()
            .iter()
            .all(|entry| entry["action"] == "held_owned")
    );
    let catalog = rusqlite::Connection::open(shared.join("catalog.sqlite3")).unwrap();
    let protected: i64 = catalog
        .query_row(
            "SELECT count(*) FROM owners WHERE generation_id=?1",
            [&original.generation_id],
            |row| row.get(0),
        )
        .unwrap();
    assert!(
        protected > 0,
        "retired active registrations must retain historical ownership"
    );
    assert!(db.verify_generation_head(&original).unwrap());
    // Missing committed ownership metadata is a fault, not an unused-blob proof.
    rusqlite::Connection::open(db.path())
        .unwrap()
        .execute(
            "UPDATE artifact_generations SET artifacts_json='[]' WHERE generation_id=?1",
            [&original.generation_id],
        )
        .unwrap();
    let held = storage(root, &["gc", "--apply"]);
    assert!(
        held["entries"]
            .as_array()
            .unwrap()
            .iter()
            .any(|entry| entry["action"] == "held_owned"
                && !entry["diagnostics"].as_array().unwrap().is_empty())
    );
}

#[test]
fn failed_project_publication_leaves_a_reconcilable_pin_and_reversible_orphan() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let shared = root.join("shared");
    configure(root, &shared);
    let original = head(&db, OutputGroup::Quant);
    rusqlite::Connection::open(db.path()).unwrap().execute_batch(
        "CREATE TRIGGER reject_shared_publication BEFORE INSERT ON artifacts
         WHEN NEW.kind='ABUNDANCE_H5' BEGIN SELECT RAISE(ABORT,'injected registration failure'); END;"
    ).unwrap();
    assert!(
        !command(root, &["quant", "--single-fragment-length", "180"])
            .status
            .success()
    );
    assert_eq!(head(&db, OutputGroup::Quant), original);
    let store = SharedStore::open(&shared).unwrap();
    let before = store.status().unwrap();
    assert_eq!(before.len(), 1);
    assert_eq!(before[0].prepared_owners, 1);
    let planned = storage(root, &["gc"]);
    assert_eq!(planned["applied"], false);
    assert_eq!(planned["entries"][0]["action"], "trash");
    assert_eq!(store.status().unwrap()[0].state, "LIVE");
    let applied = storage(root, &["gc", "--apply"]);
    assert_eq!(applied["entries"][0]["state"], "TRASHED");
    let hash = &before[0].sha256;
    storage(root, &["restore", &format!("sha256:{hash}")]);
    let reader = store.pin(hash, "integration-reader").unwrap();
    assert!(Path::new(&reader.record().storage_path).exists());
    assert_eq!(
        storage(root, &["gc", "--apply"])["entries"][0]["action"],
        "held_busy"
    );
    drop(reader);
    rusqlite::Connection::open(db.path())
        .unwrap()
        .execute_batch("DROP TRIGGER reject_shared_publication")
        .unwrap();
    success(root, &["quant"]);
    assert_ne!(head(&db, OutputGroup::Quant), original);
    assert!(
        store
            .status()
            .unwrap()
            .iter()
            .all(|blob| blob.prepared_owners == 0)
    );
    assert!(db.verify_generation_head(&original).unwrap());
}

#[test]
fn concurrent_native_projects_register_independent_owners_in_one_catalog() {
    let temp = tempfile::tempdir().unwrap();
    let shared = temp.path().join("shared");
    let first = temp.path().join("one");
    let second = temp.path().join("two");
    fs::create_dir(&first).unwrap();
    fs::create_dir(&second).unwrap();
    setup(&first, None, "both");
    setup(&second, None, "both");
    configure(&first, &shared);
    configure(&second, &shared);
    std::thread::scope(|scope| {
        scope.spawn(|| success(&first, &["quant", "--single-fragment-length", "180"]));
        scope.spawn(|| success(&second, &["quant", "--single-fragment-length", "180"]));
    });
    let store = SharedStore::open(&shared).unwrap();
    let blobs = store.status().unwrap();
    assert_eq!(blobs.iter().map(|blob| blob.owners).sum::<u64>(), 2);
    assert!(blobs.iter().all(|blob| blob.prepared_owners == 0));
    let catalog = rusqlite::Connection::open(shared.join("catalog.sqlite3")).unwrap();
    assert_eq!(
        catalog
            .query_row("SELECT count(DISTINCT project_id) FROM owners", [], |row| {
                row.get::<_, i64>(0)
            })
            .unwrap(),
        2
    );
    assert!(
        storage(&first, &["gc", "--apply"])["entries"]
            .as_array()
            .unwrap()
            .iter()
            .all(|entry| entry["action"] == "held_owned")
    );
}

#[test]
fn relative_store_configuration_follows_project_root_and_explicit_inspection_needs_no_project() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path().join("project");
    let caller = temp.path().join("caller");
    fs::create_dir(&root).unwrap();
    fs::create_dir(&caller).unwrap();
    setup(&root, None, "both");
    configure(&root, Path::new("shared"));
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_rnaa"))
        .arg("--root")
        .arg(&root)
        .args(["quant", "--single-fragment-length", "180"])
        .current_dir(&caller)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(root.join("shared/catalog.sqlite3").is_file());
    assert!(!caller.join("shared").exists());
    let result = command(
        &caller,
        &[
            "storage",
            "--shared-root",
            root.join("shared").to_str().unwrap(),
            "status",
        ],
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let status: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(status["blobs"][0]["owners"], 1);
}
