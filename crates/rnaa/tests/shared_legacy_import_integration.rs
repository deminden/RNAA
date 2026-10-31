//! Importing shared ownership must preserve real native history and old catalogs.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/recovery_fixture.rs"]
mod fixture;
use commands::{command, success};
use fixture::{Source, head, setup};
use rnaa_core::{Database, ProjectConfig, generations::OutputGroup, shared_store::SharedStore};
use rusqlite::{Connection, params};
use std::{fs, path::Path};
fn shared(root: &Path) {
    let path = root.join("rnaa.toml");
    let mut config = ProjectConfig::load(&path).unwrap();
    config.storage.shared_root = root.join("shared").display().to_string();
    config.save(&path).unwrap();
    Database::open(root)
        .unwrap()
        .set_project_config(&config)
        .unwrap();
}
fn import(root: &Path, args: &[&str]) -> serde_json::Value {
    let mut full = vec![
        "storage",
        "import-legacy",
        "--project",
        root.to_str().unwrap(),
    ];
    full.extend(args);
    let output = command(root, &full);
    assert!(
        output.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}
#[test]
fn native_history_import_keeps_original_seals_and_ownership_until_project_purge() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let source = Source::new();
    let (_paths, db) = setup(root, Some(&source), "both");
    shared(root);
    success(root, &["quant", "--single-fragment-length", "180"]);
    let old = head(&db, OutputGroup::Quant);
    success(root, &["quant", "--single-fragment-length", "170"]);
    let current = head(&db, OutputGroup::Quant);
    let original = fs::read(old.payload_path.parent().unwrap().join("seal.json")).unwrap();
    let conn = Connection::open(root.join("shared/catalog.sqlite3")).unwrap();
    conn.execute("DELETE FROM owners", []).unwrap();
    conn.execute("UPDATE blobs SET legacy_hold=1", []).unwrap();
    drop(conn);
    let before = SharedStore::open(&root.join("shared"))
        .unwrap()
        .status()
        .unwrap();
    assert!(
        before
            .iter()
            .all(|blob| blob.owners == 0 && blob.legacy_hold)
    );
    let preview = import(root, &["--complete-project-inventory"]);
    assert_eq!(preview["applied"], false);
    assert_eq!(preview["errors"], serde_json::json!([]));
    assert!(
        SharedStore::open(&root.join("shared"))
            .unwrap()
            .status()
            .unwrap()
            .iter()
            .all(|blob| blob.legacy_hold)
    );
    let applied = import(root, &["--apply", "--complete-project-inventory"]);
    assert_eq!(applied["errors"], serde_json::json!([]));
    #[cfg(unix)]
    assert_eq!(applied["complete"], true);
    assert_eq!(
        fs::read(old.payload_path.parent().unwrap().join("seal.json")).unwrap(),
        original
    );
    assert_eq!(head(&db, OutputGroup::Quant), current);
    assert!(db.verify_generation_head(&current).unwrap());
    let statuses = SharedStore::open(&root.join("shared"))
        .unwrap()
        .status()
        .unwrap();
    assert!(statuses.iter().all(|blob| blob.owners > 0));
    success(root, &["storage", "gc", "--apply"]);
    assert!(
        SharedStore::open(&root.join("shared"))
            .unwrap()
            .status()
            .unwrap()
            .iter()
            .all(|blob| blob.state == "LIVE")
    );
    success(
        root,
        &[
            "storage",
            "generations",
            "gc",
            "--apply",
            "--trash-days",
            "0",
        ],
    );
    success(
        root,
        &["storage", "generations", "gc", "--apply", "--purge-expired"],
    );
    success(root, &["storage", "gc", "--apply"]);
    #[cfg(unix)]
    assert!(
        SharedStore::open(&root.join("shared"))
            .unwrap()
            .status()
            .unwrap()
            .iter()
            .any(|blob| blob.state == "TRASHED")
    );
    assert_eq!(head(&db, OutputGroup::Quant), current);
    assert!(db.verify_generation_head(&current).unwrap());
    success(root, &["quant", "--single-fragment-length", "170"]);
    assert_eq!(head(&db, OutputGroup::Quant), current);
}
#[test]
fn explicit_store_cli_reads_an_old_project_without_upgrading_it_and_reports_busy_projects() {
    let root = tempfile::tempdir().unwrap();
    let root = fs::canonicalize(root.path()).unwrap();
    let legacy = root.join("legacy");
    fs::create_dir(&legacy).unwrap();
    let database = legacy.join("state.sqlite");
    let conn = Connection::open(&database).unwrap();
    conn.execute_batch(include_str!("../../rnaa-core/migrations/0001_init.sql"))
        .unwrap();
    conn.execute_batch(include_str!(
        "../../rnaa-core/migrations/0002_shared_store.sql"
    ))
    .unwrap();
    let id = "legacy-shared-project";
    conn.execute_batch(
        "INSERT INTO schema_migrations VALUES ('0001_init','old'),('0002_shared_store','old');",
    )
    .unwrap();
    conn.execute("INSERT INTO projects(project_id,root_dir,created_at,config_json) VALUES (?1,?2,'old','{}')",params![id,legacy.display().to_string()]).unwrap();
    drop(conn);
    let old_settings = "[obsolete_runtime]\nrscript = 'legacy-only-settings'\n";
    fs::write(legacy.join("rnaa.toml"), old_settings).unwrap();
    let before = rnaa_core::util::compute_sha256(&database).unwrap();
    let output = command(
        &legacy,
        &[
            "storage",
            "--shared-root",
            root.join("shared").to_str().unwrap(),
            "import-legacy",
            "--project",
            legacy.to_str().unwrap(),
            "--apply",
            "--complete-project-inventory",
        ],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["complete"], true);
    assert_eq!(rnaa_core::util::compute_sha256(&database).unwrap(), before);
    assert_eq!(
        fs::read_to_string(legacy.join("rnaa.toml")).unwrap(),
        old_settings
    );
    success(&root, &["init"]);
    let db = Database::open(&root).unwrap();
    let _pin =
        rnaa_core::generations::GenerationStorePin::acquire(db.path(), "shared migration CLI test")
            .unwrap();
    let output = command(
        &root,
        &[
            "storage",
            "--shared-root",
            root.join("shared").to_str().unwrap(),
            "import-legacy",
            "--project",
            root.to_str().unwrap(),
            "--apply",
        ],
    );
    assert!(!output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(!report["busy_projects"].as_array().unwrap().is_empty());
}
