//! Dependency-aware generation retirement through the native CLI.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/recovery_fixture.rs"]
mod fixture;
use commands::{command, success};
use fixture::{Source, head, replace, setup};
use rnaa_core::{
    Database, ProjectConfig,
    generations::{GenerationHead, OutputGroup},
};
use std::{fs, path::Path, sync::atomic::Ordering};

fn gc(root: &Path, args: &[&str]) -> serde_json::Value {
    let mut full = vec!["storage", "generations", "gc"];
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
fn entry<'a>(value: &'a serde_json::Value, head: &GenerationHead) -> &'a serde_json::Value {
    value["entries"]
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| entry["generation_id"] == head.generation_id)
        .unwrap()
}
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

#[test]
fn preview_keep_trash_restore_and_reuse_preserve_current_native_results() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    shared(root);
    success(root, &["quant", "--single-fragment-length", "180"]);
    let old = head(&db, OutputGroup::Quant);
    let seal = fs::read(old.payload_path.parent().unwrap().join("seal.json")).unwrap();
    success(root, &["quant", "--single-fragment-length", "170"]);
    let current = head(&db, OutputGroup::Quant);
    success(
        root,
        &[
            "storage",
            "generations",
            "keep",
            &old.generation_id,
            "--label",
            "paper",
        ],
    );
    assert_eq!(entry(&gc(root, &[]), &old)["action"], "held_reachable");
    success(
        root,
        &["storage", "generations", "unkeep", &old.generation_id],
    );
    let preview = gc(root, &[]);
    assert_eq!(entry(&preview, &old)["action"], "trash");
    assert!(old.payload_path.exists());
    let applied = gc(root, &["--apply"]);
    assert_eq!(entry(&applied, &old)["state"], "TRASHED");
    assert_eq!(head(&db, OutputGroup::Quant), current);
    assert!(db.verify_generation_head(&current).unwrap());
    assert!(!old.payload_path.exists());
    let listing = command(
        root,
        &[
            "lineage",
            "--group",
            "quant",
            "--subject",
            "RUN1",
            "--all-requests",
            "--verify",
        ],
    );
    assert!(
        listing.status.success(),
        "{}",
        String::from_utf8_lossy(&listing.stderr)
    );
    let history: serde_json::Value = serde_json::from_slice(&listing.stdout).unwrap();
    assert!(
        !history["generations"]
            .as_array()
            .unwrap()
            .iter()
            .any(|entry| entry["lineage"]["head"]["generation_id"] == old.generation_id)
    );
    success(
        root,
        &["storage", "generations", "restore", &old.generation_id],
    );
    assert_eq!(
        fs::read(old.payload_path.parent().unwrap().join("seal.json")).unwrap(),
        seal
    );
    assert_eq!(head(&db, OutputGroup::Quant), current);
    success(root, &["quant", "--single-fragment-length", "180"]);
    assert_eq!(head(&db, OutputGroup::Quant), old);
}

#[test]
fn completed_project_purge_authorizes_shared_blob_owner_release_only_after_trash_window() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    shared(root);
    success(root, &["quant", "--single-fragment-length", "180"]);
    let old = head(&db, OutputGroup::Quant);
    let old_hash = db
        .list_artifacts_for_run(Some("RUN1"))
        .unwrap()
        .into_iter()
        .find(|record| record.kind == rnaa_core::ArtifactKind::QuantAbundanceH5)
        .unwrap()
        .checksum;
    success(root, &["quant", "--single-fragment-length", "170"]);
    let current = head(&db, OutputGroup::Quant);
    let trashed = gc(root, &["--apply", "--trash-days", "0"]);
    assert_eq!(entry(&trashed, &old)["state"], "TRASHED");
    let before = command(root, &["storage", "gc", "--apply"]);
    assert!(before.status.success());
    let before: serde_json::Value = serde_json::from_slice(&before.stdout).unwrap();
    assert_eq!(
        before["entries"]
            .as_array()
            .unwrap()
            .iter()
            .find(|entry| entry["sha256"] == old_hash)
            .unwrap()["action"],
        "held_owned"
    );
    let purged = gc(root, &["--apply", "--purge-expired"]);
    assert_eq!(entry(&purged, &old)["state"], "PURGED");
    assert!(
        !command(
            root,
            &["storage", "generations", "restore", &old.generation_id]
        )
        .status
        .success()
    );
    let after = command(root, &["storage", "gc", "--apply"]);
    assert!(
        after.status.success(),
        "{}",
        String::from_utf8_lossy(&after.stderr)
    );
    let after: serde_json::Value = serde_json::from_slice(&after.stdout).unwrap();
    assert_eq!(
        after["entries"]
            .as_array()
            .unwrap()
            .iter()
            .find(|entry| entry["sha256"] == old_hash)
            .unwrap()["state"],
        "TRASHED"
    );
    assert_eq!(head(&db, OutputGroup::Quant), current);
    assert!(db.verify_generation_head(&current).unwrap());
}

#[test]
fn historical_input_bindings_and_pruning_proofs_survive_retiring_an_unneeded_fit() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "none");
    let original = head(&db, OutputGroup::Quant);
    replace(root, "both");
    let unneeded = head(&db, OutputGroup::Quant);
    let current_processing = head(&db, OutputGroup::Preprocess);
    let current_download = head(&db, OutputGroup::Download);
    source.failed.store(true, Ordering::SeqCst);
    success(
        root,
        &[
            "quant",
            "--single-fragment-length",
            "200",
            "--fastq-retention",
            "none",
        ],
    );
    assert_eq!(head(&db, OutputGroup::Quant), original);
    let report = gc(root, &["--apply"]);
    assert_eq!(entry(&report, &unneeded)["state"], "TRASHED");
    assert_eq!(
        entry(&report, &current_processing)["action"],
        "held_reachable"
    );
    assert_eq!(
        entry(&report, &current_download)["action"],
        "held_reachable"
    );
    let binding = db
        .generation_input_binding(OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    db.verify_input_binding(&binding).unwrap();
    assert!(db.verify_generation_head(&original).unwrap());
    let requests = source.requests.load(Ordering::SeqCst);
    success(root, &["quant"]);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert_eq!(head(&db, OutputGroup::Quant), original);
}
