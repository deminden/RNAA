//! Unpublished diagnostic cleanup must preserve real current scientific results.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/recovery_fixture.rs"]
mod fixture;
use commands::{command, success};
use fixture::{Source, head, setup};
use rnaa_core::{
    Database,
    generations::{Generation, GenerationHead, OutputGroup},
    jobs::{JobLease, JobOutcome},
};
use std::{fs, path::Path};

fn failed_generation(db: &Database) -> GenerationHead {
    let mut lease = JobLease::try_acquire(db, "deseq2", None, &[])
        .unwrap()
        .unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64)).unwrap();
    let mut generation = Generation::begin(&owned, OutputGroup::De, "deseq2").unwrap();
    let path = generation.head().payload_path.join("partial.tsv");
    fs::write(&path, b"partial").unwrap();
    generation.seal(&[path], &"a".repeat(64)).unwrap();
    let head = generation.head().clone();
    drop(generation);
    lease
        .finish(
            JobOutcome::Failed,
            Some("native test failure after sealing"),
        )
        .unwrap();
    drop(lease);
    head
}
fn gc(root: &Path, args: &[&str]) -> serde_json::Value {
    let mut full = vec!["storage", "generations", "unpublished", "gc"];
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
        .find(|row| row["generation_id"] == head.generation_id)
        .unwrap()
}
#[test]
fn native_preview_trash_restore_and_purge_never_publish_unfinished_results_or_damage_counts() {
    let root = tempfile::tempdir().unwrap();
    let root = root.path();
    let source = Source::new();
    let (_paths, db) = setup(root, Some(&source), "both");
    let current = head(&db, OutputGroup::Quant);
    let before = fs::read(current.payload_path.parent().unwrap().join("seal.json")).unwrap();
    let failed = failed_generation(&db);
    let original = fs::read(failed.payload_path.parent().unwrap().join("seal.json")).unwrap();
    assert_eq!(entry(&gc(root, &[]), &failed)["action"], "trash");
    assert!(failed.payload_path.exists());
    assert_eq!(
        entry(
            &gc(root, &["--apply", "--trash-days", "0", "--purge-expired"]),
            &failed
        )["state"],
        "TRASHED"
    );
    assert!(!failed.payload_path.exists());
    let output = command(root, &["storage", "usage"]);
    assert!(output.status.success());
    let usage: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        usage["categories"]["unpublished_trash"]["logical_bytes"]
            .as_u64()
            .unwrap()
            > 7
    );
    success(root, &["storage", "generations", "unpublished", "status"]);
    success(
        root,
        &[
            "storage",
            "generations",
            "unpublished",
            "restore",
            &failed.generation_id,
        ],
    );
    assert_eq!(
        fs::read(failed.payload_path.parent().unwrap().join("seal.json")).unwrap(),
        original
    );
    assert!(
        db.committed_generation_head(&failed.generation_id)
            .unwrap()
            .is_none()
    );
    assert_eq!(head(&db, OutputGroup::Quant), current);
    assert_eq!(
        entry(&gc(root, &["--apply", "--trash-days", "0"]), &failed)["state"],
        "TRASHED"
    );
    assert_eq!(
        entry(&gc(root, &["--apply", "--purge-expired"]), &failed)["state"],
        "PURGED"
    );
    assert!(
        !command(
            root,
            &[
                "storage",
                "generations",
                "unpublished",
                "restore",
                &failed.generation_id
            ]
        )
        .status
        .success()
    );
    assert_eq!(head(&db, OutputGroup::Quant), current);
    assert!(db.verify_generation_head(&current).unwrap());
    assert_eq!(
        fs::read(current.payload_path.parent().unwrap().join("seal.json")).unwrap(),
        before
    );
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), current);
}
#[test]
fn unpublished_cli_reports_live_pins_and_invalid_candidates_without_deleting_them() {
    let root = tempfile::tempdir().unwrap();
    let root = root.path();
    let source = Source::new();
    let (_paths, db) = setup(root, Some(&source), "both");
    let failed = failed_generation(&db);
    let pin =
        rnaa_core::generations::GenerationStorePin::acquire(db.path(), "CLI test reader").unwrap();
    let output = command(
        root,
        &["storage", "generations", "unpublished", "gc", "--apply"],
    );
    assert!(!output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["busy"], true);
    drop(pin);
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(root.join("rnaa.toml"), failed.payload_path.join("alias"))
            .unwrap();
        let output = command(
            root,
            &["storage", "generations", "unpublished", "gc", "--apply"],
        );
        assert!(!output.status.success());
        let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(entry(&report, &failed)["action"], "held_invalid");
        assert!(!report["errors"].as_array().unwrap().is_empty());
    }
    assert!(failed.payload_path.exists());
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1").unwrap());
}
