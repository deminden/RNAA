//! Real native manifests migrate without changing their original generation seals.
use support::native as commands;
#[path = "support/recovery_fixture.rs"]
mod fixture;
mod support;
use commands::{command, success};
use fixture::{Source, head, setup};
use rnaa_core::{Database, generations::OutputGroup};
use std::{fs, path::Path, sync::atomic::Ordering};

fn downgrade(db: &Database, omit_inputs: bool) -> Vec<(std::path::PathBuf, Vec<u8>)> {
    let mut conn = rusqlite::Connection::open(db.path()).unwrap();
    let rows = {
        let mut query=conn.prepare("SELECT generation_id,inventory_json FROM artifact_generations WHERE state='COMMITTED'").unwrap();
        query
            .query_map([], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
            })
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
    };
    let tx = conn.transaction().unwrap();
    tx.execute("DELETE FROM generation_dependencies", [])
        .unwrap();
    let mut witnesses = Vec::new();
    for (id, json) in rows {
        let mut seal: serde_json::Value = serde_json::from_str(&json).unwrap();
        seal.as_object_mut().unwrap().remove("lineage_version");
        seal.as_object_mut().unwrap().remove("dependencies");
        if omit_inputs {
            seal.as_object_mut().unwrap().remove("inputs");
        }
        let path = Path::new(seal["head"]["payload_path"].as_str().unwrap())
            .parent()
            .unwrap()
            .join("seal.json");
        let bytes = serde_json::to_vec_pretty(&seal).unwrap();
        fs::write(&path, &bytes).unwrap();
        witnesses.push((path, bytes));
        tx.execute(
            "UPDATE artifact_generations SET inventory_json=?1 WHERE generation_id=?2",
            rusqlite::params![serde_json::to_string(&seal).unwrap(), id],
        )
        .unwrap();
    }
    tx.commit().unwrap();
    witnesses
}
fn import(root: &Path, apply: bool) -> serde_json::Value {
    let mut args = vec!["storage", "generations", "import-lineage"];
    if apply {
        args.push("--apply");
    }
    let output = command(root, &args);
    assert!(
        output.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn native_legacy_ancestry_recovers_from_manifests_and_preserves_count_reuse_and_input_cleanup() {
    for omit_inputs in [false, true] {
        let source = Source::new();
        let root = tempfile::tempdir().unwrap();
        let root = root.path();
        let (_paths, db) = setup(root, Some(&source), "both");
        let counts = head(&db, OutputGroup::Quant);
        let originals = downgrade(&db, omit_inputs);
        let preview = import(root, false);
        assert_eq!(preview["complete"], true);
        assert!(!root.join(".rnaa_support/lineage_imports").exists());
        let report = import(root, true);
        assert!(
            report["entries"]
                .as_array()
                .unwrap()
                .iter()
                .all(|entry| entry["action"] == "imported")
        );
        for (path, bytes) in &originals {
            assert_eq!(fs::read(path).unwrap(), *bytes);
        }
        let lineage = db.generation_lineage(&counts).unwrap();
        assert_eq!(lineage.source, "imported");
        assert!(lineage.dependencies.len() >= 2);
        success(
            root,
            &["lineage", "--group", "quant", "--history", "--verify"],
        );
        success(root, &["quant"]);
        assert_eq!(head(&db, OutputGroup::Quant), counts);
        success(
            root,
            &["storage", "files", "clean", "--mode", "fastq", "--apply"],
        );
        assert!(db.verify_generation(OutputGroup::Download, "RUN1").unwrap());
        assert!(
            db.verify_generation(OutputGroup::Preprocess, "RUN1")
                .unwrap()
        );
        assert!(db.verify_generation_head(&counts).unwrap());
        source.failed.store(true, Ordering::SeqCst);
        let requests = source.requests.load(Ordering::SeqCst);
        success(root, &["quant"]);
        assert_eq!(head(&db, OutputGroup::Quant), counts);
        success(root, &["quant", "--single-fragment-length", "180"]);
        assert_ne!(head(&db, OutputGroup::Quant), counts);
        assert_eq!(source.requests.load(Ordering::SeqCst), requests);
        for (path, bytes) in &originals {
            assert_eq!(fs::read(path).unwrap(), *bytes);
        }
        success(root, &["storage", "generations", "gc"]);
    }
}

#[test]
fn native_normalization_and_correlation_import_keep_their_frozen_sources_and_reuse_contracts() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = support::native::setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let de = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let corr = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let originals = downgrade(&db, true);
    let unresolved = command(root, &["export", "--what", "all", "--require-current"]);
    assert!(!unresolved.status.success());
    assert!(
        String::from_utf8_lossy(&unresolved.stderr).contains("unresolved input lineage"),
        "{}",
        String::from_utf8_lossy(&unresolved.stderr)
    );
    assert_eq!(import(root, true)["complete"], true);
    let lineage = db.generation_lineage(&corr).unwrap();
    assert_eq!(lineage.source, "imported");
    assert!(
        lineage
            .dependencies
            .iter()
            .any(|edge| edge.producer_generation_id == de.generation_id)
    );
    success(root, &["export", "--what", "all", "--require-current"]);
    success(root, &["normalize"]);
    success(root, &["corr"]);
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap(),
        de
    );
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap(),
        corr
    );
    for (path, bytes) in originals {
        assert_eq!(fs::read(path).unwrap(), bytes);
    }
    success(root, &["storage", "generations", "gc"]);
}

#[test]
fn usage_json_accounts_for_native_input_trash_retired_generations_and_unpublished_bytes() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let input_bytes = db
        .list_artifacts_for_run(Some("RUN1"))
        .unwrap()
        .iter()
        .filter(|record| record.kind.as_str().starts_with("FASTQ"))
        .map(|record| record.bytes)
        .sum::<u64>();
    success(root, &["quant", "--single-fragment-length", "180"]);
    success(
        root,
        &["storage", "files", "clean", "--mode", "fastq", "--apply"],
    );
    success(root, &["storage", "generations", "gc", "--apply"]);
    let mut lease = rnaa_core::jobs::JobLease::try_acquire(&db, "deseq2", None, &[])
        .unwrap()
        .unwrap();
    let owned = lease.database();
    let generation =
        rnaa_core::generations::Generation::begin(&owned, OutputGroup::De, "deseq2").unwrap();
    let unfinished = generation.head().payload_path.join("unfinished");
    fs::write(&unfinished, b"partial").unwrap();
    drop(generation);
    lease
        .finish(rnaa_core::jobs::JobOutcome::Failed, None)
        .unwrap();
    drop(owned);
    drop(lease);
    let unknown = root.join(".rnaa_generations/unknown/payload/unindexed");
    fs::create_dir_all(unknown.parent().unwrap()).unwrap();
    fs::write(&unknown, b"unindexed").unwrap();
    let output = command(root, &["storage", "usage", "--details"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["complete"], true);
    assert_eq!(
        report["categories"]["input_trash"]["logical_bytes"],
        input_bytes
    );
    assert!(
        report["categories"]["generation_trash"]["logical_bytes"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert_eq!(
        report["categories"]["abandoned_generation"]["logical_bytes"],
        7
    );
    assert_eq!(
        report["categories"]["unindexed_generation"]["logical_bytes"],
        9
    );
    assert_eq!(fs::read(unfinished).unwrap(), b"partial");
    assert_eq!(fs::read(unknown).unwrap(), b"unindexed");
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1").unwrap());
}
