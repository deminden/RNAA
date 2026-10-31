//! Returning to an exact scientific request selects verified immutable history.
mod support;
use rnaa_core::{
    Database, ManifestArtifact,
    generations::{GenerationHead, OutputGroup},
};
use std::{fs, path::Path};
use support::native::{command, setup, success};

fn head(db: &Database, group: OutputGroup) -> GenerationHead {
    db.generation_head(group, db.project_id()).unwrap().unwrap()
}

fn history_count(root: &Path, group: &str) -> usize {
    let result = command(root, &["lineage", "--group", group, "--all-requests"]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert!(
        value["generations"]
            .as_array()
            .unwrap()
            .iter()
            .all(|entry| entry["lineage"]["head"]["output_group"] == group)
    );
    assert_eq!(
        value["generations"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|entry| entry["selected"] == true)
            .count(),
        1
    );
    value["generations"].as_array().unwrap().len()
}

fn registry(db: &Database, head: &GenerationHead) -> serde_json::Value {
    serde_json::to_value(
        db.list_artifacts()
            .unwrap()
            .into_iter()
            .filter(|record| Path::new(&record.path).starts_with(&head.payload_path))
            .collect::<Vec<_>>(),
    )
    .unwrap()
}

#[test]
fn returning_to_metadata_request_reattaches_de_and_downstream_corr_without_new_payloads() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let first_de = head(&db, OutputGroup::De);
    let first_corr = head(&db, OutputGroup::Corr);
    let original_de_records = registry(&db, &first_de);
    let original_corr_records = registry(&db, &first_corr);
    let original_sheet = fs::read_to_string(paths.samplesheet_path()).unwrap();
    let seal = first_de.payload_path.parent().unwrap().join("seal.json");
    let seal_bytes = fs::read(&seal).unwrap();
    let matrix = ManifestArtifact::snapshot(
        "VST",
        &first_de
            .paths(&paths)
            .unwrap()
            .de_project_dir(db.project_id())
            .join("vst.tsv"),
    )
    .unwrap();
    fs::write(
        paths.samplesheet_path(),
        original_sheet
            .replace("\tA\n", "\tcontrol\n")
            .replace("\tB\n", "\ttreated\n"),
    )
    .unwrap();
    success(root, &["normalize"]);
    success(root, &["corr", "--out", "edges:topk=7"]);
    let second_de = head(&db, OutputGroup::De);
    let second_corr = head(&db, OutputGroup::Corr);
    assert_ne!(second_de, first_de);
    assert_ne!(second_corr, first_corr);
    assert_eq!(history_count(root, "de"), 2);
    assert_eq!(history_count(root, "corr"), 2);

    fs::write(paths.samplesheet_path(), &original_sheet).unwrap();
    success(root, &["normalize"]);
    assert_eq!(head(&db, OutputGroup::De), first_de);
    success(root, &["corr", "--out", "edges:topk=5"]);
    assert_eq!(head(&db, OutputGroup::Corr), first_corr);
    assert_eq!(history_count(root, "de"), 2);
    assert_eq!(history_count(root, "corr"), 2);
    assert_eq!(registry(&db, &first_de), original_de_records);
    assert_eq!(registry(&db, &first_corr), original_corr_records);
    assert_eq!(fs::read(seal).unwrap(), seal_bytes);
    assert!(matrix.matches_file().unwrap());
    assert!(db.verify_generation_head(&second_de).unwrap());
    assert!(db.verify_generation_head(&second_corr).unwrap());
    for selected in [&first_de, &first_corr] {
        let event = db
            .list_events(1000)
            .unwrap()
            .into_iter()
            .find(|event| {
                event.message == "selected verified committed generation"
                    && event.context["generation_id"] == selected.generation_id
            })
            .unwrap();
        assert_eq!(event.context["producer_attempt_id"], selected.attempt_id);
        assert_ne!(event.context["reuse_attempt_id"], selected.attempt_id);
        let attempt = db
            .list_job_attempts(1000)
            .unwrap()
            .into_iter()
            .find(|attempt| {
                Some(attempt.attempt_id.as_str()) == event.context["reuse_attempt_id"].as_str()
            })
            .unwrap();
        assert_eq!(attempt.state, "SUCCEEDED");
        assert_eq!(
            attempt.request_sha256.as_deref(),
            Some(selected.request_sha256.as_str())
        );
    }
}

#[test]
fn corrupt_historical_payload_is_skipped_and_recomputed_without_overwriting_history() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let first = head(&db, OutputGroup::Corr);
    success(root, &["corr", "--out", "edges:topk=7"]);
    let second = head(&db, OutputGroup::Corr);
    let corrupt = first
        .paths(&paths)
        .unwrap()
        .corr_project_dir(db.project_id())
        .join("edges_topk.tsv");
    fs::write(&corrupt, "corruption retained for diagnosis").unwrap();
    success(root, &["corr", "--out", "edges:topk=5"]);
    let replacement = head(&db, OutputGroup::Corr);
    assert_ne!(replacement, first);
    assert_ne!(replacement, second);
    assert_eq!(replacement.request_sha256, first.request_sha256);
    assert_eq!(history_count(root, "corr"), 3);
    assert!(!db.verify_generation_head(&first).unwrap());
    assert!(db.verify_generation_head(&second).unwrap());
    assert!(db.verify_generation_head(&replacement).unwrap());
    assert_eq!(
        fs::read_to_string(corrupt).unwrap(),
        "corruption retained for diagnosis"
    );
    // The healthy selected candidate wins over the older corrupt same-request entry.
    success(root, &["corr"]);
    assert_eq!(head(&db, OutputGroup::Corr), replacement);
    assert_eq!(history_count(root, "corr"), 3);
    let inspected = command(
        root,
        &["lineage", "--group", "corr", "--all-requests", "--verify"],
    );
    assert!(!inspected.status.success());
    let report: serde_json::Value = serde_json::from_slice(&inspected.stdout).unwrap();
    let entries = report["generations"].as_array().unwrap();
    assert_eq!(
        entries
            .iter()
            .filter(|entry| entry["payload_verified"] == false)
            .count(),
        1
    );
    assert!(
        entries
            .iter()
            .any(|entry| entry["selected"] == true && entry["payload_verified"] == true)
    );
}

#[test]
fn normalization_and_contrast_results_keep_separate_request_contracts_in_shared_de_group() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root);
    success(root, &["normalize"]);
    let normalized = head(&db, OutputGroup::De);
    success(root, &["deseq2", "--contrast", "condition", "B", "A"]);
    let de = head(&db, OutputGroup::De);
    assert_ne!(normalized, de);
    success(root, &["normalize"]);
    assert_eq!(head(&db, OutputGroup::De), normalized);
    success(root, &["deseq2", "--contrast", "condition", "B", "A"]);
    assert_eq!(head(&db, OutputGroup::De), de);
    assert_eq!(history_count(root, "de"), 2);
}
