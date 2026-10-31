//! Resource failures preserve the previous published correlation generation.
mod support;
use rnaa_core::{ProjectConfig, generations::OutputGroup, util::compute_sha256};
use std::fs;
use support::native::{command, setup, success};

#[test]
fn resource_rejection_preserves_head_and_allocation_changes_reuse_history() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path());
    success(
        root.path(),
        &[
            "deseq2",
            "--design",
            "~ condition",
            "--contrast",
            "condition",
            "B",
            "A",
        ],
    );
    success(
        root.path(),
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:100",
            "--out",
            "edges:topk=5",
        ],
    );
    let head = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let published = head.paths(&paths).unwrap();
    let output = published
        .corr_project_dir(db.project_id())
        .join("edges_topk.tsv");
    let hash = compute_sha256(&output).unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.corr.max_output_mb = 1;
    config.corr.block_rows = 7;
    config.corr.threads = 2;
    config.save(&paths.config_path).unwrap();
    let failed = command(root.path(), &["corr", "--out", "dense"]);
    assert!(!failed.status.success());
    let error = String::from_utf8_lossy(&failed.stderr);
    assert!(error.contains("output estimate"), "{error}");
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        head.generation_id
    );
    assert_eq!(compute_sha256(&output).unwrap(), hash);
    success(root.path(), &["corr", "--out", "edges:topk=5"]);
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        head.generation_id
    );
    assert_eq!(compute_sha256(&output).unwrap(), hash);
    let manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(published.manifest_path("corr", db.project_id())).unwrap(),
    )
    .unwrap();
    assert_eq!(
        manifest["parameters"]["limits"]["block_rows"], 64,
        "historical producer allocation is immutable"
    );
}
