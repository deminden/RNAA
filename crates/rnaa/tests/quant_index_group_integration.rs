//! Quant worker waves share one lazy index without changing scientific reuse.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/quant.rs"]
mod quant_support;
use commands::{command, success};
use quant_support::{record_raw, setup};
use rnaa_core::{Database, ProjectConfig, ProjectPaths, RunState, generations::OutputGroup};
use std::{collections::BTreeSet, fs};

fn index_groups(paths: &ProjectPaths, db: &Database, runs: &[&str]) -> BTreeSet<String> {
    let paths = paths.with_generation_heads(db).unwrap();
    runs.iter()
        .map(|run| {
            let manifest: serde_json::Value = serde_json::from_slice(
                &fs::read(
                    paths
                        .quant_run_dir(run, "kallistors")
                        .join("rnaa_quant_manifest.json"),
                )
                .unwrap(),
            )
            .unwrap();
            let info = &manifest["parameters"]["index_execution"];
            assert_eq!(info["load_attempts"], 1);
            assert_eq!(info["successful_loads"], 1);
            assert_eq!(
                info["index_sha256"],
                db.selected_reference().unwrap().unwrap().index.checksum
            );
            info["group_id"].as_str().unwrap().to_string()
        })
        .collect()
}

#[test]
fn index_peak_limit_rejects_new_quantification_and_preserves_healthy_reuse() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path(), false);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.memory_mb = 1;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let args = ["quant", "--worker-budget", "2", "--threads", "1"];
    let result = command(root.path(), &args);
    assert!(!result.status.success());
    for run in ["RUN1", "RUN2", "RUN3", "RUN4"] {
        let failed = db.get_run(run).unwrap().unwrap();
        assert_eq!(failed.state, RunState::QuantFailed);
        assert!(
            failed
                .last_error
                .unwrap()
                .contains("resource admission failed")
        );
        assert!(
            db.generation_head(OutputGroup::Quant, run)
                .unwrap()
                .is_none()
        );
    }
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.memory_mb = 512;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    success(root.path(), &args);
    let current = paths.with_generation_heads(&db).unwrap();
    let manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(
            current
                .quant_run_dir("RUN1", "kallistors")
                .join("rnaa_quant_manifest.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(
        manifest["parameters"]["index_execution"]["layout_plan"]["transcripts"],
        1
    );
    let heads = db.generation_heads().unwrap();
    let attempts = db.list_job_attempts(100).unwrap().len();
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.memory_mb = 1;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    success(root.path(), &args);
    assert_eq!(db.generation_heads().unwrap(), heads);
    assert_eq!(db.list_job_attempts(100).unwrap().len(), attempts);
}

#[test]
fn all_worker_waves_share_one_index_and_resource_changes_preserve_reuse() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path(), false);
    success(
        root.path(),
        &[
            "quant",
            "--worker-budget",
            "2",
            "--threads",
            "1",
            "--workers",
            "2",
        ],
    );
    let runs = ["RUN1", "RUN2", "RUN3", "RUN4"];
    assert_eq!(index_groups(&paths, &db, &runs).len(), 1);
    for run in runs {
        assert_eq!(db.get_run(run).unwrap().unwrap().state, RunState::QuantDone);
    }
    let heads = db.generation_heads().unwrap();
    let attempts = db.list_job_attempts(100).unwrap().len();
    let rejected = command(root.path(), &["quant", "--match-cache-mb", "0"]);
    assert!(!rejected.status.success());
    assert_eq!(db.generation_heads().unwrap(), heads);
    assert_eq!(
        ProjectConfig::load(&paths.config_path)
            .unwrap()
            .quant
            .match_cache_mb,
        64
    );
    success(
        root.path(),
        &[
            "quant",
            "--worker-budget",
            "1",
            "--threads",
            "1",
            "--workers",
            "1",
            "--match-cache-mb",
            "4",
        ],
    );
    assert_eq!(db.generation_heads().unwrap(), heads);
    assert_eq!(db.list_job_attempts(100).unwrap().len(), attempts);
    let original = index_groups(&paths, &db, &runs);
    success(
        root.path(),
        &[
            "quant",
            "--strand",
            "forward",
            "--worker-budget",
            "2",
            "--threads",
            "1",
            "--workers",
            "2",
        ],
    );
    let changed = index_groups(&paths, &db, &runs);
    assert_eq!(changed.len(), 1);
    assert_ne!(changed, original);
    success(root.path(), &["quant", "--strand", "unstranded"]);
    assert_eq!(index_groups(&paths, &db, &runs), original);
    for head in heads
        .into_iter()
        .filter(|head| head.output_group == OutputGroup::Quant)
    {
        assert_eq!(
            db.generation_head(OutputGroup::Quant, &head.subject)
                .unwrap()
                .unwrap(),
            head
        );
    }
}

#[test]
fn one_bad_input_does_not_discard_other_workers_shared_index_or_successful_outputs() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path(), true);
    let result = command(
        root.path(),
        &[
            "quant",
            "--worker-budget",
            "2",
            "--threads",
            "1",
            "--workers",
            "2",
        ],
    );
    assert!(!result.status.success());
    let good = ["RUN1", "RUN2", "RUN3"];
    let groups = index_groups(&paths, &db, &good);
    assert_eq!(groups.len(), 1);
    let heads = good
        .iter()
        .map(|run| {
            db.generation_head(OutputGroup::Quant, run)
                .unwrap()
                .unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        db.get_run("RUN4").unwrap().unwrap().state,
        RunState::QuantFailed
    );
    let raw = paths.raw_run_dir("RUN4").join("reads.fastq");
    fs::copy(paths.raw_run_dir("RUN1").join("reads.fastq"), &raw).unwrap();
    record_raw(&db, "RUN4", &raw);
    db.set_run_state("RUN4", RunState::Verified, None).unwrap();
    success(
        root.path(),
        &["quant", "--worker-budget", "1", "--threads", "1"],
    );
    assert_eq!(
        db.get_run("RUN4").unwrap().unwrap().state,
        RunState::QuantDone
    );
    assert_eq!(index_groups(&paths, &db, &good), groups);
    for head in heads {
        assert_eq!(
            db.generation_head(OutputGroup::Quant, &head.subject)
                .unwrap()
                .unwrap(),
            head
        );
    }
}
