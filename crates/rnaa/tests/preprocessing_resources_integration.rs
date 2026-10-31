//! Actual preprocessing uses the common pool, including failures and reuse.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/quant.rs"]
mod quant_support;
use commands::{command, success};
use quant_support::setup;
use rnaa_core::{ProjectConfig, RunState, generations::OutputGroup, resources::ResourcePool};
use std::fs;
const RUN: &[&str] = &[
    "quant",
    "--preprocess",
    "--preprocess-bypass",
    "--worker-budget",
    "1",
    "--threads",
    "1",
    "--preprocess-threads",
    "1",
];

#[test]
fn insufficient_resources_cannot_bypass_and_healthy_outputs_ignore_lowered_limits() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path(), false);
    let pool_path = root.path().join("pool");
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.pool_root = pool_path.display().to_string();
    config.resources.cpu_threads = 1;
    config.resources.memory_mb = 1;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 1;
    for memory_mb in [1, 64] {
        config.resources.memory_mb = memory_mb;
        config.save(&paths.config_path).unwrap();
        db.set_project_config(&config).unwrap();
        let result = command(root.path(), RUN);
        assert!(!result.status.success());
        for run in ["RUN1", "RUN2", "RUN3", "RUN4"] {
            let state = db.get_run(run).unwrap().unwrap();
            assert_eq!(state.state, RunState::PreprocessFailed);
            assert!(
                state
                    .last_error
                    .unwrap()
                    .contains("resource admission failed")
            );
            assert!(
                db.generation_head(OutputGroup::Preprocess, run)
                    .unwrap()
                    .is_none()
            );
            assert!(
                db.generation_head(OutputGroup::Quant, run)
                    .unwrap()
                    .is_none()
            );
        }
    }
    let pool = ResourcePool::open(&pool_path).unwrap();
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.memory_mb = 512;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    success(root.path(), RUN);
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    let current = paths.with_generation_heads(&db).unwrap();
    let report: serde_json::Value = serde_json::from_slice(
        &fs::read(
            current
                .preprocess_run_dir("RUN1", "kallistors")
                .join("fasterp.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(report["mode"], "single");
    assert_eq!(
        report["execution"]["resource_reservation"]["owner"]["stage"],
        "preprocess.engine"
    );
    assert_eq!(
        report["execution"]["resource_reservation"]["requested"]["cpu_threads"],
        1
    );
    let heads = db.generation_heads().unwrap();
    let attempts = db.list_job_attempts(100).unwrap().len();
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.memory_mb = 1;
    let unavailable_pool = root.path().join("unavailable-pool");
    fs::write(&unavailable_pool, "not a directory").unwrap();
    config.resources.pool_root = unavailable_pool.display().to_string();
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    success(root.path(), RUN);
    assert_eq!(db.generation_heads().unwrap(), heads);
    assert_eq!(db.list_job_attempts(100).unwrap().len(), attempts);
    assert!(pool.snapshot().unwrap().reservations.is_empty());
}
