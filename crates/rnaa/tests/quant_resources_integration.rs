//! Real CLI processes coordinate through the same pool and recover capacity.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/quant.rs"]
mod quant_support;
use commands::success;
use quant_support::setup;
use rnaa_core::{
    Database, ProjectConfig, ProjectPaths, RunState,
    generations::OutputGroup,
    resources::{Admission, ReservationOwner, ResourcePool, ResourceRequest},
};
use std::{
    fs,
    path::Path,
    process::{Child, Command, Stdio},
    time::{Duration, Instant},
};

struct Worker(Child);
impl Drop for Worker {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
fn save(paths: &ProjectPaths, db: &Database, pool: &Path, memory: u64, cpu: u64) {
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.pool_root = pool.display().to_string();
    config.resources.memory_mb = memory;
    config.resources.cpu_threads = cpu;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 20;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
}
fn spawn(root: &Path, args: &[&str]) -> Worker {
    Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args(args)
            .env("RUST_LOG", "info")
            .stdout(Stdio::from(
                fs::File::create(root.join("worker.stdout")).unwrap(),
            ))
            .stderr(Stdio::from(
                fs::File::create(root.join("worker.stderr")).unwrap(),
            ))
            .spawn()
            .unwrap(),
    )
}
fn finish(worker: &mut Worker, root: &Path) -> std::process::ExitStatus {
    let start = Instant::now();
    loop {
        if let Some(status) = worker.0.try_wait().unwrap() {
            return status;
        }
        assert!(
            start.elapsed() < Duration::from_secs(40),
            "worker did not terminate: {}",
            fs::read_to_string(root.join("worker.stderr")).unwrap()
        );
        std::thread::sleep(Duration::from_millis(10));
    }
}
fn running(db: &Database) -> Option<String> {
    db.list_job_attempts(100)
        .unwrap()
        .into_iter()
        .find(|a| a.stage == "quant" && a.state == "RUNNING")
        .map(|a| a.attempt_id)
}
fn waiting(paths: &ProjectPaths) -> bool {
    fs::read_to_string(paths.logs_dir.join("quant.log"))
        .unwrap_or_default()
        .contains("quantification waiting for shared resources")
}

#[test]
fn independent_cli_waits_cancel_without_claims_and_remaining_project_completes() {
    let root = tempfile::tempdir().unwrap();
    let a = root.path().join("a");
    let b = root.path().join("b");
    fs::create_dir_all(&a).unwrap();
    fs::create_dir_all(&b).unwrap();
    let (paths_a, db_a) = setup(&a, false);
    let (paths_b, db_b) = setup(&b, false);
    let pool_root = root.path().join("pool");
    save(&paths_a, &db_a, &pool_root, 512, 1);
    save(&paths_b, &db_b, &pool_root, 512, 1);
    let config = ProjectConfig::load(&paths_a.config_path).unwrap();
    let pool = ResourcePool::open(&pool_root).unwrap();
    let held = match pool
        .try_reserve(
            config.resources.limits().unwrap(),
            ReservationOwner {
                project_root: root.path().into(),
                attempt_id: "outside".into(),
                stage: "test.blocker".into(),
                run_accession: None,
            },
            ResourceRequest {
                cpu_threads: 1,
                ..Default::default()
            },
        )
        .unwrap()
    {
        Admission::Granted(guard) => guard,
        Admission::Wait(_) => panic!("unexpected contention"),
    };
    let args = ["quant", "--worker-budget", "1", "--threads", "1"];
    let mut worker_a = spawn(&a, &args);
    let mut worker_b = spawn(&b, &args);
    let start = Instant::now();
    while running(&db_a).is_none()
        || running(&db_b).is_none()
        || !waiting(&paths_a)
        || !waiting(&paths_b)
    {
        assert!(worker_a.0.try_wait().unwrap().is_none());
        assert!(worker_b.0.try_wait().unwrap().is_none());
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "CLI did not reach admission"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    assert!(
        db_a.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .is_none()
    );
    assert!(
        db_b.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .is_none()
    );
    assert!(db_a.cancel_job(&running(&db_a).unwrap()).unwrap());
    assert!(!finish(&mut worker_a, &a).success());
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    drop(held);
    let status = finish(&mut worker_b, &b);
    assert!(
        status.success(),
        "{}",
        fs::read_to_string(b.join("worker.stderr")).unwrap()
    );
    let current = paths_b.with_generation_heads(&db_b).unwrap();
    let mut index_ids = std::collections::BTreeSet::new();
    for run in ["RUN1", "RUN2", "RUN3", "RUN4"] {
        assert_eq!(
            db_b.get_run(run).unwrap().unwrap().state,
            RunState::QuantDone
        );
        let manifest: serde_json::Value = serde_json::from_slice(
            &fs::read(
                current
                    .quant_run_dir(run, "kallistors")
                    .join("rnaa_quant_manifest.json"),
            )
            .unwrap(),
        )
        .unwrap();
        let claims = manifest["parameters"]["resource_execution"]["reservations"]
            .as_array()
            .unwrap();
        assert_eq!(claims.len(), 2);
        let index = claims
            .iter()
            .find(|c| c["owner"]["stage"] == "quant.index")
            .unwrap();
        index_ids.insert(index["reservation_id"].as_str().unwrap().to_string());
    }
    assert_eq!(index_ids.len(), 1);
    assert!(pool.snapshot().unwrap().reservations.is_empty());
}

#[test]
fn constrained_preprocessing_and_quantification_feed_all_waves_without_deadlock() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path(), false);
    let pool_root = root.path().join("pool");
    save(&paths, &db, &pool_root, 224, 2);
    success(
        root.path(),
        &[
            "quant",
            "--preprocess",
            "--worker-budget",
            "2",
            "--threads",
            "1",
            "--preprocess-threads",
            "1",
        ],
    );
    for run in ["RUN1", "RUN2", "RUN3", "RUN4"] {
        assert_eq!(db.get_run(run).unwrap().unwrap().state, RunState::QuantDone);
    }
    assert!(
        ResourcePool::open(&pool_root)
            .unwrap()
            .snapshot()
            .unwrap()
            .reservations
            .is_empty()
    );
}
