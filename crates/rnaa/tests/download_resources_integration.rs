//! Real download processes share capacity before HTTP and complete-copy publication.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/download.rs"]
mod download_support;
use commands::command;
use download_support::{Server, choose, download, head, setup};
use rnaa_core::{
    Database, ProjectConfig, ProjectPaths, RunState,
    generations::OutputGroup,
    resources::{Admission, ReservationOwner, ResourcePool, ResourceRequest},
};
use std::{
    fs,
    path::Path,
    process::{Child, Command, Stdio},
    sync::atomic::Ordering,
    time::{Duration, Instant},
};

struct Worker(Child);
impl Drop for Worker {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
fn spawn(root: &Path) -> Worker {
    Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args(["download", "--concurrency", "1"])
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
    let started = Instant::now();
    loop {
        if let Some(status) = worker.0.try_wait().unwrap() {
            return status;
        }
        assert!(
            started.elapsed() < Duration::from_secs(30),
            "{}",
            fs::read_to_string(root.join("worker.stderr")).unwrap()
        );
        std::thread::sleep(Duration::from_millis(10));
    }
}
fn save(paths: &ProjectPaths, db: &Database, pool: &Path, memory: u64) {
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.pool_root = pool.display().to_string();
    config.resources.memory_mb = memory;
    config.resources.cpu_threads = 1;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 15;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
}
fn attempt(db: &Database) -> Option<String> {
    db.list_job_attempts(100)
        .unwrap()
        .into_iter()
        .find(|a| a.stage == "download" && a.state == "RUNNING")
        .map(|a| a.attempt_id)
}
fn waiting(paths: &ProjectPaths) -> bool {
    fs::read_to_string(paths.logs_dir.join("download.log"))
        .unwrap_or_default()
        .contains("waiting for shared resource capacity")
}

#[test]
fn competing_cli_waits_cancel_before_http_and_survivor_publishes_under_its_claim() {
    let origin = Server::new();
    let root = tempfile::tempdir().unwrap();
    let a = root.path().join("a");
    let b = root.path().join("b");
    fs::create_dir_all(&a).unwrap();
    fs::create_dir_all(&b).unwrap();
    let (pa, da) = setup(&a);
    let (pb, db) = setup(&b);
    choose(
        &da,
        vec![origin.remote(&a, "/a.fastq", "reads.fastq", true)],
    );
    choose(
        &db,
        vec![origin.remote(&b, "/a.fastq", "reads.fastq", true)],
    );
    let pool_root = root.path().join("pool");
    save(&pa, &da, &pool_root, 256);
    save(&pb, &db, &pool_root, 256);
    let pool = ResourcePool::open(&pool_root).unwrap();
    let config = ProjectConfig::load(&pa.config_path).unwrap();
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
        Admission::Granted(g) => g,
        Admission::Wait(_) => panic!("unexpected contention"),
    };
    let mut wa = spawn(&a);
    let mut wb = spawn(&b);
    let started = Instant::now();
    while attempt(&da).is_none() || attempt(&db).is_none() || !waiting(&pa) || !waiting(&pb) {
        assert!(wa.0.try_wait().unwrap().is_none());
        assert!(wb.0.try_wait().unwrap().is_none());
        assert!(
            started.elapsed() < Duration::from_secs(10),
            "download never reached capacity wait"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
    assert_eq!(origin.count("GET"), 0);
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    assert!(da.cancel_job(&attempt(&da).unwrap()).unwrap());
    assert!(!finish(&mut wa, &a).success());
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    drop(held);
    assert!(finish(&mut wb, &b).success());
    assert_eq!(origin.count("GET"), 1);
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::Verified
    );
    let current = head(&db);
    let current_paths = current.paths(&pb).unwrap();
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(current_paths.manifest_path("download", "RUN1")).unwrap())
            .unwrap();
    assert_eq!(
        manifest["parameters"]["resource_execution"]["publication_reservation"]["owner"]["stage"],
        "download.publication"
    );
    assert_eq!(
        manifest["parameters"]["resource_execution"]["transfers"][0]["reservation"]["owner"]["stage"],
        "download.transfer"
    );
    let publication_disks = manifest["parameters"]["resource_execution"]["publication_reservation"]
        ["requested"]["disks"]
        .as_array()
        .unwrap();
    for path in [&pb.generations_dir, &pb.root] {
        let canonical = fs::canonicalize(path).unwrap();
        assert!(
            publication_disks
                .iter()
                .any(|d| d["path"].as_str() == canonical.to_str())
        );
    }
    assert!(db.verify_generation_head(&current).unwrap());
    assert!(pool.snapshot().unwrap().reservations.is_empty());
}

#[test]
fn publication_capacity_failure_preserves_previous_head_and_cached_retry_avoids_refetch() {
    let origin = Server::new();
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path());
    let pool_root = root.path().join("pool");
    save(&paths, &db, &pool_root, 1);
    let a = origin.remote(root.path(), "/a.fastq", "reads.fastq", true);
    choose(&db, vec![a.clone()]);
    assert!(
        !command(root.path(), &["download", "--concurrency", "1"])
            .status
            .success()
    );
    assert_eq!(origin.count("GET"), 0);
    assert!(
        db.generation_head(OutputGroup::Download, "RUN1")
            .unwrap()
            .is_none()
    );
    save(&paths, &db, &pool_root, 256);
    download(root.path());
    let first = head(&db);
    let b = origin.remote(root.path(), "/b.fastq", "reads.fastq", true);
    choose(&db, vec![b.clone()]);
    save(&paths, &db, &pool_root, 72);
    let failed = command(root.path(), &["download", "--concurrency", "1"]);
    assert!(!failed.status.success());
    assert_eq!(head(&db), first);
    assert!(db.verify_generation_head(&first).unwrap());
    assert_eq!(origin.count("GET"), 2);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Download, "RUN1")
            .unwrap()
            .len(),
        1
    );
    let pool = ResourcePool::open(&pool_root).unwrap();
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    save(&paths, &db, &pool_root, 256);
    download(root.path());
    let second = head(&db);
    assert_ne!(first, second);
    assert_eq!(origin.count("GET"), 2);
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    let unavailable = root.path().join("not-a-directory");
    fs::write(&unavailable, b"held").unwrap();
    save(&paths, &db, &unavailable, 1);
    origin.offline.store(true, Ordering::SeqCst);
    origin.head.store(false, Ordering::SeqCst);
    db.set_run_state("RUN1", RunState::Resolved, None).unwrap();
    download(root.path());
    assert_eq!(head(&db), second);
    assert_eq!(origin.count("GET"), 2);
    assert_eq!(origin.count("HEAD"), 0);
}
