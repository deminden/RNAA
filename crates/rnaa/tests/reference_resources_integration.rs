//! Reference admission and planning cancellation in real CLI processes.
#[path = "support/commands.rs"]
mod commands;
use commands::{command, success};
use rnaa_core::{
    Database, ProjectConfig, ProjectPaths,
    generations::OutputGroup,
    resources::{Admission, ReservationOwner, ResourcePool, ResourceRequest},
};
use std::{
    fs,
    io::Write,
    path::Path,
    process::{Child, Command, Stdio},
    time::{Duration, Instant},
};
fn sort_files(paths: &ProjectPaths) -> Vec<std::path::PathBuf> {
    let mut stack = vec![paths.generations_dir.clone()];
    let mut found = Vec::new();
    while let Some(directory) = stack.pop() {
        let Ok(entries) = fs::read_dir(directory) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if entry.file_type().is_ok_and(|t| t.is_dir()) {
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "keys") && path.is_file() {
                found.push(path);
            }
        }
    }
    found
}
struct Worker(Child);
impl Drop for Worker {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
fn setup(root: &Path, pool: &Path) -> (ProjectPaths, Database) {
    success(root, &["init"]);
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).unwrap();
    let fasta = root.join("source.fa");
    let gtf = root.join("source.gtf");
    fs::write(&fasta, format!(">tx1\n{}\n", "ACGT".repeat(150))).unwrap();
    fs::write(
        &gtf,
        "chr1\ttest\ttranscript\t1\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\";\n",
    )
    .unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.refs.custom_cdna = fasta.display().to_string();
    config.refs.custom_gtf = gtf.display().to_string();
    config.refs.index_threads = 8;
    config.resources.pool_root = pool.display().to_string();
    config.resources.cpu_threads = 1;
    config.resources.memory_mb = 256;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 15;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    (paths, db)
}
fn spawn(root: &Path) -> Worker {
    let scratch = root.join("scratch");
    fs::create_dir_all(&scratch).unwrap();
    Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args(["refs", "prepare"])
            .env("RUST_LOG", "info")
            .env("TMPDIR", scratch)
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
            start.elapsed() < Duration::from_secs(30),
            "{}",
            fs::read_to_string(root.join("worker.stderr")).unwrap()
        );
        std::thread::sleep(Duration::from_millis(10));
    }
}
fn attempt(db: &Database) -> Option<String> {
    db.list_job_attempts(100)
        .unwrap()
        .into_iter()
        .find(|a| a.stage == "refs" && a.state == "RUNNING")
        .map(|a| a.attempt_id)
}
fn waiting(paths: &ProjectPaths) -> bool {
    fs::read_to_string(paths.logs_dir.join("refs.log"))
        .unwrap_or_default()
        .contains("waiting for shared resource capacity")
}
fn large_source(root: &Path) {
    let mut file = std::io::BufWriter::new(fs::File::create(root.join("source.fa")).unwrap());
    writeln!(file, ">tx1").unwrap();
    let mut random = 59_u64;
    let mut chunk = vec![0_u8; 60000];
    for _ in 0..300 {
        for byte in &mut chunk {
            random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
            *byte = b"ACGT"[((random >> 32) & 3) as usize];
        }
        file.write_all(&chunk).unwrap();
        file.write_all(b"\n").unwrap();
    }
    file.flush().unwrap();
    drop(file);
}
#[test]
fn independent_reference_commands_wait_without_claims_and_cancellation_does_not_block_survivor() {
    let root = tempfile::tempdir().unwrap();
    let a = root.path().join("a");
    let b = root.path().join("b");
    fs::create_dir_all(&a).unwrap();
    fs::create_dir_all(&b).unwrap();
    let pool_root = root.path().join("pool");
    let (pa, da) = setup(&a, &pool_root);
    let (pb, db) = setup(&b, &pool_root);
    let pool = ResourcePool::open(&pool_root).unwrap();
    let limits = ProjectConfig::load(&pa.config_path)
        .unwrap()
        .resources
        .limits()
        .unwrap();
    let held = match pool
        .try_reserve(
            limits,
            ReservationOwner {
                project_root: root.path().into(),
                attempt_id: "blocker".into(),
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
    let (mut wa, mut wb) = (spawn(&a), spawn(&b));
    let start = Instant::now();
    while attempt(&da).is_none() || attempt(&db).is_none() || !waiting(&pa) || !waiting(&pb) {
        assert!(wa.0.try_wait().unwrap().is_none());
        assert!(wb.0.try_wait().unwrap().is_none());
        assert!(
            start.elapsed() < Duration::from_secs(10),
            "reference did not reach capacity wait"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    assert!(da.cancel_job(&attempt(&da).unwrap()).unwrap());
    assert!(!finish(&mut wa, &a).success());
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    drop(held);
    assert!(
        finish(&mut wb, &b).success(),
        "{}",
        fs::read_to_string(b.join("worker.stderr")).unwrap()
    );
    let selection = db.selected_reference().unwrap().unwrap();
    assert!(db.verify_generation_head(&selection.generation).unwrap());
    let current = selection.generation.paths(&pb).unwrap();
    let manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(current.manifest_path("refs", &selection.reference_id)).unwrap(),
    )
    .unwrap();
    assert_eq!(
        manifest["parameters"]["resource_execution"]["build_reservation"]["owner"]["stage"],
        "refs.build"
    );
    assert_eq!(manifest["parameters"]["index_threads"], 1);
    assert_eq!(manifest["parameters"]["requested_index_threads"], 8);
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    assert_eq!(fs::read_dir(b.join("scratch")).unwrap().count(), 0);
    assert!(sort_files(&pb).is_empty());
}
#[test]
fn cancellation_during_external_planning_removes_sort_runs_and_never_publishes() {
    let root = tempfile::tempdir().unwrap();
    let pool_root = root.path().join("pool");
    let (paths, db) = setup(root.path(), &pool_root);
    large_source(root.path());
    let pool = ResourcePool::open(&pool_root).unwrap();
    let mut worker = spawn(root.path());
    let start = Instant::now();
    loop {
        if pool
            .snapshot()
            .unwrap()
            .reservations
            .iter()
            .any(|r| r.owner.stage == "refs.fasta_plan")
        {
            break;
        }
        assert!(
            worker.0.try_wait().unwrap().is_none(),
            "planner exited early"
        );
        assert!(start.elapsed() < Duration::from_secs(10));
        std::thread::sleep(Duration::from_millis(5));
    }
    assert!(db.cancel_job(&attempt(&db).unwrap()).unwrap());
    assert!(!finish(&mut worker, root.path()).success());
    assert!(db.selected_reference().unwrap().is_none());
    assert!(
        db.generation_heads()
            .unwrap()
            .iter()
            .all(|h| h.output_group != OutputGroup::Refs)
    );
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    assert_eq!(fs::read_dir(paths.root.join("scratch")).unwrap().count(), 0);
    assert!(sort_files(&paths).is_empty());
}

#[test]
fn killed_planner_scratch_is_owned_and_recoverable_by_unpublished_cleanup() {
    let root = tempfile::tempdir().unwrap();
    let pool_root = root.path().join("pool");
    let (paths, db) = setup(root.path(), &pool_root);
    large_source(root.path());
    let pool = ResourcePool::open(&pool_root).unwrap();
    let mut worker = spawn(root.path());
    let start = Instant::now();
    loop {
        if !sort_files(&paths).is_empty() {
            break;
        }
        assert!(
            worker.0.try_wait().unwrap().is_none(),
            "planner exited before writing sort runs"
        );
        assert!(start.elapsed() < Duration::from_secs(15));
        std::thread::sleep(Duration::from_millis(5));
    }
    let attempt = attempt(&db).unwrap();
    worker.0.kill().unwrap();
    assert!(!finish(&mut worker, root.path()).success());
    assert!(!sort_files(&paths).is_empty());
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    // Simulate lease expiry only after the process has demonstrably exited.
    rusqlite::Connection::open(&paths.db_path)
        .unwrap()
        .execute(
            "UPDATE job_attempts SET lease_expires_ms=0 WHERE attempt_id=?1",
            [&attempt],
        )
        .unwrap();
    db.recover_jobs().unwrap();
    let gc = |extra: &[&str]| -> serde_json::Value {
        let mut args = vec!["storage", "generations", "unpublished", "gc"];
        args.extend(extra);
        let result = command(root.path(), &args);
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        serde_json::from_slice(&result.stdout).unwrap()
    };
    let preview = gc(&[]);
    let id = preview["entries"]
        .as_array()
        .unwrap()
        .iter()
        .find(|e| e["action"] == "trash")
        .expect("dead planner must be a tracked cleanup candidate")["generation_id"]
        .as_str()
        .unwrap()
        .to_string();
    let trashed = gc(&["--apply", "--trash-days", "0"]);
    assert!(
        trashed["entries"]
            .as_array()
            .unwrap()
            .iter()
            .any(|e| e["generation_id"] == id && e["state"] == "TRASHED")
    );
    assert!(sort_files(&paths).is_empty());
    let purged = gc(&["--apply", "--purge-expired"]);
    assert!(
        purged["entries"]
            .as_array()
            .unwrap()
            .iter()
            .any(|e| e["generation_id"] == id && e["state"] == "PURGED")
    );
    assert!(root.path().join("source.fa").is_file());
    assert!(db.selected_reference().unwrap().is_none());
}
