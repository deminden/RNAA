use super::*;

#[cfg(unix)]
#[test]
fn next_export_preserves_a_tracked_interrupted_directory() {
    use rnaa_core::jobs::{JobLease, JobOutcome};
    use std::os::unix::fs::MetadataExt;
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize", "--design", "~ condition"]);
    let mut lease = JobLease::try_acquire(&db, "export", None, &[])
        .unwrap()
        .unwrap();
    let attempt = lease.attempt_id().to_owned();
    let name = format!(".rnaa-export-{attempt}-fixture");
    let staged = paths.exports_dir.join(&name);
    fs::create_dir(&staged).unwrap();
    fs::write(staged.join("partial"), b"diagnostic data").unwrap();
    let meta = fs::metadata(&staged).unwrap();
    lease
        .database()
        .register_export_stage("export", &name, meta.dev(), meta.ino())
        .unwrap();
    lease
        .finish(JobOutcome::Interrupted, Some("test interruption"))
        .unwrap();
    drop(lease);

    let output = exported(root, &["export", "--what", "counts"]);
    assert!(output.join("de/gene_counts.tsv").exists());
    let catalog = rusqlite::Connection::open(&paths.db_path).unwrap();
    let (old_state, new_state): (String, i64) = catalog.query_row("SELECT (SELECT state FROM export_publications WHERE attempt_id=?1), (SELECT count(*) FROM export_publications p JOIN job_attempts j ON j.attempt_id=p.attempt_id WHERE p.state='PUBLISHED' AND j.state='SUCCEEDED')",[&attempt],|row| Ok((row.get(0)?,row.get(1)?))).unwrap();
    assert_eq!(old_state, "RECOVERED");
    assert_eq!(new_state, 1);
    let recovered_name = format!(".rnaa-recovered-{}", &name[1..]);
    assert_eq!(
        fs::read(paths.exports_dir.join(&recovered_name).join("partial")).unwrap(),
        b"diagnostic data"
    );
    assert!(!staged.exists());
    let events = db.list_events(100).unwrap();
    assert!(
        events
            .iter()
            .any(|event| event.message == "interrupted export reconciled"
                && event.context["attempt_id"] == attempt
                && event.context["exports_relative_path"] == recovered_name)
    );
}

#[test]
fn sealed_tables_export_with_provenance_and_failed_replacements_preserve_results() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize", "--design", "~ condition"]);
    success(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=5",
        ],
    );
    let current = paths.with_generation_heads(&db).unwrap();
    let source = current
        .de_project_dir(db.project_id())
        .join("gene_counts.tsv");
    let first = exported(root, &["export", "--what", "counts", "--format", "csv"]);
    let table = first.join("de/gene_counts.csv");
    assert_eq!(rows(&source, b'\t'), rows(&table, b','));
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(first.join("export-manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["format"], "csv");
    let entry = manifest["files"]
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| entry["output"] == "de/gene_counts.csv")
        .unwrap();
    assert_eq!(entry["output_sha256"], compute_sha256(&table).unwrap());
    assert_eq!(entry["rows"], 500);
    assert_eq!(entry["columns"], 9);
    assert!(manifest["files"].as_array().unwrap().iter().any(|entry| {
        entry["output"]
            .as_str()
            .unwrap()
            .starts_with("provenance/de/")
    }));
    let second = exported(
        root,
        &["export", "--what", "correlation", "--format", "tsv"],
    );
    assert!(second.join("corr/edges_topk.tsv").exists());
    let complete = exported(root, &["export", "--what", "all", "--format", "tsv"]);
    let inventory: serde_json::Value =
        serde_json::from_slice(&fs::read(complete.join("export-manifest.json")).unwrap()).unwrap();
    assert_eq!(inventory["generations"].as_array().unwrap().len(), 2);
    assert!(complete.join("de/gene_counts.tsv").exists());
    assert!(complete.join("corr/edges_topk.tsv").exists());
    let first_hash = compute_sha256(&table).unwrap();
    let before = fs::read_dir(&paths.exports_dir).unwrap().count();
    assert!(
        !command(root, &["export", "--format", "bogus"])
            .status
            .success()
    );
    assert!(
        !command(root, &["export", "--what", "bogus"])
            .status
            .success()
    );
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    let original_config = fs::read(&paths.config_path).unwrap();
    config.resources.memory_mb = 1;
    config.save(&paths.config_path).unwrap();
    assert!(
        !command(root, &["export", "--what", "counts"])
            .status
            .success()
    );
    fs::write(&paths.config_path, original_config).unwrap();
    let mut corrupted = fs::read(&source).unwrap();
    corrupted[0] ^= 1;
    fs::write(&source, corrupted).unwrap();
    assert!(
        !command(root, &["export", "--what", "counts"])
            .status
            .success()
    );
    assert_eq!(fs::read_dir(&paths.exports_dir).unwrap().count(), before);
    assert_eq!(compute_sha256(&table).unwrap(), first_hash);
    assert!(fs::read_dir(&paths.exports_dir).unwrap().all(|entry| {
        !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".rnaa-export-")
    }));
}

#[test]
fn cancellation_while_waiting_for_capacity_leaves_no_export_or_claim() {
    use rnaa_core::resources::{Admission, ReservationOwner, ResourcePool, ResourceRequest};
    use std::{
        process::{Command, Stdio},
        time::{Duration, Instant},
    };
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize", "--design", "~ condition"]);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    let pool_root = root.join("private-pool");
    config.resources.pool_root = pool_root.display().to_string();
    config.resources.cpu_threads = 1;
    config.resources.memory_mb = 512;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 20;
    config.save(&paths.config_path).unwrap();
    let pool = ResourcePool::open(&pool_root).unwrap();
    let held = match pool
        .try_reserve(
            config.resources.limits().unwrap(),
            ReservationOwner {
                project_root: root.into(),
                attempt_id: "external-blocker".into(),
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
        Admission::Granted(held) => held,
        Admission::Wait(_) => panic!("unexpected contention"),
    };
    let mut worker = Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args(["export", "--what", "counts"])
            .env("RUST_LOG", "info")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .unwrap(),
    );
    let start = Instant::now();
    let attempt = loop {
        let job = db
            .list_job_attempts(100)
            .unwrap()
            .into_iter()
            .find(|job| job.stage == "export" && job.state == "RUNNING");
        let waiting = fs::read_to_string(paths.logs_dir.join("export.log"))
            .unwrap_or_default()
            .contains("waiting for shared resource capacity");
        if let Some(job) = job
            && waiting
        {
            break job.attempt_id;
        }
        assert!(
            worker.0.try_wait().unwrap().is_none(),
            "export exited before waiting"
        );
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "export did not reach admission"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(db.cancel_job(&attempt).unwrap());
    let start = Instant::now();
    let status = loop {
        if let Some(status) = worker.0.try_wait().unwrap() {
            break status;
        }
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "cancelled export did not stop"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(!status.success());
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    assert_eq!(fs::read_dir(&paths.exports_dir).unwrap().count(), 0);
    drop(held);
    assert!(pool.snapshot().unwrap().reservations.is_empty());
}
