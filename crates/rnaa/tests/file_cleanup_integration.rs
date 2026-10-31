//! Native input cleanup, count reuse, expiration and recovery through public commands.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/recovery_fixture.rs"]
mod fixture;
use commands::{command, success};
use fixture::{head, setup};
use rnaa_core::{
    ArtifactKind, ArtifactRecord, Database, ProjectConfig,
    generations::OutputGroup,
    util::{compute_sha256, now_rfc3339},
};
use std::{fs, path::Path};

fn files(root: &Path, args: &[&str]) -> serde_json::Value {
    let mut full = vec!["storage", "files"];
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
fn fastqs(db: &Database) -> Vec<ArtifactRecord> {
    db.list_artifacts_for_run(Some("RUN1"))
        .unwrap()
        .into_iter()
        .filter(|record| {
            matches!(
                record.kind,
                ArtifactKind::FastqSingle | ArtifactKind::FastqTrimmedSingle
            )
        })
        .collect()
}

#[test]
fn recovery_origin_waits_for_http_headers_and_counts_requests_instead_of_connections() {
    use std::io::{Read, Write};
    let source = fixture::Source::new();
    let address = source
        .url
        .strip_prefix("http://")
        .unwrap()
        .split('/')
        .next()
        .unwrap();
    let mut stream = std::net::TcpStream::connect(address).unwrap();
    stream
        .set_read_timeout(Some(std::time::Duration::from_secs(2)))
        .unwrap();
    std::thread::sleep(std::time::Duration::from_millis(100));
    assert_eq!(source.requests.load(std::sync::atomic::Ordering::SeqCst), 0);
    stream.write_all(b"GET /RUN1.fastq HTTP/1.1\r\n").unwrap();
    std::thread::sleep(std::time::Duration::from_millis(100));
    assert_eq!(source.requests.load(std::sync::atomic::Ordering::SeqCst), 0);
    stream
        .write_all(b"Host: localhost\r\nConnection: close\r\n\r\n")
        .unwrap();
    let mut response = Vec::new();
    stream.read_to_end(&mut response).unwrap();
    assert!(response.starts_with(b"HTTP/1.1 200 OK\r\n"));
    assert!(response.ends_with(&fixture::reads()));
    assert_eq!(source.requests.load(std::sync::atomic::Ordering::SeqCst), 1);
}
fn config(root: &Path, db: &Database, edit: impl FnOnce(&mut ProjectConfig)) {
    let path = root.join("rnaa.toml");
    let mut config = ProjectConfig::load(&path).unwrap();
    edit(&mut config);
    config.save(&path).unwrap();
    db.set_project_config(&config).unwrap();
}

#[test]
fn preview_trash_restore_and_purge_preserve_counts_and_end_absence_authorization() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let counts = head(&db, OutputGroup::Quant);
    let pre = head(&db, OutputGroup::Preprocess);
    let reads = fastqs(&db);
    assert_eq!(reads.len(), 2);
    let preview = files(root, &["clean", "--mode", "fastq", "--trash-days", "0"]);
    assert_eq!(preview["entries"].as_array().unwrap().len(), 2);
    assert!(reads.iter().all(|record| Path::new(&record.path).exists()));
    assert!(db.retained_inputs("RUN1").unwrap().is_empty());
    let trashed = files(
        root,
        &["clean", "--mode", "fastq", "--apply", "--trash-days", "0"],
    );
    assert!(
        trashed["entries"]
            .as_array()
            .unwrap()
            .iter()
            .all(|row| row["state"] == "TRASHED")
    );
    assert!(fastqs(&db).is_empty());
    assert!(reads.iter().all(|record| !Path::new(&record.path).exists()));
    assert_eq!(db.retained_inputs("RUN1").unwrap().len(), 2);
    assert!(db.verify_generation_head(&pre).unwrap());
    assert!(db.verify_generation_head(&counts).unwrap());
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    let id = trashed["entries"][0]["trash_id"].as_str().unwrap();
    let restored = files(root, &["restore", id]);
    assert_eq!(restored["restored"].as_array().unwrap().len(), 2);
    assert_eq!(fastqs(&db).len(), 2);
    assert!(db.retained_inputs("RUN1").unwrap().is_empty());
    // Restoring ends the absence permission: an unrelated loss is damage again.
    let trimmed = reads
        .iter()
        .find(|row| row.kind == ArtifactKind::FastqTrimmedSingle)
        .unwrap();
    let body = fs::read(&trimmed.path).unwrap();
    fs::remove_file(&trimmed.path).unwrap();
    assert!(!db.verify_generation_head(&pre).unwrap());
    fs::write(&trimmed.path, body).unwrap();
    files(
        root,
        &["clean", "--mode", "fastq", "--apply", "--trash-days", "0"],
    );
    let purged = files(root, &["gc", "--apply", "--purge-expired"]);
    assert_eq!(
        purged["entries"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|row| row["state"] == "PURGED")
            .count(),
        2
    );
    assert!(db.verify_generation_head(&pre).unwrap());
    assert!(db.verify_generation_head(&counts).unwrap());
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
}

#[test]
fn changed_count_policy_restores_trashed_originals_without_a_download_source() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let old = head(&db, OutputGroup::Quant);
    let reads = fastqs(&db);
    files(root, &["clean", "--mode", "fastq", "--apply"]);
    success(root, &["quant", "--single-fragment-length", "180"]);
    assert_ne!(head(&db, OutputGroup::Quant), old);
    assert_eq!(fastqs(&db).len(), 2);
    assert!(reads.iter().all(|record| Path::new(&record.path).exists()));
}

#[test]
fn automatic_cleanup_schedule_and_reversible_retention_precedence_are_effective() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let counts = head(&db, OutputGroup::Quant);
    success(
        root,
        &["quant", "--cleanup", "fastq", "--cleanup-on", "never"],
    );
    assert_eq!(fastqs(&db).len(), 2);
    success(
        root,
        &[
            "quant",
            "--cleanup-on",
            "success",
            "--fastq-retention",
            "none",
        ],
    );
    assert!(fastqs(&db).is_empty());
    let status = files(root, &["status"]);
    assert_eq!(status["entries"].as_array().unwrap().len(), 2);
    files(
        root,
        &[
            "restore",
            status["entries"][0]["trash_id"].as_str().unwrap(),
        ],
    );
    assert_eq!(fastqs(&db).len(), 2);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
}

#[test]
fn unused_project_archive_cleanup_is_separate_and_restore_rejects_collisions() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let archive = root.join("raw/RUN1/RUN1.sra");
    fs::write(&archive, b"archive").unwrap();
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: Some("RUN1".into()),
        kind: ArtifactKind::Sra,
        path: archive.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".into(),
        checksum: compute_sha256(&archive).unwrap(),
        bytes: 7,
        created_at: now_rfc3339(),
    })
    .unwrap();
    let report = files(root, &["clean", "--mode", "sra", "--apply"]);
    assert_eq!(fastqs(&db).len(), 2);
    assert!(!archive.exists());
    fs::write(&archive, b"collision").unwrap();
    let id = report["entries"][0]["trash_id"].as_str().unwrap();
    assert!(
        !command(root, &["storage", "files", "restore", id])
            .status
            .success()
    );
    assert_eq!(fs::read(&archive).unwrap(), b"collision");
    fs::remove_file(&archive).unwrap();
    // The durable RESTORING intent is resumed before the requested operation.
    let retried = files(root, &["restore", id]);
    assert_eq!(retried["already_restored"], true);
    assert_eq!(fs::read(&archive).unwrap(), b"archive");
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1").unwrap());
}

#[test]
fn changed_source_registration_restores_physical_data_without_selecting_stale_inputs() {
    let source = fixture::Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "both");
    let reads = fastqs(&db);
    let report = files(root, &["clean", "--mode", "fastq", "--apply"]);
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.remote_files[0].url.push_str("?new_source");
    db.upsert_run(&run).unwrap();
    assert!(!db.restorable_file_trash("RUN1").unwrap());
    files(
        root,
        &[
            "restore",
            report["entries"][0]["trash_id"].as_str().unwrap(),
        ],
    );
    assert!(reads.iter().all(|record| Path::new(&record.path).exists()));
    assert!(fastqs(&db).is_empty());
    config(root, &db, |config| {
        config.quant.cleanup = rnaa_core::CleanupMode::None
    });
}

#[test]
fn always_cleanup_preserves_the_processing_error_and_success_schedule_holds_on_failure() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let counts = head(&db, OutputGroup::Quant);
    config(root, &db, |config| {
        config.quant.cleanup = rnaa_core::CleanupMode::Fastq;
        config.quant.cleanup_on = rnaa_core::CleanupWhen::Success;
        config.refs.custom_cdna = root.join("unavailable.fa").display().to_string();
    });
    let failed = command(root, &["quant"]);
    assert!(!failed.status.success());
    assert_eq!(fastqs(&db).len(), 2);
    config(root, &db, |config| {
        config.quant.cleanup_on = rnaa_core::CleanupWhen::Always
    });
    let failed = command(root, &["quant"]);
    assert!(!failed.status.success());
    assert!(
        String::from_utf8_lossy(&failed.stderr).contains("custom reference files do not exist")
    );
    assert!(fastqs(&db).is_empty());
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    assert!(db.verify_generation_head(&counts).unwrap());
    assert_eq!(
        files(root, &["status"])["entries"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn historical_count_binding_cleanup_retains_new_input_evidence_without_rewriting_old_history() {
    let source = fixture::Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "both");
    let counts = head(&db, OutputGroup::Quant);
    files(
        root,
        &["clean", "--mode", "fastq", "--apply", "--trash-days", "0"],
    );
    files(root, &["gc", "--apply", "--purge-expired"]);
    success(root, &["quant", "--single-fragment-length", "180"]);
    let new_reads = fastqs(&db);
    assert_eq!(new_reads.len(), 2);
    success(root, &["quant", "--single-fragment-length", "200"]);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    assert!(
        db.generation_input_binding(OutputGroup::Quant, "RUN1")
            .unwrap()
            .is_some()
    );
    let report = files(root, &["clean", "--mode", "fastq", "--apply"]);
    let fresh = report["entries"]
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row["action"] == "trash")
        .unwrap();
    assert!(db.verify_generation_head(&counts).unwrap());
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    files(root, &["restore", fresh["trash_id"].as_str().unwrap()]);
    assert!(new_reads.iter().all(|row| Path::new(&row.path).exists()));
    assert_eq!(fastqs(&db).len(), 2);
    let preview = command(root, &["storage", "generations", "gc"]);
    assert!(
        preview.status.success(),
        "{}",
        String::from_utf8_lossy(&preview.stderr)
    );
}

#[test]
fn finished_download_queue_exits_when_the_last_file_crosses_the_storage_threshold() {
    let source = fixture::Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "both");
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.run_accession = "RUN2".into();
    run.state = rnaa_core::RunState::Resolved;
    for file in &mut run.remote_files {
        file.url = file.url.replace("RUN1", "RUN2");
        file.basename = file.basename.replace("RUN1", "RUN2");
    }
    db.upsert_run(&run).unwrap();
    let directories = [
        "raw",
        "quant",
        "de",
        "corr",
        "refs",
        "metadata",
        ".rnaa_support",
        ".rnaa_generations",
        ".rnaa_transfers",
    ];
    let usage = || {
        directories
            .iter()
            .map(|name| rnaa_core::util::dir_size(&root.join(name)))
            .sum::<u64>()
    };
    let limit = usage() + 1024;
    config(root, &db, |config| {
        config.storage.max_active_size = format!("{limit}B")
    });
    let mut child = std::process::Command::new(env!("CARGO_BIN_EXE_rnaa"))
        .arg("--root")
        .arg(root)
        .args(["download", "--concurrency", "1"])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    let start = std::time::Instant::now();
    while child.try_wait().unwrap().is_none() {
        if start.elapsed() > std::time::Duration::from_secs(15) {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("finished queue stalled at its storage threshold");
        }
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(usage() > limit);
    assert_eq!(
        db.get_run("RUN2").unwrap().unwrap().state,
        rnaa_core::RunState::Verified
    );
}

#[test]
fn idle_forever_downloader_releases_project_pins_so_completed_run_cleanup_can_proceed() {
    struct Child(std::process::Child);
    impl Drop for Child {
        fn drop(&mut self) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
    let source = fixture::Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "both");
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.run_accession = "RUN2".into();
    run.state = rnaa_core::RunState::Resolved;
    for file in &mut run.remote_files {
        file.url = file.url.replace("RUN1", "RUN2");
        file.basename = file.basename.replace("RUN1", "RUN2");
    }
    db.upsert_run(&run).unwrap();
    let mut child = Child(
        std::process::Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args(["download", "--concurrency", "1", "--forever"])
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .unwrap(),
    );
    let start = std::time::Instant::now();
    loop {
        assert!(
            child.0.try_wait().unwrap().is_none(),
            "background downloader exited"
        );
        assert!(
            start.elapsed() < std::time::Duration::from_secs(15),
            "background downloader did not reach idle maintenance"
        );
        if db.get_run("RUN2").unwrap().unwrap().state == rnaa_core::RunState::Verified {
            let preview = db
                .cleanup_files(rnaa_core::generations::FileCleanupOptions {
                    mode: rnaa_core::CleanupMode::None,
                    run: None,
                    apply: false,
                    purge_expired: false,
                    trash_days: 7,
                    expected_policy: None,
                })
                .unwrap();
            if !preview.busy {
                break;
            }
        }
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
    files(
        root,
        &["clean", "--mode", "fastq", "--run", "RUN1", "--apply"],
    );
    assert!(fastqs(&db).is_empty());
    assert!(child.0.try_wait().unwrap().is_none());
    assert!(
        db.list_artifacts_for_run(Some("RUN2"))
            .unwrap()
            .iter()
            .any(|row| row.kind == ArtifactKind::FastqSingle && Path::new(&row.path).is_file())
    );
}
