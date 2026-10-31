use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::util::{compute_sha256, now_rfc3339};
use rnaa_core::{
    ArtifactKind, ArtifactRecord, Database, LibraryLayout, ProjectConfig, ProjectPaths, RemoteFile,
    RemoteFileKind, RunRecord, RunState,
};
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;
use std::path::Path;
use std::process::{Child, Command, Output, Stdio};
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

fn command(root: &Path, args: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_rnaa"));
    command.arg("--root").arg(root).args(args);
    command
}
fn succeeded(output: Output) {
    assert!(
        output.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}
fn setup(root: &Path, state: RunState, files: Vec<RemoteFile>) -> (ProjectPaths, Database) {
    let paths = ProjectPaths::new(root);
    paths.ensure_layout().unwrap();
    let mut config = ProjectConfig::default();
    config.download.retries = 0;
    config.save(&paths.config_path).unwrap();
    let db = Database::create(root, &config).unwrap();
    db.upsert_run(&RunRecord {
        project_id: db.project_id().into(),
        study_accession: "study".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: files,
        metadata: BTreeMap::new(),
        state,
        last_error: None,
        updated_at: now_rfc3339(),
    })
    .unwrap();
    (paths, db)
}
fn historical_raw(paths: &ProjectPaths, db: &Database) {
    let path = paths.raw_run_dir("RUN1").join("old.fastq");
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(&path, b"@old\nACGT\n+\nIIII\n").unwrap();
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: Some("RUN1".into()),
        kind: ArtifactKind::FastqSingle,
        path: path.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".into(),
        checksum: compute_sha256(&path).unwrap(),
        bytes: fs::metadata(&path).unwrap().len(),
        created_at: now_rfc3339(),
    })
    .unwrap();
}

#[test]
fn status_preserves_a_live_worker_and_cli_cancellation_fences_its_completion() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path(), RunState::Verified, vec![]);
    historical_raw(&paths, &db);
    let mut lease = JobLease::try_acquire(&db, "quant", Some("RUN1"), &[RunState::Verified])
        .unwrap()
        .unwrap();
    succeeded(command(root.path(), &["status"]).output().unwrap());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantRunning
    );
    let listing = command(root.path(), &["jobs", "list", "--json"])
        .output()
        .unwrap();
    assert!(listing.status.success());
    let jobs: serde_json::Value = serde_json::from_slice(&listing.stdout).unwrap();
    assert_eq!(jobs[0]["attempt_id"], lease.attempt_id());
    assert_eq!(jobs[0]["state"], "RUNNING");
    succeeded(
        command(root.path(), &["jobs", "cancel", lease.attempt_id()])
            .output()
            .unwrap(),
    );
    assert!(
        lease
            .database()
            .set_run_state("RUN1", RunState::QuantDone, None)
            .is_err()
    );
    assert!(lease.finish(JobOutcome::Succeeded, None).is_err());
    drop(lease);
    assert_eq!(db.list_job_attempts(1).unwrap()[0].state, "CANCELLED");
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantFailed
    );
}

#[test]
fn finite_download_failure_is_a_failed_command_and_a_durable_failed_attempt() {
    let root = tempfile::tempdir().unwrap();
    let (_paths, db) = setup(root.path(), RunState::Resolved, vec![]);
    let output = command(root.path(), &["download", "--concurrency", "1"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("download failed for 1 requested run")
    );
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::DownloadFailed
    );
    let jobs = db.list_job_attempts(10).unwrap();
    assert_eq!(jobs.len(), 1);
    assert_eq!(jobs[0].state, "FAILED");
    assert!(jobs[0].request_sha256.is_some());
}

struct KillOnDrop(Option<Child>);
impl KillOnDrop {
    fn start(mut command: Command) -> Self {
        Self(Some(
            command
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap(),
        ))
    }
    fn finish(mut self) -> Output {
        self.0.take().unwrap().wait_with_output().unwrap()
    }
}
impl Drop for KillOnDrop {
    fn drop(&mut self) {
        if let Some(child) = &mut self.0 {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

#[test]
fn competing_downloader_processes_and_status_share_one_download_attempt() {
    let root = tempfile::tempdir().unwrap();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let address = listener.local_addr().unwrap();
    let body = b"@read\nACGT\n+\nIIII\n";
    let (paths, db) = setup(
        root.path(),
        RunState::Resolved,
        vec![RemoteFile {
            url: format!("http://{address}/reads.fastq"),
            md5: None,
            bytes: Some(body.len() as u64),
            kind: RemoteFileKind::Fastq,
            basename: "reads.fastq".into(),
        }],
    );
    historical_raw(&paths, &db);
    let (arrived, request) = mpsc::channel();
    let (release, proceed) = mpsc::channel();
    let server = thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        stream
            .set_read_timeout(Some(Duration::from_secs(10)))
            .unwrap();
        {
            let mut reader = BufReader::new(&mut stream);
            let mut line = String::new();
            assert!(reader.read_line(&mut line).unwrap() > 0);
            assert!(line.starts_with("GET /reads.fastq "));
            loop {
                line.clear();
                assert!(reader.read_line(&mut line).unwrap() > 0);
                if line == "\r\n" {
                    break;
                }
            }
        }
        arrived.send(()).unwrap();
        proceed.recv_timeout(Duration::from_secs(10)).unwrap();
        write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
            body.len()
        )
        .unwrap();
        stream.write_all(body).unwrap();
    });
    let first = KillOnDrop::start(command(root.path(), &["download", "--concurrency", "1"]));
    request.recv_timeout(Duration::from_secs(10)).unwrap();
    let second = KillOnDrop::start(command(root.path(), &["download", "--concurrency", "1"]));
    succeeded(command(root.path(), &["status"]).output().unwrap());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::Downloading
    );
    assert_eq!(db.list_job_attempts(10).unwrap().len(), 1);
    release.send(()).unwrap();
    let output = first.finish();
    assert!(
        output.status.success(),
        "{}\n{}\nattempts: {:?}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
        db.list_job_attempts(10).unwrap()
    );
    succeeded(second.finish());
    server.join().unwrap();
    let jobs = db.list_job_attempts(10).unwrap();
    assert_eq!(jobs.len(), 1);
    assert_eq!(jobs[0].state, "SUCCEEDED");
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::Verified
    );
    assert_eq!(
        fs::read(
            paths
                .with_generation_heads(&db)
                .unwrap()
                .raw_run_dir("RUN1")
                .join("reads.fastq")
        )
        .unwrap(),
        body
    );
}

struct HttpReply {
    method: &'static str,
    path: &'static str,
    status: &'static str,
    body: Vec<u8>,
}
fn serve_downloads(replies: Vec<HttpReply>) -> (String, thread::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let base = format!("http://{}", listener.local_addr().unwrap());
    listener.set_nonblocking(true).unwrap();
    let handle = thread::spawn(move || {
        for reply in replies {
            let until = std::time::Instant::now() + Duration::from_secs(15);
            let mut stream = loop {
                match listener.accept() {
                    Ok((stream, _)) => break stream,
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        assert!(
                            std::time::Instant::now() < until,
                            "missing request for {}",
                            reply.path
                        );
                        thread::sleep(Duration::from_millis(5));
                    }
                    Err(error) => panic!("{error}"),
                }
            };
            stream.set_nonblocking(false).unwrap();
            stream
                .set_read_timeout(Some(Duration::from_secs(5)))
                .unwrap();
            let mut reader = BufReader::new(&mut stream);
            let mut request = String::new();
            reader.read_line(&mut request).unwrap();
            assert!(
                request.starts_with(&format!("{} {} ", reply.method, reply.path)),
                "{request}"
            );
            loop {
                let mut line = String::new();
                assert!(reader.read_line(&mut line).unwrap() > 0);
                if line == "\r\n" {
                    break;
                }
            }
            write!(
                stream,
                "HTTP/1.1 {}\r\nContent-Length: {}\r\nETag: \"v1\"\r\nConnection: close\r\n\r\n",
                reply.status,
                reply.body.len()
            )
            .unwrap();
            if reply.method != "HEAD" {
                stream.write_all(&reply.body).unwrap();
            }
        }
    });
    (base, handle)
}
fn remote(base: &str, path: &str, body: &[u8]) -> RemoteFile {
    RemoteFile {
        url: format!("{base}{path}"),
        basename: path.trim_start_matches('/').into(),
        bytes: Some(body.len() as u64),
        md5: None,
        kind: RemoteFileKind::Fastq,
    }
}

#[test]
fn failed_multifile_download_keeps_previous_generation_and_retry_publishes_one_complete_replacement()
 {
    use rnaa_core::generations::OutputGroup;
    let root = tempfile::tempdir().unwrap();
    let one = b"@a\nACGT\n+\nIIII\n";
    let two = b"@b\nTGCA\n+\nIIII\n";
    let (base, server) = serve_downloads(vec![
        HttpReply {
            method: "GET",
            path: "/old.fastq",
            status: "200 OK",
            body: one.to_vec(),
        },
        HttpReply {
            method: "GET",
            path: "/one.fastq",
            status: "200 OK",
            body: one.to_vec(),
        },
        HttpReply {
            method: "GET",
            path: "/two.fastq",
            status: "503 Unavailable",
            body: vec![],
        },
        HttpReply {
            method: "HEAD",
            path: "/one.fastq",
            status: "304 Not Modified",
            body: one.to_vec(),
        },
        HttpReply {
            method: "GET",
            path: "/two.fastq",
            status: "200 OK",
            body: two.to_vec(),
        },
    ]);
    let (paths, db) = setup(
        root.path(),
        RunState::Resolved,
        vec![remote(&base, "/old.fastq", one)],
    );
    succeeded(
        command(root.path(), &["download", "--concurrency", "1"])
            .output()
            .unwrap(),
    );
    let previous = db
        .generation_head(OutputGroup::Download, "RUN1")
        .unwrap()
        .unwrap();
    let pinned = paths.with_generation_heads(&db).unwrap();
    let records = db.list_artifacts_for_run(Some("RUN1")).unwrap();
    let old = records
        .iter()
        .map(|record| (record.path.clone(), fs::read(&record.path).unwrap()))
        .collect::<Vec<_>>();
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.remote_files = vec![
        remote(&base, "/one.fastq", one),
        remote(&base, "/two.fastq", two),
    ];
    db.upsert_run(&run).unwrap();
    let output = command(root.path(), &["download", "--concurrency", "1"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::DownloadFailed
    );
    assert_eq!(
        db.generation_head(OutputGroup::Download, "RUN1")
            .unwrap()
            .unwrap(),
        previous
    );
    assert_eq!(db.list_artifacts_for_run(Some("RUN1")).unwrap().len(), 2);
    for (path, bytes) in &old {
        assert_eq!(fs::read(path).unwrap(), *bytes);
    }
    assert!(db.verify_generation(OutputGroup::Download, "RUN1").unwrap());
    let retry = command(root.path(), &["download", "--concurrency", "1"])
        .output()
        .unwrap();
    assert!(
        retry.status.success(),
        "{}\n{:?}\n{}",
        String::from_utf8_lossy(&retry.stderr),
        db.list_job_attempts(10).unwrap(),
        fs::read_to_string(paths.logs_dir.join("download.log")).unwrap_or_default()
    );
    server.join().unwrap();
    let next = db
        .generation_head(OutputGroup::Download, "RUN1")
        .unwrap()
        .unwrap();
    assert_ne!(next.generation_id, previous.generation_id);
    let current = paths.with_generation_heads(&db).unwrap();
    assert_eq!(
        fs::read(current.raw_run_dir("RUN1").join("one.fastq")).unwrap(),
        one
    );
    assert_eq!(
        fs::read(current.raw_run_dir("RUN1").join("two.fastq")).unwrap(),
        two
    );
    assert_eq!(
        fs::read(pinned.raw_run_dir("RUN1").join("old.fastq")).unwrap(),
        one
    );
    let records = db.list_artifacts_for_run(Some("RUN1")).unwrap();
    assert_eq!(records.len(), 3);
    assert!(
        records
            .iter()
            .all(|record| Path::new(&record.path).starts_with(&next.payload_path))
    );
    assert!(db.verify_generation(OutputGroup::Download, "RUN1").unwrap());
    let attempts = db.list_job_attempts(100).unwrap();
    assert_eq!(attempts.len(), 3);
    assert_eq!(
        attempts
            .iter()
            .find(|attempt| attempt.attempt_id == next.attempt_id)
            .unwrap()
            .state,
        "SUCCEEDED"
    );
    // The listener is gone: a matching committed generation must reuse without IO.
    db.set_run_state("RUN1", RunState::DownloadFailed, None)
        .unwrap();
    succeeded(
        command(root.path(), &["download", "--concurrency", "1"])
            .output()
            .unwrap(),
    );
    assert_eq!(
        db.generation_head(OutputGroup::Download, "RUN1")
            .unwrap()
            .unwrap(),
        next
    );
    fs::write(
        current.raw_run_dir("RUN1").join("one.fastq"),
        b"@a\nTGCA\n+\nIIII\n",
    )
    .unwrap();
    let output = command(root.path(), &["download", "--concurrency", "1"])
        .output()
        .unwrap();
    assert!(
        !output.status.success(),
        "same-size corruption cannot count as reuse"
    );
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::DownloadFailed
    );
    assert_eq!(
        db.generation_head(OutputGroup::Download, "RUN1")
            .unwrap()
            .unwrap(),
        next
    );
    for (path, bytes) in old {
        assert_eq!(fs::read(path).unwrap(), bytes);
    }
}

#[test]
fn cancellation_during_an_http_body_preserves_previous_artifacts_and_never_publishes() {
    use rnaa_core::generations::OutputGroup;
    let root = tempfile::tempdir().unwrap();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let base = format!("http://{}", listener.local_addr().unwrap());
    let body = b"@a\nACGT\n+\nIIII\n";
    let (paths, db) = setup(
        root.path(),
        RunState::Resolved,
        vec![remote(&base, "/reads.fastq", body)],
    );
    historical_raw(&paths, &db);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.download.read_timeout_seconds = 1;
    config.save(&paths.config_path).unwrap();
    let (sent, received) = mpsc::channel();
    let server = thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        stream
            .set_read_timeout(Some(Duration::from_secs(5)))
            .unwrap();
        let mut reader = BufReader::new(&mut stream);
        loop {
            let mut line = String::new();
            assert!(reader.read_line(&mut line).unwrap() > 0);
            if line == "\r\n" {
                break;
            }
        }
        write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nETag: \"v1\"\r\nConnection: close\r\n\r\n",
            body.len()
        )
        .unwrap();
        stream.write_all(&body[..4]).unwrap();
        sent.send(()).unwrap();
        thread::sleep(Duration::from_secs(2));
    });
    let process = KillOnDrop::start(command(root.path(), &["download", "--concurrency", "1"]));
    received.recv_timeout(Duration::from_secs(10)).unwrap();
    let attempt = db.list_job_attempts(1).unwrap().remove(0);
    succeeded(
        command(root.path(), &["jobs", "cancel", &attempt.attempt_id])
            .output()
            .unwrap(),
    );
    assert!(!process.finish().status.success());
    server.join().unwrap();
    assert!(
        db.generation_head(OutputGroup::Download, "RUN1")
            .unwrap()
            .is_none()
    );
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::DownloadFailed
    );
    assert_eq!(db.list_job_attempts(1).unwrap()[0].state, "CANCELLED");
    assert_eq!(
        fs::read(paths.raw_run_dir("RUN1").join("old.fastq")).unwrap(),
        b"@old\nACGT\n+\nIIII\n"
    );
    assert_eq!(db.list_artifacts_for_run(Some("RUN1")).unwrap().len(), 1);
}

#[test]
fn source_metadata_changes_during_transfer_cannot_publish_a_stale_download() {
    use rnaa_core::generations::OutputGroup;
    let root = tempfile::tempdir().unwrap();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let base = format!("http://{}", listener.local_addr().unwrap());
    let body = b"@a\nACGT\n+\nIIII\n";
    let (_paths, db) = setup(
        root.path(),
        RunState::Resolved,
        vec![remote(&base, "/reads.fastq", body)],
    );
    let (started, ready) = mpsc::channel();
    let (release, resume) = mpsc::channel();
    let server = thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        stream
            .set_read_timeout(Some(Duration::from_secs(5)))
            .unwrap();
        let mut reader = BufReader::new(&mut stream);
        loop {
            let mut line = String::new();
            assert!(reader.read_line(&mut line).unwrap() > 0);
            if line == "\r\n" {
                break;
            }
        }
        started.send(()).unwrap();
        resume.recv_timeout(Duration::from_secs(10)).unwrap();
        write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
            body.len()
        )
        .unwrap();
        stream.write_all(body).unwrap();
    });
    let process = KillOnDrop::start(command(root.path(), &["download", "--concurrency", "1"]));
    ready.recv_timeout(Duration::from_secs(10)).unwrap();
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.remote_files[0].url.push_str("?source-revision=2");
    db.upsert_run(&run).unwrap();
    release.send(()).unwrap();
    assert!(!process.finish().status.success());
    server.join().unwrap();
    assert!(
        db.generation_head(OutputGroup::Download, "RUN1")
            .unwrap()
            .is_none()
    );
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::DownloadFailed
    );
    let attempt = db.list_job_attempts(1).unwrap().remove(0);
    assert_eq!(attempt.state, "FAILED");
    assert!(
        attempt
            .last_error
            .unwrap()
            .contains("source metadata changed")
    );
    assert!(db.list_artifacts_for_run(Some("RUN1")).unwrap().is_empty());
}
