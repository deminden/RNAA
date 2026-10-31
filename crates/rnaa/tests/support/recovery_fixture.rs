//! Native pipeline fixtures shared by dependency and QC recovery scenarios.
#![allow(dead_code)]
use super::commands::success;
use rnaa_core::{
    ArtifactKind, ArtifactRecord, Database, LibraryLayout, ProjectConfig, ProjectPaths, RemoteFile,
    RemoteFileKind, RunRecord, RunState,
    generations::{GenerationHead, OutputGroup},
    util::{compute_md5, compute_sha256, file_size, now_rfc3339},
};
use std::{
    collections::BTreeMap,
    fs,
    io::{BufRead, BufReader, Write},
    net::TcpListener,
    path::Path,
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, AtomicUsize, Ordering},
    },
    thread,
    time::Duration,
};
pub fn sequence() -> String {
    let mut state = 1729u64;
    (0..1024)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize] as char
        })
        .collect()
}
pub fn reads() -> Vec<u8> {
    let seq = sequence();
    (0..32)
        .map(|i| {
            format!(
                "@read{i}\n{}\n+\n{}\n",
                &seq[50 + i * 13..150 + i * 13],
                "I".repeat(100)
            )
        })
        .collect::<String>()
        .into_bytes()
}

pub fn altered_reads() -> Vec<u8> {
    let mut body = reads();
    let offset = body.iter().position(|byte| *byte == b'\n').unwrap() + 1;
    body[offset] = if body[offset] == b'A' { b'C' } else { b'A' };
    body
}

pub struct Source {
    pub url: String,
    pub failed: Arc<AtomicBool>,
    pub requests: Arc<AtomicUsize>,
    pub request_lines: Arc<Mutex<Vec<String>>>,
    pub delay_ms: Arc<AtomicUsize>,
    pub alternate: Arc<AtomicBool>,
    pub low_quality: Arc<AtomicBool>,
    stop: Arc<AtomicBool>,
    handle: Option<thread::JoinHandle<()>>,
}
impl Source {
    pub fn new() -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        listener.set_nonblocking(true).unwrap();
        let url = format!("http://{}/RUN1.fastq", listener.local_addr().unwrap());
        let failed = Arc::new(AtomicBool::new(false));
        let requests = Arc::new(AtomicUsize::new(0));
        let request_lines = Arc::new(Mutex::new(Vec::new()));
        let request_log = request_lines.clone();
        let stop = Arc::new(AtomicBool::new(false));
        let delay_ms = Arc::new(AtomicUsize::new(0));
        let alternate = Arc::new(AtomicBool::new(false));
        let low_quality = Arc::new(AtomicBool::new(false));
        let reject = low_quality.clone();
        let changed = alternate.clone();
        let (fail, seen, stopped, delay) = (
            failed.clone(),
            requests.clone(),
            stop.clone(),
            delay_ms.clone(),
        );
        let handle = thread::spawn(move || {
            'connections: while !stopped.load(Ordering::SeqCst) {
                let (mut stream, _) = match listener.accept() {
                    Ok(value) => value,
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        thread::sleep(Duration::from_millis(5));
                        continue;
                    }
                    Err(error) => panic!("HTTP fixture: {error}"),
                };
                // Accepted sockets inherit nonblocking mode on macOS. Wait for
                // actual HTTP bytes instead of counting an empty connection.
                stream.set_nonblocking(false).unwrap();
                stream
                    .set_read_timeout(Some(Duration::from_secs(10)))
                    .unwrap();
                let mut reader = BufReader::new(stream.try_clone().unwrap());
                let mut request_line = String::new();
                if reader.read_line(&mut request_line).unwrap_or_default() == 0 {
                    continue;
                }
                loop {
                    let mut line = String::new();
                    match reader.read_line(&mut line) {
                        Ok(0) | Err(_) => continue 'connections,
                        Ok(_) if line == "\r\n" => break,
                        Ok(_) => {}
                    }
                }
                request_log
                    .lock()
                    .unwrap()
                    .push(request_line.trim_end().to_string());
                seen.fetch_add(1, Ordering::SeqCst);
                thread::sleep(Duration::from_millis(delay.load(Ordering::SeqCst) as u64));
                let alternate = changed.load(Ordering::SeqCst);
                let mut body = if alternate { altered_reads() } else { reads() };
                if reject.load(Ordering::SeqCst) {
                    for byte in &mut body {
                        if *byte == b'I' {
                            *byte = b'!';
                        }
                    }
                }
                if fail.load(Ordering::SeqCst) {
                    body[0] = b'X';
                }
                let etag = if alternate {
                    "fixture-v2"
                } else {
                    "fixture-v1"
                };
                let response = format!(
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nETag: \"{etag}\"\r\nConnection: close\r\n\r\n",
                    body.len()
                );
                let _ = stream.write_all(response.as_bytes());
                let _ = stream.write_all(&body);
            }
        });
        Self {
            url,
            failed,
            requests,
            request_lines,
            delay_ms,
            alternate,
            low_quality,
            stop,
            handle: Some(handle),
        }
    }
}
impl Drop for Source {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::SeqCst);
        self.handle.take().unwrap().join().unwrap();
    }
}

pub fn setup(root: &Path, source: Option<&Source>, retention: &str) -> (ProjectPaths, Database) {
    setup_with_checksum(root, source, retention, true)
}

pub fn setup_with_checksum(
    root: &Path,
    source: Option<&Source>,
    retention: &str,
    publisher_checksum: bool,
) -> (ProjectPaths, Database) {
    success(root, &["init"]);
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).unwrap();
    let cdna = root.join("source.fa");
    let gtf = root.join("source.gtf");
    fs::write(&cdna, format!(">tx1\n{}\n", sequence())).unwrap();
    fs::write(&gtf, "chr1\tsynthetic\ttranscript\t1\t1024\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\"; gene_name \"Gene1\";\n").unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.refs.custom_cdna = cdna.display().to_string();
    config.refs.custom_gtf = gtf.display().to_string();
    config.quant.preprocess = true;
    config.quant.preprocess_threads = 1;
    config.quant.worker_budget = 1;
    config.quant.threads = 1;
    config.download.retries = 1;
    config.download.backoff_seconds = vec![0];
    // Native pipelines run concurrently in this suite; keep scheduling delays
    // from becoming transport retries. Dedicated HTTP tests cover finite stalls.
    config.download.read_timeout_seconds = 10;
    config.download.max_transfer_seconds = 60;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let raw = root.join("checksum-input.fastq");
    fs::write(&raw, reads()).unwrap();
    let md5 = compute_md5(&raw).unwrap();
    fs::remove_file(&raw).unwrap();
    let remotes = source
        .map(|source| {
            vec![
                RemoteFile {
                    url: source.url.clone(),
                    md5: publisher_checksum.then_some(md5),
                    bytes: Some(reads().len() as u64),
                    kind: RemoteFileKind::Fastq,
                    basename: "RUN1.fastq".into(),
                },
                RemoteFile {
                    url: source.url.replace("RUN1.fastq", "RUN1.sra"),
                    md5: None,
                    bytes: None,
                    kind: RemoteFileKind::Sra,
                    basename: "RUN1.sra".into(),
                },
            ]
        })
        .unwrap_or_default();
    db.upsert_run(&RunRecord {
        project_id: db.project_id().into(),
        study_accession: "synthetic".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: Some("BIO1".into()),
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: remotes,
        metadata: BTreeMap::new(),
        state: if source.is_some() {
            RunState::Resolved
        } else {
            RunState::Verified
        },
        last_error: None,
        updated_at: now_rfc3339(),
    })
    .unwrap();
    if source.is_some() {
        success(root, &["download"]);
    } else {
        let dir = paths.raw_run_dir("RUN1");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("RUN1.fastq");
        fs::write(&path, reads()).unwrap();
        db.record_artifact(&ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some("RUN1".into()),
            kind: ArtifactKind::FastqSingle,
            path: path.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".into(),
            checksum: compute_sha256(&path).unwrap(),
            bytes: file_size(&path).unwrap(),
            created_at: now_rfc3339(),
        })
        .unwrap();
    }
    success(root, &["quant", "--fastq-retention", retention]);
    (paths, db)
}
pub fn head(db: &Database, group: OutputGroup) -> GenerationHead {
    db.generation_head(group, "RUN1").unwrap().unwrap()
}
pub fn replace(root: &Path, retention: &str) {
    success(
        root,
        &[
            "quant",
            "--single-fragment-length",
            "180",
            "--fastq-retention",
            retention,
        ],
    );
}
