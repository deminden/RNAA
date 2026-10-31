//! Deterministic native download origin and project fixtures.
use super::commands::{command, success};
use rnaa_core::{
    Database, LibraryLayout, ProjectConfig, ProjectPaths, RemoteFile, RemoteFileKind, RunRecord,
    RunState,
    generations::{GenerationHead, OutputGroup},
    util::{compute_md5, now_rfc3339},
};
use std::{
    collections::BTreeMap,
    fs,
    io::{BufRead, BufReader, Write},
    net::TcpListener,
    path::Path,
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, Ordering},
    },
    thread,
    time::Duration,
};
#[derive(Clone)]
struct Origin {
    body: Vec<u8>,
    etag: Option<String>,
    status: u16,
}
pub struct Server {
    pub base: String,
    origins: Arc<Mutex<BTreeMap<String, Origin>>>,
    calls: Arc<Mutex<Vec<String>>>,
    pub offline: Arc<AtomicBool>,
    pub head: Arc<AtomicBool>,
    stop: Arc<AtomicBool>,
    handle: Option<thread::JoinHandle<()>>,
}
pub fn body(id: &str, sequence: &str) -> Vec<u8> {
    format!("@{id}\n{sequence}\n+\nIIII\n").into_bytes()
}
impl Server {
    pub fn new() -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        listener.set_nonblocking(true).unwrap();
        let base = format!("http://{}", listener.local_addr().unwrap());
        let origins = Arc::new(Mutex::new(BTreeMap::<String, Origin>::new()));
        let calls = Arc::new(Mutex::new(Vec::new()));
        let offline = Arc::new(AtomicBool::new(false));
        let head = Arc::new(AtomicBool::new(true));
        let stop = Arc::new(AtomicBool::new(false));
        let (values, seen, failed, heads, stopped) = (
            origins.clone(),
            calls.clone(),
            offline.clone(),
            head.clone(),
            stop.clone(),
        );
        let handle = thread::spawn(move || {
            while !stopped.load(Ordering::SeqCst) {
                let (mut stream, _) = match listener.accept() {
                    Ok(value) => value,
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        thread::sleep(Duration::from_millis(5));
                        continue;
                    }
                    Err(error) => panic!("{error}"),
                };
                // macOS inherits the listener's nonblocking mode on accept.
                // Header reads must wait for data rather than close a new socket.
                stream.set_nonblocking(false).unwrap();
                stream
                    .set_read_timeout(Some(Duration::from_secs(10)))
                    .unwrap();
                let mut reader = BufReader::new(stream.try_clone().unwrap());
                let mut line = String::new();
                if reader.read_line(&mut line).unwrap_or_default() == 0 {
                    continue;
                }
                let parts = line
                    .split_whitespace()
                    .map(str::to_string)
                    .collect::<Vec<_>>();
                let method = parts[0].clone();
                let path = parts[1].clone();
                let mut conditional = None;
                loop {
                    line.clear();
                    if reader.read_line(&mut line).unwrap_or_default() == 0 || line == "\r\n" {
                        break;
                    }
                    if let Some((name, value)) = line.split_once(':')
                        && name.eq_ignore_ascii_case("if-none-match")
                    {
                        conditional = Some(value.trim().to_string());
                    }
                }
                seen.lock().unwrap().push(format!("{method} {path}"));
                let value = values
                    .lock()
                    .unwrap()
                    .get(&path)
                    .cloned()
                    .unwrap_or(Origin {
                        body: vec![],
                        etag: None,
                        status: 404,
                    });
                let status = if failed.load(Ordering::SeqCst) {
                    503
                } else if method == "HEAD" && !heads.load(Ordering::SeqCst) {
                    405
                } else if method == "HEAD"
                    && value.status == 200
                    && conditional
                        .as_ref()
                        .is_some_and(|tag| value.etag.as_ref() == Some(tag))
                {
                    304
                } else {
                    value.status
                };
                let reason = match status {
                    200 => "OK",
                    304 => "Not Modified",
                    405 => "Method Not Allowed",
                    503 => "Unavailable",
                    _ => "Not Found",
                };
                let length = if status == 200 || status == 304 {
                    value.body.len()
                } else {
                    0
                };
                let mut response = format!(
                    "HTTP/1.1 {status} {reason}\r\nContent-Length: {length}\r\nConnection: close\r\n"
                );
                if let Some(tag) = value.etag {
                    response.push_str(&format!("ETag: {tag}\r\n"));
                }
                response.push_str("\r\n");
                let _ = stream.write_all(response.as_bytes());
                if method != "HEAD" && status == 200 {
                    let _ = stream.write_all(&value.body);
                }
            }
        });
        let result = Self {
            base,
            origins,
            calls,
            offline,
            head,
            stop,
            handle: Some(handle),
        };
        result.origin("/a.fastq", body("a", "ACGT"), Some("\"a1\""), 200);
        result.origin("/b.fastq", body("b", "TGCA"), Some("\"b1\""), 200);
        result
    }
    pub fn origin(&self, path: &str, body: Vec<u8>, etag: Option<&str>, status: u16) {
        self.origins.lock().unwrap().insert(
            path.into(),
            Origin {
                body,
                etag: etag.map(str::to_string),
                status,
            },
        );
    }
    pub fn count(&self, method: &str) -> usize {
        self.calls
            .lock()
            .unwrap()
            .iter()
            .filter(|call| call.starts_with(&format!("{method} ")))
            .count()
    }
    pub fn remote(&self, root: &Path, path: &str, name: &str, publisher: bool) -> RemoteFile {
        let value = self.origins.lock().unwrap()[path].clone();
        let checksum = if publisher {
            let file = root.join("checksum-input");
            fs::write(&file, &value.body).unwrap();
            let md5 = compute_md5(&file).unwrap();
            fs::remove_file(file).unwrap();
            Some(md5)
        } else {
            None
        };
        RemoteFile {
            url: format!("{}{path}", self.base),
            md5: checksum,
            bytes: Some(value.body.len() as u64),
            kind: RemoteFileKind::Fastq,
            basename: name.into(),
        }
    }
}
impl Drop for Server {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::SeqCst);
        self.handle.take().unwrap().join().unwrap();
    }
}
pub fn setup(root: &Path) -> (ProjectPaths, Database) {
    success(root, &["init"]);
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.pool_root = root.join("resources").display().to_string();
    config.download.retries = 0;
    config.download.read_timeout_seconds = 10;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    (paths, db)
}
pub fn choose(db: &Database, files: Vec<RemoteFile>) {
    if let Some(mut run) = db.get_run("RUN1").unwrap() {
        run.remote_files = files;
        db.upsert_run(&run).unwrap();
    } else {
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
            state: RunState::Resolved,
            last_error: None,
            updated_at: now_rfc3339(),
        })
        .unwrap();
    }
}
pub fn head(db: &Database) -> GenerationHead {
    db.generation_head(OutputGroup::Download, "RUN1")
        .unwrap()
        .unwrap()
}
pub fn download(root: &Path) {
    let result = command(root, &["download", "--concurrency", "1"]);
    assert!(
        result.status.success(),
        "{}\n{}\n{}\n{:?}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr),
        fs::read_to_string(root.join("logs/download.log")).unwrap_or_default(),
        Database::open(root).unwrap().list_job_attempts(20).unwrap(),
    );
}
