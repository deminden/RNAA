//! Download history requires publisher identity or current strong origin evidence.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/download.rs"]
mod download_support;
use commands::command;
use download_support::{Server, body, choose, download, head, setup};
use rnaa_core::{
    Database, RemoteFile,
    generations::{GenerationHead, OutputGroup},
};
use std::{fs, io::Write, path::Path, sync::atomic::Ordering, thread, time::Duration};
#[test]
fn origin_waits_for_late_request_bytes_on_an_accepted_nonblocking_listener() {
    use std::io::Read;
    let server = Server::new();
    let mut stream =
        std::net::TcpStream::connect(server.base.strip_prefix("http://").unwrap()).unwrap();
    stream
        .set_read_timeout(Some(Duration::from_secs(2)))
        .unwrap();
    thread::sleep(Duration::from_millis(100));
    stream.write_all(b"GET /a.fastq HTTP/1.1\r\n").unwrap();
    thread::sleep(Duration::from_millis(100));
    stream
        .write_all(b"Host: localhost\r\nConnection: close\r\n\r\n")
        .unwrap();
    let mut response = Vec::new();
    stream.read_to_end(&mut response).unwrap();
    assert!(response.starts_with(b"HTTP/1.1 200 OK\r\n"));
    assert!(response.ends_with(&body("a", "ACGT")));
    assert_eq!(server.count("GET"), 1);
}
fn original(
    root: &Path,
    server: &Server,
    db: &Database,
    publisher: bool,
) -> (RemoteFile, GenerationHead) {
    let remote = server.remote(root, "/a.fastq", "reads.fastq", publisher);
    choose(db, vec![remote.clone()]);
    download(root);
    (remote, head(db))
}
fn switch(root: &Path, server: &Server, db: &Database, publisher: bool) -> GenerationHead {
    choose(
        db,
        vec![server.remote(root, "/b.fastq", "reads.fastq", publisher)],
    );
    download(root);
    head(db)
}
#[test]
fn checksum_history_reattaches_offline_without_changing_files_or_producer() {
    let server = Server::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root);
    let (remote, first) = original(root, &server, &db, true);
    let seal = first.payload_path.parent().unwrap().join("seal.json");
    let bytes = fs::read(&seal).unwrap();
    let second = switch(root, &server, &db, true);
    server.offline.store(true, Ordering::SeqCst);
    choose(&db, vec![remote]);
    download(root);
    assert_eq!(head(&db), first);
    assert_eq!(server.count("GET"), 2);
    assert_eq!(server.count("HEAD"), 0);
    assert_eq!(fs::read(seal).unwrap(), bytes);
    assert!(db.verify_generation_head(&second).unwrap());
    assert_eq!(
        db.committed_generation_history(OutputGroup::Download, "RUN1")
            .unwrap()
            .len(),
        2
    );
    let attempts = db.list_job_attempts(100).unwrap();
    assert_eq!(
        attempts
            .iter()
            .filter(|attempt| attempt.stage == "download")
            .count(),
        3
    );
    assert!(attempts.iter().all(|attempt| attempt.state == "SUCCEEDED"));
}
#[test]
fn unchecksummed_history_requires_a_current_matching_strong_origin_validator() {
    let server = Server::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root);
    let (remote, first) = original(root, &server, &db, false);
    switch(root, &server, &db, false);
    choose(&db, vec![remote]);
    download(root);
    assert_eq!(head(&db), first);
    assert_eq!(server.count("GET"), 2);
    assert_eq!(server.count("HEAD"), 1);
    let event = db
        .list_events(100)
        .unwrap()
        .into_iter()
        .find(|event| event.message == "download selected from verified history")
        .unwrap();
    assert_eq!(
        event.context["source_validations"][0]["kind"],
        "http_strong_etag"
    );
    assert_eq!(event.context["source_validations"][0]["status"], 304);
}
#[test]
fn changed_origin_or_unavailable_validation_does_not_reattach_old_bytes() {
    for variant in ["changed", "unversioned", "weak", "head_unsupported"] {
        let server = Server::new();
        let temp = tempfile::tempdir().unwrap();
        let root = temp.path();
        let (paths, db) = setup(root);
        if variant == "unversioned" {
            server.origin("/a.fastq", body("a", "ACGT"), None, 200);
        }
        if variant == "weak" {
            server.origin("/a.fastq", body("a", "ACGT"), Some("W/\"a1\""), 200);
        }
        let (remote, first) = original(root, &server, &db, false);
        let second = switch(root, &server, &db, false);
        if variant == "changed" {
            server.origin("/a.fastq", body("c", "AAAA"), Some("\"a2\""), 200);
        }
        if variant == "head_unsupported" {
            server.head.store(false, Ordering::SeqCst);
        }
        choose(&db, vec![remote]);
        download(root);
        let current = head(&db);
        assert_ne!(current, first, "{variant}");
        assert_ne!(current, second);
        assert_eq!(server.count("GET"), 3, "{variant}");
        let expected = if variant == "changed" {
            body("c", "AAAA")
        } else {
            body("a", "ACGT")
        };
        assert_eq!(
            fs::read(
                current
                    .paths(&paths)
                    .unwrap()
                    .raw_run_dir("RUN1")
                    .join("reads.fastq")
            )
            .unwrap(),
            expected
        );
        assert!(db.verify_generation_head(&first).unwrap());
    }
}
#[test]
fn failed_origin_check_and_fetch_preserve_current_head_then_verified_retry_selects_history() {
    let server = Server::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root);
    let (remote, first) = original(root, &server, &db, false);
    let second = switch(root, &server, &db, false);
    choose(&db, vec![remote]);
    server.offline.store(true, Ordering::SeqCst);
    assert!(
        !command(root, &["download", "--concurrency", "1"])
            .status
            .success()
    );
    assert_eq!(head(&db), second);
    assert!(db.verify_generation_head(&second).unwrap());
    server.offline.store(false, Ordering::SeqCst);
    download(root);
    assert_eq!(head(&db), first);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Download, "RUN1")
            .unwrap()
            .len(),
        2
    );
}
#[test]
fn corrupt_historical_payload_is_replaced_with_a_new_generation_and_preserved_for_diagnosis() {
    let server = Server::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let (remote, first) = original(root, &server, &db, true);
    let second = switch(root, &server, &db, true);
    let corrupt = first
        .paths(&paths)
        .unwrap()
        .raw_run_dir("RUN1")
        .join("reads.fastq");
    fs::write(&corrupt, body("c", "AAAA")).unwrap();
    choose(&db, vec![remote]);
    download(root);
    let current = head(&db);
    assert_ne!(current, first);
    assert_ne!(current, second);
    assert_eq!(current.request_sha256, first.request_sha256);
    assert_eq!(fs::read(corrupt).unwrap(), body("c", "AAAA"));
    assert!(db.verify_generation_head(&current).unwrap());
    assert!(db.verify_generation_head(&second).unwrap());
}
#[test]
fn retry_of_multifile_transfer_revalidates_cached_unchecksummed_inputs() {
    let server = Server::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    server.origin("/b.fastq", body("b", "TGCA"), Some("\"b1\""), 503);
    choose(
        &db,
        vec![
            server.remote(root, "/a.fastq", "one.fastq", false),
            server.remote(root, "/b.fastq", "two.fastq", false),
        ],
    );
    assert!(
        !command(root, &["download", "--concurrency", "1"])
            .status
            .success()
    );
    assert!(
        db.generation_head(OutputGroup::Download, "RUN1")
            .unwrap()
            .is_none()
    );
    server.origin("/a.fastq", body("c", "AAAA"), Some("\"a2\""), 200);
    server.origin("/b.fastq", body("b", "TGCA"), Some("\"b1\""), 200);
    download(root);
    let current = head(&db);
    assert_eq!(
        fs::read(
            current
                .paths(&paths)
                .unwrap()
                .raw_run_dir("RUN1")
                .join("one.fastq")
        )
        .unwrap(),
        body("c", "AAAA")
    );
    assert_eq!(server.count("HEAD"), 1);
    assert_eq!(server.count("GET"), 4);
}
