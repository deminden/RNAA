use super::*;
mod admission;
use std::{
    io::{BufRead, BufReader},
    net::TcpListener,
    sync::atomic::{AtomicUsize, Ordering},
    thread,
};

struct Reply {
    status: &'static str,
    headers: String,
    body: Vec<u8>,
    expect: Vec<&'static str>,
    absent: Vec<&'static str>,
    resume: bool,
}
fn reply(status: &'static str, headers: &str, body: &[u8]) -> Reply {
    Reply {
        status,
        headers: headers.into(),
        body: body.into(),
        expect: vec![],
        absent: vec![],
        resume: false,
    }
}
fn server(replies: Vec<Reply>) -> (String, thread::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}/file", listener.local_addr().unwrap());
    listener.set_nonblocking(true).unwrap();
    let handle = thread::spawn(move || {
        for mut reply in replies {
            let until = Instant::now() + Duration::from_secs(8);
            let mut stream = loop {
                match listener.accept() {
                    Ok((stream, _)) => break stream,
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        assert!(
                            Instant::now() < until,
                            "HTTP test never received its request"
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
            let mut headers = String::new();
            let mut input = BufReader::new(&mut stream);
            loop {
                let mut line = String::new();
                assert!(input.read_line(&mut line).unwrap() > 0);
                if line == "\r\n" {
                    break;
                }
                headers.push_str(&line);
            }
            let headers = headers.to_ascii_lowercase();
            if reply.resume {
                let start: usize = headers
                    .lines()
                    .find_map(|line| line.strip_prefix("range: bytes="))
                    .unwrap()
                    .trim()
                    .trim_end_matches('-')
                    .parse()
                    .unwrap();
                assert!(start > 0 && start < reply.body.len());
                reply.headers = format!(
                    "Content-Length: {}\r\nContent-Range: bytes {}-{}/{}\r\nETag: \"v1\"\r\n",
                    reply.body.len() - start,
                    start,
                    reply.body.len() - 1,
                    reply.body.len()
                );
                reply.body = reply.body[start..].to_vec();
            }
            assert!(headers.contains("accept-encoding: identity"));
            for required in reply.expect {
                assert!(headers.contains(required), "missing {required}: {headers}");
            }
            for forbidden in reply.absent {
                assert!(
                    !headers.contains(forbidden),
                    "unexpected {forbidden}: {headers}"
                );
            }
            write!(
                stream,
                "HTTP/1.1 {}\r\nConnection: close\r\n{}\r\n",
                reply.status, reply.headers
            )
            .unwrap();
            // A cancellation test may deliberately close while data is sent.
            let _ = stream.write_all(&reply.body);
        }
    });
    (url, handle)
}
fn config() -> DownloadConfig {
    DownloadConfig {
        retries: 0,
        backoff_seconds: vec![0],
        read_timeout_seconds: 5,
        max_transfer_seconds: 10,
        ..Default::default()
    }
}
fn prefix(url: &str, path: &Path, bytes: &[u8], total: u64, etag: Option<&str>) {
    fs::write(path, bytes).unwrap();
    let hash = hash_prefix(path, bytes.len() as u64).unwrap();
    write_json_pretty(
        &sidecar(path),
        &Checkpoint {
            schema_version: 1,
            url: url.into(),
            expected_md5: None,
            expected_bytes: Some(total),
            effective_url: url.into(),
            etag: etag.map(str::to_owned),
            last_modified: None,
            total_bytes: Some(total),
            committed_bytes: bytes.len() as u64,
            prefix_sha256: digest(&hash),
        },
    )
    .unwrap();
}

#[test]
fn resumes_exact_ranges_and_restarts_when_the_server_sends_a_full_representation() -> Result<()> {
    for full in [false, true] {
        let mut response = if full {
            reply(
                "200 OK",
                "Content-Length: 8\r\nETag: \"v2\"\r\n",
                b"newbytes",
            )
        } else {
            reply(
                "206 Partial Content",
                "Content-Length: 4\r\nContent-Range: bytes 4-7/8\r\nETag: \"v1\"\r\n",
                b"efgh",
            )
        };
        response.expect = vec!["range: bytes=4-", "if-range: \"v1\""];
        let (url, server) = server(vec![response]);
        let temp = tempfile::tempdir()?;
        let part = temp.path().join("part");
        prefix(&url, &part, b"abcd", 8, Some("\"v1\""));
        download(&url, &part, None, Some(8), &config(), &|| Ok(()))?;
        assert_eq!(
            fs::read(&part)?,
            if full { b"newbytes" } else { b"abcdefgh" }
        );
        server.join().unwrap();
    }
    Ok(())
}

#[test]
fn malformed_partial_headers_never_append_to_the_verified_prefix() -> Result<()> {
    for headers in [
        "Content-Length: 4\r\nETag: \"v1\"\r\n",
        "Content-Length: 4\r\nContent-Range: bytes 3-6/8\r\nETag: \"v1\"\r\n",
        "Content-Length: 3\r\nContent-Range: bytes 4-7/8\r\nETag: \"v1\"\r\n",
        "Content-Length: 4\r\nContent-Range: bytes 4-7/*\r\nETag: \"v1\"\r\n",
    ] {
        let (url, server) = server(vec![reply("206 Partial Content", headers, b"efgh")]);
        let temp = tempfile::tempdir()?;
        let part = temp.path().join("part");
        prefix(&url, &part, b"abcd", 8, Some("\"v1\""));
        assert!(download(&url, &part, None, Some(8), &config(), &|| Ok(())).is_err());
        assert_eq!(fs::read(&part)?, b"abcd");
        server.join().unwrap();
    }
    Ok(())
}

#[test]
fn changed_entities_and_inconsistent_totals_restart_instead_of_mixing_versions() -> Result<()> {
    for bad in [
        "Content-Length: 4\r\nContent-Range: bytes 4-7/8\r\nETag: \"v2\"\r\n",
        "Content-Length: 4\r\nContent-Range: bytes 4-7/9\r\nETag: \"v1\"\r\n",
    ] {
        let mut full = reply(
            "200 OK",
            "Content-Length: 8\r\nETag: \"v2\"\r\n",
            b"newbytes",
        );
        full.absent = vec!["range:"];
        let (url, server) = server(vec![reply("206 Partial Content", bad, b"efgh"), full]);
        let temp = tempfile::tempdir()?;
        let part = temp.path().join("part");
        prefix(&url, &part, b"abcd", 8, Some("\"v1\""));
        let mut config = config();
        config.retries = 1;
        download(&url, &part, None, Some(8), &config, &|| Ok(()))?;
        assert_eq!(fs::read(&part)?, b"newbytes");
        server.join().unwrap();
    }
    Ok(())
}

#[test]
fn truncated_bodies_resume_from_the_durable_checkpoint() -> Result<()> {
    let mut resumed = reply(
        "206 Partial Content",
        "Content-Length: 4\r\nContent-Range: bytes 4-7/8\r\nETag: \"v1\"\r\n",
        b"efgh",
    );
    resumed.expect = vec!["range: bytes=4-", "if-range: \"v1\""];
    let (url, server) = server(vec![
        reply("200 OK", "Content-Length: 8\r\nETag: \"v1\"\r\n", b"abcd"),
        resumed,
    ]);
    let temp = tempfile::tempdir()?;
    let part = temp.path().join("part");
    let mut config = config();
    config.retries = 1;
    download(&url, &part, None, Some(8), &config, &|| Ok(()))?;
    assert_eq!(fs::read(&part)?, b"abcdefgh");
    server.join().unwrap();
    Ok(())
}

#[test]
fn accepts_only_proven_complete_416_responses() -> Result<()> {
    for valid in [true, false] {
        let headers = if valid {
            "Content-Range: bytes */8\r\nETag: \"v1\"\r\nContent-Length: 0\r\n"
        } else {
            "Content-Range: bytes */9\r\nETag: \"v1\"\r\nContent-Length: 0\r\n"
        };
        let (url, server) = server(vec![reply("416 Range Not Satisfiable", headers, b"")]);
        let temp = tempfile::tempdir()?;
        let part = temp.path().join("part");
        prefix(&url, &part, b"abcdefgh", 8, Some("\"v1\""));
        let result = download(&url, &part, None, Some(8), &config(), &|| Ok(()));
        assert_eq!(result.is_ok(), valid);
        assert_eq!(file_len(&part), if valid { 8 } else { 0 });
        server.join().unwrap();
    }
    Ok(())
}
fn file_len(path: &Path) -> u64 {
    fs::metadata(path).unwrap().len()
}

#[test]
fn corrupt_or_unvalidated_prefixes_are_discarded_and_uncommitted_tails_are_truncated() -> Result<()>
{
    for mode in ["corrupt", "unvalidated", "tail", "different-url"] {
        let mut response = if mode == "tail" {
            reply(
                "206 Partial Content",
                "Content-Length: 4\r\nContent-Range: bytes 4-7/8\r\nETag: \"v1\"\r\n",
                b"efgh",
            )
        } else {
            reply(
                "200 OK",
                "Content-Length: 8\r\nETag: \"v1\"\r\n",
                b"abcdefgh",
            )
        };
        if mode == "tail" {
            response.expect = vec!["range: bytes=4-"];
        } else {
            response.absent = vec!["range:"];
        }
        let (url, server) = server(vec![response]);
        let temp = tempfile::tempdir()?;
        let part = temp.path().join("part");
        prefix(
            if mode == "different-url" {
                "http://example.invalid/old"
            } else {
                &url
            },
            &part,
            b"abcd",
            8,
            if mode == "unvalidated" {
                None
            } else {
                Some("\"v1\"")
            },
        );
        if mode == "corrupt" {
            fs::write(&part, b"BAD!")?;
        }
        if mode == "tail" {
            OpenOptions::new()
                .append(true)
                .open(&part)?
                .write_all(b"uncommitted")?;
        }
        download(&url, &part, None, Some(8), &config(), &|| Ok(()))?;
        assert_eq!(fs::read(&part)?, b"abcdefgh");
        server.join().unwrap();
    }
    Ok(())
}

#[test]
fn short_valid_ranges_make_progress_and_authoritative_md5_is_enforced() -> Result<()> {
    let temp = tempfile::tempdir()?;
    let part = temp.path().join("part");
    let expected = temp.path().join("expected");
    fs::write(&expected, b"abcdefgh")?;
    let md5 = compute_md5(&expected)?.to_uppercase();
    let (url, server) = server(vec![
        reply(
            "206 Partial Content",
            "Content-Length: 4\r\nContent-Range: bytes 0-3/8\r\nETag: \"v1\"\r\n",
            b"abcd",
        ),
        reply(
            "206 Partial Content",
            "Content-Length: 4\r\nContent-Range: bytes 4-7/8\r\nETag: \"v1\"\r\n",
            b"efgh",
        ),
    ]);
    download(&url, &part, Some(&md5), Some(8), &config(), &|| Ok(()))?;
    assert_eq!(fs::read(&part)?, b"abcdefgh");
    server.join().unwrap();
    let (url, server) = self::server(vec![reply("200 OK", "Content-Length: 8\r\n", b"badbytes")]);
    assert!(download(&url, &part, Some(&md5), Some(8), &config(), &|| Ok(())).is_err());
    assert_eq!(file_len(&part), 0);
    server.join().unwrap();
    Ok(())
}

#[test]
fn encoded_http_representations_are_rejected_instead_of_silently_decompressed() -> Result<()> {
    let (url, server) = server(vec![reply(
        "200 OK",
        "Content-Length: 4\r\nContent-Encoding: gzip\r\n",
        b"fake",
    )]);
    let temp = tempfile::tempdir()?;
    assert!(
        download(
            &url,
            &temp.path().join("part"),
            None,
            Some(4),
            &config(),
            &|| Ok(())
        )
        .is_err()
    );
    server.join().unwrap();
    Ok(())
}

#[test]
fn cancellation_is_not_retried() -> Result<()> {
    let (url, server) = server(vec![reply("200 OK", "Content-Length: 8\r\n", b"abcdefgh")]);
    let temp = tempfile::tempdir()?;
    let calls = AtomicUsize::new(0);
    let mut config = config();
    config.retries = 5;
    let result = download(
        &url,
        &temp.path().join("part"),
        None,
        Some(8),
        &config,
        &|| {
            if calls.fetch_add(1, Ordering::SeqCst) >= 3 {
                bail!("cancelled");
            }
            Ok(())
        },
    );
    assert!(result.is_err());
    server.join().unwrap();
    Ok(())
}

struct KillChild(std::process::Child);
impl Drop for KillChild {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
const CRASH_BYTES: usize = 12 * 1024 * 1024;
#[test]
fn crash_worker() -> Result<()> {
    let Some(part) = std::env::var_os("RNAA_HTTP_CRASH_PART") else {
        return Ok(());
    };
    let url = std::env::var("RNAA_HTTP_CRASH_URL")?;
    let mut config = config();
    config.read_timeout_seconds = 10;
    config.max_transfer_seconds = 30;
    download(
        &url,
        Path::new(&part),
        None,
        Some(CRASH_BYTES as u64),
        &config,
        &|| Ok(()),
    )
}
#[test]
fn killed_http_process_resumes_only_its_durable_prefix() -> Result<()> {
    use std::sync::mpsc;
    let listener = TcpListener::bind("127.0.0.1:0")?;
    let url = format!("http://{}/file", listener.local_addr()?);
    let (ready, received) = mpsc::channel();
    let (resume, proceed) = mpsc::channel();
    let server = thread::spawn(move || {
        let body = vec![b'q'; CRASH_BYTES];
        let (mut stream, _) = listener.accept().unwrap();
        stream
            .set_read_timeout(Some(Duration::from_secs(10)))
            .unwrap();
        let mut input = BufReader::new(&mut stream);
        loop {
            let mut line = String::new();
            assert!(input.read_line(&mut line).unwrap() > 0);
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
        stream.write_all(&body[..9 * 1024 * 1024]).unwrap();
        ready.send(()).unwrap();
        proceed.recv_timeout(Duration::from_secs(15)).unwrap();
        drop(stream);
        let (mut stream, _) = listener.accept().unwrap();
        stream
            .set_read_timeout(Some(Duration::from_secs(10)))
            .unwrap();
        let mut input = BufReader::new(&mut stream);
        let mut offset = None;
        loop {
            let mut line = String::new();
            assert!(input.read_line(&mut line).unwrap() > 0);
            if line == "\r\n" {
                break;
            }
            if let Some(range) = line.to_ascii_lowercase().strip_prefix("range: bytes=") {
                offset = range.trim().trim_end_matches('-').parse::<usize>().ok();
            }
        }
        let offset = offset.expect("retry must resume a durable checkpoint");
        assert!((8 * 1024 * 1024..=9 * 1024 * 1024).contains(&offset));
        write!(stream,"HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {}-{}/{}\r\nETag: \"v1\"\r\nConnection: close\r\n\r\n",body.len()-offset,offset,body.len()-1,body.len()).unwrap();
        stream.write_all(&body[offset..]).unwrap();
    });
    let temp = tempfile::tempdir()?;
    let part = temp.path().join("part");
    let mut child = KillChild(
        std::process::Command::new(std::env::current_exe()?)
            .args([
                "--exact",
                "http_transfer::tests::crash_worker",
                "--nocapture",
            ])
            .env("RNAA_HTTP_CRASH_URL", &url)
            .env("RNAA_HTTP_CRASH_PART", &part)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .spawn()?,
    );
    received.recv_timeout(Duration::from_secs(15))?;
    let until = Instant::now() + Duration::from_secs(10);
    loop {
        assert!(
            child.0.try_wait()?.is_none(),
            "HTTP child exited before checkpoint"
        );
        if fs::read(sidecar(&part))
            .ok()
            .and_then(|bytes| serde_json::from_slice::<Checkpoint>(&bytes).ok())
            .is_some_and(|saved| saved.committed_bytes >= 8 * 1024 * 1024)
        {
            break;
        }
        assert!(
            Instant::now() < until,
            "HTTP child never recorded its prefix"
        );
        thread::sleep(Duration::from_millis(10));
    }
    child.0.kill()?;
    child.0.wait()?;
    resume.send(())?;
    let mut config = config();
    config.read_timeout_seconds = 5;
    config.max_transfer_seconds = 20;
    download(
        &url,
        &part,
        None,
        Some(CRASH_BYTES as u64),
        &config,
        &|| Ok(()),
    )?;
    let result = fs::read(&part)?;
    assert_eq!(result.len(), CRASH_BYTES);
    assert!(result.iter().all(|byte| *byte == b'q'));
    server.join().unwrap();
    Ok(())
}

#[test]
fn stalled_bodies_and_continuously_slow_transfers_have_finite_limits() -> Result<()> {
    for trickle in [false, true] {
        let listener = TcpListener::bind("127.0.0.1:0")?;
        let url = format!("http://{}/file", listener.local_addr()?);
        let server = thread::spawn(move || {
            let (mut stream, _) = listener.accept().unwrap();
            stream
                .set_read_timeout(Some(Duration::from_secs(5)))
                .unwrap();
            let mut input = BufReader::new(&mut stream);
            loop {
                let mut line = String::new();
                if input.read_line(&mut line).unwrap() == 0 || line == "\r\n" {
                    break;
                }
            }
            write!(stream, "HTTP/1.1 200 OK\r\nContent-Length: 100\r\nETag: \"v1\"\r\nConnection: close\r\n\r\n").unwrap();
            if trickle {
                for _ in 0..20 {
                    if stream.write_all(b"x").is_err() {
                        break;
                    }
                    thread::sleep(Duration::from_millis(100));
                }
            } else {
                thread::sleep(Duration::from_secs(2));
            }
        });
        let temp = tempfile::tempdir()?;
        let mut config = config();
        config.max_transfer_seconds = 1;
        let started = Instant::now();
        assert!(
            download(
                &url,
                &temp.path().join("part"),
                None,
                Some(100),
                &config,
                &|| Ok(())
            )
            .is_err()
        );
        assert!(started.elapsed() < Duration::from_secs(3));
        server.join().unwrap();
    }
    Ok(())
}
