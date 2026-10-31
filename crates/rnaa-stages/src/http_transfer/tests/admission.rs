use super::*;
use rnaa_core::resources::{
    ExecutionResources, ReservationOwner, ResourceConfig, ResourceRequest, ResourceReservation,
};
use std::sync::{Mutex, atomic::AtomicBool, mpsc};
const MIB: u64 = 1024 * 1024;
fn context(root: &Path, strict: bool) -> Result<ExecutionResources> {
    ExecutionResources::open(
        &ResourceConfig {
            pool_root: root.join("pool").display().to_string(),
            cpu_threads: 2,
            memory_mb: 256,
            disk_headroom_mb: if strict { u64::MAX / MIB } else { 0 },
            wait_seconds: 5,
        },
        ReservationOwner {
            project_root: root.into(),
            attempt_id: if strict { "blocker" } else { "download" }.into(),
            stage: "download".into(),
            run_accession: Some("RUN".into()),
        },
    )
}
fn blocker(context: &ExecutionResources) -> Result<ResourceReservation> {
    context.reserve(
        "blocker",
        ResourceRequest {
            cpu_threads: 1,
            ..Default::default()
        },
        &|| Ok(()),
    )
}

#[test]
fn cancellation_during_prefix_validation_preserves_payload_and_checkpoint() -> Result<()> {
    let root = tempfile::tempdir()?;
    let path = root.path().join("cancelled.part");
    prefix("http://unused", &path, b"abcde", 10, Some("\"v1\""));
    let saved = fs::read(sidecar(&path))?;
    let calls = AtomicUsize::new(0);
    let result = load_prefix("http://unused", &path, None, Some(10), &|| {
        calls.fetch_add(1, Ordering::SeqCst);
        bail!("cancelled prefix validation");
    });
    assert!(result.is_err());
    assert!(calls.load(Ordering::SeqCst) >= 2);
    assert_eq!(fs::read(&path)?, b"abcde");
    assert_eq!(fs::read(sidecar(&path))?, saved);
    Ok(())
}

#[test]
fn admitted_resume_charges_only_future_writes_and_full_restart_readmits_reclaimed_storage()
-> Result<()> {
    for full in [false, true] {
        let body = if full {
            b"abcdefghij".as_slice()
        } else {
            b"fghij".as_slice()
        };
        let headers = if full {
            "Content-Length: 10\r\nETag: \"v1\"\r\n"
        } else {
            "Content-Length: 5\r\nContent-Range: bytes 5-9/10\r\nETag: \"v1\"\r\n"
        };
        let (url, server) = server(vec![reply(
            if full {
                "200 OK"
            } else {
                "206 Partial Content"
            },
            headers,
            body,
        )]);
        let root = tempfile::tempdir()?;
        let path = root.path().join("data.part");
        prefix(&url, &path, b"abcde", 10, Some("\"v1\""));
        let context = context(root.path(), false)?;
        let execution = download_admitted(
            &url,
            &path,
            None,
            Some(10),
            &config(),
            Some(&context),
            &|| Ok(()),
        )?
        .unwrap();
        server.join().unwrap();
        assert_eq!(fs::read(&path)?, b"abcdefghij");
        let execution = serde_json::to_value(execution)?;
        assert_eq!(
            execution["max_pending_disk_bytes"],
            4 * MIB + if full { 10 } else { 5 }
        );
        assert_eq!(
            execution["reservation"]["requested"]["memory_bytes"],
            64 * MIB
        );
        assert_eq!(
            execution["reservation"]["requested"]["disks"][0]["bytes"],
            4 * MIB
        );
        assert!(context.pool().snapshot()?.reservations.is_empty());
    }
    Ok(())
}

#[test]
fn unknown_length_growth_saves_prefix_and_does_not_retry_admission_errors_as_network_failures()
-> Result<()> {
    let data = vec![b'x'; 9 * MIB as usize];
    let mut resume = reply("206 Partial Content", "", &data);
    resume.resume = true;
    resume.expect = vec!["if-range: \"v1\""];
    let (url, server) = server(vec![reply("200 OK", "ETag: \"v1\"\r\n", &data), resume]);
    let root = tempfile::tempdir()?;
    let path = root.path().join("unknown.part");
    let context = context(root.path(), false)?;
    let strict = self::context(root.path(), true)?;
    let held = Mutex::new(None);
    let mut config = config();
    config.retries = 3;
    let error = download_admitted(&url, &path, None, None, &config, Some(&context), &|| {
        if fs::metadata(&path).is_ok_and(|m| m.len() >= 7 * MIB) {
            let mut guard = held.lock().unwrap();
            if guard.is_none() {
                *guard = Some(blocker(&strict)?);
            }
        }
        Ok(())
    })
    .unwrap_err();
    assert!(
        error.is::<rnaa_core::resources::ResourceAdmissionError>(),
        "{error:#}"
    );
    assert!(format!("{error:#}").contains("write window cannot grow"));
    let saved = receipt(&path)?;
    let committed = saved["committed_bytes"].as_u64().unwrap();
    assert!((7 * MIB..=8 * MIB).contains(&committed));
    assert_eq!(fs::metadata(&path)?.len(), committed);
    assert_eq!(context.pool().snapshot()?.reservations.len(), 1);
    drop(held.into_inner().unwrap());
    download_admitted(&url, &path, None, None, &config, Some(&context), &|| Ok(()))?;
    server.join().unwrap();
    assert_eq!(fs::read(&path)?, data);
    assert!(context.pool().snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn discovered_length_drops_response_and_all_claims_before_capacity_readmission() -> Result<()> {
    let data = vec![b'x'; 16 * MIB as usize];
    let headers = format!("Content-Length: {}\r\nETag: \"v1\"\r\n", data.len());
    let (url, server) = server(vec![
        reply("200 OK", &headers, &data),
        reply("200 OK", &headers, &data),
    ]);
    let root = tempfile::tempdir()?;
    let path = root.path().join("discovery.part");
    let context = context(root.path(), false)?;
    let strict = self::context(root.path(), true)?;
    let held = Mutex::new(None);
    let inserted = AtomicBool::new(false);
    let announced = AtomicBool::new(false);
    let (tx, rx) = mpsc::channel();
    thread::scope(|scope| -> Result<()> {
        let worker = scope.spawn(|| {
            download_admitted(&url, &path, None, None, &config(), Some(&context), &|| {
                let snapshot = context.pool().snapshot()?;
                if snapshot
                    .reservations
                    .iter()
                    .any(|g| g.owner.stage == "download.transfer")
                    && !inserted.swap(true, Ordering::SeqCst)
                {
                    *held.lock().unwrap() = Some(blocker(&strict)?);
                }
                if inserted.load(Ordering::SeqCst)
                    && snapshot.reservations.len() == 1
                    && snapshot.reservations[0].owner.stage == "download.blocker"
                    && !announced.swap(true, Ordering::SeqCst)
                {
                    tx.send(())?;
                }
                Ok(())
            })
        });
        let observed = rx.recv_timeout(Duration::from_secs(5));
        drop(held.lock().unwrap().take());
        let result = worker.join().unwrap();
        observed?;
        result?;
        Ok(())
    })?;
    server.join().unwrap();
    assert_eq!(fs::read(&path)?, data);
    assert!(context.pool().snapshot()?.reservations.is_empty());
    Ok(())
}
