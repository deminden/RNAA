use super::*;
use rnaa_core::{ArtifactKind, LibraryLayout, VerifiedFile};
#[test]
fn estimates_grow_with_real_records_and_capped_metadata_never_writes() -> Result<()> {
    let root = tempfile::tempdir()?;
    let path = root.path().join("reads.fastq");
    std::fs::write(
        &path,
        format!("@read\n{}\n+\n{}\n", "ACGT".repeat(25), "I".repeat(100)),
    )?;
    let plan = FastqPlan::inspect(
        LibraryLayout::Single,
        &[VerifiedFile {
            artifact_kind: ArtifactKind::FastqSingle,
            bytes: std::fs::metadata(&path)?.len(),
            integrity_source: "fixture".into(),
            path,
            checksum_type: "none".into(),
            checksum: String::new(),
        }],
        FastqLimits::default(),
        &|| Ok(()),
    )?;
    let small = Estimate::for_plan(&plan, NativeProcessOptions::for_threads(1))?;
    let parallel = Estimate::for_plan(&plan, NativeProcessOptions::for_threads(4))?;
    assert!(parallel.memory_bytes > small.memory_bytes);
    assert_eq!(small.staging_bytes, 0);
    assert_eq!(small.read_output_bytes, plan.summary().normalized_bytes()?);
    let out = root.path().join("report.json");
    let error = write_metadata(&out, &serde_json::json!({"long":"payload"}), &mut 1).unwrap_err();
    assert!(error.downcast_ref::<ResourceAdmissionError>().is_some());
    assert!(!out.exists());
    let mut remaining = 1024;
    write_metadata(&out, &serde_json::json!({"ok":true}), &mut remaining)?;
    assert_eq!(1024 - remaining, std::fs::metadata(out)?.len());
    let _ = small.request(root.path(), root.path())?;
    Ok(())
}

#[test]
fn engine_wait_releases_scan_claim_and_cancellation_leaves_no_payloads() -> Result<()> {
    use rnaa_core::resources::{ExecutionResources, ReservationOwner, ResourceConfig};
    use std::{
        sync::{
            atomic::{AtomicBool, Ordering},
            mpsc,
        },
        time::Duration,
    };
    for cancel in [false, true] {
        let root = tempfile::tempdir()?;
        let paths = rnaa_core::ProjectPaths::new(root.path());
        paths.ensure_layout()?;
        let source = root.path().join("reads.fastq");
        std::fs::write(
            &source,
            format!("@read\n{}\n+\n{}\n", "ACGT".repeat(25), "I".repeat(100)),
        )?;
        let run = rnaa_core::RunRecord {
            project_id: "test".into(),
            study_accession: "study".into(),
            source_study_accession: None,
            run_accession: "RUN1".into(),
            sample_accession: None,
            experiment_accession: None,
            library_layout: LibraryLayout::Single,
            instrument_platform: None,
            instrument_model: None,
            remote_files: Vec::new(),
            metadata: Default::default(),
            state: rnaa_core::RunState::Verified,
            last_error: None,
            updated_at: rnaa_core::util::now_rfc3339(),
        };
        let inputs = vec![VerifiedFile {
            artifact_kind: ArtifactKind::FastqSingle,
            bytes: std::fs::metadata(&source)?.len(),
            path: source,
            checksum_type: "none".into(),
            checksum: String::new(),
            integrity_source: "fixture".into(),
        }];
        let mut config = rnaa_core::ProjectConfig::default();
        config.quant.preprocess_threads = 2;
        config.resources = ResourceConfig {
            pool_root: root.path().join("pool").display().to_string(),
            cpu_threads: 2,
            memory_mb: 128,
            disk_headroom_mb: 0,
            wait_seconds: 10,
        };
        let owner = ReservationOwner {
            project_root: root.path().into(),
            attempt_id: "current".into(),
            stage: "preprocess".into(),
            run_accession: Some("RUN1".into()),
        };
        let context = ExecutionResources::open(&config.resources, owner)?;
        let blocker = context.reserve(
            "blocker",
            ResourceRequest {
                memory_bytes: 80 * MIB,
                ..Default::default()
            },
            &|| Ok(()),
        )?;
        let cancelled = AtomicBool::new(false);
        let saw_scan = AtomicBool::new(false);
        let sent = AtomicBool::new(false);
        let (tx, rx) = mpsc::channel();
        std::thread::scope(|scope| -> Result<()> {
            let worker = scope.spawn(|| {
                super::super::FasterpInProcessPreprocessor.preprocess_with_resources(
                    &run,
                    &inputs,
                    &paths,
                    &config,
                    Some(&context),
                    &|| {
                        if cancelled.load(Ordering::SeqCst) {
                            anyhow::bail!("cancelled while waiting");
                        }
                        let snapshot = context.pool().snapshot()?;
                        if snapshot
                            .reservations
                            .iter()
                            .any(|r| r.owner.stage == "preprocess.scan")
                        {
                            saw_scan.store(true, Ordering::SeqCst);
                        } else if saw_scan.load(Ordering::SeqCst)
                            && !sent.swap(true, Ordering::SeqCst)
                        {
                            tx.send(())?;
                        }
                        Ok(())
                    },
                )
            });
            let observed = rx.recv_timeout(Duration::from_secs(5));
            // Always unblock before joining, including a failed observation.
            if observed.is_err() {
                cancelled.store(true, Ordering::SeqCst);
                drop(blocker);
                let _ = worker.join();
                observed?;
                unreachable!();
            }
            assert_eq!(context.pool().snapshot()?.reservations.len(), 1);
            assert!(!paths.preprocess_run_dir("RUN1", "kallistors").exists());
            if cancel {
                cancelled.store(true, Ordering::SeqCst);
            } else {
                drop(blocker);
            }
            let result = worker
                .join()
                .map_err(|_| anyhow::anyhow!("preprocessor panicked"))?;
            if cancel {
                let error = result
                    .err()
                    .context("cancellation unexpectedly succeeded")?;
                assert!(error.downcast_ref::<ResourceAdmissionError>().is_some());
                assert!(!paths.preprocess_run_dir("RUN1", "kallistors").exists());
            } else {
                result?;
            }
            Ok(())
        })?;
        assert!(context.pool().snapshot()?.reservations.is_empty());
    }
    Ok(())
}
