use super::*;
use crate::fastq::{FastqLimits, FastqPlan};
use crate::quant::resources::QuantResources;
use rnaa_core::resources::{ExecutionResources, ReservationOwner, ResourceConfig, ResourceRequest};
use std::{
    sync::{atomic::AtomicBool, mpsc},
    time::Duration,
};
const MIB: u64 = 1024 * 1024;
fn context(root: &Path, label: &str, memory: u64) -> Result<ExecutionResources> {
    ExecutionResources::open(
        &ResourceConfig {
            pool_root: root.join("pool").display().to_string(),
            cpu_threads: 2,
            memory_mb: memory,
            disk_headroom_mb: 0,
            wait_seconds: 10,
        },
        ReservationOwner {
            project_root: root.join(label),
            attempt_id: label.into(),
            stage: "quant".into(),
            run_accession: Some("RUN".into()),
        },
    )
}
fn estimate(
    root: &Path,
    reference: &ReferenceBundle,
    group: &QuantIndexGroup,
    sequence: &str,
) -> Result<QuantResources> {
    let files = inputs(root, sequence, false)?;
    let plan = FastqPlan::inspect(
        LibraryLayout::Single,
        &files,
        FastqLimits::default(),
        &|| Ok(()),
    )?;
    QuantResources::estimate(&group.preflight(reference, &|| Ok(()))?, &plan, 1, MIB)
}

#[test]
fn shared_graph_is_charged_once_and_survives_only_its_real_owners() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, sequence) = reference(root.path())?;
    let group = group(&reference)?;
    let estimate = estimate(root.path(), &reference, &group, &sequence)?;
    let context = context(root.path(), "owner", 512)?;
    let request = estimate.run_request(root.path(), root.path());
    let first = group.acquire(
        &reference,
        &estimate,
        request.clone(),
        Some(&context),
        &|| Ok(()),
    )?;
    let weak = Arc::downgrade(&first.index);
    let second = group.acquire(&reference, &estimate, request, Some(&context), &|| Ok(()))?;
    assert_eq!(context.pool().snapshot()?.reservations.len(), 3);
    assert_eq!(
        context.pool().snapshot()?.memory_bytes,
        estimate.index_resident_bytes + 2 * estimate.run_memory_bytes
    );
    assert!(!group.release_idle()?);
    assert_eq!(group.stats().successful_loads, 1);
    drop(first);
    assert_eq!(context.pool().snapshot()?.reservations.len(), 2);
    assert!(!group.release_idle()?);
    drop(second);
    assert_eq!(context.pool().snapshot()?.reservations.len(), 1);
    assert!(weak.upgrade().is_some());
    assert!(group.release_idle()?);
    assert!(weak.upgrade().is_none());
    assert!(context.pool().snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn idle_graph_is_evicted_before_waiting_for_run_capacity() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, sequence) = reference(root.path())?;
    let a = group(&reference)?;
    let b = group(&reference)?;
    let estimate = estimate(root.path(), &reference, &a, &sequence)?;
    let context_a = context(root.path(), "a", 256)?;
    let context_b = context(root.path(), "b", 256)?;
    let request = estimate.run_request(root.path(), root.path());
    let first = a.acquire(
        &reference,
        &estimate,
        request.clone(),
        Some(&context_a),
        &|| Ok(()),
    )?;
    let weak = Arc::downgrade(&first.index);
    drop(first);
    drop(b.acquire(
        &reference,
        &estimate,
        request.clone(),
        Some(&context_b),
        &|| Ok(()),
    )?);
    let blocker = context_b.reserve(
        "blocker",
        ResourceRequest {
            memory_bytes: 125 * MIB,
            ..Default::default()
        },
        &|| Ok(()),
    )?;
    let (tx, rx) = mpsc::channel();
    let notified = AtomicBool::new(false);
    let cancelled = AtomicBool::new(false);
    std::thread::scope(|scope| -> Result<()> {
        let worker = scope.spawn(|| {
            a.acquire(&reference, &estimate, request, Some(&context_a), &|| {
                if cancelled.load(Ordering::SeqCst) {
                    bail!("cancelled wait");
                }
                if a.stats().evictions > 0 && !notified.swap(true, Ordering::SeqCst) {
                    tx.send(())?;
                }
                Ok(())
            })
        });
        let observed = rx.recv_timeout(Duration::from_secs(5));
        if observed.is_err() {
            cancelled.store(true, Ordering::SeqCst);
            drop(blocker);
            let _ = worker.join();
            observed?;
            unreachable!();
        }
        assert!(weak.upgrade().is_none());
        assert_eq!(context_a.pool().snapshot()?.reservations.len(), 2);
        drop(blocker);
        let run = worker.join().unwrap()?;
        assert_eq!(a.stats().successful_loads, 2);
        drop(run);
        Ok(())
    })?;
    a.release_idle()?;
    b.release_idle()?;
    assert!(context_a.pool().snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn preflight_wait_can_release_the_same_groups_idle_index() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, sequence) = reference(root.path())?;
    let group = group(&reference)?;
    let estimate = estimate(root.path(), &reference, &group, &sequence)?;
    let context = context(root.path(), "owner", 256)?;
    let run = group.acquire(
        &reference,
        &estimate,
        estimate.run_request(root.path(), root.path()),
        Some(&context),
        &|| Ok(()),
    )?;
    let weak = Arc::downgrade(&run.index);
    drop(run);
    let blocker = context.reserve(
        "blocker",
        ResourceRequest {
            memory_bytes: 200 * MIB,
            ..Default::default()
        },
        &|| Ok(()),
    )?;
    let scan = group.acquire_scan(
        Some(&context),
        crate::quant::resources::scan_request(FastqLimits::default())?,
        &|| Ok(()),
    )?;
    assert!(weak.upgrade().is_none());
    assert_eq!(context.pool().snapshot()?.reservations.len(), 2);
    drop(scan);
    drop(blocker);
    assert!(context.pool().snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn cancellation_before_or_after_load_releases_both_bundled_claims() -> Result<()> {
    for after_load in [false, true] {
        let root = tempfile::tempdir()?;
        let (reference, sequence) = reference(root.path())?;
        let group = group(&reference)?;
        let estimate = estimate(root.path(), &reference, &group, &sequence)?;
        let context = context(root.path(), "owner", 256)?;
        let seen = AtomicBool::new(false);
        let result = group.acquire(
            &reference,
            &estimate,
            estimate.run_request(root.path(), root.path()),
            Some(&context),
            &|| {
                if context.pool().snapshot()?.reservations.len() == 2
                    && (!after_load || seen.swap(true, Ordering::SeqCst))
                {
                    bail!("cancelled bundled load");
                }
                Ok(())
            },
        );
        assert!(result.is_err());
        assert!(context.pool().snapshot()?.reservations.is_empty());
        assert_eq!(group.stats().load_attempts, usize::from(after_load));
        assert_eq!(group.stats().successful_loads, 0);
        assert!(!group.release_idle()?);
    }
    Ok(())
}

#[test]
fn admitted_paired_mapping_matches_the_present_engine_and_releases_run_memory() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, sequence) = reference(root.path())?;
    let files = inputs(root.path(), &sequence, true)?;
    let run = RunRecord {
        project_id: "test".into(),
        study_accession: "test".into(),
        source_study_accession: None,
        run_accession: "PAIRED".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Paired,
        instrument_platform: None,
        instrument_model: None,
        remote_files: Vec::new(),
        metadata: Default::default(),
        state: RunState::Verified,
        last_error: None,
        updated_at: "fixture".into(),
    };
    let mut config = ProjectConfig::default();
    config.quant.threads = 2;
    config.quant.match_cache_mb = 1;
    let expected = RustKallistoQuantifier.quantify(
        &run,
        &files,
        &reference,
        &ProjectPaths::new(root.path().join("baseline")),
        &config,
    )?;
    let context = context(root.path(), "paired", 512)?;
    let group = group(&reference)?;
    let actual = group.quantify_with_resources(
        &run,
        &files,
        &reference,
        &ProjectPaths::new(root.path().join("admitted")),
        &config,
        Some(&context),
        &|| Ok(()),
    )?;
    assert_eq!(
        fs::read(expected.out_dir.join("abundance.tsv"))?,
        fs::read(actual.out_dir.join("abundance.tsv"))?
    );
    let manifest: serde_json::Value = serde_json::from_slice(&fs::read(actual.manifest_path)?)?;
    assert_eq!(
        manifest["parameters"]["resource_execution"]["reservations"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(context.pool().snapshot()?.reservations.len(), 1);
    assert_eq!(context.pool().snapshot()?.cpu_threads, 0);
    group.release_idle()?;
    assert!(context.pool().snapshot()?.reservations.is_empty());
    Ok(())
}
