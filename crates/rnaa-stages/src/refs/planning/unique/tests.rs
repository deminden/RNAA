use super::*;
use rnaa_core::resources::{ReservationOwner, ResourceConfig};
fn context(root: &std::path::Path, strict: bool) -> Result<ExecutionResources> {
    ExecutionResources::open(
        &ResourceConfig {
            pool_root: root.join("pool").display().to_string(),
            cpu_threads: 2,
            memory_mb: 256,
            disk_headroom_mb: if strict { u64::MAX / (1024 * 1024) } else { 0 },
            wait_seconds: 5,
        },
        ReservationOwner {
            project_root: root.into(),
            attempt_id: "planner-test".into(),
            stage: "refs".into(),
            run_accession: None,
        },
    )
}
#[test]
fn external_merges_count_exactly_across_online_levels_and_release_all_scratch() -> Result<()> {
    let root = tempfile::tempdir()?;
    let context = context(root.path(), false)?;
    let mut words = UniqueWords::new(&context, 80 * 1024 * 1024, root.path(), &|| Ok(()))?;
    let path = words.directory.path().to_path_buf();
    let mut oracle = std::collections::BTreeSet::new();
    for group in 0..FAN_IN * FAN_IN + 7 {
        for position in 0..100 {
            let value = ((position * 37 + group % 11) % 127) as u64;
            oracle.insert(value);
            words.push(value)?;
        }
        words.flush()?;
        assert!(words.levels.iter().all(|level| level.len() < FAN_IN));
    }
    assert_eq!(words.finish()?, oracle.len() as u64);
    assert!(!path.exists());
    assert!(context.pool().snapshot()?.reservations.is_empty());
    Ok(())
}
#[test]
fn scratch_growth_refusal_releases_written_files_and_does_not_wait_with_cpu_held() -> Result<()> {
    let root = tempfile::tempdir()?;
    let context = context(root.path(), false)?;
    let mut words = UniqueWords::new(&context, 80 * 1024 * 1024, root.path(), &|| Ok(()))?;
    let path = words.directory.path().to_path_buf();
    for group in 0..4 {
        for i in 0..WORDS {
            words.push((group * WORDS + i) as u64)?;
        }
    }
    let strict = self::context(root.path(), true)?;
    let blocker = strict.reserve(
        "blocker",
        ResourceRequest {
            cpu_threads: 1,
            ..Default::default()
        },
        &|| Ok(()),
    )?;
    for i in 0..100 {
        words.push((4 * WORDS + i) as u64)?;
    }
    let error = words.flush().unwrap_err();
    assert!(error.is::<ResourceAdmissionError>(), "{error:#}");
    drop(words);
    assert!(!path.exists());
    assert_eq!(context.pool().snapshot()?.reservations.len(), 1);
    drop(blocker);
    assert!(context.pool().snapshot()?.reservations.is_empty());
    Ok(())
}
