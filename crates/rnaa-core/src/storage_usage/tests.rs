use super::*;
use crate::{
    ArtifactKind, ArtifactRecord, ProjectConfig,
    generations::{Generation, OutputGroup},
    jobs::{JobLease, JobOutcome},
};

fn setup() -> Result<(tempfile::TempDir, Database)> {
    let root = tempfile::tempdir()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    Ok((root, db))
}
fn path<'a>(report: &'a StorageUsageReport, path: &Path) -> &'a PathUsage {
    let canonical = fs::canonicalize(path).unwrap();
    report
        .paths
        .as_ref()
        .unwrap()
        .iter()
        .find(|entry| entry.path == canonical)
        .unwrap()
}

#[test]
fn accounting_distinguishes_committed_unpublished_and_unindexed_payloads_without_cleanup()
-> Result<()> {
    let (root, db) = setup()?;
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::De, "deseq2")?;
    let file = generation.head().payload_path.join("counts.tsv");
    fs::write(&file, "counts")?;
    let snapshot = crate::ManifestArtifact::snapshot("COUNTS", &file)?;
    generation.seal(std::slice::from_ref(&file), &"a".repeat(64))?;
    generation.commit(
        &[ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: None,
            kind: ArtifactKind::Counts,
            path: snapshot.path,
            blob_id: None,
            shared_path: None,
            checksum_type: snapshot.checksum_type,
            checksum: snapshot.checksum,
            bytes: snapshot.bytes,
            created_at: now_rfc3339(),
        }],
        |_| Ok(()),
    )?;
    lease.finish(JobOutcome::Succeeded, None)?;
    drop(generation);
    drop(owned);
    drop(lease);
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    let generation = Generation::begin(&owned, OutputGroup::De, "deseq2")?;
    let unfinished = generation.head().payload_path.join("partial.tsv");
    fs::write(&unfinished, "partial")?;
    let active = db.storage_usage(None, true)?;
    assert!(active.active_attempts_observed);
    assert_eq!(
        path(&active, &unfinished).category,
        "unpublished_generation"
    );
    drop(generation);
    lease.finish(JobOutcome::Failed, None)?;
    drop(owned);
    drop(lease);
    let unknown = root.path().join(".rnaa_generations/unregistered/data");
    fs::create_dir_all(unknown.parent().unwrap())?;
    fs::write(&unknown, "unknown")?;
    let report = db.storage_usage(None, true)?;
    assert!(report.complete, "{:?}", report.errors);
    assert!(!report.active_attempts_observed);
    assert_eq!(path(&report, &file).category, "committed_generation");
    assert_eq!(path(&report, &unfinished).category, "abandoned_generation");
    assert_eq!(path(&report, &unknown).category, "unindexed_generation");
    assert_eq!(
        report.logical_bytes,
        report
            .categories
            .values()
            .map(|value| value.logical_bytes)
            .sum::<u64>()
    );
    assert_eq!(fs::read(unfinished)?, b"partial");
    assert_eq!(fs::read(unknown)?, b"unknown");
    Ok(())
}

#[test]
#[cfg(unix)]
fn hardlinks_are_counted_once_across_project_and_shared_roots_without_following_symlinks()
-> Result<()> {
    let (root, db) = setup()?;
    let shared = tempfile::tempdir()?;
    let file = root.path().join("data");
    fs::write(&file, vec![b'x'; 8192])?;
    let alias = shared.path().join("alias");
    fs::hard_link(&file, &alias)?;
    let external = tempfile::tempdir()?;
    let secret = external.path().join("outside");
    fs::write(&secret, vec![b'y'; 16384])?;
    let link = root.path().join("outside-link");
    std::os::unix::fs::symlink(external.path(), &link)?;
    let report = db.storage_usage(Some(shared.path()), true)?;
    assert!(report.complete, "{:?}", report.errors);
    assert_eq!(report.observed_hardlink_aliases, Some(1));
    assert_eq!(report.unobserved_hardlinks, Some(0));
    assert_eq!(
        path(&report, &file).allocated_bytes,
        path(&report, &alias).allocated_bytes
    );
    assert_eq!(report.symlinks, 1);
    assert!(
        !report
            .paths
            .as_ref()
            .unwrap()
            .iter()
            .any(|entry| entry.path == fs::canonicalize(&secret).unwrap())
    );
    let mut expected = BTreeMap::new();
    for entry in report.paths.as_ref().unwrap() {
        let metadata = entry.path.symlink_metadata()?;
        expected
            .entry(file_key(&metadata, &entry.path))
            .or_insert(allocation(&metadata).unwrap());
    }
    assert_eq!(report.inode_allocated_bytes, Some(expected.values().sum()));
    assert!(!shared.path().join("catalog.sqlite3").exists());
    Ok(())
}

#[test]
fn nested_shared_store_is_scanned_once_and_cataloged_blobs_are_classified() -> Result<()> {
    let (root, db) = setup()?;
    let shared = root.path().join("shared");
    let store = crate::shared_store::SharedStore::open(&shared)?;
    let source = root.path().join("source");
    fs::write(&source, "blob")?;
    let sha = crate::util::compute_sha256(&source)?;
    let publication = store.prepare(&db, &source, &sha, 4)?;
    let blob = PathBuf::from(&publication.record().storage_path);
    drop(publication);
    let report = db.storage_usage(Some(&shared), true)?;
    assert_eq!(report.roots, vec![fs::canonicalize(root.path())?]);
    assert_eq!(path(&report, &blob).category, "shared_blob");
    assert_eq!(
        report
            .paths
            .as_ref()
            .unwrap()
            .iter()
            .filter(|entry| entry.path == blob)
            .count(),
        1
    );
    assert!(!report.catalog_changed);
    assert_eq!(report.shared_catalog_changed, Some(false));
    Ok(())
}
