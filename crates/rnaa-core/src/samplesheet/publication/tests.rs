use super::*;
use crate::{LibraryLayout, ProjectConfig, RunState};
use std::collections::BTreeMap;

fn run(db: &Database, condition: &str) -> RunRecord {
    RunRecord {
        project_id: db.project_id().into(),
        run_accession: "RUN1".into(),
        study_accession: "STUDY1".into(),
        sample_accession: Some("BIO1".into()),
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        source_study_accession: None,
        state: RunState::Resolved,
        last_error: None,
        metadata: BTreeMap::from([("condition".into(), condition.into())]),
        remote_files: Vec::new(),
        updated_at: now_rfc3339(),
    }
}
fn setup() -> Result<(tempfile::TempDir, ProjectPaths, Database, ColumnMap)> {
    let temp = tempfile::tempdir()?;
    let paths = ProjectPaths::new(temp.path());
    paths.ensure_layout()?;
    let db = Database::create(temp.path(), &ProjectConfig::default())?;
    let map =
        super::super::load_or_init_column_map(&paths.column_map_path(), &["condition".into()])?;
    Ok((temp, paths, db, map))
}
fn prepare(
    db: &Database,
    paths: &ProjectPaths,
    map: &ColumnMap,
    condition: &str,
    policy: ConflictPolicy,
) -> Result<PreparedPublication> {
    PreparedPublication::prepare(db, paths, &[run(db, condition)], &[], map, policy)
}
fn replace(path: &Path, old: &str, new: &str) -> Result<()> {
    let bytes = fs::read_to_string(path)?;
    assert!(bytes.contains(old));
    fs::write(path, bytes.replace(old, new))?;
    Ok(())
}

#[test]
fn edits_survive_refresh_conflicts_preserve_heads_and_explicit_policy_resolves() -> Result<()> {
    let (_temp, paths, db, map) = setup()?;
    let first = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    first.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    replace(&paths.samplesheet_path(), "control", "edited")?;
    let second = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    second.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("edited"));
    let id = db.metadata_publication()?.unwrap().publication_id;
    assert!(prepare(&db, &paths, &map, "source-update", ConflictPolicy::Fail).is_err());
    assert_eq!(db.metadata_publication()?.unwrap().publication_id, id);
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("edited"));
    let chosen = prepare(
        &db,
        &paths,
        &map,
        "source-update",
        ConflictPolicy::KeepLocal,
    )?;
    chosen.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("edited"));
    first.publication.verify()?;
    assert!(fs::read_to_string(&first.publication.merged.path)?.contains("control"));
    let source = prepare(&db, &paths, &map, "next-update", ConflictPolicy::UseSource)?;
    source.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("next-update"));
    Ok(())
}

#[test]
fn source_changes_stale_candidates_and_transaction_failures_cannot_publish() -> Result<()> {
    let (_temp, paths, db, map) = setup()?;
    let first = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    fs::write(
        paths.samplesheet_path(),
        "run_accession\tcondition\nRUN1\tedited\n",
    )?;
    assert!(first.commit_catalog(&db, &paths).is_err());
    assert!(db.metadata_publication()?.is_none());
    let second = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    let failed: Result<()> = db.atomic(|db| {
        second.commit_catalog(db, &paths)?;
        bail!("injected failure after metadata publication");
    });
    assert!(failed.is_err());
    assert!(db.metadata_publication()?.is_none());
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("edited"));
    Ok(())
}

#[test]
fn committed_views_recover_and_preserve_later_nonconflicting_edits() -> Result<()> {
    let (_temp, paths, db, map) = setup()?;
    let first = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    first.commit_catalog(&db, &paths)?;
    // Crash before the first editable view exists.
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    let mut record = run(&db, "updated");
    record.metadata.insert("donor".into(), "D1".into());
    let pending =
        PreparedPublication::prepare(&db, &paths, &[record], &[], &map, ConflictPolicy::Fail)?;
    pending.commit_catalog(&db, &paths)?;
    // A later editor adds an unrelated field before recovery.
    let text = fs::read_to_string(paths.samplesheet_path())?;
    let lines = text
        .lines()
        .enumerate()
        .map(|(i, line)| format!("{line}\t{}\n", if i == 0 { "notes" } else { "kept" }))
        .collect::<String>();
    fs::write(paths.samplesheet_path(), lines)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    let recovered = fs::read_to_string(paths.samplesheet_path())?;
    assert!(recovered.contains("updated"));
    assert!(recovered.contains("kept"));
    assert!(recovered.contains("D1"));
    // Crash after installing bytes but before recording installed_id.
    db.write(|conn| {
        conn.execute("UPDATE metadata_heads SET installed_id = NULL", [])?;
        Ok(())
    })?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    assert_eq!(fs::read_to_string(paths.samplesheet_path())?, recovered);
    Ok(())
}

#[test]
fn identity_edits_and_corrupt_witnesses_fail_without_overwriting_views() -> Result<()> {
    let (_temp, paths, db, map) = setup()?;
    let first = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    first.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    replace(&paths.samplesheet_path(), "BIO1", "BIO2")?;
    assert!(prepare(&db, &paths, &map, "control", ConflictPolicy::KeepLocal).is_err());
    replace(&paths.samplesheet_path(), "BIO2", "BIO1")?;
    fs::write(&first.publication.generated.path, b"corrupt")?;
    assert!(prepare(&db, &paths, &map, "control", ConflictPolicy::Fail).is_err());
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("control"));
    Ok(())
}

#[test]
fn narrowing_to_empty_and_reincluding_runs_preserves_their_deliberate_edits() -> Result<()> {
    let (_temp, paths, db, map) = setup()?;
    let first = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    first.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    replace(&paths.samplesheet_path(), "control", "edited")?;
    let empty = PreparedPublication::prepare(&db, &paths, &[], &[], &map, ConflictPolicy::Fail)?;
    empty.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    assert!(!fs::read_to_string(paths.samplesheet_path())?.contains("RUN1"));
    let restored = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    restored.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("edited"));
    fs::remove_file(paths.samplesheet_path())?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("edited"));
    Ok(())
}

#[test]
fn pending_conflicting_edits_are_not_overwritten_and_superseded_candidates_cannot_commit()
-> Result<()> {
    let (_temp, paths, db, map) = setup()?;
    let first = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    let superseded = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    first.commit_catalog(&db, &paths)?;
    assert!(superseded.commit_catalog(&db, &paths).is_err());
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    let pending = prepare(&db, &paths, &map, "source-update", ConflictPolicy::Fail)?;
    pending.commit_catalog(&db, &paths)?;
    replace(&paths.samplesheet_path(), "control", "later-edit")?;
    let view = fs::read(paths.samplesheet_path())?;
    assert!(recover_view(&db, &paths, ConflictPolicy::Fail).is_err());
    assert_eq!(fs::read(paths.samplesheet_path())?, view);
    recover_view(&db, &paths, ConflictPolicy::KeepLocal)?;
    assert_eq!(fs::read(paths.samplesheet_path())?, view);
    Ok(())
}

#[test]
fn override_witnesses_remain_valid_after_original_sources_change() -> Result<()> {
    let (_temp, paths, db, map) = setup()?;
    let original = paths.overrides_dir.join("source.tsv");
    fs::write(&original, "run_accession\tcondition\nRUN1\tcontrol\n")?;
    let observed = ManifestArtifact::snapshot("METADATA_OVERRIDE", &original)?;
    let first = PreparedPublication::prepare_with_sources(
        &db,
        &paths,
        &[run(&db, "control")],
        &[],
        &map,
        ConflictPolicy::Fail,
        &[observed],
    )?;
    first.commit_catalog(&db, &paths)?;
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    fs::write(&original, "run_accession\tcondition\nRUN1\tedited\n")?;
    first.publication.verify()?;
    assert!(first.verify_sources(&paths).is_err());
    let frozen = &first.publication.sources[1];
    assert!(fs::read_to_string(&frozen.path)?.contains("control"));
    fs::write(&frozen.path, b"corrupt")?;
    assert!(first.publication.verify().is_err());
    Ok(())
}

#[test]
fn metadata_head_event_and_success_roll_back_together_and_stale_owners_cannot_publish() -> Result<()>
{
    use crate::jobs::{JobLease, JobOutcome};
    let (_temp, paths, db, map) = setup()?;
    let candidate = prepare(&db, &paths, &map, "control", ConflictPolicy::Fail)?;
    let mut lease = JobLease::try_acquire(&db, "metadata", None, &[])?.unwrap();
    let owned = lease.database();
    let failed: Result<()> = owned.atomic(|db| {
        candidate.commit_catalog_preverified(db)?;
        db.complete_owned_job(JobOutcome::Succeeded, None)?;
        bail!("injected catalog failure after head/event/success");
    });
    assert!(failed.is_err());
    assert!(db.metadata_publication()?.is_none());
    assert_eq!(db.list_job_attempts(10)?[0].state, "RUNNING");
    let published_events = || {
        db.read(|conn| Ok(conn.query_row("SELECT COUNT(*) FROM events WHERE stage = 'metadata' AND message = 'published metadata revision'", [], |row| row.get::<_, i64>(0))?))
    };
    assert_eq!(published_events()?, 0);
    assert!(!paths.samplesheet_path().exists());
    owned.atomic(|db| {
        candidate.commit_catalog_preverified(db)?;
        db.complete_owned_job(JobOutcome::Succeeded, None)
    })?;
    assert_eq!(db.list_job_attempts(10)?[0].state, "SUCCEEDED");
    assert_eq!(published_events()?, 1);
    assert_eq!(
        db.metadata_publication()?.unwrap().publication_id,
        candidate.publication.publication_id
    );
    assert!(candidate.commit_catalog_preverified(&owned).is_err());
    recover_view(&db, &paths, ConflictPolicy::Fail)?;
    let later = prepare(&db, &paths, &map, "treated", ConflictPolicy::Fail)?;
    assert!(later.commit_catalog_preverified(&owned).is_err());
    assert_eq!(
        db.metadata_publication()?.unwrap().publication_id,
        candidate.publication.publication_id
    );
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(())
}
