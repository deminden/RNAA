use super::*;
use crate::{LibraryLayout, ProjectConfig, ProjectStageState};

fn run(db: &Database, id: &str, sample: Option<&str>) -> RunRecord {
    RunRecord {
        project_id: db.project_id().into(),
        study_accession: "STUDY".into(),
        source_study_accession: None,
        run_accession: id.into(),
        sample_accession: sample.map(Into::into),
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: Vec::new(),
        metadata: BTreeMap::new(),
        state: RunState::QuantDone,
        last_error: None,
        updated_at: now_rfc3339(),
    }
}

#[test]
fn resolving_narrower_and_empty_cohorts_preserves_history_and_fences_old_analysis() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let db = Database::create(dir.path(), &ProjectConfig::default())?;
    let runs = vec![
        run(&db, "RUN1", Some("BIO1")),
        run(&db, "RUN2", Some("BIO2")),
    ];
    db.upsert_runs(&runs)?;
    let legacy = CohortSnapshot::capture(&db)?;
    let selectors = db.selector_contract()?;
    let first = db.publish_resolved_cohort(&runs, &selectors, None)?;
    assert!(legacy.verify_catalog(&db).is_err());
    db.set_project_stage_state("normalize", ProjectStageState::Done, None)?;
    let snapshot = CohortSnapshot::capture(&db)?;
    let reopened = Database::open(dir.path())?;
    assert_eq!(reopened.list_active_runs()?.len(), 2);
    let second =
        reopened.publish_resolved_cohort(&runs[..1], &selectors, Some(&first.cohort_id))?;
    assert_eq!(db.list_active_runs()?.len(), 1);
    assert_eq!(db.list_runs()?.len(), 2);
    assert_eq!(db.get_run("RUN2")?.unwrap().state, RunState::QuantDone);
    assert_eq!(
        db.get_project_stage_state("normalize")?.unwrap().state,
        ProjectStageState::Pending
    );
    assert!(snapshot.verify_catalog(&db).is_err());
    assert!(
        db.publish_resolved_cohort(&runs, &selectors, Some(&first.cohort_id))
            .is_err()
    );
    db.publish_resolved_cohort(&[], &selectors, Some(&second.cohort_id))?;
    assert!(db.list_active_runs()?.is_empty());
    assert_eq!(db.list_runs()?.len(), 2);
    Ok(())
}

#[test]
fn pending_selectors_foreign_members_and_corrupt_descriptors_fail_closed() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let db = Database::create(dir.path(), &ProjectConfig::default())?;
    let runs = vec![run(&db, "RUN1", None)];
    let contract = db.selector_contract()?;
    let cohort = db.publish_resolved_cohort(&runs, &contract, None)?;
    db.add_input("condition=treated", InputType::FilterInclude)?;
    assert!(db.list_active_runs().is_err());
    assert!(
        db.publish_resolved_cohort(&runs, &contract, Some(&cohort.cohort_id))
            .is_err()
    );
    let current = db.selector_contract()?;
    let mut foreign = runs[0].clone();
    foreign.project_id = "other".into();
    assert!(
        db.publish_resolved_cohort(&[foreign], &current, Some(&cohort.cohort_id))
            .is_err()
    );
    assert!(
        db.publish_resolved_cohort(
            &[runs[0].clone(), runs[0].clone()],
            &current,
            Some(&cohort.cohort_id)
        )
        .is_err()
    );
    assert_eq!(db.get_run("RUN1")?.unwrap().project_id, db.project_id());
    db.publish_resolved_cohort(&[], &current, Some(&cohort.cohort_id))?;
    db.write(|conn| {
        conn.execute("UPDATE cohorts SET descriptor_json = '{}'", [])?;
        Ok(())
    })?;
    assert!(db.list_active_runs().is_err());
    Ok(())
}

#[test]
fn frozen_metadata_ignores_progress_but_detects_design_and_identity_changes() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let db = Database::create(dir.path(), &ProjectConfig::default())?;
    let mut value = run(&db, "RUN1", Some("BIO1"));
    value.metadata.insert("condition".into(), "control".into());
    db.upsert_run(&value)?;
    let original = CohortSnapshot::capture(&db)?;
    db.set_run_state("RUN1", RunState::QuantDone, None)?;
    original.verify_catalog(&db)?;
    value.metadata.insert("condition".into(), "treated".into());
    db.upsert_run(&value)?;
    assert!(original.verify_catalog(&db).is_err());
    assert_eq!(original.digest()?, CohortSnapshot::capture(&db)?.digest()?);
    let metadata_observed = CohortSnapshot::capture(&db)?;
    value.remote_files.push(crate::RemoteFile {
        url: "http://example.invalid/new-source.fastq".into(),
        md5: None,
        bytes: None,
        kind: crate::RemoteFileKind::Fastq,
        basename: "RUN1.fastq".into(),
    });
    db.upsert_run(&value)?;
    assert!(metadata_observed.verify_catalog(&db).is_err());
    assert_eq!(
        metadata_observed.digest()?,
        CohortSnapshot::capture(&db)?.digest()?
    );
    let cohort = db.publish_resolved_cohort(&[value.clone()], &db.selector_contract()?, None)?;
    value.sample_accession = Some("BIO2".into());
    db.upsert_run(&value)?;
    assert!(db.list_active_runs().is_err());
    assert_eq!(db.active_cohort()?.unwrap().cohort_id, cohort.cohort_id);
    Ok(())
}

#[test]
fn scientific_membership_ignores_equivalent_selectors_and_annotations_but_tracks_identity()
-> Result<()> {
    let temp = tempfile::tempdir()?;
    let db = Database::create(temp.path(), &ProjectConfig::default())?;
    let mut value = run(&db, "RUN1", Some("BIO1"));
    let first = db.publish_resolved_cohort(&[value.clone()], &db.selector_contract()?, None)?;
    let original = CohortSnapshot::capture(&db)?;
    db.add_input("RUN1", InputType::Run)?;
    value
        .metadata
        .insert("notes".into(), "updated annotation".into());
    let second = db.publish_resolved_cohort(
        &[value.clone()],
        &db.selector_contract()?,
        Some(&first.cohort_id),
    )?;
    assert_eq!(original.digest()?, CohortSnapshot::capture(&db)?.digest()?);
    assert!(original.verify_catalog(&db).is_err());
    value.sample_accession = Some("BIO2".into());
    db.publish_resolved_cohort(&[value], &db.selector_contract()?, Some(&second.cohort_id))?;
    assert_ne!(original.digest()?, CohortSnapshot::capture(&db)?.digest()?);
    Ok(())
}

#[test]
fn biological_identity_preserves_declared_accession_and_fallback_sources() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let db = Database::create(dir.path(), &ProjectConfig::default())?;
    let a = run(&db, "RUN1", Some("sam1"));
    let mut b = run(&db, "RUN2", Some(" SAM1 "));
    assert_eq!(
        BiologicalSample::for_run(&a)?,
        BiologicalSample::for_run(&b)?
    );
    assert_eq!(BiologicalSample::for_run(&a)?.id, "SAM1");
    b.metadata
        .insert("biological_sample_id".into(), "independent".into());
    assert_eq!(BiologicalSample::for_run(&b)?.id, "independent");
    assert_eq!(
        BiologicalSample::for_run(&b)?.source,
        IdentitySource::Declared
    );
    b.metadata
        .insert("biological_sample_id".into(), String::new());
    assert!(BiologicalSample::for_run(&b).is_err());
    assert_eq!(
        BiologicalSample::for_run(&run(&db, "RUN3", None))?.source,
        IdentitySource::RunFallback
    );
    Ok(())
}

#[test]
fn an_accession_can_be_a_query_and_an_exclusion_without_silently_losing_intent() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let db = Database::create(dir.path(), &ProjectConfig::default())?;
    db.add_input("SAM1", InputType::Sample)?;
    db.add_input("SAM1", InputType::SampleInclude)?;
    db.add_input("SAM1", InputType::SampleExclude)?;
    db.add_input("SAM1", InputType::SampleExclude)?;
    assert_eq!(Database::open(dir.path())?.list_inputs()?.len(), 3);
    Ok(())
}

#[test]
fn cohort_publication_rolls_back_if_member_upsert_fails() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let db = Database::create(dir.path(), &ProjectConfig::default())?;
    db.write(|conn| { conn.execute_batch("CREATE TRIGGER fail_second BEFORE INSERT ON runs WHEN NEW.run_accession = 'RUN2' BEGIN SELECT RAISE(ABORT, 'injected'); END;")?; Ok(()) })?;
    let contract = db.selector_contract()?;
    let runs = vec![run(&db, "RUN1", None), run(&db, "RUN2", None)];
    assert!(db.publish_resolved_cohort(&runs, &contract, None).is_err());
    assert!(db.list_runs()?.is_empty());
    assert!(db.active_cohort()?.is_none());
    Ok(())
}

#[test]
fn equivalent_resolutions_reuse_semantics_but_fence_superseded_producers() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let mut config = ProjectConfig::default();
    let db = Database::create(dir.path(), &config)?;
    let runs = vec![run(&db, "RUN1", None)];
    let contract = db.selector_contract()?;
    let first = db.publish_resolved_cohort(&runs, &contract, None)?;
    let snapshot = CohortSnapshot::capture(&db)?;
    db.publish_resolved_cohort(&runs, &contract, Some(&first.cohort_id))?;
    assert_eq!(snapshot.digest()?, CohortSnapshot::capture(&db)?.digest()?);
    assert!(snapshot.verify_catalog(&db).is_err());
    config.resolver.merge_strategy = MetadataMergeStrategy::PreferRemote;
    assert!(db.verify_cohort_config(&config).is_err());
    assert!(snapshot.verify_config(&config).is_err());
    db.set_project_config(&config)?;
    assert!(db.list_active_runs().is_err());
    Ok(())
}

#[test]
fn selector_change_during_a_live_fit_rolls_back_the_analysis_publication() -> Result<()> {
    use crate::{
        ArtifactKind, ArtifactRecord, ManifestArtifact,
        generations::{Generation, OutputGroup},
        jobs::{JobLease, JobOutcome},
    };
    let dir = tempfile::tempdir()?;
    let db = Database::create(dir.path(), &ProjectConfig::default())?;
    db.upsert_run(&run(&db, "RUN1", None))?;
    let mut lease = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    let owned = lease.database();
    let snapshot = CohortSnapshot::capture(&owned)?;
    owned.set_job_request(&snapshot.digest()?)?;
    let mut generation = Generation::begin(&owned, OutputGroup::De, "normalize")?;
    let file = generation.head().payload_path.join("counts.tsv");
    std::fs::write(&file, "old cohort counts")?;
    let artifact = ManifestArtifact::snapshot("COUNTS", &file)?;
    let record = ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: None,
        kind: ArtifactKind::Counts,
        path: artifact.path,
        blob_id: None,
        shared_path: None,
        checksum_type: artifact.checksum_type,
        checksum: artifact.checksum,
        bytes: artifact.bytes,
        created_at: now_rfc3339(),
    };
    generation.seal(&[file], &snapshot.digest()?)?;
    // Membership intent may change while a fit holds the analysis resource.
    Database::open(dir.path())?.add_input("new selection", InputType::FilterInclude)?;
    assert!(
        generation
            .commit(&[record], |db| snapshot.verify_catalog(db))
            .is_err()
    );
    assert!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .is_none()
    );
    assert!(db.list_artifacts()?.is_empty());
    lease.finish(JobOutcome::Failed, Some("selection changed"))?;
    Ok(())
}
