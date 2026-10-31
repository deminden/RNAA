use super::*;
use rnaa_core::{
    ArtifactKind, ArtifactRecord, LibraryLayout, ManifestArtifact, ProjectPaths, RemoteFile,
    RemoteFileKind, RunState,
    generations::{Generation, OutputGroup},
    jobs::{JobLease, JobOutcome},
    util::now_rfc3339,
};
use std::fs;

fn setup() -> Result<(tempfile::TempDir, Database, ProjectConfig, RunRecord)> {
    let temp = tempfile::tempdir()?;
    let config = ProjectConfig::default();
    config.save(&temp.path().join("rnaa.toml"))?;
    let db = Database::create(temp.path(), &config)?;
    let run = RunRecord {
        project_id: db.project_id().into(),
        study_accession: "study".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: vec![],
        metadata: Default::default(),
        state: RunState::Verified,
        last_error: None,
        updated_at: now_rfc3339(),
    };
    db.upsert_run(&run)?;
    Ok((temp, db, config, run))
}

#[test]
fn setting_observations_ignore_execution_but_detect_source_layout_and_active_policy_changes()
-> Result<()> {
    let (temp, db, expected, run) = setup()?;
    let mut current = expected.clone();
    current.quant.threads += 1;
    current.quant.worker_budget = 1;
    current.quant.preprocess_threads += 1;
    current.quant.qc_gate.min_pass_rate = Some(0.9); // Inactive preprocessing policy.
    db.set_project_config(&current)?;
    current.save(&temp.path().join("rnaa.toml"))?;
    verify_quant_current(&db, &run, &expected)?;
    current.quant.single_fragment_length += 10.0;
    db.set_project_config(&current)?;
    assert!(verify_quant_catalog(&db, &run, &expected).is_err());
    db.set_project_config(&expected)?;
    current = expected.clone();
    current.quant.preprocess = true;
    db.set_project_config(&current)?;
    assert!(verify_quant_catalog(&db, &run, &expected).is_err());
    db.set_project_config(&expected)?;
    let mut changed = run.clone();
    changed.library_layout = LibraryLayout::Paired;
    db.upsert_run(&changed)?;
    assert!(verify_quant_catalog(&db, &run, &expected).is_err());
    current = expected.clone();
    current.quant.single_fragment_length = 900.0;
    db.set_project_config(&current)?;
    verify_quant_catalog(&db, &changed, &expected)?;
    db.set_project_config(&expected)?;
    changed = run.clone();
    changed.remote_files = vec![RemoteFile {
        url: "http://example.invalid/new.fastq".into(),
        md5: None,
        bytes: None,
        kind: RemoteFileKind::Fastq,
        basename: "RUN1.fastq".into(),
    }];
    db.upsert_run(&changed)?;
    assert!(verify_quant_catalog(&db, &run, &expected).is_err());
    Ok(())
}

#[test]
fn changed_settings_cannot_commit_a_sealed_old_count_request() -> Result<()> {
    let (temp, db, config, run) = setup()?;
    let prepare = |value: &str| -> Result<(JobLease, Generation, Vec<ArtifactRecord>)> {
        let lease = JobLease::try_acquire(
            &db,
            "quant",
            Some("RUN1"),
            &[RunState::Verified, RunState::QuantDone],
        )?
        .unwrap();
        let owned = lease.database();
        owned.set_job_request(&"a".repeat(64))?;
        let generation = Generation::begin(&owned, OutputGroup::Quant, "quant")?;
        let paths = generation.paths(&ProjectPaths::new(temp.path())).unwrap();
        let path = paths
            .quant_run_dir("RUN1", "kallistors")
            .join("abundance.tsv");
        fs::create_dir_all(path.parent().unwrap())?;
        fs::write(&path, value)?;
        let input = ManifestArtifact::snapshot("ABUNDANCE_TSV", &path)?;
        let record = ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some("RUN1".into()),
            kind: ArtifactKind::QuantAbundanceTsv,
            path: input.path,
            blob_id: None,
            shared_path: None,
            checksum_type: input.checksum_type,
            checksum: input.checksum,
            bytes: input.bytes,
            created_at: now_rfc3339(),
        };
        Ok((lease, generation, vec![record]))
    };
    let (mut first_lease, mut first, records) = prepare("old complete counts")?;
    let paths = records
        .iter()
        .map(|r| std::path::PathBuf::from(&r.path))
        .collect::<Vec<_>>();
    first.seal(&paths, &"a".repeat(64))?;
    first.commit(&records, |db| {
        db.set_run_state("RUN1", RunState::QuantDone, None)
    })?;
    let previous = first.head().clone();
    first_lease.finish(JobOutcome::Succeeded, None)?;
    drop(first_lease);
    let (mut lease, mut replacement, records) = prepare("new incomplete setting observation")?;
    replacement.seal(
        &records
            .iter()
            .map(|r| std::path::PathBuf::from(&r.path))
            .collect::<Vec<_>>(),
        &"a".repeat(64),
    )?;
    let mut changed = config.clone();
    changed.quant.single_fragment_length += 10.0;
    db.set_project_config(&changed)?;
    assert!(
        replacement
            .commit(&records, |db| verify_quant_catalog(db, &run, &config))
            .is_err()
    );
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")?.unwrap(),
        previous
    );
    assert!(db.verify_generation_head(&previous)?);
    assert_eq!(
        fs::read_to_string(paths.first().unwrap())?,
        "old complete counts"
    );
    drop(replacement);
    lease.finish(JobOutcome::Failed, None)?;
    Ok(())
}
