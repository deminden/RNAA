use super::*;
use crate::{ArtifactKind, ProjectConfig, ProjectStageState, jobs::JobLease, util::compute_sha256};
use std::{
    process::{Child, Command, Stdio},
    thread,
    time::{Duration, Instant},
};
mod bindings;
mod lineage;
mod reuse;

fn setup() -> Result<(tempfile::TempDir, Database)> {
    let root = tempfile::tempdir()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    Ok((root, db))
}
fn preparing(db: &Database, value: &str) -> Result<(JobLease, Generation, Vec<ArtifactRecord>)> {
    let lease = JobLease::try_acquire(db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    owned.set_project_stage_state("deseq2", ProjectStageState::Running, None)?;
    let generation = Generation::begin(&owned, OutputGroup::De, "deseq2")?;
    let mut records = Vec::new();
    for name in ["counts.tsv", "normalized.tsv"] {
        let path = generation.head.payload_path.join(name);
        fs::write(&path, value)?;
        let snapshot = ManifestArtifact::snapshot("TEST", &path)?;
        records.push(ArtifactRecord {
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
        });
    }
    Ok((lease, generation, records))
}
fn ensure_run(db: &Database) -> Result<()> {
    use crate::{LibraryLayout, RunRecord, RunState};
    if db.get_run("RUN1")?.is_none() {
        db.upsert_run(&RunRecord {
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
        })?;
    }
    Ok(())
}

fn preparing_quant(
    db: &Database,
    value: &str,
) -> Result<(JobLease, Generation, Vec<ArtifactRecord>)> {
    use crate::RunState;
    ensure_run(db)?;
    let lease = JobLease::try_acquire(
        db,
        "quant",
        Some("RUN1"),
        &[
            RunState::Verified,
            RunState::PreprocessDone,
            RunState::QuantDone,
        ],
    )?
    .unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let generation = Generation::begin(&owned, OutputGroup::Quant, "quant")?;
    let paths = generation
        .paths(&ProjectPaths::new(db.path().parent().unwrap()))
        .unwrap();
    let output = paths.quant_run_dir("RUN1", "kallistors");
    fs::create_dir_all(&output)?;
    let mut records = Vec::new();
    for (kind, name) in [
        (ArtifactKind::QuantAbundanceH5, "abundance.h5"),
        (ArtifactKind::QuantAbundanceTsv, "abundance.tsv"),
        (ArtifactKind::QuantRunInfo, "run_info.json"),
        (ArtifactKind::Manifest, "rnaa_quant_manifest.json"),
    ] {
        let path = output.join(name);
        fs::write(&path, value)?;
        let snapshot = ManifestArtifact::snapshot(kind.as_str(), &path)?;
        records.push(ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some("RUN1".into()),
            kind,
            path: snapshot.path,
            blob_id: None,
            shared_path: None,
            checksum_type: snapshot.checksum_type,
            checksum: snapshot.checksum,
            bytes: snapshot.bytes,
            created_at: now_rfc3339(),
        });
    }
    Ok((lease, generation, records))
}

fn complete_quant(db: &Database) -> Result<()> {
    db.set_run_state("RUN1", crate::RunState::QuantDone, None)?;
    db.append_event(
        "quant",
        Some("RUN1"),
        "quant generation completed",
        serde_json::json!({}),
    )
}

fn publish_quant(db: &Database, value: &str) -> Result<GenerationHead> {
    let (mut lease, mut generation, records) = preparing_quant(db, value)?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    generation.commit(&records, complete_quant)?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(generation.head.clone())
}

#[test]
fn quant_publication_rejects_wrong_subject_and_rolls_back_shared_registration() -> Result<()> {
    let (_root, db) = setup()?;
    let previous = publish_quant(&db, "old")?;
    let (mut lease, mut generation, records) = preparing_quant(&db, "new")?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    let mut wrong_subject = records.clone();
    wrong_subject[0].run_accession = Some("RUN2".into());
    assert!(generation.commit(&wrong_subject, complete_quant).is_err());
    let blob = crate::SharedBlobRecord {
        blob_id: "prepared-blob".into(),
        storage_path: "shared".into(),
        checksum_type: "sha256".into(),
        checksum: records[0].checksum.clone(),
        bytes: records[0].bytes,
        created_at: now_rfc3339(),
    };
    let mut linked = records.clone();
    linked[0].blob_id = Some(blob.blob_id.clone());
    linked[0].shared_path = Some(blob.storage_path.clone());
    let result: Result<()> = generation.commit_prepared(
        &linked,
        |db| db.record_shared_blob(&blob),
        |db| {
            complete_quant(db)?;
            bail!("injected completion failure");
        },
    );
    assert!(result.is_err());
    assert!(db.get_shared_blob(&blob.blob_id)?.is_none());
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")?,
        Some(previous.clone())
    );
    let pinned = ProjectPaths::new(db.path().parent().unwrap()).with_generation_heads(&db)?;
    assert!(
        pinned
            .quant_run_dir("RUN1", "kallistors")
            .starts_with(&previous.payload_path)
    );
    assert!(
        !pinned
            .quant_run_dir("RUN2", "kallistors")
            .starts_with(&previous.payload_path)
    );
    assert!(
        !pinned
            .preprocess_run_dir("RUN1", "kallistors")
            .starts_with(&previous.payload_path)
    );
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1")?);
    assert_eq!(
        db.get_run("RUN1")?.unwrap().state,
        crate::RunState::QuantRunning
    );
    for record in db.list_artifacts_for_run(Some("RUN1"))? {
        assert_eq!(fs::read_to_string(record.path)?, "old");
    }
    lease.finish(JobOutcome::Failed, Some("test failure"))?;
    Ok(())
}

fn preparing_download(
    db: &Database,
    value: &str,
) -> Result<(JobLease, Generation, Vec<ArtifactRecord>)> {
    use crate::RunState;
    ensure_run(db)?;
    let lease = JobLease::try_acquire(
        db,
        "download",
        Some("RUN1"),
        &[
            RunState::Resolved,
            RunState::Verified,
            RunState::DownloadFailed,
        ],
    )?
    .unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let generation = Generation::begin(&owned, OutputGroup::Download, "download")?;
    let paths = generation
        .paths(&ProjectPaths::new(db.path().parent().unwrap()))
        .unwrap();
    fs::create_dir_all(paths.raw_run_dir("RUN1"))?;
    let mut records = Vec::new();
    for (kind, path) in [
        (
            ArtifactKind::FastqSingle,
            paths.raw_run_dir("RUN1").join("lane1.fastq"),
        ),
        (
            ArtifactKind::FastqSingle,
            paths.raw_run_dir("RUN1").join("lane2.fastq"),
        ),
        (
            ArtifactKind::Manifest,
            paths.manifest_path("download", "RUN1"),
        ),
    ] {
        fs::write(&path, value)?;
        let artifact = ManifestArtifact::snapshot(kind.as_str(), &path)?;
        records.push(ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some("RUN1".into()),
            kind,
            path: artifact.path,
            blob_id: None,
            shared_path: None,
            checksum_type: artifact.checksum_type,
            checksum: artifact.checksum,
            bytes: artifact.bytes,
            created_at: now_rfc3339(),
        });
    }
    Ok((lease, generation, records))
}
fn complete_download(db: &Database) -> Result<()> {
    db.set_run_state("RUN1", crate::RunState::Verified, None)?;
    db.append_event(
        "download",
        Some("RUN1"),
        "download generation completed",
        serde_json::json!({}),
    )
}
fn publish_download(db: &Database, value: &str) -> Result<(GenerationHead, Vec<ManifestArtifact>)> {
    let (mut lease, mut generation, records) = preparing_download(db, value)?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    generation.commit(&records, complete_download)?;
    lease.finish(JobOutcome::Succeeded, None)?;
    let snapshots = records[..2]
        .iter()
        .map(|record| ManifestArtifact::snapshot(record.kind.as_str(), Path::new(&record.path)))
        .collect::<Result<_>>()?;
    Ok((generation.head.clone(), snapshots))
}
#[test]
fn raw_generation_pruning_requires_a_committed_consumer_and_keeps_its_manifest() -> Result<()> {
    let (_root, db) = setup()?;
    let (head, raw) = publish_download(&db, "reads")?;
    let (_, manifest) = consuming_quant(&db, &[])?;
    assert!(
        db.record_fastq_retention("RUN1", &manifest, &raw[0])
            .is_err()
    );
    let (_, manifest) = consuming_quant(&db, &raw)?;
    db.record_retained_input("RUN1", &manifest, &raw[0])?;
    fs::remove_file(&raw[0].path)?;
    assert!(!db.verify_generation(OutputGroup::Download, "RUN1")?);
    fs::write(&raw[0].path, "reads")?;
    for artifact in &raw {
        db.record_fastq_retention("RUN1", &manifest, artifact)?;
        fs::remove_file(&artifact.path)?;
        db.delete_artifact(&artifact.path)?;
        assert!(db.verify_generation(OutputGroup::Download, "RUN1")?);
    }
    let manifest_path = head.payload_path.join("manifests/download/RUN1.json");
    let forged = ManifestArtifact::snapshot("FASTQ_SINGLE", &manifest_path)?;
    assert!(
        db.record_fastq_retention("RUN1", &manifest, &forged)
            .is_err()
    );
    fs::remove_file(manifest_path)?;
    assert!(!db.verify_generation(OutputGroup::Download, "RUN1")?);
    Ok(())
}

fn preparing_preprocess(
    db: &Database,
    value: &str,
) -> Result<(JobLease, Generation, Vec<ArtifactRecord>)> {
    use crate::RunState;
    ensure_run(db)?;
    let lease = JobLease::try_acquire(
        db,
        "preprocess",
        Some("RUN1"),
        &[
            RunState::Verified,
            RunState::PreprocessDone,
            RunState::PreprocessFailed,
            RunState::QuantDone,
        ],
    )?
    .unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let generation = Generation::begin(&owned, OutputGroup::Preprocess, "preprocess")?;
    let output = generation
        .paths(&ProjectPaths::new(db.path().parent().unwrap()))
        .unwrap()
        .preprocess_run_dir("RUN1", "kallistors");
    fs::create_dir_all(&output)?;
    let mut records = Vec::new();
    for (kind, name) in [
        (ArtifactKind::FastqTrimmedSingle, "trimmed.fastq"),
        (ArtifactKind::PreprocessReport, "fasterp.json"),
        (ArtifactKind::Manifest, "rnaa_preprocess_manifest.json"),
    ] {
        let path = output.join(name);
        fs::write(&path, value)?;
        let artifact = ManifestArtifact::snapshot(kind.as_str(), &path)?;
        records.push(ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some("RUN1".into()),
            kind,
            path: artifact.path,
            blob_id: None,
            shared_path: None,
            checksum_type: artifact.checksum_type,
            checksum: artifact.checksum,
            bytes: artifact.bytes,
            created_at: now_rfc3339(),
        });
    }
    Ok((lease, generation, records))
}

fn complete_preprocess(db: &Database) -> Result<()> {
    db.set_run_state("RUN1", crate::RunState::PreprocessDone, None)?;
    db.append_event(
        "preprocess",
        Some("RUN1"),
        "preprocess generation completed",
        serde_json::json!({}),
    )
}

fn publish_preprocess(db: &Database, value: &str) -> Result<(GenerationHead, ManifestArtifact)> {
    let (mut lease, mut generation, records) = preparing_preprocess(db, value)?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    generation.commit(&records, complete_preprocess)?;
    lease.finish(JobOutcome::Succeeded, None)?;
    let fastq = ManifestArtifact::snapshot(records[0].kind.as_str(), Path::new(&records[0].path))?;
    Ok((generation.head.clone(), fastq))
}

fn consuming_quant(db: &Database, inputs: &[ManifestArtifact]) -> Result<(GenerationHead, String)> {
    let (mut lease, mut generation, records) = preparing_quant(db, "counts")?;
    generation.seal_with_inputs(&required(&records), &"a".repeat(64), inputs)?;
    generation.commit(&records, complete_quant)?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok((
        generation.head.clone(),
        records
            .iter()
            .find(|r| r.kind == ArtifactKind::Manifest)
            .unwrap()
            .checksum
            .clone(),
    ))
}

#[test]
fn generation_pruning_requires_exact_committed_consumption_and_never_masks_corruption() -> Result<()>
{
    let (_root, db) = setup()?;
    let (producer, fastq) = publish_preprocess(&db, "reads")?;
    let (_, unrelated) = consuming_quant(&db, &[])?;
    assert!(
        db.record_fastq_retention("RUN1", &unrelated, &fastq)
            .is_err()
    );
    let (consumer, manifest) = consuming_quant(&db, std::slice::from_ref(&fastq))?;
    assert!(
        db.record_fastq_retention("RUN1", &"b".repeat(64), &fastq)
            .is_err()
    );
    db.record_retained_input("RUN1", &manifest, &fastq)?;
    fs::remove_file(&fastq.path)?;
    assert!(
        !db.verify_generation(OutputGroup::Preprocess, "RUN1")?,
        "old retention records alone cannot authorize generation pruning"
    );
    fs::write(&fastq.path, "reads")?;
    db.record_fastq_retention("RUN1", &manifest, &fastq)?;
    assert!(db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
    fs::remove_file(&fastq.path)?;
    assert!(db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
    db.record_fastq_retention("RUN1", &manifest, &fastq)?;
    fs::write(&fastq.path, "wrong")?;
    assert!(!db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
    assert!(
        db.record_fastq_retention("RUN1", &manifest, &fastq)
            .is_err()
    );
    fs::remove_file(&fastq.path)?;
    let report = producer
        .payload_path
        .join("preprocess/RUN1/fasterp/fasterp.json");
    let forged = ManifestArtifact::snapshot("FASTQ_TRIMMED_SINGLE", &report)?;
    assert!(
        db.record_fastq_retention("RUN1", &manifest, &forged)
            .is_err(),
        "report kind cannot be relabeled to permit pruning"
    );
    fs::remove_file(&report)?;
    assert!(!db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
    fs::write(&report, "reads")?;
    db.write(|connection| {
        connection.execute(
            "UPDATE generation_pruned_files SET quant_manifest_sha256 = ?1",
            params!["b".repeat(64)],
        )?;
        Ok(())
    })?;
    assert!(!db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
    db.write(|connection| {
        connection.execute(
            "UPDATE generation_pruned_files SET quant_manifest_sha256 = ?1",
            params![manifest],
        )?;
        Ok(())
    })?;
    assert!(db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
    let (next, _) = consuming_quant(&db, &[])?;
    assert_ne!(next.generation_id, consumer.generation_id);
    assert!(
        db.verify_generation(OutputGroup::Preprocess, "RUN1")?,
        "replacement cannot erase historical pruning evidence"
    );
    Ok(())
}

#[test]
fn rejected_qc_publication_and_failed_registration_are_atomic() -> Result<()> {
    let (_root, db) = setup()?;
    let (previous, _) = publish_preprocess(&db, "old")?;
    let (mut lease, mut generation, records) = preparing_preprocess(&db, "new")?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    let failure: Result<()> = generation.commit_outcome(
        &records,
        |_| Ok(()),
        |db| {
            db.set_run_state(
                "RUN1",
                crate::RunState::PreprocessFailed,
                Some("QC rejected"),
            )?;
            db.append_event(
                "preprocess",
                Some("RUN1"),
                "must roll back",
                serde_json::json!({}),
            )?;
            bail!("injected publication failure");
        },
        JobOutcome::Failed,
        Some("QC rejected"),
    );
    assert!(failure.is_err());
    assert_eq!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")?,
        Some(previous)
    );
    assert_eq!(
        db.get_run("RUN1")?.unwrap().state,
        crate::RunState::PreprocessRunning
    );
    assert_eq!(db.list_job_attempts(1)?[0].state, "RUNNING");
    generation.commit_outcome(
        &records,
        |_| Ok(()),
        |db| {
            db.set_run_state(
                "RUN1",
                crate::RunState::PreprocessFailed,
                Some("QC rejected"),
            )?;
            db.append_event(
                "preprocess",
                Some("RUN1"),
                "QC rejected",
                serde_json::json!({}),
            )
        },
        JobOutcome::Failed,
        Some("QC rejected"),
    )?;
    lease.finish(JobOutcome::Failed, None)?;
    assert_eq!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")?,
        Some(generation.head.clone())
    );
    assert_eq!(
        db.get_run("RUN1")?.unwrap().state,
        crate::RunState::PreprocessFailed
    );
    let attempts = db.list_job_attempts(1)?;
    assert_eq!(attempts[0].state, "FAILED");
    assert_eq!(attempts[0].last_error.as_deref(), Some("QC rejected"));
    assert!(db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
    assert!(
        !db.list_events(100)?
            .iter()
            .any(|event| event.message == "must roll back")
    );
    Ok(())
}

#[test]
fn changed_inputs_and_duplicate_inventories_cannot_be_published_or_reused() -> Result<()> {
    let (root, db) = setup()?;
    let previous = publish(&db, "old")?;
    let (mut lease, mut generation, records) = preparing(&db, "new")?;
    let input = root.path().join("source.fastq");
    fs::write(&input, "input")?;
    let snapshot = ManifestArtifact::snapshot("FASTQ_SINGLE", &input)?;
    fs::write(&input, "other")?;
    assert!(
        generation
            .seal_with_inputs(&required(&records), &"a".repeat(64), &[snapshot])
            .is_err()
    );
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?,
        Some(previous.clone())
    );
    drop(generation);
    lease.finish(JobOutcome::Failed, Some("input changed"))?;
    let mut seal = db.current_seal(OutputGroup::De, db.project_id())?.unwrap();
    seal.files.push(seal.files[0].clone());
    let json = serde_json::to_string(&seal)?;
    fs::write(
        previous.payload_path.parent().unwrap().join("seal.json"),
        &json,
    )?;
    db.write(|connection| {
        connection.execute(
            "UPDATE artifact_generations SET inventory_json = ?1 WHERE generation_id = ?2",
            params![json, previous.generation_id],
        )?;
        Ok(())
    })?;
    assert!(
        db.generation_head(OutputGroup::De, db.project_id())
            .is_err()
    );
    assert!(
        ProjectPaths::new(root.path())
            .with_generation_heads(&db)
            .is_err()
    );
    Ok(())
}

fn required(records: &[ArtifactRecord]) -> Vec<PathBuf> {
    records.iter().map(|r| PathBuf::from(&r.path)).collect()
}
fn publish(db: &Database, value: &str) -> Result<GenerationHead> {
    let (mut lease, mut generation, records) = preparing(db, value)?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    generation.commit(&records, |db| {
        db.set_project_stage_state("deseq2", ProjectStageState::Done, None)?;
        db.append_event(
            "deseq2",
            None,
            "generation completed",
            serde_json::json!({}),
        )
    })?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(generation.head.clone())
}
fn assert_head(db: &Database, expected: &GenerationHead, value: &str) -> Result<()> {
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .unwrap()
            .generation_id,
        expected.generation_id
    );
    let records = db.list_artifacts()?;
    assert_eq!(records.len(), 2);
    for record in records {
        assert!(Path::new(&record.path).starts_with(&expected.payload_path));
        assert_eq!(fs::read_to_string(&record.path)?, value);
        assert_eq!(compute_sha256(Path::new(&record.path))?, record.checksum);
    }
    Ok(())
}

#[test]
fn replacing_a_generation_preserves_old_readers_and_retires_old_registry_entries() -> Result<()> {
    let (_root, db) = setup()?;
    let first = publish(&db, "old")?;
    let pinned = ProjectPaths::new(db.path().parent().unwrap()).with_generation_heads(&db)?;
    let next = publish(&db, "new")?;
    assert_ne!(first.generation_id, next.generation_id);
    assert_head(&db, &next, "new")?;
    assert_eq!(
        fs::read_to_string(first.payload_path.join("counts.tsv"))?,
        "old"
    );
    assert_eq!(pinned.de_dir, first.payload_path.join("de"));
    assert_eq!(db.list_job_attempts(2)?[0].state, "SUCCEEDED");
    assert_eq!(
        db.get_project_stage_state("deseq2")?.unwrap().state,
        ProjectStageState::Done
    );
    Ok(())
}

#[test]
fn publication_failure_rolls_back_head_artifacts_state_event_and_attempt() -> Result<()> {
    let (_root, db) = setup()?;
    let first = publish(&db, "old")?;
    let (mut lease, mut generation, records) = preparing(&db, "new")?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    let result: Result<()> = generation.commit(&records, |db| {
        db.set_project_stage_state("deseq2", ProjectStageState::Done, None)?;
        db.append_event("deseq2", None, "must roll back", serde_json::json!({}))?;
        bail!("injected registration failure");
    });
    assert!(result.is_err());
    assert_head(&db, &first, "old")?;
    assert_eq!(
        db.get_project_stage_state("deseq2")?.unwrap().state,
        ProjectStageState::Running
    );
    assert_eq!(db.list_job_attempts(1)?[0].state, "RUNNING");
    assert!(
        !db.list_events(100)?
            .iter()
            .any(|event| event.message == "must roll back")
    );
    drop(generation);
    lease.finish(JobOutcome::Failed, Some("publication failed"))?;
    Ok(())
}

#[test]
fn missing_seals_and_incomplete_catalog_heads_never_count_as_valid_generations() -> Result<()> {
    let (_root, db) = setup()?;
    let head = publish(&db, "valid")?;
    assert!(db.verify_generation(OutputGroup::De, db.project_id())?);
    fs::write(
        head.payload_path.parent().unwrap().join("seal.json"),
        "damaged seal",
    )?;
    assert!(!db.verify_generation(OutputGroup::De, db.project_id())?);
    db.write(|connection| {
        connection.execute(
            "UPDATE artifact_generations SET state = 'SEALED' WHERE generation_id = ?1",
            params![head.generation_id],
        )?;
        Ok(())
    })?;
    assert!(
        db.generation_head(OutputGroup::De, db.project_id())
            .is_err()
    );
    assert!(
        ProjectPaths::new(db.path().parent().unwrap())
            .with_generation_heads(&db)
            .is_err()
    );
    Ok(())
}

#[test]
fn missing_outputs_mutated_seals_and_cancelled_workers_cannot_publish() -> Result<()> {
    let (_root, db) = setup()?;
    let first = publish(&db, "old")?;
    let (mut lease, mut generation, records) = preparing(&db, "new")?;
    let mut incomplete = required(&records);
    incomplete.push(generation.head.payload_path.join("missing.tsv"));
    assert!(generation.seal(&incomplete, &"a".repeat(64)).is_err());
    generation.seal(&required(&records), &"a".repeat(64))?;
    fs::write(&records[0].path, "bad")?;
    assert!(generation.commit(&records, |_| Ok(())).is_err());
    assert_head(&db, &first, "old")?;
    fs::write(&records[0].path, "new")?;
    assert!(db.cancel_job(lease.attempt_id())?);
    assert!(generation.commit(&records, |_| Ok(())).is_err());
    assert_head(&db, &first, "old")?;
    lease.finish(JobOutcome::Failed, None)?;
    drop(generation);
    assert_eq!(db.recover_generations()?, 1);
    assert_eq!(db.list_job_attempts(1)?[0].state, "CANCELLED");
    Ok(())
}

#[test]
fn atomic_contexts_close_and_panics_do_not_leak_partial_writes() -> Result<()> {
    let (_root, db) = setup()?;
    let escaped = db.atomic(|context| {
        context.set_project_stage_state("corr", ProjectStageState::Pending, None)?;
        assert_eq!(
            context.get_project_stage_state("corr")?.unwrap().state,
            ProjectStageState::Pending
        );
        Ok(context.clone())
    })?;
    assert!(
        escaped
            .set_project_stage_state("corr", ProjectStageState::Done, None)
            .is_err()
    );
    let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _: Result<()> = db.atomic(|context| {
            context.set_project_stage_state("corr", ProjectStageState::Done, None)?;
            panic!("injected atomic callback panic");
        });
    }));
    assert!(panic.is_err());
    assert_eq!(
        db.get_project_stage_state("corr")?.unwrap().state,
        ProjectStageState::Pending
    );
    db.set_project_stage_state("corr", ProjectStageState::Failed, None)?;
    Ok(())
}

#[cfg(unix)]
#[test]
fn seals_reject_output_links_outside_the_generation() -> Result<()> {
    let (root, db) = setup()?;
    let (_lease, mut generation, records) = preparing(&db, "new")?;
    let outside = root.path().join("outside");
    fs::write(&outside, "unsafe")?;
    std::os::unix::fs::symlink(&outside, generation.head.payload_path.join("link"))?;
    assert!(
        generation
            .seal(&required(&records), &"a".repeat(64))
            .is_err()
    );
    Ok(())
}

struct KillOnDrop(Child);
impl Drop for KillOnDrop {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
fn signal_and_wait(root: &Path) -> ! {
    fs::write(root.join("crash-ready"), "ready").unwrap();
    loop {
        thread::park_timeout(Duration::from_secs(1));
    }
}
fn kill_stage_worker(db: &Database, root: &Path, worker: &str, phase: &str) -> Result<()> {
    let mut child = KillOnDrop(
        Command::new(std::env::current_exe()?)
            .args(["--exact", worker, "--nocapture"])
            .env("RNAA_TEST_STAGE_GENERATION_ROOT", root)
            .env("RNAA_TEST_GENERATION_PHASE", phase)
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()?,
    );
    let deadline = Instant::now() + Duration::from_secs(10);
    while !root.join("crash-ready").exists() {
        assert!(child.0.try_wait()?.is_none(), "child exited before {phase}");
        assert!(Instant::now() < deadline, "child did not reach {phase}");
        thread::sleep(Duration::from_millis(10));
    }
    child.0.kill()?;
    child.0.wait()?;
    db.write(|connection| {
        connection.execute(
            "UPDATE job_attempts SET lease_expires_ms = 0 WHERE state = 'RUNNING'",
            [],
        )?;
        Ok(())
    })?;
    db.recover_jobs()?;
    db.recover_generations()?;
    Ok(())
}

#[test]
fn download_crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_TEST_STAGE_GENERATION_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let phase = std::env::var("RNAA_TEST_GENERATION_PHASE")?;
    let db = Database::open(root)?;
    let (_lease, mut generation, records) = preparing_download(&db, "new")?;
    if phase == "files" {
        signal_and_wait(root);
    }
    generation.seal(&required(&records), &"a".repeat(64))?;
    if phase == "sealed" {
        signal_and_wait(root);
    }
    generation.commit(&records, |db| {
        complete_download(db)?;
        if phase == "transaction" {
            signal_and_wait(root);
        }
        Ok(())
    })?;
    signal_and_wait(root);
}
#[test]
fn killed_downloads_never_expose_mixed_lane_sets_or_partial_completion() -> Result<()> {
    for phase in ["files", "sealed", "transaction", "committed"] {
        let (root, db) = setup()?;
        let (old, _) = publish_download(&db, "old")?;
        let pinned = ProjectPaths::new(root.path()).with_generation_heads(&db)?;
        kill_stage_worker(
            &db,
            root.path(),
            "generations::tests::download_crash_worker",
            phase,
        )?;
        let head = db.generation_head(OutputGroup::Download, "RUN1")?.unwrap();
        let committed = phase == "committed";
        assert_eq!(head == old, !committed);
        assert!(db.verify_generation(OutputGroup::Download, "RUN1")?);
        let records = db.list_artifacts_for_run(Some("RUN1"))?;
        assert_eq!(records.len(), 3);
        for record in records {
            assert_eq!(
                fs::read_to_string(record.path)?,
                if committed { "new" } else { "old" }
            );
        }
        assert_eq!(
            fs::read_to_string(pinned.raw_run_dir("RUN1").join("lane1.fastq"))?,
            "old"
        );
        assert_eq!(
            db.get_run("RUN1")?.unwrap().state,
            if committed {
                crate::RunState::Verified
            } else {
                crate::RunState::DownloadFailed
            }
        );
        assert_eq!(
            db.list_job_attempts(1)?[0].state,
            if committed {
                "SUCCEEDED"
            } else {
                "INTERRUPTED"
            }
        );
        assert_eq!(
            db.list_events(100)?
                .iter()
                .filter(|event| event.message == "download generation completed")
                .count(),
            if committed { 2 } else { 1 }
        );
    }
    Ok(())
}

#[test]
fn preprocess_crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_TEST_STAGE_GENERATION_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let phase = std::env::var("RNAA_TEST_GENERATION_PHASE")?;
    let db = Database::open(root)?;
    let (_lease, mut generation, records) = preparing_preprocess(&db, "new")?;
    if phase == "files" {
        signal_and_wait(root);
    }
    generation.seal(&required(&records), &"a".repeat(64))?;
    if phase == "sealed" {
        signal_and_wait(root);
    }
    generation.commit(&records, |db| {
        complete_preprocess(db)?;
        if phase == "transaction" {
            signal_and_wait(root);
        }
        Ok(())
    })?;
    signal_and_wait(root);
}

#[test]
fn killed_preprocessing_preserves_complete_generations_and_pinned_readers() -> Result<()> {
    for phase in ["files", "sealed", "transaction", "committed"] {
        let (root, db) = setup()?;
        let (previous, _) = publish_preprocess(&db, "old")?;
        let pinned = ProjectPaths::new(root.path()).with_generation_heads(&db)?;
        kill_stage_worker(
            &db,
            root.path(),
            "generations::tests::preprocess_crash_worker",
            phase,
        )?;
        let current = db
            .generation_head(OutputGroup::Preprocess, "RUN1")?
            .unwrap();
        let committed = phase == "committed";
        assert_eq!(current == previous, !committed);
        assert!(db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
        let records = db.list_artifacts_for_run(Some("RUN1"))?;
        assert_eq!(records.len(), 3);
        for record in records {
            assert_eq!(
                fs::read_to_string(record.path)?,
                if committed { "new" } else { "old" }
            );
        }
        assert_eq!(
            fs::read_to_string(
                pinned
                    .preprocess_run_dir("RUN1", "kallistors")
                    .join("trimmed.fastq")
            )?,
            "old"
        );
        assert_eq!(
            db.get_run("RUN1")?.unwrap().state,
            if committed {
                crate::RunState::PreprocessDone
            } else {
                crate::RunState::PreprocessFailed
            }
        );
        assert_eq!(
            db.list_job_attempts(1)?[0].state,
            if committed {
                "SUCCEEDED"
            } else {
                "INTERRUPTED"
            }
        );
        assert_eq!(
            db.list_events(100)?
                .iter()
                .filter(|event| event.message == "preprocess generation completed")
                .count(),
            if committed { 2 } else { 1 }
        );
    }
    Ok(())
}

#[test]
fn retention_crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_TEST_STAGE_GENERATION_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let phase = std::env::var("RNAA_TEST_GENERATION_PHASE")?;
    let db = Database::open(root)?;
    let _lease = JobLease::try_acquire(
        &db,
        "retention",
        Some("RUN1"),
        &[crate::RunState::QuantDone],
    )?
    .unwrap();
    let owned = _lease.database();
    let record = db
        .list_artifacts_for_run(Some("RUN1"))?
        .into_iter()
        .find(|record| record.kind == ArtifactKind::FastqTrimmedSingle)
        .unwrap();
    let fastq = ManifestArtifact::snapshot(record.kind.as_str(), Path::new(&record.path))?;
    let consumer = db.current_seal(OutputGroup::Quant, "RUN1")?.unwrap();
    owned.record_fastq_retention("RUN1", quant_manifest_digest(&consumer).unwrap(), &fastq)?;
    if phase == "intent" {
        signal_and_wait(root);
    }
    fs::remove_file(&fastq.path)?;
    if phase == "unlinked" {
        signal_and_wait(root);
    }
    owned.delete_artifact(&fastq.path)?;
    signal_and_wait(root);
}

#[test]
fn killed_retention_retries_each_durable_removal_boundary_without_demoting_counts() -> Result<()> {
    for phase in ["intent", "unlinked", "registry"] {
        let (root, db) = setup()?;
        let (_, fastq) = publish_preprocess(&db, "reads")?;
        let (_, manifest) = consuming_quant(&db, std::slice::from_ref(&fastq))?;
        kill_stage_worker(
            &db,
            root.path(),
            "generations::tests::retention_crash_worker",
            phase,
        )?;
        assert_eq!(
            db.get_run("RUN1")?.unwrap().state,
            crate::RunState::QuantDone
        );
        assert_eq!(db.list_job_attempts(1)?[0].state, "INTERRUPTED");
        assert!(db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
        assert!(db.verify_generation(OutputGroup::Quant, "RUN1")?);
        let mut lease = JobLease::try_acquire(
            &db,
            "retention",
            Some("RUN1"),
            &[crate::RunState::QuantDone],
        )?
        .unwrap();
        let owned = lease.database();
        owned.record_fastq_retention("RUN1", &manifest, &fastq)?;
        if Path::new(&fastq.path).exists() {
            fs::remove_file(&fastq.path)?;
        }
        owned.delete_artifact(&fastq.path)?;
        lease.finish(JobOutcome::Succeeded, None)?;
        assert!(db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
        assert!(db.verify_generation(OutputGroup::Quant, "RUN1")?);
        assert_eq!(db.retained_inputs("RUN1")?.len(), 1);
        assert_eq!(
            db.get_run("RUN1")?.unwrap().state,
            crate::RunState::QuantDone
        );
    }
    Ok(())
}

#[test]
fn pruning_authorization_rolls_back_all_evidence_when_registration_fails() -> Result<()> {
    let (_root, db) = setup()?;
    let (_, fastq) = publish_preprocess(&db, "reads")?;
    let (_, manifest) = consuming_quant(&db, std::slice::from_ref(&fastq))?;
    db.write(|connection| {
            connection.execute_batch("CREATE TRIGGER fail_pruning BEFORE INSERT ON generation_pruned_files BEGIN SELECT RAISE(ABORT, 'injected'); END;")?;
            Ok(())
        })?;
    assert!(
        db.record_fastq_retention("RUN1", &manifest, &fastq)
            .is_err()
    );
    assert!(db.retained_inputs("RUN1")?.is_empty());
    assert!(Path::new(&fastq.path).is_file());
    assert!(db.verify_generation(OutputGroup::Preprocess, "RUN1")?);
    Ok(())
}

#[test]
fn quant_crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_TEST_QUANT_GENERATION_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let phase = std::env::var("RNAA_TEST_GENERATION_PHASE")?;
    let db = Database::open(root)?;
    let (_lease, mut generation, records) = preparing_quant(&db, "new")?;
    if phase == "files" {
        signal_and_wait(root);
    }
    generation.seal(&required(&records), &"a".repeat(64))?;
    if phase == "sealed" {
        signal_and_wait(root);
    }
    generation.commit(&records, |db| {
        complete_quant(db)?;
        if phase == "transaction" {
            signal_and_wait(root);
        }
        Ok(())
    })?;
    signal_and_wait(root);
}

#[test]
fn killed_quant_producers_preserve_complete_outputs_and_atomic_run_completion() -> Result<()> {
    for phase in ["files", "sealed", "transaction", "committed"] {
        let (root, db) = setup()?;
        let old = publish_quant(&db, "old")?;
        let pinned = ProjectPaths::new(root.path()).with_generation_heads(&db)?;
        let mut child = KillOnDrop(
            Command::new(std::env::current_exe()?)
                .args([
                    "--exact",
                    "generations::tests::quant_crash_worker",
                    "--nocapture",
                ])
                .env("RNAA_TEST_QUANT_GENERATION_ROOT", root.path())
                .env("RNAA_TEST_GENERATION_PHASE", phase)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .spawn()?,
        );
        let deadline = Instant::now() + Duration::from_secs(10);
        while !root.path().join("crash-ready").exists() {
            assert!(
                child.0.try_wait()?.is_none(),
                "quant child exited before {phase}"
            );
            assert!(
                Instant::now() < deadline,
                "quant child did not reach {phase}"
            );
            thread::sleep(Duration::from_millis(10));
        }
        child.0.kill()?;
        child.0.wait()?;
        db.write(|connection| {
            connection.execute(
                "UPDATE job_attempts SET lease_expires_ms = 0 WHERE state = 'RUNNING'",
                [],
            )?;
            Ok(())
        })?;
        db.recover_jobs()?;
        db.recover_generations()?;
        let head = db.generation_head(OutputGroup::Quant, "RUN1")?.unwrap();
        let current = ProjectPaths::new(root.path()).with_generation_heads(&db)?;
        let expected = if phase == "committed" { "new" } else { "old" };
        assert_eq!(
            head.generation_id == old.generation_id,
            phase != "committed"
        );
        assert!(db.verify_generation(OutputGroup::Quant, "RUN1")?);
        let artifacts = db.list_artifacts_for_run(Some("RUN1"))?;
        assert_eq!(artifacts.len(), 4);
        for record in artifacts {
            assert_eq!(fs::read_to_string(record.path)?, expected);
        }
        assert_eq!(
            fs::read_to_string(
                current
                    .quant_run_dir("RUN1", "kallistors")
                    .join("abundance.h5")
            )?,
            expected
        );
        assert_eq!(
            fs::read_to_string(
                pinned
                    .quant_run_dir("RUN1", "kallistors")
                    .join("abundance.h5")
            )?,
            "old"
        );
        assert_eq!(
            db.get_run("RUN1")?.unwrap().state,
            if phase == "committed" {
                crate::RunState::QuantDone
            } else {
                crate::RunState::QuantFailed
            }
        );
        assert_eq!(
            db.list_job_attempts(1)?[0].state,
            if phase == "committed" {
                "SUCCEEDED"
            } else {
                "INTERRUPTED"
            }
        );
        let completed = db
            .list_events(100)?
            .iter()
            .filter(|event| event.message == "quant generation completed")
            .count();
        assert_eq!(completed, if phase == "committed" { 2 } else { 1 });
    }
    Ok(())
}

#[test]
fn crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_TEST_GENERATION_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let phase = std::env::var("RNAA_TEST_GENERATION_PHASE")?;
    let db = Database::open(root)?;
    let (_lease, mut generation, records) = preparing(&db, "new")?;
    if phase == "files" {
        signal_and_wait(root);
    }
    generation.seal(&required(&records), &"a".repeat(64))?;
    if phase == "sealed" {
        signal_and_wait(root);
    }
    generation.commit(&records, |db| {
        db.set_project_stage_state("deseq2", ProjectStageState::Done, None)?;
        db.append_event(
            "deseq2",
            None,
            "generation completed",
            serde_json::json!({}),
        )?;
        if phase == "transaction" {
            signal_and_wait(root);
        }
        Ok(())
    })?;
    signal_and_wait(root);
}

#[test]
fn killed_processes_expose_one_complete_generation_at_every_publication_boundary() -> Result<()> {
    for phase in ["files", "sealed", "transaction", "committed"] {
        let (root, db) = setup()?;
        let previous = publish(&db, "old")?;
        let mut child = KillOnDrop(
            Command::new(std::env::current_exe()?)
                .args(["--exact", "generations::tests::crash_worker", "--nocapture"])
                .env("RNAA_TEST_GENERATION_ROOT", root.path())
                .env("RNAA_TEST_GENERATION_PHASE", phase)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .spawn()?,
        );
        let deadline = Instant::now() + Duration::from_secs(10);
        while !root.path().join("crash-ready").exists() {
            assert!(
                child.0.try_wait()?.is_none(),
                "child exited before boundary {phase}"
            );
            assert!(
                Instant::now() < deadline,
                "child did not reach boundary {phase}"
            );
            thread::sleep(Duration::from_millis(10));
        }
        child.0.kill()?;
        child.0.wait()?;
        db.write(|connection| {
            connection.execute(
                "UPDATE job_attempts SET lease_expires_ms = 0 WHERE state = 'RUNNING'",
                [],
            )?;
            Ok(())
        })?;
        db.recover_jobs()?;
        db.recover_generations()?;
        if phase == "committed" {
            let head = db
                .generation_head(OutputGroup::De, db.project_id())?
                .unwrap();
            assert_ne!(head.generation_id, previous.generation_id);
            assert_head(&db, &head, "new")?;
            assert_eq!(db.list_job_attempts(1)?[0].state, "SUCCEEDED");
            assert_eq!(
                db.get_project_stage_state("deseq2")?.unwrap().state,
                ProjectStageState::Done
            );
        } else {
            assert_head(&db, &previous, "old")?;
            assert_eq!(db.list_job_attempts(1)?[0].state, "INTERRUPTED");
            let hidden: i64 = db.read(|connection| Ok(connection.query_row("SELECT COUNT(*) FROM artifact_generations WHERE state IN ('BUILDING','SEALED')", [], |row| row.get(0))?))?;
            assert_eq!(hidden, 0);
        }
        assert_eq!(
            fs::read_to_string(previous.payload_path.join("counts.tsv"))?,
            "old"
        );
    }
    Ok(())
}
fn preparing_reference(
    db: &Database,
    value: &str,
) -> Result<(JobLease, Generation, Vec<ArtifactRecord>)> {
    let lease = JobLease::try_acquire_reference(db, "reference_test")?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let generation = Generation::begin(&owned, OutputGroup::Refs, "refs")?;
    let paths = generation
        .paths(&ProjectPaths::new(db.path().parent().unwrap()))
        .unwrap();
    let dir = paths.reference_dir("reference_test");
    fs::create_dir_all(&dir)?;
    let mut outputs = vec![
        (ArtifactKind::ReferenceCdna, dir.join("transcriptome.fa")),
        (ArtifactKind::ReferenceGtf, dir.join("annotations.gtf")),
        (ArtifactKind::ReferenceIndex, dir.join("kallisto.index")),
        (ArtifactKind::ReferenceTx2Gene, dir.join("tx2gene.tsv")),
        (
            ArtifactKind::ReferenceGeneAnnotation,
            dir.join("gene_annotation.tsv"),
        ),
    ];
    outputs.push((
        ArtifactKind::Manifest,
        paths.manifest_path("refs", "reference_test"),
    ));
    let records = outputs
        .into_iter()
        .map(|(kind, path)| {
            fs::write(&path, value)?;
            let artifact = ManifestArtifact::snapshot(kind.as_str(), &path)?;
            Ok(ArtifactRecord {
                project_id: db.project_id().into(),
                run_accession: None,
                kind,
                path: artifact.path,
                blob_id: None,
                shared_path: None,
                checksum_type: artifact.checksum_type,
                checksum: artifact.checksum,
                bytes: artifact.bytes,
                created_at: now_rfc3339(),
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((lease, generation, records))
}
fn complete_reference(db: &Database) -> Result<()> {
    db.append_event(
        "refs",
        None,
        "reference generation completed",
        serde_json::json!({}),
    )
}
#[test]
fn reference_crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_TEST_STAGE_GENERATION_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let db = Database::open(root)?;
    let phase = std::env::var("RNAA_TEST_GENERATION_PHASE")?;
    let (_lease, mut generation, records) = preparing_reference(&db, "new")?;
    if phase == "files" {
        signal_and_wait(root);
    }
    generation.seal(&required(&records), &"a".repeat(64))?;
    if phase == "sealed" {
        signal_and_wait(root);
    }
    generation.commit(&records, |db| {
        complete_reference(db)?;
        if phase == "transaction" {
            signal_and_wait(root);
        }
        Ok(())
    })?;
    signal_and_wait(root);
}
#[test]
fn killed_reference_producers_keep_one_complete_bundle_and_pinned_readers() -> Result<()> {
    for phase in ["files", "sealed", "transaction", "committed"] {
        let (root, db) = setup()?;
        let (mut lease, mut generation, records) = preparing_reference(&db, "old")?;
        generation.seal(&required(&records), &"a".repeat(64))?;
        generation.commit(&records, complete_reference)?;
        let old = generation.head().clone();
        lease.finish(JobOutcome::Succeeded, None)?;
        drop(lease);
        let pinned = old.paths(&ProjectPaths::new(root.path())).unwrap();
        kill_stage_worker(
            &db,
            root.path(),
            "generations::tests::reference_crash_worker",
            phase,
        )?;
        let head = db
            .generation_head(OutputGroup::Refs, "reference_test")?
            .unwrap();
        let committed = phase == "committed";
        assert_eq!(head == old, !committed);
        assert!(db.verify_generation_head(&head)?);
        assert!(db.verify_generation_head(&old)?);
        let active = db.list_artifacts()?;
        assert_eq!(active.len(), 6);
        for record in active {
            assert_eq!(
                fs::read_to_string(record.path)?,
                if committed { "new" } else { "old" }
            );
        }
        assert_eq!(
            fs::read_to_string(
                pinned
                    .reference_dir("reference_test")
                    .join("kallisto.index")
            )?,
            "old"
        );
        assert_eq!(
            db.list_job_attempts(1)?[0].state,
            if committed {
                "SUCCEEDED"
            } else {
                "INTERRUPTED"
            }
        );
        assert_eq!(
            db.list_events(100)?
                .iter()
                .filter(|event| event.message == "reference generation completed")
                .count(),
            if committed { 2 } else { 1 }
        );
    }
    Ok(())
}
#[test]
fn failed_reference_registration_rolls_back_the_entire_replacement() -> Result<()> {
    let (_root, db) = setup()?;
    let (mut lease, mut generation, records) = preparing_reference(&db, "old")?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    generation.commit(&records, complete_reference)?;
    let old = generation.head().clone();
    lease.finish(JobOutcome::Succeeded, None)?;
    drop(lease);
    let (mut lease, mut generation, records) = preparing_reference(&db, "new")?;
    generation.seal(&required(&records), &"a".repeat(64))?;
    assert!(
        generation
            .commit(&records, |db| -> Result<()> {
                complete_reference(db)?;
                bail!("injected reference failure")
            })
            .is_err()
    );
    lease.finish(JobOutcome::Failed, Some("injected failure"))?;
    assert_eq!(
        db.generation_head(OutputGroup::Refs, "reference_test")?
            .unwrap(),
        old
    );
    assert!(db.verify_generation_head(&old)?);
    assert_eq!(db.list_artifacts()?.len(), 6);
    assert_eq!(
        db.list_events(100)?
            .iter()
            .filter(|event| event.message == "reference generation completed")
            .count(),
        1
    );
    Ok(())
}
