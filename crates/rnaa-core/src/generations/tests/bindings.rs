use super::*;
use crate::fingerprint::{InputIdentity, StageRequest};
use serde_json::{Value, json};

fn request(stage: &str, inputs: &[ManifestArtifact]) -> StageRequest<Value> {
    StageRequest {
        schema_version: 1,
        stage: stage.into(),
        algorithm_version: 1,
        engines: Default::default(),
        parameters: json!({"science":"fixture"}),
        inputs: inputs.iter().map(InputIdentity::from).collect(),
        dependencies: Default::default(),
    }
}
fn record(db: &Database, kind: ArtifactKind, path: &Path) -> Result<ArtifactRecord> {
    let snapshot = ManifestArtifact::snapshot(kind.as_str(), path)?;
    Ok(ArtifactRecord {
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
    })
}
fn publish(
    db: &Database,
    group: OutputGroup,
    inputs: &[ManifestArtifact],
    value: &str,
) -> Result<GenerationHead> {
    ensure_run(db)?;
    let run = db.get_run("RUN1")?.unwrap();
    let stage = group.as_str();
    let mut lease = JobLease::try_acquire(db, stage, Some("RUN1"), &[run.state])?.unwrap();
    let owned = lease.database();
    let request = request(stage, inputs);
    owned.set_job_request(&request.digest()?)?;
    let mut generation = Generation::begin(&owned, group, stage)?;
    let paths = generation
        .paths(&ProjectPaths::new(db.path().parent().unwrap()))
        .unwrap();
    let outputs = if group == OutputGroup::Download {
        vec![
            (
                ArtifactKind::FastqSingle,
                paths.raw_run_dir("RUN1").join("lane.fastq"),
            ),
            (
                ArtifactKind::Manifest,
                paths.manifest_path("download", "RUN1"),
            ),
        ]
    } else if group == OutputGroup::Preprocess {
        let dir = paths.preprocess_run_dir("RUN1", "kallistors");
        vec![
            (ArtifactKind::FastqTrimmedSingle, dir.join("trimmed.fastq")),
            (ArtifactKind::PreprocessReport, dir.join("fasterp.json")),
        ]
    } else {
        let dir = paths.quant_run_dir("RUN1", "kallistors");
        vec![
            (ArtifactKind::QuantAbundanceTsv, dir.join("abundance.tsv")),
            (ArtifactKind::Manifest, dir.join("rnaa_quant_manifest.json")),
        ]
    };
    let mut records = Vec::new();
    for (kind, path) in outputs {
        fs::create_dir_all(path.parent().unwrap())?;
        fs::write(&path, value)?;
        records.push(record(db, kind, &path)?);
    }
    generation.seal_request(
        &records
            .iter()
            .map(|record| PathBuf::from(&record.path))
            .collect::<Vec<_>>(),
        &request,
        inputs,
    )?;
    let head = generation.head().clone();
    generation.commit(&records, |db| {
        db.set_run_state("RUN1", crate::RunState::QuantDone, None)
    })?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(head)
}
fn input(db: &Database, head: &GenerationHead) -> Result<ManifestArtifact> {
    Ok(db
        .generation_output_artifacts(head)?
        .into_iter()
        .find(|file| file.kind == "FASTQ_SINGLE")
        .unwrap())
}

#[test]
fn full_requests_publish_atomically_and_policy_lookup_retains_exact_input_identity() -> Result<()> {
    let (_temp, db) = setup()?;
    let raw = publish(&db, OutputGroup::Download, &[], "first")?;
    let source = input(&db, &raw)?;
    let first = publish(
        &db,
        OutputGroup::Quant,
        std::slice::from_ref(&source),
        "counts",
    )?;
    let saved = db.scientific_request(&first.request_sha256)?.unwrap();
    assert_eq!(saved.digest()?, first.request_sha256);
    let witness = first.payload_path.parent().unwrap().join("seal.json");
    let seal_bytes = fs::read(&witness)?;
    db.write(|conn| {
        conn.execute(
            "DELETE FROM scientific_requests WHERE project_id=?1 AND request_sha256=?2",
            params![db.project_id(), first.request_sha256],
        )?;
        Ok(())
    })?;
    assert!(db.scientific_request(&first.request_sha256)?.is_none());
    db.index_generation_request(&first, &serde_json::to_value(&saved)?)?;
    assert_eq!(fs::read(witness)?, seal_bytes);
    let mut policy = saved.clone();
    policy.inputs.clear();
    assert_eq!(
        db.generation_policy_candidates(OutputGroup::Quant, "quant", "RUN1", &policy.digest()?)?,
        vec![first.clone()]
    );
    assert!(
        db.generation_policy_candidates(OutputGroup::Download, "quant", "RUN1", &policy.digest()?)?
            .is_empty()
    );
    let mut changed = saved.clone();
    changed.inputs[0].sha256 = "b".repeat(64);
    assert!(
        db.index_generation_request(&first, &serde_json::to_value(changed)?)
            .is_err()
    );
    let mut unknown = serde_json::to_value(&saved)?;
    unknown["unrecognized"] = json!(true);
    assert!(db.index_generation_request(&first, &unknown).is_err());
    // A failed publication leaves no request object falsely linked to a result.
    let run = db.get_run("RUN1")?.unwrap();
    let mut lease = JobLease::try_acquire(&db, "quant", Some("RUN1"), &[run.state])?.unwrap();
    let owned = lease.database();
    let mut next = request("quant", std::slice::from_ref(&source));
    next.algorithm_version = 2;
    owned.set_job_request(&next.digest()?)?;
    let mut generation = Generation::begin(&owned, OutputGroup::Quant, "quant")?;
    let path = generation.head().payload_path.join("result");
    fs::write(&path, "next")?;
    let records = vec![record(&db, ArtifactKind::QuantAbundanceTsv, &path)?];
    generation.seal_request(std::slice::from_ref(&path), &next, &[source])?;
    let failed: Result<()> = generation.commit(&records, |_| bail!("rollback request publication"));
    assert!(failed.is_err());
    assert!(db.scientific_request(&next.digest()?)?.is_none());
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")?.unwrap(),
        first
    );
    lease.finish(JobOutcome::Failed, None)?;
    Ok(())
}

#[test]
fn equivalent_new_paths_bind_and_authorize_pruning_without_editing_original_count_seals()
-> Result<()> {
    let (_temp, db) = setup()?;
    let raw_a = publish(&db, OutputGroup::Download, &[], "reads")?;
    let input_a = input(&db, &raw_a)?;
    let counts = publish(
        &db,
        OutputGroup::Quant,
        std::slice::from_ref(&input_a),
        "counts",
    )?;
    let witness = counts.payload_path.parent().unwrap().join("seal.json");
    let bytes = fs::read(&witness)?;
    let raw_b = publish(&db, OutputGroup::Download, &[], "reads")?;
    let input_b = input(&db, &raw_b)?;
    assert_eq!(InputIdentity::from(&input_a), InputIdentity::from(&input_b));
    let run = db.get_run("RUN1")?.unwrap();
    let mut lease = JobLease::try_acquire(&db, "quant", Some("RUN1"), &[run.state])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&counts.request_sha256)?;
    let verified = owned.prepare_generation_reuse(&counts)?.unwrap();
    let verified = owned.bind_generation_reuse(verified, std::slice::from_ref(&input_b))?;
    owned.reattach_generation(verified, |db| {
        db.set_run_state("RUN1", crate::RunState::QuantDone, None)
    })?;
    lease.finish(JobOutcome::Succeeded, None)?;
    drop(lease);
    let binding = db
        .generation_input_binding(OutputGroup::Quant, "RUN1")?
        .unwrap();
    assert_eq!(binding.inputs, vec![input_b.clone()]);
    let manifest = compute_sha256(
        &counts
            .paths(&ProjectPaths::new(db.path().parent().unwrap()))
            .unwrap()
            .quant_run_dir("RUN1", "kallistors")
            .join("rnaa_quant_manifest.json"),
    )?;
    db.record_fastq_retention("RUN1", &manifest, &input_b)?;
    fs::remove_file(&input_b.path)?;
    db.delete_artifact(&input_b.path)?;
    assert!(db.verify_generation_head(&raw_b)?);
    assert!(db.verified_retained_inputs("RUN1")?.contains(&input_b));
    assert!(
        db.quant_retention_evidence("RUN1")?
            .iter()
            .any(|item| item.artifact == input_b && item.generation_sha256 == manifest)
    );
    assert_eq!(fs::read(witness)?, bytes);
    assert!(db.verify_generation_head(&raw_a)?);
    // A pruned FASTQ remains absent from the restored producer's active registry.
    let run = db.get_run("RUN1")?.unwrap();
    let mut download = JobLease::try_acquire(&db, "download", Some("RUN1"), &[run.state])?.unwrap();
    let owned = download.database();
    owned.set_job_request(&raw_b.request_sha256)?;
    let verified = owned.prepare_generation_reuse(&raw_b)?.unwrap();
    assert_eq!(verified.artifact_count(), 1);
    owned.reattach_generation(verified, |_| Ok(()))?;
    download.finish(JobOutcome::Succeeded, None)?;
    assert!(
        !db.list_artifacts_for_run(Some("RUN1"))?
            .iter()
            .any(|record| record.path == input_b.path)
    );
    Ok(())
}

#[test]
fn binding_rejects_changed_content_roles_multiplicity_and_unowned_missing_inputs() -> Result<()> {
    let (_temp, db) = setup()?;
    let raw = publish(&db, OutputGroup::Download, &[], "reads")?;
    let source = input(&db, &raw)?;
    let counts = publish(
        &db,
        OutputGroup::Quant,
        std::slice::from_ref(&source),
        "counts",
    )?;
    let run = db.get_run("RUN1")?.unwrap();
    let mut lease = JobLease::try_acquire(&db, "quant", Some("RUN1"), &[run.state])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&counts.request_sha256)?;
    let mut role = source.clone();
    role.kind = "FASTQ_R1".into();
    let mut content = source.clone();
    content.checksum = "a".repeat(64);
    for inputs in [
        vec![role],
        vec![content],
        vec![source.clone(), source.clone()],
        vec![],
    ] {
        assert!(
            owned
                .bind_generation_reuse(owned.prepare_generation_reuse(&counts)?.unwrap(), &inputs)
                .is_err()
        );
    }
    fs::remove_file(&source.path)?;
    assert!(
        owned
            .bind_generation_reuse(owned.prepare_generation_reuse(&counts)?.unwrap(), &[source])
            .is_err()
    );
    assert!(
        db.generation_input_binding(OutputGroup::Quant, "RUN1")?
            .is_none()
    );
    lease.finish(JobOutcome::Failed, None)?;
    Ok(())
}

#[test]
fn failed_binding_selection_rolls_back_and_corrupt_receipts_are_not_accepted() -> Result<()> {
    let (_temp, db) = setup()?;
    let raw = publish(&db, OutputGroup::Download, &[], "reads")?;
    let source = input(&db, &raw)?;
    let counts = publish(
        &db,
        OutputGroup::Quant,
        std::slice::from_ref(&source),
        "counts",
    )?;
    let run = db.get_run("RUN1")?.unwrap();
    let mut lease = JobLease::try_acquire(&db, "quant", Some("RUN1"), &[run.state])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&counts.request_sha256)?;
    let verified = owned.bind_generation_reuse(
        owned.prepare_generation_reuse(&counts)?.unwrap(),
        std::slice::from_ref(&source),
    )?;
    let failed: Result<()> =
        owned.reattach_generation(verified, |_| bail!("injected binding rollback"));
    assert!(failed.is_err());
    assert!(
        db.generation_input_binding(OutputGroup::Quant, "RUN1")?
            .is_none()
    );
    let verified = owned
        .bind_generation_reuse(owned.prepare_generation_reuse(&counts)?.unwrap(), &[source])?;
    owned.reattach_generation(verified, |_| Ok(()))?;
    lease.finish(JobOutcome::Succeeded, None)?;
    let binding = db
        .generation_input_binding(OutputGroup::Quant, "RUN1")?
        .unwrap();
    db.write(|conn| {
        conn.execute(
            "UPDATE generation_input_bindings SET binding_sha256=?1 WHERE binding_id=?2",
            params!["0".repeat(64), binding.binding_id],
        )?;
        Ok(())
    })?;
    assert!(
        db.generation_input_binding(OutputGroup::Quant, "RUN1")
            .is_err()
    );
    Ok(())
}

#[test]
fn originals_outside_a_trimmed_only_key_require_the_exact_preprocessing_ancestry() -> Result<()> {
    let (_temp, db) = setup()?;
    let raw_a = publish(&db, OutputGroup::Download, &[], "reads")?;
    let input_a = input(&db, &raw_a)?;
    let preprocessed = publish(&db, OutputGroup::Preprocess, &[input_a], "trimmed reads")?;
    let trimmed = db
        .generation_output_artifacts(&preprocessed)?
        .into_iter()
        .find(|file| file.kind == "FASTQ_TRIMMED_SINGLE")
        .unwrap();
    let counts = publish(
        &db,
        OutputGroup::Quant,
        std::slice::from_ref(&trimmed),
        "trimmed-only fit",
    )?;
    let wrong = publish(&db, OutputGroup::Download, &[], "other reads")?;
    let wrong = input(&db, &wrong)?;
    let equivalent = publish(&db, OutputGroup::Download, &[], "reads")?;
    let equivalent_input = input(&db, &equivalent)?;
    let run = db.get_run("RUN1")?.unwrap();
    let mut lease = JobLease::try_acquire(&db, "quant", Some("RUN1"), &[run.state])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&counts.request_sha256)?;
    assert!(
        owned
            .bind_generation_reuse(
                owned.prepare_generation_reuse(&counts)?.unwrap(),
                &[wrong, trimmed.clone()]
            )
            .is_err()
    );
    let verified = owned.bind_generation_reuse(
        owned.prepare_generation_reuse(&counts)?.unwrap(),
        &[equivalent_input.clone(), trimmed],
    )?;
    owned.reattach_generation(verified, |db| {
        db.set_run_state("RUN1", crate::RunState::QuantDone, None)
    })?;
    lease.finish(JobOutcome::Succeeded, None)?;
    let manifest = compute_sha256(
        &counts
            .paths(&ProjectPaths::new(db.path().parent().unwrap()))
            .unwrap()
            .quant_run_dir("RUN1", "kallistors")
            .join("rnaa_quant_manifest.json"),
    )?;
    db.record_fastq_retention("RUN1", &manifest, &equivalent_input)?;
    fs::remove_file(&equivalent_input.path)?;
    db.delete_artifact(&equivalent_input.path)?;
    assert!(db.verify_generation_head(&equivalent)?);
    Ok(())
}
