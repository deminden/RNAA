use super::*;
use crate::{
    ArtifactKind, ArtifactRecord, ProjectPaths, ProjectStageState,
    generations::Generation,
    jobs::{JobLease, JobOutcome},
};
use std::fs;
fn setup() -> Result<(tempfile::TempDir, Database, ProjectConfig)> {
    let root = tempfile::tempdir()?;
    let mut config = ProjectConfig::default();
    config.refs.custom_cdna = root.path().join("source.fa").display().to_string();
    config.refs.custom_gtf = root.path().join("source.gtf").display().to_string();
    fs::write(&config.refs.custom_cdna, "fasta")?;
    fs::write(&config.refs.custom_gtf, "gtf")?;
    config.save(&root.path().join("rnaa.toml"))?;
    let db = Database::create(root.path(), &config)?;
    Ok((root, db, config))
}
fn candidate(
    db: &Database,
    config: &ProjectConfig,
    id: &str,
    index: &str,
    mapping: &str,
) -> Result<(
    JobLease,
    Generation,
    Vec<ArtifactRecord>,
    ReferenceSelection,
)> {
    let lease = JobLease::try_acquire_reference(db, id)?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::Refs, "refs")?;
    let paths = generation
        .paths(&ProjectPaths::new(db.path().parent().unwrap()))
        .unwrap();
    let dir = paths.reference_dir(id);
    fs::create_dir_all(&dir)?;
    let bundle = ReferenceBundle {
        id: id.into(),
        organism: "human".into(),
        ensembl_release: "custom".into(),
        cdna_path: dir.join("source.fa"),
        gtf_path: dir.join("source.gtf"),
        kallisto_index_path: dir.join("kallisto.index"),
        tx2gene_path: dir.join("tx2gene.tsv"),
        gene_annotation_path: Some(dir.join("gene_annotation.tsv")),
        manifest_path: paths.manifest_path("refs", id),
    };
    let files = [
        (ArtifactKind::ReferenceCdna, &bundle.cdna_path, "fasta"),
        (ArtifactKind::ReferenceGtf, &bundle.gtf_path, "gtf"),
        (
            ArtifactKind::ReferenceIndex,
            &bundle.kallisto_index_path,
            index,
        ),
        (
            ArtifactKind::ReferenceTx2Gene,
            &bundle.tx2gene_path,
            mapping,
        ),
        (
            ArtifactKind::ReferenceGeneAnnotation,
            bundle.gene_annotation_path.as_ref().unwrap(),
            mapping,
        ),
        (ArtifactKind::Manifest, &bundle.manifest_path, "manifest"),
    ];
    let mut records = Vec::new();
    let mut required = Vec::new();
    for (kind, path, data) in files {
        fs::write(path, data)?;
        let artifact = ManifestArtifact::snapshot(kind.as_str(), path)?;
        required.push(path.clone());
        records.push(ArtifactRecord {
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
        });
    }
    generation.seal(&required, &"a".repeat(64))?;
    let selection = ReferenceSelection::capture(&bundle, generation.head().clone(), config)?;
    Ok((lease, generation, records, selection))
}
fn publish(
    db: &Database,
    config: &ProjectConfig,
    id: &str,
    index: &str,
    mapping: &str,
    expected: Option<&ReferenceChoiceToken>,
) -> Result<ReferenceSelection> {
    let (mut lease, mut generation, records, selection) =
        candidate(db, config, id, index, mapping)?;
    generation.commit(&records, |db| db.select_reference(&selection, expected))?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(selection)
}
#[test]
fn selection_scopes_preserve_counts_for_annotation_changes_and_fence_settings_changes() -> Result<()>
{
    let (_root, db, mut config) = setup()?;
    let first = publish(&db, &config, "reference", "index A", "mapping A", None)?;
    first.verify_current(&db, &config, ReferenceScope::Analysis)?;
    let old = first.selection_token()?;
    db.set_project_stage_state("corr", ProjectStageState::Done, None)?;
    let second = publish(
        &db,
        &config,
        "reference",
        "index A",
        "mapping B",
        Some(&old),
    )?;
    assert_eq!(
        first.fingerprint(ReferenceScope::Index)?,
        second.fingerprint(ReferenceScope::Index)?
    );
    assert_ne!(
        first.fingerprint(ReferenceScope::Analysis)?,
        second.fingerprint(ReferenceScope::Analysis)?
    );
    first.verify_current(&db, &config, ReferenceScope::Index)?;
    assert!(
        first
            .verify_catalog(&db, &config, ReferenceScope::Analysis)
            .is_err()
    );
    assert_eq!(
        db.get_project_stage_state("corr")?.unwrap().state,
        ProjectStageState::Pending
    );
    config.refs.index_threads = 3;
    config.save(&db.path().parent().unwrap().join("rnaa.toml"))?;
    db.set_project_config(&config)?;
    second.verify_current(&db, &config, ReferenceScope::Analysis)?;
    config.refs.index_k = 29;
    db.set_project_config(&config)?;
    assert!(
        second
            .verify_catalog(&db, &config, ReferenceScope::Index)
            .is_err()
    );
    Ok(())
}
#[test]
fn selection_and_invalidation_roll_back_with_failed_publication() -> Result<()> {
    let (_root, db, config) = setup()?;
    let first = publish(&db, &config, "reference", "index A", "mapping A", None)?;
    db.set_project_stage_state("corr", ProjectStageState::Done, None)?;
    let events = db.list_events(100)?.len();
    let expected = first.selection_token()?;
    let (mut lease, mut generation, records, second) =
        candidate(&db, &config, "reference", "index B", "mapping B")?;
    let failed: Result<()> = generation.commit(&records, |db| {
        db.select_reference(&second, Some(&expected))?;
        bail!("injected failure")
    });
    assert!(failed.is_err());
    lease.finish(JobOutcome::Failed, Some("injected failure"))?;
    assert_eq!(db.selected_reference()?.unwrap(), first);
    assert_eq!(
        db.generation_head(OutputGroup::Refs, "reference")?.unwrap(),
        first.generation
    );
    assert_eq!(
        db.get_project_stage_state("corr")?.unwrap().state,
        ProjectStageState::Done
    );
    assert_eq!(db.list_events(100)?.len(), events + 2);
    Ok(())
}
#[test]
fn stale_producers_cannot_replace_a_newer_selection_or_forge_an_artifact_role() -> Result<()> {
    let (_root, db, config) = setup()?;
    let first = publish(&db, &config, "first", "index A", "mapping A", None)?;
    let expected = first.selection_token()?;
    let second = publish(
        &db,
        &config,
        "second",
        "index B",
        "mapping A",
        Some(&expected),
    )?;
    let (mut lease, mut generation, records, third) =
        candidate(&db, &config, "third", "index C", "mapping A")?;
    assert!(
        generation
            .commit(&records, |db| db.select_reference(&third, Some(&expected)))
            .is_err()
    );
    lease.finish(JobOutcome::Failed, Some("stale selection"))?;
    assert_eq!(db.selected_reference()?.unwrap(), second);
    assert!(db.generation_head(OutputGroup::Refs, "third")?.is_none());
    let (mut equivalent_lease, mut equivalent, records, choice) =
        candidate(&db, &config, "equivalent", "index B", "mapping A")?;
    assert!(
        equivalent
            .commit(&records, |db| db.select_reference(&choice, Some(&expected)))
            .is_err(),
        "equal data does not authorize an old producer to replace newer provenance"
    );
    equivalent_lease.finish(JobOutcome::Failed, Some("stale equivalent selection"))?;
    assert_eq!(db.selected_reference()?.unwrap(), second);
    let mut forged = second.clone();
    forged.index = forged.cdna.clone();
    assert!(
        db.select_reference(&forged, Some(&second.selection_token()?))
            .is_err()
    );
    assert_eq!(db.selected_reference()?.unwrap(), second);
    Ok(())
}
#[test]
fn source_mutation_and_corrupt_descriptor_are_not_silently_accepted() -> Result<()> {
    let (_root, db, config) = setup()?;
    let first = publish(&db, &config, "reference", "index A", "mapping A", None)?;
    let original = serde_json::to_string(&first)?;
    let mut invalid = serde_json::to_value(&first)?;
    invalid["unsupported_policy"] = serde_json::json!(true);
    db.write(|connection| {
        connection.execute(
            "UPDATE reference_selection SET descriptor_json=?1",
            [serde_json::to_string(&invalid)?],
        )?;
        Ok(())
    })?;
    assert!(db.selected_reference().is_err());
    db.write(|connection| {
        connection.execute(
            "UPDATE reference_selection SET descriptor_json=?1",
            [original],
        )?;
        Ok(())
    })?;
    fs::write(&config.refs.custom_cdna, "changed")?;
    assert!(
        first
            .verify_current(&db, &config, ReferenceScope::Index)
            .is_err()
    );
    db.write(|connection| {
        connection.execute("UPDATE reference_selection SET fingerprint='corrupt'", [])?;
        Ok(())
    })?;
    assert!(db.selected_reference().is_err());
    Ok(())
}

#[test]
fn an_analysis_cannot_publish_after_a_concurrent_reference_selection() -> Result<()> {
    let (_root, db, config) = setup()?;
    let first = publish(&db, &config, "first", "index A", "mapping A", None)?;
    let mut parent = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    let owned = parent.database();
    owned.set_job_request(&"b".repeat(64))?;
    let mut analysis = Generation::begin(&owned, OutputGroup::De, "normalize")?;
    let path = analysis.head().payload_path.join("de/result.tsv");
    fs::write(&path, "result")?;
    let snapshot = ManifestArtifact::snapshot("COUNTS", &path)?;
    let record = ArtifactRecord {
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
    };
    analysis.seal(&[path], &"b".repeat(64))?;
    let second = publish(
        &owned,
        &config,
        "second",
        "index B",
        "mapping B",
        Some(&first.selection_token()?),
    )?;
    let failed: Result<()> = analysis.commit(&[record], |db| {
        first.verify_catalog(db, &config, ReferenceScope::Analysis)
    });
    assert!(failed.is_err());
    parent.finish(JobOutcome::Failed, Some("reference changed"))?;
    assert!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .is_none()
    );
    assert_eq!(db.selected_reference()?.unwrap(), second);
    assert!(
        !db.list_artifacts()?
            .iter()
            .any(|record| record.kind == ArtifactKind::Counts)
    );
    Ok(())
}
