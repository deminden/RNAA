use super::*;

fn publish(db: &Database, value: &str, inputs: &[ManifestArtifact]) -> Result<GenerationHead> {
    let (mut lease, mut generation, records) = preparing(db, value)?;
    generation.seal_with_inputs(
        &records
            .iter()
            .map(|r| PathBuf::from(&r.path))
            .collect::<Vec<_>>(),
        &"a".repeat(64),
        inputs,
    )?;
    let head = generation.head().clone();
    generation.commit(&records, |_| Ok(()))?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(head)
}
fn input(head: &GenerationHead) -> Result<ManifestArtifact> {
    ManifestArtifact::snapshot("COUNTS", &head.payload_path.join("counts.tsv"))
}

#[test]
fn exact_edges_and_transitive_history_survive_head_replacement_and_detect_corruption() -> Result<()>
{
    let (temp, db) = setup()?;
    let first = publish(&db, "first", &[])?;
    let external = temp.path().join("external.tsv");
    fs::write(&external, "explicit external input")?;
    let external = ManifestArtifact::snapshot("SOURCE", &external)?;
    let second = publish(&db, "second", &[input(&first)?, external.clone()])?;
    let third = publish(&db, "third", &[input(&second)?])?;
    let lineage = db.generation_lineage(&second)?;
    assert!(lineage.recorded);
    assert_eq!(lineage.dependencies.len(), 1);
    assert_eq!(
        lineage.dependencies[0].producer_generation_id,
        first.generation_id
    );
    assert_eq!(lineage.external_inputs, vec![external]);
    let descendants = db.dependent_generations(&first)?;
    assert_eq!(descendants.len(), 2);
    assert!(descendants.contains(&third));
    assert!(db.verify_generation_head(&first)?);
    db.write(|conn| {
        conn.execute(
            "UPDATE generation_dependencies SET role='ALTERED' WHERE consumer_generation_id=?1",
            [&second.generation_id],
        )?;
        Ok(())
    })?;
    assert!(!db.verify_generation_head(&second)?);
    assert!(db.generation_lineage(&second).is_err());
    db.write(|conn| {
        conn.execute(
            "UPDATE generation_dependencies SET role='COUNTS' WHERE consumer_generation_id=?1",
            [&second.generation_id],
        )?;
        Ok(())
    })?;
    assert!(db.verify_generation_head(&second)?);
    let original = input(&first)?;
    db.write(|conn| {
        conn.execute(
            "UPDATE generation_files SET sha256=?1 WHERE generation_id=?2 AND path=?3",
            params!["0".repeat(64), first.generation_id, original.path],
        )?;
        Ok(())
    })?;
    assert!(!db.verify_generation_head(&first)?);
    assert!(!db.verify_generation_head(&second)?);
    db.write(|conn| {
        conn.execute(
            "UPDATE generation_files SET sha256=?1 WHERE generation_id=?2 AND path=?3",
            params![original.checksum, first.generation_id, original.path],
        )?;
        Ok(())
    })?;
    assert!(db.verify_generation_head(&first)?);
    assert!(db.verify_generation_head(&third)?);
    let original_json: String = db.read(|conn| {
        Ok(conn.query_row(
            "SELECT inventory_json FROM artifact_generations WHERE generation_id=?1",
            [&first.generation_id],
            |row| row.get(0),
        )?)
    })?;
    let mut altered: serde_json::Value = serde_json::from_str(&original_json)?;
    altered["head"]["request_sha256"] = serde_json::json!("b".repeat(64));
    db.write(|conn| {
        conn.execute(
            "UPDATE artifact_generations SET inventory_json=?1 WHERE generation_id=?2",
            params![serde_json::to_string(&altered)?, first.generation_id],
        )?;
        Ok(())
    })?;
    assert!(
        !db.verify_generation_head(&second)?,
        "edge binds the canonical producer seal as well as file bytes"
    );
    db.write(|conn| {
        conn.execute(
            "UPDATE artifact_generations SET inventory_json=?1 WHERE generation_id=?2",
            params![original_json, first.generation_id],
        )?;
        Ok(())
    })?;
    assert!(db.verify_generation_head(&second)?);
    Ok(())
}

#[test]
fn changed_dependency_catalog_rolls_back_consumer_head_registry_and_edge_publication() -> Result<()>
{
    let (_temp, db) = setup()?;
    let first = publish(&db, "first", &[])?;
    let before = db.list_artifacts()?.len();
    let (mut lease, mut candidate, records) = preparing(&db, "candidate")?;
    candidate.seal_with_inputs(
        &records
            .iter()
            .map(|r| PathBuf::from(&r.path))
            .collect::<Vec<_>>(),
        &"a".repeat(64),
        &[input(&first)?],
    )?;
    db.write(|conn| {
        conn.execute(
            "UPDATE artifact_generations SET state='ABANDONED' WHERE generation_id=?1",
            [&first.generation_id],
        )?;
        Ok(())
    })?;
    assert!(candidate.commit(&records, |_| Ok(())).is_err());
    db.write(|conn| {
        conn.execute(
            "UPDATE artifact_generations SET state='COMMITTED' WHERE generation_id=?1",
            [&first.generation_id],
        )?;
        Ok(())
    })?;
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .unwrap(),
        first
    );
    assert_eq!(db.list_artifacts()?.len(), before);
    assert_eq!(
        db.read(|conn| Ok(conn.query_row(
            "SELECT COUNT(*) FROM generation_dependencies WHERE consumer_generation_id=?1",
            [&candidate.head().generation_id],
            |row| row.get::<_, i64>(0)
        )?))?,
        0
    );
    lease.finish(JobOutcome::Failed, Some("injected stale dependency"))?;
    Ok(())
}

#[test]
fn uncommitted_managed_sources_cannot_be_reclassified_as_external_inputs() -> Result<()> {
    let (_temp, db) = setup()?;
    let (mut producer_lease, producer, records) = preparing_quant(&db, "building")?;
    let source = ManifestArtifact::snapshot("ABUNDANCE", Path::new(&records[0].path))?;
    let (mut lease, mut consumer, outputs) = preparing(&db, "consumer")?;
    let error = consumer
        .seal_with_inputs(
            &outputs
                .iter()
                .map(|r| PathBuf::from(&r.path))
                .collect::<Vec<_>>(),
            &"a".repeat(64),
            &[source],
        )
        .unwrap_err();
    assert!(error.to_string().contains("lacks a committed producer"));
    #[cfg(unix)]
    {
        let alias = db.path().parent().unwrap().join("building-source-alias");
        std::os::unix::fs::symlink(&records[0].path, &alias)?;
        let source = ManifestArtifact::snapshot("ABUNDANCE", &alias)?;
        let failed = consumer.seal_with_inputs(
            &outputs
                .iter()
                .map(|r| PathBuf::from(&r.path))
                .collect::<Vec<_>>(),
            &"a".repeat(64),
            &[source],
        );
        assert!(failed.unwrap_err().to_string().contains("canonical"));
    }
    assert!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .is_none()
    );
    lease.finish(JobOutcome::Failed, None)?;
    drop(consumer);
    producer_lease.finish(JobOutcome::Failed, None)?;
    drop(producer);
    Ok(())
}

#[test]
fn migration_backfills_file_ownership_without_rewriting_legacy_seals() -> Result<()> {
    let (temp, db) = setup()?;
    let first = publish(&db, "legacy", &[])?;
    let witness = first.payload_path.parent().unwrap().join("seal.json");
    let mut legacy: serde_json::Value = serde_json::from_slice(&fs::read(&witness)?)?;
    legacy.as_object_mut().unwrap().remove("lineage_version");
    let bytes = serde_json::to_vec(&legacy)?;
    fs::write(&witness, &bytes)?;
    db.write(|conn| {
        conn.execute("UPDATE artifact_generations SET inventory_json=?1 WHERE generation_id=?2", params![String::from_utf8(bytes.clone())?,first.generation_id])?;
        conn.execute_batch("DROP TABLE generation_dependencies; DROP TABLE generation_files; DELETE FROM schema_migrations WHERE version='0014_generation_lineage';")?;
        Ok(())
    })?;
    let upgraded = Database::open(temp.path())?;
    assert_eq!(fs::read(witness)?, bytes);
    assert!(!upgraded.generation_lineage(&first)?.recorded);
    assert!(upgraded.verify_generation_head(&first)?);
    let second = publish(&upgraded, "new", &[input(&first)?])?;
    assert_eq!(
        upgraded.generation_lineage(&second)?.dependencies[0].producer_generation_id,
        first.generation_id
    );
    assert!(upgraded.verify_generation_head(&second)?);
    Ok(())
}
