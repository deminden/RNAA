use super::*;
use crate::{
    generations::Generation,
    jobs::{JobLease, JobOutcome},
};
mod helpers;
use helpers::*;

#[test]
fn verified_manifest_import_preserves_original_seals_and_unblocks_dependency_gc() -> Result<()> {
    let (root, db) = setup()?;
    let source = root.path().join("external.tsv");
    fs::write(&source, "source")?;
    let external = ManifestArtifact::snapshot("SOURCE", &source)?;
    let first = publish(&db, "first", std::slice::from_ref(&external))?;
    let second = publish(&db, "second", &[input(&first)?])?;
    let old_first = downgrade(&db, &first, true)?;
    let old_second = downgrade(&db, &second, true)?;
    let current = publish(&db, "unrelated current", &[external])?;
    assert!(
        db.generation_gc(gc(false))?
            .entries
            .iter()
            .all(|row| row.action == "held_reachable")
    );
    let preview = db.import_generation_lineage(false)?;
    assert!(preview.complete, "{:?}", preview.errors);
    assert_eq!(
        preview
            .entries
            .iter()
            .filter(|row| row.action == "import")
            .count(),
        2
    );
    assert!(!root.path().join(".rnaa_support/lineage_imports").exists());
    assert!(!db.generation_lineage(&second)?.recorded);
    let report = db.import_generation_lineage(true)?;
    assert!(report.errors.is_empty(), "{:?}", report.errors);
    assert_eq!(
        fs::read(first.payload_path.parent().unwrap().join("seal.json"))?,
        old_first
    );
    assert_eq!(
        fs::read(second.payload_path.parent().unwrap().join("seal.json"))?,
        old_second
    );
    let lineage = db.generation_lineage(&second)?;
    assert!(lineage.recorded);
    assert_eq!(lineage.source, "imported");
    assert_eq!(
        lineage.dependencies[0].producer_generation_id,
        first.generation_id
    );
    assert_eq!(db.dependent_generations(&first)?, vec![second.clone()]);
    let report = db.generation_gc(gc(true))?;
    assert!(report.errors.is_empty(), "{:?}", report.errors);
    assert!(
        report
            .entries
            .iter()
            .filter(|row| row.generation_id != current.generation_id)
            .all(|row| row.state == "TRASHED")
    );
    db.generation_restore(&second.generation_id)?;
    assert!(db.verify_generation_head(&first)?);
    assert!(db.verify_generation_head(&second)?);
    assert_eq!(db.generation_lineage(&second)?.dependencies.len(), 1);
    assert!(
        db.import_generation_lineage(true)?
            .entries
            .iter()
            .filter(|row| row.generation_id != current.generation_id)
            .all(|row| row.action == "already_imported")
    );
    Ok(())
}

#[test]
fn corrupt_witnesses_never_authorize_legacy_collection() -> Result<()> {
    let (root, db) = setup()?;
    let source = root.path().join("source");
    fs::write(&source, "input")?;
    let first = publish(
        &db,
        "first",
        &[ManifestArtifact::snapshot("SOURCE", &source)?],
    )?;
    downgrade(&db, &first, true)?;
    let manifest = first.payload_path.join("rnaa_manifest.json");
    let original = fs::read(&manifest)?;
    fs::write(&manifest, "corruption")?;
    let report = db.import_generation_lineage(true)?;
    assert!(!report.complete);
    assert!(report.entries.iter().all(|row| row.action == "held"));
    fs::write(&manifest, original)?;
    db.import_generation_lineage(true)?;
    let imported = db
        .imported_lineage(&db.committed_seal(&first.generation_id)?.unwrap())?
        .unwrap();
    fs::write(db.lineage_import_path(&imported.import_id, false)?, b"{}")?;
    assert!(!db.verify_generation_head(&first)?);
    assert!(db.generation_lineage(&first).is_err());
    assert!(db.generation_gc(gc(true)).is_err());
    Ok(())
}

#[test]
fn importing_historical_provenance_does_not_expand_the_consumed_input_set() -> Result<()> {
    let (root, db) = setup()?;
    let source = root.path().join("source");
    fs::write(&source, "source")?;
    let external = ManifestArtifact::snapshot("SOURCE", &source)?;
    let first = publish(&db, "first", std::slice::from_ref(&external))?;
    let other = publish(&db, "other", &[external])?;
    let consumed = input(&first)?;
    let ancestry = input(&other)?;
    let result =
        publish_with_provenance(&db, "result", std::slice::from_ref(&consumed), &[ancestry])?;
    downgrade(&db, &result, false)?;
    assert!(db.import_generation_lineage(true)?.complete);
    assert_eq!(
        db.generation_inputs(OutputGroup::De, db.project_id())?
            .unwrap(),
        vec![consumed]
    );
    let lineage = db.generation_lineage(&result)?;
    assert_eq!(lineage.dependencies.len(), 2);
    assert!(
        lineage
            .dependencies
            .iter()
            .any(|edge| edge.producer_generation_id == other.generation_id)
    );
    assert!(
        db.generation_gc(gc(true))?
            .entries
            .iter()
            .all(|entry| entry.action == "held_reachable")
    );
    Ok(())
}

#[test]
fn omitted_reference_inputs_recover_the_original_local_producer_from_sealed_source_descriptors()
-> Result<()> {
    let (root, db) = setup()?;
    let external = root.path().join("source");
    fs::write(&external, "source")?;
    let producer = publish(
        &db,
        "reference source",
        &[ManifestArtifact::snapshot("SOURCE", &external)?],
    )?;
    let origin = producer.payload_path.join("counts.tsv");
    let mut lease = JobLease::try_acquire_reference(&db, "REF")?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"b".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::Refs, "refs")?;
    let cdna = generation.head().payload_path.join("transcripts.fa");
    fs::copy(&origin, &cdna)?;
    let gtf = generation.head().payload_path.join("genes.gtf");
    fs::copy(&external, &gtf)?;
    let inputs = vec![
        ManifestArtifact::snapshot("REFERENCE_CDNA", &cdna)?,
        ManifestArtifact::snapshot("REFERENCE_GTF", &gtf)?,
    ];
    let manifest = generation.head().payload_path.join("manifest.json");
    write_json_pretty(
        &manifest,
        &StageManifest {
            schema_version: 1,
            stage: "refs".into(),
            id: "REF".into(),
            started_at: now_rfc3339(),
            finished_at: now_rfc3339(),
            parameters: serde_json::json!({"transcript_source":{"kind":"LocalFile","location":origin,"sha256":compute_sha256(&origin)?},"gtf_source":{"kind":"LocalFile","location":external,"sha256":compute_sha256(&external)?}}),
            tool_versions: BTreeMap::new(),
            input_artifacts: inputs.clone(),
            output_artifacts: inputs,
            notes: Vec::new(),
        },
    )?;
    let records = [
        (ArtifactKind::ReferenceCdna, cdna.clone()),
        (ArtifactKind::ReferenceGtf, gtf.clone()),
        (ArtifactKind::Manifest, manifest.clone()),
    ]
    .into_iter()
    .map(|(kind, path)| -> Result<ArtifactRecord> {
        let file = ManifestArtifact::snapshot(kind.as_str(), &path)?;
        Ok(ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: None,
            kind,
            path: file.path,
            blob_id: None,
            shared_path: None,
            checksum_type: file.checksum_type,
            checksum: file.checksum,
            bytes: file.bytes,
            created_at: now_rfc3339(),
        })
    })
    .collect::<Result<Vec<_>>>()?;
    generation.seal_with_inputs(
        &[cdna, gtf, manifest],
        &"b".repeat(64),
        &[
            ManifestArtifact::snapshot("REFERENCE_CDNA", &origin)?,
            ManifestArtifact::snapshot("REFERENCE_GTF", &external)?,
        ],
    )?;
    let head = generation.head().clone();
    generation.commit(&records, |_| Ok(()))?;
    lease.finish(JobOutcome::Succeeded, None)?;
    drop(generation);
    drop(owned);
    drop(lease);
    downgrade(&db, &head, true)?;
    assert!(db.import_generation_lineage(true)?.complete);
    let lineage = db.generation_lineage(&head)?;
    assert!(
        lineage
            .dependencies
            .iter()
            .any(|edge| edge.producer_generation_id == producer.generation_id
                && edge.input.path == origin.display().to_string())
    );
    publish(
        &db,
        "new unrelated head",
        &[ManifestArtifact::snapshot("SOURCE", &external)?],
    )?;
    assert_eq!(
        db.generation_gc(gc(true))?
            .entries
            .iter()
            .find(|entry| entry.generation_id == producer.generation_id)
            .unwrap()
            .action,
        "held_reachable"
    );
    Ok(())
}

#[test]
fn unknown_empty_declarations_keep_history_held_and_cycles_are_rejected() -> Result<()> {
    let (root, db) = setup()?;
    let unknown = publish(&db, "unknown", &[])?;
    downgrade(&db, &unknown, true)?;
    let report = db.import_generation_lineage(true)?;
    assert!(!report.complete);
    assert!(report.entries.iter().all(|row| row.action == "held"));
    let source = root.path().join("source");
    fs::write(&source, "source")?;
    let later = publish(
        &db,
        "later",
        &[ManifestArtifact::snapshot("SOURCE", &source)?],
    )?;
    assert!(
        db.generation_gc(gc(true))?
            .entries
            .iter()
            .all(|row| row.action == "held_reachable")
    );
    assert!(later.payload_path.exists());
    let graph = BTreeMap::from([
        ("a".into(), vec!["b".into()]),
        ("b".into(), vec!["a".into()]),
    ]);
    assert!(verify_acyclic(&graph).is_err());
    Ok(())
}

#[test]
fn pins_block_import_and_transaction_failure_publishes_no_partial_graph() -> Result<()> {
    let (root, db) = setup()?;
    let source = root.path().join("source");
    fs::write(&source, "source")?;
    let first = publish(
        &db,
        "first",
        &[ManifestArtifact::snapshot("SOURCE", &source)?],
    )?;
    let second = publish(&db, "second", &[input(&first)?])?;
    downgrade(&db, &first, true)?;
    downgrade(&db, &second, true)?;
    let pin = super::super::GenerationStorePin::acquire(db.path(), "test")?;
    assert!(db.import_generation_lineage(true)?.busy);
    drop(pin);
    db.write(|conn|{conn.execute_batch("CREATE TRIGGER fail_import BEFORE INSERT ON generation_import_dependencies BEGIN SELECT RAISE(ABORT,'injected import failure'); END;")?;Ok(())})?;
    assert!(db.import_generation_lineage(true).is_err());
    let count: i64 = db.read(|conn| {
        Ok(conn.query_row(
            "SELECT count(*) FROM generation_lineage_imports",
            [],
            |row| row.get(0),
        )?)
    })?;
    assert_eq!(count, 0);
    assert!(
        db.generation_gc(gc(false))?
            .entries
            .iter()
            .all(|row| row.action == "held_reachable")
    );
    db.write(|conn| {
        conn.execute_batch("DROP TRIGGER fail_import;")?;
        Ok(())
    })?;
    assert!(db.import_generation_lineage(true)?.complete);
    Ok(())
}

#[test]
fn import_crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_LINEAGE_IMPORT_TEST_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let db = Database::open(root)?;
    fs::write(root.join("ready"), b"ready")?;
    while !root.join("go").exists() {
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    db.import_generation_lineage(true)?;
    fs::write(root.join("published"), b"published")?;
    loop {
        std::thread::sleep(std::time::Duration::from_secs(1));
    }
}

#[test]
fn killed_importers_leave_only_unpublished_witnesses_or_one_complete_committed_graph() -> Result<()>
{
    struct Child(std::process::Child);
    impl Drop for Child {
        fn drop(&mut self) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
    for before_catalog in [true, false] {
        let (root, db) = setup()?;
        let source = root.path().join("source");
        fs::write(&source, "source")?;
        let head = publish(
            &db,
            "counts",
            &[ManifestArtifact::snapshot("SOURCE", &source)?],
        )?;
        let original = downgrade(&db, &head, true)?;
        let mut child = Child(
            std::process::Command::new(std::env::current_exe()?)
                .args([
                    "--exact",
                    "generations::legacy::tests::import_crash_worker",
                    "--nocapture",
                ])
                .env("RNAA_LINEAGE_IMPORT_TEST_ROOT", root.path())
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::piped())
                .spawn()?,
        );
        let start = std::time::Instant::now();
        while !root.path().join("ready").exists() {
            assert!(child.0.try_wait()?.is_none());
            assert!(start.elapsed() < std::time::Duration::from_secs(15));
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        let connection = rusqlite::Connection::open(db.path())?;
        if before_catalog {
            connection.execute_batch("BEGIN IMMEDIATE")?;
        }
        fs::write(root.path().join("go"), b"go")?;
        let witness_directory = root.path().join(".rnaa_support/lineage_imports");
        loop {
            let reached = if before_catalog {
                witness_directory.exists() && fs::read_dir(&witness_directory)?.next().is_some()
            } else {
                root.path().join("published").exists()
            };
            if reached {
                break;
            }
            assert!(child.0.try_wait()?.is_none());
            assert!(start.elapsed() < std::time::Duration::from_secs(15));
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        if before_catalog {
            assert!(db.import_generation_lineage(false)?.busy);
        }
        child.0.kill()?;
        child.0.wait()?;
        if before_catalog {
            connection.execute_batch("ROLLBACK")?;
        }
        let report = db.import_generation_lineage(true)?;
        assert!(report.complete, "{:?}", report.errors);
        assert_eq!(
            report.entries[0].action,
            if before_catalog {
                "imported"
            } else {
                "already_imported"
            }
        );
        assert!(db.verify_generation_head(&head)?);
        assert_eq!(
            fs::read(head.payload_path.parent().unwrap().join("seal.json"))?,
            original
        );
        let count: i64 = db.read(|conn| {
            Ok(conn.query_row(
                "SELECT count(*) FROM generation_lineage_imports",
                [],
                |row| row.get(0),
            )?)
        })?;
        assert_eq!(count, 1);
    }
    Ok(())
}
