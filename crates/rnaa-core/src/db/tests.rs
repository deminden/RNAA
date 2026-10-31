use super::*;
use crate::{ArtifactRecord, LibraryLayout, ProjectStageState, RunRecord, SharedBlobRecord};
use std::collections::BTreeMap;
use std::fs;
use tempfile::TempDir;

#[test]
fn run_state_transitions_are_persisted() {
    let temp = TempDir::new().unwrap();
    let config = ProjectConfig::default();
    let db = Database::create(temp.path(), &config).unwrap();
    let run = RunRecord {
        project_id: db.project_id().to_string(),
        study_accession: "SRP1".to_string(),
        source_study_accession: Some("PRJNA1".to_string()),
        run_accession: "SRR1".to_string(),
        sample_accession: Some("SRS1".to_string()),
        experiment_accession: Some("SRX1".to_string()),
        library_layout: LibraryLayout::Paired,
        instrument_platform: Some("ILLUMINA".to_string()),
        instrument_model: Some("HiSeq".to_string()),
        remote_files: vec![],
        metadata: BTreeMap::new(),
        state: RunState::Resolved,
        last_error: None,
        updated_at: now_rfc3339(),
    };
    db.upsert_run(&run).unwrap();

    db.set_run_state("SRR1", RunState::Downloading, None)
        .unwrap();
    db.set_run_state("SRR1", RunState::Verified, None).unwrap();

    let persisted = db.get_run("SRR1").unwrap().unwrap();
    assert_eq!(persisted.state, RunState::Verified);
}

#[test]
fn upsert_run_preserves_existing_state_for_resolve_refresh() {
    let temp = TempDir::new().unwrap();
    let config = ProjectConfig::default();
    let db = Database::create(temp.path(), &config).unwrap();
    let mut run = RunRecord {
        project_id: db.project_id().to_string(),
        study_accession: "SRP1".to_string(),
        source_study_accession: Some("PRJNA1".to_string()),
        run_accession: "SRR1".to_string(),
        sample_accession: Some("SRS1".to_string()),
        experiment_accession: Some("SRX1".to_string()),
        library_layout: LibraryLayout::Paired,
        instrument_platform: Some("ILLUMINA".to_string()),
        instrument_model: Some("HiSeq".to_string()),
        remote_files: vec![],
        metadata: BTreeMap::from([("condition".to_string(), "A".to_string())]),
        state: RunState::Resolved,
        last_error: None,
        updated_at: now_rfc3339(),
    };
    db.upsert_run(&run).unwrap();
    db.set_run_state("SRR1", RunState::QuantDone, None).unwrap();

    run.metadata
        .insert("condition".to_string(), "B".to_string());
    run.state = RunState::Resolved;
    run.updated_at = now_rfc3339();
    db.upsert_run(&run).unwrap();

    let persisted = db.get_run("SRR1").unwrap().unwrap();
    assert_eq!(persisted.state, RunState::QuantDone);
    assert_eq!(
        persisted.metadata.get("condition").map(String::as_str),
        Some("B")
    );
}

#[test]
fn normalized_run_metadata_is_persisted_and_queryable() {
    let temp = TempDir::new().unwrap();
    let config = ProjectConfig::default();
    let db = Database::create(temp.path(), &config).unwrap();
    let run = RunRecord {
        project_id: db.project_id().to_string(),
        study_accession: "SRP1".to_string(),
        source_study_accession: Some("PRJNA1".to_string()),
        run_accession: "SRR1".to_string(),
        sample_accession: Some("SRS1".to_string()),
        experiment_accession: Some("SRX1".to_string()),
        library_layout: LibraryLayout::Paired,
        instrument_platform: Some("ILLUMINA".to_string()),
        instrument_model: Some("HiSeq".to_string()),
        remote_files: vec![],
        metadata: BTreeMap::from([
            ("condition".to_string(), "A".to_string()),
            ("batch".to_string(), "B1".to_string()),
            (
                "ena_library_source".to_string(),
                "TRANSCRIPTOMIC".to_string(),
            ),
        ]),
        state: RunState::Resolved,
        last_error: None,
        updated_at: now_rfc3339(),
    };
    db.upsert_run(&run).unwrap();

    let persisted = db.get_run("SRR1").unwrap().unwrap();
    assert_eq!(
        persisted.metadata.get("condition").map(String::as_str),
        Some("A")
    );
    assert_eq!(
        persisted.metadata.get("batch").map(String::as_str),
        Some("B1")
    );
    assert_eq!(
        db.list_metadata_columns().unwrap(),
        vec![
            "batch".to_string(),
            "condition".to_string(),
            "ena_library_source".to_string()
        ]
    );
}

#[test]
fn legacy_metadata_json_is_backfilled_into_normalized_table_on_open() {
    let temp = TempDir::new().unwrap();
    let config = ProjectConfig::default();
    let db = Database::create(temp.path(), &config).unwrap();
    let run = RunRecord {
        project_id: db.project_id().to_string(),
        study_accession: "SRP1".to_string(),
        source_study_accession: Some("PRJNA1".to_string()),
        run_accession: "SRR1".to_string(),
        sample_accession: Some("SRS1".to_string()),
        experiment_accession: Some("SRX1".to_string()),
        library_layout: LibraryLayout::Paired,
        instrument_platform: Some("ILLUMINA".to_string()),
        instrument_model: Some("HiSeq".to_string()),
        remote_files: vec![],
        metadata: BTreeMap::from([
            ("condition".to_string(), "A".to_string()),
            ("donor".to_string(), "D1".to_string()),
        ]),
        state: RunState::Resolved,
        last_error: None,
        updated_at: now_rfc3339(),
    };
    db.upsert_run(&run).unwrap();

    let conn = Connection::open(db.path()).unwrap();
    conn.execute(
        "DELETE FROM run_metadata WHERE project_id = ?1 AND run_accession = ?2",
        params![db.project_id(), "SRR1"],
    )
    .unwrap();
    let deleted_count: i64 = conn
        .query_row(
            "SELECT COUNT(*) FROM run_metadata WHERE run_accession = ?1",
            params!["SRR1"],
            |row| row.get(0),
        )
        .unwrap();
    assert_eq!(deleted_count, 0);

    let reopened = Database::open(temp.path()).unwrap();
    let persisted = reopened.get_run("SRR1").unwrap().unwrap();
    assert_eq!(
        persisted.metadata.get("condition").map(String::as_str),
        Some("A")
    );
    assert_eq!(
        persisted.metadata.get("donor").map(String::as_str),
        Some("D1")
    );
    assert_eq!(
        reopened.list_metadata_columns().unwrap(),
        vec!["condition".to_string(), "donor".to_string()]
    );
    let mut malformed = run.clone();
    malformed.run_accession = "SRR2".into();
    db.upsert_run(&malformed).unwrap();
    conn.execute("DELETE FROM run_metadata", []).unwrap();
    conn.execute(
        "UPDATE runs SET metadata_json = 'invalid JSON' WHERE run_accession = 'SRR2'",
        [],
    )
    .unwrap();
    assert!(Database::open(temp.path()).is_err());
    let copied: i64 = conn
        .query_row("SELECT COUNT(*) FROM run_metadata", [], |row| row.get(0))
        .unwrap();
    assert_eq!(
        copied, 0,
        "a failed backfill must roll back earlier valid rows"
    );
}

#[test]
fn artifact_blob_linkage_is_persisted() {
    let temp = TempDir::new().unwrap();
    let config = ProjectConfig::default();
    let db = Database::create(temp.path(), &config).unwrap();
    let artifact_path = temp.path().join("quant").join("SRR1").join("abundance.h5");
    fs::create_dir_all(artifact_path.parent().unwrap()).unwrap();
    fs::write(&artifact_path, b"test").unwrap();

    let blob = SharedBlobRecord {
        blob_id: "sha256:abc123".to_string(),
        storage_path: "/tmp/shared/blobs/sha256/ab/c1/abc123".to_string(),
        checksum_type: "sha256".to_string(),
        checksum: "abc123".to_string(),
        bytes: 4,
        created_at: now_rfc3339(),
    };
    db.record_shared_blob(&blob).unwrap();
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().to_string(),
        run_accession: Some("SRR1".to_string()),
        kind: crate::model::ArtifactKind::QuantAbundanceH5,
        path: artifact_path.display().to_string(),
        blob_id: Some(blob.blob_id.clone()),
        shared_path: Some(blob.storage_path.clone()),
        checksum_type: "sha256".to_string(),
        checksum: "abc123".to_string(),
        bytes: 4,
        created_at: now_rfc3339(),
    })
    .unwrap();

    let artifacts = db.list_artifacts_for_run(Some("SRR1")).unwrap();
    assert_eq!(artifacts.len(), 1);
    assert_eq!(artifacts[0].blob_id.as_deref(), Some("sha256:abc123"));
    assert_eq!(
        artifacts[0].shared_path.as_deref(),
        Some("/tmp/shared/blobs/sha256/ab/c1/abc123")
    );

    let loaded_blob = db.get_shared_blob("sha256:abc123").unwrap().unwrap();
    assert_eq!(loaded_blob.checksum, "abc123");
}

#[test]
fn project_stage_states_are_persisted() {
    let temp = TempDir::new().unwrap();
    let config = ProjectConfig::default();
    let db = Database::create(temp.path(), &config).unwrap();

    db.set_project_stage_state("normalize", ProjectStageState::Running, None)
        .unwrap();
    db.set_project_stage_state("normalize", ProjectStageState::Failed, Some("interrupted"))
        .unwrap();

    let stage = db.get_project_stage_state("normalize").unwrap().unwrap();
    assert_eq!(stage.state, ProjectStageState::Failed);
    assert_eq!(stage.last_error.as_deref(), Some("interrupted"));
}
