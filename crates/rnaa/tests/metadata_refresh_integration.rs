//! Editable metadata survives ingest/refresh and reaches real model fitting.
mod support;
use rnaa_core::{
    ProjectConfig,
    generations::OutputGroup,
    jobs::{JobLease, JobOutcome},
    metadata::{ColumnSpec, ColumnType},
};
use std::fs;
use support::native::{command, setup, success};

#[test]
fn metadata_refresh_preserves_edits_through_ingest_and_native_fitting_and_fails_closed_on_conflicts()
 {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["metadata", "refresh"]);
    let first = db.metadata_publication().unwrap().unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(paths.samplesheet_path())
        .unwrap();
    let mut headers = reader
        .headers()
        .unwrap()
        .iter()
        .map(str::to_string)
        .collect::<Vec<_>>();
    let condition = headers.iter().position(|name| name == "condition").unwrap();
    headers.push("age".into());
    let mut rows = reader
        .records()
        .map(|r| r.unwrap().iter().map(str::to_string).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    drop(reader);
    for (i, row) in rows.iter_mut().enumerate() {
        row[condition] = if i < 4 { "control" } else { "treated" }.into();
        row.push([20, 22, 21, 23, 21, 23, 20, 22][i].to_string());
    }
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(paths.samplesheet_path())
        .unwrap();
    writer.write_record(headers).unwrap();
    for row in rows {
        writer.write_record(row).unwrap();
    }
    writer.flush().unwrap();
    drop(writer);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config
        .metadata
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    config.save(&paths.config_path).unwrap();
    success(root, &["metadata", "refresh"]);
    success(root, &["download"]);
    success(
        root,
        &[
            "deseq2",
            "--design",
            "~ age + condition",
            "--contrast",
            "condition",
            "treated",
            "control",
        ],
    );
    let de = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let contract: serde_json::Value = serde_json::from_slice(
        &fs::read(
            de.paths(&paths)
                .unwrap()
                .de_project_dir(db.project_id())
                .join("model_design.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(contract["metadata"]["columns"]["age"]["type"], "numeric");
    assert_eq!(
        contract["metadata"]["columns"]["condition"]["levels"],
        serde_json::json!(["control", "treated"])
    );
    let manifest: rnaa_core::StageManifest = serde_json::from_slice(
        &fs::read(
            de.paths(&paths)
                .unwrap()
                .de_project_dir(db.project_id())
                .join("deseq2_manifest.json"),
        )
        .unwrap(),
    )
    .unwrap();
    let frozen_source = manifest
        .input_artifacts
        .iter()
        .find(|a| a.kind == "SAMPLESHEET_SOURCE")
        .unwrap();
    assert!(std::path::Path::new(&frozen_source.path).starts_with(&de.payload_path));
    let editable = fs::read_to_string(paths.samplesheet_path()).unwrap();
    let annotated = editable
        .lines()
        .enumerate()
        .map(|(i, line)| {
            format!(
                "{line}\t{}\n",
                if i == 0 {
                    "notes"
                } else {
                    "retained annotation"
                }
            )
        })
        .collect::<String>();
    fs::write(paths.samplesheet_path(), annotated).unwrap();
    assert!(
        manifest.verify_local_artifacts().unwrap(),
        "later metadata edits must not corrupt old analysis source witnesses"
    );
    let old = db.metadata_publication().unwrap().unwrap();
    let view = fs::read(paths.samplesheet_path()).unwrap();
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.metadata
        .insert("condition".into(), "source-update".into());
    db.upsert_run(&run).unwrap();
    let failed = command(root, &["metadata", "refresh"]);
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("refresh conflicts"));
    assert_eq!(
        db.metadata_publication().unwrap().unwrap().publication_id,
        old.publication_id
    );
    assert_eq!(fs::read(paths.samplesheet_path()).unwrap(), view);
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap(),
        de
    );
    success(root, &["metadata", "refresh", "--conflicts", "keep-local"]);
    assert_eq!(fs::read(paths.samplesheet_path()).unwrap(), view);
    first.verify().unwrap();
    old.verify().unwrap();
    success(root, &["metadata", "show"]);
}

#[test]
fn refresh_obeys_analysis_ownership_and_analysis_recovers_a_committed_uninstalled_view() {
    use rnaa_core::samplesheet::{
        load_or_init_column_map,
        publication::{ConflictPolicy, PreparedPublication},
    };
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["metadata", "refresh"]);
    let first = db.metadata_publication().unwrap().unwrap();
    let mut lease = JobLease::try_acquire(&db, "normalize", None, &[])
        .unwrap()
        .unwrap();
    let rejected = command(root, &["metadata", "refresh"]);
    assert!(!rejected.status.success());
    success(root, &["download"]);
    assert_eq!(
        db.metadata_publication().unwrap().unwrap().publication_id,
        first.publication_id
    );
    lease.finish(JobOutcome::Interrupted, None).unwrap();
    drop(lease);
    let mut runs = db.list_active_runs().unwrap();
    runs[0].metadata.insert("condition".into(), "B".into());
    let map = load_or_init_column_map(&paths.column_map_path(), &[]).unwrap();
    let pending = PreparedPublication::prepare(
        &db,
        &paths,
        &runs,
        &db.list_artifacts().unwrap(),
        &map,
        ConflictPolicy::Fail,
    )
    .unwrap();
    pending.commit_catalog(&db, &paths).unwrap();
    success(root, &["normalize"]);
    let de = paths
        .with_generation_heads(&db)
        .unwrap()
        .de_project_dir(db.project_id());
    let model: serde_json::Value =
        serde_json::from_slice(&fs::read(de.join("model_design.json")).unwrap()).unwrap();
    let values = model["metadata"]["columns"]["condition"]["values"]
        .as_array()
        .unwrap();
    assert_eq!(values.iter().filter(|v| *v == "B").count(), 5);
    assert_eq!(values.iter().filter(|v| *v == "A").count(), 3);
    assert!(
        fs::read_to_string(paths.samplesheet_path())
            .unwrap()
            .contains("RUN1")
    );
}

#[test]
fn adding_same_named_overrides_never_replaces_a_previously_registered_source() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path().join("project");
    success(&root, &["init"]);
    let left = temp.path().join("left");
    let right = temp.path().join("right");
    fs::create_dir_all(&left).unwrap();
    fs::create_dir_all(&right).unwrap();
    let a = left.join("metadata.tsv");
    let b = right.join("metadata.tsv");
    fs::write(&a, "run_accession\tcondition\nRUN1\tcontrol\n").unwrap();
    fs::write(&b, "run_accession\tcondition\nRUN1\ttreated\n").unwrap();
    success(
        &root,
        &[
            "add",
            "--metadata",
            a.to_str().unwrap(),
            "--metadata",
            b.to_str().unwrap(),
        ],
    );
    let db = rnaa_core::Database::open(&root).unwrap();
    let records = db.list_metadata_overrides().unwrap();
    assert_eq!(records.len(), 2);
    assert_ne!(records[0].stored_path, records[1].stored_path);
    assert_eq!(
        fs::read(&records[0].stored_path).unwrap(),
        fs::read(a).unwrap()
    );
    assert_eq!(
        fs::read(&records[1].stored_path).unwrap(),
        fs::read(b).unwrap()
    );
}
