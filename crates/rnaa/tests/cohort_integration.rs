//! Active selections must constrain every consumer without deleting old results.
mod support;
use rnaa_core::{InputType, ProjectStageState, RunState, generations::OutputGroup};
use rnaa_formats::matrix::read_matrix_tsv;
use support::native::{command, setup, success};

#[test]
fn narrower_cohort_blocks_old_matrix_and_rebuilds_only_selected_samples() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let previous = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let runs = db.list_runs().unwrap();
    let active = runs
        .into_iter()
        .filter(|run| !["RUN1", "RUN8"].contains(&run.run_accession.as_str()))
        .collect::<Vec<_>>();
    db.publish_resolved_cohort(&active, &db.selector_contract().unwrap(), None)
        .unwrap();
    for id in ["RUN1", "RUN8"] {
        db.set_run_state(id, RunState::DownloadFailed, Some("historical failure"))
            .unwrap();
    }
    assert_eq!(
        db.get_project_stage_state("normalize")
            .unwrap()
            .unwrap()
            .state,
        ProjectStageState::Pending
    );
    let rejected = command(root, &["corr", "--adjust", "none"]);
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("cohort or QC policy is stale"));
    success(root, &["normalize"]);
    let current = paths.with_generation_heads(&db).unwrap();
    let matrix = read_matrix_tsv(&current.de_project_dir(db.project_id()).join("vst.tsv")).unwrap();
    assert_eq!(
        matrix.col_ids,
        ["BIO1", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6"]
    );
    assert_eq!(db.list_runs().unwrap().len(), 8);
    assert_eq!(db.list_active_runs().unwrap().len(), 6);
    assert!(db.verify_generation_head(&previous).unwrap());
    let status = command(root, &["status"]);
    assert!(status.status.success());
    let status = String::from_utf8_lossy(&status.stdout);
    assert!(status.contains("runs_total\t6"));
    assert!(status.contains("runs_catalog_total\t8"));
    assert!(status.contains("cohort_member\tRUN2\tBIO1\tSampleAccession"));
    assert!(!status.contains("DOWNLOAD_FAILED\t2"));
    let mut historical = db.get_run("RUN1").unwrap().unwrap();
    historical.run_accession = "RUN9".into();
    historical.state = RunState::QuantDone;
    db.upsert_run(&historical).unwrap();
    success(root, &["status"]);
    assert!(
        db.eligibility("RUN9").unwrap().is_none(),
        "inactive history must not cause unrelated QC invalidation"
    );
    assert_eq!(
        db.get_project_stage_state("normalize")
            .unwrap()
            .unwrap()
            .state,
        ProjectStageState::Done
    );
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let count = db.list_job_attempts(1000).unwrap().len();
    success(root, &["download"]);
    assert_eq!(
        db.list_job_attempts(1000).unwrap().len(),
        count,
        "inactive failed downloads must not be retried"
    );
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::DownloadFailed
    );
}

#[test]
fn technical_runs_pool_into_one_sample_and_conflicting_metadata_cannot_publish() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    let previous = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let mut duplicate = db.get_run("RUN2").unwrap().unwrap();
    duplicate.sample_accession = db.get_run("RUN1").unwrap().unwrap().sample_accession;
    db.upsert_run(&duplicate).unwrap();
    let previous_counts = read_matrix_tsv(
        &previous
            .paths(&paths)
            .unwrap()
            .de_project_dir(db.project_id())
            .join("gene_counts.tsv"),
    )
    .unwrap();
    success(root, &["normalize"]);
    let pooled = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let counts = read_matrix_tsv(
        &pooled
            .paths(&paths)
            .unwrap()
            .de_project_dir(db.project_id())
            .join("gene_counts.tsv"),
    )
    .unwrap();
    assert_eq!(
        counts.col_ids,
        ["BIO0", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6", "BIO7"]
    );
    for (row, old) in counts.values.iter().zip(&previous_counts.values) {
        assert_eq!(
            row[0],
            old[0] + old[1],
            "fractional contributions must be pooled before biological-sample fitting"
        );
    }
    let descriptor: serde_json::Value = serde_json::from_slice(
        &std::fs::read(
            pooled
                .paths(&paths)
                .unwrap()
                .de_project_dir(db.project_id())
                .join("rnaa_deseq2_request.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(
        descriptor["request"]["parameters"]["sample_plan"]["samples"][0]["runs"],
        serde_json::json!(["RUN1", "RUN2"])
    );
    success(root, &["normalize"]);
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        pooled.generation_id
    );
    success(
        root,
        &[
            "corr",
            "--adjust",
            "residualize",
            "--model",
            "~ condition",
            "--out",
            "edges:topk=5",
        ],
    );
    let sheet = std::fs::read_to_string(paths.samplesheet_path()).unwrap();
    let biological_matrix = pooled
        .paths(&paths)
        .unwrap()
        .de_project_dir(db.project_id())
        .join("vst.tsv");
    let custom_biological_matrix = root.join("custom_biological_samples.tsv");
    std::fs::copy(&biological_matrix, &custom_biological_matrix).unwrap();
    success(
        root,
        &[
            "corr",
            "--matrix",
            custom_biological_matrix.to_str().unwrap(),
            "--adjust",
            "residualize",
            "--model",
            "~ condition",
        ],
    );
    let mut run_matrix = previous_counts.clone();
    run_matrix.col_ids = (1..=8).map(|i| format!("RUN{i}")).collect();
    let run_matrix_path = root.join("custom_runs.tsv");
    rnaa_formats::matrix::write_matrix_tsv(&run_matrix_path, &run_matrix).unwrap();
    let duplicate_columns = command(
        root,
        &[
            "corr",
            "--matrix",
            run_matrix_path.to_str().unwrap(),
            "--adjust",
            "none",
        ],
    );
    assert!(!duplicate_columns.status.success());
    assert!(
        String::from_utf8_lossy(&duplicate_columns.stderr)
            .contains("multiple columns for biological sample")
    );
    std::fs::write(
        paths.samplesheet_path(),
        sheet.replace("RUN2\tA", "RUN2\tB"),
    )
    .unwrap();
    let rejected = command(root, &["normalize"]);
    assert!(!rejected.status.success());
    assert!(
        String::from_utf8_lossy(&rejected.stderr).contains("conflicting or missing 'condition'")
    );
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        pooled.generation_id
    );
    std::fs::write(paths.samplesheet_path(), sheet).unwrap();
    duplicate
        .metadata
        .insert("biological_sample_id".into(), "independent RUN2".into());
    db.upsert_run(&duplicate).unwrap();
    success(root, &["normalize"]);
    db.add_input("condition=A", InputType::FilterInclude)
        .unwrap();
    let before = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    assert!(!command(root, &["normalize"]).status.success());
    let status = command(root, &["status"]);
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );
    assert!(String::from_utf8_lossy(&status.stdout).contains("cohort_ready\tno"));
    let json_status = command(root, &["status", "--json"]);
    assert!(
        json_status.status.success(),
        "{}",
        String::from_utf8_lossy(&json_status.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&json_status.stdout).unwrap();
    assert_eq!(report["cohort"]["ready"], false);
    assert!(report["cohort"]["reason"].as_str().is_some());
    assert!(report["progress"].is_null());
    assert!(report["active_runs"].is_null());
    assert!(report["qc_failed"].is_null());
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        before.generation_id
    );
}

#[test]
fn explicit_current_matrix_uses_its_proven_contributors_despite_inactive_same_sample_history() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let mut run = db.get_run("RUN2").unwrap().unwrap();
    run.sample_accession = db.get_run("RUN1").unwrap().unwrap().sample_accession;
    db.upsert_run(&run).unwrap();
    let selected = db
        .list_runs()
        .unwrap()
        .into_iter()
        .filter(|run| run.run_accession != "RUN2")
        .collect::<Vec<_>>();
    db.publish_resolved_cohort(&selected, &db.selector_contract().unwrap(), None)
        .unwrap();
    success(root, &["normalize"]);
    let matrix = paths
        .with_generation_heads(&db)
        .unwrap()
        .de_project_dir(db.project_id())
        .join("vst.tsv");
    let samples = read_matrix_tsv(&matrix).unwrap().col_ids;
    assert_eq!(samples.len(), 7);
    assert!(samples.contains(&"BIO0".into()));
    success(
        root,
        &[
            "corr",
            "--matrix",
            matrix.to_str().unwrap(),
            "--adjust",
            "residualize",
            "--model",
            "~ condition",
            "--out",
            "edges:topk=5",
        ],
    );
    assert_eq!(
        db.get_run("RUN2").unwrap().unwrap().state,
        RunState::QuantDone
    );
}
