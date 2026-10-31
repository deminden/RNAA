use super::*;
use rnaa_core::LibraryLayout;
use std::{collections::BTreeMap, path::PathBuf};

fn run(id: &str, state: RunState) -> RunRecord {
    RunRecord {
        project_id: "project".into(),
        study_accession: "study".into(),
        source_study_accession: None,
        run_accession: id.into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: vec![],
        metadata: BTreeMap::new(),
        state,
        last_error: None,
        updated_at: "2026-10-30T00:00:00Z".into(),
    }
}
fn event(run: Option<&str>, message: &str, seconds: usize) -> EventRecord {
    EventRecord {
        id: 0,
        project_id: "project".into(),
        run_accession: run.map(str::to_string),
        stage: "quant".into(),
        message: message.into(),
        ts: format!("2026-10-30T00:00:{seconds:02}Z"),
        context: serde_json::json!({}),
    }
}
fn stage(name: &str, state: ProjectStageState) -> ProjectStageRecord {
    ProjectStageRecord {
        project_id: "project".into(),
        stage: name.into(),
        state,
        last_error: None,
        updated_at: "2026-10-30T00:00:00Z".into(),
    }
}
#[test]
fn empty_cohorts_and_legacy_run_labels_do_not_claim_analysis_completion() {
    let config = ProjectConfig::default();
    let stages = [
        stage("deseq2", ProjectStageState::Done),
        stage("corr", ProjectStageState::Done),
    ];
    let empty = compute_progress_snapshot(&[], &[], &[], &stages, &config);
    assert!(!empty.normalize_done && !empty.corr_stage_done);
    assert_eq!(empty.active_stage, "waiting");
    for state in [RunState::QuantDone, RunState::CorrDone, RunState::Cleaned] {
        let snapshot = compute_progress_snapshot(&[run("RUN1", state)], &[], &[], &[], &config);
        assert_eq!(snapshot.quant_done, 1);
        assert!(!snapshot.normalize_done && !snapshot.corr_stage_done);
    }
    let snapshot = compute_progress_snapshot(
        &[run("RUN1", RunState::QuantDone)],
        &[],
        &[],
        &stages,
        &config,
    );
    assert!(snapshot.normalize_done && snapshot.corr_stage_done);
}
#[test]
fn excluded_and_unattributed_artifacts_and_ingest_events_cannot_inflate_progress() {
    let artifact = |run: Option<&str>, bytes| ArtifactRecord {
        project_id: "project".into(),
        run_accession: run.map(str::to_string),
        kind: ArtifactKind::FastqSingle,
        path: PathBuf::from("unused").display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "none".into(),
        checksum: String::new(),
        bytes,
        created_at: "2026-10-30T00:00:00Z".into(),
    };
    let artifacts = [
        artifact(Some("ACTIVE"), 10),
        artifact(Some("EXCLUDED"), 9000),
        artifact(None, 5000),
    ];
    let events = [
        event(Some("EXCLUDED"), "quantification started", 0),
        event(Some("EXCLUDED"), "quantification completed", 10),
        event(None, "quantification completed", 20),
    ];
    let snapshot = compute_progress_snapshot(
        &[run("ACTIVE", RunState::Verified)],
        &artifacts,
        &events,
        &[],
        &ProjectConfig::default(),
    );
    assert_eq!(snapshot.total_runs, 1);
    assert_eq!(snapshot.download_done, 1);
    assert_eq!(snapshot.downloaded_bytes, 10);
    assert_eq!(snapshot.quant_rate_runs_per_sec, 0.0);
}
#[test]
fn mixed_quant_completion_keeps_remaining_work_in_eta() {
    let mut config = ProjectConfig::default();
    config.quant.preprocess = false;
    let runs = [
        run("RUN1", RunState::QuantDone),
        run("RUN2", RunState::QuantDone),
        run("RUN3", RunState::Verified),
        run("RUN4", RunState::Verified),
    ];
    let events = [
        event(Some("RUN1"), "quantification started", 0),
        event(Some("RUN1"), "quantification completed", 10),
        event(Some("RUN2"), "quantification started", 10),
        event(Some("RUN2"), "quantification completed", 20),
    ];
    let snapshot = compute_progress_snapshot(&runs, &[], &events, &[], &config);
    assert_eq!(snapshot.quant_done, 2);
    assert!((snapshot.quant_rate_runs_per_sec - 0.1).abs() < 1e-12);
    assert_eq!(snapshot.processing_eta, Some(Duration::from_secs(20)));
}

#[test]
fn cohort_attempt_state_distinguishes_running_failed_and_waiting_analysis() {
    let runs = [run("RUN1", RunState::QuantDone)];
    let config = ProjectConfig::default();
    for (name, label) in [
        ("normalize", "normalization"),
        ("deseq2", "differential expression"),
        ("corr", "correlation"),
    ] {
        let running = compute_progress_snapshot(
            &runs,
            &[],
            &[],
            &[stage(name, ProjectStageState::Running)],
            &config,
        );
        assert_eq!(running.active_stage, label);
        assert!(running.active_detail.starts_with("running "));
        let failed = compute_progress_snapshot(
            &runs,
            &[],
            &[],
            &[stage(name, ProjectStageState::Failed)],
            &config,
        );
        assert_eq!(failed.active_stage, "waiting");
        assert!(failed.active_detail.contains("failed cohort analysis"));
    }
    let waiting = compute_progress_snapshot(&runs, &[], &[], &[], &config);
    assert_eq!(waiting.active_stage, "waiting");
    assert!(
        waiting
            .active_detail
            .contains("waiting for cohort analysis")
    );
}
