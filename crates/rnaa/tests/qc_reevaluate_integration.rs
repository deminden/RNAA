use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use std::process::{Command, Output};

use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::model::{
    ArtifactKind, ArtifactRecord, LibraryLayout, PreprocessQcRecord, RunRecord,
};
use rnaa_core::state::{ProjectStageState, RunState};
use rnaa_core::util::{compute_sha256, file_size, now_rfc3339};
use tempfile::TempDir;

fn rnaa_bin() -> &'static str {
    env!("CARGO_BIN_EXE_rnaa")
}

fn run_cmd(args: &[&str]) -> Output {
    Command::new(rnaa_bin())
        .args(args)
        .output()
        .expect("failed to execute rnaa")
}

fn run_ok(args: &[&str]) -> Output {
    let output = run_cmd(args);
    assert!(
        output.status.success(),
        "command failed: rnaa {}\nstdout:\n{}\nstderr:\n{}",
        args.join(" "),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn record_artifact(
    db: &Database,
    project_id: &str,
    run_accession: &str,
    kind: ArtifactKind,
    path: &Path,
) {
    db.record_artifact(&ArtifactRecord {
        project_id: project_id.to_string(),
        run_accession: Some(run_accession.to_string()),
        kind,
        path: path.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(path).expect("failed to hash artifact"),
        bytes: file_size(path).expect("failed to stat artifact"),
        created_at: now_rfc3339(),
    })
    .expect("failed to record artifact");
}

#[test]
fn qc_reevaluate_invalidates_downstream_when_gate_outcome_changes() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();

    run_ok(&["init", "--root", root.to_str().expect("invalid root path")]);

    let db = Database::open(root).expect("failed to open database");
    let project_id = db.project().expect("failed to load project").project_id;

    let run = RunRecord {
        project_id: project_id.clone(),
        study_accession: "SRPTEST".to_string(),
        source_study_accession: None,
        run_accession: "RUN1".to_string(),
        sample_accession: Some("SAMP_RUN1".to_string()),
        experiment_accession: Some("EXP_RUN1".to_string()),
        library_layout: LibraryLayout::Single,
        instrument_platform: Some("ILLUMINA".to_string()),
        instrument_model: Some("TEST".to_string()),
        remote_files: Vec::new(),
        metadata: BTreeMap::new(),
        state: RunState::CorrDone,
        last_error: None,
        updated_at: now_rfc3339(),
    };
    db.upsert_run(&run).expect("failed to upsert run");

    let report_dir = root
        .join("quant")
        .join("RUN1")
        .join("r-kallisto")
        .join("fasterp");
    fs::create_dir_all(&report_dir).expect("failed to create report dir");
    let report_path = report_dir.join("fasterp.json");
    fs::write(
        &report_path,
        r#"{
          "mode":"single",
          "threads":16,
          "summary":{"before_filtering":{"total_reads":100},"after_filtering":{"total_reads":70}},
          "filtering_result":{
            "passed_filter_reads":70,
            "low_quality_reads":20,
            "low_complexity_reads":0,
            "too_many_N_reads":0,
            "too_short_reads":10,
            "too_long_reads":0
          },
          "duplicated_reads":5,
          "duplication":{"rate":0.05},
          "adapter_cutting":{"adapter_trimmed_reads":40,"adapter_trimmed_bases":120}
        }"#,
    )
    .expect("failed to write preprocess report");
    record_artifact(
        &db,
        &project_id,
        "RUN1",
        ArtifactKind::PreprocessReport,
        &report_path,
    );

    let trimmed_fastq = report_dir.join("trimmed.fastq");
    fs::write(&trimmed_fastq, b"@r1\nACGT\n+\n!!!!\n").expect("failed to write trimmed fastq");
    record_artifact(
        &db,
        &project_id,
        "RUN1",
        ArtifactKind::FastqTrimmedSingle,
        &trimmed_fastq,
    );

    let quant_dir = root.join("quant").join("RUN1").join("r-kallisto");
    fs::create_dir_all(&quant_dir).expect("failed to create quant dir");
    let run_info = quant_dir.join("run_info.json");
    fs::write(&run_info, br#"{"run":"RUN1"}"#).expect("failed to write run info");
    record_artifact(
        &db,
        &project_id,
        "RUN1",
        ArtifactKind::QuantRunInfo,
        &run_info,
    );

    db.record_preprocess_qc(&PreprocessQcRecord {
        project_id: project_id.clone(),
        run_accession: "RUN1".to_string(),
        report_path: report_path.display().to_string(),
        mode: "single".to_string(),
        threads: 16,
        total_reads_before: 100,
        passed_reads: 70,
        failed_reads: 30,
        pass_rate: 0.70,
        failed_rate: 0.30,
        low_quality_reads: 20,
        low_complexity_reads: 0,
        too_many_n_reads: 0,
        too_short_reads: 10,
        duplicated_reads: 5,
        duplication_rate: Some(0.05),
        adapter_trimmed_reads: 40,
        adapter_trimmed_bases: 120,
        gate_passed: true,
        gate_reason: None,
        updated_at: now_rfc3339(),
    })
    .expect("failed to seed preprocess qc");

    db.set_project_stage_state("deseq2", ProjectStageState::Done, None)
        .expect("failed to seed deseq2 stage");
    db.set_project_stage_state("corr", ProjectStageState::Done, None)
        .expect("failed to seed corr stage");

    let mut config = ProjectConfig::load(&root.join("rnaa.toml")).expect("failed to load config");
    config.quant.qc_gate.min_pass_rate = Some(0.80);
    config
        .save(&root.join("rnaa.toml"))
        .expect("failed to save config");
    db.set_project_config(&config)
        .expect("failed to persist config");

    let output = run_ok(&[
        "qc",
        "reevaluate",
        "--root",
        root.to_str().expect("invalid root path"),
    ]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("qc_failed\t1"));
    assert!(stdout.contains("project_outputs_invalidated\tyes"));

    let run = db
        .list_runs()
        .expect("failed to list runs")
        .into_iter()
        .find(|item| item.run_accession == "RUN1")
        .expect("missing RUN1");
    assert_eq!(run.state, RunState::PreprocessFailed);
    assert!(
        run.last_error
            .as_deref()
            .is_some_and(|reason| reason.contains("pass_rate"))
    );

    let qc = db
        .get_preprocess_qc("RUN1")
        .expect("failed to load preprocess qc")
        .expect("missing preprocess qc");
    assert!(!qc.gate_passed);

    let deseq2_stage = db
        .get_project_stage_state("deseq2")
        .expect("failed to load deseq2 stage")
        .expect("missing deseq2 stage");
    let corr_stage = db
        .get_project_stage_state("corr")
        .expect("failed to load corr stage")
        .expect("missing corr stage");
    assert_eq!(deseq2_stage.state, ProjectStageState::Pending);
    assert_eq!(corr_stage.state, ProjectStageState::Pending);

    let status_output = run_ok(&[
        "status",
        "--root",
        root.to_str().expect("invalid root path"),
    ]);
    let status_stdout = String::from_utf8_lossy(&status_output.stdout);
    assert!(status_stdout.contains("runs_skipped_qc\t1"));
}
