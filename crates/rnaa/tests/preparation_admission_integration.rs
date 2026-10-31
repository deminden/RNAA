//! Decoded admission rejects fresh work and preserves healthy historical results.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/quant.rs"]
mod quant_support;
use commands::{command, success};
use quant_support::{record_raw, setup};
use rnaa_core::{Database, ProjectConfig, ProjectPaths, RunState, generations::OutputGroup};
use std::{fs, io::Write};

fn add_lanes(paths: &ProjectPaths, db: &Database) {
    let read = fs::read_to_string(paths.raw_run_dir("RUN1").join("reads.fastq")).unwrap();
    for lane in ["lane_a.fastq.gz", "lane_b.fastq.gz"] {
        let path = paths.raw_run_dir("RUN1").join(lane);
        let mut encoder = flate2::write::GzEncoder::new(
            fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        for _ in 0..3000 {
            encoder.write_all(read.as_bytes()).unwrap();
        }
        encoder.finish().unwrap();
        assert!(fs::metadata(&path).unwrap().len() < 1024 * 1024);
        record_raw(db, "RUN1", &path);
    }
}
fn save(paths: &ProjectPaths, db: &Database, config: &ProjectConfig) {
    config.save(&paths.config_path).unwrap();
    db.set_project_config(config).unwrap();
}
const QUANT: &[&str] = &["quant", "--worker-budget", "2", "--threads", "1"];
const PREPROCESS: &[&str] = &[
    "quant",
    "--worker-budget",
    "2",
    "--threads",
    "1",
    "--preprocess",
    "--preprocess-threads",
    "2",
];

#[test]
fn quantification_admits_total_decoded_lanes_and_reuses_without_reapplying_work_limits() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path(), false);
    add_lanes(&paths, &db);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.fastq_max_decoded_mb = 1;
    save(&paths, &db, &config);
    let result = command(root.path(), QUANT);
    assert!(!result.status.success());
    let failed = db.get_run("RUN1").unwrap().unwrap();
    assert_eq!(failed.state, RunState::QuantFailed);
    assert!(
        failed
            .last_error
            .unwrap()
            .contains("decoded FASTQ group exceeds")
    );
    assert!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .is_none()
    );
    let good = ["RUN2", "RUN3", "RUN4"].map(|run| {
        db.generation_head(OutputGroup::Quant, run)
            .unwrap()
            .unwrap()
    });
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.fastq_max_decoded_mb = 0;
    save(&paths, &db, &config);
    success(root.path(), QUANT);
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantDone
    );
    for head in good {
        assert_eq!(
            db.generation_head(OutputGroup::Quant, &head.subject)
                .unwrap()
                .unwrap(),
            head
        );
    }
    let current = paths.with_generation_heads(&db).unwrap();
    let manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(
            current
                .quant_run_dir("RUN1", "kallistors")
                .join("rnaa_quant_manifest.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(manifest["parameters"]["fastq_execution"]["lanes"], 3);
    assert!(
        manifest["parameters"]["fastq_execution"]["decoded_bytes"]
            .as_u64()
            .unwrap()
            > 1024 * 1024
    );
    let heads = db.generation_heads().unwrap();
    let attempts = db.list_job_attempts(100).unwrap().len();
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.fastq_max_decoded_mb = 1;
    save(&paths, &db, &config);
    success(root.path(), QUANT);
    assert_eq!(db.generation_heads().unwrap(), heads);
    assert_eq!(db.list_job_attempts(100).unwrap().len(), attempts);
}

#[test]
fn preprocessing_decoded_cap_is_aggregate_and_successful_outputs_survive_resource_changes() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path(), false);
    add_lanes(&paths, &db);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.preprocess_max_input_mb = 1;
    save(&paths, &db, &config);
    let result = command(root.path(), PREPROCESS);
    assert!(!result.status.success());
    let failed = db.get_run("RUN1").unwrap().unwrap();
    assert_eq!(failed.state, RunState::PreprocessFailed);
    assert!(
        failed
            .last_error
            .unwrap()
            .contains("preprocess_max_input_mb=1")
    );
    assert!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .is_none()
    );
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.preprocess_max_input_mb = 4;
    save(&paths, &db, &config);
    success(root.path(), PREPROCESS);
    let current = paths.with_generation_heads(&db).unwrap();
    let report: serde_json::Value = serde_json::from_slice(
        &fs::read(
            current
                .preprocess_run_dir("RUN1", "kallistors")
                .join("fasterp.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(report["execution"]["batch_window"]["active"], 0);
    assert!(
        report["execution"]["batch_window"]["peak"]
            .as_u64()
            .unwrap()
            <= 3
    );
    let manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(
            current
                .preprocess_run_dir("RUN1", "kallistors")
                .join("rnaa_preprocess_manifest.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(manifest["parameters"]["algorithm_version"], 3);
    let heads = db.generation_heads().unwrap();
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.preprocess_max_input_mb = 1;
    save(&paths, &db, &config);
    success(root.path(), PREPROCESS);
    assert_eq!(db.generation_heads().unwrap(), heads);
}
