use rnaa_core::traits::ReferenceManager;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;

use flate2::Compression;
use flate2::write::GzEncoder;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::model::{ArtifactKind, ArtifactRecord, LibraryLayout, RunRecord};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::{ProjectStageState, RunState};
use rnaa_core::util::{compute_sha256, file_size, now_rfc3339};
use tempfile::TempDir;

fn rnaa_bin() -> &'static str {
    env!("CARGO_BIN_EXE_rnaa")
}

fn run_cmd(args: &[&str], extra_path: Option<&Path>) -> Output {
    let mut cmd = Command::new(rnaa_bin());
    cmd.args(args);
    if let Some(path) = extra_path {
        let current_path = std::env::var("PATH").unwrap_or_default();
        cmd.env("PATH", format!("{}:{}", path.display(), current_path));
    }
    cmd.output().expect("failed to execute rnaa")
}

fn run_ok(args: &[&str], extra_path: Option<&Path>) {
    let output = run_cmd(args, extra_path);
    assert!(
        output.status.success(),
        "command failed: rnaa {}\nstdout:\n{}\nstderr:\n{}",
        args.join(" "),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn write_gzip(path: &Path, contents: &[u8]) {
    let file = fs::File::create(path).expect("failed to create gzip fixture");
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder
        .write_all(contents)
        .expect("failed to write gzip fixture");
    encoder.finish().expect("failed to finish gzip fixture");
}

fn write_minimal_inputs(root: &Path) -> (PathBuf, PathBuf) {
    let cdna = root.join("test.cdna.fa.gz");
    let gtf = root.join("test.gtf.gz");
    write_gzip(
        &cdna,
        br#">tx1
ACGTACGTACGT
>tx2
TTTTAACCGG
"#,
    );
    write_gzip(&gtf, br#"chr1	src	transcript	1	100	.	+	.	gene_id "G1"; gene_name "Gene1"; gene_biotype "protein_coding"; transcript_id "tx1";
chr1	src	transcript	101	200	.	+	.	gene_id "G2"; gene_name "Gene2"; gene_biotype "protein_coding"; transcript_id "tx2";
"#);
    (cdna, gtf)
}

fn setup_custom_reference(root: &Path) -> (ProjectPaths, Database, String) {
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).expect("failed to open database");
    let (cdna, gtf) = write_minimal_inputs(root);

    let mut config = ProjectConfig::load(&paths.config_path).expect("failed to load config");
    config.refs.custom_cdna = cdna.display().to_string();
    config.refs.custom_gtf = gtf.display().to_string();
    config
        .save(&paths.config_path)
        .expect("failed to save project config");
    db.set_project_config(&config)
        .expect("failed to persist project config");

    let cdna_sha = compute_sha256(&cdna).expect("failed to hash cdna");
    let gtf_sha = compute_sha256(&gtf).expect("failed to hash gtf");
    let reference_id = format!("custom_{}_{}", &cdna_sha[..12], &gtf_sha[..12]);
    let reference_dir = paths.reference_dir(&reference_id);
    fs::create_dir_all(&reference_dir).expect("failed to create reference dir");
    fs::write(reference_dir.join("kallisto.index"), b"index").expect("failed to write index");
    fs::write(
        reference_dir.join("tx2gene.tsv"),
        b"transcript_id\tgene_id\ntx1\tG1\ntx2\tG2\n",
    )
    .expect("failed to write tx2gene");
    fs::write(
        reference_dir.join("gene_annotation.tsv"),
        b"gene_id\tgene_name\tgene_biotype\nG1\tGene1\tprotein_coding\nG2\tGene2\tprotein_coding\n",
    )
    .expect("failed to write annotation");

    let project_id = db.project().expect("failed to load project").project_id;
    (paths, db, project_id)
}

fn record_raw_fastq(project_id: &str, run_accession: &str, path: &Path) -> ArtifactRecord {
    ArtifactRecord {
        project_id: project_id.to_string(),
        run_accession: Some(run_accession.to_string()),
        kind: ArtifactKind::FastqSingle,
        path: path.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(path).expect("failed to hash fastq"),
        bytes: file_size(path).expect("failed to stat fastq"),
        created_at: now_rfc3339(),
    }
}

fn record_quant_artifacts(
    paths: &ProjectPaths,
    db: &Database,
    project_id: &str,
    run_accession: &str,
) {
    let quant_dir = paths.quant_run_dir(run_accession, "kallistors");
    fs::create_dir_all(&quant_dir).expect("failed to create quant dir");
    let abundance_h5 = quant_dir.join("abundance.h5");
    let run_info_json = quant_dir.join("run_info.json");
    fs::write(
        quant_dir.join("abundance.tsv"),
        "target_id\tlength\teff_length\test_counts\ttpm\ntx1\t12\t10\t100\t400000\ntx2\t10\t8\t120\t600000\n",
    )
    .expect("failed to write abundance.tsv");
    fs::write(&abundance_h5, b"mock-h5").expect("failed to write abundance.h5");
    fs::write(&run_info_json, b"{\"run\":\"ok\"}").expect("failed to write run_info.json");
    db.record_artifact(&ArtifactRecord {
        project_id: project_id.to_string(),
        run_accession: Some(run_accession.to_string()),
        kind: ArtifactKind::QuantAbundanceH5,
        path: abundance_h5.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(&abundance_h5).expect("failed to hash abundance.h5"),
        bytes: file_size(&abundance_h5).expect("failed to stat abundance.h5"),
        created_at: now_rfc3339(),
    })
    .expect("failed to record abundance artifact");
    db.record_artifact(&ArtifactRecord {
        project_id: project_id.to_string(),
        run_accession: Some(run_accession.to_string()),
        kind: ArtifactKind::QuantRunInfo,
        path: run_info_json.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(&run_info_json).expect("failed to hash run_info.json"),
        bytes: file_size(&run_info_json).expect("failed to stat run_info.json"),
        created_at: now_rfc3339(),
    })
    .expect("failed to record run info artifact");
    let config = ProjectConfig::load(&paths.config_path).unwrap();
    let reference = rnaa_stages::refs::EnsemblReferenceManager::new()
        .unwrap()
        .ensure_reference(paths, db, &config)
        .unwrap();
    rnaa_stages::quant::declare_external_quant(
        run_accession,
        &quant_dir,
        &reference,
        "external restart statistical fixture",
    )
    .unwrap();
}

#[test]
fn quant_restart_recovers_preprocess_and_completes() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        None,
    );

    let (paths, db, project_id) = setup_custom_reference(root);

    let raw1 = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    let raw2 = paths.raw_run_dir("RUN2").join("RUN2.fastq");
    fs::create_dir_all(raw1.parent().expect("missing parent")).expect("failed to create raw dir");
    fs::create_dir_all(raw2.parent().expect("missing parent")).expect("failed to create raw dir");
    let sequence = "ACGTTCGATCGGATCGATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCTAGCTGACTGATCGATCGTAGCTAGCTGATCGTAGCTAGCTAGCATCGA";
    let read = format!("@read\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()));
    fs::write(&raw1, &read).expect("failed to write fastq");
    fs::write(&raw2, &read).expect("failed to write fastq");

    db.record_artifact(&record_raw_fastq(&project_id, "RUN1", &raw1))
        .expect("failed to record raw fastq");
    db.record_artifact(&record_raw_fastq(&project_id, "RUN2", &raw2))
        .expect("failed to record raw fastq");

    let runs = vec![
        RunRecord {
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
            state: RunState::PreprocessRunning,
            last_error: None,
            updated_at: now_rfc3339(),
        },
        RunRecord {
            project_id: project_id.clone(),
            study_accession: "SRPTEST".to_string(),
            source_study_accession: None,
            run_accession: "RUN2".to_string(),
            sample_accession: Some("SAMP_RUN2".to_string()),
            experiment_accession: Some("EXP_RUN2".to_string()),
            library_layout: LibraryLayout::Single,
            instrument_platform: Some("ILLUMINA".to_string()),
            instrument_model: Some("TEST".to_string()),
            remote_files: Vec::new(),
            metadata: BTreeMap::new(),
            state: RunState::Verified,
            last_error: None,
            updated_at: now_rfc3339(),
        },
    ];
    db.upsert_runs(&runs).expect("failed to upsert runs");

    // A separate downloader may own another run while this CLI invocation
    // recovers an actual orphan and quantifies the ready runs.
    let mut downloading = runs[1].clone();
    downloading.run_accession = "LIVE_DOWNLOAD".into();
    downloading.state = RunState::Resolved;
    db.upsert_run(&downloading).unwrap();
    let mut foreign_download = rnaa_core::jobs::JobLease::try_acquire(
        &db,
        "download",
        Some("LIVE_DOWNLOAD"),
        &[RunState::Resolved],
    )
    .unwrap()
    .unwrap();

    let output = run_cmd(
        &[
            "quant",
            "--root",
            root.to_str().expect("invalid root path"),
            "--preprocess",
            "--preprocess-bypass",
            "--workers",
            "1",
        ],
        None,
    );
    assert!(
        output.status.success(),
        "quant command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let reopened = Database::open(root).expect("failed to reopen db");
    let runs = reopened.list_runs().expect("failed to list runs");
    assert!(
        runs.iter()
            .filter(|run| run.run_accession != "LIVE_DOWNLOAD")
            .all(|run| run.state == RunState::QuantDone),
        "unexpected run states: {:?}\nstdout:\n{}\nstderr:\n{}",
        runs.iter()
            .map(|run| (
                run.run_accession.clone(),
                run.state.as_str().to_string(),
                run.last_error.clone()
            ))
            .collect::<Vec<_>>(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        paths
            .with_generation_heads(&db)
            .unwrap()
            .quant_run_dir("RUN1", "kallistors")
            .join("abundance.h5")
            .exists()
    );
    assert!(
        reopened
            .list_events(100)
            .expect("failed to list events")
            .iter()
            .any(|event| event.stage == "recovery" && event.message == "interrupted work recovered")
    );
    assert_eq!(
        reopened.get_run("LIVE_DOWNLOAD").unwrap().unwrap().state,
        RunState::Downloading
    );
    assert_eq!(
        reopened
            .list_job_attempts(100)
            .unwrap()
            .iter()
            .find(|job| job.attempt_id == foreign_download.attempt_id())
            .unwrap()
            .state,
        "RUNNING"
    );
    foreign_download
        .finish(rnaa_core::jobs::JobOutcome::Failed, Some("test complete"))
        .unwrap();
}

#[test]
fn normalize_restart_recovers_progress_then_rejects_confounded_design() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        None,
    );

    let (paths, db, project_id) = setup_custom_reference(root);
    let runs = ["RUN1", "RUN2"]
        .into_iter()
        .map(|run| RunRecord {
            project_id: project_id.clone(),
            study_accession: "SRPTEST".to_string(),
            source_study_accession: None,
            run_accession: run.to_string(),
            sample_accession: Some(format!("SAMP_{run}")),
            experiment_accession: Some(format!("EXP_{run}")),
            library_layout: LibraryLayout::Single,
            instrument_platform: Some("ILLUMINA".to_string()),
            instrument_model: Some("TEST".to_string()),
            remote_files: Vec::new(),
            metadata: BTreeMap::from([
                (
                    "condition".to_string(),
                    if run == "RUN1" { "A" } else { "B" }.to_string(),
                ),
                (
                    "batch".to_string(),
                    if run == "RUN1" { "B1" } else { "B2" }.to_string(),
                ),
            ]),
            state: RunState::QuantDone,
            last_error: None,
            updated_at: now_rfc3339(),
        })
        .collect::<Vec<_>>();
    db.upsert_runs(&runs).expect("failed to insert runs");
    record_quant_artifacts(&paths, &db, &project_id, "RUN1");
    record_quant_artifacts(&paths, &db, &project_id, "RUN2");
    db.set_project_stage_state("normalize", ProjectStageState::Running, None)
        .expect("failed to set project stage");

    let samplesheet = "project_id\tstudy_accession\trun_accession\tcondition\tbatch\n\
                      test\tSRPTEST\tRUN1\tA\tB1\n\
                      test\tSRPTEST\tRUN2\tB\tB2\n";
    fs::write(paths.samplesheet_path(), samplesheet).expect("failed to write samplesheet");

    let output = run_cmd(
        &[
            "normalize",
            "--root",
            root.to_str().expect("invalid root path"),
            "--design",
            "~ batch + condition",
        ],
        None,
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("rank deficient"));

    let reopened = Database::open(root).expect("failed to reopen db");
    let normalize_stage = reopened
        .get_project_stage_state("normalize")
        .expect("failed to load normalize stage")
        .expect("normalize stage missing");
    assert_eq!(normalize_stage.state, ProjectStageState::Failed);
    assert!(
        normalize_stage
            .last_error
            .as_deref()
            .unwrap_or_default()
            .contains("rank deficient")
    );
    assert!(
        reopened
            .list_runs()
            .expect("failed to list runs")
            .iter()
            .all(|run| run.state == RunState::QuantDone)
    );
    assert!(
        reopened
            .list_events(200)
            .expect("failed to list events")
            .iter()
            .any(|event| event.stage == "recovery" && event.message == "interrupted work recovered")
    );
}

#[test]
fn shared_snapshot_survives_requantification_and_preprocess_respects_budget() {
    let temp = TempDir::new().unwrap();
    let root = temp.path();
    let root_str = root.to_str().unwrap();
    run_ok(&["init", "--root", root_str], None);
    let (paths, db, project_id) = setup_custom_reference(root);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.storage.shared_root = root.join("shared").display().to_string();
    config.save(&paths.config_path).unwrap();
    let raw = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    fs::create_dir_all(raw.parent().unwrap()).unwrap();
    let sequence = "ACGTTCGATCGGATCGATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCTAGCTGACTGATCGATCGTAGCTAGCTGATCGTAGCTAGCTAGCATCGA";
    let read = format!("@r1\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()));
    fs::write(&raw, &read).unwrap();
    db.record_artifact(&record_raw_fastq(&project_id, "RUN1", &raw))
        .unwrap();
    db.upsert_run(&RunRecord {
        project_id: project_id.clone(),
        study_accession: "SRPTEST".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: Some("SAMP1".into()),
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: Vec::new(),
        metadata: BTreeMap::new(),
        state: RunState::Verified,
        last_error: None,
        updated_at: now_rfc3339(),
    })
    .unwrap();
    run_ok(&["quant", "--root", root_str, "--worker-budget", "1"], None);
    let info: serde_json::Value = serde_json::from_slice(
        &fs::read(
            paths
                .with_generation_heads(&db)
                .unwrap()
                .quant_run_dir("RUN1", "kallistors")
                .join("run_info.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(info["n_processed"], 1);
    let original = db
        .list_artifacts_for_run(Some("RUN1"))
        .unwrap()
        .into_iter()
        .find(|a| a.kind == ArtifactKind::QuantAbundanceH5)
        .unwrap();
    let shared = PathBuf::from(original.shared_path.unwrap());
    let old_checksum = compute_sha256(&shared).unwrap();
    let extra = paths.raw_run_dir("RUN1").join("extra.fastq");
    fs::write(&extra, read.repeat(2)).unwrap();
    db.record_artifact(&record_raw_fastq(&project_id, "RUN1", &extra))
        .unwrap();
    run_ok(&["quant", "--root", root_str, "--worker-budget", "1"], None);
    let info: serde_json::Value = serde_json::from_slice(
        &fs::read(
            paths
                .with_generation_heads(&db)
                .unwrap()
                .quant_run_dir("RUN1", "kallistors")
                .join("run_info.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(
        info["n_processed"], 3,
        "newly registered FASTQ must invalidate the completed request"
    );
    let manifest_path = paths
        .with_generation_heads(&db)
        .unwrap()
        .quant_run_dir("RUN1", "kallistors")
        .join("rnaa_quant_manifest.json");
    let saved_generation = compute_sha256(&manifest_path).unwrap();
    run_ok(
        &[
            "quant",
            "--root",
            root_str,
            "--worker-budget",
            "1",
            "--threads",
            "2",
        ],
        None,
    );
    assert_eq!(
        compute_sha256(&manifest_path).unwrap(),
        saved_generation,
        "CPU-only changes must reuse scientific outputs"
    );
    fs::write(&raw, read.repeat(3)).unwrap();
    db.record_artifact(&record_raw_fastq(&project_id, "RUN1", &raw))
        .unwrap();
    run_ok(&["quant", "--root", root_str, "--worker-budget", "1"], None);
    assert_eq!(compute_sha256(&shared).unwrap(), old_checksum);
    assert_ne!(
        compute_sha256(
            &paths
                .with_generation_heads(&db)
                .unwrap()
                .quant_run_dir("RUN1", "kallistors")
                .join("abundance.h5")
        )
        .unwrap(),
        old_checksum
    );
    run_ok(
        &[
            "quant",
            "--root",
            root_str,
            "--worker-budget",
            "1",
            "--preprocess",
            "--preprocess-threads",
            "16",
        ],
        None,
    );
    let report: serde_json::Value = serde_json::from_slice(
        &fs::read(
            paths
                .with_generation_heads(&db)
                .unwrap()
                .preprocess_run_dir("RUN1", "kallistors")
                .join("fasterp.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(report["threads"], 1);
    let qc = db.get_preprocess_qc("RUN1").unwrap().unwrap();
    assert_eq!(qc.total_reads_before, 5);
    assert_eq!(compute_sha256(&shared).unwrap(), old_checksum);
    // A changed QC policy must revisit eligibility even when trimmed outputs
    // can be reused. Five good reads and one rejected read pass defaults.
    let bad = format!("@bad\n{sequence}\n+\n{}\n", "!".repeat(sequence.len()));
    fs::write(&raw, format!("{}{bad}", read.repeat(3))).unwrap();
    db.record_artifact(&record_raw_fastq(&project_id, "RUN1", &raw))
        .unwrap();
    run_ok(&["quant", "--root", root_str, "--worker-budget", "1"], None);
    let mut current = ProjectConfig::load(&paths.config_path).unwrap();
    current.quant.qc_gate.min_pass_rate = Some(0.95);
    current.save(&paths.config_path).unwrap();
    let result = run_cmd(&["quant", "--root", root_str, "--worker-budget", "1"], None);
    assert!(!result.status.success());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantDone
    );
    assert!(!db.get_preprocess_qc("RUN1").unwrap().unwrap().gate_passed);
}

#[test]
fn missing_paired_mate_is_a_failed_command() {
    let temp = TempDir::new().unwrap();
    let root = temp.path().to_str().unwrap();
    run_ok(&["init", "--root", root], None);
    let (paths, db, project_id) = setup_custom_reference(temp.path());
    let raw = paths.raw_run_dir("RUN1").join("RUN1_1.fastq");
    fs::create_dir_all(raw.parent().unwrap()).unwrap();
    fs::write(&raw, "@read/1\nACGT\n+\nIIII\n").unwrap();
    let mut artifact = record_raw_fastq(&project_id, "RUN1", &raw);
    artifact.kind = ArtifactKind::FastqR1;
    db.record_artifact(&artifact).unwrap();
    db.upsert_run(&RunRecord {
        project_id,
        study_accession: "study".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Paired,
        instrument_platform: None,
        instrument_model: None,
        remote_files: vec![],
        metadata: BTreeMap::new(),
        state: RunState::Verified,
        last_error: None,
        updated_at: now_rfc3339(),
    })
    .unwrap();
    let result = run_cmd(&["quant", "--root", root, "--worker-budget", "1"], None);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stdout).contains("quant_failed\t1"));
    let run = db.get_run("RUN1").unwrap().unwrap();
    assert_eq!(run.state, RunState::QuantFailed);
    assert!(run.last_error.unwrap().contains("declared PAIRED layout"));
}
