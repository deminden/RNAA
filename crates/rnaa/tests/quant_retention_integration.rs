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
use rnaa_core::state::RunState;
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
"#,
    );
    write_gzip(&gtf, br#"chr1	src	transcript	1	100	.	+	.	gene_id "G1"; gene_name "Gene1"; gene_biotype "protein_coding"; transcript_id "tx1";
"#);
    (cdna, gtf)
}

fn setup_project(root: &Path) -> (ProjectPaths, Database, String) {
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
        b"transcript_id\tgene_id\ntx1\tG1\n",
    )
    .expect("failed to write tx2gene");
    fs::write(
        reference_dir.join("gene_annotation.tsv"),
        b"gene_id\tgene_name\tgene_biotype\nG1\tGene1\tprotein_coding\n",
    )
    .expect("failed to write annotation");

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
        state: RunState::Verified,
        last_error: None,
        updated_at: now_rfc3339(),
    };
    db.upsert_run(&run).expect("failed to upsert run");

    let raw_dir = paths.raw_run_dir("RUN1");
    fs::create_dir_all(&raw_dir).expect("failed to create raw dir");
    let raw_fastq = raw_dir.join("RUN1.fastq");
    fs::write(&raw_fastq, b"@r1\nACGT\n+\n!!!!\n").expect("failed to write fastq");
    db.record_artifact(&ArtifactRecord {
        project_id: project_id.clone(),
        run_accession: Some("RUN1".to_string()),
        kind: ArtifactKind::FastqSingle,
        path: raw_fastq.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(&raw_fastq).expect("failed to hash raw fastq"),
        bytes: file_size(&raw_fastq).expect("failed to stat raw fastq"),
        created_at: now_rfc3339(),
    })
    .expect("failed to record artifact");

    (paths, db, project_id)
}

fn seed_trimmed_fastq(paths: &ProjectPaths, db: &Database, project_id: &str) {
    use rnaa_core::{VerifiedFile, traits::Preprocessor};
    let raw = replace_raw(paths, db, &good_read());
    let run = db.get_run("RUN1").unwrap().unwrap();
    let config = ProjectConfig::load(&paths.config_path).unwrap();
    let inputs = [VerifiedFile {
        artifact_kind: ArtifactKind::FastqSingle,
        checksum_type: "sha256".into(),
        checksum: compute_sha256(&raw).unwrap(),
        bytes: file_size(&raw).unwrap(),
        integrity_source: "fixture".into(),
        path: raw,
    }];
    let artifacts = rnaa_stages::preprocess::FasterpInProcessPreprocessor
        .preprocess(&run, &inputs, paths, &config)
        .unwrap();
    for (kind, path) in [
        (
            ArtifactKind::FastqTrimmedSingle,
            artifacts.fastqs[0].path.clone(),
        ),
        (ArtifactKind::PreprocessReport, artifacts.report_json),
        (
            ArtifactKind::Manifest,
            artifacts.out_dir.join("rnaa_preprocess_manifest.json"),
        ),
    ] {
        db.record_artifact(&ArtifactRecord {
            project_id: project_id.into(),
            run_accession: Some("RUN1".into()),
            kind,
            path: path.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".into(),
            checksum: compute_sha256(&path).unwrap(),
            bytes: file_size(&path).unwrap(),
            created_at: now_rfc3339(),
        })
        .unwrap();
    }
    db.set_run_state("RUN1", RunState::PreprocessDone, None)
        .unwrap();
}

#[test]
fn quant_retention_trimmed_deletes_original_fastq() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        None,
    );
    let (paths, db, project_id) = setup_project(root);
    seed_trimmed_fastq(&paths, &db, &project_id);

    run_ok(
        &[
            "quant",
            "--root",
            root.to_str().expect("invalid root path"),
            "--preprocess",
            "--fastq-retention",
            "trimmed",
        ],
        None,
    );

    let raw_fastq = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    let trimmed_fastq = paths
        .preprocess_run_dir("RUN1", "kallistors")
        .join("trimmed.fastq");
    let artifacts = db
        .list_artifacts_for_run(Some("RUN1"))
        .expect("failed to list artifacts");
    assert!(
        !raw_fastq.exists(),
        "raw FASTQ should be removed from raw dir"
    );
    assert!(
        !paths.raw_run_dir("RUN1").exists(),
        "empty raw run dir should be pruned"
    );
    assert!(trimmed_fastq.exists(), "trimmed FASTQ should remain");
    assert!(
        artifacts
            .iter()
            .any(|item| item.kind == ArtifactKind::FastqTrimmedSingle)
    );
    assert!(
        !artifacts
            .iter()
            .any(|item| item.kind == ArtifactKind::FastqSingle)
    );
    assert!(
        !artifacts
            .iter()
            .any(|item| item.kind == ArtifactKind::Trash),
        "retention should not record trash artifacts"
    );
    let manifest = paths
        .with_generation_heads(&db)
        .unwrap()
        .quant_run_dir("RUN1", "kallistors")
        .join("rnaa_quant_manifest.json");
    let generation = compute_sha256(&manifest).unwrap();
    assert!(!db.retained_inputs("RUN1").unwrap().is_empty());
    run_ok(
        &["quant", "--root", root.to_str().unwrap(), "--preprocess"],
        None,
    );
    assert_eq!(
        compute_sha256(&manifest).unwrap(),
        generation,
        "authorized retention must preserve cache reuse"
    );
    // Updating a retention policy on already completed work must also execute it.
    run_ok(
        &[
            "quant",
            "--root",
            root.to_str().unwrap(),
            "--fastq-retention",
            "none",
        ],
        None,
    );
    assert!(!trimmed_fastq.exists());
    assert_eq!(compute_sha256(&manifest).unwrap(), generation);
}

#[test]
fn quant_retention_none_removes_both_original_and_trimmed_fastq() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        None,
    );
    let (paths, db, project_id) = setup_project(root);
    seed_trimmed_fastq(&paths, &db, &project_id);

    run_ok(
        &[
            "quant",
            "--root",
            root.to_str().expect("invalid root path"),
            "--preprocess",
            "--fastq-retention",
            "none",
        ],
        None,
    );

    let raw_fastq = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    let trimmed_fastq = paths
        .preprocess_run_dir("RUN1", "kallistors")
        .join("trimmed.fastq");

    let artifacts = db
        .list_artifacts_for_run(Some("RUN1"))
        .expect("failed to list artifacts");
    assert!(!raw_fastq.exists(), "raw FASTQ should be removed");
    assert!(!trimmed_fastq.exists(), "trimmed FASTQ should be removed");
    assert!(
        !paths.raw_run_dir("RUN1").exists(),
        "empty raw run dir should be pruned"
    );
    assert!(
        paths
            .preprocess_run_dir("RUN1", "kallistors")
            .join("fasterp.json")
            .exists(),
        "QC evidence must survive FASTQ cleanup"
    );
    assert!(!artifacts.iter().any(|item| {
        matches!(
            item.kind,
            ArtifactKind::FastqSingle | ArtifactKind::FastqTrimmedSingle
        )
    }));
    assert!(
        !artifacts
            .iter()
            .any(|item| item.kind == ArtifactKind::Trash),
        "retention should not record trash artifacts"
    );
    let manifest = paths
        .with_generation_heads(&db)
        .unwrap()
        .quant_run_dir("RUN1", "kallistors")
        .join("rnaa_quant_manifest.json");
    let generation = compute_sha256(&manifest).unwrap();
    assert!(!db.retained_inputs("RUN1").unwrap().is_empty());
    run_ok(
        &["quant", "--root", root.to_str().unwrap(), "--preprocess"],
        None,
    );
    assert_eq!(
        compute_sha256(&manifest).unwrap(),
        generation,
        "authorized retention must preserve cache reuse"
    );
}

#[test]
fn accidental_deletion_cannot_be_reconciled_as_success() {
    let temp = TempDir::new().unwrap();
    let root = temp.path();
    run_ok(&["init", "--root", root.to_str().unwrap()], None);
    let (paths, db, _) = setup_project(root);
    run_ok(&["quant", "--root", root.to_str().unwrap()], None);
    fs::remove_file(paths.raw_run_dir("RUN1").join("RUN1.fastq")).unwrap();
    let result = run_cmd(&["quant", "--root", root.to_str().unwrap()], None);
    assert!(!result.status.success());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantFailed
    );
    assert!(db.retained_inputs("RUN1").unwrap().is_empty());
}

#[test]
fn failed_native_replacement_preserves_outputs_and_retention_follows_atomic_completion() {
    use rnaa_core::generations::OutputGroup;
    let temp = TempDir::new().unwrap();
    let root = temp.path();
    let root_str = root.to_str().unwrap();
    run_ok(&["init", "--root", root_str], None);
    let (paths, db, project_id) = setup_project(root);
    run_ok(&["quant", "--root", root_str], None);
    let first = db
        .generation_head(OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    let old_files = db
        .list_artifacts_for_run(Some("RUN1"))
        .unwrap()
        .into_iter()
        .filter(|artifact| Path::new(&artifact.path).starts_with(&first.payload_path))
        .map(|artifact| {
            let contents = fs::read(&artifact.path).unwrap();
            (PathBuf::from(artifact.path), contents)
        })
        .collect::<Vec<_>>();
    assert_eq!(old_files.len(), 4);
    let assert_old_unchanged = || {
        for (path, contents) in &old_files {
            assert_eq!(fs::read(path).unwrap(), *contents);
        }
    };
    run_ok(&["quant", "--root", root_str, "--threads", "2"], None);
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        first
    );
    let raw = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    fs::write(&raw, "@broken\nACGT\n+\n!\n").unwrap();
    let result = run_cmd(
        &["quant", "--root", root_str, "--fastq-retention", "none"],
        None,
    );
    assert!(!result.status.success());
    assert!(raw.exists(), "failed replacement must not prune its inputs");
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        first
    );
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1").unwrap());
    assert_old_unchanged();
    assert!(
        !db.list_job_attempts(100)
            .unwrap()
            .iter()
            .any(|attempt| attempt.stage == "retention")
    );
    fs::write(&raw, "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n").unwrap();
    db.record_artifact(&ArtifactRecord {
        project_id,
        run_accession: Some("RUN1".into()),
        kind: ArtifactKind::FastqSingle,
        path: raw.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".into(),
        checksum: compute_sha256(&raw).unwrap(),
        bytes: file_size(&raw).unwrap(),
        created_at: now_rfc3339(),
    })
    .unwrap();
    run_ok(&["quant", "--root", root_str], None);
    let next = db
        .generation_head(OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    assert_ne!(next.generation_id, first.generation_id);
    assert!(!raw.exists());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantDone
    );
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1").unwrap());
    assert_old_unchanged();
    let artifacts = db.list_artifacts_for_run(Some("RUN1")).unwrap();
    assert_eq!(artifacts.len(), 4);
    assert!(
        artifacts
            .iter()
            .all(|artifact| Path::new(&artifact.path).starts_with(&next.payload_path))
    );
    assert!(
        artifacts
            .iter()
            .any(|artifact| artifact.kind == ArtifactKind::QuantAbundanceTsv)
    );
    let attempts = db.list_job_attempts(100).unwrap();
    assert_eq!(
        attempts
            .iter()
            .find(|attempt| attempt.attempt_id == next.attempt_id)
            .unwrap()
            .state,
        "SUCCEEDED"
    );
    let retention = attempts
        .iter()
        .find(|attempt| attempt.stage == "retention")
        .unwrap();
    assert_eq!(retention.state, "SUCCEEDED");
    assert!(retention.request_sha256.is_some());
    let current = paths.with_generation_heads(&db).unwrap();
    let manifest = current
        .quant_run_dir("RUN1", "kallistors")
        .join("rnaa_quant_manifest.json");
    let manifest_digest = compute_sha256(&manifest).unwrap();
    assert!(
        db.retained_inputs("RUN1")
            .unwrap()
            .iter()
            .all(|input| input.generation_sha256 == manifest_digest)
    );
    run_ok(&["quant", "--root", root_str], None);
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        next
    );
    fs::remove_file(next.payload_path.parent().unwrap().join("seal.json")).unwrap();
    let result = run_cmd(&["quant", "--root", root_str], None);
    assert!(
        !result.status.success(),
        "retained counts require their generation witness"
    );
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantFailed
    );
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        next
    );
    assert_old_unchanged();
}

#[test]
fn retention_from_another_generation_cannot_authorize_missing_input() {
    let temp = TempDir::new().unwrap();
    let root = temp.path();
    run_ok(&["init", "--root", root.to_str().unwrap()], None);
    let (paths, db, _) = setup_project(root);
    run_ok(&["quant", "--root", root.to_str().unwrap()], None);
    let raw = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    let snapshot = rnaa_core::ManifestArtifact::snapshot("FASTQ_SINGLE", &raw).unwrap();
    db.record_retained_input("RUN1", &"0".repeat(64), &snapshot)
        .unwrap();
    fs::remove_file(raw).unwrap();
    let result = run_cmd(&["quant", "--root", root.to_str().unwrap()], None);
    assert!(!result.status.success());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantFailed
    );
}

fn good_read() -> String {
    let sequence = "ACGTTCGATCGGATCGATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCTAGCTGACTGATCGATCGTAGCTAGCTGATCGTAGCTAGCTAGCATCGA";
    format!("@read\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()))
}

fn replace_raw(paths: &ProjectPaths, db: &Database, contents: &str) -> PathBuf {
    let raw = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    fs::create_dir_all(raw.parent().unwrap()).unwrap();
    fs::write(&raw, contents).unwrap();
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: Some("RUN1".into()),
        kind: ArtifactKind::FastqSingle,
        path: raw.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".into(),
        checksum: compute_sha256(&raw).unwrap(),
        bytes: file_size(&raw).unwrap(),
        created_at: now_rfc3339(),
    })
    .unwrap();
    raw
}

#[test]
fn preprocessing_generations_preserve_prior_files_and_authorize_only_consumed_trimmed_reads() {
    use rnaa_core::generations::OutputGroup;
    let temp = TempDir::new().unwrap();
    let root = temp.path().to_str().unwrap();
    run_ok(&["init", "--root", root], None);
    let (paths, db, _) = setup_project(temp.path());
    replace_raw(&paths, &db, &good_read().repeat(3));
    run_ok(
        &[
            "quant",
            "--root",
            root,
            "--preprocess",
            "--worker-budget",
            "1",
        ],
        None,
    );
    let first = db
        .generation_head(OutputGroup::Preprocess, "RUN1")
        .unwrap()
        .unwrap();
    let original = db
        .list_artifacts_for_run(Some("RUN1"))
        .unwrap()
        .into_iter()
        .filter(|record| Path::new(&record.path).starts_with(&first.payload_path))
        .map(|record| {
            let contents = fs::read(&record.path).unwrap();
            (PathBuf::from(record.path), contents)
        })
        .collect::<Vec<_>>();
    assert_eq!(original.len(), 3);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.qc_gate.min_pass_rate = Some(0.7);
    config.save(&paths.config_path).unwrap();
    run_ok(&["quant", "--root", root, "--worker-budget", "1"], None);
    assert_eq!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .unwrap(),
        first,
        "a policy reevaluation can reuse outputs without writing into the published generation"
    );
    replace_raw(&paths, &db, &good_read().repeat(4));
    run_ok(
        &[
            "quant",
            "--root",
            root,
            "--worker-budget",
            "1",
            "--fastq-retention",
            "original",
        ],
        None,
    );
    let next = db
        .generation_head(OutputGroup::Preprocess, "RUN1")
        .unwrap()
        .unwrap();
    assert_ne!(next.generation_id, first.generation_id);
    let current = paths.with_generation_heads(&db).unwrap();
    let trimmed = current
        .preprocess_run_dir("RUN1", "kallistors")
        .join("trimmed.fastq");
    assert!(!trimmed.exists());
    assert!(
        db.verify_generation(OutputGroup::Preprocess, "RUN1")
            .unwrap()
    );
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1").unwrap());
    assert!(
        db.retained_inputs("RUN1")
            .unwrap()
            .iter()
            .any(|record| record.artifact.path == trimmed.display().to_string())
    );
    for (path, bytes) in &original {
        assert_eq!(fs::read(path).unwrap(), *bytes);
    }
    assert!(
        !db.list_artifacts_for_run(Some("RUN1"))
            .unwrap()
            .iter()
            .any(|record| Path::new(&record.path).starts_with(&first.payload_path))
    );
    run_ok(&["quant", "--root", root], None);
    assert_eq!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .unwrap(),
        next
    );
    let report = current
        .preprocess_run_dir("RUN1", "kallistors")
        .join("fasterp.json");
    fs::write(&report, "corrupted").unwrap();
    assert!(
        !db.verify_generation(OutputGroup::Preprocess, "RUN1")
            .unwrap()
    );
    run_ok(&["quant", "--root", root, "--worker-budget", "1"], None);
    assert_ne!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .unwrap(),
        next
    );
    assert_eq!(
        fs::read_to_string(report).unwrap(),
        "corrupted",
        "repair must publish a new generation"
    );
    assert!(
        db.verify_generation(OutputGroup::Preprocess, "RUN1")
            .unwrap()
    );
    for (path, bytes) in &original {
        assert_eq!(fs::read(path).unwrap(), *bytes);
    }
}

#[test]
fn all_rejected_reads_publish_failed_qc_without_engine_error_bypass() {
    use rnaa_core::generations::OutputGroup;
    let temp = TempDir::new().unwrap();
    let root = temp.path().to_str().unwrap();
    run_ok(&["init", "--root", root], None);
    let (paths, db, _) = setup_project(temp.path());
    let args = [
        "quant",
        "--root",
        root,
        "--preprocess",
        "--preprocess-bypass",
        "--worker-budget",
        "1",
    ];
    let result = run_cmd(&args, None);
    assert!(!result.status.success());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::PreprocessDone
    );
    assert!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .is_none()
    );
    let head = db
        .generation_head(OutputGroup::Preprocess, "RUN1")
        .unwrap()
        .unwrap();
    assert!(
        db.verify_generation(OutputGroup::Preprocess, "RUN1")
            .unwrap()
    );
    let current = paths.with_generation_heads(&db).unwrap();
    assert_eq!(
        file_size(
            &current
                .preprocess_run_dir("RUN1", "kallistors")
                .join("trimmed.fastq")
        )
        .unwrap(),
        0
    );
    let qc = db.get_preprocess_qc("RUN1").unwrap().unwrap();
    assert_eq!((qc.passed_reads, qc.failed_reads), (0, 1));
    assert!(!qc.gate_passed);
    let attempt = db.list_job_attempts(1).unwrap().remove(0);
    assert_eq!(attempt.state, "FAILED");
    assert!(attempt.last_error.unwrap().contains("no reads passed"));
    let report = current
        .preprocess_run_dir("RUN1", "kallistors")
        .join("fasterp.json");
    let bytes = fs::read(&report).unwrap();
    let result = run_cmd(&args, None);
    assert!(!result.status.success());
    assert_eq!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .unwrap(),
        head
    );
    assert_eq!(fs::read(report).unwrap(), bytes);
    assert_eq!(
        db.list_job_attempts(100)
            .unwrap()
            .iter()
            .filter(|attempt| attempt.stage == "preprocess")
            .count(),
        1,
        "an unchanged QC rejection must not rerun preprocessing"
    );
    assert_eq!(
        db.eligibility("RUN1").unwrap().unwrap().descriptor.state,
        rnaa_core::eligibility::EligibilityState::Rejected
    );
    assert!(
        !db.list_events(100)
            .unwrap()
            .iter()
            .any(|event| event.message.contains("preprocessing bypassed"))
    );
}

#[test]
fn explicit_engine_error_bypass_retires_stale_trimmed_reads_and_strict_failure_preserves_prior_generations()
 {
    use rnaa_core::generations::OutputGroup;
    let temp = TempDir::new().unwrap();
    let root = temp.path().to_str().unwrap();
    run_ok(&["init", "--root", root], None);
    let (paths, db, _) = setup_project(temp.path());
    replace_raw(&paths, &db, &good_read());
    run_ok(
        &[
            "quant",
            "--root",
            root,
            "--preprocess",
            "--worker-budget",
            "1",
        ],
        None,
    );
    let first = db
        .generation_head(OutputGroup::Preprocess, "RUN1")
        .unwrap()
        .unwrap();
    let old = paths
        .with_generation_heads(&db)
        .unwrap()
        .preprocess_run_dir("RUN1", "kallistors")
        .join("trimmed.fastq");
    let old_bytes = fs::read(&old).unwrap();
    replace_raw(&paths, &db, &good_read().repeat(6000));
    run_ok(
        &[
            "quant",
            "--root",
            root,
            "--preprocess-bypass",
            "--preprocess-max-mb",
            "1",
            "--worker-budget",
            "1",
        ],
        None,
    );
    let bypass = db
        .generation_head(OutputGroup::Preprocess, "RUN1")
        .unwrap()
        .unwrap();
    assert_ne!(bypass, first);
    assert!(
        db.verify_generation(OutputGroup::Preprocess, "RUN1")
            .unwrap()
    );
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1").unwrap());
    let qc = db.get_preprocess_qc("RUN1").unwrap().unwrap();
    assert_eq!(qc.mode, "bypass");
    assert!(qc.gate_passed);
    assert!(
        !db.list_artifacts_for_run(Some("RUN1"))
            .unwrap()
            .iter()
            .any(|record| matches!(
                record.kind,
                ArtifactKind::FastqTrimmedSingle
                    | ArtifactKind::FastqTrimmedR1
                    | ArtifactKind::FastqTrimmedR2
            ))
    );
    let current = paths.with_generation_heads(&db).unwrap();
    let manifest: rnaa_core::StageManifest = serde_json::from_slice(
        &fs::read(
            current
                .quant_run_dir("RUN1", "kallistors")
                .join("rnaa_quant_manifest.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(manifest.input_artifacts.len(), 1);
    assert_eq!(manifest.input_artifacts[0].kind, "FASTQ_SINGLE");
    let info: serde_json::Value = serde_json::from_slice(
        &fs::read(
            current
                .quant_run_dir("RUN1", "kallistors")
                .join("run_info.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(info["n_processed"], 6000);
    assert_eq!(fs::read(&old).unwrap(), old_bytes);
    let quant = db
        .generation_head(OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    let result = run_cmd(
        &[
            "quant",
            "--root",
            root,
            "--preprocess-strict",
            "--worker-budget",
            "1",
        ],
        None,
    );
    assert!(!result.status.success());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::PreprocessFailed
    );
    assert_eq!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .unwrap(),
        bypass
    );
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        quant
    );
    assert_eq!(fs::read(old).unwrap(), old_bytes);
}

#[test]
fn generation_retention_preserves_reports_and_does_not_mask_report_corruption() {
    use rnaa_core::generations::OutputGroup;
    let temp = TempDir::new().unwrap();
    let root = temp.path().to_str().unwrap();
    run_ok(&["init", "--root", root], None);
    let (paths, db, _) = setup_project(temp.path());
    let raw = replace_raw(&paths, &db, &good_read().repeat(3));
    run_ok(
        &[
            "quant",
            "--root",
            root,
            "--preprocess",
            "--fastq-retention",
            "none",
            "--worker-budget",
            "1",
        ],
        None,
    );
    let preprocess = db
        .generation_head(OutputGroup::Preprocess, "RUN1")
        .unwrap()
        .unwrap();
    let quant = db
        .generation_head(OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    let current = paths.with_generation_heads(&db).unwrap();
    let report = current
        .preprocess_run_dir("RUN1", "kallistors")
        .join("fasterp.json");
    assert!(!raw.exists());
    assert!(
        !current
            .preprocess_run_dir("RUN1", "kallistors")
            .join("trimmed.fastq")
            .exists()
    );
    assert!(report.is_file());
    assert!(
        db.verify_generation(OutputGroup::Preprocess, "RUN1")
            .unwrap()
    );
    assert!(db.verify_generation(OutputGroup::Quant, "RUN1").unwrap());
    run_ok(&["quant", "--root", root], None);
    assert_eq!(
        db.generation_head(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .unwrap(),
        preprocess
    );
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        quant
    );
    run_ok(&["qc", "reevaluate", "--root", root], None);
    let before = db.get_preprocess_qc("RUN1").unwrap().unwrap();
    let mut lease = rnaa_core::jobs::JobLease::try_acquire(
        &db,
        "preprocess",
        Some("RUN1"),
        &[RunState::QuantDone],
    )
    .unwrap()
    .unwrap();
    let result = run_cmd(&["qc", "reevaluate", "--root", root], None);
    assert!(result.status.success());
    assert!(String::from_utf8_lossy(&result.stdout).contains("qc_busy_runs\t1"));
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::PreprocessRunning
    );
    assert_eq!(
        db.get_preprocess_qc("RUN1").unwrap().unwrap().updated_at,
        before.updated_at
    );
    lease
        .database()
        .set_run_state("RUN1", RunState::QuantDone, None)
        .unwrap();
    lease
        .finish(rnaa_core::jobs::JobOutcome::Succeeded, None)
        .unwrap();
    drop(lease);
    let mut analysis = rnaa_core::jobs::JobLease::try_acquire(&db, "normalize", None, &[])
        .unwrap()
        .unwrap();
    let result = run_cmd(&["qc", "reevaluate", "--root", root], None);
    assert!(!result.status.success());
    assert_eq!(
        db.get_preprocess_qc("RUN1").unwrap().unwrap().updated_at,
        before.updated_at
    );
    analysis
        .finish(rnaa_core::jobs::JobOutcome::Failed, Some("test finished"))
        .unwrap();
    drop(analysis);
    fs::write(&report, "damaged").unwrap();
    let result = run_cmd(&["qc", "reevaluate", "--root", root], None);
    assert!(!result.status.success());
    assert_eq!(
        db.get_preprocess_qc("RUN1").unwrap().unwrap().updated_at,
        before.updated_at
    );
    let result = run_cmd(&["quant", "--root", root], None);
    assert!(!result.status.success());
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantDone
    );
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        quant
    );
}

#[test]
fn downloaded_and_preprocessed_generations_remain_verified_after_all_fastqs_are_pruned() {
    use rnaa_core::{RemoteFile, RemoteFileKind, generations::OutputGroup};
    use std::{
        io::{BufRead, BufReader},
        net::TcpListener,
        thread,
    };
    let temp = TempDir::new().unwrap();
    let root = temp.path().to_str().unwrap();
    run_ok(&["init", "--root", root], None);
    let (paths, db, _) = setup_project(temp.path());
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}/download.fastq", listener.local_addr().unwrap());
    let body = good_read().repeat(3);
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.remote_files = vec![RemoteFile {
        url,
        basename: "download.fastq".into(),
        kind: RemoteFileKind::Fastq,
        bytes: Some(body.len() as u64),
        md5: None,
    }];
    db.upsert_run(&run).unwrap();
    db.set_run_state("RUN1", RunState::Resolved, None).unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.download.retries = 0;
    config.save(&paths.config_path).unwrap();
    let server = thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        stream
            .set_read_timeout(Some(std::time::Duration::from_secs(5)))
            .unwrap();
        let mut reader = BufReader::new(&mut stream);
        loop {
            let mut line = String::new();
            assert!(reader.read_line(&mut line).unwrap() > 0);
            if line == "\r\n" {
                break;
            }
        }
        write!(stream,"HTTP/1.1 200 OK\r\nContent-Length: {}\r\nETag: \"reads-v1\"\r\nConnection: close\r\n\r\n",body.len()).unwrap();
        stream.write_all(body.as_bytes()).unwrap();
    });
    run_ok(&["download", "--root", root, "--concurrency", "1"], None);
    server.join().unwrap();
    let download = db
        .generation_head(OutputGroup::Download, "RUN1")
        .unwrap()
        .unwrap();
    let downloaded = paths.with_generation_heads(&db).unwrap();
    let raw = downloaded.raw_run_dir("RUN1").join("download.fastq");
    assert!(raw.is_file());
    run_ok(
        &[
            "quant",
            "--root",
            root,
            "--preprocess",
            "--fastq-retention",
            "none",
            "--worker-budget",
            "1",
        ],
        None,
    );
    assert!(!raw.exists());
    for group in [
        OutputGroup::Download,
        OutputGroup::Preprocess,
        OutputGroup::Quant,
    ] {
        assert!(db.verify_generation(group, "RUN1").unwrap());
    }
    assert!(downloaded.manifest_path("download", "RUN1").is_file());
    let quant = db
        .generation_head(OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    // The HTTP server is gone: retention evidence must prevent a false redownload.
    run_ok(&["download", "--root", root, "--concurrency", "1"], None);
    run_ok(&["quant", "--root", root], None);
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantDone
    );
    assert_eq!(
        db.generation_head(OutputGroup::Download, "RUN1")
            .unwrap()
            .unwrap(),
        download
    );
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        quant
    );
    assert_eq!(
        db.list_job_attempts(100)
            .unwrap()
            .iter()
            .filter(|attempt| attempt.stage == "download")
            .count(),
        1
    );
}

#[test]
fn qc_policy_changes_preserve_counts_and_need_no_reads_under_every_retention_policy() {
    use rnaa_core::{eligibility::EligibilityState, generations::OutputGroup};
    for retention in ["both", "original", "trimmed", "none"] {
        let temp = TempDir::new().unwrap();
        let root = temp.path().to_str().unwrap();
        run_ok(&["init", "--root", root], None);
        let (paths, db, _) = setup_project(temp.path());
        let bad = format!("@bad\n{}\n+\n{}\n", "ACGT".repeat(25), "!".repeat(100));
        replace_raw(&paths, &db, &format!("{}{bad}", good_read().repeat(4)));
        let args = [
            "quant",
            "--root",
            root,
            "--preprocess",
            "--fastq-retention",
            retention,
            "--worker-budget",
            "1",
        ];
        run_ok(&args, None);
        let head = db
            .generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap();
        let current = paths.with_generation_heads(&db).unwrap();
        let manifest = current
            .quant_run_dir("RUN1", "kallistors")
            .join("rnaa_quant_manifest.json");
        let bytes = fs::read(&manifest).unwrap();
        let attempts = db.list_job_attempts(100).unwrap().len();
        let history = db.eligibility_history("RUN1").unwrap().len();
        let mut config = ProjectConfig::load(&paths.config_path).unwrap();
        config.quant.qc_gate.min_pass_rate = Some(0.95);
        config.save(&paths.config_path).unwrap();
        let result = run_cmd(&["quant", "--root", root], None);
        assert!(!result.status.success(), "{retention}");
        assert_eq!(
            db.get_run("RUN1").unwrap().unwrap().state,
            RunState::QuantDone
        );
        assert_eq!(
            db.eligibility("RUN1").unwrap().unwrap().descriptor.state,
            EligibilityState::Rejected
        );
        config.quant.qc_gate.min_pass_rate = Some(0.6);
        config.save(&paths.config_path).unwrap();
        run_ok(&["quant", "--root", root], None);
        assert_eq!(
            db.eligibility("RUN1").unwrap().unwrap().descriptor.state,
            EligibilityState::Eligible
        );
        assert_eq!(
            db.generation_head(OutputGroup::Quant, "RUN1")
                .unwrap()
                .unwrap(),
            head
        );
        assert_eq!(fs::read(&manifest).unwrap(), bytes);
        assert_eq!(
            db.list_job_attempts(100).unwrap().len(),
            attempts,
            "QC policy alone must not rerun an engine or cleanup"
        );
        assert_eq!(db.eligibility_history("RUN1").unwrap().len(), history + 2);
        assert!(
            db.verify_generation(OutputGroup::Preprocess, "RUN1")
                .unwrap()
        );
    }
}
