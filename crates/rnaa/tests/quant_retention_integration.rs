use std::collections::BTreeMap;
use std::fs;
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

fn make_mock_rscript(bin_dir: &Path) {
    let script = r#"#!/usr/bin/env bash
set -euo pipefail
outdir=""
manifest=""
run_id=""
reference_manifest=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) outdir="$2"; shift 2 ;;
    --manifest-out) manifest="$2"; shift 2 ;;
    --run) run_id="$2"; shift 2 ;;
    --reference-manifest) reference_manifest="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$outdir"
printf 'mock-h5' > "$outdir/abundance.h5"
cat > "$outdir/run_info.json" <<EOF
{"run":"$run_id"}
EOF
if [[ -n "${manifest}" ]]; then
  mkdir -p "$(dirname "$manifest")"
  cat > "$manifest" <<EOF
{"schema_version":1,"stage":"quant","parameters":{"reference_manifest":"$reference_manifest"}}
EOF
fi
"#;
    let path = bin_dir.join("Rscript");
    fs::write(&path, script).expect("failed to write mock Rscript");
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = fs::metadata(&path)
            .expect("failed to stat mock Rscript")
            .permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&path, perms).expect("failed to chmod mock Rscript");
    }
}

fn write_minimal_inputs(root: &Path) -> (PathBuf, PathBuf) {
    let cdna = root.join("test.cdna.fa.gz");
    let gtf = root.join("test.gtf.gz");
    fs::write(&cdna, b">tx1\nACGTACGTACGT\n").expect("failed to write cdna");
    fs::write(
        &gtf,
        b"chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; gene_biotype \"protein_coding\"; transcript_id \"tx1\";\n",
    )
    .expect("failed to write gtf");
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
    let trimmed_dir = paths.quant_run_dir("RUN1", "r-kallisto").join("fasterp");
    fs::create_dir_all(&trimmed_dir).expect("failed to create trimmed dir");
    let trimmed_fastq = trimmed_dir.join("trimmed.fastq");
    fs::write(&trimmed_fastq, b"@r1\nACGT\n+\n!!!!\n").expect("failed to write trimmed fastq");
    db.record_artifact(&ArtifactRecord {
        project_id: project_id.to_string(),
        run_accession: Some("RUN1".to_string()),
        kind: ArtifactKind::FastqTrimmedSingle,
        path: trimmed_fastq.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(&trimmed_fastq).expect("failed to hash trimmed fastq"),
        bytes: file_size(&trimmed_fastq).expect("failed to stat trimmed fastq"),
        created_at: now_rfc3339(),
    })
    .expect("failed to record trimmed artifact");
    db.set_run_state("RUN1", RunState::PreprocessDone, None)
        .expect("failed to set preprocess state");
}

#[test]
fn quant_retention_trimmed_deletes_original_fastq() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    let bin_dir = root.join("mock_bin");
    fs::create_dir_all(&bin_dir).expect("failed to create mock bin dir");
    make_mock_rscript(&bin_dir);

    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        Some(&bin_dir),
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
        Some(&bin_dir),
    );

    let raw_fastq = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    let trimmed_fastq = paths
        .quant_run_dir("RUN1", "r-kallisto")
        .join("fasterp")
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
}

#[test]
fn quant_retention_none_removes_both_original_and_trimmed_fastq() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    let bin_dir = root.join("mock_bin");
    fs::create_dir_all(&bin_dir).expect("failed to create mock bin dir");
    make_mock_rscript(&bin_dir);

    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        Some(&bin_dir),
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
        Some(&bin_dir),
    );

    let raw_fastq = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    let trimmed_fastq = paths
        .quant_run_dir("RUN1", "r-kallisto")
        .join("fasterp")
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
        !paths
            .quant_run_dir("RUN1", "r-kallisto")
            .join("fasterp")
            .exists(),
        "empty trimmed FASTQ dir should be pruned"
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
}
