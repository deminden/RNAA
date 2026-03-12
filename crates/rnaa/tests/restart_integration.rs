use std::collections::BTreeMap;
use std::fs;
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

fn make_mock_rscript(bin_dir: &Path) {
    let script = r#"#!/usr/bin/env bash
set -euo pipefail

mode=""
outdir=""
manifest=""
run_id=""
reference_manifest=""

for arg in "$@"; do
  case "$arg" in
    *quant_kallisto.R) mode="quant" ;;
    *tximport_deseq2.R) mode="deseq2" ;;
  esac
done

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

if [[ "$mode" == "quant" ]]; then
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
  exit 0
fi

cat > "$outdir/gene_counts.tsv" <<'EOF'
gene_id	RUN1	RUN2
G1	100	120
G2	50	60
EOF
cat > "$outdir/gene_counts_annotated.tsv" <<'EOF'
gene_id	gene_symbol	gene_biotype	RUN1	RUN2
G1	Gene1	protein_coding	100	120
G2	Gene2	protein_coding	50	60
EOF
cat > "$outdir/gene_norm_counts.tsv" <<'EOF'
gene_id	RUN1	RUN2
G1	10	12
G2	5	6
EOF
cat > "$outdir/gene_norm_counts_annotated.tsv" <<'EOF'
gene_id	gene_symbol	gene_biotype	RUN1	RUN2
G1	Gene1	protein_coding	10	12
G2	Gene2	protein_coding	5	6
EOF
cat > "$outdir/vst.tsv" <<'EOF'
gene_id	RUN1	RUN2
G1	2.1	2.4
G2	1.5	1.7
EOF
cat > "$outdir/vst_annotated.tsv" <<'EOF'
gene_id	gene_symbol	gene_biotype	RUN1	RUN2
G1	Gene1	protein_coding	2.1	2.4
G2	Gene2	protein_coding	1.5	1.7
EOF
printf 'placeholder' > "$outdir/vst.rds"
cat > "$outdir/size_factors.tsv" <<'EOF'
run_accession	size_factor
RUN1	1.0
RUN2	1.0
EOF
if [[ -n "${manifest}" ]]; then
  mkdir -p "$(dirname "$manifest")"
  cat > "$manifest" <<'EOF'
{"schema_version":1,"stage":"deseq2"}
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
    fs::write(&cdna, b">tx1\nACGTACGTACGT\n>tx2\nTTTTAACCGG\n").expect("failed to write cdna");
    fs::write(
        &gtf,
        b"chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; gene_biotype \"protein_coding\"; transcript_id \"tx1\";\nchr1\tsrc\tgene\t101\t200\t.\t+\t.\tgene_id \"G2\"; gene_name \"Gene2\"; gene_biotype \"protein_coding\"; transcript_id \"tx2\";\n",
    )
    .expect("failed to write gtf");
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
    let quant_dir = paths.quant_run_dir(run_accession, "r-kallisto");
    fs::create_dir_all(&quant_dir).expect("failed to create quant dir");
    let abundance_h5 = quant_dir.join("abundance.h5");
    let run_info_json = quant_dir.join("run_info.json");
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
}

#[test]
fn quant_restart_recovers_preprocess_and_completes() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    let bin_dir = root.join("mock_bin");
    fs::create_dir_all(&bin_dir).expect("failed to create mock bin dir");
    make_mock_rscript(&bin_dir);

    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        Some(&bin_dir),
    );

    let (paths, db, project_id) = setup_custom_reference(root);

    let raw1 = paths.raw_run_dir("RUN1").join("RUN1.fastq");
    let raw2 = paths.raw_run_dir("RUN2").join("RUN2.fastq");
    fs::create_dir_all(raw1.parent().expect("missing parent")).expect("failed to create raw dir");
    fs::create_dir_all(raw2.parent().expect("missing parent")).expect("failed to create raw dir");
    fs::write(&raw1, b"@r1\nACGT\n+\n!!!!\n").expect("failed to write fastq");
    fs::write(&raw2, b"@r2\nTGCA\n+\n!!!!\n").expect("failed to write fastq");

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
        Some(&bin_dir),
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
        runs.iter().all(|run| run.state == RunState::QuantDone),
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
            .quant_run_dir("RUN1", "r-kallisto")
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
}

#[test]
fn normalize_and_corr_restart_recover_project_stages() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    let bin_dir = root.join("mock_bin");
    fs::create_dir_all(&bin_dir).expect("failed to create mock bin dir");
    make_mock_rscript(&bin_dir);

    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        Some(&bin_dir),
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
            metadata: BTreeMap::new(),
            state: RunState::DeRunning,
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

    run_ok(
        &[
            "normalize",
            "--root",
            root.to_str().expect("invalid root path"),
            "--design",
            "~ batch + condition",
        ],
        Some(&bin_dir),
    );

    let reopened = Database::open(root).expect("failed to reopen db");
    let normalize_stage = reopened
        .get_project_stage_state("normalize")
        .expect("failed to load normalize stage")
        .expect("normalize stage missing");
    assert_eq!(normalize_stage.state, ProjectStageState::Done);
    assert!(
        reopened
            .list_runs()
            .expect("failed to list runs")
            .iter()
            .all(|run| run.state == RunState::DeDone)
    );

    for run in reopened.list_runs().expect("failed to list runs") {
        reopened
            .set_run_state(&run.run_accession, RunState::CorrRunning, None)
            .expect("failed to set corr state");
    }
    reopened
        .set_project_stage_state("corr", ProjectStageState::Running, None)
        .expect("failed to set corr stage");

    run_ok(
        &[
            "corr",
            "--root",
            root.to_str().expect("invalid root path"),
            "--model",
            "~ batch + condition",
            "--method",
            "pearson",
        ],
        Some(&bin_dir),
    );

    let final_db = Database::open(root).expect("failed to reopen final db");
    let corr_stage = final_db
        .get_project_stage_state("corr")
        .expect("failed to load corr stage")
        .expect("corr stage missing");
    assert_eq!(corr_stage.state, ProjectStageState::Done);
    assert!(
        final_db
            .list_runs()
            .expect("failed to list runs")
            .iter()
            .all(|run| run.state == RunState::CorrDone)
    );
    assert!(
        final_db
            .list_events(200)
            .expect("failed to list events")
            .iter()
            .any(|event| event.stage == "recovery" && event.message == "interrupted work recovered")
    );
}
