use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::model::{LibraryLayout, RunRecord};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::RunState;
use rnaa_core::util::{compute_sha256, now_rfc3339};
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
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) outdir="$2"; shift 2;;
    --manifest-out) manifest="$2"; shift 2;;
    *) shift;;
  esac
done
mkdir -p "$outdir"
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
        b"chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"Gene1\"; gene_biotype \"protein_coding\";\n",
    )
    .expect("failed to write gtf");
    (cdna, gtf)
}

fn setup_project(root: &Path, include_batch: bool) -> (ProjectPaths, String) {
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).expect("failed to open database");
    let project = db.project().expect("failed to load project");
    let (cdna, gtf) = write_minimal_inputs(root);

    let mut config = ProjectConfig::load(&paths.config_path).expect("failed to load config");
    config.refs.custom_cdna = cdna.display().to_string();
    config.refs.custom_gtf = gtf.display().to_string();
    config
        .save(&paths.config_path)
        .expect("failed to save project config");
    db.set_project_config(&config)
        .expect("failed to persist config");

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

    let runs = ["RUN1", "RUN2"]
        .into_iter()
        .map(|run| RunRecord {
            project_id: project.project_id.clone(),
            study_accession: "SRPTEST".to_string(),
            source_study_accession: None,
            run_accession: run.to_string(),
            sample_accession: Some(format!("SAMP_{run}")),
            experiment_accession: Some(format!("EXP_{run}")),
            library_layout: LibraryLayout::Single,
            instrument_platform: Some("ILLUMINA".to_string()),
            instrument_model: Some("TEST".to_string()),
            remote_files: Vec::new(),
            metadata: std::collections::BTreeMap::new(),
            state: RunState::QuantDone,
            last_error: None,
            updated_at: now_rfc3339(),
        })
        .collect::<Vec<_>>();
    db.upsert_runs(&runs).expect("failed to upsert runs");

    let samplesheet = if include_batch {
        "project_id\tstudy_accession\trun_accession\tcondition\tbatch\n\
         test\tSRPTEST\tRUN1\tA\tB1\n\
         test\tSRPTEST\tRUN2\tB\tB2\n"
    } else {
        "project_id\tstudy_accession\trun_accession\tcondition\n\
         test\tSRPTEST\tRUN1\tA\n\
         test\tSRPTEST\tRUN2\tB\n"
    };
    fs::write(paths.samplesheet_path(), samplesheet).expect("failed to write samplesheet");

    (paths, project.project_id)
}

#[test]
fn normalize_writes_outputs_and_marks_runs_done() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    let bin_dir = root.join("mock_bin");
    fs::create_dir_all(&bin_dir).expect("failed to create mock bin dir");
    make_mock_rscript(&bin_dir);

    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        Some(&bin_dir),
    );
    let (paths, project_id) = setup_project(root, true);

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

    let de_dir = paths.de_project_dir(&project_id);
    assert!(de_dir.join("gene_counts.tsv").exists());
    assert!(de_dir.join("gene_counts_annotated.tsv").exists());
    assert!(de_dir.join("gene_norm_counts.tsv").exists());
    assert!(de_dir.join("gene_norm_counts_annotated.tsv").exists());
    assert!(de_dir.join("vst.tsv").exists());
    assert!(de_dir.join("vst_annotated.tsv").exists());
    assert!(de_dir.join("deseq2_manifest.json").exists());

    let db = Database::open(root).expect("failed to reopen db");
    let runs = db.list_runs().expect("failed to list runs");
    assert!(runs.iter().all(|run| run.state == RunState::DeDone));
}

#[test]
fn normalize_fails_when_design_references_missing_column() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    let bin_dir = root.join("mock_bin");
    fs::create_dir_all(&bin_dir).expect("failed to create mock bin dir");
    make_mock_rscript(&bin_dir);

    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        Some(&bin_dir),
    );
    let _ = setup_project(root, false);

    let output = run_cmd(
        &[
            "normalize",
            "--root",
            root.to_str().expect("invalid root path"),
            "--design",
            "~ batch + condition",
        ],
        Some(&bin_dir),
    );
    assert!(!output.status.success(), "normalize unexpectedly succeeded");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("missing samplesheet columns: batch"));

    let db = Database::open(root).expect("failed to reopen db");
    let runs = db.list_runs().expect("failed to list runs");
    assert!(runs.iter().all(|run| run.state == RunState::QuantDone));
}
