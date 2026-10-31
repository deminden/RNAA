use rnaa_core::traits::ReferenceManager;
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
>tx2
TTTTAACCGG
"#,
    );
    write_gzip(&gtf, br#"chr1	src	transcript	1	100	.	+	.	gene_id "G1"; gene_name "Gene1"; gene_biotype "protein_coding"; transcript_id "tx1";
chr1	src	transcript	101	200	.	+	.	gene_id "G2"; gene_name "Gene2"; gene_biotype "protein_coding"; transcript_id "tx2";
"#);
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
            metadata: std::collections::BTreeMap::from([
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
    db.upsert_runs(&runs).expect("failed to upsert runs");

    for run in ["RUN1", "RUN2"] {
        let quant_dir = paths.quant_run_dir(run, "kallistors");
        fs::create_dir_all(&quant_dir).expect("failed to create quant dir");
        let abundance_h5 = quant_dir.join("abundance.h5");
        let run_info_json = quant_dir.join("run_info.json");
        fs::write(
            quant_dir.join("abundance.tsv"),
            format!(
                "target_id\tlength\teff_length\test_counts\ttpm\ntx1\t12\t10\t{}\t400000\ntx2\t10\t8\t{}\t600000\n",
                if run == "RUN1" { 100 } else { 120 },
                if run == "RUN1" { 50 } else { 60 },
            ),
        )
        .expect("failed to write abundance.tsv");
        fs::write(&abundance_h5, b"h5").expect("failed to write abundance.h5");
        fs::write(&run_info_json, b"{\"run\":\"ok\"}").expect("failed to write run_info.json");
        db.record_artifact(&ArtifactRecord {
            project_id: project.project_id.clone(),
            run_accession: Some(run.to_string()),
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
            project_id: project.project_id.clone(),
            run_accession: Some(run.to_string()),
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

    let reference = rnaa_stages::refs::EnsemblReferenceManager::new()
        .unwrap()
        .ensure_reference(&paths, &db, &config)
        .unwrap();
    for run in db.list_runs().unwrap() {
        rnaa_stages::quant::declare_external_quant(
            &run.run_accession,
            &paths.quant_run_dir(&run.run_accession, "kallistors"),
            &reference,
            "external underpowered statistical fixture",
        )
        .unwrap();
    }
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
fn normalize_excludes_and_demotes_runs_with_invalid_quant_outputs() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        None,
    );
    let (paths, _project_id) = setup_project(root, true);

    fs::remove_file(
        paths
            .quant_run_dir("RUN2", "kallistors")
            .join("abundance.tsv"),
    )
    .expect("failed to remove broken abundance");

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
    assert!(!output.status.success(), "normalize unexpectedly succeeded");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("need at least 2 valid quantified runs"));

    let db = Database::open(root).expect("failed to reopen db");
    let runs = db.list_runs().expect("failed to list runs");
    let run2 = runs
        .iter()
        .find(|run| run.run_accession == "RUN2")
        .expect("missing RUN2");
    assert_eq!(run2.state, RunState::QuantFailed);
    assert!(
        run2.last_error
            .as_deref()
            .unwrap_or_default()
            .contains("external abundance declaration")
    );
}

#[test]
fn normalize_rejects_confounded_design_before_advancing_runs() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        None,
    );
    let (paths, project_id) = setup_project(root, true);

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

    let db = Database::open(root).expect("failed to reopen db");
    let runs = db.list_runs().expect("failed to list runs");
    assert!(runs.iter().all(|run| run.state == RunState::QuantDone));
    assert!(
        !paths
            .de_project_dir(&project_id)
            .join("deseq2_manifest.json")
            .exists()
    );
}

#[test]
fn normalize_fails_when_design_references_missing_column() {
    let temp = TempDir::new().expect("failed to create tempdir");
    let root = temp.path();
    run_ok(
        &["init", "--root", root.to_str().expect("invalid root path")],
        None,
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
        None,
    );
    assert!(!output.status.success(), "normalize unexpectedly succeeded");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("missing samplesheet columns: batch"));

    let db = Database::open(root).expect("failed to reopen db");
    let runs = db.list_runs().expect("failed to list runs");
    assert!(runs.iter().all(|run| run.state == RunState::QuantDone));
}
