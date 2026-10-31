//! Successful native analysis from deterministic synthetic abundance inputs.
mod support;
use rnaa_core::{
    ProjectConfig, RunState,
    traits::ReferenceManager,
    util::{compute_sha256, now_rfc3339},
};
use std::{fs, path::Path};
use support::native::{command, gzip, setup, success};

#[test]
fn native_de_and_correlation_reuse_and_invalidate_real_outputs() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let current_paths = || paths.with_generation_heads(&db).unwrap();
    success(
        root,
        &[
            "deseq2",
            "--design",
            "~ condition",
            "--contrast",
            "condition",
            "B",
            "A",
        ],
    );
    let de = current_paths().de_project_dir(db.project_id());
    let manifest = || {
        current_paths()
            .de_project_dir(db.project_id())
            .join("deseq2_manifest.json")
    };
    let first = compute_sha256(&manifest()).unwrap();
    let mut results = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(de.join("de_condition_B_vs_A.tsv"))
        .unwrap();
    let rows = results
        .records()
        .map(|row| row.unwrap())
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 500);
    assert!(
        rows[..50]
            .iter()
            .filter(|row| row[2].parse::<f64>().unwrap_or_default() > 0.0)
            .count()
            >= 45
    );
    assert!(
        rows[50..100]
            .iter()
            .filter(|row| row[2].parse::<f64>().unwrap_or_default() < 0.0)
            .count()
            >= 45
    );
    success(
        root,
        &[
            "deseq2",
            "--design",
            "~ condition",
            "--contrast",
            "condition",
            "B",
            "A",
        ],
    );
    assert_eq!(compute_sha256(&manifest()).unwrap(), first);
    success(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:100",
            "--out",
            "edges:topk=5",
        ],
    );
    let corr_manifest = || current_paths().manifest_path("corr", db.project_id());
    let first_corr = compute_sha256(&corr_manifest()).unwrap();
    let selected = rnaa_formats::matrix::read_matrix_tsv(
        &current_paths().corr_dir.join("selected_matrix.tsv"),
    )
    .unwrap();
    let adjusted = rnaa_formats::matrix::read_matrix_tsv(
        &current_paths().corr_dir.join("adjusted_matrix.tsv"),
    )
    .unwrap();
    assert_eq!(
        selected.values, adjusted.values,
        "none must preserve values"
    );
    success(root, &["corr"]);
    assert_eq!(
        compute_sha256(&corr_manifest()).unwrap(),
        first_corr,
        "identical request must reuse correlation"
    );
    success(root, &["corr", "--adjust", "residualize"]);
    assert_ne!(compute_sha256(&corr_manifest()).unwrap(), first_corr);
    let request_before = compute_sha256(&manifest()).unwrap();
    let abundance = paths
        .quant_run_dir("RUN1", "kallistors")
        .join("abundance.tsv");
    let text = fs::read_to_string(&abundance).unwrap();
    fs::write(&abundance, text.replacen("\t100\t", "\t101\t", 1)).unwrap();
    let config = ProjectConfig::load(&paths.config_path).unwrap();
    let reference = rnaa_stages::refs::EnsemblReferenceManager::new()
        .unwrap()
        .ensure_reference(&paths, &db, &config)
        .unwrap();
    rnaa_stages::quant::declare_external_quant(
        "RUN1",
        abundance.parent().unwrap(),
        &reference,
        "explicitly revised simulated abundance fixture",
    )
    .unwrap();
    success(
        root,
        &[
            "deseq2",
            "--design",
            "~ condition",
            "--contrast",
            "condition",
            "B",
            "A",
        ],
    );
    assert_ne!(
        compute_sha256(&manifest()).unwrap(),
        request_before,
        "changed abundance content must recompute DE"
    );
    let edge_path = || {
        current_paths()
            .corr_project_dir(db.project_id())
            .join("edges_topk.tsv")
    };
    fs::write(edge_path(), b"corrupt output").unwrap();
    success(root, &["corr"]);
    assert!(
        fs::read_to_string(edge_path())
            .unwrap()
            .starts_with("gene_a\tgene_b\tr")
    );
    let before_metadata = compute_sha256(&corr_manifest()).unwrap();
    let sheet = fs::read_to_string(paths.samplesheet_path()).unwrap();
    fs::write(
        paths.samplesheet_path(),
        sheet.replace("RUN1\tA", "RUN1\tC"),
    )
    .unwrap();
    let stale = command(root, &["corr"]);
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("DE metadata is stale"));
    success(root, &["deseq2"]);
    let updated: serde_json::Value = serde_json::from_slice(
        &fs::read(
            current_paths()
                .de_project_dir(db.project_id())
                .join("rnaa_deseq2_request.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(
        updated["request"]["parameters"]["design_contract"]["columns"]["condition"]["values"][0],
        "C"
    );
    success(root, &["corr"]);
    assert_ne!(compute_sha256(&corr_manifest()).unwrap(), before_metadata);
    fs::write(paths.samplesheet_path(), sheet).unwrap();
    success(root, &["deseq2"]);
    let genes = root.join("selection.txt");
    fs::write(
        &genes,
        (0..100).map(|i| format!("G{i:04}\n")).collect::<String>(),
    )
    .unwrap();
    let genes_spec = format!("genes:{}", genes.display());
    success(root, &["corr", "--geneset", &genes_spec]);
    let before_selection = compute_sha256(&corr_manifest()).unwrap();
    fs::write(
        &genes,
        (10..110).map(|i| format!("G{i:04}\n")).collect::<String>(),
    )
    .unwrap();
    success(root, &["corr"]);
    assert_ne!(compute_sha256(&corr_manifest()).unwrap(), before_selection);
    let gmt = root.join("pathways.gmt");
    let pathways = format!(
        "UP\tsynthetic\t{}\nDOWN\tsynthetic\t{}\n",
        (0..50)
            .map(|i| format!("G{i:04}"))
            .collect::<Vec<_>>()
            .join("\t"),
        (50..100)
            .map(|i| format!("G{i:04}"))
            .collect::<Vec<_>>()
            .join("\t")
    );
    fs::write(&gmt, &pathways).unwrap();
    let module_args = [
        "corr",
        "--module-gmt",
        gmt.to_str().unwrap(),
        "--module-min-size",
        "2",
        "--module-top",
        "50",
        "--fgsea-permutations",
        "100",
    ];
    success(root, &module_args);
    let module_generation = compute_sha256(&corr_manifest()).unwrap();
    assert!(
        fs::read_to_string(
            current_paths()
                .corr_project_dir(db.project_id())
                .join("module_fgsea.tsv")
        )
        .unwrap()
        .lines()
        .count()
            > 1
    );
    success(root, &module_args);
    assert_eq!(compute_sha256(&corr_manifest()).unwrap(), module_generation);
    fs::write(&gmt, pathways.replace("UP\t", "RENAMED\t")).unwrap();
    success(root, &module_args);
    assert_ne!(compute_sha256(&corr_manifest()).unwrap(), module_generation);
    let jobs = db.list_job_attempts(100).unwrap();
    assert!(!jobs.is_empty());
    assert!(
        jobs.iter()
            .filter(|job| job.state == "SUCCEEDED")
            .all(|job| job.request_sha256.is_some()),
        "successful computations and reuses must both retain the scientific request identity: {jobs:?}"
    );
    let failures = jobs
        .iter()
        .filter(|job| job.state != "SUCCEEDED")
        .collect::<Vec<_>>();
    assert_eq!(
        failures.len(),
        1,
        "only the intentional stale-metadata rejection may fail"
    );
    assert_eq!(failures[0].stage, "corr");
    assert_eq!(failures[0].state, "FAILED");
    assert!(
        failures[0]
            .last_error
            .as_deref()
            .unwrap_or_default()
            .contains("DE metadata is stale")
    );
}

#[test]
fn failed_replacements_keep_the_previous_native_analysis_generation() {
    use rnaa_core::generations::OutputGroup;
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let de_args = [
        "deseq2",
        "--design",
        "~ condition",
        "--contrast",
        "condition",
        "B",
        "A",
    ];
    success(root, &de_args);
    let de = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let de_manifest = paths
        .with_generation_heads(&db)
        .unwrap()
        .de_project_dir(db.project_id())
        .join("deseq2_manifest.json");
    let de_digest = compute_sha256(&de_manifest).unwrap();
    // Contrast assembly fails after the engine has already written its matrices.
    let failure = command(
        root,
        &["deseq2", "--contrast", "condition", "MISSING_LEVEL", "A"],
    );
    assert!(!failure.status.success());
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        de.generation_id
    );
    assert_eq!(compute_sha256(&de_manifest).unwrap(), de_digest);
    assert!(
        db.list_artifacts()
            .unwrap()
            .iter()
            .filter(|a| a.run_accession.is_none() && a.path.contains("/payload/de/"))
            .all(|a| Path::new(&a.path).starts_with(&de.payload_path))
    );
    success(root, &de_args);
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        de.generation_id,
        "retrying the original request reuses the prior successful generation"
    );

    success(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:100",
            "--out",
            "edges:topk=5",
        ],
    );
    let corr = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let corr_manifest = paths
        .with_generation_heads(&db)
        .unwrap()
        .manifest_path("corr", db.project_id());
    let digest = compute_sha256(&corr_manifest).unwrap();
    // Gene selection has written its new matrix before adjustment rejects this.
    assert!(
        !command(root, &["corr", "--adjust", "invalid"])
            .status
            .success()
    );
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        corr.generation_id
    );
    assert_eq!(compute_sha256(&corr_manifest).unwrap(), digest);
    success(root, &["corr", "--adjust", "none"]);
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        corr.generation_id
    );
    let vst = paths
        .with_generation_heads(&db)
        .unwrap()
        .de_project_dir(db.project_id())
        .join("vst.tsv");
    let mut damaged = rnaa_formats::matrix::read_matrix_tsv(&vst).unwrap();
    damaged.values[0][0] += 1.0;
    rnaa_formats::matrix::write_matrix_tsv(&vst, &damaged).unwrap();
    assert!(
        !db.verify_generation(OutputGroup::De, db.project_id())
            .unwrap()
    );
    let rejected = command(root, &["corr", "--adjust", "none"]);
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("DE generation failed integrity"));
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        corr.generation_id
    );
    success(root, &de_args);
    assert!(
        db.verify_generation(OutputGroup::De, db.project_id())
            .unwrap()
    );
    assert_ne!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        de.generation_id
    );
    success(root, &["corr", "--adjust", "none"]);
}

#[test]
fn native_normalization_publishes_and_reuses_its_own_complete_generation() {
    use rnaa_core::generations::OutputGroup;
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize", "--design", "~ condition"]);
    let head = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    assert_eq!(head.stage, "normalize");
    assert!(
        db.verify_generation(OutputGroup::De, db.project_id())
            .unwrap()
    );
    let current = paths.with_generation_heads(&db).unwrap();
    let transformed = rnaa_formats::matrix::read_matrix_tsv(
        &current.de_project_dir(db.project_id()).join("vst.tsv"),
    )
    .unwrap();
    assert_eq!(
        (transformed.row_ids.len(), transformed.col_ids.len()),
        (500, 8)
    );
    assert!(
        db.list_artifacts()
            .unwrap()
            .iter()
            .any(|artifact| artifact.kind == rnaa_core::ArtifactKind::SizeFactors)
    );
    success(root, &["normalize"]);
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        head.generation_id
    );
    success(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:100",
            "--out",
            "edges:topk=5",
        ],
    );
    assert!(
        db.verify_generation(OutputGroup::Corr, db.project_id())
            .unwrap()
    );
}

#[test]
fn analysis_excludes_rejected_runs_and_refuses_a_matrix_from_an_old_qc_cohort() {
    use rnaa_core::{
        ArtifactKind, ArtifactRecord, eligibility::EligibilityState, generations::OutputGroup,
    };
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.preprocess = true;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    for run in db.list_runs().unwrap() {
        let report = root.join(format!("{}_qc.json", run.run_accession));
        let passed = if run.run_accession == "RUN1" { 80 } else { 100 };
        fs::write(&report,serde_json::to_vec(&serde_json::json!({"mode":"single","passed_reads":passed,"failed_reads":100-passed})).unwrap()).unwrap();
        db.record_artifact(&ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some(run.run_accession),
            kind: ArtifactKind::PreprocessReport,
            path: report.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".into(),
            checksum: compute_sha256(&report).unwrap(),
            bytes: fs::metadata(&report).unwrap().len(),
            created_at: now_rfc3339(),
        })
        .unwrap();
    }
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let de = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let corr = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let old_matrix = paths
        .with_generation_heads(&db)
        .unwrap()
        .de_project_dir(db.project_id())
        .join("vst.tsv");
    let abundance = paths
        .quant_run_dir("RUN1", "kallistors")
        .join("abundance.tsv");
    let counts = fs::read(&abundance).unwrap();
    config.quant.qc_gate.min_pass_rate = Some(0.95);
    config.save(&paths.config_path).unwrap();
    success(root, &["qc", "reevaluate"]);
    assert_eq!(
        db.eligibility("RUN1").unwrap().unwrap().descriptor.state,
        EligibilityState::Rejected
    );
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantDone
    );
    let result = command(root, &["corr"]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("QC policy is stale"));
    let external_matrix = root.join("externally_supplied_old_matrix.tsv");
    fs::copy(&old_matrix, &external_matrix).unwrap();
    let custom = command(
        root,
        &["corr", "--matrix", external_matrix.to_str().unwrap()],
    );
    assert!(!custom.status.success());
    assert!(String::from_utf8_lossy(&custom.stderr).contains("QC-ineligible run RUN1"));
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap(),
        de
    );
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap(),
        corr
    );
    success(root, &["normalize"]);
    let matrix = rnaa_formats::matrix::read_matrix_tsv(
        &paths
            .with_generation_heads(&db)
            .unwrap()
            .de_project_dir(db.project_id())
            .join("vst.tsv"),
    )
    .unwrap();
    assert_eq!(matrix.col_ids.len(), 7);
    assert!(!matrix.col_ids.contains(&"BIO0".to_string()));
    // --matrix was persisted by neither failed invocation; the default is current DE.
    success(root, &["corr"]);
    assert_eq!(fs::read(&abundance).unwrap(), counts);
    config.quant.qc_gate.min_pass_rate = Some(0.6);
    config.save(&paths.config_path).unwrap();
    let stale = command(root, &["corr"]);
    assert!(!stale.status.success());
    success(root, &["normalize"]);
    let matrix = rnaa_formats::matrix::read_matrix_tsv(
        &paths
            .with_generation_heads(&db)
            .unwrap()
            .de_project_dir(db.project_id())
            .join("vst.tsv"),
    )
    .unwrap();
    assert_eq!(matrix.col_ids.len(), 8);
    assert_eq!(fs::read(&abundance).unwrap(), counts);
    let before = db.get_run("RUN1").unwrap().unwrap().state;
    fs::write(root.join("RUN1_qc.json"), "corrupt report").unwrap();
    success(root, &["normalize"]);
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        before,
        "invalid QC evidence must not demote completed counts during analysis validation"
    );
    assert_eq!(
        db.eligibility("RUN1").unwrap().unwrap().descriptor.state,
        EligibilityState::InvalidEvidence
    );
    assert_eq!(fs::read(&abundance).unwrap(), counts);
}

#[test]
fn changed_annotations_refresh_analysis_and_changed_reference_rejects_supplied_counts() {
    use rnaa_core::{ProjectStageState, generations::OutputGroup};
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let first_de = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let first_corr = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let first_ref = db.selected_reference().unwrap().unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    let old_gtf = fs::File::open(&config.refs.custom_gtf).unwrap();
    let mut decoded = String::new();
    std::io::Read::read_to_string(
        &mut flate2::read::MultiGzDecoder::new(old_gtf),
        &mut decoded,
    )
    .unwrap();
    gzip(
        Path::new(&config.refs.custom_gtf),
        &decoded.replace("SymbolG", "RenamedG"),
    );
    success(root, &["refs", "prepare"]);
    let second_ref = db.selected_reference().unwrap().unwrap();
    assert_eq!(first_ref.index.checksum, second_ref.index.checksum);
    assert_ne!(
        first_ref.annotation.checksum,
        second_ref.annotation.checksum
    );
    assert_eq!(
        db.get_project_stage_state("corr").unwrap().unwrap().state,
        ProjectStageState::Pending
    );
    let stale = command(root, &["corr"]);
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("DE matrix reference is stale"));
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap(),
        first_corr
    );
    success(root, &["normalize"]);
    assert_ne!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap(),
        first_de
    );
    assert!(
        db.list_runs()
            .unwrap()
            .iter()
            .all(|run| run.state == RunState::QuantDone)
    );
    success(root, &["corr"]);
    let current_de = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let count_path = paths
        .quant_run_dir("RUN1", "kallistors")
        .join("abundance.tsv");
    let data = fs::read(&count_path).unwrap();
    config.refs.index_k = 29;
    config.save(&paths.config_path).unwrap();
    success(root, &["refs", "prepare"]);
    let stale = command(root, &["corr"]);
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("DE matrix reference is stale"));
    let incompatible = command(root, &["normalize"]);
    assert!(!incompatible.status.success());
    assert!(
        db.list_runs()
            .unwrap()
            .iter()
            .all(|run| run.state == RunState::QuantFailed)
    );
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap(),
        current_de
    );
    assert_eq!(fs::read(count_path).unwrap(), data);
}
