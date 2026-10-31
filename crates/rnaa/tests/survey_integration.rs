mod support;
use rnaa_core::{
    ProjectConfig,
    metadata::{ColumnSpec, ColumnType},
    traits::ReferenceManager,
    util::compute_sha256,
};
use std::{fs, path::PathBuf};
use support::native::{command, setup, success};

#[test]
fn survey_rejects_changed_correlation_covariates() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, _db) = setup(root);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config
        .metadata
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    config.save(&paths.config_path).unwrap();
    let write_sheet = |ages: &[u32; 8]| {
        let mut sheet = String::from("run_accession\tcondition\tage\n");
        for (i, age) in ages.iter().enumerate() {
            sheet.push_str(&format!(
                "RUN{}\t{}\t{age}\n",
                i + 1,
                if i < 4 { "A" } else { "B" }
            ));
        }
        fs::write(paths.samplesheet_path(), sheet).unwrap();
    };
    write_sheet(&[20, 22, 21, 23, 21, 23, 20, 22]);
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
    write_sheet(&[34, 20, 29, 22, 25, 31, 23, 28]);
    success(
        root,
        &[
            "corr",
            "--adjust",
            "residualize",
            "--model",
            "~ age + condition",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=5",
        ],
    );
    let output = command(root, &["survey"]);
    assert!(!output.status.success());
    let error = String::from_utf8_lossy(&output.stderr);
    assert!(
        error.contains("covariate") && error.contains("age"),
        "{error}"
    );
    assert!(
        !paths.exports_dir.exists() || fs::read_dir(&paths.exports_dir).unwrap().next().is_none()
    );
}

#[test]
fn survey_requires_matching_frozen_generations_and_summarizes_biological_samples() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize", "--design", "~ condition"]);
    success(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=5",
        ],
    );
    let incomplete = command(root, &["survey"]);
    assert!(!incomplete.status.success());
    assert!(
        String::from_utf8_lossy(&incomplete.stderr)
            .contains("requires differential-expression results")
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
    let corr_head = db
        .generation_head(rnaa_core::generations::OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let de_head = db
        .generation_head(rnaa_core::generations::OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let output = command(
        root,
        &[
            "survey",
            "--stratify",
            "condition",
            "--require-current-admission",
        ],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let directory = PathBuf::from(
        String::from_utf8(output.stdout)
            .unwrap()
            .lines()
            .find_map(|line| line.strip_prefix("survey\t"))
            .unwrap(),
    );
    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(directory.join("metadata_summary.json")).unwrap())
            .unwrap();
    assert_eq!(summary["sample_count"], 8);
    assert_eq!(summary["gene_count"], 500);
    assert_eq!(summary["correlation_gene_count"], 40);
    assert_eq!(summary["stratification"]["column"], "condition");
    assert_eq!(summary["stratification"]["levels"]["A"], 4);
    assert_eq!(summary["stratification"]["levels"]["B"], 4);
    assert_eq!(summary["stratification"]["missing"], 0);
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(directory.join("survey-manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["kind"], "survey");
    assert_eq!(
        manifest["current_cohort_validation"]["membership_sha256"]
            .as_str()
            .unwrap()
            .len(),
        64
    );
    assert_eq!(
        manifest["current_admission_validation"]["full_freshness"],
        false
    );
    assert_eq!(manifest["generations"].as_array().unwrap().len(), 2);
    assert!(
        manifest["generations"]
            .as_array()
            .unwrap()
            .iter()
            .any(|head| head["generation_id"] == corr_head.generation_id)
    );
    assert_ne!(
        manifest["survey_compatibility"]["correlation_input_generation"],
        de_head.generation_id
    );
    let current = command(root, &["survey", "--require-current"]);
    assert!(
        current.status.success(),
        "{}",
        String::from_utf8_lossy(&current.stderr)
    );
    assert!(
        manifest["files"]
            .as_array()
            .unwrap()
            .iter()
            .any(|entry| entry["output"] == "metadata_summary.json"
                && entry["output_sha256"]
                    == compute_sha256(&directory.join("metadata_summary.json")).unwrap())
    );
    assert!(directory.join("de/gene_counts.tsv").exists());
    assert!(directory.join("de/gene_norm_counts.tsv").exists());
    assert!(directory.join("de/vst.tsv").exists());
    assert!(directory.join("corr/edges_topk.tsv").exists());
    let prior = compute_sha256(&directory.join("survey-manifest.json")).unwrap();
    success(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:40",
            "--out",
            "dense",
        ],
    );
    let dense = command(root, &["survey"]);
    assert!(
        dense.status.success(),
        "{}",
        String::from_utf8_lossy(&dense.stderr)
    );
    let dense_dir = PathBuf::from(
        String::from_utf8(dense.stdout)
            .unwrap()
            .lines()
            .find_map(|line| line.strip_prefix("survey\t"))
            .unwrap(),
    );
    let dense_summary: serde_json::Value =
        serde_json::from_slice(&fs::read(dense_dir.join("metadata_summary.json")).unwrap())
            .unwrap();
    assert_eq!(dense_summary["correlation_dense_tables"], 1);
    let count = fs::read_dir(&paths.exports_dir).unwrap().count();
    let invalid = command(root, &["survey", "--stratify", "missing-column"]);
    assert!(!invalid.status.success());
    assert!(String::from_utf8_lossy(&invalid.stderr).contains("stratification column"));
    assert_eq!(fs::read_dir(&paths.exports_dir).unwrap().count(), count);
    assert_eq!(
        compute_sha256(&directory.join("survey-manifest.json")).unwrap(),
        prior
    );
    let abundance = paths
        .quant_run_dir("RUN1", "kallistors")
        .join("abundance.tsv");
    let text = fs::read_to_string(&abundance).unwrap();
    let changed = text.replacen("\t100\t", "\t101\t", 1);
    assert_ne!(text, changed);
    fs::write(&abundance, changed).unwrap();
    let config = ProjectConfig::load(&paths.config_path).unwrap();
    let reference = rnaa_stages::refs::EnsemblReferenceManager::new()
        .unwrap()
        .ensure_reference(&paths, &db, &config)
        .unwrap();
    rnaa_stages::quant::declare_external_quant(
        "RUN1",
        abundance.parent().unwrap(),
        &reference,
        "revised survey test fixture",
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
    let stale = command(root, &["survey"]);
    assert!(!stale.status.success());
    let reason = String::from_utf8_lossy(&stale.stderr);
    assert!(
        reason.contains("does not match the selected analysis VST matrix")
            || reason.contains("generation failed verification")
            || reason.contains("eligibility identities differ"),
        "{reason}"
    );
    assert_eq!(fs::read_dir(&paths.exports_dir).unwrap().count(), count);
    assert_eq!(
        compute_sha256(&directory.join("survey-manifest.json")).unwrap(),
        prior
    );
}

#[test]
fn survey_accepts_real_module_and_enrichment_outputs() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
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
    let selected = root.join("genes.txt");
    fs::write(
        &selected,
        (0..100).map(|i| format!("G{i:04}\n")).collect::<String>(),
    )
    .unwrap();
    let gmt = root.join("pathways.gmt");
    fs::write(
        &gmt,
        format!(
            "UP\tsynthetic\t{}\nDOWN\tsynthetic\t{}\n",
            (0..50)
                .map(|i| format!("G{i:04}"))
                .collect::<Vec<_>>()
                .join("\t"),
            (50..100)
                .map(|i| format!("G{i:04}"))
                .collect::<Vec<_>>()
                .join("\t")
        ),
    )
    .unwrap();
    let geneset = format!("genes:{}", selected.display());
    success(
        root,
        &[
            "corr",
            "--method",
            "pearson",
            "--adjust",
            "none",
            "--geneset",
            &geneset,
            "--out",
            "edges:topk=5",
            "--module-gmt",
            gmt.to_str().unwrap(),
            "--module-min-size",
            "2",
            "--fgsea-permutations",
            "100",
        ],
    );
    let output = command(root, &["survey"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let directory = PathBuf::from(
        String::from_utf8(output.stdout)
            .unwrap()
            .lines()
            .find_map(|line| line.strip_prefix("survey\t"))
            .unwrap(),
    );
    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(directory.join("metadata_summary.json")).unwrap())
            .unwrap();
    assert!(summary["module_tables"].as_u64().unwrap() >= 3);
    assert!(summary["enrichment_tables"].as_u64().unwrap() >= 1);
    assert_eq!(summary["correlation_gene_count"], 100);
    let strict = command(root, &["survey", "--require-current"]);
    assert!(
        strict.status.success(),
        "{}",
        String::from_utf8_lossy(&strict.stderr)
    );
    assert!(
        command(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .status
        .success()
    );
    let original_gmt = fs::read(&gmt).unwrap();
    fs::write(&gmt, b"changed\tsynthetic\tG0001\tG0002\n").unwrap();
    let changed = command(root, &["survey", "--require-current"]);
    assert!(!changed.status.success());
    assert!(
        String::from_utf8_lossy(&changed.stderr)
            .contains("reconstructed current scientific request")
    );
    fs::write(&gmt, original_gmt).unwrap();

    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    assert_eq!(
        config.corr.modules.as_ref().unwrap().gmt_path,
        gmt.canonicalize().unwrap()
    );
    let original_modules = config.corr.modules.clone();
    config.corr.modules.as_mut().unwrap().seed += 1;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let changed = command(root, &["survey", "--require-current"]);
    assert!(!changed.status.success());
    assert!(
        String::from_utf8_lossy(&changed.stderr)
            .contains("reconstructed current scientific request")
    );
    config.corr.modules = None;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let legacy = command(root, &["survey", "--require-current"]);
    assert!(!legacy.status.success());
    assert!(
        String::from_utf8_lossy(&legacy.stderr)
            .contains("older results may need correlation to be rerun")
    );
    config.corr.modules = original_modules;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    success(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            &geneset,
            "--out",
            "edges:topk=5",
        ],
    );
    assert!(
        ProjectConfig::load(&paths.config_path)
            .unwrap()
            .corr
            .modules
            .is_none()
    );
    assert!(
        command(root, &["survey", "--require-current"])
            .status
            .success()
    );
    assert!(directory.join("survey-manifest.json").exists());
}
