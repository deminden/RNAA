use super::*;

#[test]
fn optional_current_cohort_check_rejects_old_results_without_claiming_full_freshness() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let current = exported(
        root,
        &["export", "--what", "all", "--require-current-cohort"],
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(current.join("export-manifest.json")).unwrap()).unwrap();
    assert_eq!(
        manifest["current_cohort_validation"]["membership_sha256"]
            .as_str()
            .unwrap()
            .len(),
        64
    );
    assert_eq!(
        manifest["current_cohort_validation"]["full_freshness"],
        false
    );
    assert!(
        manifest["selection"]
            .as_str()
            .unwrap()
            .contains("current scientific settings are not revalidated")
    );
    let retained = compute_sha256(&current.join("export-manifest.json")).unwrap();
    let config_bytes = fs::read(&paths.config_path).unwrap();
    let mut changed_config = ProjectConfig::load(&paths.config_path).unwrap();
    changed_config.resolver.fields.push("fixture_extra".into());
    changed_config.save(&paths.config_path).unwrap();
    let unresolved = command(
        root,
        &["export", "--what", "counts", "--require-current-cohort"],
    );
    assert!(!unresolved.status.success());
    assert!(String::from_utf8_lossy(&unresolved.stderr).contains("resolver configuration changed"));
    fs::write(&paths.config_path, config_bytes).unwrap();
    let mut annotated = db.get_run("RUN1").unwrap().unwrap();
    annotated
        .metadata
        .insert("operational_note".into(), "not a membership change".into());
    db.upsert_run(&annotated).unwrap();
    assert!(
        exported(
            root,
            &["export", "--what", "counts", "--require-current-cohort"]
        )
        .exists()
    );
    let runs = db.list_runs().unwrap();
    let selected = runs
        .into_iter()
        .filter(|run| !["RUN1", "RUN8"].contains(&run.run_accession.as_str()))
        .collect::<Vec<_>>();
    db.publish_resolved_cohort(&selected, &db.selector_contract().unwrap(), None)
        .unwrap();
    let count = fs::read_dir(&paths.exports_dir).unwrap().count();
    let stale = command(
        root,
        &["export", "--what", "all", "--require-current-cohort"],
    );
    assert!(!stale.status.success());
    assert!(
        String::from_utf8_lossy(&stale.stderr).contains("does not match the current active cohort")
    );
    assert_eq!(fs::read_dir(&paths.exports_dir).unwrap().count(), count);
    assert_eq!(
        compute_sha256(&current.join("export-manifest.json")).unwrap(),
        retained
    );
    let frozen = exported(root, &["export", "--what", "counts"]);
    assert!(frozen.join("de/gene_counts.tsv").exists());
    success(root, &["normalize"]);
    let mixed = command(
        root,
        &["export", "--what", "all", "--require-current-cohort"],
    );
    assert!(!mixed.status.success());
    assert!(String::from_utf8_lossy(&mixed.stderr).contains("corr generation does not match"));
    let updated = exported(
        root,
        &["export", "--what", "counts", "--require-current-cohort"],
    );
    assert!(updated.join("de/gene_counts.tsv").exists());
}

#[test]
fn current_admission_check_rejects_changed_counts_reference_and_qc_without_refreshing_them() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let accepted = exported(
        root,
        &["export", "--what", "all", "--require-current-admission"],
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(accepted.join("export-manifest.json")).unwrap()).unwrap();
    assert_eq!(
        manifest["current_admission_validation"]["full_freshness"],
        false
    );
    assert_eq!(
        manifest["current_admission_validation"]["eligibility_sha256"]
            .as_str()
            .unwrap()
            .len(),
        64
    );
    assert_eq!(
        manifest["current_admission_validation"]["reference_sha256"]
            .as_str()
            .unwrap()
            .len(),
        64
    );
    assert_eq!(
        manifest["current_cohort_validation"]["full_freshness"],
        false
    );
    let count = fs::read_dir(&paths.exports_dir).unwrap().count();
    let abundance = paths
        .for_quant_run(&db, "RUN1")
        .unwrap()
        .quant_run_dir("RUN1", "kallistors")
        .join("abundance.tsv");
    let original = fs::read(&abundance).unwrap();
    let mut changed = original.clone();
    changed.push(b'\n');
    fs::write(&abundance, changed).unwrap();
    let stale = command(
        root,
        &["export", "--what", "all", "--require-current-admission"],
    );
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("quantified abundance inputs"));
    assert_eq!(fs::read_dir(&paths.exports_dir).unwrap().count(), count);
    let cohort_only = exported(
        root,
        &["export", "--what", "counts", "--require-current-cohort"],
    );
    assert!(cohort_only.exists());
    fs::write(&abundance, original).unwrap();

    let config = ProjectConfig::load(&paths.config_path).unwrap();
    let cdna = PathBuf::from(&config.refs.custom_cdna);
    let original_reference = fs::read(&cdna).unwrap();
    let mut modified_reference = original_reference.clone();
    modified_reference.push(0);
    fs::write(&cdna, modified_reference).unwrap();
    let stale = command(
        root,
        &["export", "--what", "counts", "--require-current-admission"],
    );
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("reference source changed"));
    fs::write(&cdna, original_reference).unwrap();

    let prior = serde_json::to_value(db.eligibility("RUN1").unwrap()).unwrap();
    let original_config = fs::read(&paths.config_path).unwrap();
    let mut changed_config = config;
    changed_config.quant.preprocess = true;
    changed_config.save(&paths.config_path).unwrap();
    let stale = command(
        root,
        &["export", "--what", "counts", "--require-current-admission"],
    );
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("QC eligibility is stale"));
    assert_eq!(
        serde_json::to_value(db.eligibility("RUN1").unwrap()).unwrap(),
        prior
    );
    fs::write(&paths.config_path, original_config).unwrap();
    assert!(
        exported(
            root,
            &["export", "--what", "counts", "--require-current-admission"]
        )
        .exists()
    );
}

#[test]
fn current_design_check_uses_only_model_metadata_and_keeps_its_scope_explicit() {
    use rnaa_core::metadata::{ColumnSpec, ColumnType};
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config
        .metadata
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let write_sheet = |ages: &[u32; 8], swap_condition: bool, note: &str| {
        let mut sheet = String::from("run_accession\tcondition\tage\tnote\n");
        for (index, age) in ages.iter().enumerate() {
            let mut condition = if index < 4 { "A" } else { "B" };
            if swap_condition && (index == 0 || index == 4) {
                condition = if condition == "A" { "B" } else { "A" };
            }
            sheet.push_str(&format!("RUN{}\t{condition}\t{age}\t{note}\n", index + 1));
        }
        fs::write(paths.samplesheet_path(), sheet).unwrap();
    };
    let ages = [20, 22, 21, 23, 21, 23, 20, 22];
    write_sheet(&ages, false, "first");
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
    let accepted = exported(
        root,
        &["export", "--what", "all", "--require-current-design"],
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(accepted.join("export-manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["current_design_validation"]["de_checked"], true);
    assert_eq!(
        manifest["current_design_validation"]["correlation_checked"],
        true
    );
    assert_eq!(
        manifest["current_design_validation"]["full_freshness"],
        false
    );
    assert!(
        command(root, &["survey", "--require-current-design"])
            .status
            .success()
    );
    write_sheet(&ages, false, "unrelated change");
    assert!(
        exported(
            root,
            &["export", "--what", "all", "--require-current-design"]
        )
        .exists()
    );
    let changed_ages = [34, 20, 29, 22, 25, 31, 23, 28];
    write_sheet(&changed_ages, false, "unrelated change");
    assert!(
        exported(
            root,
            &["export", "--what", "counts", "--require-current-design"]
        )
        .exists()
    );
    let stale_corr = command(
        root,
        &[
            "export",
            "--what",
            "correlation",
            "--require-current-design",
        ],
    );
    assert!(!stale_corr.status.success());
    assert!(
        String::from_utf8_lossy(&stale_corr.stderr)
            .contains("correlation generation does not match current design")
    );
    write_sheet(&ages, true, "unrelated change");
    let stale_de = command(
        root,
        &["export", "--what", "counts", "--require-current-design"],
    );
    assert!(!stale_de.status.success());
    assert!(
        String::from_utf8_lossy(&stale_de.stderr)
            .contains("DE generation does not match current design")
    );
    write_sheet(&ages, false, "unrelated change");
    let mut changed_config = ProjectConfig::load(&paths.config_path).unwrap();
    changed_config.corr.method = "kendall".into();
    changed_config.save(&paths.config_path).unwrap();
    db.set_project_config(&changed_config).unwrap();
    assert!(
        exported(
            root,
            &["export", "--what", "all", "--require-current-design"]
        )
        .exists(),
        "correlation method is outside the model-design check"
    );
    changed_config.deseq2.design = "~ 1".into();
    changed_config.save(&paths.config_path).unwrap();
    db.set_project_config(&changed_config).unwrap();
    let stale_design = command(
        root,
        &["export", "--what", "counts", "--require-current-design"],
    );
    assert!(!stale_design.status.success());
    assert!(
        String::from_utf8_lossy(&stale_design.stderr)
            .contains("DE generation does not match current design")
    );
}

#[test]
fn current_request_check_reconstructs_scientific_settings_and_selected_inputs() {
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
    let genes = root.join("selected-genes.txt");
    fs::write(
        &genes,
        (0..60).map(|i| format!("G{i:04}\n")).collect::<String>(),
    )
    .unwrap();
    let geneset = format!("genes:{}", genes.display());
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
    let accepted = exported(root, &["export", "--what", "all", "--require-current"]);
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(accepted.join("export-manifest.json")).unwrap()).unwrap();
    assert_eq!(
        manifest["current_request_validation"]["requests"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        manifest["current_request_validation"]["external_files_locked"],
        false
    );
    assert!(
        command(root, &["survey", "--require-current"])
            .status
            .success()
    );
    assert!(
        exported(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .exists()
    );

    let original_genes = fs::read(&genes).unwrap();
    fs::write(&genes, b"G0001\nG0002\n").unwrap();
    let stale = command(root, &["export", "--what", "all", "--require-current"]);
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains(
        "correlation generation does not match the reconstructed current scientific request"
    ));
    fs::write(&genes, original_genes).unwrap();

    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.corr.method = "kendall".into();
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let stale = command(root, &["export", "--what", "all", "--require-current"]);
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains(
        "correlation generation does not match the reconstructed current scientific request"
    ));
    config.corr.method = "pearson".into();
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();

    let original_contrasts = db.list_contrasts().unwrap();
    db.replace_contrasts(&[rnaa_core::ContrastSpec::new(
        None,
        "condition".into(),
        "A".into(),
        "B".into(),
    )])
    .unwrap();
    let stale = command(root, &["export", "--what", "de", "--require-current"]);
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("current contrast selection"));
    assert!(
        exported(
            root,
            &["export", "--what", "de", "--require-current-design"]
        )
        .exists()
    );
    db.replace_contrasts(&original_contrasts).unwrap();
    let current_paths = paths.with_generation_heads(&db).unwrap();
    let vst = current_paths
        .de_project_dir(db.project_id())
        .join("vst.tsv");
    success(
        root,
        &[
            "corr",
            "--matrix",
            vst.to_str().unwrap(),
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
            .matrix
            .is_none()
    );
    assert!(
        exported(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .exists()
    );
    let external = root.join("external-matrix.tsv");
    let mut input = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(vst)
        .unwrap();
    let headers = input.byte_headers().unwrap().clone();
    let mut output = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&external)
        .unwrap();
    output.write_byte_record(&headers).unwrap();
    for (index, row) in input.records().enumerate() {
        let row = row.unwrap();
        let fields = row
            .iter()
            .enumerate()
            .map(|(column, value)| {
                if index == 0 && column == 1 {
                    (value.parse::<f64>().unwrap() + 1.0).to_string()
                } else {
                    value.to_owned()
                }
            })
            .collect::<Vec<_>>();
        output.write_record(fields).unwrap();
    }
    output.flush().unwrap();
    success(
        root,
        &[
            "corr",
            "--matrix",
            external.to_str().unwrap(),
            "--adjust",
            "none",
            "--geneset",
            &geneset,
            "--out",
            "edges:topk=5",
        ],
    );
    let custom_export = exported(
        root,
        &["export", "--what", "correlation", "--require-current"],
    );
    assert!(custom_export.exists());
    let contents = fs::read(&external).unwrap();
    fs::write(&external, b"gene_id\tRUN1\nG0001\t1\n").unwrap();
    let changed = command(
        root,
        &["export", "--what", "correlation", "--require-current"],
    );
    assert!(!changed.status.success());
    assert!(
        String::from_utf8_lossy(&changed.stderr)
            .contains("reconstructed current scientific request")
    );
    fs::write(&external, contents).unwrap();
    assert!(
        exported(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .exists()
    );
    let mut matrix_config = ProjectConfig::load(&paths.config_path).unwrap();
    let selection = matrix_config.corr.matrix.take();
    matrix_config.save(&paths.config_path).unwrap();
    db.set_project_config(&matrix_config).unwrap();
    let switched = command(
        root,
        &["export", "--what", "correlation", "--require-current"],
    );
    assert!(!switched.status.success());
    matrix_config.corr.matrix = selection;
    matrix_config.save(&paths.config_path).unwrap();
    db.set_project_config(&matrix_config).unwrap();
    fs::remove_file(&external).unwrap();
    assert!(
        !command(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .status
        .success()
    );
    assert!(custom_export.join("export-manifest.json").exists());
}
