//! Metadata types, references and design admission through the native workflow.
mod support;
use rnaa_core::{
    ProjectConfig, ProjectStageState,
    generations::OutputGroup,
    metadata::{ColumnSpec, ColumnType},
};
use std::{fs, path::Path};
use support::native::{command, setup, success};

fn json(path: &Path) -> serde_json::Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}

#[test]
fn native_de_and_correlation_use_the_same_typed_model_and_declared_reference() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let ages = [20, 22, 21, 23, 21, 23, 20, 22];
    let mut sheet = "run_accession\tcondition\tage\n".to_string();
    for (i, age) in ages.into_iter().enumerate() {
        sheet.push_str(&format!(
            "RUN{}\t{}\t{age}\n",
            i + 1,
            if i < 4 { "1" } else { "2" }
        ));
    }
    fs::write(paths.samplesheet_path(), sheet).unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.metadata.columns.insert(
        "condition".into(),
        ColumnSpec {
            kind: ColumnType::Factor,
            levels: vec!["2".into(), "1".into()],
            reference: Some("2".into()),
        },
    );
    config
        .metadata
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    config.save(&paths.config_path).unwrap();
    success(
        root,
        &[
            "deseq2",
            "--design",
            "~ age + condition",
            "--contrast",
            "condition",
            "1",
            "2",
            "--contrast",
            "condition",
            "2",
            "1",
        ],
    );
    let de = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let de_dir = de.paths(&paths).unwrap().de_project_dir(db.project_id());
    let contract = json(&de_dir.join("model_design.json"));
    assert_eq!(
        contract["metadata"]["columns"]["condition"]["type"],
        "factor"
    );
    assert_eq!(
        contract["metadata"]["columns"]["condition"]["levels"],
        serde_json::json!(["2", "1"])
    );
    assert_eq!(contract["metadata"]["columns"]["age"]["type"], "numeric");
    assert_eq!(contract["factor_references"]["condition"], "2");
    assert_eq!(contract["rank"], 3);
    assert_eq!(contract["residual_degrees_of_freedom"], 5);
    let mut forward = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(de_dir.join("de_condition_1_vs_2.tsv"))
        .unwrap();
    let mut reverse = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(de_dir.join("de_condition_2_vs_1.tsv"))
        .unwrap();
    let column = forward
        .headers()
        .unwrap()
        .iter()
        .position(|name| name == "log2FoldChange")
        .unwrap();
    let mut checked = 0;
    for (a, b) in forward.records().zip(reverse.records()) {
        let (a, b) = (a.unwrap(), b.unwrap());
        assert_eq!(&a[0], &b[0]);
        if let (Ok(a), Ok(b)) = (a[column].parse::<f64>(), b[column].parse::<f64>()) {
            assert!((a + b).abs() < 1e-10);
            checked += 1;
        }
    }
    assert!(
        checked > 400,
        "opposite resolved contrasts must yield opposite fitted effects"
    );
    success(root, &["deseq2"]);
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        de.generation_id
    );
    success(
        root,
        &[
            "corr",
            "--adjust",
            "residualize",
            "--model",
            "~ age + condition",
            "--out",
            "edges:topk=5",
        ],
    );
    let corr_contract = json(
        &paths
            .with_generation_heads(&db)
            .unwrap()
            .corr_dir
            .join("model_design.json"),
    );
    assert_eq!(contract["matrix"], corr_contract["matrix"]);
    assert_eq!(
        contract["coefficient_names"],
        corr_contract["coefficient_names"]
    );
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config
        .metadata
        .columns
        .insert("condition".into(), ColumnSpec::new(ColumnType::Numeric));
    config.save(&paths.config_path).unwrap();
    let stale = command(root, &["corr"]);
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("metadata schema are stale"));
    let states = db
        .list_runs()
        .unwrap()
        .into_iter()
        .map(|run| run.state)
        .collect::<Vec<_>>();
    let failed = command(root, &["deseq2"]);
    assert!(!failed.status.success());
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        de.generation_id
    );
    assert_eq!(
        states,
        db.list_runs()
            .unwrap()
            .into_iter()
            .map(|run| run.state)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        db.get_project_stage_state("deseq2").unwrap().unwrap().state,
        ProjectStageState::Failed
    );
}

#[test]
fn preflight_rejects_confounding_missing_values_and_insufficient_biological_replication_without_a_generation()
 {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    let head = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    let before = fs::read_dir(root.join(".rnaa_generations"))
        .unwrap()
        .count();
    let original = fs::read_to_string(paths.samplesheet_path()).unwrap();
    let mut confounded = "run_accession\tcondition\tbatch\n".to_string();
    for i in 1..=8 {
        let level = if i <= 4 { "A" } else { "B" };
        confounded.push_str(&format!("RUN{i}\t{level}\t{level}\n"));
    }
    fs::write(paths.samplesheet_path(), confounded).unwrap();
    let result = command(root, &["normalize", "--design", "~ batch + condition"]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("rank deficient"));
    fs::write(
        paths.samplesheet_path(),
        original.replace("RUN1\tA", "RUN1\tNA"),
    )
    .unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.metadata.missing_values = vec!["NA".into()];
    config.save(&paths.config_path).unwrap();
    let result = command(root, &["normalize", "--design", "~ condition"]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("missing values"));
    fs::write(paths.samplesheet_path(), original).unwrap();
    config = ProjectConfig::load(&paths.config_path).unwrap();
    config.metadata.minimum_replicates = 5;
    config.save(&paths.config_path).unwrap();
    let result = command(
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
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("minimum_replicates"));
    assert_eq!(
        fs::read_dir(root.join(".rnaa_generations"))
            .unwrap()
            .count(),
        before,
        "preflight must fail before staging a replacement generation"
    );
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        head.generation_id
    );
    assert_eq!(
        db.get_project_stage_state("deseq2").unwrap().unwrap().state,
        ProjectStageState::Failed
    );
}

#[test]
fn reference_transforms_and_full_indicators_keep_the_same_fitted_group_effect() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let mut baseline: Option<Vec<Option<f64>>> = None;
    for formula in [
        "~ condition",
        "~ relevel(condition, ref='B')",
        "~ 0 + condition",
    ] {
        success(
            root,
            &[
                "deseq2",
                "--design",
                formula,
                "--contrast",
                "condition",
                "B",
                "A",
            ],
        );
        let current = paths.with_generation_heads(&db).unwrap();
        let dir = current.de_project_dir(db.project_id());
        let contract = json(&dir.join("model_design.json"));
        if formula.contains("0 +") {
            assert!(contract["contrasts"][0]["reference"].is_null());
            assert_eq!(
                contract["contrasts"][0]["weights"],
                serde_json::json!([-1.0, 1.0])
            );
        }
        if formula.contains("relevel") {
            assert_eq!(contract["contrasts"][0]["reference"], "B");
            assert!(
                contract["contrasts"][0]["factor"]
                    .as_str()
                    .unwrap()
                    .contains("relevel")
            );
        }
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(dir.join("de_condition_B_vs_A.tsv"))
            .unwrap();
        let column = reader
            .headers()
            .unwrap()
            .iter()
            .position(|name| name == "log2FoldChange")
            .unwrap();
        let effects = reader
            .records()
            .map(|row| {
                row.unwrap()[column]
                    .parse::<f64>()
                    .ok()
                    .filter(|value| value.is_finite())
            })
            .collect::<Vec<_>>();
        if let Some(reference) = &baseline {
            let comparable = effects
                .iter()
                .zip(reference)
                .filter_map(|(a, b)| a.zip(*b))
                .collect::<Vec<(f64, f64)>>();
            assert!(comparable.len() > 400);
            let error = comparable
                .iter()
                .map(|(a, b)| (a - b).abs() / (1.0 + b.abs()))
                .fold(0.0_f64, f64::max);
            assert!(
                error < 1e-3,
                "equivalent coding changed the fitted contrast: {formula}, max relative error {error}"
            );
        } else {
            baseline = Some(effects);
        }
    }
}

#[test]
fn categorical_contrasts_preserve_the_two_group_zero_rule_with_a_third_expressed_group() {
    use rnaa_core::traits::ReferenceManager;
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let config = ProjectConfig::load(&paths.config_path).unwrap();
    let reference = rnaa_stages::refs::EnsemblReferenceManager::new()
        .unwrap()
        .ensure_reference(&paths, &db, &config)
        .unwrap();
    let mut sheet = "run_accession\tcondition\n".to_string();
    for i in 1..=8 {
        let run = format!("RUN{i}");
        let condition = if i <= 3 {
            "A"
        } else if i <= 6 {
            "B"
        } else {
            "C"
        };
        sheet.push_str(&format!("{run}\t{condition}\n"));
        let dir = paths.quant_run_dir(&run, "kallistors");
        let path = dir.join("abundance.tsv");
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)
            .unwrap();
        let header = reader.headers().unwrap().clone();
        let mut rows = reader
            .records()
            .map(|row| row.unwrap().iter().map(str::to_string).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        drop(reader);
        rows[0][3] = if condition == "C" { "1000" } else { "0" }.into();
        let mass = rows
            .iter()
            .map(|row| row[3].parse::<f64>().unwrap() / row[2].parse::<f64>().unwrap())
            .sum::<f64>();
        for row in &mut rows {
            row[4] = (1e6 * row[3].parse::<f64>().unwrap() / row[2].parse::<f64>().unwrap() / mass)
                .to_string();
        }
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)
            .unwrap();
        writer.write_record(&header).unwrap();
        for row in rows {
            writer.write_record(row).unwrap();
        }
        writer.flush().unwrap();
        drop(writer);
        rnaa_stages::quant::declare_external_quant(
            &run,
            &dir,
            &reference,
            "synthetic three-group character-contrast zero-rule fixture",
        )
        .unwrap();
    }
    fs::write(paths.samplesheet_path(), sheet).unwrap();
    success(
        root,
        &[
            "deseq2",
            "--design",
            "~ 0 + condition",
            "--contrast",
            "condition",
            "B",
            "A",
        ],
    );
    let dir = paths
        .with_generation_heads(&db)
        .unwrap()
        .de_project_dir(db.project_id());
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(dir.join("de_condition_B_vs_A.tsv"))
        .unwrap();
    let headers = reader.headers().unwrap().clone();
    let effect = headers
        .iter()
        .position(|name| name == "log2FoldChange")
        .unwrap();
    let pvalue = headers.iter().position(|name| name == "pvalue").unwrap();
    let row = reader
        .records()
        .map(|row| row.unwrap())
        .find(|row| &row[0] == "G0000")
        .unwrap();
    assert_eq!(row[effect].parse::<f64>().unwrap(), 0.0);
    assert_eq!(row[pvalue].parse::<f64>().unwrap(), 1.0);
}
