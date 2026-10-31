//! Scientific reuse ignores unrelated metadata and consumes current used values.
mod support;
use rnaa_core::{
    Database, ManifestArtifact, ProjectConfig, ProjectPaths,
    generations::{GenerationHead, OutputGroup},
    metadata::{ColumnSpec, ColumnType},
};
use rnaa_formats::matrix::read_matrix_tsv;
use std::fs;
use support::native::{command, setup, success};

fn head(db: &Database, group: OutputGroup) -> GenerationHead {
    db.generation_head(group, db.project_id()).unwrap().unwrap()
}
fn sheet(paths: &ProjectPaths, ages: &[&str; 8], reverse: bool) {
    let mut text = "notes\tage\trun_accession\tcondition\tread1_path\tunused_score\n".to_string();
    let order = if reverse {
        (0..8).rev().collect::<Vec<_>>()
    } else {
        (0..8).collect()
    };
    for i in order {
        text.push_str(&format!(
            "annotation {i}\t{}\tRUN{}\t{}\t/new/location/{i}.fastq\tunavailable\n",
            ages[i],
            i + 1,
            if i < 4 { "A" } else { "B" }
        ));
    }
    fs::write(paths.samplesheet_path(), text).unwrap();
}
fn configure_age(paths: &ProjectPaths) {
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config
        .metadata
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    config.save(&paths.config_path).unwrap();
}

#[test]
fn operational_and_unused_metadata_reuse_results_and_retain_each_attempts_source_witness() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let de = head(&db, OutputGroup::De);
    let corr = head(&db, OutputGroup::Corr);
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.metadata
        .insert("operational_note".into(), "new catalog annotation".into());
    db.upsert_run(&run).unwrap();
    sheet(
        &paths,
        &["21", "22", "20", "23", "21", "23", "20", "22"],
        true,
    );
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config
        .metadata
        .columns
        .insert("unused_score".into(), ColumnSpec::new(ColumnType::Numeric));
    config
        .metadata
        .columns
        .insert("notes".into(), ColumnSpec::new(ColumnType::Ignore));
    config.metadata.missing_values = vec!["not collected".into()];
    config.metadata.minimum_replicates = 3;
    config.corr.model = "unused when adjustment is none".into();
    config.save(&paths.config_path).unwrap();
    success(root, &["corr"]);
    success(root, &["normalize"]);
    assert_eq!(head(&db, OutputGroup::De), de);
    assert_eq!(head(&db, OutputGroup::Corr), corr);
    assert!(db.verify_generation_head(&de).unwrap());
    assert!(db.verify_generation_head(&corr).unwrap());
    let source =
        ManifestArtifact::snapshot("SAMPLESHEET_SOURCE", &paths.samplesheet_path()).unwrap();
    let events = db.list_events(1000).unwrap();
    let reused = events
        .iter()
        .filter(|event| event.message.contains("reconciled from existing outputs"))
        .collect::<Vec<_>>();
    assert_eq!(reused.len(), 2);
    for event in reused {
        let witnesses: Vec<ManifestArtifact> =
            serde_json::from_value(event.context["metadata_witnesses"].clone()).unwrap();
        let witness = witnesses
            .iter()
            .find(|w| w.kind == "SAMPLESHEET_SOURCE")
            .unwrap();
        assert_eq!(witness.checksum, source.checksum);
        assert_ne!(witness.path, source.path);
        assert!(witness.matches_file().unwrap());
        assert_eq!(
            event.context["configuration_observation"]["schema"]["columns"]["unused_score"]["type"],
            "numeric"
        );
    }
    let descriptor: serde_json::Value = serde_json::from_slice(
        &fs::read(
            de.paths(&paths)
                .unwrap()
                .de_project_dir(db.project_id())
                .join("rnaa_deseq2_request.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert!(
        descriptor["request"]["inputs"]
            .as_array()
            .unwrap()
            .iter()
            .all(|input| !["SAMPLESHEET", "SAMPLESHEET_SOURCE", "METADATA"]
                .contains(&input["role"].as_str().unwrap()))
    );
    assert!(
        descriptor["request"]["parameters"]["design_contract"]["columns"]
            .get("notes")
            .is_none()
    );
    assert!(
        descriptor["input_witnesses"]
            .as_array()
            .unwrap()
            .iter()
            .any(|w| w["kind"] == "SAMPLESHEET_SOURCE")
    );
}

#[test]
fn correlation_only_covariate_changes_reuse_de_and_fit_current_values() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    configure_age(&paths);
    sheet(
        &paths,
        &["20", "22", "21", "23", "21", "23", "20", "22"],
        false,
    );
    success(root, &["normalize"]);
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
    let de = head(&db, OutputGroup::De);
    let previous = head(&db, OutputGroup::Corr);
    let old_matrix = read_matrix_tsv(
        &previous
            .paths(&paths)
            .unwrap()
            .corr_dir
            .join("adjusted_matrix.tsv"),
    )
    .unwrap();
    sheet(
        &paths,
        &["34", "20", "29", "22", "25", "31", "23", "28"],
        true,
    );
    // No DE rerun: its model never used age, but the consumer does.
    success(root, &["corr"]);
    let current = head(&db, OutputGroup::Corr);
    assert_ne!(current.generation_id, previous.generation_id);
    assert_eq!(head(&db, OutputGroup::De), de);
    let model: serde_json::Value = serde_json::from_slice(
        &fs::read(
            current
                .paths(&paths)
                .unwrap()
                .corr_dir
                .join("model_design.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(
        model["metadata"]["columns"]["age"]["values"],
        serde_json::json!([34.0, 20.0, 29.0, 22.0, 25.0, 31.0, 23.0, 28.0])
    );
    let adjusted = read_matrix_tsv(
        &current
            .paths(&paths)
            .unwrap()
            .corr_dir
            .join("adjusted_matrix.tsv"),
    )
    .unwrap();
    let difference = old_matrix
        .values
        .iter()
        .flatten()
        .zip(adjusted.values.iter().flatten())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(difference > 1e-4);
    let ages = [34.0, 20.0, 29.0, 22.0, 25.0, 31.0, 23.0, 28.0];
    for row in adjusted.values.iter().take(20) {
        let dot = row.iter().zip(ages).map(|(x, age)| x * age).sum::<f64>();
        let scale = 1.0
            + row
                .iter()
                .zip(ages)
                .map(|(x, age)| (x * age).abs())
                .sum::<f64>();
        assert!(
            dot.abs() < scale * 1e-8,
            "residuals must remove the current covariate"
        );
    }
    sheet(
        &paths,
        &[
            "34.0", "20.0", "29.0", "22.0", "25.0", "31.0", "23.0", "28.0",
        ],
        false,
    );
    success(root, &["corr"]);
    success(root, &["normalize"]);
    assert_eq!(head(&db, OutputGroup::Corr), current);
    assert_eq!(head(&db, OutputGroup::De), de);
    // Restoring the original covariate values selects its exact historical fit.
    sheet(
        &paths,
        &["20", "22", "21", "23", "21", "23", "20", "22"],
        true,
    );
    success(root, &["corr"]);
    assert_eq!(head(&db, OutputGroup::Corr), previous);
    assert_eq!(head(&db, OutputGroup::De), de);
    assert!(db.verify_generation_head(&current).unwrap());
}

#[test]
fn changes_to_producer_model_values_reject_stale_matrices_then_recompute() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    success(root, &["normalize"]);
    success(root, &["corr", "--adjust", "none", "--out", "edges:topk=5"]);
    let de = head(&db, OutputGroup::De);
    let corr = head(&db, OutputGroup::Corr);
    let text = fs::read_to_string(paths.samplesheet_path()).unwrap();
    fs::write(
        paths.samplesheet_path(),
        text.replace("\tA\n", "\tcontrol\n")
            .replace("\tB\n", "\ttreated\n"),
    )
    .unwrap();
    let rejected = command(root, &["corr"]);
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("DE metadata is stale"));
    assert_eq!(head(&db, OutputGroup::De), de);
    assert_eq!(head(&db, OutputGroup::Corr), corr);
    success(root, &["normalize"]);
    success(root, &["corr"]);
    assert_ne!(head(&db, OutputGroup::De).generation_id, de.generation_id);
    let updated_de = head(&db, OutputGroup::De);
    let old_matrix = rnaa_core::util::compute_sha256(
        &de.paths(&paths)
            .unwrap()
            .de_project_dir(db.project_id())
            .join("vst.tsv"),
    )
    .unwrap();
    let new_matrix = rnaa_core::util::compute_sha256(
        &updated_de
            .paths(&paths)
            .unwrap()
            .de_project_dir(db.project_id())
            .join("vst.tsv"),
    )
    .unwrap();
    if old_matrix == new_matrix {
        assert_eq!(head(&db, OutputGroup::Corr), corr);
    } else {
        assert_ne!(
            head(&db, OutputGroup::Corr).generation_id,
            corr.generation_id
        );
    }
    assert!(db.verify_generation_head(&de).unwrap());
    assert!(db.verify_generation_head(&corr).unwrap());
}
