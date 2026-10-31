//! Module report limits, resource-only reuse and failed admission through the CLI.
mod support;
use rnaa_core::{
    Database, ProjectConfig, ProjectPaths, generations::OutputGroup, util::compute_sha256,
};
use std::{fs, path::Path};
use support::native::{command, setup, success};

fn prepare(root: &Path) -> (ProjectPaths, Database, Vec<String>) {
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
    let selection = root.join("genes.txt");
    fs::write(
        &selection,
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
    let args = vec![
        "corr".into(),
        "--method".into(),
        "pearson".into(),
        "--adjust".into(),
        "none".into(),
        "--geneset".into(),
        format!("genes:{}", selection.display()),
        "--out".into(),
        "edges:topk=5".into(),
        "--module-gmt".into(),
        gmt.display().to_string(),
        "--module-min-size".into(),
        "2".into(),
        "--fgsea-permutations".into(),
        "100".into(),
    ];
    (paths, db, args)
}
fn with_top(args: &[String], top: &str) -> Vec<String> {
    let mut args = args.to_vec();
    args.extend(["--module-top".into(), top.into()]);
    args
}

#[test]
fn report_limits_change_the_report_and_preserve_full_scores_and_resource_only_reuse() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db, base) = prepare(root.path());
    let args = with_top(&base, "1");
    success(
        root.path(),
        &args.iter().map(String::as_str).collect::<Vec<_>>(),
    );
    let first = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let first_paths = first.paths(&paths).unwrap();
    let network_manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(first_paths.manifest_path("corr", db.project_id())).unwrap(),
    )
    .unwrap();
    assert_eq!(
        network_manifest["parameters"]["method"], "pearson",
        "module options must not override the selected ordinary network method"
    );
    let directory = first_paths.corr_project_dir(db.project_id());
    let scores = compute_sha256(&directory.join("module_eigengenes.tsv")).unwrap();
    let members = compute_sha256(&directory.join("module_membership.tsv")).unwrap();
    assert!(
        fs::read_to_string(directory.join("module_membership.tsv"))
            .unwrap()
            .lines()
            .count()
            > fs::read_to_string(directory.join("modules_top_genes.tsv"))
                .unwrap()
                .lines()
                .count()
    );
    let args = with_top(&base, "3");
    success(
        root.path(),
        &args.iter().map(String::as_str).collect::<Vec<_>>(),
    );
    let second = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    assert_ne!(first.generation_id, second.generation_id);
    let selected = second
        .paths(&paths)
        .unwrap()
        .corr_project_dir(db.project_id());
    assert_eq!(
        compute_sha256(&selected.join("module_eigengenes.tsv")).unwrap(),
        scores
    );
    assert_eq!(
        compute_sha256(&selected.join("module_membership.tsv")).unwrap(),
        members
    );
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.corr.block_rows = 7;
    config.corr.threads = 2;
    config.save(&paths.config_path).unwrap();
    success(
        root.path(),
        &args.iter().map(String::as_str).collect::<Vec<_>>(),
    );
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        second.generation_id
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(selected.join("module_fgsea_manifest.json")).unwrap())
            .unwrap();
    assert_eq!(manifest["parameters"]["limits"]["block_rows"], 64);
    assert_eq!(manifest["parameters"]["algorithm_version"], 2);
    assert_eq!(
        manifest["parameters"]["method"], "spearman",
        "module construction has its own declared method"
    );
}

#[test]
fn bad_pathways_and_module_resource_limits_preserve_the_published_generation() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db, base) = prepare(root.path());
    let args = with_top(&base, "1");
    success(
        root.path(),
        &args.iter().map(String::as_str).collect::<Vec<_>>(),
    );
    let head = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let directory = head
        .paths(&paths)
        .unwrap()
        .corr_project_dir(db.project_id());
    let hash = compute_sha256(&directory.join("module_eigengenes.tsv")).unwrap();
    let gmt = root.path().join("pathways.gmt");
    let original = fs::read(&gmt).unwrap();
    fs::write(&gmt, "DUP\tx\tG0000\tG0001\nDUP\tx\tG0002\tG0003\n").unwrap();
    let failure = command(
        root.path(),
        &args.iter().map(String::as_str).collect::<Vec<_>>(),
    );
    assert!(!failure.status.success());
    assert!(String::from_utf8_lossy(&failure.stderr).contains("duplicate GMT"));
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        head.generation_id
    );
    assert_eq!(
        compute_sha256(&directory.join("module_eigengenes.tsv")).unwrap(),
        hash
    );
    fs::write(gmt, original).unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.corr.max_output_mb = 1;
    config.save(&paths.config_path).unwrap();
    let changed = with_top(&base, "2");
    let failure = command(
        root.path(),
        &changed.iter().map(String::as_str).collect::<Vec<_>>(),
    );
    assert!(!failure.status.success());
    assert!(
        String::from_utf8_lossy(&failure.stderr).contains("module output estimate"),
        "{}",
        String::from_utf8_lossy(&failure.stderr)
    );
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        head.generation_id
    );
    assert_eq!(
        compute_sha256(&directory.join("module_eigengenes.tsv")).unwrap(),
        hash
    );
}
