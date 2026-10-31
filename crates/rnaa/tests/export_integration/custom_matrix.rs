use super::*;

#[test]
fn strict_custom_matrix_export_does_not_require_a_de_result() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root);
    let input =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/native_workflow_counts.tsv");
    success(
        root,
        &[
            "corr",
            "--matrix",
            input.to_str().unwrap(),
            "--adjust",
            "none",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=5",
        ],
    );
    assert!(
        db.generation_head(rnaa_core::generations::OutputGroup::De, db.project_id())
            .unwrap()
            .is_none()
    );
    let output = exported(
        root,
        &["export", "--what", "correlation", "--require-current"],
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(output.join("export-manifest.json")).unwrap()).unwrap();
    assert_eq!(
        manifest["current_request_validation"]["requests"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        manifest["current_request_validation"]["external_matrix"]["bytes_match_request"],
        true
    );
    assert_eq!(
        manifest["current_request_validation"]["external_matrix"]["upstream_provenance_verified"],
        false
    );
}

#[test]
fn strict_custom_matrix_rebuilds_residualization_against_current_metadata() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, _db) = setup(root);
    let input =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/native_workflow_counts.tsv");
    let subset = root.join("subset.tsv");
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input)
        .unwrap();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&subset)
        .unwrap();
    writer
        .write_record(reader.headers().unwrap().iter().take(7))
        .unwrap();
    for row in reader.records() {
        writer.write_record(row.unwrap().iter().take(7)).unwrap();
    }
    writer.flush().unwrap();
    success(
        root,
        &[
            "corr",
            "--matrix",
            subset.to_str().unwrap(),
            "--adjust",
            "residualize",
            "--model",
            "~ condition",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=5",
        ],
    );
    assert!(
        exported(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .exists()
    );
    let original = fs::read(&subset).unwrap();
    let mut changed_values = original.clone();
    let body = changed_values
        .iter()
        .position(|byte| *byte == b'\n')
        .unwrap()
        + 1;
    let first_value = body
        + changed_values[body..]
            .iter()
            .position(|byte| *byte == b'\t')
            .unwrap()
        + 1;
    changed_values[first_value] = if changed_values[first_value] == b'9' {
        b'8'
    } else {
        b'9'
    };
    fs::write(&subset, changed_values).unwrap();
    assert!(
        exported(
            root,
            &[
                "export",
                "--what",
                "correlation",
                "--require-current-design"
            ]
        )
        .exists()
    );
    assert!(
        !command(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .status
        .success()
    );
    fs::write(&subset, original).unwrap();
    let sheet = fs::read_to_string(paths.samplesheet_path()).unwrap();
    fs::write(
        paths.samplesheet_path(),
        sheet.replace("RUN1\tA", "RUN1\tB"),
    )
    .unwrap();
    let changed = command(
        root,
        &["export", "--what", "correlation", "--require-current"],
    );
    assert!(!changed.status.success());
    assert!(
        String::from_utf8_lossy(&changed.stderr).contains("current design or used sample metadata")
    );
}

#[test]
fn strict_custom_matrix_with_module_analysis_checks_both_external_inputs() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, _db) = setup(root);
    let source =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/native_workflow_counts.tsv");
    let matrix = root.join("matrix.tsv");
    fs::copy(source, &matrix).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&matrix)
        .unwrap();
    let genes = reader
        .records()
        .map(|row| row.unwrap()[0].to_owned())
        .collect::<Vec<_>>();
    let gmt = root.join("pathways.gmt");
    fs::write(
        &gmt,
        format!(
            "FIRST\tfixture\t{}\nSECOND\tfixture\t{}\n",
            genes[..50].join("\t"),
            genes[50..100].join("\t")
        ),
    )
    .unwrap();
    let selected = root.join("selected.txt");
    fs::write(&selected, format!("{}\n", genes[..100].join("\n"))).unwrap();
    let geneset = format!("genes:{}", selected.display());
    success(
        root,
        &[
            "corr",
            "--matrix",
            matrix.to_str().unwrap(),
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
    assert!(
        exported(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .exists()
    );
    fs::write(&gmt, b"CHANGED\tfixture\tG0001\tG0002\n").unwrap();
    assert!(
        !command(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .status
        .success()
    );
}

#[test]
fn strict_custom_matrix_accounts_for_its_size_before_loading() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let source =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/native_workflow_counts.tsv");
    let matrix = root.join("matrix.tsv");
    fs::copy(source, &matrix).unwrap();
    success(
        root,
        &[
            "corr",
            "--matrix",
            matrix.to_str().unwrap(),
            "--adjust",
            "none",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=5",
        ],
    );
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.pool_root = root.join("limited-pool").display().to_string();
    config.resources.memory_mb = 128;
    config.resources.disk_headroom_mb = 0;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    fs::OpenOptions::new()
        .write(true)
        .open(&matrix)
        .unwrap()
        .set_len(8 * 1024 * 1024)
        .unwrap();
    let rejected = command(
        root,
        &["export", "--what", "correlation", "--require-current"],
    );
    assert!(!rejected.status.success());
    assert!(
        String::from_utf8_lossy(&rejected.stderr)
            .contains("resource request exceeds its own CPU or memory capacity")
    );
    assert_eq!(fs::read_dir(&paths.exports_dir).unwrap().count(), 0);
}

#[test]
fn growing_external_matrix_while_export_waits_is_rejected_before_publication() {
    use rnaa_core::resources::{Admission, ReservationOwner, ResourcePool, ResourceRequest};
    use std::{
        io::Write,
        process::{Command, Stdio},
        time::{Duration, Instant},
    };
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    let source =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/native_workflow_counts.tsv");
    let matrix = root.join("matrix.tsv");
    fs::copy(source, &matrix).unwrap();
    success(
        root,
        &[
            "corr",
            "--matrix",
            matrix.to_str().unwrap(),
            "--adjust",
            "none",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=5",
        ],
    );
    let original = fs::read(&matrix).unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    let pool_root = root.join("private-pool");
    config.resources.pool_root = pool_root.display().to_string();
    config.resources.cpu_threads = 1;
    config.resources.memory_mb = 512;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 30;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let pool = ResourcePool::open(&pool_root).unwrap();
    let held = match pool
        .try_reserve(
            config.resources.limits().unwrap(),
            ReservationOwner {
                project_root: root.into(),
                attempt_id: "matrix-budget-blocker".into(),
                stage: "test.blocker".into(),
                run_accession: None,
            },
            ResourceRequest {
                memory_bytes: 447 * 1024 * 1024,
                ..Default::default()
            },
        )
        .unwrap()
    {
        Admission::Granted(held) => held,
        Admission::Wait(_) => panic!("unexpected contention"),
    };
    let mut worker = Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args(["export", "--what", "correlation", "--require-current"])
            .env("RUST_LOG", "info")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .unwrap(),
    );
    let started = Instant::now();
    loop {
        let log = fs::read_to_string(paths.logs_dir.join("export.log")).unwrap_or_default();
        if log.lines().any(|line| {
            line.contains("export.write") && line.contains("waiting for shared resource capacity")
        }) {
            break;
        }
        assert!(
            worker.0.try_wait().unwrap().is_none(),
            "export exited before write admission"
        );
        assert!(
            started.elapsed() < Duration::from_secs(20),
            "export did not wait at write admission: {log}"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
    fs::OpenOptions::new()
        .append(true)
        .open(&matrix)
        .unwrap()
        .write_all(b"\nextra\t1\n")
        .unwrap();
    drop(held);
    let started = Instant::now();
    let status = loop {
        if let Some(status) = worker.0.try_wait().unwrap() {
            break status;
        }
        assert!(
            started.elapsed() < Duration::from_secs(20),
            "export did not reject the grown matrix"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(!status.success());
    assert_eq!(fs::read_dir(&paths.exports_dir).unwrap().count(), 0);
    fs::write(&matrix, original).unwrap();
    assert!(
        exported(
            root,
            &["export", "--what", "correlation", "--require-current"]
        )
        .exists()
    );
}
