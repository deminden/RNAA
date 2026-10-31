use std::{fs, process::Command};

fn run(root: &std::path::Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rnaa"))
        .arg("--root")
        .arg(root)
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn doctor_is_read_only_and_reports_project_errors_as_json() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path().join("not-yet-created");
    let output = run(&root, &["doctor", "--json"]);
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["checks"][0]["status"], "uninitialized");
    assert!(!root.exists(), "doctor must not initialize a project");

    assert!(run(&root, &["init"]).status.success());
    let before = fs::read(root.join("rnaa.toml")).unwrap();
    let output = run(&root, &["doctor", "--json"]);
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["ok"], true);
    assert!(report["project_id"].as_str().is_some());
    assert_eq!(report["engines"]["rsdeseq2"], "0.3.1");
    assert_eq!(fs::read(root.join("rnaa.toml")).unwrap(), before);

    fs::remove_file(root.join("rnaa.toml")).unwrap();
    let output = run(&root, &["doctor", "--json"]);
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        report["checks"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c["name"] == "configuration"
                && c["status"] == "ok"
                && c["detail"].as_str().unwrap().contains("catalog"))
    );
    let db = rusqlite::Connection::open(root.join("state.sqlite")).unwrap();
    db.execute("UPDATE projects SET config_json = 'not-json'", [])
        .unwrap();
    drop(db);
    let output = run(&root, &["doctor", "--json"]);
    assert!(!output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        report["checks"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c["name"] == "configuration" && c["status"] == "error")
    );
    let db = rusqlite::Connection::open(root.join("state.sqlite")).unwrap();
    db.execute(
        "UPDATE projects SET config_json = ?1",
        [r#"{"project":{"name":"bad","schema_version":2}}"#],
    )
    .unwrap();
    drop(db);
    let output = run(&root, &["doctor", "--json"]);
    assert!(!output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        report["checks"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c["name"] == "configuration"
                && c["status"] == "error"
                && c["detail"].as_str().unwrap().contains("schema_version"))
    );
    // A valid file takes precedence over a broken catalog copy, as in normal startup.
    fs::write(root.join("rnaa.toml"), &before).unwrap();
    assert!(run(&root, &["doctor", "--json"]).status.success());

    fs::write(root.join("rnaa.toml"), "not = [valid").unwrap();
    let output = run(&root, &["doctor", "--json"]);
    assert!(!output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["ok"], false);
    assert!(
        report["checks"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c["name"] == "configuration" && c["status"] == "error")
    );

    fs::write(root.join("rnaa.toml"), before).unwrap();
    fs::write(root.join("state.sqlite"), b"not a database").unwrap();
    let output = run(&root, &["doctor", "--json"]);
    assert!(!output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        report["checks"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c["name"] == "catalog" && c["status"] == "error")
    );
}
