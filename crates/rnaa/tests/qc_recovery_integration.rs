//! Reconstruct missing QC evidence independently of valid count payloads.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/recovery_fixture.rs"]
mod fixture;
use commands::{command, success};
use fixture::{Source, head, setup, setup_with_checksum};
use rnaa_core::{Database, ProjectPaths, generations::OutputGroup};
use std::{fs, path::Path, sync::atomic::Ordering};

fn report(paths: &ProjectPaths, db: &Database) -> std::path::PathBuf {
    head(db, OutputGroup::Preprocess)
        .paths(paths)
        .unwrap()
        .preprocess_run_dir("RUN1", "kallistors")
        .join("fasterp.json")
}
fn count_generations(db: &Database) -> usize {
    db.committed_generation_history(OutputGroup::Quant, "RUN1")
        .unwrap()
        .len()
}
fn assert_admitted(db: &Database) {
    assert!(
        db.eligibility("RUN1")
            .unwrap()
            .unwrap()
            .descriptor
            .eligible()
    );
}
fn recover(root: &Path) {
    let output = command(root, &["quant"]);
    assert!(
        output.status.success(),
        "{}\n{}\n{:?}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
        Database::open(root).unwrap().list_job_attempts(50).unwrap(),
    );
}

#[test]
fn missing_or_corrupt_only_report_is_rebuilt_after_all_reads_are_pruned_without_refitting_counts() {
    for corrupt in [false, true] {
        let source = Source::new();
        let temp = tempfile::tempdir().unwrap();
        let root = temp.path();
        let (paths, db) = setup(root, Some(&source), "none");
        let counts = head(&db, OutputGroup::Quant);
        let original = head(&db, OutputGroup::Preprocess);
        let damaged = report(&paths, &db);
        let seal = counts.payload_path.parent().unwrap().join("seal.json");
        let sealed = fs::read(&seal).unwrap();
        if corrupt {
            fs::write(&damaged, "damaged report").unwrap();
        } else {
            fs::remove_file(&damaged).unwrap();
        }
        recover(root);
        assert_eq!(head(&db, OutputGroup::Quant), counts);
        assert_eq!(count_generations(&db), 1);
        assert_eq!(fs::read(&seal).unwrap(), sealed);
        assert_ne!(head(&db, OutputGroup::Preprocess), original);
        assert!(
            db.verify_generation(OutputGroup::Preprocess, "RUN1")
                .unwrap()
        );
        assert!(db.verify_generation_head(&counts).unwrap());
        assert_admitted(&db);
        assert_eq!(source.requests.load(Ordering::SeqCst), 2);
        let plan = db
            .list_events(100)
            .unwrap()
            .into_iter()
            .find(|event| event.message == "planned QC report recovery")
            .unwrap();
        assert_eq!(plan.context["purpose"], "qc_report");
        assert_eq!(plan.context["action"], "download");
        assert_eq!(damaged.exists(), corrupt);
        if corrupt {
            assert_eq!(fs::read_to_string(&damaged).unwrap(), "damaged report");
        }
        success(root, &["quant"]);
        assert_eq!(head(&db, OutputGroup::Quant), counts);
        assert_eq!(source.requests.load(Ordering::SeqCst), 2);
    }
}

#[test]
fn retained_originals_rebuild_report_locally_and_preserve_original_counts() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, None, "both");
    let counts = head(&db, OutputGroup::Quant);
    let original = head(&db, OutputGroup::Preprocess);
    fs::write(report(&paths, &db), "damaged report").unwrap();
    recover(root);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    assert_eq!(count_generations(&db), 1);
    assert_ne!(head(&db, OutputGroup::Preprocess), original);
    assert_admitted(&db);
}

#[test]
fn retained_trimmed_reads_do_not_substitute_for_original_qc_inputs() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, None, "trimmed");
    let counts = head(&db, OutputGroup::Quant);
    let processed = head(&db, OutputGroup::Preprocess);
    fs::remove_file(report(&paths, &db)).unwrap();
    assert!(!command(root, &["quant"]).status.success());
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    assert_eq!(head(&db, OutputGroup::Preprocess), processed);
    assert_eq!(count_generations(&db), 1);
    assert!(
        !db.eligibility("RUN1")
            .unwrap()
            .unwrap()
            .descriptor
            .eligible()
    );
    assert!(
        db.get_run("RUN1")
            .unwrap()
            .unwrap()
            .last_error
            .unwrap()
            .contains("restore/register original reads")
    );
}

#[test]
fn failed_report_input_recovery_preserves_count_and_report_heads_then_retry_repairs_evidence() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "none");
    let counts = head(&db, OutputGroup::Quant);
    let processed = head(&db, OutputGroup::Preprocess);
    let downloaded = head(&db, OutputGroup::Download);
    fs::write(report(&paths, &db), "damaged report").unwrap();
    source.failed.store(true, Ordering::SeqCst);
    assert!(!command(root, &["quant"]).status.success());
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    assert_eq!(head(&db, OutputGroup::Preprocess), processed);
    assert_eq!(head(&db, OutputGroup::Download), downloaded);
    assert!(db.verify_generation_head(&counts).unwrap());
    source.failed.store(false, Ordering::SeqCst);
    recover(root);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    assert_eq!(count_generations(&db), 1);
    assert_admitted(&db);
}

#[test]
fn changed_unchecksummed_origin_during_report_recovery_invalidates_old_counts() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup_with_checksum(root, Some(&source), "none", false);
    let counts = head(&db, OutputGroup::Quant);
    fs::remove_file(report(&paths, &db)).unwrap();
    source.alternate.store(true, Ordering::SeqCst);
    recover(root);
    assert_ne!(head(&db, OutputGroup::Quant), counts);
    assert_eq!(count_generations(&db), 2);
    assert!(db.verify_generation_head(&counts).unwrap());
    assert_admitted(&db);
    assert_eq!(source.requests.load(Ordering::SeqCst), 2);
}

#[test]
fn reconstructed_report_applies_current_gate_and_never_bypasses_rejected_reads() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup_with_checksum(root, Some(&source), "none", false);
    let counts = head(&db, OutputGroup::Quant);
    let processed = head(&db, OutputGroup::Preprocess);
    fs::remove_file(report(&paths, &db)).unwrap();
    source.low_quality.store(true, Ordering::SeqCst);
    assert!(
        !command(root, &["quant", "--preprocess-bypass"])
            .status
            .success()
    );
    let rejected = head(&db, OutputGroup::Preprocess);
    assert_ne!(rejected, processed);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    assert_eq!(count_generations(&db), 1);
    assert_eq!(
        db.eligibility("RUN1").unwrap().unwrap().descriptor.state,
        rnaa_core::eligibility::EligibilityState::Rejected,
    );
    assert_eq!(
        db.get_preprocess_qc("RUN1").unwrap().unwrap().passed_reads,
        0
    );
    let requests = source.requests.load(Ordering::SeqCst);
    assert!(!command(root, &["quant"]).status.success());
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert_eq!(head(&db, OutputGroup::Preprocess), rejected);
    assert_eq!(head(&db, OutputGroup::Quant), counts);
}

#[test]
fn competing_report_recoveries_share_one_download_and_one_reconstructed_report() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "none");
    let counts = head(&db, OutputGroup::Quant);
    fs::remove_file(report(&paths, &db)).unwrap();
    source.delay_ms.store(500, Ordering::SeqCst);
    let mut first = std::process::Command::new(env!("CARGO_BIN_EXE_rnaa"))
        .args(["--root", root.to_str().unwrap(), "quant"])
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    let start = std::time::Instant::now();
    while source.requests.load(Ordering::SeqCst) < 2 {
        assert!(start.elapsed() < std::time::Duration::from_secs(15));
        assert!(first.try_wait().unwrap().is_none());
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    recover(root);
    let result = first.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(head(&db, OutputGroup::Quant), counts);
    assert_eq!(count_generations(&db), 1);
    assert_eq!(
        source.requests.load(Ordering::SeqCst),
        2,
        "HTTP: {:?}; Recovery events: {:#?}; attempts: {:#?}",
        source.request_lines.lock().unwrap(),
        db.list_events(150)
            .unwrap()
            .into_iter()
            .filter(|event| matches!(
                event.stage.as_str(),
                "download" | "preprocess" | "quant" | "recovery"
            ))
            .collect::<Vec<_>>(),
        db.list_job_attempts(30).unwrap()
    );
    assert_eq!(
        db.committed_generation_history(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .len(),
        2,
    );
    assert_admitted(&db);
}
