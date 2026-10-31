//! Pruned-input recovery through real HTTP, preprocessing and quantification.
#[path = "support/commands.rs"]
mod commands;
#[path = "support/recovery_fixture.rs"]
mod fixture;
use commands::{command, success};
use fixture::{Source, altered_reads, head, replace, setup};
use rnaa_core::{
    ArtifactKind, ProjectConfig, RunState, generations::OutputGroup, util::compute_md5,
};
use std::{
    fs,
    path::Path,
    sync::{Arc, atomic::Ordering},
    thread,
    time::{Duration, Instant},
};

#[test]
fn historical_counts_return_without_redownload_after_all_inputs_were_pruned() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "none");
    let original = head(&db, OutputGroup::Quant);
    let seal = original.payload_path.parent().unwrap().join("seal.json");
    let original_seal = fs::read(&seal).unwrap();
    replace(root, "none");
    let replacement = head(&db, OutputGroup::Quant);
    assert_ne!(original, replacement);
    let requests = source.requests.load(Ordering::SeqCst);
    source.failed.store(true, Ordering::SeqCst);
    success(root, &["quant", "--single-fragment-length", "200"]);
    assert_eq!(head(&db, OutputGroup::Quant), original);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert_eq!(fs::read(seal).unwrap(), original_seal);
    assert!(db.verify_generation_head(&replacement).unwrap());
    let binding = db
        .generation_input_binding(OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    assert_eq!(binding.generation_id, original.generation_id);
    assert!(
        binding
            .inputs
            .iter()
            .all(|input| !Path::new(&input.path).exists())
    );
    assert!(
        db.eligibility("RUN1")
            .unwrap()
            .unwrap()
            .descriptor
            .eligible()
    );
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), original);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    let history = db
        .committed_generation_history(OutputGroup::Quant, "RUN1")
        .unwrap();
    assert_eq!(history.len(), 2);
    let inspected = command(
        root,
        &[
            "lineage",
            "--group",
            "quant",
            "--all-requests",
            "--history",
            "--verify",
        ],
    );
    assert!(
        inspected.status.success(),
        "{}",
        String::from_utf8_lossy(&inspected.stderr)
    );
    assert!(
        ProjectConfig::load(&paths.config_path)
            .unwrap()
            .quant
            .preprocess
    );
}

#[test]
fn historical_selection_retains_newer_equivalent_inputs_with_durable_binding_receipts() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "none");
    let original = head(&db, OutputGroup::Quant);
    replace(root, "both");
    let replacement = head(&db, OutputGroup::Quant);
    let download = head(&db, OutputGroup::Download);
    let preprocessing = head(&db, OutputGroup::Preprocess);
    assert!(
        !db.list_artifacts_for_run(Some("RUN1"))
            .unwrap()
            .iter()
            .filter(|item| item.kind.as_str().starts_with("FASTQ"))
            .collect::<Vec<_>>()
            .is_empty()
    );
    let requests = source.requests.load(Ordering::SeqCst);
    success(
        root,
        &[
            "quant",
            "--single-fragment-length",
            "200",
            "--fastq-retention",
            "none",
        ],
    );
    assert_eq!(head(&db, OutputGroup::Quant), original);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert!(
        db.list_artifacts_for_run(Some("RUN1"))
            .unwrap()
            .iter()
            .all(|item| !item.kind.as_str().starts_with("FASTQ"))
    );
    assert!(db.verify_generation_head(&download).unwrap());
    assert!(db.verify_generation_head(&preprocessing).unwrap());
    assert!(
        db.eligibility("RUN1")
            .unwrap()
            .unwrap()
            .descriptor
            .eligible()
    );
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), original);
    replace(root, "none");
    assert_eq!(head(&db, OutputGroup::Quant), replacement);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert!(db.verify_generation_head(&download).unwrap());
    assert!(db.verify_generation_head(&preprocessing).unwrap());
    assert_eq!(
        db.committed_generation_history(OutputGroup::Quant, "RUN1")
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn restoring_preprocessing_reuses_old_counts_and_rebinds_its_pruned_report_sources() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "none");
    let original = head(&db, OutputGroup::Quant);
    let preprocessing = head(&db, OutputGroup::Preprocess);
    success(
        root,
        &["quant", "--no-preprocess", "--fastq-retention", "both"],
    );
    let raw_counts = head(&db, OutputGroup::Quant);
    assert_ne!(raw_counts, original);
    let requests = source.requests.load(Ordering::SeqCst);
    success(
        root,
        &["quant", "--preprocess", "--fastq-retention", "none"],
    );
    assert_eq!(head(&db, OutputGroup::Quant), original);
    assert_eq!(head(&db, OutputGroup::Preprocess), preprocessing);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert!(
        db.generation_input_binding(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .is_some()
    );
    assert!(
        db.eligibility("RUN1")
            .unwrap()
            .unwrap()
            .descriptor
            .eligible()
    );
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), original);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Quant, "RUN1")
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn corrupt_historical_counts_cannot_replace_current_counts_when_recovery_fails() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "none");
    let original = head(&db, OutputGroup::Quant);
    replace(root, "none");
    let replacement = head(&db, OutputGroup::Quant);
    let corrupt = original
        .paths(&paths)
        .unwrap()
        .quant_run_dir("RUN1", "kallistors")
        .join("abundance.tsv");
    fs::write(&corrupt, "corrupt historical table").unwrap();
    source.failed.store(true, Ordering::SeqCst);
    let failed = command(root, &["quant", "--single-fragment-length", "200"]);
    assert!(!failed.status.success());
    assert_eq!(head(&db, OutputGroup::Quant), replacement);
    assert!(db.verify_generation_head(&replacement).unwrap());
    source.failed.store(false, Ordering::SeqCst);
    success(root, &["quant"]);
    let rebuilt = head(&db, OutputGroup::Quant);
    assert_ne!(rebuilt, original);
    assert_ne!(rebuilt, replacement);
    assert_eq!(rebuilt.request_sha256, original.request_sha256);
    assert_eq!(
        fs::read_to_string(corrupt).unwrap(),
        "corrupt historical table"
    );
    assert!(db.verify_generation_head(&rebuilt).unwrap());
}

#[test]
fn an_older_healthy_report_repairs_current_qc_evidence_without_reads_or_network() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "none");
    let old_preprocessing = head(&db, OutputGroup::Preprocess);
    replace(root, "none");
    let current_counts = head(&db, OutputGroup::Quant);
    let bad_preprocessing = head(&db, OutputGroup::Preprocess);
    assert_ne!(old_preprocessing, bad_preprocessing);
    let corrupt = bad_preprocessing
        .paths(&paths)
        .unwrap()
        .preprocess_run_dir("RUN1", "kallistors")
        .join("fasterp.json");
    fs::write(&corrupt, "corrupt current report").unwrap();
    let requests = source.requests.load(Ordering::SeqCst);
    source.failed.store(true, Ordering::SeqCst);
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Preprocess), old_preprocessing);
    assert_eq!(head(&db, OutputGroup::Quant), current_counts);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert!(
        db.eligibility("RUN1")
            .unwrap()
            .unwrap()
            .descriptor
            .eligible()
    );
    assert_eq!(
        fs::read_to_string(corrupt).unwrap(),
        "corrupt current report"
    );
}

#[test]
fn changed_local_read_bytes_do_not_reuse_an_earlier_request_with_the_same_options() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "both");
    let original = head(&db, OutputGroup::Quant);
    replace(root, "both");
    let raw = db
        .list_artifacts_for_run(Some("RUN1"))
        .unwrap()
        .into_iter()
        .find(|item| item.kind == ArtifactKind::FastqSingle)
        .unwrap();
    let mut changed = fs::read(&raw.path).unwrap();
    let offset = changed.iter().position(|byte| *byte == b'\n').unwrap() + 1;
    changed[offset] = if changed[offset] == b'A' { b'C' } else { b'A' };
    fs::write(&raw.path, changed).unwrap();
    success(root, &["quant", "--single-fragment-length", "200"]);
    let current = head(&db, OutputGroup::Quant);
    assert_ne!(current, original);
    assert_ne!(current.request_sha256, original.request_sha256);
    assert!(db.verify_generation_head(&current).unwrap());
}

#[test]
fn historical_reuse_preserves_requests_computed_from_retained_trimmed_inputs_only() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "trimmed");
    let original = head(&db, OutputGroup::Quant);
    let requests = source.requests.load(Ordering::SeqCst);
    replace(root, "trimmed");
    let trimmed_only = head(&db, OutputGroup::Quant);
    let manifest: rnaa_core::StageManifest = serde_json::from_slice(
        &fs::read(
            trimmed_only
                .paths(&paths)
                .unwrap()
                .quant_run_dir("RUN1", "kallistors")
                .join("rnaa_quant_manifest.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert!(
        manifest
            .input_artifacts
            .iter()
            .all(|input| input.kind.starts_with("FASTQ_TRIMMED_"))
    );
    success(root, &["quant", "--single-fragment-length", "200"]);
    assert_eq!(head(&db, OutputGroup::Quant), original);
    replace(root, "none");
    assert_eq!(head(&db, OutputGroup::Quant), trimmed_only);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), trimmed_only);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Quant, "RUN1")
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn competing_history_requests_select_one_result_without_recovery_or_new_payloads() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "none");
    let original = head(&db, OutputGroup::Quant);
    replace(root, "none");
    let jobs = db
        .list_job_attempts(1000)
        .unwrap()
        .iter()
        .filter(|attempt| attempt.stage == "quant")
        .count();
    let requests = source.requests.load(Ordering::SeqCst);
    source.failed.store(true, Ordering::SeqCst);
    let barrier = Arc::new(std::sync::Barrier::new(2));
    let mut handles = Vec::new();
    for _ in 0..2 {
        let root = root.to_path_buf();
        let barrier = barrier.clone();
        handles.push(thread::spawn(move || {
            barrier.wait();
            command(&root, &["quant", "--single-fragment-length", "200"])
        }));
    }
    for handle in handles {
        let result = handle.join().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    assert_eq!(head(&db, OutputGroup::Quant), original);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert_eq!(
        db.list_job_attempts(1000)
            .unwrap()
            .iter()
            .filter(|attempt| attempt.stage == "quant")
            .count(),
        jobs + 1
    );
    assert_eq!(
        db.committed_generation_history(OutputGroup::Quant, "RUN1")
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn local_trimmed_only_history_keeps_its_original_source_proof_after_all_reads_are_pruned() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "trimmed");
    let original = head(&db, OutputGroup::Quant);
    replace(root, "none");
    let trimmed_only = head(&db, OutputGroup::Quant);
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), trimmed_only);
    success(root, &["quant", "--single-fragment-length", "200"]);
    assert_eq!(head(&db, OutputGroup::Quant), original);
    replace(root, "none");
    assert_eq!(head(&db, OutputGroup::Quant), trimmed_only);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Quant, "RUN1")
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn current_source_content_is_checked_against_pruned_count_history_not_only_source_selection() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "none");
    let original = head(&db, OutputGroup::Quant);
    source.alternate.store(true, Ordering::SeqCst);
    let checksum = root.join("updated-source.fastq");
    fs::write(&checksum, altered_reads()).unwrap();
    let md5 = compute_md5(&checksum).unwrap();
    fs::remove_file(checksum).unwrap();
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.remote_files[0].md5 = Some(md5);
    db.upsert_run(&run).unwrap();
    success(root, &["quant"]);
    let current = head(&db, OutputGroup::Quant);
    assert_ne!(current, original);
    assert_ne!(current.request_sha256, original.request_sha256);
    // Simulate a retained old head from an earlier version after a new source
    // has been selected and both generations' reads have been lawfully pruned.
    let run = db.get_run("RUN1").unwrap().unwrap();
    let mut lease =
        rnaa_core::jobs::JobLease::try_acquire(&db, "quant", Some("RUN1"), &[run.state])
            .unwrap()
            .unwrap();
    let owned = lease.database();
    owned.set_job_request(&original.request_sha256).unwrap();
    let verified = owned.prepare_generation_reuse(&original).unwrap().unwrap();
    owned
        .reattach_generation(verified, |db| {
            db.set_run_state("RUN1", RunState::QuantDone, None)
        })
        .unwrap();
    lease
        .finish(rnaa_core::jobs::JobOutcome::Succeeded, None)
        .unwrap();
    drop(lease);
    assert_eq!(head(&db, OutputGroup::Quant), original);
    let requests = source.requests.load(Ordering::SeqCst);
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), current);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Quant, "RUN1")
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn restored_download_source_reselects_the_original_preprocessing_and_counts_offline() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "both");
    let original_download = head(&db, OutputGroup::Download);
    let original_preprocessing = head(&db, OutputGroup::Preprocess);
    let original_counts = head(&db, OutputGroup::Quant);
    let original_sources = db.get_run("RUN1").unwrap().unwrap().remote_files;
    source.alternate.store(true, Ordering::SeqCst);
    let checksum = root.join("updated-source.fastq");
    fs::write(&checksum, altered_reads()).unwrap();
    let md5 = compute_md5(&checksum).unwrap();
    fs::remove_file(checksum).unwrap();
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.remote_files[0].md5 = Some(md5);
    db.upsert_run(&run).unwrap();
    success(root, &["quant"]);
    let replacement = head(&db, OutputGroup::Quant);
    assert_ne!(replacement, original_counts);
    let requests = source.requests.load(Ordering::SeqCst);
    source.failed.store(true, Ordering::SeqCst);
    let mut run = db.get_run("RUN1").unwrap().unwrap();
    run.remote_files = original_sources;
    db.upsert_run(&run).unwrap();
    success(root, &["download"]);
    assert_eq!(head(&db, OutputGroup::Download), original_download);
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Preprocess), original_preprocessing);
    assert_eq!(head(&db, OutputGroup::Quant), original_counts);
    assert_eq!(source.requests.load(Ordering::SeqCst), requests);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Download, "RUN1")
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        db.committed_generation_history(OutputGroup::Preprocess, "RUN1")
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        db.committed_generation_history(OutputGroup::Quant, "RUN1")
            .unwrap()
            .len(),
        2
    );
    assert!(db.verify_generation_head(&replacement).unwrap());
}

#[test]
fn pruned_reads_recover_from_a_durable_pending_state_and_publish_exact_dependency_edges() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "none");
    let old = head(&db, OutputGroup::Quant);
    let raw = head(&db, OutputGroup::Download);
    assert!(
        db.list_artifacts_for_run(Some("RUN1"))
            .unwrap()
            .iter()
            .all(|a| !a.kind.as_str().starts_with("FASTQ"))
    );
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.single_fragment_length = 180.0;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    db.set_run_state(
        "RUN1",
        RunState::Resolved,
        Some("simulated stop after recovery planning"),
    )
    .unwrap();
    success(root, &["quant"]);
    let new = head(&db, OutputGroup::Quant);
    let download = head(&db, OutputGroup::Download);
    let preprocess = head(&db, OutputGroup::Preprocess);
    assert_ne!(new.generation_id, old.generation_id);
    assert_ne!(download.generation_id, raw.generation_id);
    assert_eq!(source.requests.load(Ordering::SeqCst), 2);
    let lineage = db.generation_lineage(&new).unwrap();
    assert!(lineage.recorded);
    assert!(
        lineage
            .dependencies
            .iter()
            .any(|d| d.producer_generation_id == download.generation_id)
    );
    assert!(
        lineage
            .dependencies
            .iter()
            .any(|d| d.producer_generation_id == preprocess.generation_id)
    );
    assert!(
        lineage
            .dependencies
            .iter()
            .any(|d| d.input.kind == "REFERENCE_INDEX")
    );
    assert!(
        db.generation_lineage(&preprocess)
            .unwrap()
            .dependencies
            .iter()
            .any(|d| d.producer_generation_id == download.generation_id)
    );
    assert!(
        db.dependent_generations(&download)
            .unwrap()
            .iter()
            .any(|h| h.generation_id == new.generation_id)
    );
    assert!(db.verify_generation_head(&old).unwrap());
    assert!(db.verify_generation_head(&raw).unwrap());
    assert!(db.verify_generation_head(&new).unwrap());
    let report = command(
        root,
        &[
            "lineage",
            "--group",
            "quant",
            "--subject",
            "RUN1",
            "--history",
            "--verify",
        ],
    );
    assert!(
        report.status.success(),
        "{}",
        String::from_utf8_lossy(&report.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&report.stdout).unwrap();
    let entries = report["generations"].as_array().unwrap();
    assert!(entries.len() >= 4);
    assert!(
        entries
            .iter()
            .all(|entry| entry["payload_verified"] == true)
    );
    assert!(
        entries
            .iter()
            .any(|entry| entry["lineage"]["head"]["generation_id"] == new.generation_id)
    );
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), new);
    assert_eq!(source.requests.load(Ordering::SeqCst), 2);
}

#[test]
fn retained_trimmed_reads_avoid_redownload_and_reprocessing() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "trimmed");
    let old = head(&db, OutputGroup::Quant);
    let processed = head(&db, OutputGroup::Preprocess);
    source.failed.store(true, Ordering::SeqCst);
    replace(root, "trimmed");
    let new = head(&db, OutputGroup::Quant);
    assert_ne!(new.generation_id, old.generation_id);
    assert_eq!(head(&db, OutputGroup::Preprocess), processed);
    assert_eq!(source.requests.load(Ordering::SeqCst), 1);
    assert!(
        db.generation_lineage(&new)
            .unwrap()
            .dependencies
            .iter()
            .any(|d| d.producer_generation_id == processed.generation_id)
    );
    assert!(db.verify_generation_head(&old).unwrap());
    success(root, &["quant"]);
    assert_eq!(head(&db, OutputGroup::Quant), new);
    let report_path = processed
        .paths(&paths)
        .unwrap()
        .preprocess_run_dir("RUN1", "kallistors")
        .join("fasterp.json");
    fs::write(report_path, b"corrupt QC report").unwrap();
    let rejected = command(
        root,
        &[
            "lineage",
            "--group",
            "quant",
            "--subject",
            "RUN1",
            "--history",
            "--verify",
        ],
    );
    assert!(!rejected.status.success());
    let report: serde_json::Value = serde_json::from_slice(&rejected.stdout).unwrap();
    assert!(
        report["generations"]
            .as_array()
            .unwrap()
            .iter()
            .any(|entry| entry["payload_verified"] == false)
    );
    assert!(
        db.verify_generation_head(&new).unwrap(),
        "QC report admission is independent of count payload validity"
    );
}

#[test]
fn disabling_preprocessing_recovers_originals_instead_of_treating_trimmed_reads_as_raw() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root, Some(&source), "trimmed");
    let previous = head(&db, OutputGroup::Quant);
    let processed = head(&db, OutputGroup::Preprocess);
    success(
        root,
        &["quant", "--no-preprocess", "--fastq-retention", "both"],
    );
    let current = head(&db, OutputGroup::Quant);
    assert_ne!(current.generation_id, previous.generation_id);
    assert_eq!(
        source.requests.load(Ordering::SeqCst),
        2,
        "{}",
        fs::read_to_string(paths.logs_dir.join("quant.log")).unwrap_or_default()
    );
    assert_eq!(head(&db, OutputGroup::Preprocess), processed);
    let manifest: rnaa_core::StageManifest = serde_json::from_slice(
        &fs::read(
            current
                .paths(&paths)
                .unwrap()
                .quant_run_dir("RUN1", "kallistors")
                .join("rnaa_quant_manifest.json"),
        )
        .unwrap(),
    )
    .unwrap();
    assert!(
        manifest
            .input_artifacts
            .iter()
            .all(|file| file.kind == "FASTQ_SINGLE")
    );
    assert_eq!(
        manifest.parameters["request"]["parameters"]["preprocess"],
        false
    );
    assert!(db.verify_generation_head(&previous).unwrap());
}

#[test]
fn failed_restaging_preserves_prior_results_and_retry_recovers() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "none");
    let old = head(&db, OutputGroup::Quant);
    let download = head(&db, OutputGroup::Download);
    source.failed.store(true, Ordering::SeqCst);
    let failed = command(root, &["quant", "--single-fragment-length", "180"]);
    assert!(!failed.status.success());
    assert_eq!(head(&db, OutputGroup::Quant), old);
    assert_eq!(head(&db, OutputGroup::Download), download);
    assert!(db.verify_generation_head(&old).unwrap());
    assert!(
        db.list_job_attempts(100)
            .unwrap()
            .iter()
            .any(|a| a.stage == "download" && a.state == "FAILED")
    );
    source.failed.store(false, Ordering::SeqCst);
    success(root, &["quant"]);
    assert_ne!(
        head(&db, OutputGroup::Quant).generation_id,
        old.generation_id
    );
    assert!(db.verify_generation_head(&old).unwrap());
}

#[test]
fn missing_local_inputs_block_without_erasing_counts_or_running_a_false_replacement() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, None, "none");
    let old = head(&db, OutputGroup::Quant);
    let jobs = db
        .list_job_attempts(100)
        .unwrap()
        .iter()
        .filter(|a| a.stage == "quant")
        .count();
    let failed = command(root, &["quant", "--single-fragment-length", "180"]);
    assert!(!failed.status.success());
    assert_eq!(head(&db, OutputGroup::Quant), old);
    assert!(db.verify_generation_head(&old).unwrap());
    assert_eq!(
        db.list_job_attempts(100)
            .unwrap()
            .iter()
            .filter(|a| a.stage == "quant")
            .count(),
        jobs
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
fn competing_recovery_commands_wait_for_one_replacement_producer() {
    let source = Source::new();
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (_paths, db) = setup(root, Some(&source), "none");
    source.delay_ms.store(500, Ordering::SeqCst);
    let mut first = std::process::Command::new(env!("CARGO_BIN_EXE_rnaa"))
        .args([
            "--root",
            root.to_str().unwrap(),
            "quant",
            "--single-fragment-length",
            "180",
        ])
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    let start = Instant::now();
    while source.requests.load(Ordering::SeqCst) < 2 {
        assert!(start.elapsed() < Duration::from_secs(10));
        assert!(first.try_wait().unwrap().is_none());
        thread::sleep(Duration::from_millis(10));
    }
    success(root, &["quant", "--single-fragment-length", "180"]);
    let result = first.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
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
        db.list_job_attempts(100)
            .unwrap()
            .iter()
            .filter(|a| a.stage == "quant")
            .count(),
        2
    );
    assert_eq!(
        db.get_run("RUN1").unwrap().unwrap().state,
        RunState::QuantDone
    );
    assert_eq!(
        db.list_job_attempts(100)
            .unwrap()
            .iter()
            .filter(|attempt| attempt.stage == "download")
            .count(),
        2
    );
}
