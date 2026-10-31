//! Cohort attempts never replace per-run ingest progress or diagnostics.
mod support;
use rnaa_core::{
    Database, ProjectStageState, RunState,
    generations::OutputGroup,
    model::{ArtifactKind, ArtifactRecord},
    util::{compute_sha256, now_rfc3339},
};
use rusqlite::OptionalExtension;
use std::{
    fs,
    path::Path,
    process::{Child, Command, Stdio},
    time::{Duration, Instant},
};
use support::native::{command, setup, success};

fn snapshot(db: &Database) -> serde_json::Value {
    serde_json::to_value(db.list_runs().unwrap()).unwrap()
}
struct Worker(Child);
impl Drop for Worker {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
fn de(root: &Path) {
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
}
fn corr(root: &Path) {
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
}

#[test]
fn correlation_reuses_without_capacity_and_cancels_while_waiting() {
    use rnaa_core::resources::{Admission, ReservationOwner, ResourcePool, ResourceRequest};
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    de(root);
    corr(root);
    let previous = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let mut config = rnaa_core::ProjectConfig::load(&paths.config_path).unwrap();
    let pool_root = root.join("private-pool");
    config.resources.pool_root = pool_root.display().to_string();
    config.resources.cpu_threads = 1;
    config.resources.memory_mb = 512;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 20;
    config.save(&paths.config_path).unwrap();
    let pool = ResourcePool::open(&pool_root).unwrap();
    let held = match pool
        .try_reserve(
            config.resources.limits().unwrap(),
            ReservationOwner {
                project_root: root.into(),
                attempt_id: "external-blocker".into(),
                stage: "test.blocker".into(),
                run_accession: None,
            },
            ResourceRequest {
                cpu_threads: 1,
                ..Default::default()
            },
        )
        .unwrap()
    {
        Admission::Granted(held) => held,
        Admission::Wait(_) => panic!("unexpected contention"),
    };
    corr(root);
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        previous.generation_id
    );
    let mut worker = Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args([
                "corr",
                "--adjust",
                "none",
                "--geneset",
                "topvar:40",
                "--out",
                "edges:topk=6",
            ])
            .env("RUST_LOG", "info")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .unwrap(),
    );
    let start = Instant::now();
    let attempt = loop {
        let job = db
            .list_job_attempts(100)
            .unwrap()
            .into_iter()
            .find(|job| job.stage == "corr" && job.state == "RUNNING");
        let waiting = fs::read_to_string(paths.logs_dir.join("corr.log"))
            .unwrap_or_default()
            .contains("corr.compute")
            && fs::read_to_string(paths.logs_dir.join("corr.log"))
                .unwrap_or_default()
                .contains("waiting for shared resource capacity");
        if let Some(job) = job
            && waiting
        {
            break job.attempt_id;
        }
        assert!(
            worker.0.try_wait().unwrap().is_none(),
            "correlation exited before waiting"
        );
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "correlation did not reach admission"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(db.cancel_job(&attempt).unwrap());
    let start = Instant::now();
    let status = loop {
        if let Some(status) = worker.0.try_wait().unwrap() {
            break status;
        }
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "cancelled correlation did not stop"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(!status.success());
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        previous.generation_id
    );
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    drop(held);
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    success(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=6",
        ],
    );
    assert_ne!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        previous.generation_id
    );
    assert!(pool.snapshot().unwrap().reservations.is_empty());
}

#[test]
fn correlation_network_waits_for_memory_without_holding_cpu_and_cancels() {
    use rnaa_core::resources::{Admission, ReservationOwner, ResourcePool, ResourceRequest};
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    de(root);
    corr(root);
    let previous = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let mut config = rnaa_core::ProjectConfig::load(&paths.config_path).unwrap();
    let pool_root = root.join("private-pool");
    config.resources.pool_root = pool_root.display().to_string();
    config.resources.cpu_threads = 1;
    config.resources.memory_mb = 512;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 20;
    config.corr.block_rows = 1024;
    config.save(&paths.config_path).unwrap();
    let pool = ResourcePool::open(&pool_root).unwrap();
    let held = match pool
        .try_reserve(
            config.resources.limits().unwrap(),
            ReservationOwner {
                project_root: root.into(),
                attempt_id: "memory-blocker".into(),
                stage: "test.blocker".into(),
                run_accession: None,
            },
            ResourceRequest {
                memory_bytes: 444 * 1024 * 1024,
                ..Default::default()
            },
        )
        .unwrap()
    {
        Admission::Granted(held) => held,
        Admission::Wait(_) => panic!("unexpected contention"),
    };
    let args = [
        "corr",
        "--adjust",
        "none",
        "--geneset",
        "topvar:400",
        "--out",
        "dense",
    ];
    let mut worker = Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args(args)
            .env("RUST_LOG", "info")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .unwrap(),
    );
    let start = Instant::now();
    let attempt = loop {
        let job = db
            .list_job_attempts(100)
            .unwrap()
            .into_iter()
            .find(|job| job.stage == "corr" && job.state == "RUNNING");
        let log = fs::read_to_string(paths.logs_dir.join("corr.log")).unwrap_or_default();
        if let Some(job) = job
            && log.lines().any(|line| {
                line.contains("corr.network")
                    && line.contains("waiting for shared resource capacity")
            })
        {
            break job.attempt_id;
        }
        assert!(
            worker.0.try_wait().unwrap().is_none(),
            "correlation exited before waiting for memory"
        );
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "correlation did not wait for memory"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    let waiting = pool.snapshot().unwrap();
    assert_eq!(waiting.reservations.len(), 1);
    assert_eq!(waiting.cpu_threads, 0, "waiting correlation retained CPU");
    assert_eq!(waiting.memory_bytes, 444 * 1024 * 1024);
    assert!(db.cancel_job(&attempt).unwrap());
    let start = Instant::now();
    let status = loop {
        if let Some(status) = worker.0.try_wait().unwrap() {
            break status;
        }
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "cancelled correlation did not stop"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(!status.success());
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        previous.generation_id
    );
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    drop(held);
    success(root, &args);
    let fresh = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    assert_ne!(fresh.generation_id, previous.generation_id);
    assert!(pool.snapshot().unwrap().reservations.is_empty());
    let mut config = rnaa_core::ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.memory_mb = 1;
    config.save(&paths.config_path).unwrap();
    success(root, &args);
    let rejected = command(
        root,
        &[
            "corr",
            "--adjust",
            "none",
            "--geneset",
            "topvar:40",
            "--out",
            "edges:topk=7",
        ],
    );
    assert!(!rejected.status.success());
    assert!(
        String::from_utf8_lossy(&rejected.stderr)
            .contains("resource request exceeds its own CPU or memory capacity")
    );
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        fresh.generation_id
    );
    assert!(pool.snapshot().unwrap().reservations.is_empty());
}

#[test]
fn correlation_adjustment_waits_for_memory_without_holding_cpu_and_cancels() {
    use rnaa_core::resources::{Admission, ReservationOwner, ResourcePool, ResourceRequest};
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    de(root);
    corr(root);
    let previous = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let mut config = rnaa_core::ProjectConfig::load(&paths.config_path).unwrap();
    let pool_root = root.join("private-pool");
    config.resources.pool_root = pool_root.display().to_string();
    config.resources.cpu_threads = 1;
    config.resources.memory_mb = 512;
    config.resources.disk_headroom_mb = 0;
    config.resources.wait_seconds = 20;
    config.save(&paths.config_path).unwrap();
    let pool = ResourcePool::open(&pool_root).unwrap();
    let held = match pool
        .try_reserve(
            config.resources.limits().unwrap(),
            ReservationOwner {
                project_root: root.into(),
                attempt_id: "memory-blocker".into(),
                stage: "test.blocker".into(),
                run_accession: None,
            },
            ResourceRequest {
                memory_bytes: 500 * 1024 * 1024,
                ..Default::default()
            },
        )
        .unwrap()
    {
        Admission::Granted(held) => held,
        Admission::Wait(_) => panic!("unexpected contention"),
    };
    let args = [
        "corr",
        "--adjust",
        "residualize",
        "--model",
        "~ condition",
        "--geneset",
        "topvar:40",
        "--out",
        "edges:topk=6",
    ];
    let mut worker = Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args(args)
            .env("RUST_LOG", "info")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .unwrap(),
    );
    let start = Instant::now();
    let attempt = loop {
        let job = db
            .list_job_attempts(100)
            .unwrap()
            .into_iter()
            .find(|job| job.stage == "corr" && job.state == "RUNNING");
        let log = fs::read_to_string(paths.logs_dir.join("corr.log")).unwrap_or_default();
        if let Some(job) = job
            && log.lines().any(|line| {
                line.contains("corr.adjust")
                    && line.contains("waiting for shared resource capacity")
            })
        {
            break job.attempt_id;
        }
        assert!(
            worker.0.try_wait().unwrap().is_none(),
            "correlation exited before waiting for adjustment memory"
        );
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "correlation did not wait for adjustment memory"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    let waiting = pool.snapshot().unwrap();
    assert_eq!(waiting.reservations.len(), 1);
    assert_eq!(waiting.cpu_threads, 0);
    assert!(db.cancel_job(&attempt).unwrap());
    let start = Instant::now();
    let status = loop {
        if let Some(status) = worker.0.try_wait().unwrap() {
            break status;
        }
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "cancelled adjustment did not stop"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(!status.success());
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        previous.generation_id
    );
    assert_eq!(pool.snapshot().unwrap().reservations.len(), 1);
    drop(held);
    success(root, &args);
    assert_ne!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        previous.generation_id
    );
    assert!(pool.snapshot().unwrap().reservations.is_empty());
}

#[test]
fn cancelling_running_correlation_preserves_the_published_head_and_releases_claims() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let (paths, db) = setup(root);
    de(root);
    corr(root);
    let previous = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    let mut config = rnaa_core::ProjectConfig::load(&paths.config_path).unwrap();
    let pool_root = root.join("private-pool");
    config.resources.pool_root = pool_root.display().to_string();
    config.resources.cpu_threads = 1;
    config.resources.memory_mb = 512;
    config.resources.disk_headroom_mb = 0;
    config.corr.threads = 1;
    config.corr.block_rows = 64;
    config.save(&paths.config_path).unwrap();
    let input = root.join("long-running.tsv");
    let mut matrix = String::from("gene_id\tRUN1\tRUN2\tRUN3\tRUN4\tRUN5\tRUN6\tRUN7\tRUN8\n");
    for gene in 0..4000 {
        matrix.push_str(&format!("G{gene}"));
        for sample in 0..8 {
            matrix.push_str(&format!(
                "\t{}",
                (gene * 13 + sample * sample + gene * sample) % 97
            ));
        }
        matrix.push('\n');
    }
    fs::write(&input, matrix).unwrap();
    let pool = rnaa_core::resources::ResourcePool::open(&pool_root).unwrap();
    let mut worker = Worker(
        Command::new(env!("CARGO_BIN_EXE_rnaa"))
            .arg("--root")
            .arg(root)
            .args([
                "corr",
                "--matrix",
                input.to_str().unwrap(),
                "--adjust",
                "none",
                "--geneset",
                "all",
                "--method",
                "kendall",
                "--out",
                "edges:topk=5",
            ])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .unwrap(),
    );
    let start = Instant::now();
    let attempt = loop {
        let job = db
            .list_job_attempts(100)
            .unwrap()
            .into_iter()
            .find(|job| job.stage == "corr" && job.state == "RUNNING");
        if let Some(job) = job {
            let connection = rusqlite::Connection::open(db.path()).unwrap();
            let payload: Option<String> = connection.query_row(
                "SELECT payload_path FROM artifact_generations WHERE attempt_id=?1 AND state='BUILDING'",
                [&job.attempt_id], |row| row.get(0),
            ).optional().unwrap();
            let active = pool.snapshot().unwrap().reservations.iter().any(|claim| {
                claim.owner.attempt_id == job.attempt_id && claim.owner.stage == "corr.network"
            });
            let writing = payload.is_some_and(|path| {
                std::fs::metadata(
                    Path::new(&path)
                        .join("corr")
                        .join(db.project_id())
                        .join("edges_topk.tsv"),
                )
                .is_ok_and(|metadata| metadata.len() > 0)
            });
            if active && writing {
                break job.attempt_id;
            }
        }
        assert!(
            worker.0.try_wait().unwrap().is_none(),
            "correlation exited before producing a network tile"
        );
        assert!(
            start.elapsed() < Duration::from_secs(30),
            "correlation did not begin writing network output"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(db.cancel_job(&attempt).unwrap());
    let start = Instant::now();
    let status = loop {
        if let Some(status) = worker.0.try_wait().unwrap() {
            break status;
        }
        assert!(
            start.elapsed() < Duration::from_secs(15),
            "running correlation did not observe cancellation"
        );
        std::thread::sleep(Duration::from_millis(10));
    };
    assert!(!status.success());
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        previous.generation_id
    );
    assert!(pool.snapshot().unwrap().reservations.is_empty());
}

#[test]
fn successful_analysis_and_historical_reuse_preserve_complete_run_records() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path());
    db.set_run_state(
        "RUN1",
        RunState::QuantDone,
        Some("retained ingest diagnostic"),
    )
    .unwrap();
    let original = snapshot(&db);
    de(root.path());
    assert_eq!(snapshot(&db), original);
    let first = db
        .generation_head(OutputGroup::De, db.project_id())
        .unwrap()
        .unwrap();
    de(root.path());
    assert_eq!(snapshot(&db), original);
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        first.generation_id
    );
    corr(root.path());
    assert_eq!(snapshot(&db), original);
    let first = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    corr(root.path());
    assert_eq!(snapshot(&db), original);
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        first.generation_id
    );
    let status = command(root.path(), &["status"]);
    assert!(status.status.success());
    let text = String::from_utf8_lossy(&status.stdout);
    assert!(text.contains("run_progress_model\tingest_and_quant_only"));
    assert!(text.contains("progress_quant\t8/8"));
    assert!(text.contains("cohort_stage\tdeseq2\tDONE"));
    assert!(text.contains("cohort_stage\tcorr\tDONE"));
    let json_status = command(root.path(), &["status", "--json"]);
    assert!(
        json_status.status.success(),
        "{}",
        String::from_utf8_lossy(&json_status.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&json_status.stdout).unwrap();
    assert_eq!(report["schema_version"], 1);
    assert_eq!(report["cohort"]["ready"], true);
    assert_eq!(report["progress"]["total_runs"], 8);
    assert_eq!(report["progress"]["quant_done"], 8);
    assert_eq!(report["progress"]["normalize_done"], true);
    assert_eq!(report["progress"]["correlation_done"], true);
    assert!(
        report["stages"]
            .as_array()
            .unwrap()
            .iter()
            .any(|s| s["stage"] == "corr" && s["state"] == "DONE")
    );
    assert_eq!(snapshot(&db), original);
    assert!(
        db.list_job_attempts(100)
            .unwrap()
            .iter()
            .filter(|job| matches!(job.stage.as_str(), "deseq2" | "corr"))
            .all(|job| job.run_accession.is_none() && job.request_sha256.is_some())
    );
    assert!(
        paths
            .with_generation_heads(&db)
            .unwrap()
            .de_project_dir(db.project_id())
            .join("vst.tsv")
            .is_file()
    );
}

#[test]
fn publication_failures_preserve_runs_previous_heads_and_ingest_diagnostics() {
    let root = tempfile::tempdir().unwrap();
    let (paths, db) = setup(root.path());
    db.set_run_state("RUN1", RunState::QuantDone, Some("ingest note"))
        .unwrap();
    de(root.path());
    corr(root.path());
    let original = snapshot(&db);
    let heads = db.generation_heads().unwrap();
    let published = paths.with_generation_heads(&db).unwrap();
    let matrix = published.de_project_dir(db.project_id()).join("vst.tsv");
    let hash = compute_sha256(&matrix).unwrap();
    for (group, args, stage) in [
        ("de", vec!["normalize"], "normalize"),
        ("corr", vec!["corr", "--out", "dense"], "corr"),
    ] {
        let connection = rusqlite::Connection::open(db.path()).unwrap();
        connection.execute_batch(&format!("CREATE TRIGGER reject_analysis_head BEFORE INSERT ON generation_heads WHEN NEW.output_group='{group}' BEGIN SELECT RAISE(ABORT,'injected analysis publication failure'); END;")).unwrap();
        let result = command(root.path(), &args);
        assert!(!result.status.success());
        assert!(
            String::from_utf8_lossy(&result.stderr)
                .contains("injected analysis publication failure")
        );
        assert_eq!(snapshot(&db), original);
        assert_eq!(db.generation_heads().unwrap(), heads);
        assert_eq!(compute_sha256(&matrix).unwrap(), hash);
        assert_eq!(
            db.get_project_stage_state(stage).unwrap().unwrap().state,
            ProjectStageState::Failed
        );
        assert!(
            db.list_job_attempts(100)
                .unwrap()
                .iter()
                .any(|job| job.stage == stage && job.state == "FAILED")
        );
        connection
            .execute_batch("DROP TRIGGER reject_analysis_head")
            .unwrap();
    }
    corr(root.path());
    assert_eq!(snapshot(&db), original);
}

#[test]
fn explicit_matrix_uses_quantified_members_without_fake_de_progress() {
    let root = tempfile::tempdir().unwrap();
    let (_paths, db) = setup(root.path());
    let original = snapshot(&db);
    let input =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/native_workflow_counts.tsv");
    success(
        root.path(),
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
    assert_eq!(snapshot(&db), original);
    assert!(
        db.generation_head(OutputGroup::De, db.project_id())
            .unwrap()
            .is_none()
    );
    assert_eq!(
        db.get_project_stage_state("corr").unwrap().unwrap().state,
        ProjectStageState::Done
    );
}

#[test]
fn custom_matrix_reuse_and_invalid_replacement_preserve_the_previous_result() {
    let root = tempfile::tempdir().unwrap();
    let (_, db) = setup(root.path());
    let matrix = root.path().join("custom.tsv");
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/native_workflow_counts.tsv"),
        &matrix,
    )
    .unwrap();
    let args = [
        "corr",
        "--matrix",
        matrix.to_str().unwrap(),
        "--adjust",
        "none",
        "--geneset",
        "all",
        "--out",
        "edges:topk=5",
    ];
    success(root.path(), &args);
    let first = db
        .generation_head(OutputGroup::Corr, db.project_id())
        .unwrap()
        .unwrap();
    success(root.path(), &args);
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        first.generation_id
    );
    let contents = fs::read_to_string(&matrix).unwrap();
    let (header, data) = contents.split_once('\n').unwrap();
    let (gene, rest) = data.split_once('\t').unwrap();
    let (_, remaining) = rest.split_once('\t').unwrap();
    fs::write(&matrix, format!("{header}\n{gene}\tinvalid\t{remaining}")).unwrap();
    let result = command(root.path(), &args);
    assert!(!result.status.success());
    assert_eq!(
        db.generation_head(OutputGroup::Corr, db.project_id())
            .unwrap()
            .unwrap()
            .generation_id,
        first.generation_id
    );
}

#[test]
fn narrower_cohort_progress_excludes_history_and_does_not_mark_analysis_complete() {
    let root = tempfile::tempdir().unwrap();
    let (_paths, db) = setup(root.path());
    de(root.path());
    corr(root.path());
    let original = snapshot(&db);
    let runs = db.list_runs().unwrap();
    let selectors = db.selector_contract().unwrap();
    let previous = db.active_cohort().unwrap();
    db.publish_resolved_cohort(
        &runs[..6],
        &selectors,
        previous.as_ref().map(|cohort| cohort.cohort_id.as_str()),
    )
    .unwrap();
    let raw = root.path().join("excluded.fastq");
    fs::write(&raw, b"excluded history").unwrap();
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: Some("RUN8".into()),
        kind: ArtifactKind::FastqSingle,
        path: raw.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".into(),
        checksum: compute_sha256(&raw).unwrap(),
        bytes: 16,
        created_at: now_rfc3339(),
    })
    .unwrap();
    let status = command(root.path(), &["status"]);
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let text = String::from_utf8_lossy(&status.stdout);
    for value in [
        "runs_catalog_total\t8",
        "runs_total\t6",
        "progress_download\t6/6",
        "progress_quant\t6/6",
        "downloaded_bytes\t0",
        "progress_deseq2\t0/1",
        "progress_corr\t0/1",
    ] {
        assert!(text.contains(value), "missing {value}: {text}");
    }
    let output = command(root.path(), &["status", "--json"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["runs_catalog_total"], 8);
    assert_eq!(report["progress"]["total_runs"], 6);
    assert_eq!(report["progress"]["quant_done"], 6);
    assert_eq!(report["progress"]["normalize_done"], false);
    assert_eq!(report["active_runs"].as_array().unwrap().len(), 6);
    assert_eq!(snapshot(&db), original);
}

#[test]
fn analysis_progress_crash_worker() {
    let Some(root) = std::env::var_os("RNAA_ANALYSIS_PROGRESS_ROOT") else {
        return;
    };
    let stage = std::env::var("RNAA_ANALYSIS_PROGRESS_STAGE").unwrap();
    let root = Path::new(&root);
    let db = Database::open(root).unwrap();
    let lease = rnaa_core::jobs::JobLease::try_acquire(&db, &stage, None, &[])
        .unwrap()
        .unwrap();
    lease
        .database()
        .set_project_stage_state(&stage, ProjectStageState::Running, None)
        .unwrap();
    rnaa_core::util::publish_file(&root.join("analysis-ready"), |staged| {
        fs::write(staged, lease.attempt_id())?;
        Ok(())
    })
    .unwrap();
    loop {
        std::thread::sleep(Duration::from_millis(50));
    }
}
#[test]
fn killed_cohort_workers_recover_without_resetting_runs_or_live_ingest() {
    let root = tempfile::tempdir().unwrap();
    let (_paths, db) = setup(root.path());
    let mut ingest = rnaa_core::jobs::JobLease::try_acquire(
        &db,
        "download",
        Some("RUN8"),
        &[RunState::QuantDone],
    )
    .unwrap()
    .unwrap();
    let original = snapshot(&db);
    for stage in ["normalize", "deseq2", "corr"] {
        let marker = root.path().join("analysis-ready");
        let _ = fs::remove_file(&marker);
        let mut child = Worker(
            Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "analysis_progress_crash_worker", "--nocapture"])
                .env("RNAA_ANALYSIS_PROGRESS_ROOT", root.path())
                .env("RNAA_ANALYSIS_PROGRESS_STAGE", stage)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .spawn()
                .unwrap(),
        );
        let deadline = Instant::now() + Duration::from_secs(10);
        while !marker.exists() {
            assert!(child.0.try_wait().unwrap().is_none(), "worker exited");
            if Instant::now() >= deadline {
                panic!("worker did not claim analysis");
            }
            std::thread::sleep(Duration::from_millis(10));
        }
        assert_eq!(snapshot(&db), original);
        assert_eq!(
            db.get_project_stage_state(stage).unwrap().unwrap().state,
            ProjectStageState::Running
        );
        child.0.kill().unwrap();
        child.0.wait().unwrap();
        let attempt = fs::read_to_string(&marker).unwrap();
        assert!(!attempt.is_empty(), "worker published an empty attempt ID");
        rusqlite::Connection::open(db.path())
            .unwrap()
            .execute(
                "UPDATE job_attempts SET lease_expires_ms=0 WHERE attempt_id=?1",
                [attempt],
            )
            .unwrap();
        assert!(db.recover_jobs().unwrap() > 0);
        assert_eq!(snapshot(&db), original);
        assert_eq!(
            db.get_project_stage_state(stage).unwrap().unwrap().state,
            ProjectStageState::Failed
        );
        assert_eq!(
            db.get_run("RUN8").unwrap().unwrap().state,
            RunState::Downloading
        );
    }
    ingest
        .finish(rnaa_core::jobs::JobOutcome::Failed, Some("test finished"))
        .unwrap();
}
