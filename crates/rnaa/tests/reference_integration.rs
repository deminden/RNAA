use rnaa_core::{Database, ProjectConfig, ProjectPaths, generations::OutputGroup};
use std::{
    fs,
    path::Path,
    process::{Command, Output, Stdio},
};
fn command(root: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rnaa"))
        .arg("--root")
        .arg(root)
        .args(args)
        .output()
        .unwrap()
}
fn setup(root: &Path) -> (Database, ProjectPaths) {
    let init = command(root, &["init"]);
    assert!(init.status.success());
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).unwrap();
    let cdna = root.join("source.fa");
    let gtf = root.join("source.gtf");
    fs::write(&cdna,format!(">tx1\n{}\n","ACGTTCGATCGGATCGATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCTAGCTGACTGATCGATCGTAGCTAGCTGATCGTAGCTAGCTAGCATCGA")).unwrap();
    fs::write(&gtf,"chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\"; gene_name \"Gene1\";\n").unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.refs.custom_cdna = cdna.display().to_string();
    config.refs.custom_gtf = gtf.display().to_string();
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    (db, paths)
}
#[test]
fn competing_reference_processes_publish_one_complete_bundle_and_reuse_it() {
    let root = tempfile::tempdir().unwrap();
    let (db, _paths) = setup(root.path());
    let children = (0..2)
        .map(|_| {
            Command::new(env!("CARGO_BIN_EXE_rnaa"))
                .arg("--root")
                .arg(root.path())
                .args(["refs", "prepare"])
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap()
        })
        .collect::<Vec<_>>();
    let results = children
        .into_iter()
        .map(|child| child.wait_with_output().unwrap())
        .collect::<Vec<_>>();
    for result in &results {
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    let manifests = results
        .iter()
        .map(|result| {
            String::from_utf8_lossy(&result.stdout)
                .lines()
                .find_map(|line| line.strip_prefix("manifest\t").map(str::to_string))
                .unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(manifests[0], manifests[1]);
    let heads = db.generation_heads().unwrap();
    assert_eq!(heads.len(), 1);
    assert_eq!(heads[0].output_group, OutputGroup::Refs);
    assert!(db.verify_generation_head(&heads[0]).unwrap());
    assert_eq!(db.list_artifacts().unwrap().len(), 6);
    assert_eq!(
        db.list_events(100)
            .unwrap()
            .iter()
            .filter(|event| event.message == "reference bundle prepared")
            .count(),
        1
    );
    let before = fs::read(&manifests[0]).unwrap();
    let jobs = db.list_job_attempts(100).unwrap().len();
    let cached = command(root.path(), &["refs", "prepare"]);
    assert!(cached.status.success());
    assert_eq!(fs::read(&manifests[0]).unwrap(), before);
    assert_eq!(db.list_job_attempts(100).unwrap().len(), jobs);
}
#[test]
fn a_bad_reference_source_fails_without_replacing_registered_good_outputs() {
    let root = tempfile::tempdir().unwrap();
    let (db, paths) = setup(root.path());
    let valid = command(root.path(), &["refs", "prepare"]);
    assert!(valid.status.success());
    let previous = db.generation_heads().unwrap().remove(0);
    let index = previous
        .paths(&paths)
        .unwrap()
        .reference_dir(&previous.subject)
        .join("kallisto.index");
    let bytes = fs::read(&index).unwrap();
    fs::write(root.path().join("source.gtf"), "invalid annotation").unwrap();
    let invalid = command(root.path(), &["refs", "prepare"]);
    assert!(!invalid.status.success());
    assert_eq!(db.generation_heads().unwrap(), vec![previous.clone()]);
    assert!(db.verify_generation_head(&previous).unwrap());
    assert_eq!(fs::read(index).unwrap(), bytes);
    assert_eq!(db.list_artifacts().unwrap().len(), 6);
    assert_eq!(db.list_job_attempts(1).unwrap()[0].state, "FAILED");
}

#[test]
fn annotation_changes_reuse_pruned_native_counts_but_index_changes_block_analysis() {
    use rnaa_core::{
        ArtifactKind, ArtifactRecord, LibraryLayout, RunRecord, RunState,
        util::{compute_sha256, file_size, now_rfc3339},
    };
    use std::collections::BTreeMap;
    let root = tempfile::tempdir().unwrap();
    let (db, paths) = setup(root.path());
    db.upsert_run(&RunRecord {
        project_id: db.project_id().into(),
        study_accession: "synthetic".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: Vec::new(),
        metadata: BTreeMap::from([("condition".into(), "A".into())]),
        state: RunState::Verified,
        last_error: None,
        updated_at: now_rfc3339(),
    })
    .unwrap();
    let raw = paths.raw_run_dir("RUN1").join("run.fastq");
    fs::create_dir_all(raw.parent().unwrap()).unwrap();
    let sequence = fs::read_to_string(root.path().join("source.fa"))
        .unwrap()
        .lines()
        .nth(1)
        .unwrap()
        .to_string();
    fs::write(
        &raw,
        format!("@read\n{sequence}\n+\n{}\n", "I".repeat(sequence.len())).repeat(3),
    )
    .unwrap();
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: Some("RUN1".into()),
        kind: ArtifactKind::FastqSingle,
        path: raw.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".into(),
        checksum: compute_sha256(&raw).unwrap(),
        bytes: file_size(&raw).unwrap(),
        created_at: now_rfc3339(),
    })
    .unwrap();
    fs::write(
        paths.samplesheet_path(),
        "run_accession\tcondition\nRUN1\tA\n",
    )
    .unwrap();
    let result = command(
        root.path(),
        &["quant", "--worker-budget", "1", "--fastq-retention", "none"],
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!raw.exists());
    let count_head = db
        .generation_head(OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    let manifest = paths
        .with_generation_heads(&db)
        .unwrap()
        .quant_run_dir("RUN1", "kallistors")
        .join("rnaa_quant_manifest.json");
    let bytes = fs::read(&manifest).unwrap();
    let before = db.selected_reference().unwrap().unwrap();
    let quant_jobs = db
        .list_job_attempts(100)
        .unwrap()
        .iter()
        .filter(|job| job.stage == "quant")
        .count();
    let gtf = root.path().join("source.gtf");
    let annotation = fs::read_to_string(&gtf).unwrap();
    fs::write(&gtf, annotation.replace("Gene1", "NewName")).unwrap();
    let prepared = command(root.path(), &["refs", "prepare"]);
    assert!(
        prepared.status.success(),
        "{}",
        String::from_utf8_lossy(&prepared.stderr)
    );
    let after = db.selected_reference().unwrap().unwrap();
    assert_eq!(before.index.checksum, after.index.checksum);
    assert_ne!(before.annotation.checksum, after.annotation.checksum);
    let reused = command(root.path(), &["quant"]);
    assert!(
        reused.status.success(),
        "{}",
        String::from_utf8_lossy(&reused.stderr)
    );
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        count_head
    );
    assert_eq!(fs::read(&manifest).unwrap(), bytes);
    assert_eq!(
        db.list_job_attempts(100)
            .unwrap()
            .iter()
            .filter(|job| job.stage == "quant")
            .count(),
        quant_jobs
    );
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.refs.index_k = 29;
    config.save(&paths.config_path).unwrap();
    db.set_project_stage_state("normalize", rnaa_core::ProjectStageState::Done, None)
        .unwrap();
    let status = command(root.path(), &["status"]);
    assert!(status.status.success());
    assert!(String::from_utf8_lossy(&status.stdout).contains("reference_ready\tno"));
    assert_eq!(
        db.get_project_stage_state("normalize")
            .unwrap()
            .unwrap()
            .state,
        rnaa_core::ProjectStageState::Pending
    );
    let prepared = command(root.path(), &["refs", "prepare"]);
    assert!(
        prepared.status.success(),
        "{}",
        String::from_utf8_lossy(&prepared.stderr)
    );
    assert_ne!(
        db.selected_reference().unwrap().unwrap().index.checksum,
        after.index.checksum
    );
    let analysis = command(root.path(), &["normalize", "--design", "~ condition"]);
    assert!(!analysis.status.success());
    let run = db.get_run("RUN1").unwrap().unwrap();
    assert_eq!(run.state, RunState::QuantFailed);
    assert!(run.last_error.unwrap().contains("current reference"));
    assert_eq!(
        db.generation_head(OutputGroup::Quant, "RUN1")
            .unwrap()
            .unwrap(),
        count_head
    );
    assert_eq!(fs::read(manifest).unwrap(), bytes);
}

fn prepare_reference(root: &Path) {
    let result = command(root, &["refs", "prepare"]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}
fn selected(db: &Database) -> rnaa_core::generations::GenerationHead {
    db.selected_reference().unwrap().unwrap().generation
}
fn index_k(paths: &ProjectPaths, db: &Database, k: usize) {
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.refs.index_k = k;
    config.refs.index_threads = 1;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
}
#[test]
fn restored_index_options_select_original_reference_without_building_another_payload() {
    let root = tempfile::tempdir().unwrap();
    let (db, paths) = setup(root.path());
    index_k(&paths, &db, 31);
    prepare_reference(root.path());
    let first = selected(&db);
    let seal = first.payload_path.parent().unwrap().join("seal.json");
    let bytes = fs::read(&seal).unwrap();
    index_k(&paths, &db, 29);
    prepare_reference(root.path());
    let second = selected(&db);
    assert_ne!(first, second);
    index_k(&paths, &db, 31);
    prepare_reference(root.path());
    assert_eq!(selected(&db), first);
    assert_eq!(
        db.generation_head(OutputGroup::Refs, &first.subject)
            .unwrap()
            .unwrap(),
        first
    );
    assert_eq!(
        db.committed_generation_history(OutputGroup::Refs, &first.subject)
            .unwrap()
            .len(),
        2
    );
    assert_eq!(fs::read(seal).unwrap(), bytes);
    assert!(db.verify_generation_head(&second).unwrap());
    let attempts = db.list_job_attempts(100).unwrap();
    assert_eq!(
        attempts
            .iter()
            .filter(|attempt| attempt.stage == "refs")
            .count(),
        3
    );
    assert!(attempts.iter().all(|attempt| attempt.state == "SUCCEEDED"));
    prepare_reference(root.path());
    assert_eq!(db.list_job_attempts(100).unwrap().len(), attempts.len());
}
#[test]
fn returning_to_another_retained_source_records_an_owned_reference_selection() {
    let root = tempfile::tempdir().unwrap();
    let (db, paths) = setup(root.path());
    prepare_reference(root.path());
    let first = selected(&db);
    let gtf = root.path().join("source.gtf");
    let original = fs::read_to_string(&gtf).unwrap();
    fs::write(&gtf, original.replace("Gene1", "GeneTwo")).unwrap();
    prepare_reference(root.path());
    let second = selected(&db);
    assert_ne!(first.subject, second.subject);
    fs::write(&gtf, &original).unwrap();
    prepare_reference(root.path());
    assert_eq!(selected(&db), first);
    let latest = db.list_job_attempts(1).unwrap().remove(0);
    assert_eq!(latest.stage, "refs");
    assert_eq!(latest.state, "SUCCEEDED");
    assert_ne!(latest.attempt_id, first.attempt_id);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Refs, &first.subject)
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        db.committed_generation_history(OutputGroup::Refs, &second.subject)
            .unwrap()
            .len(),
        1
    );
    assert!(db.verify_generation_head(&second).unwrap());
    assert!(
        db.selected_reference()
            .unwrap()
            .unwrap()
            .generation
            .payload_path
            .starts_with(fs::canonicalize(&paths.generations_dir).unwrap())
    );
}
#[test]
fn corrupt_historical_reference_is_preserved_and_rebuilt_in_a_fresh_generation() {
    let root = tempfile::tempdir().unwrap();
    let (db, paths) = setup(root.path());
    index_k(&paths, &db, 31);
    prepare_reference(root.path());
    let first = selected(&db);
    index_k(&paths, &db, 29);
    prepare_reference(root.path());
    let second = selected(&db);
    let corrupt = first
        .paths(&paths)
        .unwrap()
        .reference_dir(&first.subject)
        .join("kallisto.index");
    fs::write(&corrupt, "corrupt historical index").unwrap();
    index_k(&paths, &db, 31);
    prepare_reference(root.path());
    let replacement = selected(&db);
    assert_ne!(replacement, first);
    assert_ne!(replacement, second);
    assert_eq!(replacement.request_sha256, first.request_sha256);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Refs, &first.subject)
            .unwrap()
            .len(),
        3
    );
    assert!(db.verify_generation_head(&replacement).unwrap());
    assert!(db.verify_generation_head(&second).unwrap());
    assert_eq!(
        fs::read_to_string(corrupt).unwrap(),
        "corrupt historical index"
    );
}
#[test]
fn moved_identical_reference_sources_retain_the_payload_and_record_current_source_observations() {
    let root = tempfile::tempdir().unwrap();
    let (db, paths) = setup(root.path());
    prepare_reference(root.path());
    let first = selected(&db);
    let seal = first.payload_path.parent().unwrap().join("seal.json");
    let bytes = fs::read(&seal).unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    let cdna = root.path().join("renamed.fa");
    let gtf = root.path().join("renamed.gtf");
    fs::rename(&config.refs.custom_cdna, &cdna).unwrap();
    fs::rename(&config.refs.custom_gtf, &gtf).unwrap();
    config.refs.custom_cdna = cdna.display().to_string();
    config.refs.custom_gtf = gtf.display().to_string();
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    prepare_reference(root.path());
    assert_eq!(selected(&db), first);
    assert_eq!(fs::read(seal).unwrap(), bytes);
    assert_eq!(
        db.committed_generation_history(OutputGroup::Refs, &first.subject)
            .unwrap()
            .len(),
        1
    );
    let event = db
        .list_events(100)
        .unwrap()
        .into_iter()
        .find(|event| event.message == "reference bundle selected from verified history")
        .unwrap();
    assert_eq!(
        event.context["source_observations"]["transcriptome"]["location"],
        config.refs.custom_cdna
    );
    db.selected_reference()
        .unwrap()
        .unwrap()
        .verify_current(
            &db,
            &config,
            rnaa_core::references::ReferenceScope::Analysis,
        )
        .unwrap();
}
