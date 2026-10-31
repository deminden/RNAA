//! Ready-reference raw-count fixture shared by index and preparation admission checks.
use super::commands::success;
use rnaa_core::{
    Database, LibraryLayout, ManifestArtifact, ProjectConfig, ProjectPaths, RunRecord, RunState,
    model::{ArtifactKind, ArtifactRecord},
    util::now_rfc3339,
};
use std::{collections::BTreeMap, fs, path::Path};
pub fn setup(root: &Path, malformed_last: bool) -> (ProjectPaths, Database) {
    success(root, &["init"]);
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).unwrap();
    let mut random = 59_u64;
    let sequence = (0..1200)
        .map(|_| {
            random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[((random >> 32) & 3) as usize] as char
        })
        .collect::<String>();
    let fasta = root.join("transcripts.fa");
    let annotation = root.join("transcripts.gtf");
    fs::write(&fasta, format!(">tx1\n{sequence}\n")).unwrap();
    fs::write(
        &annotation,
        "chr1\tsynthetic\ttranscript\t1\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\";\n",
    )
    .unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.resources.pool_root = root.join("resources").display().to_string();
    config.refs.custom_cdna = fasta.display().to_string();
    config.refs.custom_gtf = annotation.display().to_string();
    config.refs.index_k = 21;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let read = format!("@read\n{}\n+\n{}\n", &sequence[30..130], "I".repeat(100));
    for run in ["RUN1", "RUN2", "RUN3", "RUN4"] {
        let raw = paths.raw_run_dir(run).join("reads.fastq");
        fs::create_dir_all(raw.parent().unwrap()).unwrap();
        fs::write(
            &raw,
            if malformed_last && run == "RUN4" {
                "@broken\nACGT\n+\nI\n"
            } else {
                &read
            },
        )
        .unwrap();
        db.upsert_run(&RunRecord {
            project_id: db.project_id().into(),
            study_accession: "fixture".into(),
            source_study_accession: None,
            run_accession: run.into(),
            sample_accession: Some(format!("BIO_{run}")),
            experiment_accession: None,
            library_layout: LibraryLayout::Single,
            instrument_platform: None,
            instrument_model: None,
            remote_files: vec![],
            metadata: BTreeMap::new(),
            state: RunState::Verified,
            last_error: None,
            updated_at: now_rfc3339(),
        })
        .unwrap();
        record_raw(&db, run, &raw);
    }
    // These tests constrain quant/preprocessing phases. Reference preparation has
    // its own admission tests and must finish before they hold the shared pool.
    success(root, &["refs", "prepare"]);
    (paths, db)
}
pub fn record_raw(db: &Database, run: &str, path: &Path) {
    let input = ManifestArtifact::snapshot("FASTQ_SINGLE", path).unwrap();
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: Some(run.into()),
        kind: ArtifactKind::FastqSingle,
        path: input.path,
        blob_id: None,
        shared_path: None,
        checksum_type: input.checksum_type,
        checksum: input.checksum,
        bytes: input.bytes,
        created_at: now_rfc3339(),
    })
    .unwrap();
}
