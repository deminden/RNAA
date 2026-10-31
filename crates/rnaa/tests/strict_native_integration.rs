//! Strict request validation must check native quantification, not just supplied tables.
#[path = "support/commands.rs"]
mod commands;

use rnaa_core::{
    ArtifactKind, ArtifactRecord, Database, LibraryLayout, ProjectConfig, ProjectPaths, RunRecord,
    RunState,
    util::{compute_sha256, file_size, now_rfc3339},
};
use std::{collections::BTreeMap, fs, io::Write};

#[test]
fn current_request_rejects_changed_native_quantification_settings() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    commands::success(root, &["init"]);
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).unwrap();
    let mut input = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(include_bytes!("../../../fixtures/native_workflow_counts.tsv").as_slice());
    let samples = input
        .headers()
        .unwrap()
        .iter()
        .skip(1)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let rows = input.records().map(Result::unwrap).collect::<Vec<_>>();
    let mut fasta = String::new();
    let mut gtf = String::new();
    let mut reads = vec![Vec::new(); samples.len()];
    for (index, row) in rows.iter().enumerate() {
        let mut state = (index as u64 + 1).wrapping_mul(0x9e3779b97f4a7c15);
        let sequence = (0..1000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize] as char
            })
            .collect::<String>();
        let gene = &row[0];
        fasta.push_str(&format!(">t{gene}\n{sequence}\n"));
        gtf.push_str(&format!("chr1\tsynthetic\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"{gene}\"; transcript_id \"t{gene}\"; gene_name \"{gene}\";\n"));
        for (column, output) in reads.iter_mut().enumerate() {
            let count = 1 + (row[column + 1].parse::<f64>().unwrap() as usize % 12);
            for copy in 0..count {
                writeln!(
                    output,
                    "@g{index}r{copy}\n{}\n+\n{}",
                    &sequence[..100],
                    "I".repeat(100)
                )
                .unwrap();
            }
        }
    }
    let fasta_path = root.join("reference.fa");
    let gtf_path = root.join("reference.gtf");
    fs::write(&fasta_path, fasta).unwrap();
    fs::write(&gtf_path, gtf).unwrap();
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.refs.custom_cdna = fasta_path.display().to_string();
    config.refs.custom_gtf = gtf_path.display().to_string();
    config.quant.preprocess = false;
    config.quant.threads = 1;
    config.quant.worker_budget = 2;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    for (index, sample) in samples.iter().enumerate() {
        db.upsert_run(&RunRecord {
            project_id: db.project_id().into(),
            study_accession: "synthetic".into(),
            source_study_accession: None,
            run_accession: sample.clone(),
            sample_accession: Some(format!("BIO{index}")),
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
        let directory = paths.raw_run_dir(sample);
        fs::create_dir_all(&directory).unwrap();
        let fastq = directory.join(format!("{sample}.fastq"));
        fs::write(&fastq, &reads[index]).unwrap();
        db.record_artifact(&ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some(sample.clone()),
            kind: ArtifactKind::FastqSingle,
            path: fastq.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".into(),
            checksum: compute_sha256(&fastq).unwrap(),
            bytes: file_size(&fastq).unwrap(),
            created_at: now_rfc3339(),
        })
        .unwrap();
    }
    commands::success(root, &["quant", "--fastq-retention", "both"]);
    fs::write(
        paths.samplesheet_path(),
        format!("run_accession\n{}\n", samples.join("\n")),
    )
    .unwrap();
    commands::success(root, &["normalize", "--design", "~ 1"]);
    commands::success(root, &["export", "--what", "counts", "--require-current"]);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.quant.single_fragment_length += 20.0;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let stale = commands::command(root, &["export", "--what", "counts", "--require-current"]);
    assert!(!stale.status.success());
    assert!(
        String::from_utf8_lossy(&stale.stderr)
            .contains("does not match current scientific inputs and settings"),
        "{}",
        String::from_utf8_lossy(&stale.stderr)
    );
    config.quant.single_fragment_length -= 20.0;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let quant = db
        .generation_head(rnaa_core::generations::OutputGroup::Quant, "RUN1")
        .unwrap()
        .unwrap();
    let mut conn = rusqlite::Connection::open(db.path()).unwrap();
    let tx = conn.transaction().unwrap();
    let json: String = tx
        .query_row(
            "SELECT inventory_json FROM artifact_generations WHERE generation_id=?1",
            [&quant.generation_id],
            |row| row.get(0),
        )
        .unwrap();
    let mut seal: serde_json::Value = serde_json::from_str(&json).unwrap();
    seal.as_object_mut().unwrap().remove("lineage_version");
    seal.as_object_mut().unwrap().remove("dependencies");
    fs::write(
        quant.payload_path.parent().unwrap().join("seal.json"),
        serde_json::to_vec_pretty(&seal).unwrap(),
    )
    .unwrap();
    tx.execute(
        "DELETE FROM generation_dependencies WHERE consumer_generation_id=?1",
        [&quant.generation_id],
    )
    .unwrap();
    tx.execute(
        "UPDATE artifact_generations SET inventory_json=?1 WHERE generation_id=?2",
        rusqlite::params![serde_json::to_string(&seal).unwrap(), quant.generation_id],
    )
    .unwrap();
    tx.commit().unwrap();
    let unresolved = commands::command(root, &["export", "--what", "counts", "--require-current"]);
    assert!(!unresolved.status.success());
    let reason = String::from_utf8_lossy(&unresolved.stderr);
    assert!(
        reason.contains("quant generation for RUN1 has unresolved input lineage")
            || reason.contains("selected de generation failed verification"),
        "{reason}"
    );
}
