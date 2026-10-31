//! Successful native analysis from deterministic synthetic abundance inputs.
use flate2::{Compression, write::GzEncoder};
use rnaa_core::state::RunState;
use rnaa_core::traits::ReferenceManager;
use rnaa_core::util::now_rfc3339;
use rnaa_core::{Database, LibraryLayout, ProjectConfig, ProjectPaths, RunRecord};
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::Path;

pub use super::commands::{command, success};
pub fn gzip(path: &Path, text: &str) {
    let mut encoder = GzEncoder::new(fs::File::create(path).unwrap(), Compression::default());
    encoder.write_all(text.as_bytes()).unwrap();
    encoder.finish().unwrap();
}
pub fn setup(root: &Path) -> (ProjectPaths, Database) {
    success(root, &["init"]);
    let paths = ProjectPaths::new(root);
    let db = Database::open(root).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(include_bytes!("../../../../fixtures/native_workflow_counts.tsv").as_slice());
    let samples = reader
        .headers()
        .unwrap()
        .iter()
        .skip(1)
        .map(str::to_string)
        .collect::<Vec<_>>();
    let rows = reader.records().map(|row| row.unwrap()).collect::<Vec<_>>();
    let names = rows
        .iter()
        .map(|row| format!("t{}.1|{}.1|label|", &row[0], &row[0]))
        .collect::<Vec<_>>();
    let mut fasta = String::new();
    let mut gtf = String::new();
    for row in &rows {
        let gene = &row[0];
        fasta.push_str(&format!(
            ">t{gene}.1|{gene}.1|label|\n{}\n",
            "ACGT".repeat(250)
        ));
        gtf.push_str(&format!("chr1\tsynthetic\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"{gene}.1\"; transcript_id \"t{gene}.1\"; gene_name \"Symbol{gene}\"; gene_type \"protein_coding\";\n"));
    }
    let cdna = root.join("synthetic.fa.gz");
    let annotation = root.join("synthetic.gtf.gz");
    gzip(&cdna, &fasta);
    gzip(&annotation, &gtf);
    let mut config = ProjectConfig::load(&paths.config_path).unwrap();
    config.refs.custom_cdna = cdna.display().to_string();
    config.refs.custom_gtf = annotation.display().to_string();
    config.deseq2.design = "~ condition".into();
    config.corr.model = "~ condition".into();
    config.corr.threads = 1;
    config.save(&paths.config_path).unwrap();
    db.set_project_config(&config).unwrap();
    let mut sheet = "run_accession\tcondition\n".to_string();
    for (column, sample) in samples.iter().enumerate() {
        let condition = if column < 4 { "A" } else { "B" };
        sheet.push_str(&format!("{sample}\t{condition}\n"));
        db.upsert_run(&RunRecord {
            project_id: db.project_id().into(),
            study_accession: "synthetic".into(),
            source_study_accession: None,
            run_accession: sample.clone(),
            sample_accession: Some(format!("BIO{column}")),
            experiment_accession: None,
            library_layout: LibraryLayout::Single,
            instrument_platform: None,
            instrument_model: None,
            remote_files: vec![],
            metadata: BTreeMap::from([("condition".into(), condition.into())]),
            state: RunState::QuantDone,
            last_error: None,
            updated_at: now_rfc3339(),
        })
        .unwrap();
        let directory = paths.quant_run_dir(sample, "kallistors");
        fs::create_dir_all(&directory).unwrap();
        let counts = rows
            .iter()
            .map(|row| row[column + 1].parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        let library = counts.iter().sum::<f64>();
        let result = kallistors::quant::QuantResult {
            tpm: counts.iter().map(|value| value / library * 1e6).collect(),
            est_counts: counts,
            eff_lengths: vec![100.0; rows.len()],
            post_bias: None,
        };
        let target_names = names.iter().map(String::as_str).collect::<Vec<_>>();
        let lengths = vec![1000_u32; rows.len()];
        kallistors::quant::write_abundance_tsv_from_parts(
            &directory.join("abundance.tsv"),
            &target_names,
            &lengths,
            &result,
        )
        .unwrap();
        kallistors::quant::write_abundance_h5(
            &directory.join("abundance.h5"),
            kallistors::quant::H5QuantOutput {
                target_names: &target_names,
                target_lengths: &lengths,
                result: &result,
                bootstrap_est_counts: &[],
                reads_processed: library as u64,
                fragment_length_hist: None,
                bias_observed: None,
                index_version: 1,
                kallisto_version: "synthetic-fixture",
                call: "native analysis integration",
                start_time: "0",
            },
        )
        .unwrap();
        fs::write(
            directory.join("run_info.json"),
            serde_json::to_vec(&serde_json::json!({"n_targets":rows.len(), "n_processed":library}))
                .unwrap(),
        )
        .unwrap();
    }
    let reference = rnaa_stages::refs::EnsemblReferenceManager::new()
        .unwrap()
        .ensure_reference(&paths, &db, &config)
        .unwrap();
    for run in db.list_runs().unwrap() {
        rnaa_stages::quant::declare_external_quant(
            &run.run_accession,
            &paths.quant_run_dir(&run.run_accession, "kallistors"),
            &reference,
            "deterministic simulated abundance fixture",
        )
        .unwrap();
    }
    fs::write(paths.samplesheet_path(), sheet).unwrap();
    (paths, db)
}
