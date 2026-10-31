//! Mandatory offline numerical gates against pinned independent scientific output.
#[path = "scientific_reference/correlation.rs"]
mod correlation;
#[path = "scientific_reference/import.rs"]
mod import;
#[path = "scientific_reference/modules.rs"]
mod modules;
mod support;
#[path = "scientific_reference/workflows.rs"]
mod workflows;
use anyhow::{Context, Result};
use rnaa_core::{
    ContrastSpec, ProjectConfig, ProjectPaths, ReferenceBundle, traits::DifferentialExpression,
};
use rnaa_stages::deseq2::RustDeseq2Runner;
use std::fs;
use support::{assert_close, read_table, reference};

#[test]
fn default_de_workflow_matches_r_deseq2_including_parametric_failure_local_fallback() -> Result<()>
{
    let oracle = reference("deseq2")?;
    let expected = read_table(&oracle.join("results_wald_ratio.tsv"))?;
    let counts = read_table(&oracle.join("counts.tsv"))?;
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let samples = (1..=4).map(|i| format!("sample{i}")).collect::<Vec<_>>();
    let map = root.path().join("tx2gene.tsv");
    let mut mapping = String::from("transcript\tgene\n");
    for row in &counts {
        mapping.push_str(&format!("t{}\t{}\n", row["gene"], row["gene"]));
    }
    fs::write(&map, mapping)?;
    let sheet = root.path().join("samples.tsv");
    let mut metadata = String::from("run_accession\tcondition\n");
    for (i, sample) in samples.iter().enumerate() {
        metadata.push_str(&format!("{sample}\t{}\n", if i < 2 { "A" } else { "B" }));
        let directory = paths.quant_run_dir(sample, "kallistors");
        fs::create_dir_all(&directory)?;
        let total = counts
            .iter()
            .map(|row| row[sample].parse::<f64>())
            .collect::<std::result::Result<Vec<_>, _>>()?
            .iter()
            .sum::<f64>();
        let mut abundance = String::from("target_id\tlength\teff_length\test_counts\ttpm\n");
        for row in &counts {
            let count = row[sample].parse::<f64>()?;
            abundance.push_str(&format!(
                "t{}\t1000\t1000\t{count}\t{}\n",
                row["gene"],
                count / total * 1e6
            ));
        }
        fs::write(directory.join("abundance.tsv"), abundance)?;
    }
    fs::write(&sheet, metadata)?;
    let reference = ReferenceBundle {
        id: "oracle".into(),
        organism: "synthetic".into(),
        ensembl_release: "pinned".into(),
        cdna_path: map.clone(),
        gtf_path: map.clone(),
        kallisto_index_path: map.clone(),
        tx2gene_path: map.clone(),
        gene_annotation_path: None,
        manifest_path: map,
    };
    let mut config = ProjectConfig::default();
    config.deseq2.design = "~ condition".into();
    let contrast = ContrastSpec::new(
        Some("B_vs_A".into()),
        "condition".into(),
        "B".into(),
        "A".into(),
    );
    RustDeseq2Runner.deseq2(
        "oracle",
        &sheet,
        &reference,
        "~ condition",
        &[contrast],
        &paths,
        &config,
    )?;
    let directory = paths.de_project_dir("oracle");
    let actual = read_table(&directory.join("de_B_vs_A.tsv"))?;
    assert_eq!(actual.len(), expected.len());
    for (row, expected) in actual.iter().zip(&expected) {
        assert_eq!(row["gene_id"], expected["gene"]);
        for column in [
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
        ] {
            assert_close(
                row[column].parse()?,
                expected[column].parse()?,
                1e-8,
                1e-8,
                &format!("DESeq2 1.52.0 {} {column}", row["gene_id"]),
            );
        }
    }
    let actual = read_table(&directory.join("gene_norm_counts.tsv"))?;
    let expected = read_table(&oracle.join("normalized_counts_ratio.tsv"))?;
    for (row, expected) in actual.iter().zip(expected) {
        assert_eq!(row["gene_id"], expected["gene"]);
        for sample in &samples {
            assert_close(
                row[sample].parse()?,
                expected[sample].parse()?,
                1e-10,
                1e-11,
                &format!("normalized counts {} {sample}", row["gene_id"]),
            );
        }
    }
    let actual = read_table(&directory.join("size_factors.tsv"))?;
    let expected = read_table(&oracle.join("size_factors_ratio.tsv"))?;
    for (row, expected) in actual.iter().zip(expected) {
        assert_eq!(row["sample"], expected["sample"]);
        assert_close(
            row["size_factor"].parse()?,
            expected["size_factor"].parse()?,
            1e-12,
            1e-12,
            &format!("size factor {}", row["sample"]),
        );
    }
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(directory.join("deseq2_manifest.json"))?)?;
    assert_eq!(
        manifest["parameters"]["dispersion_fit"]["selected"],
        "local"
    );
    assert!(
        manifest["parameters"]["dispersion_fit"]["fallback_reason"]
            .as_str()
            .context("missing fallback evidence")?
            .contains("parametric")
    );
    Ok(())
}
