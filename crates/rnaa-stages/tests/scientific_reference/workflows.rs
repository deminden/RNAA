//! Public RNAA workflow versus independent R fits, including filtering/missingness.
use super::*;
use rnaa_core::metadata::{ColumnSpec, ColumnType};
use std::path::Path;

fn compare_tables(actual: &Path, expected: &Path, columns: &[&str], label: &str) -> Result<()> {
    let actual = read_table(actual)?;
    let expected = read_table(expected)?;
    assert_eq!(actual.len(), expected.len(), "{label}: row count");
    let actual = actual
        .iter()
        .map(|row| (&row["gene_id"], row))
        .collect::<std::collections::BTreeMap<_, _>>();
    for row in &expected {
        let gene = &row["gene_id"];
        let observed = actual.get(gene).context("missing oracle gene")?;
        for &column in columns {
            let expected = &row[column];
            let actual = &observed[column];
            let label = format!("{label} {gene} {column}");
            if expected == "NA" || expected == "NaN" {
                assert!(
                    actual.is_empty() || actual == "NA" || actual == "NaN",
                    "{label}: expected undefined, got {actual}"
                );
            } else {
                let tolerance = if columns.contains(&"log2FoldChange") {
                    1e-7
                } else {
                    1e-9
                };
                assert_close(
                    actual.parse()?,
                    expected.parse()?,
                    tolerance,
                    tolerance,
                    &label,
                );
            }
        }
    }
    Ok(())
}

fn run_case(case: &str, design: &str) -> Result<()> {
    let oracle = reference("de_workflows")?.join(case);
    let counts = read_table(&oracle.join("counts.tsv"))?;
    let metadata = read_table(&oracle.join("metadata.tsv"))?;
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let mapping = root.path().join("tx2gene.tsv");
    let mut map = String::from("transcript\tgene\n");
    for row in &counts {
        map.push_str(&format!("t{}\t{}\n", row["gene_id"], row["gene_id"]));
    }
    fs::write(&mapping, map)?;
    for sample in &metadata {
        let id = &sample["run_accession"];
        let total = counts
            .iter()
            .map(|row| row[id].parse::<f64>().unwrap())
            .sum::<f64>();
        let directory = paths.quant_run_dir(id, "kallistors");
        fs::create_dir_all(&directory)?;
        let mut abundance = String::from("target_id\tlength\teff_length\test_counts\ttpm\n");
        for row in &counts {
            let count = row[id].parse::<f64>()?;
            abundance.push_str(&format!(
                "t{}\t1000\t1000\t{count}\t{}\n",
                row["gene_id"],
                count / total * 1e6
            ));
        }
        fs::write(directory.join("abundance.tsv"), abundance)?;
    }
    let reference = ReferenceBundle {
        id: "oracle".into(),
        organism: "synthetic".into(),
        ensembl_release: "pinned".into(),
        cdna_path: mapping.clone(),
        gtf_path: mapping.clone(),
        kallisto_index_path: mapping.clone(),
        tx2gene_path: mapping.clone(),
        gene_annotation_path: None,
        manifest_path: mapping,
    };
    let mut config = ProjectConfig::default();
    config.deseq2.design = design.into();
    config
        .metadata
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    let contrasts = [("B_vs_A", "B", "A"), ("A_vs_B", "A", "B")]
        .into_iter()
        .map(|(name, a, b)| {
            ContrastSpec::new(Some(name.into()), "condition".into(), a.into(), b.into())
        })
        .collect::<Vec<_>>();
    RustDeseq2Runner.deseq2(
        "oracle",
        &oracle.join("metadata.tsv"),
        &reference,
        design,
        &contrasts,
        &paths,
        &config,
    )?;
    let directory = paths.de_project_dir("oracle");
    for direction in ["B_vs_A", "A_vs_B"] {
        compare_tables(
            &directory.join(format!("de_{direction}.tsv")),
            &oracle.join(format!("{direction}.tsv")),
            &[
                "baseMean",
                "log2FoldChange",
                "lfcSE",
                "stat",
                "pvalue",
                "padj",
            ],
            &format!("{case} {direction}"),
        )?;
    }
    let samples = metadata
        .iter()
        .map(|row| row["run_accession"].as_str())
        .collect::<Vec<_>>();
    compare_tables(
        &directory.join("gene_norm_counts.tsv"),
        &oracle.join("normalized_counts.tsv"),
        &samples,
        &format!("{case} normalized counts"),
    )?;
    compare_tables(
        &directory.join("vst.tsv"),
        &oracle.join("vst.tsv"),
        &samples,
        &format!("{case} VST"),
    )?;
    let actual = read_table(&directory.join("size_factors.tsv"))?;
    let expected = read_table(&oracle.join("size_factors.tsv"))?;
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(&expected) {
        assert_eq!(actual["sample"], expected["sample"]);
        assert_close(
            actual["size_factor"].parse()?,
            expected["size_factor"].parse()?,
            1e-12,
            1e-12,
            &format!("{case} size factor {}", actual["sample"]),
        );
    }
    Ok(())
}
#[test]
fn categorical_and_reverse_contrasts_match_r() -> Result<()> {
    run_case("categorical", "~ condition")
}
#[test]
fn numeric_covariate_matches_r() -> Result<()> {
    run_case("numeric", "~ age + condition")
}
#[test]
fn numeric_interaction_matches_r() -> Result<()> {
    run_case("numeric_interaction", "~ age * condition")
}
#[test]
fn factor_interaction_matches_r() -> Result<()> {
    run_case("factor_interaction", "~ batch * condition")
}
#[test]
fn three_levels_character_contrast_zero_handling_matches_r() -> Result<()> {
    run_case("three_levels", "~ condition")
}
