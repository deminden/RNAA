//! Independent R graph/PCA/fgsea reference values through the native module API.
use super::*;
use rnaa_core::model::AdjustedMatrix;
use rnaa_stages::corr::{
    CorrelationLimits, ModuleGseaParams, run_module_gsea_with_limits, standardized_module_score,
};
use std::{collections::BTreeMap, path::Path};

#[test]
fn standardized_pca_matches_r_including_scales_zero_average_and_ties() -> Result<()> {
    let oracle = reference("modules")?;
    for case in [
        "unequal_scales",
        "zero_average",
        "tied_leading",
        "three_samples",
    ] {
        let directory = oracle.join(case);
        let matrix = rnaa_formats::matrix::read_matrix_tsv(&directory.join("input.tsv"))?;
        let result = standardized_module_score(&matrix.values)?;
        if case == "unequal_scales" {
            let mut mean = (0..matrix.col_ids.len())
                .map(|sample| {
                    matrix.values.iter().map(|row| row[sample]).sum::<f64>()
                        / matrix.values.len() as f64
                })
                .collect::<Vec<_>>();
            let center = mean.iter().sum::<f64>() / mean.len() as f64;
            for value in &mut mean {
                *value -= center;
            }
            let sd = (mean.iter().map(|v| v * v).sum::<f64>() / (mean.len() - 1) as f64).sqrt();
            assert!(
                mean.iter()
                    .zip(&result.values)
                    .any(|(a, b)| (a / sd - b).abs() > 0.1),
                "reference must distinguish PCA from the former arithmetic mean"
            );
        }
        let expected = read_table(&directory.join("scores.tsv"))?;
        assert_eq!(result.values.len(), expected.len());
        for (sample, (value, expected)) in result.values.iter().zip(&expected).enumerate() {
            assert_eq!(matrix.col_ids[sample], expected["sample"]);
            assert_close(
                *value,
                expected["score"].parse()?,
                1e-9,
                1e-9,
                &format!("R PCA {case} {}", expected["sample"]),
            );
        }
        let diagnostics = read_table(&directory.join("diagnostics.tsv"))?;
        assert_close(
            result.variance_explained,
            diagnostics[0]["variance_explained"].parse()?,
            1e-9,
            1e-9,
            &format!("R PCA {case} variance explained"),
        );
        assert_eq!(
            result.leading_multiplicity,
            diagnostics[0]["leading_multiplicity"].parse::<usize>()?
        );
        let wgcna = directory.join("wgcna.tsv");
        if wgcna.exists() {
            for (value, expected) in result.values.iter().zip(read_table(&wgcna)?) {
                assert_close(
                    *value,
                    expected["normalized_wgcna"].parse()?,
                    1e-9,
                    1e-9,
                    &format!("WGCNA {case} score"),
                );
            }
        }
        let mut reversed = matrix.values;
        reversed.reverse();
        for (left, right) in result
            .values
            .iter()
            .zip(standardized_module_score(&reversed)?.values)
        {
            assert_close(
                *left,
                right,
                1e-9,
                1e-9,
                &format!("{case} membership-order invariance"),
            );
        }
    }
    Ok(())
}

fn keyed(path: &Path, keys: &[&str]) -> Result<BTreeMap<Vec<String>, BTreeMap<String, String>>> {
    let mut result = BTreeMap::new();
    for row in read_table(path)? {
        let key = keys.iter().map(|key| row[*key].clone()).collect::<Vec<_>>();
        anyhow::ensure!(
            !result.contains_key(&key),
            "duplicate reference/output identity {key:?} in {}",
            path.display()
        );
        result.insert(key, row);
    }
    Ok(result)
}
fn compare(actual: &Path, expected: &Path, keys: &[&str], values: &[&str]) -> Result<()> {
    let actual = keyed(actual, keys)?;
    let expected = keyed(expected, keys)?;
    assert_eq!(
        actual.keys().collect::<Vec<_>>(),
        expected.keys().collect::<Vec<_>>(),
        "{} reference identities",
        actual.len()
    );
    for (key, row) in expected {
        for &column in values {
            assert_close(
                actual[&key][column].parse()?,
                row[column].parse()?,
                1e-9,
                1e-9,
                &format!("R module {key:?} {column}"),
            );
        }
    }
    Ok(())
}
#[test]
fn full_modules_scores_rankings_and_enrichment_match_r_across_allocations_and_report_limits()
-> Result<()> {
    let oracle = reference("modules")?.join("workflow");
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let matrix = AdjustedMatrix {
        matrix_path: oracle.join("input.tsv"),
        genes: vec![],
        samples: vec![],
    };
    let mut original_scores = None;
    for (block_rows, threads, top_genes) in [(1, 1, 1), (7, 2, 5), (64, 3, 500)] {
        let params = ModuleGseaParams {
            gmt_path: oracle.join("pathways.gmt"),
            top_genes,
            min_module_size: 2,
            permutations: 100,
            seed: 42,
            pathway_min_size: 10,
            pathway_max_size: None,
        };
        run_module_gsea_with_limits(
            &matrix,
            &paths,
            "oracle",
            &params,
            CorrelationLimits {
                block_rows,
                ..Default::default()
            },
            threads,
        )?;
        let directory = paths.corr_project_dir("oracle");
        compare(
            &directory.join("module_membership.tsv"),
            &oracle.join("membership.tsv"),
            &["module_id", "gene"],
            &["module_size", "connectivity", "rank_in_module"],
        )?;
        compare(
            &directory.join("module_eigengenes.tsv"),
            &oracle.join("scores.tsv"),
            &["module_id"],
            &(1..=8)
                .map(|i| format!("sample{i}"))
                .collect::<Vec<_>>()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
        )?;
        compare(
            &directory.join("module_t_values.tsv"),
            &oracle.join("rankings.tsv"),
            &["module_id", "gene"],
            &["spearman_r", "t_value"],
        )?;
        compare(
            &directory.join("module_fgsea.tsv"),
            &oracle.join("enrichment.tsv"),
            &["module_id", "pathway"],
            &[
                "es",
                "nes",
                "p_value",
                "padj",
                "pathway_size",
                "leading_edge_size",
            ],
        )?;
        let scores = fs::read(directory.join("module_eigengenes.tsv"))?;
        if let Some(expected) = &original_scores {
            assert_eq!(
                &scores, expected,
                "report limits and thread/tile allocation cannot alter PCA scores"
            );
        } else {
            original_scores = Some(scores);
        }
        let top = read_table(&directory.join("modules_top_genes.tsv"))?;
        assert_eq!(top.len(), 2 * top_genes.min(15));
        assert!(top.iter().all(|row| row["module_size"] == "15"));
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(directory.join("module_fgsea_manifest.json"))?)?;
        assert_eq!(manifest["parameters"]["admission"]["eligible_genes"], 30);
        assert_eq!(manifest["parameters"]["admission"]["filtered_pathways"], 1);
        assert_eq!(
            manifest["parameters"]["admission"]["unmatched_gene_entries"],
            2
        );
        assert!(
            manifest["parameters"]["graph_statistics"]["peak_tile_entries"]
                .as_u64()
                .unwrap()
                <= (block_rows * block_rows) as u64
        );
        assert!(
            manifest["parameters"]["graph_statistics"]["retained_neighbors"]
                .as_u64()
                .unwrap()
                <= 31 * 20
        );
    }
    Ok(())
}
