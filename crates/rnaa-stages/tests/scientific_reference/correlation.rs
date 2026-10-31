use super::*;
use rnaa_core::{
    model::{AdjustedMatrix, CorrelationMethod, OutputMode},
    traits::Correlator,
};
use rnaa_stages::corr::{CorrelationLimits, MinCorrCorrelator};
fn compare(method: CorrelationMethod) -> Result<()> {
    let oracle = reference("correlation")?;
    for case in [
        "ties_and_outliers",
        "minimum_two_samples",
        "three_samples",
        "near_perfect",
    ] {
        let input = oracle.join(case).join("input.tsv");
        let data = read_table(&input)?;
        let genes = data
            .iter()
            .map(|row| row["gene_id"].clone())
            .collect::<Vec<_>>();
        let samples = (1..data[0].len()).map(|i| format!("sample{i}")).collect();
        let matrix = AdjustedMatrix {
            matrix_path: input,
            genes: genes.clone(),
            samples,
        };
        let root = tempfile::tempdir()?;
        let paths = ProjectPaths::new(root.path());
        paths.ensure_layout()?;
        MinCorrCorrelator.correlate(&matrix, method, &OutputMode::Dense, 1, &paths, "oracle")?;
        let actual = read_table(&paths.corr_project_dir("oracle").join("corr_dense.tsv"))?;
        let expected = read_table(&oracle.join(case).join(format!("{method}.tsv")))?;
        assert_eq!(actual.len(), expected.len());
        for (row, expected) in actual.iter().zip(expected) {
            assert_eq!(row["gene_id"], expected["gene_id"]);
            for gene in &genes {
                let label = format!("{method} {case} {}->{gene}", row["gene_id"]);
                let actual = row[gene].parse::<f64>()?;
                if expected[gene] == "NA" || expected[gene] == "NaN" {
                    assert!(actual.is_nan(), "{label}: expected undefined, got {actual}");
                } else {
                    assert_close(actual, expected[gene].parse()?, 1e-9, 1e-9, &label);
                }
            }
        }
        for block_rows in [1, 2, 3] {
            let correlator = MinCorrCorrelator.with_limits(CorrelationLimits {
                block_rows,
                ..Default::default()
            });
            for (mode, name) in [
                (OutputMode::TopK { k: 2 }, "topk2"),
                (OutputMode::Threshold { min_abs_r: 0.0 }, "threshold0"),
                (OutputMode::Threshold { min_abs_r: 0.75 }, "threshold0.75"),
                (OutputMode::Threshold { min_abs_r: 1.0 }, "threshold1"),
            ] {
                let outputs = correlator.correlate(&matrix, method, &mode, 2, &paths, "oracle")?;
                let actual = read_table(&outputs[0])?;
                let expected = read_table(&oracle.join(case).join(format!("{method}_{name}.tsv")))?;
                assert_eq!(
                    actual.len(),
                    expected.len(),
                    "{method} {case} {name} block={block_rows}: edge count"
                );
                for (actual, expected) in actual.iter().zip(&expected) {
                    let label = format!("{method} {case} {name} block={block_rows}");
                    assert_eq!(
                        actual["gene_a"], expected["gene_a"],
                        "{label}: source order"
                    );
                    assert_eq!(
                        actual["gene_b"], expected["gene_b"],
                        "{label}: target order"
                    );
                    assert_close(
                        actual["r"].parse()?,
                        expected["r"].parse()?,
                        1e-9,
                        1e-9,
                        &label,
                    );
                }
            }
        }
    }
    Ok(())
}
#[test]
fn pearson_matches_pinned_r() -> Result<()> {
    compare(CorrelationMethod::Pearson)
}
#[test]
fn spearman_matches_pinned_r() -> Result<()> {
    compare(CorrelationMethod::Spearman)
}
#[test]
fn kendall_matches_pinned_r() -> Result<()> {
    compare(CorrelationMethod::Kendall)
}
#[test]
fn bicor_matches_pinned_wgcna() -> Result<()> {
    compare(CorrelationMethod::Bicor)
}
#[test]
fn hellcor_matches_pinned_hellcor() -> Result<()> {
    compare(CorrelationMethod::Hellcor)
}
