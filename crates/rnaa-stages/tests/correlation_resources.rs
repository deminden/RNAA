//! Resource rejection and bounded selection through the public correlation API.
use anyhow::Result;
use rnaa_core::{
    ProjectPaths,
    model::{AdjustedMatrix, CorrelationMethod, OutputMode},
    traits::Correlator,
};
use rnaa_stages::corr::{CorrelationLimits, MinCorrCorrelator};
use std::fs;

fn input(root: &std::path::Path, genes: usize) -> Result<AdjustedMatrix> {
    let path = root.join("matrix.tsv");
    let mut data = String::from("gene_id\tsample1\tsample2\tsample3\tsample4\n");
    for i in 0..genes {
        data.push_str(&format!(
            "gene{i:05}\t{}\t{}\t{}\t{}\n",
            i % 11,
            (i * 3 + 7) % 17,
            (i + 5) % 19,
            (i * 7) % 23
        ));
    }
    fs::write(&path, data)?;
    Ok(AdjustedMatrix {
        matrix_path: path,
        genes: vec![],
        samples: vec![],
    })
}
#[test]
fn bounded_topk_matches_across_tiles_and_threads() -> Result<()> {
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let matrix = input(root.path(), 257)?;
    let mut expected = None;
    for (block_rows, threads) in [(1, 1), (7, 2), (64, 3)] {
        let outputs = MinCorrCorrelator
            .with_limits(CorrelationLimits {
                block_rows,
                ..Default::default()
            })
            .correlate(
                &matrix,
                CorrelationMethod::Spearman,
                &OutputMode::TopK { k: 5 },
                threads,
                &paths,
                "resource",
            )?;
        let bytes = fs::read(&outputs[0])?;
        if let Some(expected) = &expected {
            assert_eq!(
                &bytes, expected,
                "output must preserve exact ties across resource changes"
            );
        } else {
            expected = Some(bytes);
        }
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(outputs.last().unwrap())?)?;
        let stats = &manifest["parameters"]["statistics"];
        assert!(stats["peak_block_entries"].as_u64().unwrap() <= (block_rows * block_rows) as u64);
        assert!(stats["peak_retained_edges"].as_u64().unwrap() <= (block_rows * 5) as u64);
        assert_eq!(stats["peak_stripe_entries"], 0);
        assert_eq!(
            fs::read_to_string(&outputs[0])?.lines().count(),
            1 + 257 * 5
        );
    }
    Ok(())
}
#[test]
fn oversized_output_or_working_set_is_rejected_before_output() -> Result<()> {
    for (limits, reason) in [
        (
            CorrelationLimits {
                max_output_mb: 1,
                ..Default::default()
            },
            "output estimate",
        ),
        (
            CorrelationLimits {
                max_memory_mb: 1,
                ..Default::default()
            },
            "working estimate",
        ),
    ] {
        let root = tempfile::tempdir()?;
        let paths = ProjectPaths::new(root.path());
        paths.ensure_layout()?;
        let matrix = input(root.path(), 1000)?;
        let error = MinCorrCorrelator
            .with_limits(limits)
            .correlate(
                &matrix,
                CorrelationMethod::Pearson,
                &OutputMode::Dense,
                1,
                &paths,
                "resource",
            )
            .unwrap_err();
        assert!(error.to_string().contains(reason), "{error:#}");
        assert!(!paths.corr_project_dir("resource").exists());
    }
    Ok(())
}
#[test]
fn invalid_modes_limits_and_single_sample_fail_before_output() -> Result<()> {
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let matrix = input(root.path(), 3)?;
    for mode in [
        OutputMode::TopK { k: 0 },
        OutputMode::Threshold {
            min_abs_r: f64::NAN,
        },
        OutputMode::Threshold { min_abs_r: 1.1 },
    ] {
        assert!(
            MinCorrCorrelator
                .correlate(&matrix, CorrelationMethod::Pearson, &mode, 1, &paths, "bad")
                .is_err()
        );
    }
    assert!(
        MinCorrCorrelator
            .with_limits(CorrelationLimits {
                block_rows: 0,
                ..Default::default()
            })
            .correlate(
                &matrix,
                CorrelationMethod::Pearson,
                &OutputMode::Dense,
                1,
                &paths,
                "bad"
            )
            .is_err()
    );
    fs::write(&matrix.matrix_path, "gene_id\tone\na\t1\nb\t2\n")?;
    assert!(
        MinCorrCorrelator
            .correlate(
                &matrix,
                CorrelationMethod::Pearson,
                &OutputMode::Dense,
                1,
                &paths,
                "bad"
            )
            .unwrap_err()
            .to_string()
            .contains("two samples")
    );
    assert!(!paths.corr_project_dir("bad").exists());
    Ok(())
}
