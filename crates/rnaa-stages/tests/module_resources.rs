//! Module admission rejects ambiguous inputs and unbounded requests before output.
use anyhow::Result;
use rnaa_core::{ProjectPaths, model::AdjustedMatrix};
use rnaa_stages::corr::{
    CorrelationLimits, ModuleGseaParams, run_module_gsea_with_limits, standardized_module_score,
};
use std::{fs, path::Path};

fn setup(root: &Path, genes: usize) -> Result<(ProjectPaths, AdjustedMatrix, ModuleGseaParams)> {
    let paths = ProjectPaths::new(root);
    paths.ensure_layout()?;
    let input = root.join("input.tsv");
    let gmt = root.join("pathways.gmt");
    let mut data = String::from("gene_id\tsample1\tsample2\tsample3\tsample4\n");
    for i in 0..genes {
        data.push_str(&format!(
            "G{i}\t{}\t{}\t{}\t{}\n",
            i % 7,
            (i * 3 + 1) % 11,
            (i * 5 + 7) % 13,
            (i * 11) % 17
        ));
    }
    fs::write(&input, data)?;
    fs::write(&gmt, "PATH\tsynthetic\tG0\tG1\tG2\n")?;
    Ok((
        paths,
        AdjustedMatrix {
            matrix_path: input,
            genes: vec![],
            samples: vec![],
        },
        ModuleGseaParams {
            gmt_path: gmt,
            top_genes: 1,
            min_module_size: 2,
            permutations: 100,
            seed: 42,
            pathway_min_size: 2,
            pathway_max_size: None,
        },
    ))
}
#[test]
fn malformed_or_nonoverlapping_pathways_fail_before_outputs() -> Result<()> {
    for (gmt, reason) in [
        ("bad\tline\n", "invalid"),
        ("\tx\tG0\n", "invalid"),
        ("PATH\tx\tG0\tG1\nPATH\tx\tG1\tG2\n", "duplicate"),
        ("PATH\tx\tG0\tG0\n", "duplicate"),
        ("PATH\tx\tG0\t\n", "empty"),
        ("\n", "no pathways"),
        ("PATH\tx\tUNKNOWN\n", "no GMT pathways"),
    ] {
        let root = tempfile::tempdir()?;
        let (paths, matrix, params) = setup(root.path(), 30)?;
        fs::write(&params.gmt_path, gmt)?;
        let error = run_module_gsea_with_limits(
            &matrix,
            &paths,
            "module",
            &params,
            CorrelationLimits::default(),
            1,
        )
        .unwrap_err();
        assert!(error.to_string().contains(reason), "{error:#}");
        assert!(!paths.corr_project_dir("module").exists());
    }
    Ok(())
}
#[test]
fn options_are_rejected_instead_of_silently_clamped() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (paths, matrix, base) = setup(root.path(), 30)?;
    let mut invalid = Vec::new();
    let mut p = base.clone();
    p.top_genes = 0;
    invalid.push(p);
    let mut p = base.clone();
    p.min_module_size = 1;
    invalid.push(p);
    let mut p = base.clone();
    p.min_module_size = 31;
    invalid.push(p);
    let mut p = base.clone();
    p.permutations = 99;
    invalid.push(p);
    let mut p = base.clone();
    p.pathway_min_size = 0;
    invalid.push(p);
    let mut p = base.clone();
    p.pathway_max_size = Some(1);
    invalid.push(p);
    let mut p = base.clone();
    p.pathway_max_size = Some(30);
    invalid.push(p);
    for params in invalid {
        assert!(
            run_module_gsea_with_limits(
                &matrix,
                &paths,
                "module",
                &params,
                CorrelationLimits::default(),
                1
            )
            .is_err()
        );
        assert!(!paths.corr_project_dir("module").exists());
    }
    Ok(())
}
#[test]
fn graph_pca_permutation_and_output_estimates_are_checked_before_work() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (paths, matrix, mut params) = setup(root.path(), 1000)?;
    for (limits, reason) in [
        (
            CorrelationLimits {
                max_memory_mb: 1,
                ..Default::default()
            },
            "working estimate",
        ),
        (
            CorrelationLimits {
                max_output_mb: 1,
                ..Default::default()
            },
            "output estimate",
        ),
        (
            CorrelationLimits {
                block_rows: 0,
                ..Default::default()
            },
            "invalid module resource",
        ),
    ] {
        let error =
            run_module_gsea_with_limits(&matrix, &paths, "module", &params, limits, 1).unwrap_err();
        assert!(error.to_string().contains(reason), "{error:#}");
        assert!(!paths.corr_project_dir("module").exists());
    }
    params.permutations = usize::MAX;
    assert!(
        run_module_gsea_with_limits(
            &matrix,
            &paths,
            "module",
            &params,
            CorrelationLimits::default(),
            2
        )
        .unwrap_err()
        .to_string()
        .contains("working estimate")
    );
    assert!(!paths.corr_project_dir("module").exists());
    Ok(())
}
#[test]
fn bounded_graph_keeps_only_neighbors_and_reports_full_membership() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (paths, matrix, params) = setup(root.path(), 257)?;
    let outputs = run_module_gsea_with_limits(
        &matrix,
        &paths,
        "module",
        &params,
        CorrelationLimits {
            block_rows: 7,
            ..Default::default()
        },
        2,
    )?;
    let manifest: serde_json::Value = serde_json::from_slice(&fs::read(outputs.last().unwrap())?)?;
    let stats = &manifest["parameters"]["graph_statistics"];
    assert!(stats["peak_tile_entries"].as_u64().unwrap() <= 49);
    assert!(stats["retained_neighbors"].as_u64().unwrap() <= 257 * 20);
    assert_eq!(manifest["parameters"]["threads"], 2);
    assert!(
        fs::read_to_string(
            paths
                .corr_project_dir("module")
                .join("module_membership.tsv")
        )?
        .lines()
        .count()
            > fs::read_to_string(
                paths
                    .corr_project_dir("module")
                    .join("modules_top_genes.tsv")
            )?
            .lines()
            .count()
    );
    Ok(())
}
#[test]
fn pca_rejects_constants_nonfinite_and_wrong_dimensions() {
    for rows in [
        vec![],
        vec![vec![1.0, 2.0, 3.0]],
        vec![vec![1.0, 2.0], vec![2.0, 3.0]],
        vec![vec![1.0, 2.0, 3.0], vec![1.0, 1.0, 1.0]],
        vec![vec![1.0, 2.0, 3.0], vec![1.0, f64::NAN, 2.0]],
        vec![vec![1.0, 2.0, 3.0], vec![1.0, 2.0]],
    ] {
        assert!(standardized_module_score(&rows).is_err());
    }
}
