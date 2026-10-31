//! Reproducible bounded module graph, PCA and enrichment process-level scale probe.
mod support;
use anyhow::{Context, Result, ensure};
use rnaa_core::ProjectPaths;
use rnaa_stages::corr::{CorrelationLimits, ModuleGseaParams, run_module_gsea_with_limits};
use std::{fs, time::Instant};

fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let number = |i: usize, default| -> Result<usize> {
        args.get(i)
            .map(|s| s.parse().context("expected positive integer"))
            .unwrap_or(Ok(default))
    };
    let (genes, samples, top, block, threads) = (
        number(0, 2000)?,
        number(1, 16)?,
        number(2, 10)?,
        number(3, 64)?,
        number(4, 2)?,
    );
    ensure!(
        args.get(5).is_none_or(|method| method == "spearman"),
        "module scale probe uses Spearman"
    );
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let matrix = support::input(root.path(), genes, samples)?;
    let gmt = root.path().join("pathways.gmt");
    fs::write(
        &gmt,
        format!(
            "PROBE\tsynthetic\t{}\n",
            (0..100.min(genes / 2))
                .map(|i| format!("gene{i:07}"))
                .collect::<Vec<_>>()
                .join("\t")
        ),
    )?;
    let params = ModuleGseaParams {
        gmt_path: gmt,
        top_genes: top,
        min_module_size: 2,
        permutations: 100,
        seed: 42,
        pathway_min_size: 10,
        pathway_max_size: None,
    };
    let start = Instant::now();
    let outputs = run_module_gsea_with_limits(
        &matrix,
        &paths,
        "scale",
        &params,
        CorrelationLimits {
            block_rows: block,
            ..Default::default()
        },
        threads,
    )?;
    let elapsed = start.elapsed();
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(outputs.last().context("missing manifest")?)?)?;
    let output_bytes = outputs
        .iter()
        .map(|path| fs::metadata(path).map(|m| m.len()))
        .collect::<std::io::Result<Vec<_>>>()?
        .iter()
        .sum::<u64>();
    println!(
        "{}",
        serde_json::to_string_pretty(
            &serde_json::json!({"workflow":"module_enrichment","genes":genes,"samples":samples,"top_reported":top,"block_rows":block,"threads":threads,"method":"spearman","permutations":100,"seed":42,
        "module_wall_seconds":elapsed.as_secs_f64(),"input_bytes":fs::metadata(&matrix.matrix_path)?.len(),"output_bytes":output_bytes,
        "statistics":manifest["parameters"]["graph_statistics"],"modules":manifest["parameters"]["modules"],"resource_estimate":manifest["parameters"]["admission"],
        "profile":if cfg!(debug_assertions) { "debug" } else { "release" }})
        )?
    );
    Ok(())
}
