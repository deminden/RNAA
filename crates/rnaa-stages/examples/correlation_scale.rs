//! Reproducible process-level scale probe; measure RSS externally with time -l/-v.
use anyhow::{Context, Result};
use rnaa_core::{
    ProjectPaths,
    model::{CorrelationMethod, OutputMode},
    traits::Correlator,
};
use rnaa_stages::corr::{CorrelationLimits, MinCorrCorrelator};
use std::{fs, time::Instant};
mod support;

fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let number = |i: usize, default: usize| -> Result<usize> {
        args.get(i)
            .map(|s| s.parse().context("expected positive integer"))
            .unwrap_or(Ok(default))
    };
    let (genes, samples, k, block, threads) = (
        number(0, 4000)?,
        number(1, 16)?,
        number(2, 10)?,
        number(3, 64)?,
        number(4, 2)?,
    );
    let method = args
        .get(5)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(CorrelationMethod::Pearson);
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let matrix = support::input(root.path(), genes, samples)?;
    let path = matrix.matrix_path.clone();
    let start = Instant::now();
    let outputs = MinCorrCorrelator
        .with_limits(CorrelationLimits {
            block_rows: block,
            ..Default::default()
        })
        .correlate(
            &matrix,
            method,
            &OutputMode::TopK { k },
            threads,
            &paths,
            "scale",
        )?;
    let elapsed = start.elapsed();
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(outputs.last().context("missing manifest")?)?)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({
            "genes":genes,"samples":samples,"k":k,"block_rows":block,"threads":threads,"method":method.to_string(),
            "correlation_wall_seconds":elapsed.as_secs_f64(),"input_bytes":fs::metadata(path)?.len(),
            "output_bytes":fs::metadata(&outputs[0])?.len(),"statistics":manifest["parameters"]["statistics"],
            "resource_estimate":manifest["parameters"]["resource_estimate"],
            "profile":if cfg!(debug_assertions) { "debug" } else { "release" }
        }))?
    );
    Ok(())
}
