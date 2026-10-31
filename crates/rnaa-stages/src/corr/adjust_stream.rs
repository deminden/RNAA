//! Row-wise matrix adjustment with bounded input and output streams.
use super::limits::AdjustmentEstimate;
use super::residualize::{prepare_samples, residual_basis, residualize_row};
use super::stream::OutputBudget;
use anyhow::{Context, Result, bail};
use nalgebra::DMatrix;
use rnaa_core::{metadata::MetadataConfig, model::AdjustedMatrix, paths::ProjectPaths};
use std::{
    collections::BTreeSet,
    fs::{self, File},
    io::Read,
    path::Path,
};

/// Adjust a validated matrix one row at a time, respecting admitted byte counts.
#[allow(clippy::too_many_arguments)]
pub fn adjust_matrix_streaming(
    matrix_path: &Path,
    samplesheet: &Path,
    adjustment: &str,
    model: &str,
    paths: &ProjectPaths,
    schema: &MetadataConfig,
    estimate: AdjustmentEstimate,
    check: &impl Fn() -> Result<()>,
) -> Result<AdjustedMatrix> {
    check()?;
    let file = File::open(matrix_path)?;
    let before = file.metadata()?;
    if !before.is_file() || before.len() != estimate.input_bytes {
        bail!("matrix size changed after adjustment admission");
    }
    let cap = estimate
        .input_bytes
        .checked_add(1)
        .context("matrix size overflow")?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.take(cap));
    let samples = reader
        .headers()?
        .iter()
        .skip(1)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    if samples.len() != estimate.samples
        || samples.is_empty()
        || samples.iter().any(|id| id.trim().is_empty())
        || samples.iter().collect::<BTreeSet<_>>().len() != samples.len()
    {
        bail!("matrix contains duplicate or empty sample IDs");
    }
    let basis = match adjustment {
        "none" => None,
        "residualize" => {
            let prepared = prepare_samples(&samples, samplesheet, model, schema)?;
            let design = DMatrix::from_row_slice(
                samples.len(),
                prepared.contract.coefficient_names.len(),
                &prepared.contract.matrix,
            );
            rnaa_core::util::write_json_pretty(
                &paths.corr_dir.join("model_design.json"),
                &prepared.contract,
            )?;
            residual_basis(&design)?
        }
        other => bail!("unsupported adjustment '{other}'; expected none or residualize"),
    };
    check()?;
    let output = paths.corr_dir.join("adjusted_matrix.tsv");
    let file = File::create(&output)?;
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(OutputBudget::new(file, estimate.output_bytes));
    writer.write_record(std::iter::once("gene_id").chain(samples.iter().map(String::as_str)))?;
    let mut genes = Vec::new();
    let mut seen = BTreeSet::new();
    let mut record = csv::StringRecord::new();
    let mut gene_id_bytes = 0_u64;
    while reader.read_record(&mut record)? {
        if genes.len().is_multiple_of(256) {
            check()?;
        }
        if record.is_empty() {
            continue;
        }
        let gene = record.get(0).unwrap_or_default();
        gene_id_bytes = gene_id_bytes
            .checked_add(u64::try_from(gene.len())?)
            .context("matrix gene identifier size overflow")?;
        if genes.len() >= estimate.genes || gene_id_bytes > estimate.gene_id_bytes {
            bail!("matrix shape changed during adjustment");
        }
        if gene.trim().is_empty() || !seen.insert(gene.to_owned()) {
            bail!("matrix contains duplicate or empty gene IDs");
        }
        let values = record
            .iter()
            .skip(1)
            .map(|field| {
                field
                    .parse::<f64>()
                    .with_context(|| format!("failed to parse float '{field}'"))
            })
            .collect::<Result<Vec<_>>>()?;
        if values.len() != samples.len() || values.iter().any(|value| !value.is_finite()) {
            bail!("matrix values must be finite and match the sample dimensions");
        }
        let adjusted = residualize_row(&values, basis.as_ref());
        if adjusted.iter().any(|value| !value.is_finite()) {
            bail!("matrix values must be finite and match the sample dimensions");
        }
        let mut fields = Vec::with_capacity(adjusted.len() + 1);
        fields.push(gene.to_owned());
        fields.extend(adjusted.iter().map(|value| value.to_string()));
        writer.write_record(&fields)?;
        genes.push(gene.to_owned());
    }
    if genes.len() != estimate.genes || gene_id_bytes != estimate.gene_id_bytes {
        bail!("matrix shape changed during adjustment");
    }
    if reader.position().byte() != estimate.input_bytes
        || reader.get_ref().get_ref().metadata()?.len() != estimate.input_bytes
        || fs::metadata(matrix_path)?.len() != estimate.input_bytes
    {
        bail!("matrix size changed during adjustment");
    }
    check()?;
    writer.flush()?;
    Ok(AdjustedMatrix {
        matrix_path: output,
        genes,
        samples,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::corr::{Residualizer, estimate_adjustment};
    use rnaa_formats::matrix::{read_matrix_tsv, write_matrix_tsv};

    #[test]
    fn streamed_adjustment_matches_previous_outputs_byte_for_byte() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let input = temp.path().join("input.tsv");
        let sheet = temp.path().join("samples.tsv");
        fs::write(
            &input,
            "gene_id\tS1\tS2\tS3\tS4\n\"a\tquoted\"\t1\t2\t3\t4\nb\t0\t1\t9\t10\nc\t1e-300\t2\t4\t8\n\"d\nquoted\"\t1\t3\t5\t7\n",
        )?;
        fs::write(&sheet, "sample_id\tcondition\nS1\tA\nS2\tA\nS3\tB\nS4\tB\n")?;
        let previous = ProjectPaths::new(temp.path().join("previous"));
        let current = ProjectPaths::new(temp.path().join("current"));
        fs::create_dir_all(&previous.corr_dir)?;
        fs::create_dir_all(&current.corr_dir)?;
        let schema = MetadataConfig::default();
        for mode in ["none", "residualize"] {
            let old_path = previous.corr_dir.join("adjusted_matrix.tsv");
            if mode == "none" {
                write_matrix_tsv(&old_path, &read_matrix_tsv(&input)?)?;
            } else {
                Residualizer.adjust_typed(&input, &sheet, "~ condition", &previous, &schema)?;
            }
            let adjusted = adjust_matrix_streaming(
                &input,
                &sheet,
                mode,
                "~ condition",
                &current,
                &schema,
                estimate_adjustment(&input, &sheet, mode == "residualize")?,
                &|| Ok(()),
            )?;
            assert_eq!(
                fs::read(&adjusted.matrix_path)?,
                fs::read(&old_path)?,
                "{mode}"
            );
            assert_eq!(adjusted.genes, ["a\tquoted", "b", "c", "d\nquoted"]);
            assert_eq!(adjusted.samples, ["S1", "S2", "S3", "S4"]);
        }
        Ok(())
    }

    #[test]
    fn streamed_adjustment_bounds_output_checks_input_and_honors_cancellation() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let input = temp.path().join("input.tsv");
        let sheet = temp.path().join("samples.tsv");
        fs::write(&input, "gene_id\tS1\tS2\ng1\t1\t2\ng2\t3\t4\n")?;
        fs::write(&sheet, "sample_id\nS1\nS2\n")?;
        let paths = ProjectPaths::new(temp.path());
        fs::create_dir_all(&paths.corr_dir)?;
        let estimate = estimate_adjustment(&input, &sheet, false)?;
        let schema = MetadataConfig::default();
        let result = adjust_matrix_streaming(
            &input,
            &sheet,
            "none",
            "~ 1",
            &paths,
            &schema,
            AdjustmentEstimate {
                output_bytes: 1,
                ..estimate
            },
            &|| Ok(()),
        );
        assert!(result.is_err());
        assert!(fs::metadata(paths.corr_dir.join("adjusted_matrix.tsv"))?.len() <= 1);
        let checks = std::cell::Cell::new(0);
        let result = adjust_matrix_streaming(
            &input,
            &sheet,
            "none",
            "~ 1",
            &paths,
            &schema,
            estimate,
            &|| {
                checks.set(checks.get() + 1);
                if checks.get() == 3 {
                    bail!("cancelled in adjustment");
                }
                Ok(())
            },
        );
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("cancelled in adjustment")
        );
        assert!(
            adjust_matrix_streaming(
                &input,
                &sheet,
                "none",
                "~ 1",
                &paths,
                &schema,
                AdjustmentEstimate {
                    input_bytes: estimate.input_bytes - 1,
                    ..estimate
                },
                &|| Ok(())
            )
            .is_err()
        );
        fs::write(&input, "gene_id\tS1\tS2\ng1\t1\t2\ng1\t3\t4\n")?;
        assert!(
            adjust_matrix_streaming(
                &input,
                &sheet,
                "none",
                "~ 1",
                &paths,
                &schema,
                estimate_adjustment(&input, &sheet, false)?,
                &|| Ok(())
            )
            .is_err()
        );
        fs::write(&input, "gene_id\tS1\tS2\ng1\t100\t2\ng2\t3\t4\n")?;
        let earlier = estimate_adjustment(&input, &sheet, false)?;
        fs::write(&input, "gene_id\tS1\tS2\ng123\t1\t2\ng2\t3\t4\n")?;
        assert_eq!(fs::metadata(&input)?.len(), earlier.input_bytes);
        assert!(
            adjust_matrix_streaming(
                &input,
                &sheet,
                "none",
                "~ 1",
                &paths,
                &schema,
                earlier,
                &|| Ok(())
            )
            .unwrap_err()
            .to_string()
            .contains("shape changed")
        );
        Ok(())
    }

    #[test]
    fn adjustment_checks_cancellation_after_processing_rows() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let input = temp.path().join("input.tsv");
        let sheet = temp.path().join("samples.tsv");
        let mut data = String::from("gene_id\tS1\tS2\n");
        for index in 0..600 {
            data.push_str(&format!("gene{index}\t{index}\t{}\n", index + 1));
        }
        fs::write(&input, data)?;
        fs::write(&sheet, "sample_id\nS1\nS2\n")?;
        let paths = ProjectPaths::new(temp.path());
        fs::create_dir_all(&paths.corr_dir)?;
        let estimate = estimate_adjustment(&input, &sheet, false)?;
        let checks = std::cell::Cell::new(0);
        let result = adjust_matrix_streaming(
            &input,
            &sheet,
            "none",
            "~ 1",
            &paths,
            &MetadataConfig::default(),
            estimate,
            &|| {
                checks.set(checks.get() + 1);
                if checks.get() == 4 {
                    bail!("cancelled after rows");
                }
                Ok(())
            },
        );
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("cancelled after rows")
        );
        assert_eq!(checks.get(), 4);
        Ok(())
    }
}
