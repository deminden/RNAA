use super::*;
use rnaa_stages::tximport::{complete_length_offsets, import_kallisto_abundances};

#[test]
fn split_technical_runs_preserve_r_import_counts_and_lengths() -> Result<()> {
    use rnaa_stages::tximport::{
        QuantifiedSample, RunAbundance, import_grouped_kallisto_abundances,
    };
    let oracle = reference("tximport")?;
    for case in ["fractional_and_zero_abundance", "zero_length_without_mass"] {
        let root = tempfile::tempdir()?;
        let mut samples = Vec::new();
        for sample in 1..=4 {
            let rows = read_table(&oracle.join(case).join(format!("sample{sample}.tsv")))?;
            let mut runs = Vec::new();
            for lane in 1..=2 {
                let run_id = format!("sample{sample}_lane{lane}");
                let path = root.path().join(format!("{run_id}.tsv"));
                let mut writer = csv::WriterBuilder::new()
                    .delimiter(b'\t')
                    .from_path(&path)?;
                writer.write_record(["target_id", "length", "eff_length", "est_counts", "tpm"])?;
                for row in &rows {
                    writer.write_record([
                        row["target_id"].clone(),
                        row["length"].clone(),
                        row["eff_length"].clone(),
                        (row["est_counts"].parse::<f64>()? / 2.0).to_string(),
                        row["tpm"].clone(),
                    ])?;
                }
                writer.flush()?;
                runs.push(RunAbundance { run_id, path });
            }
            samples.push(QuantifiedSample {
                sample_id: format!("sample{sample}"),
                runs,
            });
        }
        for mode in ["no", "scaledTPM", "lengthScaledTPM"] {
            let actual = import_grouped_kallisto_abundances(
                &samples,
                &oracle.join(case).join("tx2gene.tsv"),
                mode,
            )?;
            let expected = read_table(&oracle.join(case).join(format!("{mode}_counts.tsv")))?;
            assert_eq!(actual.gene_ids.len(), expected.len());
            for (gene, row) in expected.iter().enumerate() {
                assert_eq!(actual.gene_ids[gene], row["gene_id"]);
                for (sample, id) in actual.sample_ids.iter().enumerate() {
                    assert_eq!(
                        actual.counts[gene * 4 + sample],
                        row[id].parse::<u32>()?,
                        "{case} pooled {mode} {gene} {id}"
                    );
                }
            }
            if let Some(factors) = actual.normalization_factors {
                for (gene, row) in read_table(&oracle.join(case).join("length_factors.tsv"))?
                    .iter()
                    .enumerate()
                {
                    for (sample, id) in actual.sample_ids.iter().enumerate() {
                        assert_close(
                            factors[gene * 4 + sample],
                            row[id].parse()?,
                            1e-10,
                            1e-10,
                            &format!("{case} pooled length {gene} {id}"),
                        );
                    }
                }
            } else {
                assert_ne!(mode, "no");
            }
        }
    }
    Ok(())
}
#[test]
fn all_gene_import_modes_rounding_and_offsets_match_pinned_tximport_deseq2() -> Result<()> {
    let oracle = reference("tximport")?;
    for case in ["fractional_and_zero_abundance", "zero_length_without_mass"] {
        let samples = (1..=4)
            .map(|i| {
                (
                    format!("sample{i}"),
                    oracle.join(case).join(format!("sample{i}.tsv")),
                )
            })
            .collect::<Vec<_>>();
        for mode in ["no", "scaledTPM", "lengthScaledTPM"] {
            let imported =
                import_kallisto_abundances(&samples, &oracle.join(case).join("tx2gene.tsv"), mode)?;
            let expected = read_table(&oracle.join(case).join(format!("{mode}_counts.tsv")))?;
            assert_eq!(imported.gene_ids.len(), expected.len());
            for (gene, row) in imported.gene_ids.iter().zip(&expected).enumerate() {
                assert_eq!(row.0, row.1.get("gene_id").unwrap());
                for (sample, id) in imported.sample_ids.iter().enumerate() {
                    assert_eq!(
                        imported.counts[gene * 4 + sample],
                        row.1[id].parse::<u32>()?,
                        "{case} {mode} {} {id}",
                        row.0
                    );
                }
            }
            if mode == "no" {
                let mut factors = imported
                    .normalization_factors
                    .context("missing length factors")?;
                for (gene, row) in read_table(&oracle.join(case).join("length_factors.tsv"))?
                    .iter()
                    .enumerate()
                {
                    for (sample, id) in imported.sample_ids.iter().enumerate() {
                        assert_close(
                            factors[gene * 4 + sample],
                            row[id].parse()?,
                            1e-10,
                            1e-10,
                            &format!("{case} length factors {gene} {id}"),
                        );
                    }
                }
                complete_length_offsets(&imported.counts, &mut factors, 4)?;
                let normalized = read_table(&oracle.join(case).join("normalized_counts.tsv"))?;
                for (gene, row) in read_table(&oracle.join(case).join("normalization_factors.tsv"))?
                    .iter()
                    .enumerate()
                {
                    for (sample, id) in imported.sample_ids.iter().enumerate() {
                        assert_close(
                            factors[gene * 4 + sample],
                            row[id].parse()?,
                            1e-10,
                            1e-10,
                            &format!("{case} complete offsets {gene} {id}"),
                        );
                        assert_close(
                            imported.counts[gene * 4 + sample] as f64 / factors[gene * 4 + sample],
                            normalized[gene][id].parse()?,
                            1e-9,
                            1e-10,
                            &format!("{case} normalized counts {gene} {id}"),
                        );
                    }
                }
            } else {
                assert!(imported.normalization_factors.is_none());
            }
        }
    }
    Ok(())
}
