use super::*;

fn table(path: &Path, rows: &str) -> Result<()> {
    std::fs::write(
        path,
        format!("target_id\tlength\teff_length\test_counts\ttpm\n{rows}"),
    )?;
    Ok(())
}
fn contribution(id: &str, path: &Path) -> RunAbundance {
    RunAbundance {
        run_id: id.into(),
        path: path.into(),
    }
}

#[test]
fn pools_fractional_counts_and_length_exposure_before_rounding_and_scaling() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let a = dir.path().join("a.tsv");
    let b = dir.path().join("b.tsv");
    let c = dir.path().join("c.tsv");
    table(
        &a,
        "t1\t1000\t100\t0.25\t1\nt2\t1000\t200\t1.25\t1\nt3\t1000\t100\t10\t1\n",
    )?;
    table(
        &b,
        "t1\t1000\t50\t0.25\t1\nt2\t1000\t300\t0.75\t1\nt3\t1000\t200\t20\t1\n",
    )?;
    table(
        &c,
        "t1\t1000\t100\t5\t181818.18181818182\nt2\t1000\t200\t5\t90909.09090909091\nt3\t1000\t150\t30\t727272.7272727273\n",
    )?;
    let mapping = dir.path().join("mapping.tsv");
    std::fs::write(&mapping, "transcript_id\tgene_id\nt1\tG1\nt2\tG1\nt3\tG2\n")?;
    let samples = vec![
        QuantifiedSample {
            sample_id: "A".into(),
            runs: vec![contribution("R2", &b), contribution("R1", &a)],
        },
        QuantifiedSample {
            sample_id: "B".into(),
            runs: vec![contribution("R3", &c)],
        },
    ];
    let raw = import_grouped_kallisto_abundances(&samples, &mapping, "no")?;
    assert_eq!(raw.sample_ids, ["A", "B"]);
    assert_eq!(
        raw.counts,
        [2, 10, 30, 30],
        "round only after pooling fractional transcript and gene counts"
    );
    let factors = raw.normalization_factors.unwrap();
    assert!((factors[0] - (15.0_f64 / 13.0).sqrt()).abs() < 1e-10);
    assert!((factors[1] - (13.0_f64 / 15.0).sqrt()).abs() < 1e-10);
    assert!(factors[2..].iter().all(|value| (value - 1.0).abs() < 1e-12));
    assert_eq!(
        import_grouped_kallisto_abundances(&samples, &mapping, "scaledTPM")?.counts,
        [2, 11, 30, 29]
    );
    assert_eq!(
        import_grouped_kallisto_abundances(&samples, &mapping, "lengthScaledTPM")?.counts,
        [2, 11, 30, 29]
    );
    // Different run traversal order cannot change the result.
    let mut reordered = samples.clone();
    reordered.reverse();
    reordered[1].runs.reverse();
    assert_eq!(
        import_grouped_kallisto_abundances(&reordered, &mapping, "lengthScaledTPM")?.counts,
        import_grouped_kallisto_abundances(&samples, &mapping, "lengthScaledTPM")?.counts
    );
    Ok(())
}

#[test]
fn pooling_rejects_duplicate_runs_mismatched_references_and_invalid_abundances() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let a = dir.path().join("a.tsv");
    let b = dir.path().join("b.tsv");
    let mapping = dir.path().join("mapping.tsv");
    std::fs::write(&mapping, "transcript_id\tgene_id\nt1\tG1\n")?;
    table(&a, "t1\t1000\t100\t1\t1000000\n")?;
    let sample = QuantifiedSample {
        sample_id: "A".into(),
        runs: vec![contribution("R1", &a), contribution("R2", &b)],
    };
    for rows in [
        "other\t1000\t100\t1\t1000000\n",
        "t1\t999\t100\t1\t1000000\n",
        "t1\t1000\t100\tNaN\t1\n",
        "t1\t1000\t0\t1\t1\n",
    ] {
        table(&b, rows)?;
        assert!(pool(&sample, &dir.path().join("output.tsv")).is_err());
    }
    let mut duplicate = sample.clone();
    duplicate.sample_id = "B".into();
    assert!(
        import_grouped_kallisto_abundances(&[sample, duplicate], &mapping, "no")
            .unwrap_err()
            .to_string()
            .contains("exactly one")
    );
    Ok(())
}

#[test]
fn splitting_a_library_across_runs_does_not_change_the_biological_import() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let full = dir.path().join("full.tsv");
    let a = dir.path().join("a.tsv");
    let b = dir.path().join("b.tsv");
    let other = dir.path().join("other.tsv");
    table(
        &full,
        "t1\t1000\t100\t2\t333333.3333333333\nt2\t1000\t200\t8\t666666.6666666666\n",
    )?;
    table(
        &a,
        "t1\t1000\t100\t0.5\t333333.3333333333\nt2\t1000\t200\t2\t666666.6666666666\n",
    )?;
    table(
        &b,
        "t1\t1000\t100\t1.5\t333333.3333333333\nt2\t1000\t200\t6\t666666.6666666666\n",
    )?;
    table(
        &other,
        "t1\t1000\t100\t6\t750000\nt2\t1000\t200\t4\t250000\n",
    )?;
    let mapping = dir.path().join("mapping.tsv");
    std::fs::write(&mapping, "transcript_id\tgene_id\nt1\tG1\nt2\tG2\n")?;
    let reference = vec![
        QuantifiedSample {
            sample_id: "A".into(),
            runs: vec![contribution("FULL", &full)],
        },
        QuantifiedSample {
            sample_id: "B".into(),
            runs: vec![contribution("OTHER", &other)],
        },
    ];
    let mut split = reference.clone();
    split[0].runs = vec![contribution("R1", &a), contribution("R2", &b)];
    for mode in ["no", "scaledTPM", "lengthScaledTPM"] {
        let expected = import_grouped_kallisto_abundances(&reference, &mapping, mode)?;
        let actual = import_grouped_kallisto_abundances(&split, &mapping, mode)?;
        assert_eq!(actual.sample_ids, expected.sample_ids);
        assert_eq!(actual.counts, expected.counts);
        if let (Some(a), Some(b)) = (actual.normalization_factors, expected.normalization_factors) {
            assert!(a.iter().zip(b).all(|(a, b)| (a - b).abs() < 1e-12));
        }
    }
    Ok(())
}

#[test]
fn unexpressed_targets_keep_positive_effective_lengths_without_inventing_mass() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let a = dir.path().join("a.tsv");
    let b = dir.path().join("b.tsv");
    let output = dir.path().join("pooled.tsv");
    table(&a, "t1\t1000\t100\t0\t0\n")?;
    table(&b, "t1\t1000\t300\t0\t0\n")?;
    pool(
        &QuantifiedSample {
            sample_id: "A".into(),
            runs: vec![contribution("R1", &a), contribution("R2", &b)],
        },
        &output,
    )?;
    let row: Row = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(output)?
        .deserialize()
        .next()
        .unwrap()?;
    assert_eq!(row.eff_length, 200.0);
    assert_eq!(row.est_counts, 0.0);
    assert_eq!(row.tpm, 0.0);
    Ok(())
}
