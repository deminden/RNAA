use super::*;

#[test]
fn groups_runs_in_deterministic_sample_order_and_preserves_membership() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("metadata.tsv");
    std::fs::write(
        &path,
        "run_accession\tsample_accession\tcondition\tbatch\tlane\nR2\tBIO1\tA\tB1\t2\nR3\tBIO2\tB\tB1\t1\nR1\tBIO1\tA\tB1\t1\n",
    )?;
    let plan = SamplePlan::from_sheet(&path, "~ batch + condition", &[])?;
    assert_eq!(plan.sample_ids(), ["BIO1", "BIO2"]);
    assert_eq!(plan.samples[0].runs, ["R1", "R2"]);
    assert_eq!(plan.samples[0].metadata["lane"], "");
    assert_eq!(plan.samples[0].metadata["condition"], "A");
    let output = dir.path().join("sample_metadata.tsv");
    plan.write_sheet(&output)?;
    let text = std::fs::read_to_string(output)?;
    assert!(text.contains("sample_id\tbiological_sample_id\tidentity_source\trun_accessions"));
    let encoded = serde_json::to_vec(&plan)?;
    assert_eq!(plan, serde_json::from_slice::<SamplePlan>(&encoded)?);
    Ok(())
}

#[test]
fn conflicting_model_values_and_single_biological_samples_are_rejected() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("metadata.tsv");
    std::fs::write(
        &path,
        "run_accession\tsample_accession\tcondition\nR1\tBIO1\tA\nR2\tBIO1\tB\nR3\tBIO2\tB\n",
    )?;
    assert!(
        SamplePlan::from_sheet(&path, "~ condition", &[])
            .unwrap_err()
            .to_string()
            .contains("conflicting")
    );
    std::fs::write(
        &path,
        "run_accession\tsample_accession\tcondition\nR1\tBIO1\tA\nR2\tBIO1\tA\n",
    )?;
    assert!(
        SamplePlan::from_sheet(&path, "~ condition", &[])
            .unwrap_err()
            .to_string()
            .contains("found 1")
    );
    Ok(())
}

#[test]
fn duplicate_rows_headers_and_identity_covariates_are_rejected() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("metadata.tsv");
    for text in [
        "run_accession\tcondition\tcondition\nR1\tA\tB\n",
        "run_accession\tcondition\nR1\tA\nR1\tB\n",
    ] {
        std::fs::write(&path, text)?;
        assert!(SamplePlan::from_sheet(&path, "~ condition", &[]).is_err());
    }
    std::fs::write(&path, "run_accession\tcondition\nR1\tA\nR2\tB\n")?;
    assert!(SamplePlan::from_sheet(&path, "~ run_accession", &[]).is_err());
    let plan = SamplePlan::from_sheet(&path, "~ condition", &[])?;
    assert_eq!(plan.samples[0].identity_source, IdentitySource::RunFallback);
    Ok(())
}

#[test]
fn typed_numeric_values_agree_across_technical_runs_without_merging_factor_labels() -> Result<()> {
    use crate::metadata::{ColumnSpec, ColumnType, MetadataConfig};
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("metadata.tsv");
    std::fs::write(
        &path,
        "run_accession\tsample_accession\tage\tcondition\nR1\tBIO1\t21\t01\nR2\tBIO1\t21.0\t01\nR3\tBIO2\t22\t1\n",
    )?;
    let mut schema = MetadataConfig::default();
    schema
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    let plan = SamplePlan::from_sheet_with_metadata(&path, "~ age + condition", &[], &schema)?;
    assert_eq!(plan.samples[0].metadata["age"], "21");
    assert_eq!(plan.samples[0].metadata["condition"], "01");
    assert_eq!(plan.samples[1].metadata["condition"], "1");
    Ok(())
}
