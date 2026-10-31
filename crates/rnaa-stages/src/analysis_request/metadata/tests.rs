use super::*;
use rnaa_core::{
    ContrastSpec, ManifestArtifact,
    fingerprint::StageRequest,
    metadata::{ColumnSpec, ColumnType},
};
use std::fs;

fn sheet(path: &Path, decimal: bool, changed_note: bool) -> Result<()> {
    let mut text = "run_accession\tcondition\tage\tnotes\tread1_path\n".to_string();
    for i in 0..8 {
        text.push_str(&format!(
            "RUN{i}\t{}\t{}{}\t{}\t/location/{i}\n",
            if i < 4 { "01" } else { "1" },
            20 + i % 4,
            if decimal { ".0" } else { "" },
            if changed_note { "changed" } else { "original" }
        ));
    }
    fs::write(path, text)?;
    Ok(())
}
fn config() -> ProjectConfig {
    let mut config = ProjectConfig::default();
    config
        .metadata
        .columns
        .insert("age".into(), ColumnSpec::new(ColumnType::Numeric));
    config
}

#[test]
fn scientific_projection_keeps_normalized_used_values_and_effective_factor_coding() -> Result<()> {
    let temp = tempfile::tempdir()?;
    let path = temp.path().join("metadata.tsv");
    let mut config = config();
    let contrasts = vec![ContrastSpec::new(
        None,
        "condition".into(),
        "1".into(),
        "01".into(),
    )];
    sheet(&path, false, false)?;
    let (plan, original) = prepare_de(&path, "~ age + condition", &contrasts, &config)?;
    sheet(&path, true, true)?;
    config
        .metadata
        .columns
        .insert("notes".into(), ColumnSpec::new(ColumnType::Ignore));
    config
        .metadata
        .columns
        .insert("unused".into(), ColumnSpec::new(ColumnType::Numeric));
    config.metadata.minimum_replicates = 3;
    config.metadata.missing_values.push("not collected".into());
    let (changed_plan, equivalent) = prepare_de(&path, "~ age + condition", &contrasts, &config)?;
    assert_eq!(original, equivalent);
    assert_eq!(plan.membership(), changed_plan.membership());
    assert_ne!(
        plan, changed_plan,
        "exact source metadata remains distinct provenance"
    );
    config.metadata.columns.insert(
        "condition".into(),
        ColumnSpec {
            kind: ColumnType::Factor,
            levels: vec!["01".into(), "1".into()],
            reference: Some("1".into()),
        },
    );
    let (_, relevelled) = prepare_de(&path, "~ age + condition", &contrasts, &config)?;
    assert_ne!(original, relevelled);
    config.metadata.minimum_replicates = 5;
    assert!(
        prepare_de(&path, "~ age + condition", &contrasts, &config).is_err(),
        "a reused request must still pass current replication policy"
    );
    Ok(())
}

#[test]
fn fixed_formulas_ignore_unused_covariates_but_dot_expansion_consumes_them() -> Result<()> {
    let temp = tempfile::tempdir()?;
    let path = temp.path().join("metadata.tsv");
    sheet(&path, false, false)?;
    let config = config();
    let (_, fixed) = prepare_de(&path, "~ condition", &[], &config)?;
    let (_, expanded) = prepare_de(&path, "~ .", &[], &config)?;
    assert_eq!(fixed.columns.len(), 1);
    assert_eq!(expanded.columns.len(), 2);
    let text = fs::read_to_string(&path)?.replace("\t20\t", "\t21\t");
    fs::write(&path, text)?;
    assert_eq!(fixed, prepare_de(&path, "~ condition", &[], &config)?.1);
    assert_ne!(expanded, prepare_de(&path, "~ .", &[], &config)?.1);
    Ok(())
}

#[test]
fn raw_witnesses_are_retained_without_changing_scientific_identity_or_overwriting_history()
-> Result<()> {
    let temp = tempfile::tempdir()?;
    let path = temp.path().join("samplesheet.tsv");
    fs::write(&path, "original source")?;
    let source = ManifestArtifact::snapshot("SAMPLESHEET_SOURCE", &path)?;
    let make = |input| AnalysisRequest {
        request: StageRequest {
            schema_version: 1,
            stage: "test".into(),
            algorithm_version: 1,
            engines: Default::default(),
            parameters: true,
            inputs: scientific_inputs(std::slice::from_ref(&input)),
            dependencies: Default::default(),
        },
        inputs: vec![input],
        policy: None,
        observations: Vec::new(),
    };
    let old = make(source);
    let root = temp.path().join("observations");
    let retained = old.retain_metadata_witnesses(&root)?;
    fs::write(&path, "a later source with unchanged model values")?;
    let new = make(ManifestArtifact::snapshot("SAMPLESHEET_SOURCE", &path)?);
    assert_eq!(old.request.digest()?, new.request.digest()?);
    assert!(
        old.verify_inputs().is_err(),
        "raw source changes must still fence a running attempt"
    );
    assert!(retained[0].matches_file()?);
    let later = new.retain_metadata_witnesses(&root)?;
    assert_ne!(later[0].path, retained[0].path);
    fs::write(&path, "original source")?;
    fs::write(&retained[0].path, "corrupt observation")?;
    assert!(old.retain_metadata_witnesses(&root).is_err());
    assert_eq!(
        fs::read_to_string(&retained[0].path)?,
        "corrupt observation"
    );
    Ok(())
}
