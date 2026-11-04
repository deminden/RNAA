use super::*;

#[test]
fn upgraded_engine_cannot_reuse_prior_preprocessing_identity() -> anyhow::Result<()> {
    let run = RunRecord {
        project_id: "test".into(),
        study_accession: "study".into(),
        source_study_accession: None,
        run_accession: "RUN1".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: Vec::new(),
        metadata: Default::default(),
        state: rnaa_core::RunState::Verified,
        last_error: None,
        updated_at: String::new(),
    };
    let current = request_from_inputs(&run, &[], &ProjectConfig::default())?;
    let mut old = current.clone();
    old.algorithm_version = 3;
    old.engines.insert("fasterp".into(), "0.2.1".into());
    assert_ne!(current.digest()?, old.digest()?);
    assert_eq!(current.engines["fasterp"], "0.2.1+git.dcca6ad25e23");

    let mut manifest = StageManifest {
        schema_version: 1,
        stage: "preprocess".into(),
        id: run.run_accession.clone(),
        started_at: String::new(),
        finished_at: String::new(),
        parameters: scientific_parameters(&run),
        tool_versions: Default::default(),
        input_artifacts: Vec::new(),
        output_artifacts: Vec::new(),
        notes: Vec::new(),
    };
    assert!(processing_identity_matches(&run, &manifest));
    manifest.parameters["algorithm_version"] = serde_json::json!(3);
    manifest.parameters["engine_version"] = serde_json::json!("0.2.1");
    assert!(!processing_identity_matches(&run, &manifest));
    Ok(())
}
