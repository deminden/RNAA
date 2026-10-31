use super::*;
use rnaa_core::{
    ArtifactKind, ArtifactRecord, ManifestArtifact,
    fingerprint::StageRequest,
    generations::{Generation, OutputGroup},
    jobs::{JobLease, JobOutcome},
    util::now_rfc3339,
};

fn plan(config: &ProjectConfig) -> AnalysisRequest<bool> {
    AnalysisRequest {
        observations: Vec::new(),
        request: StageRequest {
            schema_version: 1,
            stage: "normalize".into(),
            algorithm_version: 1,
            engines: Default::default(),
            parameters: true,
            inputs: Vec::new(),
            dependencies: Default::default(),
        },
        inputs: Vec::new(),
        policy: Some(ObservedConfig::de(config, &config.deseq2.design)),
    }
}

#[test]
fn upstream_quantification_policy_edits_fence_both_de_and_correlation_publication() -> Result<()> {
    let temp = tempfile::tempdir()?;
    let mut config = ProjectConfig::default();
    let db = Database::create(temp.path(), &config)?;
    config.save(&temp.path().join("rnaa.toml"))?;
    let de = plan(&config);
    let mut corr = plan(&config);
    corr.policy = Some(ObservedConfig::corr(&config));
    config.quant.threads += 1;
    db.set_project_config(&config)?;
    config.save(&temp.path().join("rnaa.toml"))?;
    de.verify_current_settings(&db)?;
    corr.verify_current_settings(&db)?;
    config.quant.single_fragment_length += 10.0;
    db.set_project_config(&config)?;
    config.save(&temp.path().join("rnaa.toml"))?;
    assert!(de.verify_current_settings(&db).is_err());
    assert!(corr.verify_current_settings(&db).is_err());
    Ok(())
}
#[test]
fn model_policy_changes_cannot_publish_a_sealed_old_design_but_threads_can_change() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let mut config = ProjectConfig::default();
    config.deseq2.design = "~ condition".into();
    let db = Database::create(dir.path(), &config)?;
    config.save(&dir.path().join("rnaa.toml"))?;
    let request = plan(&config);
    request.verify_current_settings(&db)?;
    config.quant.threads += 1;
    config.corr.threads += 1;
    db.set_project_config(&config)?;
    config.save(&dir.path().join("rnaa.toml"))?;
    request.verify_current_settings(&db)?;
    let mut lease = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    let owned = lease.database();
    let digest = request.request.digest()?;
    owned.set_job_request(&digest)?;
    let mut generation = Generation::begin(&owned, OutputGroup::De, "normalize")?;
    let path = generation.head().payload_path.join("counts.tsv");
    std::fs::write(&path, "old design result")?;
    let file = ManifestArtifact::snapshot("COUNTS", &path)?;
    let record = ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: None,
        kind: ArtifactKind::Counts,
        path: file.path,
        blob_id: None,
        shared_path: None,
        checksum_type: file.checksum_type,
        checksum: file.checksum,
        bytes: file.bytes,
        created_at: now_rfc3339(),
    };
    generation.seal(&[path], &digest)?;
    config.metadata.minimum_replicates = 3;
    db.set_project_config(&config)?;
    assert!(
        generation
            .commit(&[record], |db| request.verify_settings_catalog(db))
            .is_err()
    );
    assert!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .is_none()
    );
    assert!(db.list_artifacts()?.is_empty());
    lease.finish(JobOutcome::Failed, Some("model policy changed"))?;
    Ok(())
}
#[test]
fn disk_edits_and_corr_input_path_changes_are_observed_without_changing_thread_policy() -> Result<()>
{
    let dir = tempfile::tempdir()?;
    let mut config = ProjectConfig::default();
    let db = Database::create(dir.path(), &config)?;
    config.save(&dir.path().join("rnaa.toml"))?;
    let mut request = plan(&config);
    request.policy = Some(ObservedConfig::corr(&config));
    request.verify_current_settings(&db)?;
    config.corr.threads = 7;
    config.save(&dir.path().join("rnaa.toml"))?;
    request.verify_current_settings(&db)?;
    config.corr.geneset = "genes:new-input.tsv".into();
    config.save(&dir.path().join("rnaa.toml"))?;
    assert!(request.verify_current_settings(&db).is_err());
    request.verify_settings_catalog(&db)?;
    Ok(())
}

#[test]
fn module_options_fence_correlation_publication_but_not_de() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let mut config = ProjectConfig::default();
    config.corr.modules = Some(rnaa_core::config::CorrModuleConfig {
        gmt_path: dir.path().join("pathways.gmt"),
        top_genes: 500,
        min_module_size: 25,
        permutations: 1000,
        seed: 42,
        pathway_min_size: 10,
        pathway_max_size: None,
    });
    let db = Database::create(dir.path(), &config)?;
    config.save(&dir.path().join("rnaa.toml"))?;
    let de = plan(&config);
    let mut corr = plan(&config);
    corr.policy = Some(ObservedConfig::corr(&config));
    corr.verify_current_settings(&db)?;
    config.corr.modules.as_mut().unwrap().seed += 1;
    db.set_project_config(&config)?;
    config.save(&dir.path().join("rnaa.toml"))?;
    assert!(corr.verify_current_settings(&db).is_err());
    de.verify_current_settings(&db)?;
    Ok(())
}

#[test]
fn changing_external_matrix_selection_fences_correlation_publication() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let mut config = ProjectConfig::default();
    config.corr.matrix = Some(dir.path().join("first.tsv"));
    let db = Database::create(dir.path(), &config)?;
    config.save(&dir.path().join("rnaa.toml"))?;
    let mut request = plan(&config);
    request.policy = Some(ObservedConfig::corr(&config));
    request.verify_current_settings(&db)?;
    config.corr.matrix = Some(dir.path().join("second.tsv"));
    db.set_project_config(&config)?;
    config.save(&dir.path().join("rnaa.toml"))?;
    assert!(request.verify_current_settings(&db).is_err());
    Ok(())
}
