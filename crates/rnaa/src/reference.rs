//! Reference admission shared by quantification and cohort analysis.
use anyhow::{Context, Result, bail};
use rnaa_core::{
    Database, ProjectConfig, ReferenceBundle,
    references::{ReferenceScope, ReferenceSelection},
};

pub(crate) fn pin(
    db: &Database,
    bundle: &ReferenceBundle,
    config: &ProjectConfig,
    scope: ReferenceScope,
) -> Result<ReferenceSelection> {
    let selected = db
        .selected_reference()?
        .context("no verified reference is selected; prepare references first")?;
    if bundle.id != selected.reference_id
        || bundle.kallisto_index_path.to_string_lossy() != selected.index.path
        || bundle.tx2gene_path.to_string_lossy() != selected.mapping.path
        || bundle.cdna_path.to_string_lossy() != selected.cdna.path
        || bundle.gtf_path.to_string_lossy() != selected.gtf.path
        || bundle
            .gene_annotation_path
            .as_ref()
            .map(|path| path.to_string_lossy().into_owned())
            .as_deref()
            != Some(selected.annotation.path.as_str())
    {
        bail!("selected reference changed while preparing the workflow; retry");
    }
    selected.verify_catalog(db, config, scope)?;
    Ok(selected)
}
pub(crate) fn bind<P>(
    request: &mut rnaa_stages::analysis_request::AnalysisRequest<P>,
    selection: &ReferenceSelection,
) -> Result<()> {
    request.request.dependencies.insert(
        "reference".into(),
        selection.fingerprint(ReferenceScope::Analysis)?,
    );
    Ok(())
}

pub(crate) struct Status {
    pub(crate) selection: Option<ReferenceSelection>,
    pub(crate) reason: Option<String>,
}
pub(crate) fn refresh_status(db: &Database, config: &ProjectConfig) -> Result<Status> {
    let selection = db.selected_reference()?;
    let reason = match &selection {
        Some(value) => value
            .verify_current(db, config, ReferenceScope::Analysis)
            .err()
            .map(|error| format!("{error:#}")),
        None => Some("no verified reference selected".into()),
    };
    if let Some(reason) = &reason {
        db.with_idle_analysis(|db| {
            let disk = ProjectConfig::load(
                &std::path::Path::new(&db.project()?.root_dir).join("rnaa.toml"),
            )?;
            if rnaa_core::references::ReferencePolicy::from_config(&disk)
                != rnaa_core::references::ReferencePolicy::from_config(config)
            {
                return Ok(());
            }
            db.atomic(|db| {
                if db
                    .selected_reference()?
                    .map(|value| value.generation.generation_id)
                    != selection
                        .as_ref()
                        .map(|value| value.generation.generation_id.clone())
                {
                    return Ok(());
                }
                let mut changed = false;
                for stage in ["normalize", "deseq2", "corr"] {
                    if db
                        .get_project_stage_state(stage)?
                        .is_some_and(|value| value.state == rnaa_core::ProjectStageState::Done)
                    {
                        db.set_project_stage_state(
                            stage,
                            rnaa_core::ProjectStageState::Pending,
                            Some(reason),
                        )?;
                        changed = true;
                    }
                }
                if changed {
                    db.append_event(
                        "refs",
                        None,
                        "reference validity changed",
                        serde_json::json!({"reason":reason}),
                    )?;
                }
                Ok(())
            })
        })?;
    }
    Ok(Status { selection, reason })
}
