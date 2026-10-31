//! Admit the selected DE matrix and prepare current metadata for its consumer.
use super::inputs::write_filtered_deseq2_samplesheet;
use crate::quality;
use anyhow::{Context, Result, bail};
use rnaa_core::{
    Database, ManifestArtifact, ProjectConfig, ProjectPaths, RunRecord,
    references::{ReferenceScope, ReferenceSelection},
    sample_plan::SamplePlan,
};
use rnaa_stages::analysis_request::DeOptions;
use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

pub(crate) fn load_descriptor(
    paths: &ProjectPaths,
    db: &Database,
    reference: &ReferenceSelection,
) -> Result<serde_json::Value> {
    let descriptor: serde_json::Value = serde_json::from_slice(&fs::read(
        paths
            .de_project_dir(db.project_id())
            .join("rnaa_deseq2_request.json"),
    )?)?;
    if descriptor["request"]["dependencies"]["reference"]
        != reference.fingerprint(ReferenceScope::Analysis)?
    {
        bail!("DE matrix reference is stale; rerun normalize or deseq2 before correlation");
    }
    Ok(descriptor)
}

pub(crate) fn current_metadata(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    eligibility: &quality::Snapshot,
    descriptor: serde_json::Value,
    runs: &[RunRecord],
    staging: &Path,
) -> Result<(PathBuf, ManifestArtifact)> {
    let options: DeOptions = serde_json::from_value(descriptor["request"]["parameters"].clone())
        .context("DE model descriptor is stale; rerun normalize or deseq2 before correlation")?;
    let samples = descriptor["request"]["parameters"]["samples"]
        .as_array()
        .context("DE cohort descriptor missing; rerun normalize or deseq2")?
        .iter()
        .map(|sample| {
            sample
                .as_str()
                .map(str::to_string)
                .context("invalid DE cohort descriptor")
        })
        .collect::<Result<BTreeSet<_>>>()?;
    let expected = db
        .list_active_runs()?
        .iter()
        .filter(|run| eligibility.eligible.contains(&run.run_accession))
        .map(|run| Ok(rnaa_core::cohorts::BiologicalSample::for_run(run)?.id))
        .collect::<Result<BTreeSet<_>>>()?;
    if descriptor["request"]["dependencies"]["eligibility"] != eligibility.digest()?
        || samples != expected
    {
        bail!(
            "DE matrix cohort or QC policy is stale; rerun normalize or deseq2 before correlation"
        );
    }
    let filtered = staging.join("samplesheet_quantified.tsv");
    let source = write_filtered_deseq2_samplesheet(paths, runs, &filtered)?;
    let stale = "DE metadata is stale; DE model settings or metadata schema are stale or invalid; rerun normalize or deseq2 before correlation";
    if !options.matches_metadata(&filtered, config).context(stale)? {
        bail!("{stale}");
    }
    // The producer only needs its own design variables. Adjustment may consume
    // different covariates, so build it from current metadata after producer proof.
    let plan = SamplePlan::from_sheet_with_metadata(
        &filtered,
        if config.corr.adjust == "none" {
            "~ 1"
        } else {
            &config.corr.model
        },
        &[],
        &config.metadata,
    )?;
    let sheet = staging.join("sample_metadata.tsv");
    plan.write_sheet(&sheet)?;
    Ok((sheet, source))
}
