//! Scientific setting observations fence count publication independently of threads.
use super::quant_request;
use anyhow::{Context, Result, bail};
use rnaa_core::{Database, ProjectConfig, ReferenceBundle, RunRecord, eligibility::QcPolicy};

pub fn quant_policy_fingerprint(
    run: &RunRecord,
    reference: &ReferenceBundle,
    config: &ProjectConfig,
) -> Result<String> {
    quant_request(run, reference, config, &[])?.digest()
}

fn matches(
    expected_run: &RunRecord,
    expected: &ProjectConfig,
    current_run: &RunRecord,
    current: &ProjectConfig,
) -> bool {
    let single = expected_run.library_layout == rnaa_core::LibraryLayout::Single;
    expected_run.library_layout == current_run.library_layout
        && expected.download.prefer == current.download.prefer
        && expected.quant.engine == current.quant.engine
        && expected.quant.strand == current.quant.strand
        && expected.quant.preprocess == current.quant.preprocess
        && (!single
            || expected.quant.single_fragment_length == current.quant.single_fragment_length
                && expected.quant.single_fragment_sd == current.quant.single_fragment_sd)
        && QcPolicy::from_config(expected) == QcPolicy::from_config(current)
}

pub fn verify_quant_catalog(db: &Database, run: &RunRecord, config: &ProjectConfig) -> Result<()> {
    let current: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
    let latest = db
        .get_run(&run.run_accession)?
        .context("quantification run disappeared")?;
    if !matches(run, config, &latest, &current)
        || crate::download::ReqwestDownloader.source_request_digest(run, config.download.prefer)?
            != crate::download::ReqwestDownloader
                .source_request_digest(&latest, current.download.prefer)?
    {
        bail!(
            "quantification settings, source selection, layout or QC admission changed; retry with current settings"
        );
    }
    Ok(())
}

pub fn source_binding_current(
    db: &Database,
    run: &RunRecord,
    config: &ProjectConfig,
) -> Result<bool> {
    match db.generation_head(
        rnaa_core::generations::OutputGroup::Download,
        &run.run_accession,
    )? {
        Some(head) => Ok(head.request_sha256
            == crate::download::ReqwestDownloader
                .source_request_digest(run, config.download.prefer)?),
        None => Ok(run.remote_files.is_empty()),
    }
}
pub fn verify_quant_current(db: &Database, run: &RunRecord, config: &ProjectConfig) -> Result<()> {
    verify_quant_catalog(db, run, config)?;
    let current =
        ProjectConfig::load(&std::path::Path::new(&db.project()?.root_dir).join("rnaa.toml"))?;
    if !matches(run, config, run, &current) {
        bail!(
            "quantification settings or QC admission changed on disk; retry with current settings"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests;
