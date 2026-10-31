//! Typed model/membership projection; raw tables remain immutable witnesses.
use super::{AnalysisRequest, DeOptions};
use crate::design::{PreparedDesign, Purpose, ScientificDesign};
use anyhow::Result;
use rnaa_core::{ContrastSpec, ProjectConfig, sample_plan::SamplePlan};
use std::path::Path;

pub(super) fn prepare_de(
    sheet: &Path,
    design: &str,
    contrasts: &[ContrastSpec],
    config: &ProjectConfig,
) -> Result<(SamplePlan, ScientificDesign)> {
    let plan = SamplePlan::from_sheet_with_metadata(sheet, design, contrasts, &config.metadata)?;
    let (headers, rows) = plan.model_table()?;
    let prepared = PreparedDesign::new(
        &headers,
        &rows,
        &plan.sample_ids(),
        design,
        &config.metadata,
        contrasts,
        Purpose::DifferentialExpression,
    )?;
    Ok((plan, prepared.contract.scientific()))
}

impl DeOptions {
    /// Revalidate actual current metadata before admitting an existing DE matrix.
    /// This also checks current schema constraints, even when they do not alter X.
    pub fn matches_metadata(&self, sheet: &Path, config: &ProjectConfig) -> Result<bool> {
        if self.design != config.deseq2.design
            || self.transform != config.deseq2.transform
            || self.counts_from_abundance != config.deseq2.counts_from_abundance
        {
            return Ok(false);
        }
        let (plan, design) = prepare_de(sheet, &self.design, &self.contrasts, config)?;
        Ok(self.sample_plan == plan.membership() && self.design_contract == design)
    }
}

pub(super) fn scientific_inputs(
    inputs: &[rnaa_core::ManifestArtifact],
) -> Vec<rnaa_core::fingerprint::InputIdentity> {
    inputs
        .iter()
        .filter(|input| !metadata_role(&input.kind))
        .map(rnaa_core::fingerprint::InputIdentity::from)
        .collect()
}

fn metadata_role(role: &str) -> bool {
    matches!(role, "SAMPLESHEET" | "METADATA" | "SAMPLESHEET_SOURCE")
}

impl<P> AnalysisRequest<P> {
    /// Retain current-attempt evidence without rewriting the reused generation.
    pub fn retain_metadata_witnesses(
        &self,
        root: &Path,
    ) -> Result<Vec<rnaa_core::ManifestArtifact>> {
        use anyhow::{Context, bail};
        use rnaa_core::{ManifestArtifact, fingerprint::InputIdentity};
        let mut retained = Vec::new();
        for source in self
            .inputs
            .iter()
            .filter(|input| metadata_role(&input.kind))
        {
            if !source.matches_file()? {
                bail!("metadata changed before retaining reuse evidence");
            }
            let dir = root.join(&source.checksum).join(&source.kind);
            std::fs::create_dir_all(&dir)?;
            let path = dir.join(
                Path::new(&source.path)
                    .file_name()
                    .context("metadata witness has no file name")?,
            );
            if !path.exists() {
                let staged = tempfile::NamedTempFile::new_in(&dir)?;
                std::fs::copy(&source.path, staged.path())?;
                let copied = ManifestArtifact::snapshot(&source.kind, staged.path())?;
                if copied.checksum != source.checksum
                    || copied.bytes != source.bytes
                    || !source.matches_file()?
                {
                    bail!("metadata changed while retaining reuse evidence");
                }
                staged.as_file().sync_all()?;
                if let Err(error) = staged.persist_noclobber(&path)
                    && error.error.kind() != std::io::ErrorKind::AlreadyExists
                {
                    return Err(error.into());
                }
                #[cfg(unix)]
                std::fs::File::open(&dir)?.sync_all()?;
            }
            let witness = ManifestArtifact::snapshot(&source.kind, &path)?;
            if InputIdentity::from(&witness) != InputIdentity::from(source) {
                bail!(
                    "retained metadata observation is corrupt: {}",
                    path.display()
                );
            }
            retained.push(witness);
        }
        for source in self
            .inputs
            .iter()
            .filter(|input| metadata_role(&input.kind))
        {
            if !source.matches_file()? {
                bail!("metadata changed while retaining reuse evidence");
            }
        }
        Ok(retained)
    }
}

#[cfg(test)]
mod tests;
