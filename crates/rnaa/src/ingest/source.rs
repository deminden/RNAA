//! Current read identity, including exact evidence for intentionally pruned files.
use crate::ingest::inputs::load_raw_fastq_artifacts;
use anyhow::{Result, bail};
use rnaa_core::{
    Database, ManifestArtifact, ProjectConfig, RunRecord,
    fingerprint::InputIdentity,
    generations::{GenerationHead, OutputGroup},
};
use std::path::Path;
pub(super) fn raw(input: &ManifestArtifact) -> bool {
    matches!(
        input.kind.as_str(),
        "FASTQ_SINGLE" | "FASTQ_R1" | "FASTQ_R2"
    )
}
pub(super) fn trimmed(input: &ManifestArtifact) -> bool {
    matches!(
        input.kind.as_str(),
        "FASTQ_TRIMMED_SINGLE" | "FASTQ_TRIMMED_R1" | "FASTQ_TRIMMED_R2"
    )
}
pub(super) fn identities(inputs: &[ManifestArtifact]) -> Vec<InputIdentity> {
    let mut values = inputs.iter().map(InputIdentity::from).collect::<Vec<_>>();
    values.sort();
    values
}
pub(super) fn ordered(mut inputs: Vec<ManifestArtifact>) -> Vec<ManifestArtifact> {
    inputs.sort_by(|a, b| (&a.kind, &a.path).cmp(&(&b.kind, &b.path)));
    inputs
}
fn registry(db: &Database, run: &str) -> Result<Vec<ManifestArtifact>> {
    Ok(ordered(
        load_raw_fastq_artifacts(db, run)?
            .into_iter()
            .map(|file| ManifestArtifact {
                kind: file.artifact_kind.as_str().into(),
                path: file.path.display().to_string(),
                checksum_type: file.checksum_type,
                checksum: file.checksum,
                bytes: file.bytes,
            })
            .collect(),
    ))
}
pub(super) struct Source {
    pub raw: Vec<ManifestArtifact>,
    download: Option<GenerationHead>,
    registered: Vec<ManifestArtifact>,
}
impl Source {
    pub fn capture(db: &Database, run: &RunRecord, config: &ProjectConfig) -> Result<Option<Self>> {
        if !rnaa_stages::quant::source_binding_current(db, run, config)? {
            return Ok(None);
        }
        let download = db.generation_head(OutputGroup::Download, &run.run_accession)?;
        let registered = registry(db, &run.run_accession)?;
        let inputs = if let Some(head) = &download {
            if !db.verify_generation_head(head)? {
                return Ok(None);
            }
            let inputs = ordered(
                db.generation_output_artifacts(head)?
                    .into_iter()
                    .filter(raw)
                    .collect(),
            );
            if registered.iter().any(|current| {
                !inputs
                    .iter()
                    .any(|source| source.path == current.path && source.kind == current.kind)
            }) {
                return Ok(None);
            }
            inputs
        } else if !registered.is_empty() {
            if registered
                .iter()
                .any(|file| !Path::new(&file.path).is_file())
            {
                return Ok(None);
            }
            let mut inputs = Vec::new();
            for item in &registered {
                inputs.push(ManifestArtifact::snapshot(
                    &item.kind,
                    Path::new(&item.path),
                )?);
            }
            ordered(inputs)
        } else {
            let mut inputs = ordered(
                db.generation_inputs(OutputGroup::Quant, &run.run_accession)?
                    .unwrap_or_default()
                    .into_iter()
                    .filter(raw)
                    .collect(),
            );
            if inputs.is_empty() {
                inputs = ordered(
                    db.generation_inputs(OutputGroup::Preprocess, &run.run_accession)?
                        .unwrap_or_default()
                        .into_iter()
                        .filter(raw)
                        .collect(),
                );
            }
            let retained = db.verified_retained_inputs(&run.run_accession)?;
            if inputs.iter().any(|input| !retained.contains(input)) {
                return Ok(None);
            }
            inputs
        };
        let valid = match run.library_layout {
            rnaa_core::LibraryLayout::Single => {
                !inputs.is_empty() && inputs.iter().all(|file| file.kind == "FASTQ_SINGLE")
            }
            rnaa_core::LibraryLayout::Paired => {
                !inputs.is_empty()
                    && inputs
                        .iter()
                        .all(|file| matches!(file.kind.as_str(), "FASTQ_R1" | "FASTQ_R2"))
                    && inputs.iter().filter(|file| file.kind == "FASTQ_R1").count()
                        == inputs.iter().filter(|file| file.kind == "FASTQ_R2").count()
            }
            rnaa_core::LibraryLayout::Unknown => false,
        };
        Ok(valid.then_some(Self {
            raw: inputs,
            download,
            registered,
        }))
    }
    pub fn verify_catalog(
        &self,
        db: &Database,
        run: &RunRecord,
        config: &ProjectConfig,
    ) -> Result<()> {
        rnaa_stages::quant::verify_quant_catalog(db, run, config)?;
        if db.generation_head(OutputGroup::Download, &run.run_accession)? != self.download
            || registry(db, &run.run_accession)? != self.registered
        {
            bail!("current read source changed during historical selection");
        }
        Ok(())
    }
    pub fn verify_current(
        &self,
        db: &Database,
        run: &RunRecord,
        config: &ProjectConfig,
    ) -> Result<()> {
        rnaa_stages::quant::verify_quant_current(db, run, config)?;
        self.verify_catalog(db, run, config)?;
        if let Some(head) = &self.download {
            if !db.verify_generation_head(head)? {
                bail!("current read producer failed verification");
            }
        } else {
            let retained = db.verified_retained_inputs(&run.run_accession)?;
            for input in &self.raw {
                if Path::new(&input.path).exists() {
                    if !input.matches_file()? {
                        bail!("current local reads changed during historical selection");
                    }
                } else if !retained.contains(input) {
                    bail!("current local reads disappeared without retention evidence");
                }
            }
        }
        Ok(())
    }
}
