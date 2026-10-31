//! Shared publication contracts for cohort analysis outputs.
pub(crate) mod annotations;
pub(crate) mod corr;
pub(crate) mod custom_matrix;
pub(crate) mod de;
pub(crate) mod geneset;
pub(crate) mod inputs;
pub(crate) mod producer;
mod reuse;
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::generations::OutputGroup;
use rnaa_core::model::{ArtifactKind, ArtifactRecord};
use rnaa_core::util::now_rfc3339;
use std::path::PathBuf;

pub(crate) fn prepare_analysis_artifacts(
    db: &Database,
    outputs: &[PathBuf],
    group: OutputGroup,
) -> Result<Vec<ArtifactRecord>> {
    let mut records = std::collections::BTreeMap::new();
    for path in outputs {
        let name = path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or_default();
        let kind = match group {
            OutputGroup::Refs
            | OutputGroup::Download
            | OutputGroup::Preprocess
            | OutputGroup::Quant => {
                bail!("ingest and references require its own output contract")
            }
            OutputGroup::De => match name {
                "sample_metadata.tsv" => ArtifactKind::Samplesheet,
                "gene_counts.tsv" | "gene_counts_annotated.tsv" => ArtifactKind::Counts,
                "gene_norm_counts.tsv" | "gene_norm_counts_annotated.tsv" => {
                    ArtifactKind::NormalizedCounts
                }
                "size_factors.tsv" => ArtifactKind::SizeFactors,
                "vst.tsv" | "vst_annotated.tsv" => ArtifactKind::Vst,
                "vst.rds" => ArtifactKind::VstRds,
                name if name.ends_with(".json") => ArtifactKind::Manifest,
                name if name.starts_with("de_") && name.ends_with(".tsv") => ArtifactKind::DeTable,
                _ => bail!("unknown DE output contract: {}", path.display()),
            },
            OutputGroup::Corr => {
                if name.ends_with(".json") {
                    ArtifactKind::Manifest
                } else if name.contains("edges") {
                    ArtifactKind::CorrEdges
                } else if name.contains("dense") {
                    ArtifactKind::CorrDense
                } else if name.contains("modules_")
                    || matches!(name, "module_membership.tsv" | "module_eigengenes.tsv")
                {
                    ArtifactKind::CorrModules
                } else if name.contains("module_fgsea") || name.contains("module_t_values") {
                    ArtifactKind::EnrichmentTable
                } else if name.contains("adjusted_matrix") {
                    ArtifactKind::AdjustedMatrix
                } else {
                    bail!("unknown correlation output contract: {}", path.display());
                }
            }
        };
        let snapshot = rnaa_core::ManifestArtifact::snapshot(kind.as_str(), path)?;
        records.insert(
            path.clone(),
            ArtifactRecord {
                project_id: db.project_id().into(),
                run_accession: None,
                kind,
                path: snapshot.path,
                blob_id: None,
                shared_path: None,
                checksum_type: snapshot.checksum_type,
                checksum: snapshot.checksum,
                bytes: snapshot.bytes,
                created_at: now_rfc3339(),
            },
        );
    }
    Ok(records.into_values().collect())
}
