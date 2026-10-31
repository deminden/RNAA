//! Shared registry retirement for fresh publication and historical selection.
use super::{GenerationHead, OutputGroup};
use crate::{ArtifactRecord, Database, ManifestArtifact};
use anyhow::{Context, Result, bail};
use rusqlite::params;
use std::path::Path;

impl Database {
    /// Typed original output witnesses, including intentionally pruned FASTQs.
    /// This validates catalog membership; verify_generation_head checks payloads.
    pub fn generation_output_artifacts(
        &self,
        head: &GenerationHead,
    ) -> Result<Vec<ManifestArtifact>> {
        let seal = self
            .committed_seal(&head.generation_id)?
            .context("output generation is not committed")?;
        if seal.head != *head {
            bail!("output generation does not match its sealed head");
        }
        let json: String = self.read(|conn| {
            Ok(conn.query_row(
                "SELECT artifacts_json FROM artifact_generations
            WHERE project_id=?1 AND generation_id=?2 AND state='COMMITTED'",
                params![self.project_id(), head.generation_id],
                |row| row.get(0),
            )?)
        })?;
        let records: Vec<ArtifactRecord> = serde_json::from_str(&json)?;
        let mut outputs = Vec::new();
        for record in records {
            if record.project_id != self.project_id()
                || !seal.files.iter().any(|file| {
                    file.path == record.path
                        && file.checksum_type == record.checksum_type
                        && file.checksum == record.checksum
                        && file.bytes == record.bytes
                })
            {
                bail!("output registration differs from its immutable inventory");
            }
            outputs.push(ManifestArtifact {
                kind: record.kind.as_str().into(),
                path: record.path,
                checksum_type: record.checksum_type,
                checksum: record.checksum,
                bytes: record.bytes,
            });
        }
        Ok(outputs)
    }

    /// Retire registrations only. Immutable payloads and legacy files stay intact.
    /// The caller holds the publication/selection transaction.
    pub(super) fn retire_generation_artifacts(
        &self,
        group: OutputGroup,
        subject: &str,
    ) -> Result<()> {
        let previous = self.generation_artifacts(group, subject)?;
        if let Some(previous) = previous {
            for record in previous {
                self.delete_artifact(&record.path)?;
            }
        } else {
            // Retire legacy registry entries; their physical files are kept.
            let root = self
                .path()
                .parent()
                .context("project database has no root")?;
            let legacy = root.join(group.as_str());
            let legacy_manifest = root
                .join("manifests/corr")
                .join(format!("{}.json", subject));
            let legacy_records = match group {
                OutputGroup::Download | OutputGroup::Preprocess | OutputGroup::Quant => {
                    self.list_artifacts_for_run(Some(subject))?
                }
                OutputGroup::Refs | OutputGroup::De | OutputGroup::Corr => self.list_artifacts()?,
            };
            for record in legacy_records {
                let retire = match group {
                    OutputGroup::Refs => {
                        record.run_accession.is_none()
                            && (Path::new(&record.path).starts_with(legacy.join(subject))
                                || Path::new(&record.path)
                                    == root
                                        .join("manifests/refs")
                                        .join(format!("{}.json", subject)))
                    }
                    OutputGroup::Download => {
                        record.run_accession.as_deref() == Some(subject)
                            && (matches!(
                                record.kind,
                                crate::ArtifactKind::FastqR1
                                    | crate::ArtifactKind::FastqR2
                                    | crate::ArtifactKind::FastqSingle
                                    | crate::ArtifactKind::Sra
                            ) || record.kind == crate::ArtifactKind::Manifest
                                && Path::new(&record.path)
                                    == root
                                        .join("manifests/download")
                                        .join(format!("{}.json", subject)))
                    }
                    OutputGroup::Preprocess => {
                        record.run_accession.as_deref() == Some(subject)
                            && (matches!(
                                record.kind,
                                crate::ArtifactKind::FastqTrimmedR1
                                    | crate::ArtifactKind::FastqTrimmedR2
                                    | crate::ArtifactKind::FastqTrimmedSingle
                                    | crate::ArtifactKind::PreprocessReport
                            ) || record.kind == crate::ArtifactKind::Manifest
                                && Path::new(&record.path)
                                    .file_name()
                                    .is_some_and(|name| name == "rnaa_preprocess_manifest.json"))
                    }
                    OutputGroup::Quant => {
                        record.run_accession.as_deref() == Some(subject)
                            && matches!(
                                record.kind,
                                crate::ArtifactKind::QuantDir
                                    | crate::ArtifactKind::QuantAbundanceH5
                                    | crate::ArtifactKind::QuantAbundanceTsv
                                    | crate::ArtifactKind::QuantRunInfo
                            )
                            || record.run_accession.as_deref() == Some(subject)
                                && record.kind == crate::ArtifactKind::Manifest
                                && Path::new(&record.path)
                                    .file_name()
                                    .is_some_and(|name| name == "rnaa_quant_manifest.json")
                    }
                    OutputGroup::De | OutputGroup::Corr => {
                        record.run_accession.is_none()
                            && (Path::new(&record.path).starts_with(&legacy)
                                || (group == OutputGroup::Corr
                                    && Path::new(&record.path) == legacy_manifest))
                    }
                };
                if retire {
                    self.delete_artifact(&record.path)?;
                }
            }
        }
        Ok(())
    }
}
