use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result, bail};
use rnaa_core::config::ProjectConfig;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{QuantArtifacts, ReferenceBundle, RunRecord, VerifiedFile};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::Quantifier;
use rnaa_core::util::{command_version, file_size, now_rfc3339, write_json_pretty};
use serde_json::json;

#[derive(Debug, Clone, Default)]
pub struct RKallistoQuantifier;

impl Quantifier for RKallistoQuantifier {
    fn quantify(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        reference: &ReferenceBundle,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<QuantArtifacts> {
        let read1 = fastqs
            .iter()
            .find(|item| {
                item.artifact_kind == rnaa_core::ArtifactKind::FastqR1
                    || item.artifact_kind == rnaa_core::ArtifactKind::FastqSingle
                    || item.artifact_kind == rnaa_core::ArtifactKind::FastqTrimmedR1
                    || item.artifact_kind == rnaa_core::ArtifactKind::FastqTrimmedSingle
            })
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "no read1 FASTQ artifact found for run {}",
                    run.run_accession
                )
            })?;
        let read2 = fastqs.iter().find(|item| {
            item.artifact_kind == rnaa_core::ArtifactKind::FastqR2
                || item.artifact_kind == rnaa_core::ArtifactKind::FastqTrimmedR2
        });

        let out_dir = paths.quant_run_dir(&run.run_accession, &config.quant.engine);
        fs::create_dir_all(&out_dir)
            .with_context(|| format!("failed to create {}", out_dir.display()))?;
        let manifest_path = out_dir.join("rnaa_quant_manifest.json");
        let started_at = now_rfc3339();

        let script_path = find_quant_script(paths)?;
        if config.quant.engine != "r-kallisto" {
            bail!(
                "unsupported quant engine '{}'; currently only 'r-kallisto' is implemented",
                config.quant.engine
            );
        }

        let mut cmd = Command::new("Rscript");
        cmd.arg("--vanilla")
            .arg(&script_path)
            .arg("--run")
            .arg(&run.run_accession)
            .arg("--fastq1")
            .arg(&read1.path)
            .arg("--index")
            .arg(&reference.kallisto_index_path)
            .arg("--threads")
            .arg(config.quant.threads.to_string())
            .arg("--outdir")
            .arg(&out_dir)
            .arg("--reference-manifest")
            .arg(&reference.manifest_path)
            .arg("--manifest-out")
            .arg(&manifest_path);
        if let Some(read2) = read2 {
            cmd.arg("--fastq2").arg(&read2.path);
        }
        let status = cmd
            .status()
            .with_context(|| format!("failed to spawn Rscript for run {}", run.run_accession))?;
        if !status.success() {
            bail!(
                "Rscript quantification failed for {} with status {}",
                run.run_accession,
                status
            );
        }

        let abundance_h5 = out_dir.join("abundance.h5");
        let run_info_json = out_dir.join("run_info.json");
        if !abundance_h5.exists() || file_size(&abundance_h5).unwrap_or_default() == 0 {
            bail!(
                "quant output missing or empty abundance.h5 for {}",
                run.run_accession
            );
        }
        if !run_info_json.exists() || file_size(&run_info_json).unwrap_or_default() == 0 {
            bail!(
                "quant output missing or empty run_info.json for {}",
                run.run_accession
            );
        }
        let abundance_tsv = out_dir.join("abundance.tsv");
        if abundance_tsv.exists() {
            fs::remove_file(&abundance_tsv).with_context(|| {
                format!(
                    "failed to remove non-canonical abundance.tsv for {}",
                    run.run_accession
                )
            })?;
        }

        // Fallback manifest writer if the R script did not generate it.
        if !manifest_path.exists() {
            let manifest = StageManifest {
                schema_version: 1,
                stage: "quant".to_string(),
                id: run.run_accession.clone(),
                started_at,
                finished_at: now_rfc3339(),
                parameters: json!({
                    "engine": config.quant.engine,
                    "threads": config.quant.threads,
                    "index": reference.kallisto_index_path.display().to_string(),
                    "outdir": out_dir.display().to_string()
                }),
                tool_versions: BTreeMap::from_iter(
                    [
                        (
                            "Rscript".to_string(),
                            command_version("Rscript", &["--version"]),
                        ),
                        (
                            "kallisto".to_string(),
                            command_version("kallisto", &["version"]),
                        ),
                    ]
                    .into_iter()
                    .filter_map(|(name, version)| version.map(|version| (name, version))),
                ),
                input_artifacts: fastqs
                    .iter()
                    .map(|item| ManifestArtifact {
                        kind: item.artifact_kind.as_str().to_string(),
                        path: item.path.display().to_string(),
                        checksum_type: item.checksum_type.clone(),
                        checksum: item.checksum.clone(),
                        bytes: item.bytes,
                    })
                    .collect(),
                output_artifacts: vec![
                    ManifestArtifact {
                        kind: "KALLISTO_DIR".to_string(),
                        path: out_dir.display().to_string(),
                        checksum_type: "none".to_string(),
                        checksum: String::new(),
                        bytes: 0,
                    },
                    ManifestArtifact {
                        kind: "ABUNDANCE_H5".to_string(),
                        path: abundance_h5.display().to_string(),
                        checksum_type: "none".to_string(),
                        checksum: String::new(),
                        bytes: file_size(&abundance_h5).unwrap_or_default(),
                    },
                    ManifestArtifact {
                        kind: "RUN_INFO_JSON".to_string(),
                        path: run_info_json.display().to_string(),
                        checksum_type: "none".to_string(),
                        checksum: String::new(),
                        bytes: file_size(&run_info_json).unwrap_or_default(),
                    },
                ],
                notes: vec![
                    "Fallback manifest emitted by Rust wrapper because R manifest was absent."
                        .to_string(),
                ],
            };
            write_json_pretty(&manifest_path, &manifest)?;
        }

        Ok(QuantArtifacts {
            run_accession: run.run_accession.clone(),
            out_dir,
            abundance_h5,
            run_info_json,
            manifest_path,
        })
    }
}

pub fn reconcile_quant_artifacts(
    run: &RunRecord,
    reference: &ReferenceBundle,
    paths: &ProjectPaths,
    config: &ProjectConfig,
) -> Result<Option<QuantArtifacts>> {
    let out_dir = paths.quant_run_dir(&run.run_accession, &config.quant.engine);
    let abundance_h5 = out_dir.join("abundance.h5");
    let run_info_json = out_dir.join("run_info.json");
    let manifest_path = out_dir.join("rnaa_quant_manifest.json");

    for path in [&abundance_h5, &run_info_json, &manifest_path] {
        if !path.exists() || file_size(path).unwrap_or_default() == 0 {
            return Ok(None);
        }
    }

    if !quant_manifest_matches_reference(&manifest_path, &reference.manifest_path)? {
        return Ok(None);
    }

    Ok(Some(QuantArtifacts {
        run_accession: run.run_accession.clone(),
        out_dir,
        abundance_h5,
        run_info_json,
        manifest_path,
    }))
}

fn quant_manifest_matches_reference(
    manifest_path: &Path,
    expected_reference: &Path,
) -> Result<bool> {
    let raw = fs::read_to_string(manifest_path)
        .with_context(|| format!("failed to read quant manifest {}", manifest_path.display()))?;
    let value: serde_json::Value = serde_json::from_str(&raw)
        .with_context(|| format!("failed to parse quant manifest {}", manifest_path.display()))?;
    let manifest_reference = value
        .get("parameters")
        .and_then(|params| params.get("reference_manifest"))
        .and_then(|value| value.as_str())
        .map(PathBuf::from);
    Ok(manifest_reference.as_deref() == Some(expected_reference))
}

#[allow(dead_code)]
fn _exists_nonempty(path: &Path) -> bool {
    path.exists() && file_size(path).unwrap_or_default() > 0
}

fn find_quant_script(paths: &ProjectPaths) -> Result<std::path::PathBuf> {
    let mut candidates = Vec::new();
    candidates.push(paths.root.join("r").join("quant_kallisto.R"));
    if let Ok(dir) = std::env::var("RNAA_R_SCRIPTS_DIR") {
        candidates.push(std::path::PathBuf::from(dir).join("quant_kallisto.R"));
    }
    candidates.push(
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("r")
            .join("quant_kallisto.R"),
    );
    for candidate in candidates {
        if candidate.exists() {
            return Ok(candidate);
        }
    }
    bail!(
        "quant script not found; checked project r/quant_kallisto.R, RNAA_R_SCRIPTS_DIR, and bundled source path"
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use rnaa_core::config::ProjectConfig;
    use rnaa_core::model::{LibraryLayout, RunRecord};
    use rnaa_core::paths::ProjectPaths;
    use rnaa_core::state::RunState;
    use rnaa_core::util::{now_rfc3339, write_json_pretty};
    use serde_json::json;
    use std::collections::BTreeMap;

    #[test]
    fn reconcile_quant_artifacts_accepts_matching_reference_manifest() {
        let temp = tempfile::tempdir().unwrap();
        let paths = ProjectPaths::new(temp.path());
        paths.ensure_layout().unwrap();
        let config = ProjectConfig::default();
        let run = RunRecord {
            project_id: "p".to_string(),
            study_accession: "SRP1".to_string(),
            source_study_accession: None,
            run_accession: "SRR1".to_string(),
            sample_accession: None,
            experiment_accession: None,
            library_layout: LibraryLayout::Single,
            instrument_platform: None,
            instrument_model: None,
            remote_files: Vec::new(),
            metadata: BTreeMap::new(),
            state: RunState::Verified,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        let reference = ReferenceBundle {
            id: "human".to_string(),
            organism: "human".to_string(),
            ensembl_release: "115".to_string(),
            cdna_path: temp.path().join("cdna.fa.gz"),
            gtf_path: temp.path().join("genes.gtf.gz"),
            kallisto_index_path: temp.path().join("kallisto.index"),
            tx2gene_path: temp.path().join("tx2gene.tsv"),
            gene_annotation_path: None,
            manifest_path: temp.path().join("ref_manifest.json"),
        };

        fs::write(&reference.manifest_path, "{}").unwrap();
        let out_dir = paths.quant_run_dir(&run.run_accession, &config.quant.engine);
        fs::create_dir_all(&out_dir).unwrap();
        fs::write(out_dir.join("abundance.h5"), b"abc").unwrap();
        fs::write(out_dir.join("run_info.json"), b"{}").unwrap();
        write_json_pretty(
            &out_dir.join("rnaa_quant_manifest.json"),
            &json!({
                "parameters": {
                    "reference_manifest": reference.manifest_path.display().to_string()
                }
            }),
        )
        .unwrap();

        let artifacts = reconcile_quant_artifacts(&run, &reference, &paths, &config).unwrap();
        assert!(artifacts.is_some());
    }
}
