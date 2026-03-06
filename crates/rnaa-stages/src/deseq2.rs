use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, Result, anyhow, bail};
use rnaa_core::config::ProjectConfig;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{ContrastSpec, ReferenceBundle};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::DifferentialExpression;
use rnaa_core::util::{command_version, file_size, now_rfc3339, write_json_pretty};
use serde_json::json;

#[derive(Debug, Clone, Default)]
pub struct Deseq2Runner;

impl DifferentialExpression for Deseq2Runner {
    fn deseq2(
        &self,
        project_id: &str,
        reference: &ReferenceBundle,
        design: &str,
        contrasts: &[ContrastSpec],
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<Vec<PathBuf>> {
        let script_path = find_deseq2_script(paths)?;
        let samplesheet_path = paths.samplesheet_path();
        if !samplesheet_path.exists() {
            bail!(
                "samplesheet missing at {}; run `rnaa resolve` first",
                samplesheet_path.display()
            );
        }
        let de_dir = paths.de_project_dir(project_id);
        fs::create_dir_all(&de_dir)
            .with_context(|| format!("failed to create {}", de_dir.display()))?;
        let manifest_path = de_dir.join("deseq2_manifest.json");

        let mut cmd = Command::new("Rscript");
        cmd.arg("--vanilla")
            .arg(&script_path)
            .arg("--samplesheet")
            .arg(&samplesheet_path)
            .arg("--quant-root")
            .arg(&paths.quant_dir)
            .arg("--engine")
            .arg(&config.quant.engine)
            .arg("--tx2gene")
            .arg(&reference.tx2gene_path)
            .arg("--gene-annotation")
            .arg(
                reference
                    .gene_annotation_path
                    .as_ref()
                    .ok_or_else(|| anyhow!("reference missing gene annotation table"))?,
            )
            .arg("--design")
            .arg(design)
            .arg("--transform")
            .arg(&config.deseq2.transform)
            .arg("--counts-from-abundance")
            .arg(&config.deseq2.counts_from_abundance)
            .arg("--ignore-tx-version")
            .arg("true")
            .arg("--ignore-after-bar")
            .arg("false")
            .arg("--outdir")
            .arg(&de_dir)
            .arg("--manifest-out")
            .arg(&manifest_path);
        for contrast in contrasts {
            let value = format!(
                "{}|{}|{}|{}",
                contrast.factor, contrast.level_a, contrast.level_b, contrast.name
            );
            cmd.arg("--contrast").arg(value);
        }

        let status = cmd
            .status()
            .context("failed to spawn Rscript for DESeq2 stage")?;
        if !status.success() {
            bail!("DESeq2 R script failed with status {status}");
        }

        let mut outputs = vec![
            de_dir.join("gene_counts.tsv"),
            de_dir.join("gene_counts_annotated.tsv"),
            de_dir.join("gene_norm_counts.tsv"),
            de_dir.join("gene_norm_counts_annotated.tsv"),
            de_dir.join("vst.tsv"),
            de_dir.join("vst_annotated.tsv"),
            manifest_path.clone(),
        ];
        let optional = [de_dir.join("vst.rds"), de_dir.join("size_factors.tsv")];
        for path in optional {
            if path.exists() {
                outputs.push(path);
            }
        }
        for contrast in contrasts {
            let path = de_dir.join(format!("de_{}.tsv", contrast.name));
            if path.exists() {
                outputs.push(path);
            }
            let annotated = de_dir.join(format!("de_{}_annotated.tsv", contrast.name));
            if annotated.exists() {
                outputs.push(annotated);
            }
        }

        for required in ["gene_counts.tsv", "gene_norm_counts.tsv", "vst.tsv"] {
            let path = de_dir.join(required);
            if !path.exists() || file_size(&path).unwrap_or_default() == 0 {
                bail!("DESeq2 output missing or empty {}", path.display());
            }
        }

        if !manifest_path.exists() {
            let manifest = StageManifest {
                schema_version: 1,
                stage: "deseq2".to_string(),
                id: project_id.to_string(),
                started_at: now_rfc3339(),
                finished_at: now_rfc3339(),
                parameters: json!({
                    "design": design,
                    "transform": config.deseq2.transform,
                    "counts_from_abundance": config.deseq2.counts_from_abundance,
                    "engine": config.quant.engine
                }),
                tool_versions: BTreeMap::from_iter(
                    [(
                        "Rscript".to_string(),
                        command_version("Rscript", &["--version"]),
                    )]
                    .into_iter()
                    .filter_map(|(name, value)| value.map(|value| (name, value))),
                ),
                input_artifacts: vec![ManifestArtifact {
                    kind: "SAMPLESHEET".to_string(),
                    path: samplesheet_path.display().to_string(),
                    checksum_type: "none".to_string(),
                    checksum: String::new(),
                    bytes: file_size(&samplesheet_path).unwrap_or_default(),
                }],
                output_artifacts: outputs
                    .iter()
                    .map(|path| ManifestArtifact {
                        kind: "DE_OUTPUT".to_string(),
                        path: path.display().to_string(),
                        checksum_type: "none".to_string(),
                        checksum: String::new(),
                        bytes: file_size(path).unwrap_or_default(),
                    })
                    .collect(),
                notes: vec!["Fallback DESeq2 manifest written by Rust wrapper.".to_string()],
            };
            write_json_pretty(&manifest_path, &manifest)?;
            if !outputs.iter().any(|path| path == &manifest_path) {
                outputs.push(manifest_path);
            }
        }

        Ok(outputs)
    }
}

fn find_deseq2_script(paths: &ProjectPaths) -> Result<PathBuf> {
    let mut candidates = Vec::new();
    candidates.push(paths.root.join("r").join("tximport_deseq2.R"));
    if let Ok(dir) = std::env::var("RNAA_R_SCRIPTS_DIR") {
        candidates.push(PathBuf::from(dir).join("tximport_deseq2.R"));
    }
    candidates.push(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("r")
            .join("tximport_deseq2.R"),
    );
    for candidate in candidates {
        if candidate.exists() {
            return Ok(candidate);
        }
    }
    bail!(
        "DESeq2 script not found; checked project r/tximport_deseq2.R, RNAA_R_SCRIPTS_DIR, and bundled source path"
    )
}
