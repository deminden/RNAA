//! Cohort-stage requests captured before fitting and verified before publication.
use crate::corr::ModuleGseaParams;
use anyhow::{Result, bail};
use rnaa_core::config::ProjectConfig;
use rnaa_core::fingerprint::{InputIdentity, StageRequest};
use rnaa_core::util::write_json_pretty;
use rnaa_core::{ContrastSpec, ManifestArtifact, ProjectPaths, ReferenceBundle, StageManifest};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
mod metadata;
mod policy;
use policy::ObservedConfig;

#[derive(Debug, Clone)]
pub struct AnalysisRequest<P> {
    pub request: StageRequest<P>,
    pub inputs: Vec<ManifestArtifact>,
    policy: Option<ObservedConfig>,
    observations: Vec<ManifestArtifact>,
}
impl<P> AnalysisRequest<P> {
    /// Mutable source observations fence this attempt but are not retained inputs.
    pub fn observe_input(&mut self, input: ManifestArtifact) {
        self.observations.push(input);
    }
}
impl<P> AnalysisRequest<P> {
    pub fn relocate_input(&mut self, role: &str, path: &Path) -> Result<()> {
        let matches = self
            .inputs
            .iter()
            .enumerate()
            .filter(|(_, input)| input.kind == role)
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        if matches.len() != 1 {
            bail!("relocation requires one unambiguous input role '{role}'");
        }
        let replacement = ManifestArtifact::snapshot(role, path)?;
        let index = matches[0];
        if InputIdentity::from(&replacement) != InputIdentity::from(&self.inputs[index]) {
            bail!("frozen input changed while relocating '{role}'");
        }
        self.inputs[index] = replacement;
        Ok(())
    }
    pub fn verify_inputs(&self) -> Result<()> {
        for input in self.inputs.iter().chain(&self.observations) {
            if !input.matches_file()? {
                bail!("analysis input changed during computation: {}", input.path);
            }
        }
        Ok(())
    }
}
impl<P: Serialize> AnalysisRequest<P> {
    pub fn seal(&self, path: &Path, outputs: &[PathBuf]) -> Result<()> {
        self.verify_inputs()?;
        let mut manifest: StageManifest = serde_json::from_slice(&fs::read(path)?)?;
        manifest.input_artifacts = self.inputs.clone();
        manifest.output_artifacts = outputs
            .iter()
            .filter(|output| output.as_path() != path)
            .map(|output| ManifestArtifact::snapshot("ANALYSIS_OUTPUT", output))
            .collect::<Result<_>>()?;
        manifest.parameters["request_sha256"] = serde_json::to_value(self.request.digest()?)?;
        manifest.parameters["request"] = serde_json::to_value(&self.request)?;
        manifest.parameters["configuration_observation"] = self.configuration_observation()?;
        write_json_pretty(path, &manifest)
    }
    pub fn write_descriptor(&self, path: &Path) -> Result<()> {
        write_json_pretty(
            path,
            &serde_json::json!({"schema_version":1,"request_sha256":self.request.digest()?,"request": &self.request,"input_witnesses":&self.inputs,"configuration_observation":self.configuration_observation()?}),
        )
    }
    pub fn matches_descriptor(&self, path: &Path) -> Result<bool> {
        let Some(value) = fs::read(path)
            .ok()
            .and_then(|bytes| serde_json::from_slice::<serde_json::Value>(&bytes).ok())
        else {
            return Ok(false);
        };
        Ok(value["schema_version"] == 1
            && value["request_sha256"] == self.request.digest()?
            && value["request"] == serde_json::to_value(&self.request)?)
    }

    pub fn matches(
        &self,
        path: &Path,
        stage: &str,
        id: &str,
        required: &[PathBuf],
    ) -> Result<bool> {
        let Some(manifest) = fs::read(path)
            .ok()
            .and_then(|bytes| serde_json::from_slice::<StageManifest>(&bytes).ok())
        else {
            return Ok(false);
        };
        if !required
            .iter()
            .filter(|output| output.as_path() != path)
            .all(|output| {
                manifest
                    .output_artifacts
                    .iter()
                    .any(|artifact| Path::new(&artifact.path) == output)
            })
        {
            return Ok(false);
        }
        Ok(manifest.stage == stage
            && manifest.id == id
            && manifest.parameters["request_sha256"] == self.request.digest()?
            && manifest.parameters["request"] == serde_json::to_value(&self.request)?
            // Physical input paths may belong to a previous immutable generation.
            // The typed request already binds exact roles, names, hashes and sizes.
            && self.inputs.iter().map(ManifestArtifact::matches_file).collect::<Result<Vec<_>>>()?.into_iter().all(|valid| valid)
            && manifest.verify_local_artifacts()?)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeOptions {
    design: String,
    transform: String,
    counts_from_abundance: String,
    contrasts: Vec<ContrastSpec>,
    samples: Vec<String>,
    sample_plan: rnaa_core::sample_plan::SampleMembershipPlan,
    design_contract: crate::design::ScientificDesign,
}
#[allow(clippy::too_many_arguments)]
pub fn deseq_request(
    paths: &ProjectPaths,
    sheet: &Path,
    reference: &ReferenceBundle,
    design: &str,
    contrasts: &[ContrastSpec],
    config: &ProjectConfig,
    stage: &str,
) -> Result<AnalysisRequest<DeOptions>> {
    config.validate()?;
    let mut inputs = vec![
        ManifestArtifact::snapshot("SAMPLESHEET", sheet)?,
        ManifestArtifact::snapshot("TX2GENE", &reference.tx2gene_path)?,
    ];
    if let Some(annotation) = &reference.gene_annotation_path {
        inputs.push(ManifestArtifact::snapshot("ANNOTATION", annotation)?);
    }
    let (sample_plan, prepared) = metadata::prepare_de(sheet, design, contrasts, config)?;
    let samples = sample_plan.sample_ids();
    for run in sample_plan.samples.iter().flat_map(|sample| &sample.runs) {
        inputs.push(ManifestArtifact::snapshot(
            &format!("ABUNDANCE:{run}"),
            &paths
                .quant_run_dir(run, &config.quant.engine)
                .join("abundance.tsv"),
        )?);
        let dir = paths.quant_run_dir(run, &config.quant.engine);
        let provenance = if dir.join("rnaa_quant_import.json").exists() {
            dir.join("rnaa_quant_import.json")
        } else {
            dir.join("rnaa_quant_manifest.json")
        };
        inputs.push(ManifestArtifact::snapshot(
            &format!("QUANT_PROVENANCE:{run}"),
            &provenance,
        )?);
    }
    let request = StageRequest {
        schema_version: 1,
        stage: stage.into(),
        algorithm_version: crate::deseq2::DE_ALGORITHM_VERSION,
        engines: BTreeMap::from([
            ("rsdeseq2".into(), crate::engines::RSDESEQ2.into()),
            (
                "dependency_lock".into(),
                crate::engines::DEPENDENCY_LOCK_SHA256.into(),
            ),
        ]),
        parameters: DeOptions {
            design: design.into(),
            transform: config.deseq2.transform.clone(),
            counts_from_abundance: config.deseq2.counts_from_abundance.clone(),
            contrasts: contrasts.to_vec(),
            samples,
            sample_plan: sample_plan.membership(),
            design_contract: prepared,
        },
        inputs: metadata::scientific_inputs(&inputs),
        dependencies: BTreeMap::new(),
    };
    Ok(AnalysisRequest {
        request,
        inputs,
        policy: Some(ObservedConfig::de(config, design)),
        observations: Vec::new(),
    })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorrOptions {
    adjust: String,
    model: Option<String>,
    method: String,
    geneset: String,
    output: String,
    modules: Option<ModuleOptions>,
    design_contract: Option<crate::design::ScientificDesign>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
struct ModuleOptions {
    #[serde(default)]
    algorithm_version: u32,
    top: usize,
    minimum: usize,
    permutations: usize,
    seed: u64,
    #[serde(default)]
    pathway_minimum: usize,
    pathway_maximum: Option<usize>,
}
pub fn correlation_request(
    matrix: &Path,
    sheet: &Path,
    annotation: Option<&Path>,
    config: &ProjectConfig,
    modules: Option<&ModuleGseaParams>,
) -> Result<AnalysisRequest<CorrOptions>> {
    correlation_request_inner(matrix, sheet, annotation, config, modules, None)
}

/// Reconstruct a request without reading beyond the matrix size admitted by the caller.
pub fn correlation_request_exact(
    matrix: &Path,
    sheet: &Path,
    annotation: Option<&Path>,
    config: &ProjectConfig,
    modules: Option<&ModuleGseaParams>,
    matrix_bytes: u64,
) -> Result<AnalysisRequest<CorrOptions>> {
    correlation_request_inner(
        matrix,
        sheet,
        annotation,
        config,
        modules,
        Some(matrix_bytes),
    )
}

fn correlation_request_inner(
    matrix: &Path,
    sheet: &Path,
    annotation: Option<&Path>,
    config: &ProjectConfig,
    modules: Option<&ModuleGseaParams>,
    matrix_bytes: Option<u64>,
) -> Result<AnalysisRequest<CorrOptions>> {
    config.validate()?;
    let design_contract = if config.corr.adjust == "residualize" {
        let prepared = if let Some(bytes) = matrix_bytes {
            crate::corr::prepare_residualization_exact(
                matrix,
                bytes,
                sheet,
                &config.corr.model,
                &config.metadata,
            )?
        } else {
            crate::corr::prepare_residualization(
                matrix,
                sheet,
                &config.corr.model,
                &config.metadata,
            )?
        };
        Some(prepared.contract.scientific())
    } else if config.corr.adjust == "none" {
        None
    } else {
        bail!(
            "unsupported matrix adjustment '{}'; expected none or residualize",
            config.corr.adjust
        );
    };
    let matrix_witness = if let Some(bytes) = matrix_bytes {
        ManifestArtifact::snapshot_exact("MATRIX", matrix, bytes)?
    } else {
        ManifestArtifact::snapshot("MATRIX", matrix)?
    };
    let mut inputs = vec![
        matrix_witness,
        ManifestArtifact::snapshot("METADATA", sheet)?,
    ];
    if let Some(annotation) = annotation {
        inputs.push(ManifestArtifact::snapshot("ANNOTATION", annotation)?);
    }
    let geneset = if let Some(path) = config.corr.geneset.strip_prefix("genes:") {
        inputs.push(ManifestArtifact::snapshot(
            "GENE_SELECTION",
            Path::new(path),
        )?);
        "genes:file".to_string()
    } else {
        config.corr.geneset.clone()
    };
    if let Some(module) = modules {
        module.validate()?;
        inputs.push(ManifestArtifact::snapshot("PATHWAYS", &module.gmt_path)?);
    }
    let request = StageRequest {
        schema_version: 1,
        stage: "corr".into(),
        algorithm_version: crate::corr::CORRELATION_ALGORITHM_VERSION,
        engines: BTreeMap::from([
            ("mincorr".into(), crate::engines::MINCORR_REVISION.into()),
            ("rsfgsea".into(), crate::engines::RSFGSEA.into()),
            ("rsdeseq2_design".into(), crate::engines::RSDESEQ2.into()),
            (
                "dependency_lock".into(),
                crate::engines::DEPENDENCY_LOCK_SHA256.into(),
            ),
        ]),
        parameters: CorrOptions {
            adjust: config.corr.adjust.clone(),
            model: (config.corr.adjust != "none").then(|| config.corr.model.clone()),
            method: config.corr.method.clone(),
            geneset,
            output: config.corr.output.clone(),
            modules: modules.map(|p| ModuleOptions {
                algorithm_version: crate::corr::MODULE_ALGORITHM_VERSION,
                top: p.top_genes,
                minimum: p.min_module_size,
                permutations: p.permutations,
                seed: p.seed,
                pathway_minimum: p.pathway_min_size,
                pathway_maximum: p.pathway_max_size,
            }),
            design_contract,
        },
        inputs: metadata::scientific_inputs(&inputs),
        dependencies: BTreeMap::new(),
    };
    Ok(AnalysisRequest {
        request,
        inputs,
        policy: Some(ObservedConfig::corr(config)),
        observations: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_matrix_requests_match_normal_requests_and_reject_size_changes() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let matrix = temp.path().join("matrix.tsv");
        let sheet = temp.path().join("samples.tsv");
        fs::write(&matrix, b"gene_id\tS1\tS2\tS3\tS4\nG1\t1\t2\t3\t4\n")?;
        fs::write(
            &sheet,
            b"sample_id\tcondition\nS1\tA\nS2\tA\nS3\tB\nS4\tB\n",
        )?;
        let mut config = ProjectConfig::default();
        config.corr.model = "~ condition".into();
        let admitted = fs::metadata(&matrix)?.len();
        for adjust in ["residualize", "none"] {
            config.corr.adjust = adjust.into();
            let ordinary = correlation_request(&matrix, &sheet, None, &config, None)?;
            let exact = correlation_request_exact(&matrix, &sheet, None, &config, None, admitted)?;
            assert_eq!(ordinary.request.digest()?, exact.request.digest()?);
            assert!(
                correlation_request_exact(&matrix, &sheet, None, &config, None, admitted - 1)
                    .is_err()
            );
            let original = fs::read(&matrix)?;
            fs::write(&matrix, [original.as_slice(), b"extra"].concat())?;
            assert!(
                correlation_request_exact(&matrix, &sheet, None, &config, None, admitted).is_err()
            );
            fs::write(&matrix, original)?;
        }
        Ok(())
    }
    #[test]
    fn manifests_require_complete_output_coverage_and_unchanged_inputs() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let source = temp.path().join("source");
        let output = temp.path().join("output");
        let manifest_path = temp.path().join("manifest.json");
        let descriptor = temp.path().join("request.json");
        fs::write(&source, b"source")?;
        fs::write(&output, b"result")?;
        let inputs = vec![ManifestArtifact::snapshot("MATRIX", &source)?];
        let plan = AnalysisRequest {
            request: StageRequest {
                schema_version: 1,
                stage: "corr".into(),
                algorithm_version: 2,
                engines: BTreeMap::new(),
                parameters: true,
                inputs: inputs.iter().map(InputIdentity::from).collect(),
                dependencies: BTreeMap::new(),
            },
            inputs,
            policy: None,
            observations: Vec::new(),
        };
        write_json_pretty(
            &manifest_path,
            &StageManifest {
                schema_version: 1,
                stage: "corr".into(),
                id: "project".into(),
                started_at: String::new(),
                finished_at: String::new(),
                parameters: serde_json::json!({}),
                tool_versions: BTreeMap::new(),
                input_artifacts: vec![],
                output_artifacts: vec![],
                notes: vec![],
            },
        )?;
        let required = vec![output.clone(), manifest_path.clone()];
        plan.seal(&manifest_path, &required)?;
        plan.write_descriptor(&descriptor)?;
        assert!(plan.matches(&manifest_path, "corr", "project", &required)?);
        assert!(plan.matches_descriptor(&descriptor)?);
        let missing = temp.path().join("unlisted");
        fs::write(&missing, b"unlisted")?;
        let mut expanded = required.clone();
        expanded.push(missing);
        assert!(!plan.matches(&manifest_path, "corr", "project", &expanded)?);
        fs::write(&descriptor, b"{partial")?;
        assert!(!plan.matches_descriptor(&descriptor)?);
        fs::write(&source, b"edited")?;
        assert!(!plan.matches(&manifest_path, "corr", "project", &required)?);
        assert!(plan.seal(&manifest_path, &required).is_err());
        Ok(())
    }
}
