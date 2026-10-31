//! Payload/scientific validity is distinct from current physical input admission.
use super::{quant_request, request_matches};
use anyhow::Result;
use rnaa_core::{
    ManifestArtifact, ProjectConfig, ProjectPaths, QuantArtifacts, ReferenceBundle, RunRecord,
    StageManifest, VerifiedFile,
    fingerprint::{InputIdentity, RetainedInput},
    util::compute_sha256,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
};

fn identities(inputs: &[ManifestArtifact]) -> Vec<InputIdentity> {
    let mut values = inputs.iter().map(InputIdentity::from).collect::<Vec<_>>();
    values.sort();
    values
}

fn load_payload(
    run: &RunRecord,
    paths: &ProjectPaths,
    config: &ProjectConfig,
) -> Result<Option<(QuantArtifacts, StageManifest)>> {
    let out_dir = paths.quant_run_dir(&run.run_accession, &config.quant.engine);
    let artifacts = QuantArtifacts {
        run_accession: run.run_accession.clone(),
        abundance_h5: out_dir.join("abundance.h5"),
        run_info_json: out_dir.join("run_info.json"),
        manifest_path: out_dir.join("rnaa_quant_manifest.json"),
        out_dir,
    };
    let Some(manifest) = fs::read(&artifacts.manifest_path)
        .ok()
        .and_then(|bytes| serde_json::from_slice::<StageManifest>(&bytes).ok())
    else {
        return Ok(None);
    };
    if manifest.schema_version != 1
        || manifest.stage != "quant"
        || manifest.id != run.run_accession
        || manifest.input_artifacts.is_empty()
    {
        return Ok(None);
    }
    let required = [
        artifacts.abundance_h5.clone(),
        artifacts.out_dir.join("abundance.tsv"),
        artifacts.run_info_json.clone(),
    ];
    if manifest.output_artifacts.len() != required.len() {
        return Ok(None);
    }
    for path in required {
        let Some(output) = manifest
            .output_artifacts
            .iter()
            .find(|file| Path::new(&file.path) == path)
        else {
            return Ok(None);
        };
        if !output.matches_file()? {
            return Ok(None);
        }
    }
    Ok(Some((artifacts, manifest)))
}

fn scientific_match(
    run: &RunRecord,
    reference: &ReferenceBundle,
    config: &ProjectConfig,
    manifest: &StageManifest,
    inputs: &[ManifestArtifact],
) -> Result<bool> {
    let originals_omitted = manifest.input_artifacts.iter().all(|input| {
        !matches!(
            input.kind.as_str(),
            "FASTQ_SINGLE" | "FASTQ_R1" | "FASTQ_R2"
        )
    }) && manifest
        .input_artifacts
        .iter()
        .any(|input| input.kind.starts_with("FASTQ_TRIMMED_"));
    // Additional originals are provenance observations whose ancestry the caller
    // verifies through the immutable preprocessing producer, not extra engine input.
    let projected = inputs
        .iter()
        .filter(|input| {
            !originals_omitted
                || !matches!(
                    input.kind.as_str(),
                    "FASTQ_SINGLE" | "FASTQ_R1" | "FASTQ_R2"
                )
        })
        .cloned()
        .collect::<Vec<_>>();
    let inputs = projected.as_slice();
    if identities(inputs) != identities(&manifest.input_artifacts) {
        return Ok(false);
    }
    let Some(saved) = manifest.parameters.get("request") else {
        return Ok(false);
    };
    request_matches(
        saved,
        &manifest.parameters["request_sha256"],
        &quant_request(run, reference, config, inputs)?,
    )
}

/// Validate an immutable result's science and complete outputs. The caller must
/// separately prove current source/processing/QC admission and input ownership or
/// authorized pruning before selecting it.
pub fn historical_quant_artifacts(
    run: &RunRecord,
    reference: &ReferenceBundle,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    inputs: &[ManifestArtifact],
) -> Result<Option<QuantArtifacts>> {
    let Some((artifacts, manifest)) = load_payload(run, paths, config)? else {
        return Ok(None);
    };
    if !scientific_match(run, reference, config, &manifest, inputs)? {
        return Ok(None);
    }
    Ok(Some(artifacts))
}

pub fn reconcile_quant_artifacts(
    run: &RunRecord,
    reference: &ReferenceBundle,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    registered_inputs: &[VerifiedFile],
    retained: &[RetainedInput],
) -> Result<Option<QuantArtifacts>> {
    reconcile_quant_artifacts_with_inputs(
        run,
        reference,
        paths,
        config,
        registered_inputs,
        retained,
        None,
    )
}

pub fn reconcile_quant_artifacts_with_inputs(
    run: &RunRecord,
    reference: &ReferenceBundle,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    registered_inputs: &[VerifiedFile],
    retained: &[RetainedInput],
    bound_inputs: Option<&[ManifestArtifact]>,
) -> Result<Option<QuantArtifacts>> {
    let Some((artifacts, manifest)) = load_payload(run, paths, config)? else {
        return Ok(None);
    };
    let inputs = bound_inputs.unwrap_or(&manifest.input_artifacts);
    if !scientific_match(run, reference, config, &manifest, inputs)? {
        return Ok(None);
    }
    let generation = compute_sha256(&artifacts.manifest_path)?;
    let mut registered = BTreeMap::new();
    for file in registered_inputs {
        if registered
            .insert(file.path.display().to_string(), file)
            .is_some()
        {
            return Ok(None);
        }
    }
    let mut paths = BTreeSet::new();
    for input in inputs {
        if !paths.insert(&input.path) {
            return Ok(None);
        }
        let authorized = retained
            .iter()
            .any(|record| record.generation_sha256 == generation && record.artifact == *input);
        if let Some(file) = registered.remove(&input.path) {
            if file.artifact_kind.as_str() != input.kind {
                return Ok(None);
            }
        } else if !authorized {
            return Ok(None);
        }
        if Path::new(&input.path).exists() {
            if !input.matches_file()? {
                return Ok(None);
            }
        } else if !authorized {
            return Ok(None);
        }
    }
    if !registered.is_empty() {
        return Ok(None);
    }
    Ok(Some(artifacts))
}
