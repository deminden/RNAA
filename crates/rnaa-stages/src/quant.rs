mod external;
mod index_group;
mod index_plan;
mod input;
mod native;
mod policy;
mod reconcile;
mod resources;
pub use external::{declare_external_quant, external_quant_matches};
pub use index_group::{QuantIndexGroup, QuantIndexStats};
pub use index_plan::QuantIndexPlan;
pub use policy::{
    quant_policy_fingerprint, source_binding_current, verify_quant_catalog, verify_quant_current,
};
pub use reconcile::{
    historical_quant_artifacts, reconcile_quant_artifacts, reconcile_quant_artifacts_with_inputs,
};
use std::collections::BTreeMap;
use std::fs;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result};
use rnaa_core::config::ProjectConfig;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{QuantArtifacts, ReferenceBundle, RunRecord, VerifiedFile};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::Quantifier;
use rnaa_core::util::{compute_sha256, file_size, now_rfc3339, write_json_pretty};
use serde_json::json;

#[derive(Debug, Clone, Default)]
pub struct RustKallistoQuantifier;

impl Quantifier for RustKallistoQuantifier {
    fn quantify(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        reference: &ReferenceBundle,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<QuantArtifacts> {
        config.quant.validate()?;
        let group = QuantIndexGroup::new(ManifestArtifact::snapshot(
            "REFERENCE_INDEX",
            &reference.kallisto_index_path,
        )?)?;
        group.quantify(run, fastqs, reference, paths, config)
    }
}

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct QuantSettings {
    layout: rnaa_core::LibraryLayout,
    strand: rnaa_core::config::QuantStrand,
    single_fragment_length: Option<f64>,
    single_fragment_sd: Option<f64>,
    preprocess: bool,
}

fn quant_request(
    run: &RunRecord,
    reference: &ReferenceBundle,
    config: &ProjectConfig,
    inputs: &[ManifestArtifact],
) -> Result<rnaa_core::fingerprint::StageRequest<QuantSettings>> {
    use rnaa_core::fingerprint::{InputIdentity, StageRequest};
    let single = run.library_layout == rnaa_core::LibraryLayout::Single;
    Ok(StageRequest {
        schema_version: 1,
        stage: "quant".into(),
        algorithm_version: 4,
        engines: BTreeMap::from([
            ("kallistors".into(), crate::engines::KALLISTORS.into()),
            (
                "dependency_lock".into(),
                crate::engines::DEPENDENCY_LOCK_SHA256.into(),
            ),
        ]),
        parameters: QuantSettings {
            layout: run.library_layout,
            strand: config.quant.strand,
            single_fragment_length: single.then_some(config.quant.single_fragment_length),
            single_fragment_sd: single.then_some(config.quant.single_fragment_sd),
            preprocess: config.quant.preprocess,
        },
        inputs: inputs.iter().map(InputIdentity::from).collect(),
        dependencies: BTreeMap::from([(
            "reference_index".into(),
            compute_sha256(&reference.kallisto_index_path)?,
        )]),
    })
}

/// Revision 3 coupled gate policy to quantification. Validate its original
/// digest before projecting only those two fields out of the scientific key.
/// Keep the original manifest untouched: retention receipts bind its exact bytes.
fn request_matches(
    saved: &serde_json::Value,
    digest: &serde_json::Value,
    expected: &rnaa_core::fingerprint::StageRequest<QuantSettings>,
) -> Result<bool> {
    use rnaa_core::fingerprint::StageRequest;
    let Ok(mut request) = serde_json::from_value::<StageRequest<serde_json::Value>>(saved.clone())
    else {
        return Ok(false);
    };
    // A typed parse must not silently discard unrecognized top-level fields.
    if serde_json::to_value(&request)? != *saved || request.digest()? != *digest {
        return Ok(false);
    }
    if request.algorithm_version == 3 {
        let Some(parameters) = request.parameters.as_object_mut() else {
            return Ok(false);
        };
        let strict = parameters.remove("preprocess_strict");
        let gate = parameters.remove("qc_gate");
        if !matches!(
            strict,
            Some(serde_json::Value::Null | serde_json::Value::Bool(_))
        ) {
            return Ok(false);
        }
        if !gate.is_some_and(|value| {
            value.is_null() || serde_json::from_value::<rnaa_core::PreprocessQcGate>(value).is_ok()
        }) {
            return Ok(false);
        }
        request.algorithm_version = 4;
    }
    if serde_json::from_value::<QuantSettings>(request.parameters.clone()).is_err() {
        return Ok(false);
    }
    Ok(request.digest()? == expected.digest()?)
}

pub fn quant_request_digest(
    run: &RunRecord,
    reference: &ReferenceBundle,
    config: &ProjectConfig,
    inputs: &[ManifestArtifact],
) -> Result<String> {
    quant_request(run, reference, config, inputs)?.digest()
}

pub fn quant_request_value(
    run: &RunRecord,
    reference: &ReferenceBundle,
    config: &ProjectConfig,
    inputs: &[ManifestArtifact],
) -> Result<rnaa_core::fingerprint::StageRequest<serde_json::Value>> {
    Ok(serde_json::from_value(serde_json::to_value(
        quant_request(run, reference, config, inputs)?,
    )?)?)
}

/// Capture registry inputs before execution. Read digests from payloads rather
/// than trusting potentially obsolete database checksums or sizes.
pub fn snapshot_quant_inputs(inputs: &[VerifiedFile]) -> Result<Vec<ManifestArtifact>> {
    inputs
        .iter()
        .map(|file| ManifestArtifact::snapshot(file.artifact_kind.as_str(), &file.path))
        .collect()
}

pub fn attach_quant_request(
    artifacts: &QuantArtifacts,
    run: &RunRecord,
    reference: &ReferenceBundle,
    config: &ProjectConfig,
    inputs: Vec<ManifestArtifact>,
) -> Result<()> {
    for input in &inputs {
        if !input.matches_file()? {
            anyhow::bail!(
                "quantification input changed while computing: {}",
                input.path
            );
        }
    }
    let mut manifest: StageManifest = serde_json::from_slice(&fs::read(&artifacts.manifest_path)?)?;
    let request = quant_request(run, reference, config, &inputs)?;
    manifest.parameters["request_sha256"] = serde_json::to_value(request.digest()?)?;
    manifest.parameters["request"] = serde_json::to_value(request)?;
    // Preserve historical preprocessing provenance when quantification is repeated
    // from retained trimmed reads after the original FASTQs have been removed.
    manifest.parameters["input_lineage"] = serde_json::to_value(&manifest.input_artifacts)?;
    manifest.input_artifacts = inputs;
    write_json_pretty(&artifacts.manifest_path, &manifest)
}

#[cfg(test)]
mod request_tests {
    use super::*;
    use rnaa_core::fingerprint::StageRequest;
    #[test]
    fn legacy_qc_policy_can_change_without_rewriting_or_trusting_a_tampered_request() -> Result<()>
    {
        let expected = StageRequest {
            schema_version: 1,
            stage: "quant".into(),
            algorithm_version: 4,
            engines: BTreeMap::new(),
            parameters: QuantSettings {
                layout: rnaa_core::LibraryLayout::Single,
                strand: rnaa_core::config::QuantStrand::Unstranded,
                single_fragment_length: Some(200.0),
                single_fragment_sd: Some(20.0),
                preprocess: true,
            },
            inputs: Vec::new(),
            dependencies: BTreeMap::new(),
        };
        let mut legacy = serde_json::to_value(&expected)?;
        legacy["algorithm_version"] = serde_json::json!(3);
        legacy["parameters"]["preprocess_strict"] = serde_json::json!(true);
        legacy["parameters"]["qc_gate"] =
            serde_json::to_value(rnaa_core::PreprocessQcGate::default())?;
        let digest = |value: &serde_json::Value| -> Result<serde_json::Value> {
            Ok(serde_json::json!(
                serde_json::from_value::<StageRequest<serde_json::Value>>(value.clone())?
                    .digest()?
            ))
        };
        let original = legacy.clone();
        let original_digest = digest(&legacy)?;
        assert!(request_matches(&legacy, &original_digest, &expected)?);
        assert_eq!(legacy, original);
        legacy["parameters"]["preprocess_strict"] = serde_json::json!(false);
        legacy["parameters"]["qc_gate"]["min_pass_rate"] = serde_json::json!(0.99);
        assert!(
            !request_matches(&legacy, &original_digest, &expected)?,
            "validate the saved original digest before projection"
        );
        assert!(request_matches(&legacy, &digest(&legacy)?, &expected)?);
        legacy["parameters"]["single_fragment_length"] = serde_json::json!(250.0);
        assert!(!request_matches(&legacy, &digest(&legacy)?, &expected)?);
        legacy = original.clone();
        legacy["parameters"]["unknown_count_option"] = serde_json::json!(true);
        assert!(!request_matches(&legacy, &digest(&legacy)?, &expected)?);
        legacy = original;
        legacy["unexpected_top_level"] = serde_json::json!(true);
        assert!(!request_matches(&legacy, &digest(&legacy)?, &expected)?);
        let current = serde_json::to_value(&expected)?;
        assert!(request_matches(&current, &digest(&current)?, &expected)?);
        Ok(())
    }
}
