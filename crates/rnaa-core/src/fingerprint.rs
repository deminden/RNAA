//! Semantic request identities and explicit evidence for intentionally removed inputs.
use crate::manifest::ManifestArtifact;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct InputIdentity {
    pub role: String,
    pub name: String,
    pub sha256: String,
    pub bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageRequest<P> {
    pub schema_version: u32,
    pub stage: String,
    pub algorithm_version: u32,
    pub engines: BTreeMap<String, String>,
    pub parameters: P,
    pub inputs: Vec<InputIdentity>,
    pub dependencies: BTreeMap<String, String>,
}

impl<P: Serialize> StageRequest<P> {
    /// Sort inputs without deduplicating: two identical lanes still contain twice
    /// the reads. Paths and execution resources are provenance, not semantic inputs.
    pub fn digest(&self) -> Result<String> {
        let mut value = serde_json::to_value(self)?;
        let mut inputs = self.inputs.clone();
        inputs.sort();
        value["inputs"] = serde_json::to_value(inputs)?;
        Ok(hex::encode(Sha256::digest(serde_json::to_vec(&value)?)))
    }
}

impl From<&ManifestArtifact> for InputIdentity {
    fn from(artifact: &ManifestArtifact) -> Self {
        Self {
            role: artifact.kind.clone(),
            name: std::path::Path::new(&artifact.path)
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .into_owned(),
            sha256: artifact.checksum.clone(),
            bytes: artifact.bytes,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RetainedInput {
    pub generation_sha256: String,
    pub artifact: ManifestArtifact,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fingerprints_preserve_multiplicity_and_ignore_execution_location_and_order() -> Result<()> {
        let file = |path: &str| ManifestArtifact {
            kind: "FASTQ_R1".into(),
            path: path.into(),
            checksum_type: "sha256".into(),
            checksum: "a".repeat(64),
            bytes: 100,
        };
        let a = InputIdentity::from(&file("/project/a.fastq"));
        assert_eq!(a, InputIdentity::from(&file("/another/project/a.fastq")));
        let mut request = StageRequest {
            schema_version: 1,
            stage: "quant".into(),
            algorithm_version: 1,
            engines: BTreeMap::new(),
            parameters: true,
            inputs: vec![a.clone()],
            dependencies: BTreeMap::new(),
        };
        let single = request.digest()?;
        request.inputs.push(a);
        assert_ne!(request.digest()?, single);
        let mut b = request.inputs[0].clone();
        b.name = "b.fastq".into();
        request.inputs.push(b);
        let before = request.digest()?;
        request.inputs.reverse();
        assert_eq!(request.digest()?, before);
        request.dependencies.insert("index".into(), "new".into());
        assert_ne!(request.digest()?, before);
        Ok(())
    }
}
