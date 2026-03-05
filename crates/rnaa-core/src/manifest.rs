use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ManifestArtifact {
    pub kind: String,
    pub path: String,
    pub checksum_type: String,
    pub checksum: String,
    pub bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageManifest {
    pub schema_version: u32,
    pub stage: String,
    pub id: String,
    pub started_at: String,
    pub finished_at: String,
    pub parameters: serde_json::Value,
    pub tool_versions: BTreeMap<String, String>,
    pub input_artifacts: Vec<ManifestArtifact>,
    pub output_artifacts: Vec<ManifestArtifact>,
    pub notes: Vec<String>,
}
