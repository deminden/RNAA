use std::collections::BTreeMap;
use std::fs::{self, File, Metadata};
use std::io::Read;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

#[cfg(unix)]
fn same_observation(left: &Metadata, right: &Metadata) -> bool {
    use std::os::unix::fs::MetadataExt;
    left.is_file()
        && right.is_file()
        && left.dev() == right.dev()
        && left.ino() == right.ino()
        && left.len() == right.len()
        && left.mtime() == right.mtime()
        && left.mtime_nsec() == right.mtime_nsec()
        && left.ctime() == right.ctime()
        && left.ctime_nsec() == right.ctime_nsec()
}

#[cfg(not(unix))]
fn same_observation(left: &Metadata, right: &Metadata) -> bool {
    left.is_file()
        && right.is_file()
        && left.len() == right.len()
        && left.modified().ok() == right.modified().ok()
        && left.created().ok() == right.created().ok()
}

// Read no more than the original length plus one byte. The descriptor and
// pathname must still refer to the observed file after hashing. This detects
// ordinary concurrent changes but does not lock out external writers.
fn stable_digest<F: FnMut(u64)>(
    path: &Path,
    expected: Option<u64>,
    mut after_read: F,
) -> Result<Option<(String, u64)>> {
    let before_path = match fs::metadata(path) {
        Ok(metadata) if metadata.is_file() => metadata,
        Ok(_) => return Ok(None),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => {
            return Err(error).with_context(|| format!("failed to inspect {}", path.display()));
        }
    };
    if expected.is_some_and(|size| size != before_path.len()) {
        return Ok(None);
    }
    let mut input =
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let before = input.metadata()?;
    if !same_observation(&before_path, &before) {
        return Ok(None);
    }
    let mut remaining = before.len();
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 64 * 1024];
    while remaining > 0 {
        let limit = usize::try_from(remaining.min(buffer.len() as u64))?;
        let read = input
            .read(&mut buffer[..limit])
            .with_context(|| format!("failed to read {}", path.display()))?;
        if read == 0 {
            return Ok(None);
        }
        hasher.update(&buffer[..read]);
        remaining -= read as u64;
        after_read(before.len() - remaining);
    }
    if input.read(&mut buffer[..1])? != 0 {
        return Ok(None);
    }
    let after = input.metadata()?;
    let after_path = match fs::metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => {
            return Err(error).with_context(|| format!("failed to inspect {}", path.display()));
        }
    };
    if !same_observation(&before, &after) || !same_observation(&before, &after_path) {
        return Ok(None);
    }
    Ok(Some((hex::encode(hasher.finalize()), before.len())))
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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

impl StageManifest {
    /// Strict local-stage validation; remote download/source manifests are not
    /// eligible because a local digest alone cannot attest their provenance.
    pub fn verify_local_artifacts(&self) -> anyhow::Result<bool> {
        if self.schema_version != 1
            || self.input_artifacts.is_empty()
            || self.output_artifacts.is_empty()
        {
            return Ok(false);
        }
        for artifact in self.input_artifacts.iter().chain(&self.output_artifacts) {
            if !artifact.matches_file()? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

impl ManifestArtifact {
    pub fn snapshot(kind: &str, path: &std::path::Path) -> anyhow::Result<Self> {
        Self::snapshot_with_size(kind, path, None)
    }

    /// Capture a regular file only if it still has the size admitted by the caller.
    pub fn snapshot_exact(kind: &str, path: &std::path::Path, bytes: u64) -> anyhow::Result<Self> {
        Self::snapshot_with_size(kind, path, Some(bytes))
    }

    fn snapshot_with_size(
        kind: &str,
        path: &std::path::Path,
        expected: Option<u64>,
    ) -> anyhow::Result<Self> {
        let (checksum, bytes) = stable_digest(path, expected, |_| {})?.with_context(|| {
            format!(
                "artifact changed, differs from its admitted size, or is not a regular file: {}",
                path.display()
            )
        })?;
        Ok(Self {
            kind: kind.into(),
            path: path.display().to_string(),
            checksum_type: "sha256".into(),
            checksum,
            bytes,
        })
    }
    pub fn matches_file(&self) -> anyhow::Result<bool> {
        let path = std::path::Path::new(&self.path);
        if self.checksum_type != "sha256" {
            return Ok(false);
        }
        Ok(stable_digest(path, Some(self.bytes), |_| {})?
            .is_some_and(|(checksum, bytes)| checksum == self.checksum && bytes == self.bytes))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[cfg(unix)]
    fn snapshots_reject_growth_truncation_and_path_replacement_during_reads() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let path = temp.path().join("input");
        let original = vec![b'a'; 128 * 1024];
        for change in ["grow", "truncate", "replace", "rewrite"] {
            fs::write(&path, &original)?;
            let mut changed = false;
            let mut observed = 0;
            let result = stable_digest(&path, None, |bytes| {
                observed = bytes;
                if changed {
                    return;
                }
                changed = true;
                match change {
                    "grow" => {
                        use std::io::Write;
                        fs::OpenOptions::new()
                            .append(true)
                            .open(&path)
                            .unwrap()
                            .write_all(b"extra")
                            .unwrap();
                    }
                    "truncate" => fs::write(&path, b"short").unwrap(),
                    "replace" => {
                        let replacement = temp.path().join("replacement");
                        fs::write(&replacement, &original).unwrap();
                        fs::rename(replacement, &path).unwrap();
                    }
                    "rewrite" => fs::write(&path, vec![b'b'; original.len()]).unwrap(),
                    _ => unreachable!(),
                }
            })?;
            assert!(changed, "{change} was not exercised");
            assert!(result.is_none(), "{change} was not detected");
            assert!(
                observed <= original.len() as u64,
                "read beyond the original bound"
            );
        }
        fs::write(&path, &original)?;
        let witness = ManifestArtifact::snapshot("TEST", &path)?;
        assert_eq!(witness.bytes, original.len() as u64);
        assert_eq!(witness.checksum, crate::util::compute_sha256(&path)?);
        assert_eq!(
            ManifestArtifact::snapshot_exact("TEST", &path, witness.bytes)?,
            witness
        );
        assert!(ManifestArtifact::snapshot_exact("TEST", &path, witness.bytes - 1).is_err());
        assert!(witness.matches_file()?);
        fs::remove_file(&path)?;
        assert!(!witness.matches_file()?);
        assert!(ManifestArtifact::snapshot("TEST", &path).is_err());
        Ok(())
    }

    #[test]
    fn verifies_inputs_and_outputs_including_same_size_corruption() -> anyhow::Result<()> {
        let root = tempfile::TempDir::new()?;
        let input = root.path().join("input");
        let output = root.path().join("output");
        std::fs::write(&input, b"source")?;
        std::fs::write(&output, b"result")?;
        let artifact = |path: &std::path::Path| -> anyhow::Result<ManifestArtifact> {
            Ok(ManifestArtifact {
                kind: "TEST".into(),
                path: path.display().to_string(),
                checksum_type: "sha256".into(),
                checksum: crate::util::compute_sha256(path)?,
                bytes: crate::util::file_size(path)?,
            })
        };
        let manifest = StageManifest {
            schema_version: 1,
            stage: "test".into(),
            id: "test".into(),
            started_at: String::new(),
            finished_at: String::new(),
            parameters: serde_json::Value::Null,
            tool_versions: BTreeMap::new(),
            input_artifacts: vec![artifact(&input)?],
            output_artifacts: vec![artifact(&output)?],
            notes: vec![],
        };
        assert!(manifest.verify_local_artifacts()?);
        std::fs::write(&output, b"broken")?;
        assert!(!manifest.verify_local_artifacts()?);
        std::fs::write(&output, b"result")?;
        std::fs::write(&input, b"edited")?;
        assert!(!manifest.verify_local_artifacts()?);
        Ok(())
    }
}
