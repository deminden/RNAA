//! Exact immutable producer edges, kept distinct from external source witnesses.
use super::{GenerationHead, Seal};
use crate::{Database, ManifestArtifact};
use anyhow::{Context, Result, bail};
use rusqlite::{OptionalExtension, params};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenerationDependency {
    pub producer_generation_id: String,
    pub producer_seal_sha256: String,
    pub input: ManifestArtifact,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenerationLineage {
    pub schema_version: u32,
    pub recorded: bool,
    pub source: String,
    pub head: GenerationHead,
    pub dependencies: Vec<GenerationDependency>,
    pub external_inputs: Vec<ManifestArtifact>,
}

impl Database {
    pub fn committed_generation_head(&self, id: &str) -> Result<Option<GenerationHead>> {
        Ok(self.committed_seal(id)?.map(|seal| seal.head))
    }
    pub(super) fn resolve_generation_dependencies(
        &self,
        inputs: &[ManifestArtifact],
        own_payload: &Path,
    ) -> Result<Vec<GenerationDependency>> {
        let mut dependencies = Vec::new();
        let mut seen = BTreeSet::new();
        for input in inputs {
            let canonical = if Path::new(&input.path).exists() {
                std::fs::canonicalize(&input.path)?
            } else {
                // Reuse may bind a sealed producer whose FASTQ was deliberately
                // pruned. Its exact indexed ownership is still required below.
                Path::new(&input.path).to_path_buf()
            };
            if canonical
                .components()
                .any(|part| part.as_os_str() == ".rnaa_generations")
                && canonical != Path::new(&input.path)
            {
                bail!(
                    "managed generation inputs must use their sealed canonical paths: {}",
                    input.path
                );
            }
            if !seen.insert((&input.kind, &input.path)) {
                bail!("duplicate generation input role/path");
            }
            let owners = self.read(|conn| {
                let mut query = conn.prepare("SELECT f.generation_id FROM generation_files f JOIN artifact_generations g ON g.generation_id=f.generation_id WHERE g.project_id=?1 AND g.state='COMMITTED' AND f.path=?2 AND f.sha256=?3 AND f.bytes=?4 AND NOT EXISTS(SELECT 1 FROM generation_retirements r WHERE r.generation_id=g.generation_id)")?;
                Ok(query.query_map(params![self.project_id(),input.path,input.checksum,i64::try_from(input.bytes)?], |row| row.get::<_,String>(0))?.collect::<std::result::Result<Vec<_>,_>>()?)
            })?;
            if owners.len() > 1 {
                bail!("generation input has ambiguous producer ownership");
            }
            if let Some(owner) = owners.first() {
                let producer = self
                    .committed_seal(owner)?
                    .context("generation input producer disappeared")?;
                if !producer.files.iter().any(|file| same_file(file, input)) {
                    bail!("producer inventory differs from its indexed input");
                }
                dependencies.push(GenerationDependency {
                    producer_generation_id: owner.clone(),
                    producer_seal_sha256: seal_digest(&producer)?,
                    input: input.clone(),
                });
            } else if Path::new(&input.path)
                .components()
                .any(|component| component.as_os_str() == ".rnaa_generations")
            {
                // A known generation path must never silently become an unmanaged source.
                // Consumer-owned frozen metadata has not been committed yet.
                if !Path::new(&input.path).starts_with(own_payload) {
                    bail!(
                        "generation input lacks a committed producer: {}",
                        input.path
                    );
                }
            }
        }
        dependencies.sort_by(|a, b| {
            (&a.producer_generation_id, &a.input.kind, &a.input.path).cmp(&(
                &b.producer_generation_id,
                &b.input.kind,
                &b.input.path,
            ))
        });
        Ok(dependencies)
    }

    pub(super) fn verify_dependency_catalog(
        &self,
        dependencies: &[GenerationDependency],
    ) -> Result<()> {
        self.read(|conn| {
            for dependency in dependencies {
                let json: Option<String> = conn.query_row("SELECT g.inventory_json FROM generation_files f JOIN artifact_generations g ON g.generation_id=f.generation_id WHERE g.project_id=?1 AND g.state='COMMITTED' AND f.generation_id=?2 AND f.path=?3 AND f.sha256=?4 AND f.bytes=?5",
                    params![self.project_id(),dependency.producer_generation_id,dependency.input.path,dependency.input.checksum,i64::try_from(dependency.input.bytes)?], |row| row.get(0)).optional()?;
                let producer: Seal = serde_json::from_str(&json.context("generation dependency changed before publication")?)?;
                if seal_digest(&producer)? != dependency.producer_seal_sha256 { bail!("generation dependency seal changed before publication"); }
            }
            Ok(())
        })
    }

    pub(super) fn lineage_matches_catalog(&self, seal: &Seal) -> Result<bool> {
        let indexed = self.read(|conn| {
            let mut query = conn.prepare("SELECT path,sha256,bytes FROM generation_files WHERE generation_id=?1 ORDER BY path")?;
            Ok(query.query_map([&seal.head.generation_id], |row| Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,i64>(2)?)))?.collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        let mut expected = seal
            .files
            .iter()
            .map(|file| {
                Ok((
                    file.path.clone(),
                    file.checksum.clone(),
                    i64::try_from(file.bytes)?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        expected.sort();
        if indexed != expected {
            return Ok(false);
        }
        if seal.lineage_version == 0 {
            return Ok(seal.dependencies.is_empty() && self.imported_lineage(seal).is_ok());
        }
        if seal.lineage_version != 1 {
            return Ok(false);
        }
        let Ok(recorded) = self.recorded_dependencies(&seal.head.generation_id) else {
            return Ok(false);
        };
        if recorded != seal.dependencies || self.verify_dependency_catalog(&recorded).is_err() {
            return Ok(false);
        }
        Ok(recorded.iter().all(|dependency| {
            seal.inputs.contains(&dependency.input)
                && dependency.producer_generation_id != seal.head.generation_id
        }))
    }

    fn recorded_dependencies(&self, id: &str) -> Result<Vec<GenerationDependency>> {
        self.read(|conn| {
            let mut query = conn.prepare("SELECT producer_generation_id,artifact_json,role,path,producer_seal_sha256 FROM generation_dependencies WHERE consumer_generation_id=?1 ORDER BY producer_generation_id,role,path")?;
            query.query_map([id], |row| Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,String>(2)?,row.get::<_,String>(3)?,row.get::<_,String>(4)?)))?
                .map(|row| {let (producer,json,role,path,digest)=row?; let input: ManifestArtifact=serde_json::from_str(&json)?;
                    if input.kind != role || input.path != path { bail!("generation dependency index differs from its input witness"); }
                    Ok(GenerationDependency {producer_generation_id:producer,producer_seal_sha256:digest,input})}).collect()
        })
    }

    pub fn generation_lineage(&self, head: &GenerationHead) -> Result<GenerationLineage> {
        let seal = self
            .committed_seal(&head.generation_id)?
            .context("generation is not committed in this project")?;
        if seal.head != *head || !self.lineage_matches_catalog(&seal)? {
            bail!("generation lineage does not match its sealed catalog");
        }
        let (recorded, dependencies) = self.effective_dependencies(&seal)?;
        let inputs = self.ancestry_inputs(&seal)?;
        let external_inputs = inputs
            .into_iter()
            .filter(|input| {
                !Path::new(&input.path).starts_with(&head.payload_path)
                    && !dependencies
                        .iter()
                        .any(|dependency| dependency.input == *input)
            })
            .collect();
        Ok(GenerationLineage {
            schema_version: 1,
            recorded,
            source: if seal.lineage_version == 1 {
                "native"
            } else if recorded {
                "imported"
            } else {
                "unresolved_legacy"
            }
            .into(),
            head: head.clone(),
            dependencies,
            external_inputs,
        })
    }

    /// Immutable downstream history; current heads and settings determine admission separately.
    pub fn dependent_generations(&self, producer: &GenerationHead) -> Result<Vec<GenerationHead>> {
        self.generation_lineage(producer)?;
        let ids = self.read(|conn| {
            let mut query = conn.prepare("WITH RECURSIVE all_edges AS (SELECT consumer_generation_id,producer_generation_id FROM generation_dependencies UNION SELECT consumer_generation_id,producer_generation_id FROM generation_import_dependencies), dependents(id) AS (SELECT consumer_generation_id FROM all_edges WHERE producer_generation_id=?1 UNION SELECT d.consumer_generation_id FROM all_edges d JOIN dependents p ON d.producer_generation_id=p.id) SELECT id FROM dependents ORDER BY id")?;
            Ok(query.query_map([&producer.generation_id], |row| row.get::<_,String>(0))?.collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        ids.iter()
            .map(|id| {
                self.committed_seal(id)?
                    .map(|seal| seal.head)
                    .context("dependent generation is not committed")
            })
            .collect()
    }
}
fn same_file(a: &ManifestArtifact, b: &ManifestArtifact) -> bool {
    a.path == b.path
        && a.checksum_type == b.checksum_type
        && a.checksum == b.checksum
        && a.bytes == b.bytes
}
fn seal_digest(seal: &Seal) -> Result<String> {
    Ok(hex::encode(Sha256::digest(serde_json::to_vec(seal)?)))
}

/// Resolve aliases even after the final file (or a symlink target) was pruned.
/// The nearest existing prefix is canonical; absent tail names retain identity.
pub(crate) fn canonical_missing_path(path: &Path) -> Result<PathBuf> {
    fn resolve(path: &Path, depth: usize) -> Result<PathBuf> {
        if depth > 40 {
            bail!("too many source path aliases");
        }
        match fs::canonicalize(path) {
            Ok(path) => return Ok(path),
            Err(error)
                if matches!(
                    error.kind(),
                    std::io::ErrorKind::NotFound | std::io::ErrorKind::NotADirectory
                ) => {}
            Err(error) => return Err(error.into()),
        }
        let mut ancestor = path.to_path_buf();
        let mut tail = Vec::new();
        loop {
            if ancestor.as_os_str().is_empty() {
                ancestor = PathBuf::from(".");
            }
            let metadata = match ancestor.symlink_metadata() {
                Ok(metadata) => Some(metadata),
                Err(error)
                    if matches!(
                        error.kind(),
                        std::io::ErrorKind::NotFound | std::io::ErrorKind::NotADirectory
                    ) =>
                {
                    None
                }
                Err(error) => return Err(error.into()),
            };
            if metadata.is_some_and(|metadata| metadata.file_type().is_symlink()) {
                let target = fs::read_link(&ancestor)?;
                let mut target = if target.is_absolute() {
                    target
                } else {
                    ancestor
                        .parent()
                        .context("source alias has no parent")?
                        .join(target)
                };
                for name in tail.iter().rev() {
                    target.push(name);
                }
                return resolve(&target, depth + 1);
            }
            match fs::canonicalize(&ancestor) {
                Ok(mut prefix) => {
                    for name in tail.iter().rev() {
                        prefix.push(name);
                    }
                    return Ok(prefix);
                }
                Err(error)
                    if matches!(
                        error.kind(),
                        std::io::ErrorKind::NotFound | std::io::ErrorKind::NotADirectory
                    ) => {}
                Err(error) => return Err(error.into()),
            }
            tail.push(
                ancestor
                    .file_name()
                    .context("cannot resolve source path prefix")?
                    .to_os_string(),
            );
            ancestor = ancestor
                .parent()
                .context("source path has no existing prefix")?
                .to_path_buf();
        }
    }
    resolve(path, 0)
}
