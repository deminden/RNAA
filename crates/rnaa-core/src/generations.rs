//! Immutable stage output sets. SQLite is the publication boundary; payload
//! paths never change, and readers pin a committed generation before opening files.
use crate::{
    ArtifactRecord, Database, ManifestArtifact, ProjectPaths,
    jobs::JobOutcome,
    util::{now_rfc3339, write_json_pretty},
};
use anyhow::{Context, Result, bail};
use rusqlite::{OptionalExtension, params};
use serde::{Deserialize, Serialize};
use std::{
    fs,
    path::{Path, PathBuf},
};
mod bindings;
mod legacy;
mod lineage;
pub use legacy::{LineageImportEntry, LineageImportReport};
pub(crate) mod pins;
mod registry;
mod requests;
mod reuse;
mod verification;
pub use bindings::GenerationInputBinding;
pub(crate) use lineage::canonical_missing_path;
pub use lineage::{GenerationDependency, GenerationLineage};
pub use pins::GenerationStorePin;
mod gc;
pub use gc::{GenerationGcEntry, GenerationGcOptions, GenerationGcReport};
pub use gc::{UnpublishedCleanupEntry, UnpublishedCleanupReport};
mod fastq_retention;
mod file_trash;
pub use file_trash::{CleanupPolicy, FileCleanupEntry, FileCleanupOptions, FileCleanupReport};
pub use reuse::VerifiedGenerationReuse;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum OutputGroup {
    Refs,
    Download,
    Preprocess,
    Quant,
    De,
    Corr,
}
impl OutputGroup {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Refs => "refs",
            Self::Download => "download",
            Self::Preprocess => "preprocess",
            Self::Quant => "quant",
            Self::De => "de",
            Self::Corr => "corr",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenerationHead {
    pub generation_id: String,
    pub output_group: OutputGroup,
    pub stage: String,
    pub subject: String,
    pub attempt_id: String,
    pub payload_path: PathBuf,
    pub request_sha256: String,
}

impl GenerationHead {
    pub fn paths(&self, base: &ProjectPaths) -> Result<ProjectPaths> {
        base.for_generation(self)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Seal {
    schema_version: u32,
    head: GenerationHead,
    files: Vec<ManifestArtifact>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    inputs: Vec<ManifestArtifact>,
    #[serde(default, skip_serializing_if = "is_zero")]
    lineage_version: u32,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    dependencies: Vec<GenerationDependency>,
}
fn is_zero(value: &u32) -> bool {
    *value == 0
}

pub struct Generation {
    db: Database,
    head: GenerationHead,
    seal: Option<Seal>,
    committed: bool,
    scientific_request: Option<serde_json::Value>,
}

impl Generation {
    pub fn begin(db: &Database, group: OutputGroup, stage: &str) -> Result<Self> {
        let identity = db
            .job
            .as_ref()
            .context("generation creation requires an owned job")?;
        if !matches!(
            (group, stage),
            (OutputGroup::Refs, "refs")
                | (OutputGroup::Download, "download")
                | (OutputGroup::Preprocess, "preprocess")
                | (OutputGroup::Quant, "quant")
                | (OutputGroup::De, "normalize" | "deseq2")
                | (OutputGroup::Corr, "corr")
        ) {
            bail!(
                "stage {stage} cannot publish the {} output group",
                group.as_str()
            );
        }
        let subject = db.read(|connection| {
            let (actual, run, resource): (String, Option<String>, String) = connection.query_row(
                "SELECT stage, run_accession, resource FROM job_attempts WHERE attempt_id = ?1",
                params![identity.attempt_id],
                |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
            )?;
            if actual != stage {
                bail!("generation stage does not match its owning attempt");
            }
            match group {
                OutputGroup::Refs => {
                    if run.is_some() {
                        bail!("reference generation cannot belong to a run");
                    }
                    Ok(resource
                        .strip_prefix("reference:")
                        .context("reference job requires a reference resource")?
                        .to_string())
                }
                OutputGroup::Download | OutputGroup::Preprocess | OutputGroup::Quant => {
                    run.context("ingest generation requires a run subject")
                }
                OutputGroup::De | OutputGroup::Corr => {
                    if run.is_some() {
                        bail!("analysis generation requires a project subject");
                    }
                    Ok(db.project_id().into())
                }
            }
        })?;
        if matches!(
            group,
            OutputGroup::Refs
                | OutputGroup::Download
                | OutputGroup::Preprocess
                | OutputGroup::Quant
        ) && (Path::new(&subject).components().count() != 1
            || !matches!(
                Path::new(&subject).components().next(),
                Some(std::path::Component::Normal(_))
            ))
        {
            bail!("run subject must be a single path component");
        }
        let id = uuid::Uuid::new_v4().to_string();
        let root = fs::canonicalize(db.path().parent().context("database has no project root")?)?;
        let payload = root.join(".rnaa_generations").join(&id).join("payload");
        fs::create_dir_all(&payload)?;
        fs::create_dir_all(payload.join(group.as_str()))?;
        if group == OutputGroup::Refs {
            fs::create_dir_all(payload.join("manifests/refs"))?;
        }
        if group == OutputGroup::Download {
            fs::create_dir_all(payload.join("manifests/download"))?;
        }
        if group == OutputGroup::Corr {
            fs::create_dir_all(payload.join("manifests/corr"))?;
        }
        let head = GenerationHead {
            generation_id: id,
            output_group: group,
            stage: stage.into(),
            subject,
            attempt_id: identity.attempt_id.clone(),
            payload_path: payload,
            request_sha256: String::new(),
        };
        db.write(|connection| {
            connection.execute("INSERT INTO artifact_generations(generation_id, project_id, output_group, stage,
                subject, attempt_id, payload_path, state, created_at) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, 'BUILDING', ?8)",
                params![head.generation_id, db.project_id(), group.as_str(), stage, head.subject, head.attempt_id, head.payload_path.to_string_lossy(), now_rfc3339()])?;
            Ok(())
        })?;
        Ok(Self {
            db: db.clone(),
            head,
            seal: None,
            committed: false,
            scientific_request: None,
        })
    }

    pub fn head(&self) -> &GenerationHead {
        &self.head
    }

    pub fn paths(&self, base: &ProjectPaths) -> Result<ProjectPaths> {
        base.for_generation(&self.head)
    }

    /// Verify complete required output coverage, then flush every payload file and
    /// directory before recording the sealed inventory. The seal is outside payload
    /// so it never includes a digest of itself.
    pub fn seal(&mut self, required: &[PathBuf], request: &str) -> Result<()> {
        self.seal_with_inputs(required, request, &[])
    }

    /// Store the full scientific request with the committed result, using its
    /// digest as the immutable identity. Physical input witnesses remain separate.
    pub fn seal_request<P: Serialize>(
        &mut self,
        required: &[PathBuf],
        request: &crate::fingerprint::StageRequest<P>,
        inputs: &[ManifestArtifact],
    ) -> Result<()> {
        if request.stage != self.head.stage {
            bail!("scientific request stage differs from its generation");
        }
        let value = serde_json::to_value(request)?;
        requests::canonical_request(&value)?;
        self.seal_with_inputs(required, &request.digest()?, inputs)?;
        self.scientific_request = Some(value);
        Ok(())
    }

    pub fn seal_with_inputs(
        &mut self,
        required: &[PathBuf],
        request: &str,
        inputs: &[ManifestArtifact],
    ) -> Result<()> {
        self.scientific_request = None;
        self.db.check_job()?;
        for input in inputs {
            if !input.matches_file()? {
                bail!("generation input changed before sealing: {}", input.path);
            }
        }
        let files = inventory(&self.head.payload_path, true)?;
        if files.is_empty() || required.is_empty() {
            bail!("cannot publish an empty generation");
        }
        for required in required {
            if !files
                .iter()
                .any(|file| Path::new(&file.path) == required && file.bytes > 0)
            {
                bail!(
                    "required output missing from generation: {}",
                    required.display()
                );
            }
        }
        self.head.request_sha256 = request.into();
        let seal = Seal {
            schema_version: 1,
            head: self.head.clone(),
            files,
            inputs: inputs.to_vec(),
            lineage_version: 1,
            dependencies: self
                .db
                .resolve_generation_dependencies(inputs, &self.head.payload_path)?,
        };
        let root = self
            .head
            .payload_path
            .parent()
            .context("generation has no directory")?;
        write_json_pretty(&root.join("seal.json"), &seal)?;
        fs::File::open(root)?.sync_all()?;
        fs::File::open(root.parent().context("generation store has no root")?)?.sync_all()?;
        fs::File::open(
            root.parent()
                .and_then(Path::parent)
                .context("generation store has no project root")?,
        )?
        .sync_all()?;
        self.db.write(|connection| {
            let expected: Option<String> = connection.query_row("SELECT request_sha256 FROM job_attempts WHERE attempt_id = ?1", params![self.head.attempt_id], |row| row.get(0))?;
            if expected.as_deref() != Some(request) { bail!("generation does not match the attempt's frozen request"); }
            let changed = connection.execute("UPDATE artifact_generations SET state = 'SEALED', request_sha256 = ?1, inventory_json = ?2
                WHERE generation_id = ?3 AND state = 'BUILDING'", params![request, serde_json::to_string(&seal)?, self.head.generation_id])?;
            if changed != 1 { bail!("generation is no longer building"); }
            Ok(())
        })?;
        self.seal = Some(seal);
        Ok(())
    }

    /// Hash and validate before beginning the short publication transaction. The
    /// callback performs only state/event changes; artifact registration and attempt
    /// success belong to the same commit as the authoritative generation pointer.
    pub fn commit<T>(
        &mut self,
        records: &[ArtifactRecord],
        complete: impl FnOnce(&Database) -> Result<T>,
    ) -> Result<T> {
        self.commit_prepared(records, |_| Ok(()), complete)
    }

    /// Register prepared auxiliary database records in the same transaction, before
    /// artifacts that refer to them. Both callbacks must avoid filesystem work.
    pub fn commit_prepared<T>(
        &mut self,
        records: &[ArtifactRecord],
        register_dependencies: impl FnOnce(&Database) -> Result<()>,
        complete: impl FnOnce(&Database) -> Result<T>,
    ) -> Result<T> {
        self.commit_outcome(
            records,
            register_dependencies,
            complete,
            JobOutcome::Succeeded,
            None,
        )
    }

    /// A complete diagnostic generation may record a scientifically rejected QC
    /// result while the attempt itself fails. Publication remains one transaction.
    pub fn commit_outcome<T>(
        &mut self,
        records: &[ArtifactRecord],
        register_dependencies: impl FnOnce(&Database) -> Result<()>,
        complete: impl FnOnce(&Database) -> Result<T>,
        outcome: JobOutcome,
        message: Option<&str>,
    ) -> Result<T> {
        if matches!(outcome, JobOutcome::Interrupted) {
            bail!("interrupted work cannot publish a generation");
        }
        if records.is_empty() {
            bail!("generation must register its output artifacts");
        }
        let seal = self
            .seal
            .as_ref()
            .context("generation must be sealed before publication")?;
        if inventory(&self.head.payload_path, false)? != seal.files {
            bail!("generation changed after sealing");
        }
        for input in &seal.inputs {
            if !input.matches_file()? {
                bail!("generation input changed after sealing: {}", input.path);
            }
        }
        for record in records {
            if record.project_id != self.db.project_id()
                || record.run_accession.as_deref()
                    != match self.head.output_group {
                        OutputGroup::Download | OutputGroup::Preprocess | OutputGroup::Quant => {
                            Some(self.head.subject.as_str())
                        }
                        OutputGroup::Refs | OutputGroup::De | OutputGroup::Corr => None,
                    }
                || !seal.files.iter().any(|file| {
                    file.path == record.path
                        && file.checksum_type == record.checksum_type
                        && file.checksum == record.checksum
                        && file.bytes == record.bytes
                })
            {
                bail!(
                    "artifact registration does not match the sealed payload: {}",
                    record.path
                );
            }
        }
        let result = self.db.atomic(|db| {
            if let Some(request) = &self.scientific_request
                && db.store_scientific_request(request)? != self.head.request_sha256 {
                bail!("scientific request changed after sealing");
            }
            db.verify_dependency_catalog(&seal.dependencies)?;
            db.retire_generation_artifacts(self.head.output_group, &self.head.subject)?;
            register_dependencies(db)?;
            for record in records { db.record_artifact(record)?; }
            db.write(|connection| {
                let changed = connection.execute("UPDATE artifact_generations SET state = 'COMMITTED', artifacts_json = ?1, committed_at = ?2
                    WHERE generation_id = ?3 AND state = 'SEALED' AND inventory_json = ?4", params![serde_json::to_string(records)?, now_rfc3339(), self.head.generation_id, serde_json::to_string(seal)?])?;
                if changed != 1 { bail!("sealed generation no longer matches its inventory"); }
                for file in &seal.files {
                    connection.execute("INSERT INTO generation_files(generation_id,path,sha256,bytes) VALUES (?1,?2,?3,?4)", params![self.head.generation_id,file.path,file.checksum,i64::try_from(file.bytes)?])?;
                }
                for dependency in &seal.dependencies {
                    connection.execute("INSERT INTO generation_dependencies(consumer_generation_id,producer_generation_id,producer_seal_sha256,role,path,artifact_json) VALUES (?1,?2,?3,?4,?5,?6)",
                        params![self.head.generation_id,dependency.producer_generation_id,dependency.producer_seal_sha256,dependency.input.kind,dependency.input.path,serde_json::to_string(&dependency.input)?])?;
                }
                connection.execute("INSERT INTO generation_heads(project_id, output_group, subject, generation_id) VALUES (?1, ?2, ?3, ?4)
                    ON CONFLICT(project_id, output_group, subject) DO UPDATE SET generation_id = excluded.generation_id",
                    params![db.project_id(), self.head.output_group.as_str(), self.head.subject, self.head.generation_id])?;
                connection.execute("DELETE FROM generation_binding_heads WHERE project_id=?1 AND output_group=?2 AND subject=?3",
                    params![db.project_id(), self.head.output_group.as_str(), self.head.subject])?;
                Ok(())
            })?;
            let value = complete(db)?;
            db.complete_owned_job(outcome, message)?;
            Ok(value)
        })?;
        self.committed = true;
        Ok(result)
    }
}

impl Drop for Generation {
    fn drop(&mut self) {
        if !self.committed {
            let _ = self.db.write(|connection| {
                connection.execute("UPDATE artifact_generations SET state = 'ABANDONED' WHERE generation_id = ?1 AND state IN ('BUILDING','SEALED')", params![self.head.generation_id])?;
                Ok(())
            });
        }
    }
}

fn valid_sha256(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

fn same_file(left: &ManifestArtifact, right: &ManifestArtifact) -> bool {
    left.path == right.path
        && left.checksum_type == "sha256"
        && left.checksum_type == right.checksum_type
        && left.checksum == right.checksum
        && left.bytes == right.bytes
}

fn is_trimmed_fastq(kind: &str) -> bool {
    matches!(
        kind,
        "FASTQ_TRIMMED_R1" | "FASTQ_TRIMMED_R2" | "FASTQ_TRIMMED_SINGLE"
    )
}

fn prunable_fastq(group: OutputGroup, kind: &str) -> bool {
    match group {
        OutputGroup::Download => matches!(kind, "FASTQ_R1" | "FASTQ_R2" | "FASTQ_SINGLE"),
        OutputGroup::Preprocess => is_trimmed_fastq(kind),
        _ => false,
    }
}

fn quant_manifest_digest(seal: &Seal) -> Option<&str> {
    seal.files
        .iter()
        .find(|file| {
            Path::new(&file.path)
                .file_name()
                .is_some_and(|name| name == "rnaa_quant_manifest.json")
        })
        .map(|file| file.checksum.as_str())
}

fn inventory(root: &Path, flush: bool) -> Result<Vec<ManifestArtifact>> {
    let mut files = Vec::new();
    let mut directories = Vec::new();
    for entry in walkdir::WalkDir::new(root).follow_links(false) {
        let entry = entry?;
        if entry.file_type().is_symlink() {
            bail!(
                "generation cannot contain symlinks: {}",
                entry.path().display()
            );
        }
        if !entry.file_type().is_file() && !entry.file_type().is_dir() {
            bail!("generation can contain only regular files and directories");
        }
        if entry.file_type().is_file() {
            if flush {
                fs::File::open(entry.path())?.sync_all()?;
            }
            files.push(ManifestArtifact::snapshot("GENERATION_FILE", entry.path())?);
        } else if entry.file_type().is_dir() {
            directories.push(entry.path().to_path_buf());
        } else {
            bail!(
                "unsupported generation file type: {}",
                entry.path().display()
            );
        }
    }
    if flush {
        for directory in directories.into_iter().rev() {
            fs::File::open(directory)?.sync_all()?;
        }
    }
    files.sort_by(|a, b| a.path.cmp(&b.path));
    Ok(files)
}

#[cfg(test)]
mod tests;
