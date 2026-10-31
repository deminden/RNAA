//! Reversible cleanup of exclusively project-owned, count-backed input files.
mod files;
mod io;
mod proof;
mod restore;
#[cfg(test)]
mod tests;
use super::{OutputGroup, Seal, canonical_missing_path, quant_manifest_digest};
use crate::{
    ArtifactKind, ArtifactRecord, CleanupMode, CleanupWhen, Database, LibraryLayout,
    ManifestArtifact, ProjectConfig,
    util::{compute_sha256, write_json_pretty},
};
use anyhow::{Context, Result, bail};
use rusqlite::params;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};
fn digest(bytes: impl AsRef<[u8]>) -> String {
    hex::encode(Sha256::digest(bytes.as_ref()))
}
fn now_ms() -> Result<i64> {
    Ok(i64::try_from(
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_millis(),
    )?)
}
pub(super) fn fastq(kind: ArtifactKind) -> bool {
    matches!(
        kind,
        ArtifactKind::FastqR1
            | ArtifactKind::FastqR2
            | ArtifactKind::FastqSingle
            | ArtifactKind::FastqTrimmedR1
            | ArtifactKind::FastqTrimmedR2
            | ArtifactKind::FastqTrimmedSingle
    )
}
fn selected(mode: CleanupMode, kind: ArtifactKind) -> bool {
    match mode {
        CleanupMode::None => false,
        CleanupMode::Fastq => fastq(kind),
        CleanupMode::Sra => kind == ArtifactKind::Sra,
        CleanupMode::All => fastq(kind) || kind == ArtifactKind::Sra,
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CleanupPolicy {
    pub mode: CleanupMode,
    pub when: CleanupWhen,
    pub trash_days: u64,
}
impl CleanupPolicy {
    pub fn from_config(config: &ProjectConfig) -> Self {
        Self {
            mode: config.quant.cleanup,
            when: config.quant.cleanup_on,
            trash_days: config.quant.trash_days,
        }
    }
}
#[derive(Clone)]
pub struct FileCleanupOptions {
    pub mode: CleanupMode,
    pub run: Option<String>,
    pub apply: bool,
    pub purge_expired: bool,
    pub trash_days: u64,
    pub expected_policy: Option<CleanupPolicy>,
}
#[derive(Debug, Serialize)]
pub struct FileCleanupEntry {
    pub trash_id: Option<String>,
    pub operation_id: Option<String>,
    pub run: String,
    pub path: String,
    pub state: String,
    pub action: String,
    pub bytes: u64,
    pub reason: Option<String>,
}
#[derive(Debug, Serialize)]
pub struct FileCleanupReport {
    pub applied: bool,
    pub busy: bool,
    pub entries: Vec<FileCleanupEntry>,
    pub errors: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Descriptor {
    schema_version: u32,
    trash_id: String,
    operation_id: String,
    project_id: String,
    run: String,
    layout: LibraryLayout,
    source_selection: String,
    record: ArtifactRecord,
    source: ManifestArtifact,
    canonical_source: String,
    producer: Option<String>,
    producer_seal: Option<String>,
    consumer: String,
    consumer_seal: String,
    consumer_manifest: String,
    binding: Option<String>,
    binding_sha256: Option<String>,
    mode: CleanupMode,
    created_at_ms: i64,
    purge_after_ms: i64,
}
pub(super) struct GenerationLink {
    pub(super) producer: Option<String>,
    pub(super) consumer: String,
    pub(super) restorable: bool,
    pub(super) absence: bool,
}
#[derive(Clone)]
pub(super) struct Stored {
    pub descriptor: Descriptor,
    pub state: String,
}
impl Stored {
    pub fn active_absence(&self) -> bool {
        matches!(
            self.state.as_str(),
            "MOVING" | "TRASHED" | "PURGING" | "PURGED"
        )
    }
    fn entry(&self, action: &str) -> FileCleanupEntry {
        FileCleanupEntry {
            trash_id: Some(self.descriptor.trash_id.clone()),
            operation_id: Some(self.descriptor.operation_id.clone()),
            run: self.descriptor.run.clone(),
            path: self.descriptor.record.path.clone(),
            state: self.state.clone(),
            action: action.into(),
            bytes: self.descriptor.source.bytes,
            reason: None,
        }
    }
}

impl Database {
    /// Finish only previously authorized moves; no new cleanup or expiration.
    /// Returns false when live readers prevent recovering a pending operation.
    pub fn recover_input_cleanup(&self) -> Result<bool> {
        let pending: bool = self.read(|conn| Ok(conn.query_row("SELECT EXISTS(SELECT 1 FROM file_trash WHERE project_id=?1 AND state IN ('MOVING','RESTORING','PURGING'))", [self.project_id()], |row| row.get(0))?))?;
        if !pending {
            return Ok(true);
        }
        let Some(_lock) = super::pins::maintenance_lock(self)? else {
            return Ok(false);
        };
        let mut db = self.clone();
        db.generation_maintenance = true;
        db.recover_file_trash()?;
        Ok(true)
    }
    pub(super) fn file_trash_rows(&self) -> Result<Vec<Stored>> {
        self.file_trash_rows_filtered(None, None)
    }
    fn file_trash_rows_filtered(
        &self,
        run: Option<&str>,
        source: Option<&str>,
    ) -> Result<Vec<Stored>> {
        let rows=self.read(|conn| {
            let run_filter = if run.is_some() { "run_accession=?2" } else { "?2 IS NULL" };
            let source_filter = if source.is_some() { "canonical_source=?3" } else { "?3 IS NULL" };
            let sql = format!("SELECT trash_id,run_accession,canonical_source,operation_id,state,descriptor_json,descriptor_sha256,created_at_ms,purge_after_ms FROM file_trash WHERE project_id=?1 AND {run_filter} AND {source_filter} ORDER BY created_at_ms,trash_id");
            let mut query=conn.prepare(&sql)?;
            Ok(query.query_map(params![self.project_id(),run,source],|row|Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,String>(2)?,row.get::<_,String>(3)?,row.get::<_,String>(4)?,row.get::<_,String>(5)?,row.get::<_,String>(6)?,row.get::<_,i64>(7)?,row.get::<_,i64>(8)?)))?.collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        rows.into_iter()
            .map(
                |(id, run, source, operation, state, json, hash, created, expiry)| {
                    let descriptor: Descriptor = serde_json::from_str(&json)?;
                    if serde_json::to_string(&descriptor)? != json
                        || digest(json.as_bytes()) != hash
                        || descriptor.trash_id != id
                        || descriptor.run != run
                        || descriptor.canonical_source != source
                        || descriptor.operation_id != operation
                        || descriptor.created_at_ms != created
                        || descriptor.purge_after_ms != expiry
                    {
                        bail!("file trash catalog differs from its immutable descriptor");
                    }
                    self.verify_trash_descriptor(&descriptor)?;
                    let witness = self.file_receipt_path(&id, false)?;
                    let bytes = fs::read(witness)?;
                    if serde_json::to_value(&descriptor)?
                        != serde_json::from_slice::<serde_json::Value>(&bytes)?
                    {
                        bail!("file trash receipt witness changed");
                    }
                    Ok(Stored { descriptor, state })
                },
            )
            .collect()
    }
    pub fn file_trash_status(&self) -> Result<serde_json::Value> {
        let _pin = super::pins::read_lock(self.path())?;
        let entries=self.file_trash_rows()?.into_iter().map(|row| {
            serde_json::json!({"trash_id":row.descriptor.trash_id,"operation_id":row.descriptor.operation_id,
                "run":row.descriptor.run,"kind":row.descriptor.record.kind,"source":row.descriptor.canonical_source,
                "state":row.state,"bytes":row.descriptor.source.bytes,"producer_generation":row.descriptor.producer,
                "consumer_generation":row.descriptor.consumer,"created_at_ms":row.descriptor.created_at_ms,"purge_after_ms":row.descriptor.purge_after_ms})
        }).collect::<Vec<_>>();
        Ok(serde_json::json!({"schema_version":1,"entries":entries}))
    }
    fn verify_cleanup_policy(&self, expected: Option<&CleanupPolicy>) -> Result<()> {
        if let Some(expected) = expected {
            let catalog: ProjectConfig = serde_json::from_str(&self.project()?.config_json)?;
            let disk = ProjectConfig::load(
                &self
                    .path()
                    .parent()
                    .context("project has no root")?
                    .join("rnaa.toml"),
            )?;
            if CleanupPolicy::from_config(&catalog) != *expected
                || CleanupPolicy::from_config(&disk) != *expected
            {
                bail!("cleanup policy changed during file cleanup; retry");
            }
        }
        Ok(())
    }
    pub fn cleanup_files(&self, options: FileCleanupOptions) -> Result<FileCleanupReport> {
        let mut report = FileCleanupReport {
            applied: options.apply,
            busy: false,
            entries: Vec::new(),
            errors: Vec::new(),
        };
        let Some(_lock) = super::pins::maintenance_lock(self)? else {
            report.busy = true;
            return Ok(report);
        };
        let mut db = self.clone();
        db.generation_maintenance = true;
        if options.apply {
            db.recover_file_trash()?;
        }
        db.verify_cleanup_policy(options.expected_policy.as_ref())?;
        let prior = db.file_trash_rows()?;
        let artifacts = db.list_artifacts()?;
        let mut multiplicity = BTreeMap::<String, usize>::new();
        for artifact in &artifacts {
            let path = canonical_missing_path(Path::new(&artifact.path))?
                .display()
                .to_string();
            *multiplicity.entry(path).or_default() += 1;
        }
        let mut operations = BTreeMap::new();
        let mut prepared = Vec::new();
        for artifact in artifacts {
            let Some(run) = artifact.run_accession.clone() else {
                continue;
            };
            if !selected(options.mode, artifact.kind)
                || options.run.as_ref().is_some_and(|wanted| *wanted != run)
            {
                continue;
            }
            let operation = operations
                .entry(run.clone())
                .or_insert_with(|| uuid::Uuid::new_v4().to_string());
            match db.prepare_trash_descriptor(
                &artifact,
                operation,
                options.mode,
                options.trash_days,
                &multiplicity,
            ) {
                Ok(descriptor) => {
                    report.entries.push(
                        Stored {
                            descriptor: descriptor.clone(),
                            state: "LIVE".into(),
                        }
                        .entry("trash"),
                    );
                    prepared.push(descriptor);
                }
                Err(error) => {
                    report.errors.push(format!("{}: {error:#}", artifact.path));
                    report.entries.push(FileCleanupEntry {
                        trash_id: None,
                        operation_id: None,
                        run,
                        path: artifact.path,
                        state: "LIVE".into(),
                        action: "held".into(),
                        bytes: artifact.bytes,
                        reason: Some(format!("{error:#}")),
                    });
                }
            }
        }
        // An operation is a complete run group; one rejected member holds all.
        for (run, reason) in db.incomplete_cleanup_groups(&prepared)? {
            report.errors.push(format!("{run}: {reason}"));
            for entry in &mut report.entries {
                if entry.run == run {
                    entry.action = "held".into();
                    entry.reason = Some(reason.clone());
                }
            }
        }
        let held = report
            .entries
            .iter()
            .filter(|entry| entry.action == "held")
            .map(|entry| entry.run.clone())
            .collect::<std::collections::BTreeSet<_>>();
        prepared.retain(|descriptor| !held.contains(&descriptor.run));
        for entry in &mut report.entries {
            if entry.action == "trash" && held.contains(&entry.run) {
                entry.action = "held_run_group".into();
            }
        }
        let timestamp = now_ms()?;
        let unexpired_groups = prior
            .iter()
            .filter(|row| row.state != "TRASHED" || row.descriptor.purge_after_ms > timestamp)
            .map(|row| row.descriptor.operation_id.clone())
            .collect::<std::collections::BTreeSet<_>>();
        let expired = prior
            .iter()
            .filter(|row| {
                row.state == "TRASHED"
                    && !unexpired_groups.contains(&row.descriptor.operation_id)
                    && options
                        .run
                        .as_ref()
                        .is_none_or(|run| *run == row.descriptor.run)
            })
            .cloned()
            .collect::<Vec<_>>();
        for row in &prior {
            let purge = options.purge_expired
                && expired
                    .iter()
                    .any(|item| item.descriptor.trash_id == row.descriptor.trash_id);
            report
                .entries
                .push(row.entry(if purge { "purge" } else { "retained" }));
        }
        if options.apply {
            db.verify_cleanup_policy(options.expected_policy.as_ref())?;
            db.mark_file_trash(&prepared)?;
            db.recover_file_trash()?;
            if options.purge_expired {
                db.verify_cleanup_policy(options.expected_policy.as_ref())?;
                db.mark_file_purge(&expired)?;
                db.recover_file_trash()?;
            }
            let states = db
                .file_trash_rows()?
                .into_iter()
                .map(|row| (row.descriptor.trash_id, row.state))
                .collect::<BTreeMap<_, _>>();
            for entry in &mut report.entries {
                if let Some(id) = &entry.trash_id
                    && let Some(state) = states.get(id)
                {
                    entry.state = state.clone();
                }
            }
        }
        Ok(report)
    }
}
