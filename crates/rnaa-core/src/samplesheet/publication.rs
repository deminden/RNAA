//! Immutable metadata witnesses and a transactional head with replayable view installation.
use super::{ColumnMap, render_samplesheet, table::Table};
use crate::util::{now_rfc3339, publish_file, write_json_pretty};
use crate::{ArtifactRecord, Database, ManifestArtifact, ProjectPaths, RunRecord};
use anyhow::{Context, Result, bail};
use rusqlite::{OptionalExtension, params};
use serde::{Deserialize, Serialize};
use std::{fs, path::Path};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ConflictPolicy {
    #[default]
    Fail,
    KeepLocal,
    UseSource,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Publication {
    pub schema_version: u32,
    pub publication_id: String,
    pub project_id: String,
    pub parent: Option<String>,
    pub policy: ConflictPolicy,
    pub generated: ManifestArtifact,
    pub history_generated: ManifestArtifact,
    pub history_merged: ManifestArtifact,
    pub local: Option<ManifestArtifact>,
    pub merged: ManifestArtifact,
    pub column_map: ManifestArtifact,
    pub sources: Vec<ManifestArtifact>,
    pub source_observations: Vec<ManifestArtifact>,
}
impl Publication {
    pub fn verify(&self) -> Result<()> {
        if self.schema_version != 1 || uuid::Uuid::parse_str(&self.publication_id).is_err() {
            bail!("invalid metadata publication identity");
        }
        for input in [
            &self.generated,
            &self.merged,
            &self.column_map,
            &self.history_generated,
            &self.history_merged,
        ]
        .into_iter()
        .chain(self.local.iter())
        .chain(&self.sources)
        {
            if !input.matches_file()? {
                bail!(
                    "metadata publication witness is missing or corrupt: {}",
                    input.path
                );
            }
        }
        let receipt = Path::new(&self.generated.path)
            .parent()
            .context("metadata witness has no directory")?
            .join("receipt.json");
        if serde_json::from_slice::<Self>(&fs::read(receipt)?)? != *self {
            bail!("metadata receipt differs from its committed catalog descriptor");
        }
        if self.sources.len() != self.source_observations.len() + 1
            || self
                .sources
                .iter()
                .skip(1)
                .zip(&self.source_observations)
                .any(|(copy, original)| {
                    copy.checksum != original.checksum
                        || copy.bytes != original.bytes
                        || copy.kind != original.kind
                })
        {
            bail!("metadata source witness inventory is inconsistent");
        }
        Ok(())
    }
}

/// Prepared on disk before publication; failed/superseded candidates remain inspectable.
pub struct PreparedPublication {
    pub publication: Publication,
    observed_local: Option<ManifestArtifact>,
    observed_map: ManifestArtifact,
    observed_sources: Vec<ManifestArtifact>,
}

impl PreparedPublication {
    pub fn prepare(
        db: &Database,
        paths: &ProjectPaths,
        runs: &[RunRecord],
        artifacts: &[ArtifactRecord],
        map: &ColumnMap,
        policy: ConflictPolicy,
    ) -> Result<Self> {
        Self::prepare_with_sources(db, paths, runs, artifacts, map, policy, &[])
    }

    pub fn prepare_with_sources(
        db: &Database,
        paths: &ProjectPaths,
        runs: &[RunRecord],
        artifacts: &[ArtifactRecord],
        map: &ColumnMap,
        policy: ConflictPolicy,
        source_observations: &[ManifestArtifact],
    ) -> Result<Self> {
        recover_view(db, paths, policy)?;
        let previous = db.metadata_publication()?;
        let observed_local = paths
            .samplesheet_path()
            .exists()
            .then(|| ManifestArtifact::snapshot("EDITABLE_METADATA", &paths.samplesheet_path()))
            .transpose()?;
        let observed_map = ManifestArtifact::snapshot("COLUMN_MAP", &paths.column_map_path())?;
        let id = uuid::Uuid::new_v4().to_string();
        let dir = paths.metadata_dir.join("revisions").join(&id);
        fs::create_dir_all(&dir)?;
        let generated_bytes = render_samplesheet(db.project_id(), runs, artifacts, map)?.0;
        let generated = snapshot_bytes(
            "GENERATED_METADATA",
            &dir.join("generated.tsv"),
            &generated_bytes,
        )?;
        let mut sources = vec![snapshot_bytes(
            "RESOLVED_RUNS",
            &dir.join("runs.json"),
            &serde_json::to_vec_pretty(runs)?,
        )?];
        for (i, source) in source_observations.iter().enumerate() {
            let copy = snapshot_bytes(
                &source.kind,
                &dir.join(format!("source-{i}.tsv")),
                &fs::read(&source.path)?,
            )?;
            if source.checksum != copy.checksum
                || source.bytes != copy.bytes
                || !source.matches_file()?
            {
                bail!("metadata override changed while freezing its source witness; retry");
            }
            sources.push(copy);
        }
        let local = observed_local
            .as_ref()
            .map(|source| {
                snapshot_bytes(
                    "EDITED_METADATA",
                    &dir.join("local.tsv"),
                    &fs::read(&source.path)?,
                )
            })
            .transpose()?;
        let column_map = snapshot_bytes(
            "COLUMN_MAP",
            &dir.join("column_map.tsv"),
            &fs::read(&observed_map.path)?,
        )?;
        if local.as_ref().map(|p| (&p.checksum, p.bytes))
            != observed_local.as_ref().map(|p| (&p.checksum, p.bytes))
            || column_map.checksum != observed_map.checksum
            || column_map.bytes != observed_map.bytes
        {
            bail!("metadata sources changed while their witnesses were copied; retry");
        }
        let remote = Table::parse(&generated_bytes)?;
        let baseline = match &previous {
            Some(previous) => {
                previous.verify()?;
                Table::parse(&fs::read(&previous.history_generated.path)?)?
            }
            None => remote.clone(),
        };
        let editable = match &local {
            Some(local) => Table::parse(&fs::read(&local.path)?)?,
            None => remote.clone(),
        };
        let local_table = match &previous {
            Some(previous) => {
                Table::parse(&fs::read(&previous.history_merged.path)?)?.overlay(&editable)
            }
            None => editable,
        };
        let all_source = baseline.overlay(&remote);
        let history_bytes = Table::merge(&baseline, &local_table, &all_source, policy)
            .with_context(|| format!("refresh candidate versions retained at {}", dir.display()))?;
        let history_merged = snapshot_bytes(
            "MERGED_METADATA_HISTORY",
            &dir.join("history_merged.tsv"),
            &history_bytes,
        )?;
        let history_generated = snapshot_bytes(
            "GENERATED_METADATA_HISTORY",
            &dir.join("history_generated.tsv"),
            &all_source.project(&all_source)?,
        )?;
        let bytes = Table::parse(&history_bytes)?.project(&remote)?;
        let merged = snapshot_bytes("PUBLISHED_METADATA", &dir.join("merged.tsv"), &bytes)?;
        let publication = Publication {
            schema_version: 1,
            publication_id: id,
            project_id: db.project_id().into(),
            parent: previous.map(|p| p.publication_id),
            policy,
            generated,
            history_generated,
            history_merged,
            local,
            merged,
            column_map,
            sources,
            source_observations: source_observations.to_vec(),
        };
        write_json_pretty(&dir.join("receipt.json"), &publication)?;
        let prepared = Self {
            publication,
            observed_local,
            observed_map,
            observed_sources: source_observations.to_vec(),
        };
        prepared.verify_sources(paths)?;
        prepared.publication.verify()?;
        Ok(prepared)
    }

    pub fn verify_sources(&self, paths: &ProjectPaths) -> Result<()> {
        if !matches_observation(&self.observed_local, &paths.samplesheet_path())?
            || !self.observed_map.matches_file()?
            || self
                .observed_sources
                .iter()
                .map(ManifestArtifact::matches_file)
                .collect::<Result<Vec<_>>>()?
                .contains(&false)
        {
            bail!(
                "editable metadata or column map changed during refresh; retry without overwriting the edit"
            );
        }
        Ok(())
    }

    /// May run inside the same database transaction as resolved cohort publication.
    pub fn commit_catalog(&self, db: &Database, paths: &ProjectPaths) -> Result<()> {
        self.verify_sources(paths)?;
        self.publication.verify()?;
        self.commit_catalog_preverified(db)
    }
    /// Call after verifying witnesses/sources outside the writer transaction.
    /// Immutable payloads are already prepared; only the catalog head is checked here.
    pub fn commit_catalog_preverified(&self, db: &Database) -> Result<()> {
        db.atomic_or_current(|db| self.commit_in_transaction(db))
    }
    fn commit_in_transaction(&self, db: &Database) -> Result<()> {
        if db.metadata_publication()?.map(|p| p.publication_id) != self.publication.parent {
            bail!("metadata refresh was superseded; retry");
        }
        db.write(|conn| {
                conn.execute("INSERT INTO metadata_publications(publication_id, project_id, descriptor_json, created_at) VALUES (?1, ?2, ?3, ?4)",
                    params![self.publication.publication_id, db.project_id(), serde_json::to_string(&self.publication)?, now_rfc3339()])?;
                conn.execute("INSERT INTO metadata_heads(project_id, publication_id) VALUES (?1, ?2) ON CONFLICT(project_id) DO UPDATE SET publication_id = excluded.publication_id",
                    params![db.project_id(), self.publication.publication_id])?;
                Ok(())
            })?;
        db.append_event("metadata", None, "published metadata revision", serde_json::json!({"publication_id":self.publication.publication_id,"parent":self.publication.parent,"policy":self.publication.policy}))?;
        Ok(())
    }
}

impl Database {
    pub fn metadata_publication(&self) -> Result<Option<Publication>> {
        self.read(|conn| {
            let row: Option<(String, String)> = conn.query_row(
                "SELECT p.publication_id, p.descriptor_json FROM metadata_heads h JOIN metadata_publications p ON p.publication_id = h.publication_id WHERE h.project_id = ?1",
                [self.project_id()], |row| Ok((row.get(0)?, row.get(1)?))).optional()?;
            row.map(|(id, json)| {
                let publication: Publication = serde_json::from_str(&json)?;
                if publication.project_id != self.project_id() || publication.publication_id != id {
                    bail!("metadata publication belongs to another project or identity");
                }
                Ok(publication)
            }).transpose()
        })
    }
}

/// Replay a committed but uninstalled view. A later editor's changes are merged
/// against the exact pre-publication view, never mistaken for source-side edits.
pub fn recover_view(db: &Database, paths: &ProjectPaths, policy: ConflictPolicy) -> Result<()> {
    let Some(head) = db.metadata_publication()? else {
        return Ok(());
    };
    let installed: Option<String> = db.read(|conn| {
        Ok(conn.query_row(
            "SELECT installed_id FROM metadata_heads WHERE project_id = ?1",
            [db.project_id()],
            |row| row.get(0),
        )?)
    })?;
    if installed.as_deref() == Some(&head.publication_id) && paths.samplesheet_path().exists() {
        return Ok(());
    }
    head.verify()?;
    let path = paths.samplesheet_path();
    let expected = fs::read(&head.merged.path)?;
    let current = path
        .exists()
        .then(|| ManifestArtifact::snapshot("EDITABLE_METADATA", &path))
        .transpose()?;
    let bytes = if current
        .as_ref()
        .is_some_and(|p| p.checksum == head.merged.checksum && p.bytes == head.merged.bytes)
        || matches_observation(&head.local, &path)?
    {
        expected
    } else if let Some(current) = &current {
        let remote = Table::parse(&expected)?;
        let base = match &head.local {
            Some(previous) => Table::parse(&fs::read(&previous.path)?)?,
            None => remote.clone(),
        };
        Table::merge(&base, &Table::parse(&fs::read(&current.path)?)?, &remote, policy)
            .context("pending metadata publication conflicts with a later edit; retained receipt records both original versions")?
    } else {
        expected
    };
    if !matches_observation(&current, &path)? {
        bail!("metadata changed while recovering its view; retry");
    }
    publish_file(&path, |staged| {
        fs::write(staged, &bytes)?;
        Ok(())
    })?;
    db.write(|conn| {
        let changed = conn.execute("UPDATE metadata_heads SET installed_id = ?1 WHERE project_id = ?2 AND publication_id = ?1", params![head.publication_id, db.project_id()])?;
        if changed != 1 { bail!("metadata publication was superseded during view recovery"); }
        Ok(())
    })
}

fn matches_observation(observed: &Option<ManifestArtifact>, path: &Path) -> Result<bool> {
    match observed {
        Some(previous) if path.exists() => {
            Ok(ManifestArtifact::snapshot("VIEW", path)?.checksum == previous.checksum)
        }
        None => Ok(!path.exists()),
        _ => Ok(false),
    }
}
fn snapshot_bytes(kind: &str, path: &Path, bytes: &[u8]) -> Result<ManifestArtifact> {
    publish_file(path, |staged| {
        fs::write(staged, bytes)?;
        Ok(())
    })?;
    ManifestArtifact::snapshot(kind, path)
}

#[cfg(test)]
mod tests;
