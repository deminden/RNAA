//! Artifact, shared-blob and retained-input catalog access.
use super::Database;
use crate::model::{ArtifactRecord, SharedBlobRecord};
use crate::util::now_rfc3339;
use anyhow::{Context, Result};
use rusqlite::{OptionalExtension, Row, params};
use std::str::FromStr;

impl Database {
    pub fn record_artifact(&self, artifact: &ArtifactRecord) -> Result<()> {
        let canonical =
            crate::generations::canonical_missing_path(std::path::Path::new(&artifact.path))?;
        self.write(|conn| {
            let retired:bool=conn.query_row("SELECT EXISTS(SELECT 1 FROM generation_files f JOIN generation_retirements r ON r.generation_id=f.generation_id WHERE f.path=?1)",[canonical.to_string_lossy().as_ref()],|row|row.get(0))?;
            if retired {anyhow::bail!("restore a retired generation before registering its artifacts");}
            conn.execute(
                "INSERT OR REPLACE INTO artifacts (project_id, run_accession, kind, path, blob_id, shared_path, checksum_type, checksum, bytes, created_at) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)",
                params![
                    artifact.project_id,
                    artifact.run_accession,
                    artifact.kind.as_str(),
                    artifact.path,
                    artifact.blob_id,
                    artifact.shared_path,
                    artifact.checksum_type,
                    artifact.checksum,
                    artifact.bytes as i64,
                    artifact.created_at
                ],
            )
            .with_context(|| format!("failed to record artifact {}", artifact.path))?;
            if let Some(blob_id) = artifact.blob_id.as_deref() {
                conn.execute(
                    "INSERT OR REPLACE INTO shared_blob_refs (blob_id, project_id, artifact_path, referenced_at) VALUES (?1, ?2, ?3, ?4)",
                    params![blob_id, artifact.project_id, artifact.path, now_rfc3339()],
                )
                .with_context(|| {
                    format!(
                        "failed to link artifact {} to shared blob {blob_id}",
                        artifact.path
                    )
                })?;
            }
            Ok(())
        })
    }

    pub fn record_shared_blob(&self, blob: &SharedBlobRecord) -> Result<()> {
        self.write(|conn| {
            conn.execute(
                "INSERT OR REPLACE INTO shared_blobs (blob_id, storage_path, checksum_type, checksum, bytes, created_at) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
                params![
                    blob.blob_id,
                    blob.storage_path,
                    blob.checksum_type,
                    blob.checksum,
                    blob.bytes as i64,
                    blob.created_at
                ],
            )
            .with_context(|| format!("failed to record shared blob {}", blob.blob_id))?;
            Ok(())
        })
    }

    pub fn get_shared_blob(&self, blob_id: &str) -> Result<Option<SharedBlobRecord>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT blob_id, storage_path, checksum_type, checksum, bytes, created_at FROM shared_blobs WHERE blob_id = ?1",
            )?;
            stmt.query_row(params![blob_id], parse_shared_blob)
                .optional()
                .context("failed to load shared blob")
        })
    }

    pub fn list_artifacts(&self) -> Result<Vec<ArtifactRecord>> {
        self.list_artifacts_for_run(None)
    }

    pub fn list_artifacts_for_run(
        &self,
        run_accession: Option<&str>,
    ) -> Result<Vec<ArtifactRecord>> {
        self.read(|conn| {
            let sql = if run_accession.is_some() {
                "SELECT project_id, run_accession, kind, path, blob_id, shared_path, checksum_type, checksum, bytes, created_at FROM artifacts WHERE run_accession = ?1 ORDER BY kind, path"
            } else {
                "SELECT project_id, run_accession, kind, path, blob_id, shared_path, checksum_type, checksum, bytes, created_at FROM artifacts ORDER BY kind, path"
            };
            let mut stmt = conn.prepare(sql)?;
            let rows = if let Some(run_accession) = run_accession {
                stmt.query_map(params![run_accession], parse_artifact)?
            } else {
                stmt.query_map([], parse_artifact)?
            };
            let mut artifacts = Vec::new();
            for row in rows {
                artifacts.push(row?);
            }
            Ok(artifacts)
        })
    }

    pub fn delete_artifact(&self, path: &str) -> Result<()> {
        self.write(|conn| {
            conn.execute("DELETE FROM artifacts WHERE path = ?1", params![path])
                .with_context(|| format!("failed to delete artifact record {path}"))?;
            Ok(())
        })
    }

    /// Authorize absence before unlinking. Existing files must still verify against
    /// their snapshot, so a failed unlink cannot hide subsequent modifications.
    pub fn record_retained_input(
        &self,
        run: &str,
        generation: &str,
        artifact: &crate::manifest::ManifestArtifact,
    ) -> Result<()> {
        self.write(|conn| {
            conn.execute("INSERT OR REPLACE INTO retained_inputs (project_id, run_accession, generation_sha256, path, artifact_json, recorded_at) VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
                params![self.project_id, run, generation, artifact.path, serde_json::to_string(artifact)?, now_rfc3339()])?;
            Ok(())
        })
    }

    pub fn retained_inputs(&self, run: &str) -> Result<Vec<crate::fingerprint::RetainedInput>> {
        let mut evidence: Vec<crate::fingerprint::RetainedInput> = self.read(|conn| {
            let mut stmt = conn.prepare("SELECT generation_sha256, artifact_json FROM retained_inputs WHERE project_id = ?1 AND run_accession = ?2")?;
            let rows = stmt.query_map(params![self.project_id, run], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
            })?;
            rows.map(|row| {
                let (generation_sha256, json) = row?;
                Ok(crate::fingerprint::RetainedInput {
                    generation_sha256,
                    artifact: serde_json::from_str(&json)?,
                })
            })
            .collect()
        })?;
        for item in self.file_trash_retained_inputs(run)? {
            if !evidence.iter().any(|prior| {
                prior.generation_sha256 == item.generation_sha256 && prior.artifact == item.artifact
            }) {
                evidence.push(item);
            }
        }
        Ok(evidence)
    }

    pub fn disk_usage_bytes(&self) -> Result<u64> {
        self.read(|conn| {
            let total: Option<i64> =
                conn.query_row("SELECT SUM(bytes) FROM artifacts", [], |row| row.get(0))?;
            Ok(total.unwrap_or_default() as u64)
        })
    }
}

fn parse_artifact(row: &Row<'_>) -> rusqlite::Result<ArtifactRecord> {
    let kind_text: String = row.get(2)?;
    let bytes: i64 = row.get(8)?;
    let kind = crate::model::ArtifactKind::from_str(&kind_text).map_err(|err| {
        rusqlite::Error::FromSqlConversionFailure(
            2,
            rusqlite::types::Type::Text,
            Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                err.to_string(),
            )),
        )
    })?;
    Ok(ArtifactRecord {
        project_id: row.get(0)?,
        run_accession: row.get(1)?,
        kind,
        path: row.get(3)?,
        blob_id: row.get(4)?,
        shared_path: row.get(5)?,
        checksum_type: row.get(6)?,
        checksum: row.get(7)?,
        bytes: bytes as u64,
        created_at: row.get(9)?,
    })
}

fn parse_shared_blob(row: &Row<'_>) -> rusqlite::Result<SharedBlobRecord> {
    let bytes: i64 = row.get(4)?;
    Ok(SharedBlobRecord {
        blob_id: row.get(0)?,
        storage_path: row.get(1)?,
        checksum_type: row.get(2)?,
        checksum: row.get(3)?,
        bytes: bytes as u64,
        created_at: row.get(5)?,
    })
}
