//! Full typed scientific requests, indexed independently of physical file paths.
use super::{GenerationHead, OutputGroup};
use crate::{Database, fingerprint::StageRequest};
use anyhow::{Context, Result, bail};
use rusqlite::{OptionalExtension, params};
use serde_json::Value;

pub(super) fn canonical_request(value: &Value) -> Result<StageRequest<Value>> {
    let mut request: StageRequest<Value> = serde_json::from_value(value.clone())?;
    if request.schema_version != 1 || serde_json::to_value(&request)? != *value {
        bail!("unsupported or incomplete scientific request");
    }
    request.inputs.sort();
    Ok(request)
}

impl Database {
    pub(super) fn store_scientific_request(&self, value: &Value) -> Result<String> {
        let request = canonical_request(value)?;
        let digest = request.digest()?;
        let json = serde_json::to_string(&request)?;
        let mut policy = request.clone();
        policy.inputs.clear();
        let policy_digest = policy.digest()?;
        self.write(|conn| {
            conn.execute("INSERT INTO scientific_requests(project_id,request_sha256,stage,policy_sha256,request_json)
                VALUES (?1,?2,?3,?4,?5) ON CONFLICT DO NOTHING",
                params![self.project_id(),digest,request.stage,policy_digest,json])?;
            let stored: (String,String,String) = conn.query_row("SELECT stage,policy_sha256,request_json FROM scientific_requests
                WHERE project_id=?1 AND request_sha256=?2",params![self.project_id(),digest], |row| Ok((row.get(0)?,row.get(1)?,row.get(2)?)))?;
            if stored != (request.stage.clone(),policy_digest.clone(),json.clone()) {
                bail!("scientific request catalog differs from its content identity");
            }
            Ok(())
        })?;
        Ok(digest)
    }

    /// Index an older request without rewriting its seal. The complete request
    /// must hash to the identity already committed by its original producer.
    pub fn index_generation_request(&self, head: &GenerationHead, value: &Value) -> Result<()> {
        let request = canonical_request(value)?;
        if request.digest()? != head.request_sha256 || request.stage != head.stage {
            bail!("historical scientific request differs from its committed identity");
        }
        self.atomic_or_current(|db| {
            if db
                .committed_seal(&head.generation_id)?
                .as_ref()
                .map(|seal| &seal.head)
                != Some(head)
            {
                bail!("scientific request requires its exact committed generation");
            }
            db.store_scientific_request(value)?;
            Ok(())
        })
    }

    pub fn scientific_request(&self, digest: &str) -> Result<Option<StageRequest<Value>>> {
        let stored: Option<(String,String,String)> = self.read(|conn| Ok(conn.query_row(
            "SELECT stage,policy_sha256,request_json FROM scientific_requests WHERE project_id=?1 AND request_sha256=?2",
            params![self.project_id(),digest], |row| Ok((row.get(0)?,row.get(1)?,row.get(2)?))).optional()?))?;
        stored
            .map(|(stage, policy_digest, json)| {
                let request = canonical_request(&serde_json::from_str(&json)?)?;
                let mut policy = request.clone();
                policy.inputs.clear();
                if request.digest()? != digest
                    || request.stage != stage
                    || policy.digest()? != policy_digest
                {
                    bail!("scientific request catalog failed identity verification");
                }
                Ok(request)
            })
            .transpose()
    }

    pub fn generation_policy_candidates(
        &self,
        group: OutputGroup,
        stage: &str,
        subject: &str,
        policy: &str,
    ) -> Result<Vec<GenerationHead>> {
        let ids = self.read(|conn| {
            let mut query = conn.prepare("SELECT g.generation_id FROM artifact_generations g JOIN scientific_requests r
                ON r.project_id=g.project_id AND r.request_sha256=g.request_sha256
                LEFT JOIN generation_heads h ON h.project_id=g.project_id AND h.output_group=g.output_group
                AND h.subject=g.subject AND h.generation_id=g.generation_id
                WHERE g.project_id=?1 AND g.output_group=?2 AND g.stage=?3 AND g.subject=?4 AND g.state='COMMITTED'
                AND NOT EXISTS(SELECT 1 FROM generation_retirements x WHERE x.generation_id=g.generation_id)
                AND r.stage=?3 AND r.policy_sha256=?5 ORDER BY (h.generation_id IS NOT NULL) DESC,g.committed_at DESC,g.generation_id DESC")?;
            Ok(query.query_map(params![self.project_id(),group.as_str(),stage,subject,policy], |row| row.get::<_,String>(0))?
                .collect::<std::result::Result<Vec<_>,_>>()?)
        })?;
        ids.iter()
            .map(|id| {
                let head = self
                    .committed_seal(id)?
                    .context("request candidate disappeared")?
                    .head;
                let mut request = self
                    .scientific_request(&head.request_sha256)?
                    .context("indexed request disappeared")?;
                request.inputs.clear();
                if head.stage != stage || request.digest()? != policy {
                    bail!("request policy differs from its committed index");
                }
                Ok(head)
            })
            .collect()
    }
}
