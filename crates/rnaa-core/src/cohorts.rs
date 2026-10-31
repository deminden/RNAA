//! Durable selection and biological identity, separate from retained run history.
use crate::fingerprint::StageRequest;
use crate::util::now_rfc3339;
use crate::{Database, InputType, MetadataMergeStrategy, RunRecord, RunState};
use anyhow::{Context, Result, bail};
use rusqlite::{OptionalExtension, params};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SelectorContract {
    inputs: Vec<(String, InputType)>,
    overrides: Vec<(i64, String, MetadataMergeStrategy)>,
    resolver: crate::ResolverConfig,
}
impl SelectorContract {
    pub fn verify_config(&self, config: &crate::ProjectConfig) -> Result<()> {
        if self.resolver != config.resolver {
            bail!("resolver configuration changed; resolve the cohort again");
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum IdentitySource {
    Declared,
    SampleAccession,
    RunFallback,
}
impl IdentitySource {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Declared => "declared",
            Self::SampleAccession => "sample_accession",
            Self::RunFallback => "run_fallback",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BiologicalSample {
    pub id: String,
    pub source: IdentitySource,
}
impl BiologicalSample {
    pub fn for_run(run: &RunRecord) -> Result<Self> {
        if let Some(id) = run.metadata.get("biological_sample_id") {
            if id.trim().is_empty() {
                bail!("empty biological_sample_id for {}", run.run_accession);
            }
            return Ok(Self {
                id: id.trim().into(),
                source: IdentitySource::Declared,
            });
        }
        if let Some(id) = run
            .sample_accession
            .as_deref()
            .map(str::trim)
            .filter(|id| !id.is_empty())
        {
            return Ok(Self {
                id: id.to_ascii_uppercase(),
                source: IdentitySource::SampleAccession,
            });
        }
        Ok(Self {
            id: run.run_accession.clone(),
            source: IdentitySource::RunFallback,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CohortMember {
    pub run: String,
    pub sample: BiologicalSample,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ActiveCohort {
    schema_version: u32,
    pub cohort_id: String,
    project_id: String,
    selectors: SelectorContract,
    pub members: Vec<CohortMember>,
}

/// Exact catalog projection used by an analysis. Progress and timestamps are
/// deliberately excluded; metadata and sample identity remain scientific inputs.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MemberMetadata {
    member: CohortMember,
    study: String,
    experiment: Option<String>,
    sample_accession: Option<String>,
    layout: String,
    platform: Option<String>,
    instrument: Option<String>,
    remote_files: Vec<crate::RemoteFile>,
    metadata: BTreeMap<String, String>,
}
impl MemberMetadata {
    fn from_run(run: &RunRecord) -> Result<Self> {
        Ok(Self {
            member: CohortMember {
                run: run.run_accession.clone(),
                sample: BiologicalSample::for_run(run)?,
            },
            study: run.study_accession.clone(),
            experiment: run.experiment_accession.clone(),
            sample_accession: run.sample_accession.clone(),
            layout: run.library_layout.as_str().into(),
            platform: run.instrument_platform.clone(),
            instrument: run.instrument_model.clone(),
            remote_files: run.remote_files.clone(),
            metadata: run.metadata.clone(),
        })
    }
}

#[derive(Debug, Clone)]
pub struct CohortSnapshot {
    head: Option<String>,
    selectors: SelectorContract,
    members: Vec<MemberMetadata>,
}
impl CohortSnapshot {
    pub fn capture(db: &Database) -> Result<Self> {
        let head = db.active_cohort()?.map(|cohort| cohort.cohort_id);
        let selectors = db.selector_contract()?;
        let members = db
            .list_active_runs()?
            .iter()
            .map(MemberMetadata::from_run)
            .collect::<Result<_>>()?;
        let snapshot = Self {
            head,
            selectors,
            members,
        };
        snapshot.verify_catalog(db)?;
        Ok(snapshot)
    }
    pub fn digest(&self) -> Result<String> {
        // Membership determines scientific admission. Full source metadata,
        // selectors and the concrete head remain transient publication guards.
        StageRequest {
            schema_version: 1,
            stage: "active_cohort".into(),
            algorithm_version: 2,
            engines: BTreeMap::new(),
            parameters: self
                .members
                .iter()
                .map(|row| &row.member)
                .collect::<Vec<_>>(),
            inputs: Vec::new(),
            dependencies: BTreeMap::new(),
        }
        .digest()
    }
    pub fn verify_catalog(&self, db: &Database) -> Result<()> {
        let members = db
            .list_active_runs()?
            .iter()
            .map(MemberMetadata::from_run)
            .collect::<Result<Vec<_>>>()?;
        if self.head != db.active_cohort()?.map(|cohort| cohort.cohort_id)
            || self.selectors != db.selector_contract()?
            || self.members != members
        {
            bail!(
                "active cohort or biological sample metadata changed during analysis; retry analysis"
            );
        }
        Ok(())
    }
    pub fn verify_config(&self, config: &crate::ProjectConfig) -> Result<()> {
        self.selectors.verify_config(config)
    }
}

impl Database {
    pub fn selector_contract(&self) -> Result<SelectorContract> {
        self.read(|conn| {
            let mut stmt = conn.prepare("SELECT input_id, input_type FROM inputs WHERE project_id = ?1 ORDER BY input_id, input_type")?;
            let inputs = stmt.query_map([self.project_id()], |row| Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?)))?
                .map(|row| { let (id, kind) = row?; Ok((id, kind.parse()?)) }).collect::<Result<_>>()?;
            let mut stmt = conn.prepare("SELECT id, stored_path, merge_strategy FROM metadata_overrides WHERE project_id = ?1 ORDER BY id")?;
            let overrides = stmt.query_map([self.project_id()], |row| Ok((row.get::<_, i64>(0)?, row.get::<_, String>(1)?, row.get::<_, String>(2)?)))?
                .map(|row| { let (id, path, strategy) = row?; Ok((id, path, serde_json::from_str(&strategy)?)) }).collect::<Result<_>>()?;
            let config: String = conn.query_row("SELECT config_json FROM projects WHERE project_id = ?1", [self.project_id()], |row| row.get(0))?;
            let config: crate::ProjectConfig = serde_json::from_str(&config)?;
            Ok(SelectorContract { inputs, overrides, resolver: config.resolver })
        })
    }

    pub fn verify_cohort_config(&self, config: &crate::ProjectConfig) -> Result<()> {
        if let Some(cohort) = self.active_cohort()?
            && cohort.selectors.resolver != config.resolver
        {
            bail!("resolver configuration changed; run `rnaa resolve` before processing");
        }
        Ok(())
    }

    pub fn active_cohort(&self) -> Result<Option<ActiveCohort>> {
        self.read(|conn| {
            let row: Option<(String, String, String)> = conn.query_row(
                "SELECT c.cohort_id, c.project_id, c.descriptor_json FROM cohort_heads h JOIN cohorts c ON c.cohort_id = h.cohort_id WHERE h.project_id = ?1",
                [self.project_id()], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
            ).optional()?;
            let Some((id, project, json)) = row else { return Ok(None); };
            let descriptor: ActiveCohort = serde_json::from_str(&json).context("invalid active cohort descriptor")?;
            let ids = descriptor.members.iter().map(|member| &member.run).collect::<BTreeSet<_>>();
            if descriptor.schema_version != 1 || descriptor.cohort_id != id
                || descriptor.project_id != project || project != self.project_id()
                || ids.len() != descriptor.members.len()
                || ids.iter().any(|id| id.trim().is_empty()) {
                bail!("invalid active cohort identity or member list");
            }
            Ok(Some(descriptor))
        })
    }

    pub fn list_active_runs(&self) -> Result<Vec<RunRecord>> {
        let selection = self.selector_contract()?;
        let Some(cohort) = self.active_cohort()? else {
            if !selection.inputs.is_empty() || !selection.overrides.is_empty() {
                bail!("cohort has not been resolved; run `rnaa resolve` before processing");
            }
            // Existing projects and explicitly imported runs have no selectors.
            return self.list_runs();
        };
        if cohort.selectors != selection {
            bail!("cohort selectors changed; run `rnaa resolve` before processing");
        }
        let mut catalog = self
            .list_runs()?
            .into_iter()
            .map(|run| (run.run_accession.clone(), run))
            .collect::<BTreeMap<_, _>>();
        cohort
            .members
            .iter()
            .map(|member| {
                let run = catalog
                    .remove(&member.run)
                    .with_context(|| format!("active cohort member {} is missing", member.run))?;
                if BiologicalSample::for_run(&run)? != member.sample {
                    bail!(
                        "biological identity changed for {}; resolve the cohort again",
                        member.run
                    );
                }
                Ok(run)
            })
            .collect()
    }

    pub fn list_active_runs_in_states(&self, states: &[RunState]) -> Result<Vec<RunRecord>> {
        Ok(self
            .list_active_runs()?
            .into_iter()
            .filter(|run| states.contains(&run.state))
            .collect())
    }

    pub fn publish_resolved_cohort(
        &self,
        runs: &[RunRecord],
        selectors: &SelectorContract,
        expected_head: Option<&str>,
    ) -> Result<ActiveCohort> {
        self.publish_resolved_cohort_with(runs, selectors, expected_head, |_| Ok(()))
    }

    pub fn publish_resolved_cohort_with(
        &self,
        runs: &[RunRecord],
        selectors: &SelectorContract,
        expected_head: Option<&str>,
        publish: impl FnOnce(&Database) -> Result<()>,
    ) -> Result<ActiveCohort> {
        let ids = runs
            .iter()
            .map(|run| &run.run_accession)
            .collect::<BTreeSet<_>>();
        if ids.len() != runs.len()
            || runs.iter().any(|run| {
                run.project_id != self.project_id() || run.run_accession.trim().is_empty()
            })
        {
            bail!("resolved cohort contains duplicate, empty or foreign run IDs");
        }
        let mut members = runs
            .iter()
            .map(|run| {
                Ok(CohortMember {
                    run: run.run_accession.clone(),
                    sample: BiologicalSample::for_run(run)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        members.sort_by(|a, b| a.run.cmp(&b.run));
        let cohort = ActiveCohort {
            schema_version: 1,
            cohort_id: uuid::Uuid::new_v4().to_string(),
            project_id: self.project_id().into(),
            selectors: selectors.clone(),
            members,
        };
        self.atomic(|db| {
            if db.selector_contract()? != *selectors || db.active_cohort()?.as_ref().map(|cohort| cohort.cohort_id.as_str()) != expected_head {
                bail!("resolution was superseded or selectors changed; retry `rnaa resolve`");
            }
            db.upsert_runs(runs)?;
            db.write(|conn| {
                conn.execute("INSERT INTO cohorts (cohort_id, project_id, descriptor_json, created_at) VALUES (?1, ?2, ?3, ?4)", params![cohort.cohort_id, self.project_id(), serde_json::to_string(&cohort)?, now_rfc3339()])?;
                conn.execute("INSERT INTO cohort_heads (project_id, cohort_id) VALUES (?1, ?2) ON CONFLICT(project_id) DO UPDATE SET cohort_id = excluded.cohort_id", params![self.project_id(), cohort.cohort_id])?;
                invalidate_analysis(conn, self.project_id(), "active cohort was resolved again")
            })?;
            db.append_event("resolve", None, "published active cohort", serde_json::json!({"cohort_id":cohort.cohort_id,"run_count":cohort.members.len()}))?;
            publish(db)?;
            Ok(())
        })?;
        Ok(cohort)
    }
}

pub(crate) fn invalidate_analysis(
    conn: &rusqlite::Connection,
    project: &str,
    reason: &str,
) -> Result<()> {
    conn.execute("UPDATE project_stages SET state = 'PENDING', last_error = ?1, updated_at = ?2 WHERE project_id = ?3 AND stage IN ('normalize', 'deseq2', 'corr')", params![reason, now_rfc3339(), project])?;
    Ok(())
}

#[cfg(test)]
mod tests;
