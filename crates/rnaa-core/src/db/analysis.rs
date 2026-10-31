//! Analysis stage state and contrast catalog access.
use super::Database;
use crate::model::{ContrastSpec, ProjectStageRecord};
use crate::state::ProjectStageState;
use crate::util::now_rfc3339;
use anyhow::{Context, Result};
use rusqlite::{OptionalExtension, Row, params};
use std::str::FromStr;

impl Database {
    pub fn replace_contrasts(&self, contrasts: &[ContrastSpec]) -> Result<()> {
        self.write(|conn| {
            conn.execute(
                "DELETE FROM contrasts WHERE project_id = ?1",
                params![self.project_id],
            )?;
            for contrast in contrasts {
                conn.execute(
                    "INSERT INTO contrasts (project_id, name, factor, level_a, level_b) VALUES (?1, ?2, ?3, ?4, ?5)",
                    params![
                        self.project_id,
                        contrast.name,
                        contrast.factor,
                        contrast.level_a,
                        contrast.level_b
                    ],
                )?;
            }
            Ok(())
        })
    }

    pub fn list_contrasts(&self) -> Result<Vec<ContrastSpec>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT name, factor, level_a, level_b FROM contrasts WHERE project_id = ?1 ORDER BY name",
            )?;
            let rows = stmt.query_map(params![self.project_id], |row| {
                Ok(ContrastSpec {
                    name: row.get(0)?,
                    factor: row.get(1)?,
                    level_a: row.get(2)?,
                    level_b: row.get(3)?,
                })
            })?;
            let mut contrasts = Vec::new();
            for row in rows {
                contrasts.push(row?);
            }
            Ok(contrasts)
        })
    }

    pub fn set_project_stage_state(
        &self,
        stage: &str,
        state: ProjectStageState,
        last_error: Option<&str>,
    ) -> Result<()> {
        self.write(|conn| {
            conn.execute(
                "INSERT INTO project_stages (project_id, stage, state, last_error, updated_at)
                 VALUES (?1, ?2, ?3, ?4, ?5)
                 ON CONFLICT(project_id, stage) DO UPDATE SET
                   state = excluded.state,
                   last_error = excluded.last_error,
                   updated_at = excluded.updated_at",
                params![
                    self.project_id,
                    stage,
                    state.as_str(),
                    last_error,
                    now_rfc3339()
                ],
            )
            .with_context(|| format!("failed to set project stage state for {stage}"))?;
            Ok(())
        })
    }

    pub fn get_project_stage_state(&self, stage: &str) -> Result<Option<ProjectStageRecord>> {
        self.read(|conn| {
            conn.query_row(
                "SELECT project_id, stage, state, last_error, updated_at FROM project_stages WHERE project_id = ?1 AND stage = ?2",
                params![self.project_id, stage],
                parse_project_stage,
            )
            .optional()
            .context("failed to load project stage")
        })
    }

    pub fn list_project_stage_states(&self) -> Result<Vec<ProjectStageRecord>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT project_id, stage, state, last_error, updated_at FROM project_stages WHERE project_id = ?1 ORDER BY stage",
            )?;
            let rows = stmt.query_map(params![self.project_id], parse_project_stage)?;
            let mut stages = Vec::new();
            for row in rows {
                stages.push(row?);
            }
            Ok(stages)
        })
    }
}

fn parse_project_stage(row: &Row<'_>) -> rusqlite::Result<ProjectStageRecord> {
    let state_text: String = row.get(2)?;
    let state = ProjectStageState::from_str(&state_text).map_err(|err| {
        rusqlite::Error::FromSqlConversionFailure(
            2,
            rusqlite::types::Type::Text,
            Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                err.to_string(),
            )),
        )
    })?;
    Ok(ProjectStageRecord {
        project_id: row.get(0)?,
        stage: row.get(1)?,
        state,
        last_error: row.get(3)?,
        updated_at: row.get(4)?,
    })
}
