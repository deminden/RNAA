//! Audited legacy projection and analysis-writer exclusion for schema upgrades.
use super::*;
use std::{
    fs::File,
    time::{Duration, Instant},
};

#[derive(Debug, serde::Serialize)]
pub struct LegacyRunAnalysisState {
    pub run_accession: String,
    pub project_id: String,
    pub analysis_state: RunState,
    pub analysis_error: Option<String>,
    pub source_updated_at: String,
    pub projected_at: String,
}

impl Database {
    pub(super) fn run_progress_upgrade_guard(conn: &Connection) -> Result<Option<File>> {
        let Some(path) = conn.path().filter(|path| !path.is_empty()) else {
            return Ok(None);
        };
        let deadline = Instant::now() + Duration::from_secs(2);
        loop {
            let exists: bool = conn.query_row("SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name='schema_migrations')", [], |row| row.get(0))?;
            if exists && conn.query_row("SELECT EXISTS(SELECT 1 FROM schema_migrations WHERE version='0021_run_progress')", [], |row| row.get::<_, bool>(0))? { return Ok(None); }
            if let Some(lock) = crate::jobs::analysis_upgrade_lock(Path::new(path))? {
                return Ok(Some(lock));
            }
            if Instant::now() >= deadline {
                bail!(
                    "run progress upgrade requires an idle cohort analysis worker; retry when it finishes"
                );
            }
            std::thread::sleep(Duration::from_millis(25));
        }
    }

    pub fn legacy_run_analysis_states(&self) -> Result<Vec<LegacyRunAnalysisState>> {
        self.read(|connection| {
            let mut query = connection.prepare("SELECT run_accession,project_id,analysis_state,analysis_error,source_updated_at,projected_at FROM legacy_run_analysis_states ORDER BY run_accession")?;
            let values = query.query_map([], |row| Ok((row.get::<_, String>(0)?,row.get::<_, String>(1)?,row.get::<_, String>(2)?,row.get::<_, Option<String>>(3)?,row.get::<_, String>(4)?,row.get::<_, String>(5)?)))?.collect::<std::result::Result<Vec<_>, _>>()?;
            values.into_iter().map(|(run_accession,project_id,state,analysis_error,source_updated_at,projected_at)| Ok(LegacyRunAnalysisState { run_accession,project_id,analysis_state:RunState::from_str(&state)?,analysis_error,source_updated_at,projected_at })).collect()
        })
    }
}

#[cfg(test)]
mod tests;
