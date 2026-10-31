//! Automatic cleanup delegates proof, ownership and recoverable moves to core.
use anyhow::{Result, bail};
use rnaa_core::{
    CleanupMode, CleanupWhen, Database, ProjectConfig, ProjectPaths, RunRecord,
    generations::{CleanupPolicy, FileCleanupOptions},
    traits::Cleaner,
};

#[derive(Debug, Clone, Default)]
pub struct TrashCleaner;

impl TrashCleaner {
    pub fn after_command(
        &self,
        db: &Database,
        config: &ProjectConfig,
        succeeded: bool,
    ) -> Result<()> {
        if !self.apply(db, config, succeeded, None)? {
            bail!(
                "input cleanup is deferred because project readers or workers are active; retry storage files clean --apply"
            );
        }
        Ok(())
    }

    /// A busy maintenance boundary is retried at a later idle scheduler tick.
    pub fn try_completed_run(
        &self,
        db: &Database,
        config: &ProjectConfig,
        run: &str,
    ) -> Result<bool> {
        self.apply(db, config, true, Some(run.into()))
    }

    fn apply(
        &self,
        db: &Database,
        config: &ProjectConfig,
        succeeded: bool,
        run: Option<String>,
    ) -> Result<bool> {
        if config.quant.cleanup == CleanupMode::None
            || config.quant.cleanup_on == CleanupWhen::Never
            || (config.quant.cleanup_on == CleanupWhen::Success && !succeeded)
        {
            return Ok(true);
        }
        let report = db.cleanup_files(FileCleanupOptions {
            mode: config.quant.cleanup,
            run,
            apply: true,
            purge_expired: false,
            trash_days: config.quant.trash_days,
            expected_policy: Some(CleanupPolicy::from_config(config)),
        })?;
        if report.busy {
            return Ok(false);
        }
        if !report.errors.is_empty() {
            bail!("input cleanup held files: {}", report.errors.join("; "));
        }
        Ok(true)
    }
}

impl Cleaner for TrashCleaner {
    fn cleanup(&self, run: &RunRecord, paths: &ProjectPaths, config: &ProjectConfig) -> Result<()> {
        if !self.apply(
            &Database::open(&paths.root)?,
            config,
            true,
            Some(run.run_accession.clone()),
        )? {
            bail!("project inputs are in use; retry cleanup");
        }
        Ok(())
    }
}
