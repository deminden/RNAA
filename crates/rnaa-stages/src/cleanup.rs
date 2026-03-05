use anyhow::Result;
use rnaa_core::config::ProjectConfig;
use rnaa_core::model::RunRecord;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::Cleaner;

#[derive(Debug, Clone, Default)]
pub struct TrashCleaner;

impl Cleaner for TrashCleaner {
    fn cleanup(
        &self,
        _run: &RunRecord,
        _paths: &ProjectPaths,
        _config: &ProjectConfig,
    ) -> Result<()> {
        Ok(())
    }
}
