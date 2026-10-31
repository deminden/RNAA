//! Configuration observations fence fresh and historical reference publication.
use super::*;
use rnaa_core::references::ReferencePolicy;

pub(super) struct ReferenceObservation {
    policy: ReferencePolicy,
    disk: Option<ReferencePolicy>,
}
impl ReferenceObservation {
    pub(super) fn capture(db: &Database, config: &ProjectConfig) -> Result<Self> {
        let value = Self {
            policy: ReferencePolicy::from_config(config),
            disk: Self::disk_policy(db)?,
        };
        value.verify_catalog(db)?;
        Ok(value)
    }
    fn disk_policy(db: &Database) -> Result<Option<ReferencePolicy>> {
        let path = ProjectPaths::new(
            db.path()
                .parent()
                .context("reference project has no root")?,
        )
        .config_path;
        if path.exists() {
            Ok(Some(ReferencePolicy::from_config(&ProjectConfig::load(
                &path,
            )?)))
        } else {
            Ok(None)
        }
    }
    pub(super) fn verify_catalog(&self, db: &Database) -> Result<()> {
        let config: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
        if ReferencePolicy::from_config(&config) != self.policy {
            bail!("reference settings changed during preparation; retry");
        }
        Ok(())
    }
    pub(super) fn verify_current(&self, db: &Database) -> Result<()> {
        self.verify_catalog(db)?;
        if Self::disk_policy(db)? != self.disk {
            bail!("reference settings changed on disk during preparation; retry");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reference_observation_fences_scientific_policy_but_allows_execution_changes() -> Result<()> {
        let root = tempfile::tempdir()?;
        let paths = ProjectPaths::new(root.path());
        paths.ensure_layout()?;
        let mut config = ProjectConfig::default();
        config.save(&paths.config_path)?;
        let db = Database::create(root.path(), &config)?;
        let snapshot = ReferenceObservation::capture(&db, &config)?;
        config.refs.index_threads += 1;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
        snapshot.verify_current(&db)?;
        config.refs.index_k = 29;
        config.save(&paths.config_path)?;
        assert!(snapshot.verify_current(&db).is_err());
        config.refs.index_k = 31;
        config.save(&paths.config_path)?;
        config.refs.custom_gtf = "changed.gtf".into();
        db.set_project_config(&config)?;
        assert!(snapshot.verify_catalog(&db).is_err());
        Ok(())
    }
}
