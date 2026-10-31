//! Optional, deliberately narrow current-cohort validation for frozen outputs.
use anyhow::{Context, Result, bail};
use rnaa_core::{
    Database, ProjectConfig, ProjectPaths, cohorts::CohortSnapshot, generations::GenerationHead,
};

pub(super) struct Guard {
    snapshot: CohortSnapshot,
    digest: String,
}

impl Guard {
    pub(super) fn capture(
        paths: &ProjectPaths,
        db: &Database,
        config: &ProjectConfig,
        heads: &[GenerationHead],
    ) -> Result<Self> {
        let snapshot = CohortSnapshot::capture(db)?;
        let digest = snapshot.digest()?;
        for head in heads {
            let request = db.scientific_request(&head.request_sha256)?.context("selected generation lacks a verified scientific request; cannot check current cohort")?;
            if request.dependencies.get("active_cohort") != Some(&digest) {
                bail!(
                    "selected {} generation does not match the current active cohort",
                    head.output_group.as_str()
                );
            }
        }
        let guard = Self { snapshot, digest };
        guard.verify(paths, db, config)?;
        Ok(guard)
    }

    pub(super) fn digest(&self) -> &str {
        &self.digest
    }

    pub(super) fn verify(
        &self,
        paths: &ProjectPaths,
        db: &Database,
        config: &ProjectConfig,
    ) -> Result<()> {
        self.snapshot.verify_catalog(db)?;
        let persisted: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
        let disk = ProjectConfig::load(&paths.config_path)?;
        for observed in [config, &persisted, &disk] {
            self.snapshot.verify_config(observed)?;
            db.verify_cohort_config(observed)?;
        }
        Ok(())
    }
}

#[cfg(all(test, any(target_os = "macos", target_os = "linux")))]
mod tests {
    use super::*;
    use rnaa_core::jobs::{JobLease, JobOutcome};
    use std::fs;

    #[test]
    fn changed_catalog_at_publication_restores_staging() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let paths = ProjectPaths::new(temp.path());
        let config = ProjectConfig::default();
        config.save(&paths.config_path)?;
        let db = Database::create(&paths.root, &config)?;
        let mut lease =
            JobLease::try_acquire(&db, "export", None, &[])?.context("missing export lease")?;
        let owned = lease.database();
        let guard = Guard::capture(&paths, &owned, &config, &[])?;
        let staged = paths.root.join("staged");
        let destination = paths.root.join("published");
        fs::create_dir(&staged)?;
        fs::write(staged.join("data"), b"unchanged")?;
        let mut changed = config.clone();
        changed.resolver.fields.push("changed".into());
        db.set_project_config(&changed)?;
        let result =
            super::super::publication::publish(&owned, &staged, &destination, |transaction| {
                guard.verify(&paths, transaction, &config)
            });
        assert!(result.is_err());
        assert!(!destination.exists());
        assert_eq!(fs::read(staged.join("data"))?, b"unchanged");
        lease.finish(JobOutcome::Failed, Some("test"))?;
        Ok(())
    }
}
