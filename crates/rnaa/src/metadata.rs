//! Canonical samplesheet projection from stored run metadata.
use anyhow::Result;
use rnaa_core::Database;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::samplesheet::load_or_init_column_map;

pub(crate) fn cmd_refresh(
    paths: &ProjectPaths,
    db: &Database,
    policy: rnaa_core::samplesheet::publication::ConflictPolicy,
) -> Result<()> {
    use rnaa_core::jobs::{JobLease, JobOutcome};
    let mut lease = JobLease::try_acquire(db, "metadata", None, &[])?.ok_or_else(|| {
        anyhow::anyhow!("metadata refresh requires the analysis/resolver to finish")
    })?;
    match refresh_owned(paths, &lease.database(), policy, true) {
        Ok(()) => {
            lease.finish(JobOutcome::Succeeded, None)?;
            rnaa_core::samplesheet::publication::recover_view(db, paths, policy)
        }
        Err(error) => {
            lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))?;
            Err(error)
        }
    }
}

pub(crate) fn refresh_samplesheet(paths: &ProjectPaths, db: &Database) -> Result<()> {
    if db
        .with_idle_analysis(|db| refresh_owned(paths, db, Default::default(), false))?
        .is_none()
    {
        tracing::info!("metadata refresh deferred while analysis/resolution is active");
    }
    Ok(())
}

fn refresh_owned(
    paths: &ProjectPaths,
    db: &Database,
    policy: rnaa_core::samplesheet::publication::ConflictPolicy,
    complete_attempt: bool,
) -> Result<()> {
    let observed = rnaa_core::cohorts::CohortSnapshot::capture(db)?;
    let runs = db.list_active_runs()?;
    let artifacts = db.list_artifacts()?;
    let metadata_columns = db.list_metadata_columns()?;
    let column_map = load_or_init_column_map(&paths.column_map_path(), &metadata_columns)?;
    let publication = rnaa_core::samplesheet::publication::PreparedPublication::prepare(
        db,
        paths,
        &runs,
        &artifacts,
        &column_map,
        policy,
    )?;
    publication.verify_sources(paths)?;
    db.atomic(|db| {
        observed.verify_catalog(db)?;
        publication.commit_catalog_preverified(db)?;
        if complete_attempt {
            db.complete_owned_job(rnaa_core::jobs::JobOutcome::Succeeded, None)?;
        }
        Ok(())
    })?;
    if !complete_attempt {
        rnaa_core::samplesheet::publication::recover_view(db, paths, policy)?;
    }
    Ok(())
}
