//! Preserve catalog-owned leftovers from interrupted exports without deleting data.
use anyhow::{Context, Result, bail};
use rnaa_core::{Database, ProjectPaths, export_recovery::PendingExport};
use std::{fs, path::Path};

#[cfg(unix)]
fn identity(metadata: &fs::Metadata) -> (u64, u64) {
    use std::os::unix::fs::MetadataExt;
    (metadata.dev(), metadata.ino())
}

#[cfg(not(unix))]
fn identity(_metadata: &fs::Metadata) -> (u64, u64) {
    (0, 0)
}

pub(super) fn register(db: &Database, stage: &str, path: &Path) -> Result<()> {
    let metadata = path.symlink_metadata()?;
    if !metadata.is_dir() || metadata.file_type().is_symlink() {
        bail!("export staging path is not a real directory");
    }
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .context("invalid export staging name")?;
    let (device, inode) = identity(&metadata);
    db.register_export_stage(stage, name, device, inode)
}

fn probe(path: &Path, record: &PendingExport) -> Result<Option<bool>> {
    match path.symlink_metadata() {
        Ok(metadata) => Ok(Some(
            metadata.is_dir()
                && !metadata.file_type().is_symlink()
                && identity(&metadata) == (record.device, record.inode),
        )),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(error) => Err(error.into()),
    }
}

pub(super) fn reconcile(paths: &ProjectPaths, db: &Database) -> Result<()> {
    match paths.exports_dir.symlink_metadata() {
        Ok(meta) if !meta.is_dir() || meta.file_type().is_symlink() => {
            bail!("export directory must be a real directory")
        }
        Ok(_) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error.into()),
    }
    for record in db.pending_export_recovery()? {
        db.check_job()?;
        let staged = paths.exports_dir.join(&record.staged_name);
        let published = paths.exports_dir.join(
            record
                .staged_name
                .strip_prefix('.')
                .context("invalid staged name")?,
        );
        let recovered = paths
            .exports_dir
            .join(format!(".rnaa-recovered-{}", &record.staged_name[1..]));
        let result = match (
            probe(&staged, &record)?,
            probe(&recovered, &record)?,
            probe(&published, &record)?,
        ) {
            (Some(true), None, None) => {
                super::publication::rename_exclusive(&staged, &recovered)?;
                tracing::warn!(attempt = %record.attempt_id, stage = %record.stage, path = %recovered.display(), "interrupted export preserved for inspection");
                "RECOVERED"
            }
            (None, Some(true), None) => "RECOVERED",
            (None, None, Some(true)) => {
                tracing::warn!(attempt = %record.attempt_id, stage = %record.stage, path = %published.display(), "complete export remained after interrupted publication; inspect before use");
                "ORPHANED"
            }
            (None, None, None) => "MISSING",
            _ => {
                tracing::warn!(attempt = %record.attempt_id, stage = %record.stage, "export recovery held because directory identity or destination changed");
                continue;
            }
        };
        db.finish_export_recovery(&record.attempt_id, result)?;
    }
    Ok(())
}

#[cfg(all(test, any(target_os = "macos", target_os = "linux")))]
mod tests {
    use super::*;
    use rnaa_core::{
        ProjectConfig,
        jobs::{JobLease, JobOutcome},
    };

    #[test]
    fn interrupted_directories_are_preserved_and_conflicts_are_held() -> Result<()> {
        const CHILD_ROOT: &str = "RNAA_TEST_EXPORT_RECOVERY_ROOT";
        if let Some(root) = std::env::var_os(CHILD_ROOT) {
            let paths = ProjectPaths::new(root);
            let db = Database::open(&paths.root)?;
            let lease =
                JobLease::try_acquire(&db, "export", None, &[])?.context("missing child lease")?;
            let name = format!(".rnaa-export-{}-fixture", lease.attempt_id());
            let staged = paths.exports_dir.join(name);
            fs::create_dir(&staged)?;
            fs::write(staged.join("data"), b"preserve me")?;
            register(&lease.database(), "export", &staged)?;
            if std::env::var_os("RNAA_TEST_EXPORT_RENAMED").is_some() {
                super::super::publication::rename_exclusive(
                    &staged,
                    &paths.exports_dir.join(
                        staged
                            .file_name()
                            .unwrap()
                            .to_str()
                            .unwrap()
                            .trim_start_matches('.'),
                    ),
                )?;
            }
            std::process::exit(73);
        }

        for scenario in ["staged", "recovery_move", "published"] {
            let renamed = scenario == "published";
            let temp = tempfile::tempdir()?;
            let paths = ProjectPaths::new(temp.path());
            fs::create_dir(&paths.exports_dir)?;
            let db = Database::create(&paths.root, &ProjectConfig::default())?;
            let mut child = std::process::Command::new(std::env::current_exe()?);
            child.arg("--exact").arg("export::recovery::tests::interrupted_directories_are_preserved_and_conflicts_are_held").env(CHILD_ROOT, &paths.root);
            if renamed {
                child.env("RNAA_TEST_EXPORT_RENAMED", "1");
            }
            assert_eq!(child.status()?.code(), Some(73));
            let previous = db
                .list_job_attempts(10)?
                .pop()
                .context("missing child attempt")?;
            let conn = rusqlite::Connection::open(&paths.db_path)?;
            conn.execute(
                "UPDATE job_attempts SET lease_expires_ms=0 WHERE attempt_id=?1",
                [&previous.attempt_id],
            )?;
            let mut lease = JobLease::try_acquire(&db, "survey", None, &[])?
                .context("missing recovery lease")?;
            let owned = lease.database();
            let staged_name = format!(".rnaa-export-{}-fixture", previous.attempt_id);
            let staged = paths.exports_dir.join(&staged_name);
            let recovered = paths
                .exports_dir
                .join(format!(".rnaa-recovered-{}", &staged_name[1..]));
            let published = paths.exports_dir.join(&staged_name[1..]);
            let unknown = paths
                .exports_dir
                .join(".rnaa-export-00000000-0000-4000-8000-000000000001-unknown");
            fs::create_dir(&unknown)?;
            if !renamed {
                fs::create_dir(&recovered)?;
                fs::write(recovered.join("other"), b"unrelated")?;
                reconcile(&paths, &owned)?;
                assert_eq!(fs::read(staged.join("data"))?, b"preserve me");
                assert_eq!(fs::read(recovered.join("other"))?, b"unrelated");
                assert_eq!(owned.pending_export_recovery()?.len(), 1);
                fs::remove_dir_all(&recovered)?;
                let displaced = paths.exports_dir.join("user-preserved-original");
                fs::rename(&staged, &displaced)?;
                fs::create_dir(&staged)?;
                fs::write(staged.join("other"), b"replacement")?;
                reconcile(&paths, &owned)?;
                assert_eq!(fs::read(staged.join("other"))?, b"replacement");
                assert_eq!(owned.pending_export_recovery()?.len(), 1);
                fs::remove_dir_all(&staged)?;
                fs::rename(&displaced, &staged)?;
                if scenario == "recovery_move" {
                    super::super::publication::rename_exclusive(&staged, &recovered)?;
                }
            }
            reconcile(&paths, &owned)?;
            let result = if renamed { &published } else { &recovered };
            assert_eq!(fs::read(result.join("data"))?, b"preserve me");
            assert!(unknown.is_dir());
            assert!(!staged.exists());
            assert!(owned.pending_export_recovery()?.is_empty());
            reconcile(&paths, &owned)?;
            lease.finish(JobOutcome::Succeeded, None)?;
        }
        Ok(())
    }
}
