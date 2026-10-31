//! Project opening, root selection and logging setup.
use anyhow::Context;
use anyhow::Result;
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::paths::ProjectPaths;
use std::fs;
use std::path::{Path, PathBuf};
use tracing_subscriber::EnvFilter;
use tracing_subscriber::prelude::*;

pub(crate) fn resolve_root(root: Option<PathBuf>) -> Result<PathBuf> {
    match root {
        Some(root) => Ok(root),
        None => std::env::current_dir().context("failed to determine current directory"),
    }
}

pub(crate) fn open_project(root: &Path) -> Result<(ProjectPaths, Database, ProjectConfig)> {
    let (paths, db, config) = open_ingest_project(root)?;
    let paths = paths.with_generation_heads(&db)?;
    Ok((paths, db, config))
}

/// Ingest resolves pinned snapshots per worker. Between workers it may restore
/// or trash inputs, so its command shell must not retain a project reader pin.
pub(crate) fn open_ingest_project(root: &Path) -> Result<(ProjectPaths, Database, ProjectConfig)> {
    let paths = ProjectPaths::new(root);
    let db = Database::open(root)?;
    let config = if paths.config_path.exists() {
        ProjectConfig::load(&paths.config_path)?
    } else {
        let project = db.project()?;
        serde_json::from_str(&project.config_json)
            .context("failed to parse config_json from database")?
    };
    Ok((paths, db, config))
}

pub(crate) fn init_logging(
    paths: &ProjectPaths,
    command: &str,
    log_file: Option<&Path>,
) -> Result<tracing_appender::non_blocking::WorkerGuard> {
    fs::create_dir_all(&paths.logs_dir)
        .with_context(|| format!("failed to create {}", paths.logs_dir.display()))?;
    let target_path = log_file
        .map(PathBuf::from)
        .unwrap_or_else(|| paths.logs_dir.join(format!("{command}.log")));
    if let Some(parent) = target_path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let file_name = target_path
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or_else(|| anyhow::anyhow!("invalid log file path: {}", target_path.display()))?;
    let file_dir = target_path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("invalid log file path: {}", target_path.display()))?;
    let file_appender = tracing_appender::rolling::never(file_dir, file_name);
    let (non_blocking, guard) = tracing_appender::non_blocking(file_appender);
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    let fmt_layer = tracing_subscriber::fmt::layer()
        .with_writer(non_blocking)
        .with_ansi(false);

    let _ = tracing_subscriber::registry()
        .with(filter)
        .with(fmt_layer)
        .try_init();
    Ok(guard)
}
