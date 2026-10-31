//! Read-only diagnostics for project configuration and the SQLite catalog.
use anyhow::{Context, Result, bail};
use rnaa_core::{config::ProjectConfig, paths::ProjectPaths};
use rusqlite::{Connection, OpenFlags};
use serde::Serialize;
use std::path::Path;

#[derive(Serialize)]
struct Check {
    name: &'static str,
    status: &'static str,
    detail: String,
}

#[derive(Serialize)]
struct Report {
    schema_version: u32,
    root: String,
    project_id: Option<String>,
    checks: Vec<Check>,
    engines: serde_json::Value,
    ok: bool,
}

fn check(name: &'static str, result: Result<String>) -> Check {
    match result {
        Ok(detail) => Check {
            name,
            status: "ok",
            detail,
        },
        Err(error) => Check {
            name,
            status: "error",
            detail: format!("{error:#}"),
        },
    }
}

fn inspect_catalog(path: &Path) -> Result<(String, String)> {
    let db = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    let result: String = db.query_row("PRAGMA quick_check(1)", [], |row| row.get(0))?;
    if result != "ok" {
        bail!("SQLite quick_check: {result}");
    }
    let mut statement = db.prepare("SELECT project_id, config_json FROM projects LIMIT 2")?;
    let ids = statement
        .query_map([], |row| {
            Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
        })?
        .collect::<rusqlite::Result<Vec<_>>>()?;
    if ids.len() != 1 {
        bail!("expected exactly one project record");
    }
    Ok(ids[0].clone())
}

pub(crate) fn run(root: &Path, json: bool) -> Result<()> {
    let paths = ProjectPaths::new(root);
    let mut checks = Vec::new();
    let catalog = if paths.db_path.exists() {
        match inspect_catalog(&paths.db_path) {
            Ok((id, config)) => {
                checks.push(check("catalog", Ok(id.clone())));
                Some((id, config))
            }
            Err(error) => {
                checks.push(check("catalog", Err(error)));
                None
            }
        }
    } else {
        checks.push(Check {
            name: "catalog",
            status: "uninitialized",
            detail: "No project catalog found".into(),
        });
        None
    };
    let project_id = catalog.as_ref().map(|(id, _)| id.clone());
    let config = if paths.config_path.exists() {
        Some(("file", ProjectConfig::load(&paths.config_path)))
    } else if let Some((_, stored)) = &catalog {
        Some((
            "catalog",
            (|| {
                let config: ProjectConfig = serde_json::from_str(stored)
                    .context("failed to parse catalog configuration")?;
                config.validate().context("invalid catalog configuration")?;
                Ok(config)
            })(),
        ))
    } else {
        None
    };
    if let Some((source, config)) = config {
        match config {
            Ok(config) => {
                checks.push(Check {
                    name: "configuration",
                    status: "ok",
                    detail: format!("Configuration from {source} is valid"),
                });
                checks.push(check(
                    "resources",
                    config.resources.limits().map(|limits| {
                        format!(
                            "{} CPUs, {} bytes memory, {} bytes disk headroom",
                            limits.cpu_threads, limits.memory_bytes, limits.disk_headroom_bytes
                        )
                    }),
                ));
            }
            Err(error) => checks.push(check("configuration", Err(error))),
        }
    } else {
        checks.push(Check {
            name: "configuration",
            status: "missing",
            detail: "No configuration file or readable project catalog".into(),
        });
    }
    let ok = checks.iter().all(|item| item.status != "error");
    let report = Report {
        schema_version: 1,
        root: root.display().to_string(),
        project_id,
        checks,
        engines: serde_json::json!({"kallistors": rnaa_stages::engines::KALLISTORS, "rsdeseq2": rnaa_stages::engines::RSDESEQ2, "rsfgsea": rnaa_stages::engines::RSFGSEA, "mincorr": rnaa_stages::engines::MINCORR, "mincorr_revision": rnaa_stages::engines::MINCORR_REVISION}),
        ok,
    };
    if json {
        println!("{}", serde_json::to_string_pretty(&report)?);
    } else {
        println!(
            "project_id\t{}",
            report.project_id.as_deref().unwrap_or("uninitialized")
        );
        println!("root\t{}", report.root);
        for item in &report.checks {
            println!(
                "check\t{}\t{}\t{}",
                item.name,
                item.status,
                item.detail.replace(['\n', '\r', '\t'], " ")
            );
        }
        println!("tool\tnative_http\tOK\treqwest");
        println!(
            "analysis_engines\tOK\tkallistors {} + rsdeseq2 {} (in-process Rust)",
            rnaa_stages::engines::KALLISTORS,
            rnaa_stages::engines::RSDESEQ2
        );
    }
    if !ok {
        bail!("project diagnostics found errors");
    }
    Ok(())
}
