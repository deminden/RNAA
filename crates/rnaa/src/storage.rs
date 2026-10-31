//! Storage admission accounting and shared-blob materialization.
use crate::cli::{
    FileStorageCommand, GenerationStorageCommand, StorageArgs, StorageCommand, parse_cleanup_mode,
};
use anyhow::{Context, Result, bail};
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::shared_store::{BlobPublication, GcOptions, SharedStore};
use rnaa_core::util::{dir_size, parse_size_bytes};
use std::path::{Path, PathBuf};

pub(crate) fn configured_storage_limit_bytes(config: &ProjectConfig) -> Result<Option<u64>> {
    let raw = config.storage.max_active_size.trim();
    if raw.is_empty() {
        return Ok(None);
    }
    Ok(Some(parse_size_bytes(raw)?))
}

pub(crate) fn project_active_storage_bytes(paths: &ProjectPaths) -> u64 {
    [
        "raw",
        "quant",
        "de",
        "corr",
        "refs",
        "metadata",
        ".rnaa_support",
        ".rnaa_generations",
        ".rnaa_transfers",
    ]
    .into_iter()
    .map(|directory| dir_size(&paths.root.join(directory)))
    .sum()
}

pub(crate) fn prepare_shared_blob(
    db: &Database,
    config: &ProjectConfig,
    source_path: &Path,
    checksum: &str,
    bytes: u64,
) -> Result<Option<BlobPublication>> {
    let Some(shared_root) = resolve_shared_root(
        config,
        db.path().parent().context("project database has no root")?,
    ) else {
        return Ok(None);
    };
    SharedStore::open(&shared_root)?
        .prepare(db, source_path, checksum, bytes)
        .map(Some)
}

fn resolve_shared_root(config: &ProjectConfig, project_root: &Path) -> Option<PathBuf> {
    let configured = config.storage.shared_root.trim();
    if !configured.is_empty() {
        let path = PathBuf::from(configured);
        return Some(if path.is_absolute() {
            path
        } else {
            project_root.join(path)
        });
    }
    std::env::var("RNAA_SHARED_ROOT").ok().and_then(|value| {
        let trimmed = value.trim();
        if trimmed.is_empty() {
            None
        } else {
            let path = PathBuf::from(trimmed);
            Some(if path.is_absolute() {
                path
            } else {
                project_root.join(path)
            })
        }
    })
}

pub(crate) fn cmd_storage(project_root: &Path, args: StorageArgs) -> Result<()> {
    let config_path = ProjectPaths::new(project_root).config_path;
    // Explicit legacy-store inspection needs only the declared catalogs. An old
    // project's unrelated configuration must not force a settings migration.
    let standalone_import =
        args.shared_root.is_some() && matches!(&args.command, StorageCommand::ImportLegacy { .. });
    let config = if config_path.exists() && !standalone_import {
        ProjectConfig::load(&config_path)?
    } else {
        ProjectConfig::default()
    };
    let open_store = || -> Result<SharedStore> {
        let root = args
            .shared_root
            .clone()
            .or_else(|| resolve_shared_root(&config, project_root))
            .context("no shared store configured; provide --shared-root or storage.shared_root")?;
        SharedStore::open(&root)
    };
    match args.command {
        StorageCommand::Usage {
            include_shared,
            details,
        } => {
            let shared = if include_shared {
                Some(
                    args.shared_root
                        .clone()
                        .or_else(|| resolve_shared_root(&config, project_root))
                        .context(
                            "--include-shared requires --shared-root or storage.shared_root",
                        )?,
                )
            } else {
                if args.shared_root.is_some() {
                    bail!("use --include-shared with an explicit --shared-root for accounting");
                }
                None
            };
            let report = Database::open(project_root)?.storage_usage(shared.as_deref(), details)?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            if !report.complete {
                bail!("storage accounting could not inspect all paths; inspect the JSON errors");
            }
        }
        StorageCommand::Files { command } => {
            if args.shared_root.is_some() {
                bail!("project file operations use --root, not --shared-root");
            }
            cmd_files(&Database::open(project_root)?, &config, command)?;
        }
        StorageCommand::Generations { command } => {
            if args.shared_root.is_some() {
                bail!("project generation operations use --root, not --shared-root");
            }
            cmd_generations(&Database::open(project_root)?, &config, command)?;
        }
        StorageCommand::ImportLegacy {
            project,
            apply,
            complete_project_inventory,
        } => {
            let report = open_store()?.import_legacy(
                &project,
                rnaa_core::shared_store::LegacyImportOptions {
                    apply,
                    complete_project_inventory,
                },
            )?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            if !report.busy_projects.is_empty() {
                bail!(
                    "legacy projects are busy; retry shared import after their workers and readers finish"
                );
            }
            if !report.errors.is_empty() {
                bail!(
                    "legacy shared import held invalid or incomplete evidence; inspect the JSON diagnostics"
                );
            }
        }
        StorageCommand::Status => {
            let store = open_store()?;
            println!(
                "{}",
                serde_json::to_string_pretty(
                    &serde_json::json!({"schema_version":1,"shared_root":store.root(),"blobs":store.status()?,"unindexed_files":"preserved"})
                )?
            );
        }
        StorageCommand::Gc {
            apply,
            purge_expired,
            trash_days,
        } => {
            let store = open_store()?;
            let report = store.gc(GcOptions {
                apply,
                purge_expired,
                trash_days: trash_days.unwrap_or(config.quant.trash_days),
            })?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            if !report.errors.is_empty() {
                bail!(
                    "shared collection retained {} blob(s) with errors; inspect the JSON diagnostics",
                    report.errors.len()
                );
            }
        }
        StorageCommand::Restore { blob } => {
            let store = open_store()?;
            let hash = blob.strip_prefix("sha256:").unwrap_or(&blob);
            store.restore(hash)?;
            println!(
                "{}",
                serde_json::json!({"blob_id":format!("sha256:{hash}"),"state":"LIVE"})
            );
        }
    }
    Ok(())
}

fn cmd_files(db: &Database, config: &ProjectConfig, command: FileStorageCommand) -> Result<()> {
    use rnaa_core::generations::FileCleanupOptions;
    let options = match command {
        FileStorageCommand::Status => {
            println!(
                "{}",
                serde_json::to_string_pretty(&db.file_trash_status()?)?
            );
            return Ok(());
        }
        FileStorageCommand::Restore { trash_id } => {
            println!(
                "{}",
                serde_json::to_string_pretty(&db.restore_file_trash(&trash_id)?)?
            );
            return Ok(());
        }
        FileStorageCommand::Clean {
            mode,
            run,
            apply,
            trash_days,
        } => FileCleanupOptions {
            mode: mode
                .as_deref()
                .map(parse_cleanup_mode)
                .transpose()?
                .unwrap_or(config.quant.cleanup),
            run,
            apply,
            purge_expired: false,
            trash_days: trash_days.unwrap_or(config.quant.trash_days),
            expected_policy: None,
        },
        FileStorageCommand::Gc {
            apply,
            purge_expired,
            run,
        } => FileCleanupOptions {
            mode: rnaa_core::CleanupMode::None,
            run,
            apply,
            purge_expired,
            trash_days: config.quant.trash_days,
            expected_policy: None,
        },
    };
    let report = db.cleanup_files(options)?;
    println!("{}", serde_json::to_string_pretty(&report)?);
    if report.busy {
        bail!("project inputs are in use; retry after readers and workers finish");
    }
    if !report.errors.is_empty() {
        bail!(
            "input cleanup held {} file(s); inspect the JSON diagnostics",
            report.errors.len()
        );
    }
    Ok(())
}

fn cmd_generations(
    db: &Database,
    config: &ProjectConfig,
    command: GenerationStorageCommand,
) -> Result<()> {
    match command {
        GenerationStorageCommand::Unpublished { command } => {
            cmd_unpublished(db, config, command)?;
        }
        GenerationStorageCommand::ImportLineage { apply } => {
            let report = db.import_generation_lineage(apply)?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            if report.busy {
                bail!("project generations are in use; retry ancestry import after readers finish");
            }
            if !report.errors.is_empty() {
                bail!(
                    "ancestry import held {} legacy generation(s); inspect the JSON reasons",
                    report.errors.len()
                );
            }
        }
        GenerationStorageCommand::Status => println!(
            "{}",
            serde_json::to_string_pretty(&db.generation_storage_status()?)?
        ),
        GenerationStorageCommand::Gc {
            apply,
            purge_expired,
            trash_days,
        } => {
            let report = db.generation_gc(rnaa_core::generations::GenerationGcOptions {
                apply,
                purge_expired,
                trash_days: trash_days.unwrap_or(config.quant.trash_days),
            })?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            if report.busy {
                bail!(
                    "project generations are in use; retry after workers and pinned readers finish"
                );
            }
            if !report.errors.is_empty() {
                bail!(
                    "generation cleanup preserved {} invalid result(s); inspect the JSON errors",
                    report.errors.len()
                );
            }
        }
        GenerationStorageCommand::Restore { generation } => {
            db.generation_restore(&generation)?;
            println!(
                "{}",
                serde_json::json!({"generation_id":generation,"state":"LIVE"})
            );
        }
        GenerationStorageCommand::Keep { generation, label } => {
            db.generation_keep(&generation, Some(&label))?;
            println!(
                "{}",
                serde_json::json!({"generation_id":generation,"keep":label})
            );
        }
        GenerationStorageCommand::Unkeep { generation } => {
            db.generation_keep(&generation, None)?;
            println!(
                "{}",
                serde_json::json!({"generation_id":generation,"keep":null})
            );
        }
    }
    Ok(())
}

fn cmd_unpublished(
    db: &Database,
    config: &ProjectConfig,
    command: crate::cli::UnpublishedStorageCommand,
) -> Result<()> {
    use crate::cli::UnpublishedStorageCommand;
    match command {
        UnpublishedStorageCommand::Status => {
            println!(
                "{}",
                serde_json::to_string_pretty(&db.unpublished_generation_status()?)?
            );
        }
        UnpublishedStorageCommand::Restore { generation } => {
            db.restore_unpublished_generation(&generation)?;
            println!(
                "{}",
                serde_json::json!({"generation_id":generation,"state":"ABANDONED","published":false})
            );
        }
        UnpublishedStorageCommand::Gc {
            apply,
            purge_expired,
            trash_days,
        } => {
            let report =
                db.unpublished_generation_gc(rnaa_core::generations::GenerationGcOptions {
                    apply,
                    purge_expired,
                    trash_days: trash_days.unwrap_or(config.quant.trash_days),
                })?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            if report.busy {
                bail!(
                    "project is busy; retry unpublished cleanup after readers and workers finish"
                );
            }
            if !report.errors.is_empty() {
                bail!("unpublished cleanup held invalid or pending data; inspect the JSON errors");
            }
        }
    }
    Ok(())
}
