//! Command dispatch; workflow and storage logic live in dedicated modules.
use clap::Parser;
mod analysis;
mod bootstrap;
mod cli;
mod doctor;
mod export;
mod ingest;
mod lineage;
mod metadata;
mod progress;
mod project;
mod quality;
mod reference;
mod resources;
mod run_state;
mod status;
mod storage;
use crate::analysis::corr::cmd_corr;
use crate::analysis::de::{cmd_deseq2, cmd_normalize};
use crate::bootstrap::{init_logging, open_ingest_project, open_project, resolve_root};
use crate::cli::{Cli, Commands, JobsCommand, parse_download_preference};
use crate::ingest::{cmd_download, cmd_quant, cmd_run};
use crate::project::{cmd_add, cmd_init, cmd_refs, cmd_resolve};
use crate::quality::cmd_qc;
use crate::status::cmd_status;
use anyhow::{Result, bail};

fn main() {
    if let Err(err) = run() {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    let root = resolve_root(cli.root.clone())?;
    let log_file = cli.log_file.clone();

    match cli.command {
        Commands::Init(args) => cmd_init(&root, args, log_file.as_deref()),
        Commands::Add(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "add", log_file.as_deref())?;
            cmd_add(&paths, &db, &config, args)
        }
        Commands::Resolve(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "resolve", log_file.as_deref())?;
            cmd_resolve(&paths, &db, &config, args)
        }
        Commands::Metadata(args) => {
            let (paths, db, _config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "metadata", log_file.as_deref())?;
            match args.command {
                cli::MetadataCommand::Refresh { conflicts } => {
                    metadata::cmd_refresh(&paths, &db, cli::metadata_conflict_policy(&conflicts)?)?;
                    println!("metadata\trefreshed");
                    Ok(())
                }
                cli::MetadataCommand::Show => {
                    let publication = db.metadata_publication()?;
                    if let Some(publication) = &publication {
                        publication.verify()?;
                    }
                    println!("{}", serde_json::to_string_pretty(&publication)?);
                    Ok(())
                }
            }
        }
        Commands::Lineage(args) => {
            let (_paths, db, _config) = open_project(&root)?;
            lineage::inspect(&db, args)
        }
        Commands::Status { json } => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "status", log_file.as_deref())?;
            cmd_status(&paths, &db, &config, json)
        }
        Commands::Storage(args) => storage::cmd_storage(&root, args),
        Commands::Jobs(args) => {
            let (_paths, db, _config) = open_project(&root)?;
            match args.command {
                JobsCommand::List { json } => {
                    let jobs = db.list_job_attempts(100)?;
                    if json {
                        println!("{}", serde_json::to_string_pretty(&jobs)?);
                    } else {
                        println!(
                            "attempt_id\tstage\trun\tstate\tfencing_token\towner\tcancel_requested"
                        );
                        for job in jobs {
                            println!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                job.attempt_id,
                                job.stage,
                                job.run_accession.as_deref().unwrap_or("project"),
                                job.state,
                                job.fencing_token,
                                job.owner,
                                job.cancel_requested
                            );
                        }
                    }
                    Ok(())
                }
                JobsCommand::Cancel { attempt } => {
                    if !db.cancel_job(&attempt)? {
                        bail!("job {attempt} is not running in this project");
                    }
                    println!("cancel_requested\t{attempt}");
                    Ok(())
                }
            }
        }
        Commands::Doctor { json } => doctor::run(&root, json),
        Commands::Download(args) => {
            let (paths, db, mut config) = open_ingest_project(&root)?;
            let _log_guard = init_logging(&paths, "download", log_file.as_deref())?;
            if let Some(prefer) = &args.prefer {
                config.download.prefer = parse_download_preference(prefer)?;
            }
            if let Some(concurrency) = args.concurrency {
                config.download.concurrency = concurrency;
            }
            cmd_download(&paths, db, config, args)
        }
        Commands::Refs(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "refs", log_file.as_deref())?;
            cmd_refs(&paths, &db, config, args)
        }
        Commands::Quant(args) => {
            let (paths, db, config) = open_ingest_project(&root)?;
            let _log_guard = init_logging(&paths, "quant", log_file.as_deref())?;
            cmd_quant(&paths, &db, config, args)
        }
        Commands::Normalize(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "normalize", log_file.as_deref())?;
            cmd_normalize(&paths, &db, config, args)
        }
        Commands::Deseq2(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "deseq2", log_file.as_deref())?;
            cmd_deseq2(&paths, &db, config, args)
        }
        Commands::Corr(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "corr", log_file.as_deref())?;
            cmd_corr(&paths, &db, config, args)
        }
        Commands::Qc(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "qc", log_file.as_deref())?;
            cmd_qc(&paths, &db, &config, args)
        }
        Commands::Run(args) => {
            let (paths, db, config) = open_ingest_project(&root)?;
            let _log_guard = init_logging(&paths, "run", log_file.as_deref())?;
            cmd_run(&paths, &db, config, args)
        }
        Commands::Export(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "export", log_file.as_deref())?;
            export::run(&paths, &db, &config, args)
        }
        Commands::Survey(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "survey", log_file.as_deref())?;
            export::survey(&paths, &db, &config, args)
        }
    }
}
