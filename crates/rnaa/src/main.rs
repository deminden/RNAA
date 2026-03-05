use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result, bail};
use chrono::{DateTime, Utc};
use clap::{Args, Parser, Subcommand};
use regex::Regex;
use rnaa_core::config::{CleanupMode, DownloadPreference, MetadataMergeStrategy, ProjectConfig};
use rnaa_core::model::{
    ArtifactKind, ArtifactRecord, ContrastSpec, CorrelationMethod, EventRecord, InputType,
    OutputMode, RemoteFile, RemoteFileKind, RunRecord, VerifiedFile,
};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::samplesheet::{load_or_init_column_map, write_samplesheet};
use rnaa_core::state::RunState;
use rnaa_core::traits::{
    Correlator, DifferentialExpression, MatrixAdjuster, Quantifier, ReferenceManager,
};
use rnaa_core::util::{compute_sha256, dir_size, file_size, now_rfc3339, sanitize_basename};
use rnaa_core::{Database, ProjectRecord};
use rnaa_formats::matrix::{MatrixData, read_matrix_tsv, write_matrix_tsv};
use rnaa_stages::corr::{MinCorrCorrelator, Residualizer};
use rnaa_stages::deseq2::Deseq2Runner;
use rnaa_stages::download::ShellDownloader;
use rnaa_stages::quant::RKallistoQuantifier;
use rnaa_stages::refs::EnsemblReferenceManager;
use rnaa_stages::resolver::EnaResolver;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::prelude::*;

fn main() {
    if let Err(err) = run() {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    let root = resolve_root(cli.root)?;

    match cli.command {
        Commands::Init(args) => cmd_init(&root, args),
        Commands::Add(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "add")?;
            cmd_add(&paths, &db, &config, args)
        }
        Commands::Resolve(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "resolve")?;
            cmd_resolve(&paths, &db, &config, args)
        }
        Commands::Status => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "status")?;
            cmd_status(&paths, &db, &config)
        }
        Commands::Doctor => {
            let paths = ProjectPaths::new(&root);
            paths.ensure_layout()?;
            let _log_guard = init_logging(&paths, "doctor")?;
            let project = Database::open(&root).ok().and_then(|db| db.project().ok());
            cmd_doctor(&paths, project.as_ref())
        }
        Commands::Download(args) => {
            let (paths, db, mut config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "download")?;
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
            let _log_guard = init_logging(&paths, "refs")?;
            cmd_refs(&paths, &db, config, args)
        }
        Commands::Quant(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "quant")?;
            cmd_quant(&paths, &db, config, args)
        }
        Commands::Deseq2(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "deseq2")?;
            cmd_deseq2(&paths, &db, config, args)
        }
        Commands::Corr(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "corr")?;
            cmd_corr(&paths, &db, config, args)
        }
        Commands::Run(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "run")?;
            cmd_run(&paths, &db, config, args)
        }
        Commands::Export(_args) => cmd_stub("export", &root, "export"),
        Commands::Survey(_args) => cmd_stub("survey", &root, "survey"),
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "rnaa",
    version,
    about = "RNAA - Rust Nucleotide Alignment and Analytics"
)]
struct Cli {
    #[arg(long, global = true)]
    root: Option<PathBuf>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Init(InitArgs),
    Add(AddArgs),
    Resolve(ResolveArgs),
    Refs(RefsArgs),
    Download(DownloadArgs),
    Quant(QuantArgs),
    Deseq2(Deseq2Args),
    Corr(CorrArgs),
    Run(RunArgs),
    Status,
    Doctor,
    Export(ExportArgs),
    Survey(SurveyArgs),
}

#[derive(Args, Debug)]
struct InitArgs {
    #[arg(long)]
    config: Option<PathBuf>,
}

#[derive(Args, Debug)]
struct AddArgs {
    #[arg(long = "id", required = true)]
    ids: Vec<String>,

    #[arg(long = "metadata")]
    metadata: Vec<PathBuf>,

    #[arg(long, default_value = "prefer_local")]
    merge_strategy: String,
}

#[derive(Args, Debug)]
struct ResolveArgs {
    #[arg(long)]
    force: bool,
}

#[derive(Args, Debug)]
struct RefsArgs {
    #[command(subcommand)]
    command: RefsCommand,
}

#[derive(Subcommand, Debug)]
enum RefsCommand {
    Prepare(RefsPrepareArgs),
}

#[derive(Args, Debug)]
struct RefsPrepareArgs {
    #[arg(long)]
    organism: Option<String>,
    #[arg(long)]
    ensembl: Option<String>,
    #[arg(long)]
    gtf: Option<PathBuf>,
    #[arg(long)]
    cdna: Option<PathBuf>,
}

#[derive(Args, Debug)]
struct DownloadArgs {
    #[arg(long)]
    concurrency: Option<usize>,

    #[arg(long)]
    forever: bool,

    #[arg(long)]
    prefer: Option<String>,
}

#[derive(Args, Debug)]
struct QuantArgs {
    #[arg(long)]
    engine: Option<String>,
    #[arg(long)]
    threads: Option<usize>,
    #[arg(long)]
    cleanup: Option<String>,
    #[arg(long = "trash-days")]
    trash_days: Option<u64>,
}

#[derive(Args, Debug)]
struct Deseq2Args {
    #[arg(long)]
    design: Option<String>,
    #[arg(long = "contrast", num_args = 3)]
    contrast: Vec<String>,
}

#[derive(Args, Debug)]
struct CorrArgs {
    #[arg(long)]
    matrix: Option<String>,
    #[arg(long)]
    adjust: Option<String>,
    #[arg(long)]
    model: Option<String>,
    #[arg(long)]
    method: Option<String>,
    #[arg(long)]
    geneset: Option<String>,
    #[arg(long = "out")]
    out_spec: Option<String>,
}

#[derive(Args, Debug)]
struct RunArgs {
    #[arg(long, default_value_t = false)]
    no_download: bool,
    #[arg(long, default_value_t = false)]
    no_corr: bool,
}

#[derive(Args, Debug)]
struct ExportArgs {
    #[arg(long)]
    format: Option<String>,
    #[arg(long)]
    what: Option<String>,
}

#[derive(Args, Debug)]
struct SurveyArgs {
    #[arg(long)]
    stratify: Option<String>,
}

fn cmd_init(root: &Path, args: InitArgs) -> Result<()> {
    let paths = ProjectPaths::new(root);
    if paths.db_path.exists() {
        bail!("RNAA project already exists at {}", root.display());
    }
    paths.ensure_layout()?;
    let _log_guard = init_logging(&paths, "init")?;

    let mut config = if let Some(path) = args.config {
        ProjectConfig::load(&path)?
    } else {
        ProjectConfig::default()
    };
    if config.project.name == "rnaa-project"
        && let Some(name) = root.file_name().and_then(|value| value.to_str())
    {
        config.project.name = name.to_string();
    }
    config.save(&paths.config_path)?;
    let db = Database::create(root, &config)?;

    println!("project_id\t{}", db.project_id());
    println!("root\t{}", root.display());
    println!("config\t{}", paths.config_path.display());
    println!("database\t{}", paths.db_path.display());
    Ok(())
}

fn cmd_add(
    paths: &ProjectPaths,
    db: &Database,
    _config: &ProjectConfig,
    args: AddArgs,
) -> Result<()> {
    let merge_strategy = parse_merge_strategy(&args.merge_strategy)?;
    for input_id in &args.ids {
        let input_type = InputType::from_accession(input_id)?;
        db.add_input(input_id, input_type)?;
    }

    for metadata_path in &args.metadata {
        let file_name = metadata_path
            .file_name()
            .and_then(|value| value.to_str())
            .ok_or_else(|| {
                anyhow::anyhow!("invalid metadata file name {}", metadata_path.display())
            })?;
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .context("system clock before UNIX_EPOCH")?
            .as_secs();
        let dest = paths
            .overrides_dir
            .join(format!("{stamp}_{}", sanitize_basename(file_name)));
        fs::copy(metadata_path, &dest).with_context(|| {
            format!(
                "failed to copy metadata override {} to {}",
                metadata_path.display(),
                dest.display()
            )
        })?;
        db.add_metadata_override(metadata_path, &dest, merge_strategy)?;
    }

    println!("added_inputs\t{}", args.ids.len());
    println!("metadata_overrides\t{}", args.metadata.len());
    Ok(())
}

fn cmd_resolve(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    args: ResolveArgs,
) -> Result<()> {
    let resolver = EnaResolver::new(config.resolver.fields.clone())?;
    if args.force {
        tracing::info!(
            "force flag set; resolver will refresh ENA metadata and upsert existing runs"
        );
    }
    let resolved = resolver.resolve_project(db, paths, config)?;
    println!("runs\t{}", resolved.runs.len());
    println!("samplesheet\t{}", paths.samplesheet_path().display());
    println!("column_map\t{}", paths.column_map_path().display());
    Ok(())
}

fn cmd_status(paths: &ProjectPaths, db: &Database, config: &ProjectConfig) -> Result<()> {
    let project = db.project()?;
    let counts = db.state_counts()?;
    let inputs = db.list_inputs()?;
    let overrides = db.list_metadata_overrides()?;
    let disk_usage = db.disk_usage_bytes()?;
    let runs = db.list_runs()?;
    let artifacts = db.list_artifacts()?;
    let events = db.list_events(10_000)?;

    let total_runs = runs.len() as u64;
    let download_done = runs
        .iter()
        .filter(|run| is_download_done_state(run.state))
        .count() as u64;
    let quant_done = runs
        .iter()
        .filter(|run| is_quant_done_state(run.state))
        .count() as u64;
    let de_done = runs
        .iter()
        .filter(|run| is_de_done_state(run.state))
        .count() as u64;
    let corr_done = runs
        .iter()
        .filter(|run| is_corr_done_state(run.state))
        .count() as u64;

    let quant_remaining_runs = (runs
        .iter()
        .filter(|run| {
            matches!(
                run.state,
                RunState::Verified
                    | RunState::QuantRunning
                    | RunState::QuantFailed
                    | RunState::DeRunning
                    | RunState::DeDone
                    | RunState::DeFailed
                    | RunState::CorrRunning
                    | RunState::CorrDone
                    | RunState::CorrFailed
                    | RunState::Cleaned
            )
        })
        .count() as u64)
        .saturating_sub(quant_done);

    let downloaded_bytes = artifacts
        .iter()
        .filter(|artifact| {
            matches!(
                artifact.kind,
                ArtifactKind::FastqR1
                    | ArtifactKind::FastqR2
                    | ArtifactKind::FastqSingle
                    | ArtifactKind::Sra
            )
        })
        .map(|artifact| artifact.bytes)
        .sum::<u64>();
    let expected_total_bytes =
        estimate_total_download_bytes(&runs, config.download.prefer).max(downloaded_bytes);
    let remaining_download_bytes = expected_total_bytes.saturating_sub(downloaded_bytes);

    let download_rate_bps = estimate_download_speed_bps(&events);
    let quant_rate_runs_per_sec = estimate_stage_runs_per_sec(
        &events,
        "quant",
        "quantification started",
        "quantification completed",
    );
    let de_avg_secs = estimate_stage_avg_duration_secs(
        &events,
        "deseq2",
        "deseq2 stage started",
        "deseq2 stage completed",
    );
    let corr_avg_secs = estimate_stage_avg_duration_secs(
        &events,
        "corr",
        "corr stage started",
        "corr stage completed",
    );

    let download_eta = duration_from_rate(remaining_download_bytes as f64, download_rate_bps);
    let quant_eta = duration_from_rate(quant_remaining_runs as f64, quant_rate_runs_per_sec);
    let de_eta = if de_done < total_runs {
        de_avg_secs.map(std::time::Duration::from_secs_f64)
    } else {
        Some(std::time::Duration::from_secs(0))
    };
    let corr_eta = if corr_done < total_runs {
        corr_avg_secs.map(std::time::Duration::from_secs_f64)
    } else {
        Some(std::time::Duration::from_secs(0))
    };
    let processing_eta = match (quant_eta, de_eta, corr_eta) {
        (Some(q), Some(d), Some(c)) => Some(q + d + c),
        (Some(q), _, _) => Some(q),
        _ => None,
    };
    let total_eta = match (download_eta, processing_eta) {
        (Some(d), Some(p)) => Some(std::cmp::max(d, p)),
        (Some(d), None) => Some(d),
        (None, Some(p)) => Some(p),
        (None, None) => None,
    };
    let limiting_stage = match (download_eta, processing_eta) {
        (Some(d), Some(p)) if d >= p => "download",
        (Some(_), Some(_)) => "processing",
        (Some(_), None) => "download",
        (None, Some(_)) => "processing",
        (None, None) => "unknown",
    };

    println!("project_id\t{}", project.project_id);
    println!("root\t{}", project.root_dir);
    println!("inputs\t{}", inputs.len());
    println!("runs_total\t{}", total_runs);
    println!("metadata_overrides\t{}", overrides.len());
    println!("samplesheet\t{}", paths.samplesheet_path().display());
    println!(
        "disk_usage_bytes\t{}",
        disk_usage.max(dir_size(&paths.raw_dir))
    );
    println!(
        "progress_download\t{}/{} ({:.1}%)",
        download_done,
        total_runs,
        ratio_pct(download_done, total_runs)
    );
    println!(
        "progress_quant\t{}/{} ({:.1}%)",
        quant_done,
        total_runs,
        ratio_pct(quant_done, total_runs)
    );
    println!(
        "progress_deseq2\t{}/{} ({:.1}%)",
        de_done,
        total_runs,
        ratio_pct(de_done, total_runs)
    );
    println!(
        "progress_corr\t{}/{} ({:.1}%)",
        corr_done,
        total_runs,
        ratio_pct(corr_done, total_runs)
    );
    println!("downloaded_bytes\t{}", downloaded_bytes);
    println!("expected_total_download_bytes\t{}", expected_total_bytes);
    if download_rate_bps > 0.0 {
        println!(
            "download_speed_mbps\t{:.2}",
            download_rate_bps / (1024.0 * 1024.0)
        );
    } else {
        println!("download_speed_mbps\tunknown");
    }
    if quant_rate_runs_per_sec > 0.0 {
        println!(
            "quant_rate_runs_per_hour\t{:.2}",
            quant_rate_runs_per_sec * 3600.0
        );
    } else {
        println!("quant_rate_runs_per_hour\tunknown");
    }
    println!(
        "eta_download\t{}",
        download_eta
            .map(format_duration)
            .unwrap_or_else(|| "unknown".to_string())
    );
    println!(
        "eta_processing\t{}",
        processing_eta
            .map(format_duration)
            .unwrap_or_else(|| "unknown".to_string())
    );
    println!(
        "eta_total\t{}",
        total_eta
            .map(format_duration)
            .unwrap_or_else(|| "unknown".to_string())
    );
    println!("limiting_stage\t{}", limiting_stage);

    if counts.is_empty() {
        println!("states\t0");
    } else {
        println!("state_counts");
        for (state, count) in counts {
            println!("{}\t{}", state.as_str(), count);
        }
    }

    let errors = db.last_errors(10)?;
    if !errors.is_empty() {
        println!("last_errors");
        for (run, state, error) in errors {
            println!("{run}\t{state}\t{error}");
        }
    }
    Ok(())
}

fn cmd_doctor(paths: &ProjectPaths, project: Option<&ProjectRecord>) -> Result<()> {
    if let Some(project) = project {
        println!("project_id\t{}", project.project_id);
    } else {
        println!("project_id\tuninitialized");
    }
    println!("root\t{}", paths.root.display());

    for binary in ["curl", "wget", "Rscript", "kallisto", "xsra", "mincorr"] {
        match which::which(binary) {
            Ok(path) => println!("tool\t{binary}\tOK\t{}", path.display()),
            Err(_) => println!("tool\t{binary}\tMISSING"),
        }
    }

    if which::which("Rscript").is_ok() {
        let status = Command::new("Rscript")
            .args([
                "--vanilla",
                "-e",
                "pkgs <- c('jsonlite','tximport','DESeq2'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) { cat(paste(missing, collapse=',')); quit(status=1) }",
            ])
            .output()
            .context("failed to execute Rscript for package check")?;
        if status.status.success() {
            println!("r_packages\tOK");
        } else {
            println!(
                "r_packages\tMISSING\t{}",
                String::from_utf8_lossy(&status.stdout).trim()
            );
        }
    }

    if let Ok(output) = Command::new("df")
        .args(["-Pk", paths.root.to_string_lossy().as_ref()])
        .output()
    {
        let text = String::from_utf8_lossy(&output.stdout);
        if let Some(line) = text.lines().nth(1) {
            println!("disk\t{}", line.trim());
        }
    }
    Ok(())
}

fn cmd_download(
    paths: &ProjectPaths,
    db: Database,
    config: ProjectConfig,
    args: DownloadArgs,
) -> Result<()> {
    let downloader = ShellDownloader;
    let prefer = args
        .prefer
        .as_deref()
        .map(parse_download_preference)
        .transpose()?;
    downloader.run_worker(
        db.clone(),
        paths.clone(),
        config.clone(),
        args.concurrency.unwrap_or(config.download.concurrency),
        args.forever,
        prefer,
    )?;
    if !args.forever {
        refresh_samplesheet(paths, &db)?;
    }
    println!("download\tcomplete");
    Ok(())
}

fn cmd_refs(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: RefsArgs,
) -> Result<()> {
    match args.command {
        RefsCommand::Prepare(prepare) => cmd_refs_prepare(paths, db, config, prepare),
    }
}

fn cmd_refs_prepare(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    args: RefsPrepareArgs,
) -> Result<()> {
    let mut changed = false;
    if let Some(organism) = args.organism {
        config.refs.organism = organism;
        changed = true;
    }
    if let Some(ensembl) = args.ensembl {
        config.refs.ensembl = ensembl;
        changed = true;
    }
    match (args.gtf, args.cdna) {
        (Some(gtf), Some(cdna)) => {
            config.refs.custom_gtf = gtf.display().to_string();
            config.refs.custom_cdna = cdna.display().to_string();
            changed = true;
        }
        (None, None) => {}
        _ => bail!("--gtf and --cdna must be provided together"),
    }

    if changed {
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }

    let manager = EnsemblReferenceManager::new()?;
    let reference = manager.ensure_reference(paths, &config)?;
    let outputs = vec![
        (ArtifactKind::ReferenceCdna, reference.cdna_path.clone()),
        (ArtifactKind::ReferenceGtf, reference.gtf_path.clone()),
        (
            ArtifactKind::ReferenceIndex,
            reference.kallisto_index_path.clone(),
        ),
        (
            ArtifactKind::ReferenceTx2Gene,
            reference.tx2gene_path.clone(),
        ),
        (ArtifactKind::Manifest, reference.manifest_path.clone()),
    ];
    for (kind, path) in outputs {
        let artifact = ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: None,
            kind,
            path: path.display().to_string(),
            checksum_type: "sha256".to_string(),
            checksum: compute_sha256(&path)?,
            bytes: file_size(&path)?,
            created_at: now_rfc3339(),
        };
        db.record_artifact(&artifact)?;
    }
    if let Some(annotation_path) = &reference.gene_annotation_path {
        let artifact = ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: None,
            kind: ArtifactKind::ReferenceGeneAnnotation,
            path: annotation_path.display().to_string(),
            checksum_type: "sha256".to_string(),
            checksum: compute_sha256(annotation_path)?,
            bytes: file_size(annotation_path)?,
            created_at: now_rfc3339(),
        };
        db.record_artifact(&artifact)?;
    }
    db.append_event(
        "refs",
        None,
        "reference bundle prepared",
        serde_json::json!({
            "reference_id": reference.id,
            "organism": reference.organism,
            "ensembl_release": reference.ensembl_release
        }),
    )?;

    println!("reference_id\t{}", reference.id);
    println!("cdna\t{}", reference.cdna_path.display());
    println!("gtf\t{}", reference.gtf_path.display());
    println!("index\t{}", reference.kallisto_index_path.display());
    println!("manifest\t{}", reference.manifest_path.display());
    Ok(())
}

fn cmd_quant(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    args: QuantArgs,
) -> Result<()> {
    if let Some(engine) = args.engine {
        config.quant.engine = engine;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(threads) = args.threads {
        config.quant.threads = threads;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(cleanup) = args.cleanup {
        config.quant.cleanup = parse_cleanup_mode(&cleanup)?;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(trash_days) = args.trash_days {
        config.quant.trash_days = trash_days;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }

    let runs = db.list_runs_in_states(&[RunState::Verified, RunState::QuantFailed])?;
    if runs.is_empty() {
        println!("quant\tno_verified_runs");
        return Ok(());
    }
    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, &config)?;
    let quantifier = RKallistoQuantifier;

    let mut completed = 0_u64;
    let mut failed = 0_u64;
    for run in runs {
        db.set_run_state(&run.run_accession, RunState::QuantRunning, None)?;
        db.append_event(
            "quant",
            Some(&run.run_accession),
            "quantification started",
            serde_json::json!({
                "engine": config.quant.engine,
                "threads": config.quant.threads
            }),
        )?;
        let fastqs = load_fastq_artifacts(db, &run.run_accession)?;
        if fastqs.is_empty() {
            let message = format!("no FASTQ artifacts found for {}", run.run_accession);
            db.set_run_state(&run.run_accession, RunState::QuantFailed, Some(&message))?;
            db.append_event(
                "quant",
                Some(&run.run_accession),
                "quantification failed",
                serde_json::json!({ "error": message }),
            )?;
            failed += 1;
            continue;
        }

        match quantifier.quantify(&run, &fastqs, &reference, paths, &config) {
            Ok(artifacts) => {
                let records = vec![
                    ArtifactRecord {
                        project_id: db.project_id().to_string(),
                        run_accession: Some(run.run_accession.clone()),
                        kind: ArtifactKind::QuantDir,
                        path: artifacts.out_dir.display().to_string(),
                        checksum_type: "none".to_string(),
                        checksum: String::new(),
                        bytes: dir_size(&artifacts.out_dir),
                        created_at: now_rfc3339(),
                    },
                    ArtifactRecord {
                        project_id: db.project_id().to_string(),
                        run_accession: Some(run.run_accession.clone()),
                        kind: ArtifactKind::Manifest,
                        path: artifacts.manifest_path.display().to_string(),
                        checksum_type: "sha256".to_string(),
                        checksum: compute_sha256(&artifacts.manifest_path)?,
                        bytes: file_size(&artifacts.manifest_path)?,
                        created_at: now_rfc3339(),
                    },
                ];
                for record in records {
                    db.record_artifact(&record)?;
                }
                db.set_run_state(&run.run_accession, RunState::QuantDone, None)?;
                db.append_event(
                    "quant",
                    Some(&run.run_accession),
                    "quantification completed",
                    serde_json::json!({
                        "out_dir": artifacts.out_dir.display().to_string(),
                        "abundance_tsv": artifacts.abundance_tsv.display().to_string(),
                        "run_info_json": artifacts.run_info_json.display().to_string()
                    }),
                )?;
                completed += 1;
            }
            Err(err) => {
                let message = format!("{err:#}");
                db.set_run_state(&run.run_accession, RunState::QuantFailed, Some(&message))?;
                db.append_event(
                    "quant",
                    Some(&run.run_accession),
                    "quantification failed",
                    serde_json::json!({ "error": message }),
                )?;
                failed += 1;
            }
        }
    }

    println!("quant_completed\t{}", completed);
    println!("quant_failed\t{}", failed);
    Ok(())
}

fn cmd_deseq2(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    args: Deseq2Args,
) -> Result<()> {
    let runs = db.list_runs_in_states(&[RunState::QuantDone, RunState::DeFailed])?;
    if runs.is_empty() {
        println!("deseq2\tno_quantified_runs");
        return Ok(());
    }

    let samplesheet_path = paths.samplesheet_path();
    if !samplesheet_path.exists() {
        bail!(
            "samplesheet missing at {}; run `rnaa resolve` first",
            samplesheet_path.display()
        );
    }
    let header = read_samplesheet_header(&samplesheet_path)?;

    let design = args
        .design
        .clone()
        .unwrap_or_else(|| config.deseq2.design.clone());
    validate_design_formula(&header, &design)?;
    config.deseq2.design = design.clone();

    let contrasts = if args.contrast.is_empty() {
        db.list_contrasts()?
    } else {
        parse_cli_contrasts(&args.contrast)?
    };
    if !contrasts.is_empty() {
        for contrast in &contrasts {
            if !header.contains(&contrast.factor) {
                bail!(
                    "contrast factor '{}' not found in samplesheet columns",
                    contrast.factor
                );
            }
        }
    }

    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    db.replace_contrasts(&contrasts)?;

    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, &config)?;
    let runner = Deseq2Runner;
    let target_runs = runs
        .iter()
        .map(|run| run.run_accession.clone())
        .collect::<Vec<_>>();
    for run_accession in &target_runs {
        db.set_run_state(run_accession, RunState::DeRunning, None)?;
    }
    db.append_event(
        "deseq2",
        None,
        "deseq2 stage started",
        serde_json::json!({
            "design": design,
            "contrasts": contrasts,
            "runs": target_runs.len(),
        }),
    )?;

    let outputs = match runner.deseq2(
        db.project_id(),
        &reference,
        &design,
        &contrasts,
        paths,
        &config,
    ) {
        Ok(outputs) => outputs,
        Err(err) => {
            let message = format!("{err:#}");
            for run_accession in &target_runs {
                db.set_run_state(run_accession, RunState::DeFailed, Some(&message))?;
            }
            db.append_event(
                "deseq2",
                None,
                "deseq2 stage failed",
                serde_json::json!({ "error": message }),
            )?;
            return Err(err);
        }
    };

    let mut recorded = 0_u64;
    for output in &outputs {
        if !output.exists() {
            continue;
        }
        let file_name = output
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or_default();
        let kind = if file_name == "gene_counts.tsv" {
            ArtifactKind::Counts
        } else if file_name == "gene_norm_counts.tsv" {
            ArtifactKind::NormalizedCounts
        } else if file_name == "vst.tsv" {
            ArtifactKind::Vst
        } else if file_name == "vst.rds" {
            ArtifactKind::VstRds
        } else if file_name.starts_with("de_") && file_name.ends_with(".tsv") {
            ArtifactKind::DeTable
        } else if file_name.ends_with("_manifest.json") || file_name == "deseq2_manifest.json" {
            ArtifactKind::Manifest
        } else {
            continue;
        };
        db.record_artifact(&ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: None,
            kind,
            path: output.display().to_string(),
            checksum_type: "sha256".to_string(),
            checksum: compute_sha256(output)?,
            bytes: file_size(output)?,
            created_at: now_rfc3339(),
        })?;
        recorded += 1;
    }

    for run_accession in &target_runs {
        db.set_run_state(run_accession, RunState::DeDone, None)?;
    }
    db.append_event(
        "deseq2",
        None,
        "deseq2 stage completed",
        serde_json::json!({
            "outputs_recorded": recorded,
            "runs_marked": target_runs.len(),
        }),
    )?;

    println!("deseq2_outputs\t{}", outputs.len());
    println!("deseq2_recorded\t{}", recorded);
    Ok(())
}

fn cmd_corr(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    args: CorrArgs,
) -> Result<()> {
    if let Some(adjust) = args.adjust {
        config.corr.adjust = adjust;
    }
    if let Some(model) = args.model {
        config.corr.model = model;
    }
    if let Some(method) = args.method {
        config.corr.method = method;
    }
    if let Some(geneset) = args.geneset {
        config.corr.geneset = geneset;
    }
    if let Some(out_spec) = args.out_spec {
        config.corr.output = out_spec;
    }
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;

    let runs = db.list_runs_in_states(&[
        RunState::DeDone,
        RunState::CorrFailed,
        RunState::CorrRunning,
        RunState::CorrDone,
    ])?;
    if runs.is_empty() {
        bail!("no DE-complete runs found; run `rnaa deseq2` first");
    }

    let matrix_path = args
        .matrix
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or_else(|| paths.de_project_dir(db.project_id()).join("vst.tsv"));
    if !matrix_path.exists() {
        bail!(
            "missing matrix at {}; run `rnaa deseq2` first or pass --matrix",
            matrix_path.display()
        );
    }

    let selected_matrix = materialize_geneset_matrix(paths, &matrix_path, &config.corr.geneset)?;
    let adjuster = Residualizer;
    let adjusted = adjuster.adjust(
        &selected_matrix,
        &paths.samplesheet_path(),
        &config.corr.model,
        paths,
    )?;

    for run in &runs {
        db.set_run_state(&run.run_accession, RunState::CorrRunning, None)?;
    }
    db.append_event(
        "corr",
        None,
        "corr stage started",
        serde_json::json!({
            "matrix": adjusted.matrix_path.display().to_string(),
            "model": config.corr.model,
            "method": config.corr.method,
            "output": config.corr.output,
            "geneset": config.corr.geneset
        }),
    )?;

    let method = config.corr.method.parse::<CorrelationMethod>()?;
    let output_mode = parse_output_mode(&config.corr.output)?;
    let correlator = MinCorrCorrelator;
    let outputs =
        match correlator.correlate(&adjusted, method, &output_mode, paths, db.project_id()) {
            Ok(outputs) => outputs,
            Err(err) => {
                let message = format!("{err:#}");
                for run in &runs {
                    db.set_run_state(&run.run_accession, RunState::CorrFailed, Some(&message))?;
                }
                db.append_event(
                    "corr",
                    None,
                    "corr stage failed",
                    serde_json::json!({ "error": message }),
                )?;
                return Err(err);
            }
        };

    let mut recorded = 0_u64;
    for output in &outputs {
        if !output.exists() {
            continue;
        }
        let file_name = output
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or_default();
        let kind = if file_name.contains("edges") {
            ArtifactKind::CorrEdges
        } else if file_name.contains("dense") {
            ArtifactKind::CorrDense
        } else if file_name.contains("adjusted_matrix") {
            ArtifactKind::AdjustedMatrix
        } else if file_name.ends_with(".json") {
            ArtifactKind::Manifest
        } else {
            continue;
        };
        db.record_artifact(&ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: None,
            kind,
            path: output.display().to_string(),
            checksum_type: "sha256".to_string(),
            checksum: compute_sha256(output)?,
            bytes: file_size(output)?,
            created_at: now_rfc3339(),
        })?;
        recorded += 1;
    }
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().to_string(),
        run_accession: None,
        kind: ArtifactKind::AdjustedMatrix,
        path: adjusted.matrix_path.display().to_string(),
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(&adjusted.matrix_path)?,
        bytes: file_size(&adjusted.matrix_path)?,
        created_at: now_rfc3339(),
    })?;
    recorded += 1;

    for run in &runs {
        db.set_run_state(&run.run_accession, RunState::CorrDone, None)?;
    }
    db.append_event(
        "corr",
        None,
        "corr stage completed",
        serde_json::json!({
            "outputs_recorded": recorded,
            "runs_marked": runs.len()
        }),
    )?;
    println!("corr_outputs\t{}", outputs.len());
    println!("corr_recorded\t{}", recorded);
    Ok(())
}

fn cmd_run(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: RunArgs,
) -> Result<()> {
    cmd_resolve(paths, db, &config, ResolveArgs { force: false })?;
    if !args.no_download {
        cmd_download(
            paths,
            db.clone(),
            config.clone(),
            DownloadArgs {
                concurrency: None,
                forever: false,
                prefer: None,
            },
        )?;
    }
    cmd_quant(
        paths,
        db,
        config.clone(),
        QuantArgs {
            engine: None,
            threads: None,
            cleanup: None,
            trash_days: None,
        },
    )?;
    cmd_deseq2(
        paths,
        db,
        config.clone(),
        Deseq2Args {
            design: None,
            contrast: Vec::new(),
        },
    )?;
    if !args.no_corr {
        cmd_corr(
            paths,
            db,
            config,
            CorrArgs {
                matrix: None,
                adjust: None,
                model: None,
                method: None,
                geneset: None,
                out_spec: None,
            },
        )?;
    }
    println!("run\tcompleted");
    Ok(())
}

fn load_fastq_artifacts(db: &Database, run_accession: &str) -> Result<Vec<VerifiedFile>> {
    let artifacts = db.list_artifacts_for_run(Some(run_accession))?;
    let mut verified = artifacts
        .into_iter()
        .filter(|artifact| {
            matches!(
                artifact.kind,
                ArtifactKind::FastqR1 | ArtifactKind::FastqR2 | ArtifactKind::FastqSingle
            )
        })
        .map(|artifact| VerifiedFile {
            artifact_kind: artifact.kind,
            path: PathBuf::from(artifact.path),
            checksum_type: artifact.checksum_type,
            checksum: artifact.checksum,
            bytes: artifact.bytes,
            integrity_source: "artifact_db".to_string(),
        })
        .collect::<Vec<_>>();
    verified.sort_by(|left, right| left.path.cmp(&right.path));
    Ok(verified)
}

fn refresh_samplesheet(paths: &ProjectPaths, db: &Database) -> Result<()> {
    let runs = db.list_runs()?;
    let artifacts = db.list_artifacts()?;
    let metadata_columns = runs
        .iter()
        .flat_map(|run| run.metadata.keys().cloned())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let column_map = load_or_init_column_map(&paths.column_map_path(), &metadata_columns)?;
    write_samplesheet(
        &paths.samplesheet_path(),
        db.project_id(),
        &runs,
        &artifacts,
        &column_map,
    )?;
    Ok(())
}

fn is_download_done_state(state: RunState) -> bool {
    matches!(
        state,
        RunState::Verified
            | RunState::QuantRunning
            | RunState::QuantDone
            | RunState::QuantFailed
            | RunState::DeRunning
            | RunState::DeDone
            | RunState::DeFailed
            | RunState::CorrRunning
            | RunState::CorrDone
            | RunState::CorrFailed
            | RunState::Cleaned
    )
}

fn is_quant_done_state(state: RunState) -> bool {
    matches!(
        state,
        RunState::QuantDone
            | RunState::DeRunning
            | RunState::DeDone
            | RunState::DeFailed
            | RunState::CorrRunning
            | RunState::CorrDone
            | RunState::CorrFailed
            | RunState::Cleaned
    )
}

fn is_de_done_state(state: RunState) -> bool {
    matches!(
        state,
        RunState::DeDone | RunState::CorrRunning | RunState::CorrDone | RunState::CorrFailed
    )
}

fn is_corr_done_state(state: RunState) -> bool {
    matches!(state, RunState::CorrDone | RunState::Cleaned)
}

fn ratio_pct(done: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        (done as f64) * 100.0 / (total as f64)
    }
}

fn estimate_total_download_bytes(runs: &[RunRecord], prefer: DownloadPreference) -> u64 {
    runs.iter()
        .map(|run| estimate_run_download_bytes(run, prefer))
        .sum()
}

fn estimate_run_download_bytes(run: &RunRecord, prefer: DownloadPreference) -> u64 {
    let fastq_bytes = run
        .remote_files
        .iter()
        .filter(|file| {
            matches!(
                file.kind,
                RemoteFileKind::Fastq | RemoteFileKind::FastqR1 | RemoteFileKind::FastqR2
            )
        })
        .map(|file| file.bytes.unwrap_or_default())
        .sum::<u64>();
    let sra_bytes = run
        .remote_files
        .iter()
        .filter(|file| file.kind == RemoteFileKind::Sra)
        .map(|file| file.bytes.unwrap_or_default())
        .sum::<u64>();
    let preferred = match prefer {
        DownloadPreference::Fastq if fastq_bytes > 0 => fastq_bytes,
        DownloadPreference::Sra if sra_bytes > 0 => sra_bytes,
        DownloadPreference::Fastq if sra_bytes > 0 => sra_bytes,
        DownloadPreference::Sra if fastq_bytes > 0 => fastq_bytes,
        _ => 0,
    };
    if preferred > 0 {
        return preferred;
    }
    // If remote sizes are missing, fall back to a conservative estimate from any listed files.
    run.remote_files
        .iter()
        .filter_map(|file: &RemoteFile| file.bytes)
        .sum()
}

fn estimate_download_speed_bps(events: &[EventRecord]) -> f64 {
    let mut total_bytes = 0_f64;
    let mut total_seconds = 0_f64;
    for run in distinct_event_runs(events, "download", "download started", "download verified") {
        if let (Some(start), Some(done)) = (
            find_event_ts(events, "download", "download started", Some(&run)),
            find_event_ts(events, "download", "download verified", Some(&run)),
        ) && done > start
        {
            let secs = (done - start).num_milliseconds() as f64 / 1000.0;
            if secs > 0.0 {
                let run_bytes = events
                    .iter()
                    .find(|event| {
                        event.stage == "download"
                            && event.message == "download verified"
                            && event.run_accession.as_deref() == Some(run.as_str())
                    })
                    .and_then(|event| event.context.get("bytes"))
                    .and_then(|value| value.as_f64())
                    .unwrap_or(0.0);
                total_seconds += secs;
                total_bytes += run_bytes;
            }
        }
    }

    if total_bytes > 0.0 && total_seconds > 0.0 {
        return total_bytes / total_seconds;
    }

    // Fallback: coarse stage-level estimate from first start to last completion when per-run bytes are unavailable.
    let first_start = earliest_event_ts(events, "download", "download started");
    let last_done = latest_event_ts(events, "download", "download verified");
    let completed = events
        .iter()
        .filter(|event| event.stage == "download" && event.message == "download verified")
        .count();
    if completed == 0 {
        return 0.0;
    }
    let avg_bytes_per_run = events
        .iter()
        .filter(|event| event.stage == "download" && event.message == "download verified")
        .filter_map(|event| event.context.get("bytes").and_then(|value| value.as_f64()))
        .sum::<f64>()
        / (completed as f64);
    match (first_start, last_done) {
        (Some(start), Some(done)) if done > start => {
            let secs = (done - start).num_milliseconds() as f64 / 1000.0;
            if secs > 0.0 {
                (avg_bytes_per_run * completed as f64) / secs
            } else {
                0.0
            }
        }
        _ => 0.0,
    }
}

fn estimate_stage_runs_per_sec(
    events: &[EventRecord],
    stage: &str,
    started_message: &str,
    completed_message: &str,
) -> f64 {
    let starts = events
        .iter()
        .filter(|event| event.stage == stage && event.message == started_message)
        .count() as f64;
    let completed = events
        .iter()
        .filter(|event| event.stage == stage && event.message == completed_message)
        .count() as f64;
    if starts == 0.0 || completed == 0.0 {
        return 0.0;
    }
    let first_start = earliest_event_ts(events, stage, started_message);
    let last_done = latest_event_ts(events, stage, completed_message);
    match (first_start, last_done) {
        (Some(start), Some(done)) if done > start => {
            let secs = (done - start).num_milliseconds() as f64 / 1000.0;
            if secs > 0.0 { completed / secs } else { 0.0 }
        }
        _ => 0.0,
    }
}

fn estimate_stage_avg_duration_secs(
    events: &[EventRecord],
    stage: &str,
    started_message: &str,
    completed_message: &str,
) -> Option<f64> {
    let start = latest_event_ts(events, stage, started_message)?;
    let done = latest_event_ts(events, stage, completed_message)?;
    if done <= start {
        return None;
    }
    Some((done - start).num_milliseconds() as f64 / 1000.0)
}

fn duration_from_rate(work_left: f64, rate_per_sec: f64) -> Option<std::time::Duration> {
    if work_left <= 0.0 {
        return Some(std::time::Duration::from_secs(0));
    }
    if rate_per_sec <= 0.0 {
        return None;
    }
    Some(std::time::Duration::from_secs_f64(work_left / rate_per_sec))
}

fn format_duration(duration: std::time::Duration) -> String {
    let secs = duration.as_secs();
    let hours = secs / 3600;
    let minutes = (secs % 3600) / 60;
    let seconds = secs % 60;
    format!("{hours:02}:{minutes:02}:{seconds:02}")
}

fn distinct_event_runs(
    events: &[EventRecord],
    stage: &str,
    started_message: &str,
    completed_message: &str,
) -> std::collections::BTreeSet<String> {
    events
        .iter()
        .filter(|event| {
            event.stage == stage
                && (event.message == started_message || event.message == completed_message)
        })
        .filter_map(|event| event.run_accession.clone())
        .collect()
}

fn parse_event_ts(ts: &str) -> Option<DateTime<Utc>> {
    DateTime::parse_from_rfc3339(ts)
        .ok()
        .map(|value| value.with_timezone(&Utc))
}

fn find_event_ts(
    events: &[EventRecord],
    stage: &str,
    message: &str,
    run_accession: Option<&str>,
) -> Option<DateTime<Utc>> {
    events
        .iter()
        .find(|event| {
            event.stage == stage
                && event.message == message
                && event.run_accession.as_deref() == run_accession
        })
        .and_then(|event| parse_event_ts(&event.ts))
}

fn earliest_event_ts(events: &[EventRecord], stage: &str, message: &str) -> Option<DateTime<Utc>> {
    events
        .iter()
        .filter(|event| event.stage == stage && event.message == message)
        .filter_map(|event| parse_event_ts(&event.ts))
        .min()
}

fn latest_event_ts(events: &[EventRecord], stage: &str, message: &str) -> Option<DateTime<Utc>> {
    events
        .iter()
        .filter(|event| event.stage == stage && event.message == message)
        .filter_map(|event| parse_event_ts(&event.ts))
        .max()
}

fn parse_output_mode(spec: &str) -> Result<OutputMode> {
    let normalized = spec.trim().to_ascii_lowercase();
    if normalized == "dense" {
        return Ok(OutputMode::Dense);
    }
    if let Some(value) = normalized.strip_prefix("edges:topk=") {
        let k = value
            .parse::<usize>()
            .with_context(|| format!("invalid topk value in output spec '{spec}'"))?;
        return Ok(OutputMode::TopK { k: k.max(1) });
    }
    if let Some(value) = normalized.strip_prefix("edges:threshold=") {
        let min_abs_r = value
            .parse::<f64>()
            .with_context(|| format!("invalid threshold value in output spec '{spec}'"))?;
        return Ok(OutputMode::Threshold { min_abs_r });
    }
    bail!(
        "unsupported output spec '{}'; expected dense | edges:topk=N | edges:threshold=R",
        spec
    )
}

fn materialize_geneset_matrix(
    paths: &ProjectPaths,
    input_matrix: &Path,
    geneset: &str,
) -> Result<PathBuf> {
    let matrix = read_matrix_tsv(input_matrix)?;
    if matrix.row_ids.len() < 2 {
        return Ok(input_matrix.to_path_buf());
    }
    let selected_indices = parse_geneset(&matrix, geneset)?;
    if selected_indices.len() < 2 {
        bail!(
            "geneset '{}' selected fewer than 2 genes ({})",
            geneset,
            selected_indices.len()
        );
    }
    if selected_indices.len() == matrix.row_ids.len() {
        return Ok(input_matrix.to_path_buf());
    }
    let selected = MatrixData {
        row_ids: selected_indices
            .iter()
            .map(|&idx| matrix.row_ids[idx].clone())
            .collect(),
        col_ids: matrix.col_ids.clone(),
        values: selected_indices
            .iter()
            .map(|&idx| matrix.values[idx].clone())
            .collect(),
    };
    let out_path = paths.corr_dir.join("selected_matrix.tsv");
    write_matrix_tsv(&out_path, &selected)?;
    Ok(out_path)
}

fn parse_geneset(matrix: &MatrixData, geneset: &str) -> Result<Vec<usize>> {
    let spec = geneset.trim();
    if spec.is_empty() || spec.eq_ignore_ascii_case("all") {
        return Ok((0..matrix.row_ids.len()).collect());
    }
    if let Some(value) = spec.strip_prefix("topvar:") {
        let n = value
            .trim()
            .parse::<usize>()
            .with_context(|| format!("invalid topvar geneset: '{geneset}'"))?;
        return Ok(select_topvar(matrix, n));
    }
    if let Some(value) = spec.strip_prefix("mean_gt:") {
        let threshold = value
            .trim()
            .parse::<f64>()
            .with_context(|| format!("invalid mean_gt geneset: '{geneset}'"))?;
        return Ok(matrix
            .values
            .iter()
            .enumerate()
            .filter_map(|(idx, row)| {
                let mean = row.iter().copied().sum::<f64>() / row.len() as f64;
                if mean > threshold { Some(idx) } else { None }
            })
            .collect());
    }
    if let Some(value) = spec.strip_prefix("genes:") {
        let gene_file = PathBuf::from(value.trim());
        let content = fs::read_to_string(&gene_file)
            .with_context(|| format!("failed to read geneset file {}", gene_file.display()))?;
        let gene_set = content
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .collect::<std::collections::BTreeSet<_>>();
        return Ok(matrix
            .row_ids
            .iter()
            .enumerate()
            .filter_map(|(idx, gene)| {
                if gene_set.contains(gene.as_str()) {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect());
    }
    if let Some(value) = spec.strip_prefix("regex:") {
        let pattern = Regex::new(value.trim())
            .with_context(|| format!("invalid regex geneset: '{}'", geneset))?;
        return Ok(matrix
            .row_ids
            .iter()
            .enumerate()
            .filter_map(|(idx, gene)| {
                if pattern.is_match(gene) {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect());
    }
    bail!("unsupported geneset spec '{}'", geneset)
}

fn select_topvar(matrix: &MatrixData, n: usize) -> Vec<usize> {
    let mut idx_var = matrix
        .values
        .iter()
        .enumerate()
        .map(|(idx, row)| {
            let mean = row.iter().copied().sum::<f64>() / row.len() as f64;
            let var = row
                .iter()
                .map(|value| {
                    let d = *value - mean;
                    d * d
                })
                .sum::<f64>()
                / row.len() as f64;
            (idx, var)
        })
        .collect::<Vec<_>>();
    idx_var.sort_by(|left, right| right.1.total_cmp(&left.1));
    idx_var
        .into_iter()
        .take(n.min(matrix.row_ids.len()))
        .map(|(idx, _)| idx)
        .collect()
}

fn cmd_stub(command: &str, root: &Path, detail: &str) -> Result<()> {
    let paths = ProjectPaths::new(root);
    if paths.db_path.exists() {
        let _ = paths.ensure_layout();
        let _log_guard = init_logging(&paths, command)?;
        tracing::warn!("{detail} is not implemented yet");
    }
    bail!("{detail} is scaffolded but not implemented yet")
}

fn resolve_root(root: Option<PathBuf>) -> Result<PathBuf> {
    match root {
        Some(root) => Ok(root),
        None => std::env::current_dir().context("failed to determine current directory"),
    }
}

fn open_project(root: &Path) -> Result<(ProjectPaths, Database, ProjectConfig)> {
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

fn init_logging(
    paths: &ProjectPaths,
    command: &str,
) -> Result<tracing_appender::non_blocking::WorkerGuard> {
    fs::create_dir_all(&paths.logs_dir)
        .with_context(|| format!("failed to create {}", paths.logs_dir.display()))?;
    let file_appender = tracing_appender::rolling::never(&paths.logs_dir, format!("{command}.log"));
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

fn parse_merge_strategy(value: &str) -> Result<MetadataMergeStrategy> {
    match value {
        "override" => Ok(MetadataMergeStrategy::Override),
        "prefer_remote" => Ok(MetadataMergeStrategy::PreferRemote),
        "prefer_local" => Ok(MetadataMergeStrategy::PreferLocal),
        "union-with-prefix" | "union_with_prefix" => Ok(MetadataMergeStrategy::UnionWithPrefix),
        _ => bail!("unsupported merge strategy: {value}"),
    }
}

fn parse_download_preference(value: &str) -> Result<DownloadPreference> {
    match value {
        "fastq" | "ena_fastq" => Ok(DownloadPreference::Fastq),
        "sra" => Ok(DownloadPreference::Sra),
        _ => bail!("unsupported download preference: {value}"),
    }
}

fn parse_cleanup_mode(value: &str) -> Result<CleanupMode> {
    match value {
        "none" => Ok(CleanupMode::None),
        "fastq" => Ok(CleanupMode::Fastq),
        "sra" => Ok(CleanupMode::Sra),
        "all" => Ok(CleanupMode::All),
        _ => bail!("unsupported cleanup mode: {value}"),
    }
}

fn parse_cli_contrasts(items: &[String]) -> Result<Vec<ContrastSpec>> {
    if !items.len().is_multiple_of(3) {
        bail!(
            "contrast values must be provided in groups of 3: --contrast <factor> <level_a> <level_b>"
        );
    }
    let mut contrasts = Vec::new();
    for chunk in items.chunks(3) {
        contrasts.push(ContrastSpec::new(
            None,
            chunk[0].clone(),
            chunk[1].clone(),
            chunk[2].clone(),
        ));
    }
    Ok(contrasts)
}

fn read_samplesheet_header(path: &Path) -> Result<std::collections::BTreeSet<String>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read samplesheet {}", path.display()))?;
    let headers = reader
        .headers()
        .with_context(|| format!("failed to parse samplesheet header {}", path.display()))?;
    Ok(headers.iter().map(ToString::to_string).collect())
}

fn validate_design_formula(
    header: &std::collections::BTreeSet<String>,
    design: &str,
) -> Result<()> {
    let mut missing = Vec::new();
    for token in design
        .split(|ch: char| !ch.is_ascii_alphanumeric() && ch != '_')
        .filter(|token| !token.is_empty())
    {
        if token == "batch" || token == "condition" {
            if !header.contains(token) {
                missing.push(token.to_string());
            }
            continue;
        }
        if token.chars().all(|ch| ch.is_ascii_digit()) {
            continue;
        }
        if token == "Intercept" {
            continue;
        }
        if token != "batch" && token != "condition" && !header.contains(token) {
            missing.push(token.to_string());
        }
    }
    missing.sort();
    missing.dedup();
    if !missing.is_empty() {
        bail!(
            "design formula '{}' references missing samplesheet columns: {}",
            design,
            missing.join(", ")
        );
    }
    Ok(())
}
