use std::collections::VecDeque;
use std::fs;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex, mpsc};
use std::thread;
use std::time::Duration;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result, bail};
use chrono::{DateTime, Utc};
use clap::{Args, Parser, Subcommand};
use regex::Regex;
use rnaa_core::config::{CleanupMode, DownloadPreference, MetadataMergeStrategy, ProjectConfig};
use rnaa_core::model::{
    ArtifactKind, ArtifactRecord, ContrastSpec, CorrelationMethod, EventRecord, InputType,
    OutputMode, QuantArtifacts, RemoteFile, RemoteFileKind, RunRecord, SharedBlobRecord,
    VerifiedFile,
};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::samplesheet::{load_or_init_column_map, write_samplesheet};
use rnaa_core::state::RunState;
use rnaa_core::traits::{
    Correlator, DifferentialExpression, MatrixAdjuster, Preprocessor, Quantifier, ReferenceManager,
};
use rnaa_core::util::{
    compute_sha256, dir_size, file_size, now_rfc3339, sanitize_basename, write_json_pretty,
};
use rnaa_core::{Database, ProjectRecord};
use rnaa_formats::matrix::{MatrixData, read_matrix_tsv, write_matrix_tsv};
use rnaa_stages::corr::{MinCorrCorrelator, ModuleGseaParams, Residualizer, run_module_gsea};
use rnaa_stages::deseq2::Deseq2Runner;
use rnaa_stages::download::ShellDownloader;
use rnaa_stages::preprocess::FasterpInProcessPreprocessor;
use rnaa_stages::quant::{RKallistoQuantifier, reconcile_quant_artifacts};
use rnaa_stages::refs::EnsemblReferenceManager;
use rnaa_stages::resolver::EnaResolver;
use serde_json::json;
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
        Commands::Status => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "status", log_file.as_deref())?;
            cmd_status(&paths, &db, &config)
        }
        Commands::Doctor => {
            let paths = ProjectPaths::new(&root);
            paths.ensure_layout()?;
            let _log_guard = init_logging(&paths, "doctor", log_file.as_deref())?;
            let project = Database::open(&root).ok().and_then(|db| db.project().ok());
            cmd_doctor(&paths, project.as_ref())
        }
        Commands::Download(args) => {
            let (paths, db, mut config) = open_project(&root)?;
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
            let (paths, db, config) = open_project(&root)?;
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
        Commands::Run(args) => {
            let (paths, db, config) = open_project(&root)?;
            let _log_guard = init_logging(&paths, "run", log_file.as_deref())?;
            cmd_run(&paths, &db, config, args)
        }
        Commands::Export(_args) => cmd_stub("export", &root, "export", log_file.as_deref()),
        Commands::Survey(_args) => cmd_stub("survey", &root, "survey", log_file.as_deref()),
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

    #[arg(long, global = true)]
    log_file: Option<PathBuf>,

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
    Normalize(NormalizeArgs),
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
    #[arg(long = "id")]
    ids: Vec<String>,

    #[arg(long = "sample")]
    samples: Vec<String>,

    #[arg(long = "include-sample")]
    include_samples: Vec<String>,

    #[arg(long = "exclude-sample")]
    exclude_samples: Vec<String>,

    #[arg(long = "where")]
    where_predicates: Vec<String>,

    #[arg(long = "exclude-where")]
    exclude_where_predicates: Vec<String>,

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

#[derive(Args, Debug, Clone)]
struct QuantArgs {
    #[arg(long)]
    engine: Option<String>,
    #[arg(long)]
    threads: Option<usize>,
    #[arg(long)]
    workers: Option<usize>,
    #[arg(long, default_value_t = false, conflicts_with = "no_preprocess")]
    preprocess: bool,
    #[arg(long, default_value_t = false, conflicts_with = "preprocess")]
    no_preprocess: bool,
    #[arg(long, default_value_t = false, conflicts_with = "preprocess_bypass")]
    preprocess_strict: bool,
    #[arg(long, default_value_t = false, conflicts_with = "preprocess_strict")]
    preprocess_bypass: bool,
    #[arg(long = "preprocess-max-mb")]
    preprocess_max_mb: Option<usize>,
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
struct NormalizeArgs {
    #[arg(long)]
    design: Option<String>,
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
    #[arg(long = "module-gmt")]
    module_gmt: Option<PathBuf>,
    #[arg(long = "module-top", default_value_t = 500)]
    module_top: usize,
    #[arg(long = "module-min-size", default_value_t = 25)]
    module_min_size: usize,
    #[arg(long = "fgsea-permutations", default_value_t = 1000)]
    fgsea_permutations: usize,
    #[arg(long = "fgsea-seed", default_value_t = 42)]
    fgsea_seed: u64,
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

fn cmd_init(root: &Path, args: InitArgs, log_file: Option<&Path>) -> Result<()> {
    let paths = ProjectPaths::new(root);
    if paths.db_path.exists() {
        bail!("RNAA project already exists at {}", root.display());
    }
    paths.ensure_layout()?;
    let _log_guard = init_logging(&paths, "init", log_file)?;

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
    let mut total_inputs = 0_usize;
    for input_id in &args.ids {
        let input_type = InputType::from_accession(input_id)?;
        db.add_input(input_id, input_type)?;
        total_inputs += 1;
    }
    for sample in &args.samples {
        db.add_input(sample, InputType::Sample)?;
        total_inputs += 1;
    }
    for sample in &args.include_samples {
        db.add_input(sample, InputType::SampleInclude)?;
        total_inputs += 1;
    }
    for sample in &args.exclude_samples {
        db.add_input(sample, InputType::SampleExclude)?;
        total_inputs += 1;
    }
    for predicate in &args.where_predicates {
        validate_where_predicate(predicate)?;
        db.add_input(predicate, InputType::FilterInclude)?;
        total_inputs += 1;
    }
    for predicate in &args.exclude_where_predicates {
        validate_where_predicate(predicate)?;
        db.add_input(predicate, InputType::FilterExclude)?;
        total_inputs += 1;
    }
    if total_inputs == 0 && args.metadata.is_empty() {
        bail!(
            "nothing to add; provide at least one --id/--sample/--include-sample/--exclude-sample/--where/--exclude-where or --metadata"
        );
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

    println!("added_inputs\t{}", total_inputs);
    println!("added_accessions\t{}", args.ids.len());
    println!("added_samples\t{}", args.samples.len());
    println!("added_include_samples\t{}", args.include_samples.len());
    println!("added_exclude_samples\t{}", args.exclude_samples.len());
    println!("added_where_predicates\t{}", args.where_predicates.len());
    println!(
        "added_exclude_where_predicates\t{}",
        args.exclude_where_predicates.len()
    );
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
    let _resolved = resolver.resolve_project(db, paths, config)?;
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
    let progress = compute_progress_snapshot(&runs, &artifacts, &events, config);

    println!("project_id\t{}", project.project_id);
    println!("root\t{}", project.root_dir);
    println!("inputs\t{}", inputs.len());
    println!("runs_total\t{}", progress.total_runs);
    println!("metadata_overrides\t{}", overrides.len());
    println!("samplesheet\t{}", paths.samplesheet_path().display());
    println!(
        "disk_usage_bytes\t{}",
        disk_usage.max(dir_size(&paths.raw_dir))
    );
    println!(
        "progress_download\t{}/{} ({:.1}%)",
        progress.download_done,
        progress.total_runs,
        ratio_pct(progress.download_done, progress.total_runs)
    );
    println!(
        "progress_quant\t{}/{} ({:.1}%)",
        progress.quant_done,
        progress.total_runs,
        ratio_pct(progress.quant_done, progress.total_runs)
    );
    println!(
        "progress_deseq2\t{}/1 ({:.1}%)",
        progress.normalize_done as u8,
        if progress.normalize_done { 100.0 } else { 0.0 }
    );
    println!(
        "progress_corr\t{}/1 ({:.1}%)",
        progress.corr_stage_done as u8,
        if progress.corr_stage_done { 100.0 } else { 0.0 }
    );
    println!("downloaded_bytes\t{}", progress.downloaded_bytes);
    println!(
        "expected_total_download_bytes\t{}",
        progress.expected_total_bytes
    );
    if progress.download_rate_bps > 0.0 {
        println!(
            "download_speed_mbps\t{:.2}",
            progress.download_rate_bps / (1024.0 * 1024.0)
        );
    } else {
        println!("download_speed_mbps\tunknown");
    }
    if progress.quant_rate_runs_per_sec > 0.0 {
        println!(
            "quant_rate_runs_per_hour\t{:.2}",
            progress.quant_rate_runs_per_sec * 3600.0
        );
    } else {
        println!("quant_rate_runs_per_hour\tunknown");
    }
    println!(
        "eta_download\t{}",
        format_eta(
            progress.download_eta,
            progress.download_done >= progress.total_runs
        )
    );
    println!(
        "eta_processing\t{}",
        format_eta(progress.processing_eta, progress.active_stage == "complete")
    );
    println!(
        "eta_total\t{}",
        format_eta(progress.total_eta, progress.active_stage == "complete")
    );
    println!("limiting_stage\t{}", progress.limiting_stage);
    println!("active\t{}", progress.active_detail);

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
            blob_id: None,
            shared_path: None,
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
            blob_id: None,
            shared_path: None,
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
        if threads == 0 {
            bail!("--threads must be >= 1");
        }
        config.quant.threads = threads;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(workers) = args.workers {
        config.quant.workers = workers;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if args.preprocess {
        config.quant.preprocess = true;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if args.no_preprocess {
        config.quant.preprocess = false;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if args.preprocess_strict {
        config.quant.preprocess_strict = true;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if args.preprocess_bypass {
        config.quant.preprocess_strict = false;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
    }
    if let Some(preprocess_max_mb) = args.preprocess_max_mb {
        config.quant.preprocess_max_input_mb = preprocess_max_mb;
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

    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, &config)?;
    reconcile_quant_state(db, paths, &config, &reference)?;
    let runs = db.list_runs_in_states(&[RunState::Verified, RunState::QuantFailed])?;
    if runs.is_empty() {
        println!("quant\tno_verified_runs");
        return Ok(());
    }
    let worker_count = resolve_quant_worker_count(&config, runs.len());
    let mut completed = 0_u64;
    let mut failed = 0_u64;
    if worker_count <= 1 {
        for run in runs {
            if process_quant_run(db, paths, &config, &reference, run)? {
                completed += 1;
            } else {
                failed += 1;
            }
        }
    } else {
        println!("quant_workers\t{}", worker_count);
        let queue = Arc::new(Mutex::new(VecDeque::from(runs)));
        let (tx, rx) = mpsc::channel::<Result<bool>>();
        let mut handles = Vec::new();

        for _ in 0..worker_count {
            let db = db.clone();
            let paths = paths.clone();
            let config = config.clone();
            let reference = reference.clone();
            let queue = Arc::clone(&queue);
            let tx = tx.clone();
            handles.push(thread::spawn(move || {
                loop {
                    let maybe_run = {
                        let mut guard = queue.lock().expect("quant queue mutex poisoned");
                        guard.pop_front()
                    };
                    let Some(run) = maybe_run else {
                        break;
                    };
                    let result = process_quant_run(&db, &paths, &config, &reference, run);
                    let _ = tx.send(result);
                }
            }));
        }
        drop(tx);

        for result in rx {
            match result {
                Ok(true) => completed += 1,
                Ok(false) => failed += 1,
                Err(err) => {
                    failed += 1;
                    tracing::error!("quant worker failed: {err:#}");
                }
            }
        }
        for handle in handles {
            if let Err(err) = handle.join() {
                bail!("quant worker thread panicked: {err:?}");
            }
        }
    }

    println!("quant_completed\t{}", completed);
    println!("quant_failed\t{}", failed);
    Ok(())
}

fn process_quant_run(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::model::ReferenceBundle,
    run: RunRecord,
) -> Result<bool> {
    let quantifier = RKallistoQuantifier;
    let preprocessor = FasterpInProcessPreprocessor;

    if !config.quant.preprocess
        && let Some(artifacts) = reconcile_quant_artifacts(&run, reference, paths, config)?
    {
        finalize_quant_success(
            db,
            config,
            &run.run_accession,
            artifacts,
            false,
            "quantification reused from existing outputs",
        )?;
        return Ok(true);
    }

    db.set_run_state(&run.run_accession, RunState::QuantRunning, None)?;
    db.append_event(
        "quant",
        Some(&run.run_accession),
        "quantification started",
        serde_json::json!({
            "engine": config.quant.engine,
            "threads": config.quant.threads,
            "workers": config.quant.workers,
            "preprocess": config.quant.preprocess,
            "preprocess_strict": config.quant.preprocess_strict
        }),
    )?;
    let raw_fastqs = load_fastq_artifacts(db, &run.run_accession)?;
    if raw_fastqs.is_empty() {
        let message = format!("no FASTQ artifacts found for {}", run.run_accession);
        db.set_run_state(&run.run_accession, RunState::QuantFailed, Some(&message))?;
        db.append_event(
            "quant",
            Some(&run.run_accession),
            "quantification failed",
            serde_json::json!({ "error": message }),
        )?;
        return Ok(false);
    }

    let mut quant_fastqs = raw_fastqs.clone();
    if config.quant.preprocess {
        match preprocessor.preprocess(&run, &raw_fastqs, paths, config) {
            Ok(preprocessed) => {
                for item in &preprocessed.fastqs {
                    db.record_artifact(&ArtifactRecord {
                        project_id: db.project_id().to_string(),
                        run_accession: Some(run.run_accession.clone()),
                        kind: item.artifact_kind,
                        path: item.path.display().to_string(),
                        blob_id: None,
                        shared_path: None,
                        checksum_type: item.checksum_type.clone(),
                        checksum: item.checksum.clone(),
                        bytes: item.bytes,
                        created_at: now_rfc3339(),
                    })?;
                }
                db.record_artifact(&ArtifactRecord {
                    project_id: db.project_id().to_string(),
                    run_accession: Some(run.run_accession.clone()),
                    kind: ArtifactKind::PreprocessReport,
                    path: preprocessed.report_json.display().to_string(),
                    blob_id: None,
                    shared_path: None,
                    checksum_type: "sha256".to_string(),
                    checksum: compute_sha256(&preprocessed.report_json)?,
                    bytes: file_size(&preprocessed.report_json)?,
                    created_at: now_rfc3339(),
                })?;
                db.append_event(
                    "quant",
                    Some(&run.run_accession),
                    "preprocessing completed",
                    serde_json::json!({
                        "tool": preprocessed.tool_name,
                        "version": preprocessed.tool_version,
                        "reused": preprocessed.reused,
                        "passed_reads": preprocessed.passed_reads,
                        "failed_reads": preprocessed.failed_reads,
                        "report_json": preprocessed.report_json.display().to_string(),
                        "outputs": preprocessed
                            .fastqs
                            .iter()
                            .map(|item| item.path.display().to_string())
                            .collect::<Vec<_>>()
                    }),
                )?;
                quant_fastqs = preprocessed.fastqs;
            }
            Err(err) if !config.quant.preprocess_strict => {
                db.append_event(
                    "quant",
                    Some(&run.run_accession),
                    "preprocessing failed; using raw FASTQ",
                    serde_json::json!({ "error": format!("{err:#}") }),
                )?;
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
                return Ok(false);
            }
        }
    }

    match quantifier.quantify(&run, &quant_fastqs, reference, paths, config) {
        Ok(artifacts) => {
            finalize_quant_success(
                db,
                config,
                &run.run_accession,
                artifacts,
                config.quant.preprocess,
                "quantification completed",
            )?;
            Ok(true)
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
            Ok(false)
        }
    }
}

fn reconcile_quant_state(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::model::ReferenceBundle,
) -> Result<()> {
    if config.quant.preprocess {
        return Ok(());
    }

    let runs = db.list_runs_in_states(&[
        RunState::Verified,
        RunState::QuantRunning,
        RunState::QuantFailed,
    ])?;
    for run in runs {
        if let Some(artifacts) = reconcile_quant_artifacts(&run, reference, paths, config)? {
            finalize_quant_success(
                db,
                config,
                &run.run_accession,
                artifacts,
                false,
                "quantification reconciled from existing outputs",
            )?;
        }
    }
    Ok(())
}

fn finalize_quant_success(
    db: &Database,
    config: &ProjectConfig,
    run_accession: &str,
    artifacts: QuantArtifacts,
    preprocessed: bool,
    event_message: &str,
) -> Result<()> {
    let abundance_checksum = compute_sha256(&artifacts.abundance_h5)?;
    let abundance_bytes = file_size(&artifacts.abundance_h5)?;
    let shared_link = register_shared_blob(
        db,
        config,
        &artifacts.abundance_h5,
        &abundance_checksum,
        abundance_bytes,
    )?;
    let records = vec![
        ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: Some(run_accession.to_string()),
            kind: ArtifactKind::QuantDir,
            path: artifacts.out_dir.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "none".to_string(),
            checksum: String::new(),
            bytes: dir_size(&artifacts.out_dir),
            created_at: now_rfc3339(),
        },
        ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: Some(run_accession.to_string()),
            kind: ArtifactKind::QuantAbundanceH5,
            path: artifacts.abundance_h5.display().to_string(),
            blob_id: shared_link.as_ref().map(|link| link.blob_id.clone()),
            shared_path: shared_link
                .as_ref()
                .map(|link| link.shared_path.display().to_string()),
            checksum_type: "sha256".to_string(),
            checksum: abundance_checksum,
            bytes: abundance_bytes,
            created_at: now_rfc3339(),
        },
        ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: Some(run_accession.to_string()),
            kind: ArtifactKind::QuantRunInfo,
            path: artifacts.run_info_json.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".to_string(),
            checksum: compute_sha256(&artifacts.run_info_json)?,
            bytes: file_size(&artifacts.run_info_json)?,
            created_at: now_rfc3339(),
        },
        ArtifactRecord {
            project_id: db.project_id().to_string(),
            run_accession: Some(run_accession.to_string()),
            kind: ArtifactKind::Manifest,
            path: artifacts.manifest_path.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".to_string(),
            checksum: compute_sha256(&artifacts.manifest_path)?,
            bytes: file_size(&artifacts.manifest_path)?,
            created_at: now_rfc3339(),
        },
    ];
    for record in records {
        db.record_artifact(&record)?;
    }
    db.set_run_state(run_accession, RunState::QuantDone, None)?;
    db.append_event(
        "quant",
        Some(run_accession),
        event_message,
        serde_json::json!({
            "out_dir": artifacts.out_dir.display().to_string(),
            "abundance_h5": artifacts.abundance_h5.display().to_string(),
            "run_info_json": artifacts.run_info_json.display().to_string(),
            "shared_blob_id": shared_link.as_ref().map(|link| link.blob_id.clone()),
            "shared_blob_path": shared_link
                .as_ref()
                .map(|link| link.shared_path.display().to_string()),
            "preprocessed": preprocessed
        }),
    )?;
    Ok(())
}

fn resolve_quant_worker_count(config: &ProjectConfig, runs_total: usize) -> usize {
    if runs_total == 0 {
        return 1;
    }
    let per_run_threads = config.quant.threads.max(1);
    let auto = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .saturating_div(per_run_threads)
        .max(1);
    let requested = if config.quant.workers == 0 {
        auto
    } else {
        config.quant.workers
    };
    requested.max(1).min(runs_total)
}

fn cmd_deseq2(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: Deseq2Args,
) -> Result<()> {
    let contrasts = if args.contrast.is_empty() {
        db.list_contrasts()?
    } else {
        parse_cli_contrasts(&args.contrast)?
    };
    let (_outputs, _recorded) =
        run_deseq2_stage(paths, db, config, args.design, contrasts, "deseq2", true)?;
    Ok(())
}

fn cmd_normalize(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: NormalizeArgs,
) -> Result<()> {
    let (outputs, recorded) = run_deseq2_stage(
        paths,
        db,
        config,
        args.design,
        Vec::new(),
        "normalize",
        false,
    )?;
    println!("normalize_outputs\t{}", outputs);
    println!("normalize_recorded\t{}", recorded);
    Ok(())
}

fn run_deseq2_stage(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    design_override: Option<String>,
    contrasts: Vec<ContrastSpec>,
    stage_name: &str,
    store_contrasts: bool,
) -> Result<(u64, u64)> {
    let runs = db.list_runs_in_states(&[
        RunState::QuantDone,
        RunState::DeRunning,
        RunState::DeFailed,
        RunState::DeDone,
        RunState::CorrRunning,
        RunState::CorrFailed,
        RunState::CorrDone,
        RunState::Cleaned,
    ])?;
    if runs.is_empty() {
        println!("{stage_name}\tno_quantified_runs");
        return Ok((0, 0));
    }

    let samplesheet_path = paths.samplesheet_path();
    if !samplesheet_path.exists() {
        bail!(
            "samplesheet missing at {}; run `rnaa resolve` first",
            samplesheet_path.display()
        );
    }
    let header = read_samplesheet_header(&samplesheet_path)?;
    let design = design_override.unwrap_or_else(|| config.deseq2.design.clone());
    validate_design_formula(&header, &design)?;
    for contrast in &contrasts {
        if !header.contains(&contrast.factor) {
            bail!(
                "contrast factor '{}' not found in samplesheet columns",
                contrast.factor
            );
        }
    }

    config.deseq2.design = design.clone();
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    if store_contrasts {
        db.replace_contrasts(&contrasts)?;
    }

    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, &config)?;
    let de_request_manifest = paths
        .de_project_dir(db.project_id())
        .join("rnaa_deseq2_request.json");
    if let Some((outputs, recorded)) = reconcile_deseq2_stage(
        paths,
        db,
        &runs,
        &reference,
        &design,
        &contrasts,
        &config,
        store_contrasts,
        &de_request_manifest,
    )? {
        return Ok((outputs, recorded));
    }

    let runner = Deseq2Runner;
    let target_runs = runs
        .iter()
        .map(|run| run.run_accession.clone())
        .collect::<Vec<_>>();
    for run_accession in &target_runs {
        db.set_run_state(run_accession, RunState::DeRunning, None)?;
    }
    db.append_event(
        stage_name,
        None,
        "stage started",
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
                stage_name,
                None,
                "stage failed",
                serde_json::json!({ "error": message }),
            )?;
            return Err(err);
        }
    };

    let mut recorded = 0_u64;
    write_json_pretty(
        &de_request_manifest,
        &json!({
            "schema_version": 1,
            "stage": stage_name,
            "design": design,
            "transform": config.deseq2.transform,
            "counts_from_abundance": config.deseq2.counts_from_abundance,
            "tx2gene": reference.tx2gene_path.display().to_string(),
            "gene_annotation": reference
                .gene_annotation_path
                .as_ref()
                .map(|path| path.display().to_string()),
            "store_contrasts": store_contrasts,
            "contrasts": contrasts,
        }),
    )?;
    recorded += register_deseq2_outputs(db, &outputs)?;
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().to_string(),
        run_accession: None,
        kind: ArtifactKind::Manifest,
        path: de_request_manifest.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(&de_request_manifest)?,
        bytes: file_size(&de_request_manifest)?,
        created_at: now_rfc3339(),
    })?;
    recorded += 1;

    mark_runs_de_done(db, &runs)?;
    db.append_event(
        stage_name,
        None,
        "stage completed",
        serde_json::json!({
            "outputs_recorded": recorded,
            "runs_marked": target_runs.len(),
        }),
    )?;

    Ok((outputs.len() as u64, recorded))
}

fn reconcile_deseq2_stage(
    paths: &ProjectPaths,
    db: &Database,
    runs: &[RunRecord],
    reference: &rnaa_core::model::ReferenceBundle,
    design: &str,
    contrasts: &[ContrastSpec],
    config: &ProjectConfig,
    store_contrasts: bool,
    request_manifest_path: &Path,
) -> Result<Option<(u64, u64)>> {
    if !request_manifest_path.exists() || file_size(request_manifest_path).unwrap_or_default() == 0
    {
        return Ok(None);
    }
    let expected = json!({
        "stage": if store_contrasts { "deseq2" } else { "normalize" },
        "design": design,
        "transform": config.deseq2.transform,
        "counts_from_abundance": config.deseq2.counts_from_abundance,
        "tx2gene": reference.tx2gene_path.display().to_string(),
        "gene_annotation": reference
            .gene_annotation_path
            .as_ref()
            .map(|path| path.display().to_string()),
        "store_contrasts": store_contrasts,
        "contrasts": contrasts,
    });
    let actual: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(request_manifest_path)
            .with_context(|| format!("failed to read {}", request_manifest_path.display()))?,
    )
    .with_context(|| format!("failed to parse {}", request_manifest_path.display()))?;
    for key in [
        "stage",
        "design",
        "transform",
        "counts_from_abundance",
        "tx2gene",
        "gene_annotation",
        "store_contrasts",
        "contrasts",
    ] {
        if actual.get(key) != expected.get(key) {
            return Ok(None);
        }
    }

    let mut required = vec![
        paths
            .de_project_dir(db.project_id())
            .join("gene_counts.tsv"),
        paths
            .de_project_dir(db.project_id())
            .join("gene_counts_annotated.tsv"),
        paths
            .de_project_dir(db.project_id())
            .join("gene_norm_counts.tsv"),
        paths
            .de_project_dir(db.project_id())
            .join("gene_norm_counts_annotated.tsv"),
        paths.de_project_dir(db.project_id()).join("vst.tsv"),
        paths
            .de_project_dir(db.project_id())
            .join("vst_annotated.tsv"),
        paths
            .de_project_dir(db.project_id())
            .join("deseq2_manifest.json"),
    ];
    if paths
        .de_project_dir(db.project_id())
        .join("vst.rds")
        .exists()
    {
        required.push(paths.de_project_dir(db.project_id()).join("vst.rds"));
    }
    for contrast in contrasts {
        required.push(
            paths
                .de_project_dir(db.project_id())
                .join(format!("de_{}.tsv", contrast.name)),
        );
        required.push(
            paths
                .de_project_dir(db.project_id())
                .join(format!("de_{}_annotated.tsv", contrast.name)),
        );
    }
    if required
        .iter()
        .any(|path| !path.exists() || file_size(path).unwrap_or_default() == 0)
    {
        return Ok(None);
    }

    let mut outputs = required.clone();
    outputs.push(request_manifest_path.to_path_buf());
    let mut recorded = register_deseq2_outputs(db, &outputs)?;
    recorded += register_single_artifact(db, ArtifactKind::Manifest, request_manifest_path)?;
    mark_runs_de_done(db, runs)?;
    db.append_event(
        expected["stage"].as_str().unwrap_or("deseq2"),
        None,
        "stage reconciled from existing outputs",
        json!({
            "request_manifest": request_manifest_path.display().to_string(),
            "outputs_recorded": recorded,
            "runs_marked": runs.len(),
        }),
    )?;
    Ok(Some((outputs.len() as u64, recorded)))
}

fn register_deseq2_outputs(db: &Database, outputs: &[PathBuf]) -> Result<u64> {
    let mut recorded = 0_u64;
    for output in outputs {
        if !output.exists() {
            continue;
        }
        let file_name = output
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or_default();
        let kind = if file_name == "gene_counts.tsv" {
            ArtifactKind::Counts
        } else if file_name == "gene_counts_annotated.tsv" {
            ArtifactKind::Counts
        } else if file_name == "gene_norm_counts.tsv" {
            ArtifactKind::NormalizedCounts
        } else if file_name == "gene_norm_counts_annotated.tsv" {
            ArtifactKind::NormalizedCounts
        } else if file_name == "vst.tsv" {
            ArtifactKind::Vst
        } else if file_name == "vst_annotated.tsv" {
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
        recorded += register_single_artifact(db, kind, output)?;
    }
    Ok(recorded)
}

fn register_corr_outputs(db: &Database, outputs: &[PathBuf]) -> Result<u64> {
    let mut recorded = 0_u64;
    for output in outputs {
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
        } else if file_name.contains("modules_") {
            ArtifactKind::CorrModules
        } else if file_name.contains("module_fgsea") || file_name.contains("module_t_values") {
            ArtifactKind::EnrichmentTable
        } else if file_name.contains("adjusted_matrix") {
            ArtifactKind::AdjustedMatrix
        } else if file_name.ends_with(".json") {
            ArtifactKind::Manifest
        } else {
            continue;
        };
        recorded += register_single_artifact(db, kind, output)?;
    }
    Ok(recorded)
}

fn register_single_artifact(db: &Database, kind: ArtifactKind, path: &Path) -> Result<u64> {
    if !path.exists() {
        return Ok(0);
    }
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().to_string(),
        run_accession: None,
        kind,
        path: path.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(path)?,
        bytes: file_size(path)?,
        created_at: now_rfc3339(),
    })?;
    Ok(1)
}

fn mark_runs_de_done(db: &Database, runs: &[RunRecord]) -> Result<()> {
    for run in runs {
        db.set_run_state(&run.run_accession, RunState::DeDone, None)?;
    }
    Ok(())
}

fn resolve_reference_gene_annotation(
    paths: &ProjectPaths,
    config: &ProjectConfig,
) -> Result<Option<PathBuf>> {
    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, config)?;
    Ok(reference.gene_annotation_path)
}

fn load_gene_annotation_map(
    path: &Path,
) -> Result<std::collections::HashMap<String, (String, String)>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read gene annotation {}", path.display()))?;
    let headers = reader
        .headers()
        .with_context(|| format!("failed to parse gene annotation header {}", path.display()))?
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let gene_id_idx = headers
        .iter()
        .position(|col| col == "gene_id")
        .ok_or_else(|| anyhow::anyhow!("gene annotation missing gene_id column"))?;
    let gene_name_idx = headers
        .iter()
        .position(|col| col == "gene_name")
        .ok_or_else(|| anyhow::anyhow!("gene annotation missing gene_name column"))?;
    let biotype_idx = headers
        .iter()
        .position(|col| col == "gene_biotype")
        .ok_or_else(|| anyhow::anyhow!("gene annotation missing gene_biotype column"))?;
    let mut map = std::collections::HashMap::new();
    for row in reader.records() {
        let row = row?;
        let gene_id = normalize_stable_gene_id(row.get(gene_id_idx).unwrap_or_default());
        if gene_id.is_empty() {
            continue;
        }
        let gene_name = row.get(gene_name_idx).unwrap_or_default().to_string();
        let biotype = row.get(biotype_idx).unwrap_or_default().to_string();
        map.entry(gene_id).or_insert((gene_name, biotype));
    }
    Ok(map)
}

fn normalize_stable_gene_id(value: &str) -> String {
    value.split('.').next().unwrap_or(value).to_string()
}

fn annotate_corr_outputs(outputs: &[PathBuf], annotation_path: &Path) -> Result<Vec<PathBuf>> {
    let annotations = load_gene_annotation_map(annotation_path)?;
    let mut annotated_outputs = Vec::new();
    for output in outputs {
        let Some(file_name) = output.file_name().and_then(|v| v.to_str()) else {
            continue;
        };
        if file_name == "edges_topk.tsv" || file_name == "edges_threshold.tsv" {
            let annotated_path = output.with_file_name(file_name.replace(".tsv", "_annotated.tsv"));
            annotate_edges_table(output, &annotated_path, &annotations)?;
            annotated_outputs.push(annotated_path);
        } else if file_name == "modules_top_genes.tsv" {
            let annotated_path = output.with_file_name("modules_top_genes_annotated.tsv");
            annotate_single_gene_table(
                output,
                &annotated_path,
                "gene",
                &annotations,
                Some(("gene_symbol", "gene_biotype")),
            )?;
            annotated_outputs.push(annotated_path);
        } else if file_name == "module_t_values.tsv" {
            let annotated_path = output.with_file_name("module_t_values_annotated.tsv");
            annotate_single_gene_table(
                output,
                &annotated_path,
                "gene",
                &annotations,
                Some(("gene_symbol", "gene_biotype")),
            )?;
            annotated_outputs.push(annotated_path);
        }
    }
    Ok(annotated_outputs)
}

fn annotate_edges_table(
    input: &Path,
    output: &Path,
    annotations: &std::collections::HashMap<String, (String, String)>,
) -> Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input)
        .with_context(|| format!("failed to read {}", input.display()))?;
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    writer.write_record([
        "gene_a_id",
        "gene_a_symbol",
        "gene_a_biotype",
        "gene_b_id",
        "gene_b_symbol",
        "gene_b_biotype",
        "r",
    ])?;
    for row in reader.records() {
        let row = row?;
        let gene_a = row.get(0).unwrap_or_default();
        let gene_b = row.get(1).unwrap_or_default();
        let r = row.get(2).unwrap_or_default();
        let (gene_a_symbol, gene_a_biotype) = annotations
            .get(gene_a)
            .cloned()
            .unwrap_or_else(|| (String::new(), String::new()));
        let (gene_b_symbol, gene_b_biotype) = annotations
            .get(gene_b)
            .cloned()
            .unwrap_or_else(|| (String::new(), String::new()));
        writer.write_record([
            gene_a,
            gene_a_symbol.as_str(),
            gene_a_biotype.as_str(),
            gene_b,
            gene_b_symbol.as_str(),
            gene_b_biotype.as_str(),
            r,
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn annotate_single_gene_table(
    input: &Path,
    output: &Path,
    gene_column: &str,
    annotations: &std::collections::HashMap<String, (String, String)>,
    extra_columns: Option<(&str, &str)>,
) -> Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input)
        .with_context(|| format!("failed to read {}", input.display()))?;
    let headers = reader.headers()?.clone();
    let gene_idx = headers
        .iter()
        .position(|col| col == gene_column)
        .ok_or_else(|| anyhow::anyhow!("{} missing {} column", input.display(), gene_column))?;
    let mut out_headers = headers.iter().map(ToString::to_string).collect::<Vec<_>>();
    if let Some((symbol_col, biotype_col)) = extra_columns {
        out_headers.insert(gene_idx + 1, symbol_col.to_string());
        out_headers.insert(gene_idx + 2, biotype_col.to_string());
    }
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    writer.write_record(out_headers)?;
    for row in reader.records() {
        let row = row?;
        let gene_id = row.get(gene_idx).unwrap_or_default();
        let (symbol, biotype) = annotations
            .get(gene_id)
            .cloned()
            .unwrap_or_else(|| (String::new(), String::new()));
        let mut out = row.iter().map(ToString::to_string).collect::<Vec<_>>();
        if extra_columns.is_some() {
            out.insert(gene_idx + 1, symbol);
            out.insert(gene_idx + 2, biotype);
        }
        writer.write_record(out)?;
    }
    writer.flush()?;
    Ok(())
}

fn reconcile_corr_stage(
    paths: &ProjectPaths,
    db: &Database,
    runs: &[RunRecord],
    matrix_path: &Path,
    config: &ProjectConfig,
    args: &CorrArgs,
    request_manifest_path: &Path,
) -> Result<Option<(u64, u64)>> {
    if !request_manifest_path.exists() || file_size(request_manifest_path).unwrap_or_default() == 0
    {
        return Ok(None);
    }
    let expected = json!({
        "stage": "corr",
        "matrix": matrix_path.display().to_string(),
        "model": config.corr.model,
        "method": config.corr.method,
        "output": config.corr.output,
        "geneset": config.corr.geneset,
        "gene_annotation": resolve_reference_gene_annotation(paths, config)?
            .map(|path| path.display().to_string()),
        "module_gmt": args.module_gmt.as_ref().map(|path| path.display().to_string()),
        "module_top": args.module_top,
        "module_min_size": args.module_min_size,
        "fgsea_permutations": args.fgsea_permutations,
        "fgsea_seed": args.fgsea_seed,
    });
    let actual: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(request_manifest_path)
            .with_context(|| format!("failed to read {}", request_manifest_path.display()))?,
    )
    .with_context(|| format!("failed to parse {}", request_manifest_path.display()))?;
    for key in [
        "stage",
        "matrix",
        "model",
        "method",
        "output",
        "geneset",
        "gene_annotation",
        "module_gmt",
        "module_top",
        "module_min_size",
        "fgsea_permutations",
        "fgsea_seed",
    ] {
        if actual.get(key) != expected.get(key) {
            return Ok(None);
        }
    }

    let corr_dir = paths.corr_project_dir(db.project_id());
    let mut outputs = vec![
        corr_dir.join("adjusted_matrix.tsv"),
        paths.manifest_path("corr", db.project_id()),
    ];
    match parse_output_mode(&config.corr.output)? {
        OutputMode::TopK { .. } => {
            outputs.push(corr_dir.join("edges_topk.tsv"));
            outputs.push(corr_dir.join("edges_topk_annotated.tsv"));
        }
        OutputMode::Threshold { .. } => {
            outputs.push(corr_dir.join("edges_threshold.tsv"));
            outputs.push(corr_dir.join("edges_threshold_annotated.tsv"));
        }
        OutputMode::Dense => outputs.push(corr_dir.join("corr_dense.tsv")),
    }
    if args.module_gmt.is_some() {
        outputs.push(corr_dir.join("modules_top_genes.tsv"));
        outputs.push(corr_dir.join("modules_top_genes_annotated.tsv"));
        outputs.push(corr_dir.join("module_t_values.tsv"));
        outputs.push(corr_dir.join("module_t_values_annotated.tsv"));
        outputs.push(corr_dir.join("module_fgsea.tsv"));
        outputs.push(corr_dir.join("module_fgsea_manifest.json"));
    }
    if outputs
        .iter()
        .any(|path| !path.exists() || file_size(path).unwrap_or_default() == 0)
    {
        return Ok(None);
    }

    let mut recorded = register_corr_outputs(db, &outputs)?;
    recorded += register_single_artifact(db, ArtifactKind::Manifest, request_manifest_path)?;
    for run in runs {
        db.set_run_state(&run.run_accession, RunState::CorrDone, None)?;
    }
    db.append_event(
        "corr",
        None,
        "corr stage reconciled from existing outputs",
        json!({
            "request_manifest": request_manifest_path.display().to_string(),
            "outputs_recorded": recorded,
            "runs_marked": runs.len(),
        }),
    )?;
    Ok(Some((outputs.len() as u64, recorded)))
}

fn cmd_corr(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    args: CorrArgs,
) -> Result<()> {
    if let Some(adjust) = &args.adjust {
        config.corr.adjust = adjust.clone();
    }
    if let Some(model) = &args.model {
        config.corr.model = model.clone();
    }
    if let Some(method) = &args.method {
        config.corr.method = method.clone();
    }
    if let Some(geneset) = &args.geneset {
        config.corr.geneset = geneset.clone();
    }
    if let Some(out_spec) = &args.out_spec {
        config.corr.output = out_spec.clone();
    }
    if args.module_gmt.is_some() {
        config.corr.method = "spearman".to_string();
    }
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    let reference = EnsemblReferenceManager::new()?.ensure_reference(paths, &config)?;

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

    let corr_dir = paths.corr_project_dir(db.project_id());
    let corr_request_manifest = corr_dir.join("rnaa_corr_request.json");
    if let Some((_outputs, _recorded)) = reconcile_corr_stage(
        paths,
        db,
        &runs,
        &matrix_path,
        &config,
        &args,
        &corr_request_manifest,
    )? {
        return Ok(());
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
            "geneset": config.corr.geneset,
            "gene_annotation": reference
                .gene_annotation_path
                .as_ref()
                .map(|path| path.display().to_string()),
            "module_gmt": args.module_gmt.as_ref().map(|path| path.display().to_string()),
            "module_top": args.module_top,
            "module_min_size": args.module_min_size,
            "fgsea_permutations": args.fgsea_permutations,
        }),
    )?;

    let method = config.corr.method.parse::<CorrelationMethod>()?;
    let output_mode = parse_output_mode(&config.corr.output)?;
    let correlator = MinCorrCorrelator;
    let mut outputs =
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
    if let Some(module_gmt) = args.module_gmt.as_ref() {
        let params = ModuleGseaParams {
            gmt_path: module_gmt.clone(),
            top_genes: args.module_top,
            min_module_size: args.module_min_size,
            permutations: args.fgsea_permutations,
            seed: args.fgsea_seed,
        };
        let module_outputs = run_module_gsea(&adjusted, paths, db.project_id(), &params)
            .with_context(|| {
                format!("module enrichment failed for gmt {}", module_gmt.display())
            })?;
        outputs.extend(module_outputs);
    }
    if let Some(annotation_path) = &reference.gene_annotation_path {
        let annotated = annotate_corr_outputs(&outputs, annotation_path)?;
        outputs.extend(annotated);
    }

    write_json_pretty(
        &corr_request_manifest,
        &json!({
            "schema_version": 1,
            "stage": "corr",
            "matrix": matrix_path.display().to_string(),
            "model": config.corr.model,
            "method": config.corr.method,
            "output": config.corr.output,
            "geneset": config.corr.geneset,
            "gene_annotation": reference
                .gene_annotation_path
                .as_ref()
                .map(|path| path.display().to_string()),
            "module_gmt": args.module_gmt.as_ref().map(|path| path.display().to_string()),
            "module_top": args.module_top,
            "module_min_size": args.module_min_size,
            "fgsea_permutations": args.fgsea_permutations,
            "fgsea_seed": args.fgsea_seed,
        }),
    )?;

    let mut recorded = register_corr_outputs(db, &outputs)?;
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().to_string(),
        run_accession: None,
        kind: ArtifactKind::AdjustedMatrix,
        path: adjusted.matrix_path.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(&adjusted.matrix_path)?,
        bytes: file_size(&adjusted.matrix_path)?,
        created_at: now_rfc3339(),
    })?;
    recorded += 1;
    db.record_artifact(&ArtifactRecord {
        project_id: db.project_id().to_string(),
        run_accession: None,
        kind: ArtifactKind::Manifest,
        path: corr_request_manifest.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(&corr_request_manifest)?,
        bytes: file_size(&corr_request_manifest)?,
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
    Ok(())
}

fn cmd_run(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: RunArgs,
) -> Result<()> {
    cmd_resolve(paths, db, &config, ResolveArgs { force: false })?;
    let quant_args = QuantArgs {
        engine: None,
        threads: None,
        workers: None,
        preprocess: false,
        no_preprocess: false,
        preprocess_strict: false,
        preprocess_bypass: false,
        preprocess_max_mb: None,
        cleanup: None,
        trash_days: None,
    };

    // Pipeline mode: run download and quant together so verified runs are consumed immediately.
    let mut download_handle = if args.no_download {
        None
    } else {
        let downloader = ShellDownloader;
        let db_for_download = db.clone();
        let paths_for_download = paths.clone();
        let config_for_download = config.clone();
        Some(thread::spawn(move || {
            downloader.run_worker(
                db_for_download,
                paths_for_download,
                config_for_download.clone(),
                config_for_download.download.concurrency,
                false,
                None,
            )
        }))
    };

    let started_at = std::time::Instant::now();
    let mut last_progress = started_at;
    let mut last_progress_snapshot = None;
    print_run_progress(
        paths,
        db,
        &config,
        started_at.elapsed(),
        &mut last_progress_snapshot,
    )?;

    loop {
        if last_progress.elapsed() >= Duration::from_secs(10) {
            print_run_progress(
                paths,
                db,
                &config,
                started_at.elapsed(),
                &mut last_progress_snapshot,
            )?;
            last_progress = std::time::Instant::now();
        }
        let quant_ready = db.list_runs_in_states(&[RunState::Verified, RunState::QuantFailed])?;
        let mut did_work = false;
        if !quant_ready.is_empty() {
            cmd_quant(paths, db, config.clone(), quant_args.clone())?;
            last_progress = std::time::Instant::now();
            print_run_progress(
                paths,
                db,
                &config,
                started_at.elapsed(),
                &mut last_progress_snapshot,
            )?;
            did_work = true;
        }

        if let Some(handle) = &download_handle
            && handle.is_finished()
        {
            let join_result = download_handle.take().expect("download handle present");
            let worker_result = join_result
                .join()
                .map_err(|_| anyhow::anyhow!("download worker thread panicked"))?;
            worker_result?;
            refresh_samplesheet(paths, db)?;
            last_progress = std::time::Instant::now();
            print_run_progress(
                paths,
                db,
                &config,
                started_at.elapsed(),
                &mut last_progress_snapshot,
            )?;
            did_work = true;
        }

        let pending_download = if args.no_download {
            false
        } else {
            has_pending_download_runs(db)?
        };
        let pending_quant = !db
            .list_runs_in_states(&[RunState::Verified, RunState::QuantFailed])?
            .is_empty();

        if download_handle.is_none() && !pending_download && !pending_quant {
            break;
        }
        if !did_work {
            thread::sleep(Duration::from_secs(2));
        }
    }

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
            config.clone(),
            CorrArgs {
                matrix: None,
                adjust: None,
                model: None,
                method: None,
                geneset: None,
                out_spec: None,
                module_gmt: None,
                module_top: 500,
                module_min_size: 25,
                fgsea_permutations: 1000,
                fgsea_seed: 42,
            },
        )?;
    }
    print_run_progress(
        paths,
        db,
        &config,
        started_at.elapsed(),
        &mut last_progress_snapshot,
    )?;
    println!("run\tcompleted");
    Ok(())
}

fn print_run_progress(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    elapsed: Duration,
    last_snapshot: &mut Option<String>,
) -> Result<()> {
    let runs = db.list_runs()?;
    let artifacts = db.list_artifacts()?;
    let events = db.list_events(5_000)?;
    let progress = compute_progress_snapshot(&runs, &artifacts, &events, config);
    let eta_total_text = format_eta(progress.total_eta, progress.active_stage == "complete");
    let eta_download_text = format_eta(
        progress.download_eta,
        progress.download_done >= progress.total_runs,
    );
    let eta_processing_text =
        format_eta(progress.processing_eta, progress.active_stage == "complete");
    let mut lines = vec![
        format!("Project\t{}", db.project_id()),
        format!("Elapsed\t{}", format_duration(elapsed)),
        format!(
            "Stage\t{} | download {}/{} | quant {}/{} | normalize {} | corr {}",
            progress.active_stage,
            progress.download_done,
            progress.total_runs,
            progress.quant_done,
            progress.total_runs,
            if progress.normalize_done {
                "1/1"
            } else {
                "0/1"
            },
            if progress.corr_stage_done {
                "1/1"
            } else {
                "0/1"
            }
        ),
        format!(
            "ETA\t{} total | {} download | {} processing | limiting {}",
            eta_total_text, eta_download_text, eta_processing_text, progress.limiting_stage
        ),
        format!("Active\t{}", progress.active_detail),
    ];
    if progress.download_rate_bps > 0.0 {
        lines.push(format!(
            "Speed\t{:.2} MiB/s download",
            progress.download_rate_bps / (1024.0 * 1024.0)
        ));
    }
    lines.push(format!(
        "Paths\t{} | {}",
        paths.samplesheet_path().display(),
        paths.logs_dir.display()
    ));
    let snapshot = lines.join("\n");
    if last_snapshot.as_deref() != Some(snapshot.as_str()) {
        println!("{snapshot}");
        *last_snapshot = Some(snapshot);
    }
    Ok(())
}

fn format_eta(value: Option<Duration>, completed: bool) -> String {
    match value {
        Some(duration) if completed => format_duration(duration),
        Some(duration) if duration > Duration::from_secs(0) => format_duration(duration),
        _ if completed => "00:00:00".to_string(),
        _ => "estimating".to_string(),
    }
}

struct ProgressSnapshot {
    total_runs: u64,
    download_done: u64,
    quant_done: u64,
    normalize_done: bool,
    corr_stage_done: bool,
    downloaded_bytes: u64,
    expected_total_bytes: u64,
    download_rate_bps: f64,
    quant_rate_runs_per_sec: f64,
    download_eta: Option<Duration>,
    processing_eta: Option<Duration>,
    total_eta: Option<Duration>,
    limiting_stage: &'static str,
    active_stage: &'static str,
    active_detail: String,
}

fn compute_progress_snapshot(
    runs: &[RunRecord],
    artifacts: &[ArtifactRecord],
    events: &[EventRecord],
    config: &ProjectConfig,
) -> ProgressSnapshot {
    let total_runs = runs.len() as u64;
    let download_active = runs
        .iter()
        .filter(|run| matches!(run.state, RunState::Downloading))
        .count() as u64;
    let quant_running = runs
        .iter()
        .filter(|run| matches!(run.state, RunState::QuantRunning))
        .count() as u64;
    let quant_done = runs
        .iter()
        .filter(|run| is_quant_done_state(run.state))
        .count() as u64;
    let normalize_done = runs.iter().all(|run| is_de_done_state(run.state));
    let corr_stage_done = runs.iter().all(|run| is_corr_done_state(run.state));
    let download_done = artifacts
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
        .filter_map(|artifact| artifact.run_accession.clone())
        .collect::<std::collections::BTreeSet<_>>()
        .len() as u64;
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
        estimate_total_download_bytes(runs, config.download.prefer).max(downloaded_bytes);
    let remaining_download_bytes = expected_total_bytes.saturating_sub(downloaded_bytes);
    let download_rate_bps = estimate_download_speed_bps(events);
    let quant_rate_runs_per_sec = estimate_stage_runs_per_sec(
        events,
        "quant",
        "quantification started",
        "quantification completed",
    );
    let de_avg_secs = estimate_stage_avg_duration_secs(
        events,
        "deseq2",
        "deseq2 stage started",
        "deseq2 stage completed",
    );
    let corr_avg_secs = estimate_stage_avg_duration_secs(
        events,
        "corr",
        "corr stage started",
        "corr stage completed",
    );
    let download_eta = duration_from_rate(remaining_download_bytes as f64, download_rate_bps);
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
    let quant_eta = duration_from_rate(quant_remaining_runs as f64, quant_rate_runs_per_sec);
    let de_eta = if normalize_done {
        Some(Duration::from_secs(0))
    } else {
        de_avg_secs.map(Duration::from_secs_f64)
    };
    let corr_eta = if corr_stage_done {
        Some(Duration::from_secs(0))
    } else {
        corr_avg_secs.map(Duration::from_secs_f64)
    };
    let processing_eta = match (quant_eta, de_eta, corr_eta) {
        (Some(q), Some(d), Some(c)) => Some(q + d + c),
        (Some(q), _, _) => Some(q),
        (_, Some(d), Some(c)) => Some(d + c),
        (_, Some(d), _) => Some(d),
        (_, _, Some(c)) => Some(c),
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
    let active_stage = if quant_running > 0 {
        "quantification"
    } else if download_active > 0 {
        "download"
    } else if !corr_stage_done && normalize_done {
        "correlation"
    } else if !normalize_done && quant_done == total_runs {
        "normalization"
    } else if corr_stage_done {
        "complete"
    } else {
        "waiting"
    };
    let active_detail = describe_active_work(
        runs,
        total_runs,
        download_active,
        quant_running,
        quant_done,
        normalize_done,
        corr_stage_done,
        active_stage,
    );

    ProgressSnapshot {
        total_runs,
        download_done,
        quant_done,
        normalize_done,
        corr_stage_done,
        downloaded_bytes,
        expected_total_bytes,
        download_rate_bps,
        quant_rate_runs_per_sec,
        download_eta,
        processing_eta,
        total_eta,
        limiting_stage,
        active_stage,
        active_detail,
    }
}

fn describe_active_work(
    runs: &[RunRecord],
    total_runs: u64,
    download_active: u64,
    quant_running: u64,
    quant_done: u64,
    normalize_done: bool,
    corr_done: bool,
    active_stage: &str,
) -> String {
    match active_stage {
        "download" => format!(
            "downloading {download_active} run(s); {}/{} completed",
            runs.iter()
                .filter(|run| matches!(run.state, RunState::Downloaded | RunState::Verified))
                .count(),
            total_runs
        ),
        "quantification" => {
            format!("quantifying {quant_running} run(s); {quant_done}/{total_runs} complete")
        }
        "normalization" => {
            let samples = runs.len();
            if normalize_done {
                format!("normalization complete for {samples} samples")
            } else {
                format!("running DESeq2 normalization on {samples} samples")
            }
        }
        "correlation" => {
            let samples = runs.len();
            if corr_done {
                format!("correlation complete for {samples} samples")
            } else {
                format!("running adjusted correlation on {samples} samples")
            }
        }
        "complete" => format!("all stages complete for {} runs", total_runs),
        _ => {
            let failed = runs
                .iter()
                .filter(|run| {
                    matches!(
                        run.state,
                        RunState::DownloadFailed
                            | RunState::VerifyFailed
                            | RunState::QuantFailed
                            | RunState::DeFailed
                            | RunState::CorrFailed
                    )
                })
                .count();
            if failed > 0 {
                format!("waiting with {failed} failed run(s) needing retry")
            } else {
                "waiting for next eligible work item".to_string()
            }
        }
    }
}

fn has_pending_download_runs(db: &Database) -> Result<bool> {
    Ok(!db
        .list_runs_in_states(&[
            RunState::Resolved,
            RunState::Downloading,
            RunState::Downloaded,
            RunState::DownloadFailed,
            RunState::VerifyFailed,
        ])?
        .is_empty())
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

fn cmd_stub(command: &str, root: &Path, detail: &str, log_file: Option<&Path>) -> Result<()> {
    let paths = ProjectPaths::new(root);
    if paths.db_path.exists() {
        let _ = paths.ensure_layout();
        let _log_guard = init_logging(&paths, command, log_file)?;
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

#[derive(Debug, Clone)]
struct SharedBlobLink {
    blob_id: String,
    shared_path: PathBuf,
}

fn register_shared_blob(
    db: &Database,
    config: &ProjectConfig,
    source_path: &Path,
    checksum: &str,
    bytes: u64,
) -> Result<Option<SharedBlobLink>> {
    let Some(shared_root) = resolve_shared_root(config) else {
        return Ok(None);
    };
    if checksum.len() < 4 {
        bail!("invalid sha256 checksum for shared blob: {checksum}");
    }

    let blob_id = format!("sha256:{checksum}");
    let shared_path = shared_root
        .join("blobs")
        .join("sha256")
        .join(&checksum[0..2])
        .join(&checksum[2..4])
        .join(checksum);

    materialize_shared_blob(source_path, &shared_path)?;

    let shared_checksum = compute_sha256(&shared_path)?;
    if shared_checksum != checksum {
        bail!(
            "shared blob checksum mismatch for {} (expected {}, got {})",
            shared_path.display(),
            checksum,
            shared_checksum
        );
    }
    let shared_bytes = file_size(&shared_path)?;
    if shared_bytes != bytes {
        bail!(
            "shared blob size mismatch for {} (expected {}, got {})",
            shared_path.display(),
            bytes,
            shared_bytes
        );
    }

    db.record_shared_blob(&SharedBlobRecord {
        blob_id: blob_id.clone(),
        storage_path: shared_path.display().to_string(),
        checksum_type: "sha256".to_string(),
        checksum: checksum.to_string(),
        bytes,
        created_at: now_rfc3339(),
    })?;

    Ok(Some(SharedBlobLink {
        blob_id,
        shared_path,
    }))
}

fn resolve_shared_root(config: &ProjectConfig) -> Option<PathBuf> {
    let configured = config.storage.shared_root.trim();
    if !configured.is_empty() {
        return Some(PathBuf::from(configured));
    }
    std::env::var("RNAA_SHARED_ROOT").ok().and_then(|value| {
        let trimmed = value.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(PathBuf::from(trimmed))
        }
    })
}

fn materialize_shared_blob(source_path: &Path, shared_path: &Path) -> Result<()> {
    if shared_path.exists() {
        return Ok(());
    }
    let parent = shared_path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("invalid shared blob path {}", shared_path.display()))?;
    fs::create_dir_all(parent).with_context(|| format!("failed to create {}", parent.display()))?;

    match fs::hard_link(source_path, shared_path) {
        Ok(()) => return Ok(()),
        Err(err) if err.kind() == ErrorKind::AlreadyExists => return Ok(()),
        Err(_) => {}
    }

    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system clock before UNIX_EPOCH")?
        .as_nanos();
    let temp_path = shared_path.with_extension(format!("tmp-{stamp}"));
    fs::copy(source_path, &temp_path).with_context(|| {
        format!(
            "failed to copy {} into shared store temp {}",
            source_path.display(),
            temp_path.display()
        )
    })?;
    match fs::rename(&temp_path, shared_path) {
        Ok(()) => Ok(()),
        Err(err) if err.kind() == ErrorKind::AlreadyExists => {
            let _ = fs::remove_file(&temp_path);
            Ok(())
        }
        Err(err) => {
            let _ = fs::remove_file(&temp_path);
            Err(err).with_context(|| {
                format!(
                    "failed to finalize shared blob {} -> {}",
                    temp_path.display(),
                    shared_path.display()
                )
            })
        }
    }
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

fn validate_where_predicate(value: &str) -> Result<()> {
    let trimmed = value.trim();
    let valid = trimmed.contains("!=") || trimmed.contains('=');
    if !valid {
        bail!(
            "invalid where predicate '{}'; expected 'field=value' or 'field!=value'",
            value
        );
    }
    Ok(())
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
