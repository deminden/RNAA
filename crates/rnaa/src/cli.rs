//! CLI arguments and conversion into typed workflow settings.
use anyhow::{Result, bail};
use clap::{Args, Parser, Subcommand};
use rnaa_core::config::{
    CleanupMode, CleanupWhen, DownloadPreference, FastqRetention, MetadataMergeStrategy,
};
use rnaa_core::model::ContrastSpec;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "rnaa",
    version,
    about = "RNAA - Rust Nucleotide Alignment and Analytics"
)]
pub(crate) struct Cli {
    #[arg(long, global = true)]
    pub(crate) root: Option<PathBuf>,

    #[arg(long, global = true)]
    pub(crate) log_file: Option<PathBuf>,

    #[command(subcommand)]
    pub(crate) command: Commands,
}

#[derive(Subcommand, Debug)]
pub(crate) enum Commands {
    Init(InitArgs),
    Add(AddArgs),
    Resolve(ResolveArgs),
    Metadata(MetadataArgs),
    Lineage(LineageArgs),
    Refs(RefsArgs),
    Download(DownloadArgs),
    Quant(QuantArgs),
    Normalize(NormalizeArgs),
    Deseq2(Deseq2Args),
    Corr(CorrArgs),
    Qc(QcArgs),
    Run(RunArgs),
    Status {
        /// Emit a machine-readable project status report.
        #[arg(long)]
        json: bool,
    },
    Storage(StorageArgs),
    Jobs(JobsArgs),
    Doctor {
        /// Emit a machine-readable diagnostic report.
        #[arg(long)]
        json: bool,
    },
    Export(ExportArgs),
    Survey(SurveyArgs),
}

#[derive(Args, Debug)]
pub(crate) struct StorageArgs {
    /// Override the configured shared store location.
    #[arg(long, global = true)]
    pub(crate) shared_root: Option<PathBuf>,
    #[command(subcommand)]
    pub(crate) command: StorageCommand,
}
#[derive(Subcommand, Debug)]
pub(crate) enum StorageCommand {
    /// Measure project bytes and distinguish incomplete generations and trash.
    Usage {
        #[arg(long)]
        include_shared: bool,
        #[arg(long)]
        details: bool,
    },
    /// Inspect, trash, restore or expire count-backed FASTQ/SRA inputs.
    Files {
        #[command(subcommand)]
        command: FileStorageCommand,
    },
    /// Inspect and collect unreachable project generations.
    Generations {
        #[command(subcommand)]
        command: GenerationStorageCommand,
    },
    /// Show store-wide catalog ownership and durable pins as JSON.
    Status,
    /// Import ownership from an explicit inventory of idle legacy projects.
    ImportLegacy {
        #[arg(long, required = true)]
        project: Vec<PathBuf>,
        #[arg(long)]
        apply: bool,
        /// Assert that every legacy project using this store is supplied.
        #[arg(long)]
        complete_project_inventory: bool,
    },
    /// Preview collection; --apply moves unowned blobs to reversible trash.
    Gc {
        #[arg(long)]
        apply: bool,
        /// Also delete trash whose original retention deadline has passed.
        #[arg(long)]
        purge_expired: bool,
        /// Retention for newly trashed blobs; cannot shorten existing deadlines.
        #[arg(long)]
        trash_days: Option<u64>,
    },
    /// Restore a blob by SHA-256 or sha256:<digest> identity.
    Restore { blob: String },
}

#[derive(Subcommand, Debug)]
pub(crate) enum FileStorageCommand {
    Status,
    /// Preview input cleanup. Explicit --apply moves complete run groups to trash.
    Clean {
        #[arg(long, value_parser = ["none", "fastq", "sra", "all"])]
        mode: Option<String>,
        #[arg(long)]
        run: Option<String>,
        #[arg(long)]
        apply: bool,
        #[arg(long)]
        trash_days: Option<u64>,
    },
    /// Preview expiration of existing trash without selecting new input files.
    Gc {
        #[arg(long)]
        apply: bool,
        #[arg(long)]
        purge_expired: bool,
        #[arg(long)]
        run: Option<String>,
    },
    /// Restore the complete run operation containing this file trash identity.
    Restore {
        trash_id: String,
    },
}

#[derive(Subcommand, Debug)]
pub(crate) enum GenerationStorageCommand {
    Status,
    /// Reversible cleanup of tracked unpublished data from stopped jobs.
    Unpublished {
        #[command(subcommand)]
        command: UnpublishedStorageCommand,
    },
    /// Preview rebuilding missing ancestry from immutable legacy evidence.
    ImportLineage {
        #[arg(long)]
        apply: bool,
    },
    Gc {
        #[arg(long)]
        apply: bool,
        #[arg(long)]
        purge_expired: bool,
        #[arg(long)]
        trash_days: Option<u64>,
    },
    Restore {
        generation: String,
    },
    Keep {
        generation: String,
        #[arg(long, default_value = "manual")]
        label: String,
    },
    Unkeep {
        generation: String,
    },
}

#[derive(Subcommand, Debug)]
pub(crate) enum UnpublishedStorageCommand {
    Status,
    Gc {
        #[arg(long)]
        apply: bool,
        #[arg(long)]
        purge_expired: bool,
        #[arg(long)]
        trash_days: Option<u64>,
    },
    /// Restore diagnostic files; the generation remains unpublished.
    Restore {
        generation: String,
    },
}

#[derive(Args, Debug)]
pub(crate) struct JobsArgs {
    #[command(subcommand)]
    pub(crate) command: JobsCommand,
}

#[derive(Subcommand, Debug)]
pub(crate) enum JobsCommand {
    List {
        #[arg(long)]
        json: bool,
    },
    Cancel {
        attempt: String,
    },
}

#[derive(Args, Debug)]
pub(crate) struct InitArgs {
    #[arg(long)]
    pub(crate) config: Option<PathBuf>,
}

#[derive(Args, Debug)]
pub(crate) struct AddArgs {
    #[arg(long = "id")]
    pub(crate) ids: Vec<String>,

    #[arg(long = "sample")]
    pub(crate) samples: Vec<String>,

    #[arg(long = "include-sample")]
    pub(crate) include_samples: Vec<String>,

    #[arg(long = "exclude-sample")]
    pub(crate) exclude_samples: Vec<String>,

    #[arg(long = "where")]
    pub(crate) where_predicates: Vec<String>,

    #[arg(long = "exclude-where")]
    pub(crate) exclude_where_predicates: Vec<String>,

    #[arg(long = "metadata")]
    pub(crate) metadata: Vec<PathBuf>,

    #[arg(long, default_value = "prefer_local")]
    pub(crate) merge_strategy: String,
}

#[derive(Args, Debug)]
pub(crate) struct ResolveArgs {
    #[arg(long)]
    pub(crate) force: bool,
    #[arg(long, default_value = "fail", value_parser = ["fail", "keep-local", "use-source"])]
    pub(crate) metadata_conflicts: String,
}

#[derive(Args, Debug)]
pub(crate) struct MetadataArgs {
    #[command(subcommand)]
    pub(crate) command: MetadataCommand,
}

#[derive(Args, Debug)]
pub(crate) struct LineageArgs {
    #[arg(long, value_parser = ["refs", "download", "preprocess", "quant", "de", "corr"])]
    pub(crate) group: Option<String>,
    #[arg(long)]
    pub(crate) subject: Option<String>,
    /// Include the selected generations' immutable upstream history.
    #[arg(long)]
    pub(crate) history: bool,
    /// Include all committed request variants for the selected output subjects.
    #[arg(long)]
    pub(crate) all_requests: bool,
    /// Rehash payloads and validate authorized pruning; may read large files.
    #[arg(long)]
    pub(crate) verify: bool,
}

#[derive(Subcommand, Debug)]
pub(crate) enum MetadataCommand {
    /// Refresh generated fields while preserving deliberate sample metadata edits.
    Refresh {
        #[arg(long, default_value = "fail", value_parser = ["fail", "keep-local", "use-source"])]
        conflicts: String,
    },
    /// Show and verify the current immutable metadata publication.
    Show,
}

pub(crate) fn metadata_conflict_policy(
    value: &str,
) -> Result<rnaa_core::samplesheet::publication::ConflictPolicy> {
    use rnaa_core::samplesheet::publication::ConflictPolicy;
    match value {
        "fail" => Ok(ConflictPolicy::Fail),
        "keep-local" => Ok(ConflictPolicy::KeepLocal),
        "use-source" => Ok(ConflictPolicy::UseSource),
        _ => bail!("unknown metadata conflict policy '{value}'"),
    }
}

#[derive(Args, Debug)]
pub(crate) struct RefsArgs {
    #[command(subcommand)]
    pub(crate) command: RefsCommand,
}

#[derive(Subcommand, Debug)]
pub(crate) enum RefsCommand {
    Prepare(RefsPrepareArgs),
}

#[derive(Args, Debug)]
pub(crate) struct RefsPrepareArgs {
    #[arg(long)]
    pub(crate) organism: Option<String>,
    #[arg(long)]
    pub(crate) ensembl: Option<String>,
    #[arg(long)]
    pub(crate) gtf: Option<PathBuf>,
    #[arg(long)]
    pub(crate) cdna: Option<PathBuf>,
}

#[derive(Args, Debug)]
pub(crate) struct DownloadArgs {
    #[arg(long)]
    pub(crate) concurrency: Option<usize>,
    #[arg(long = "storage-limit")]
    pub(crate) storage_limit: Option<String>,

    #[arg(long)]
    pub(crate) forever: bool,

    #[arg(long)]
    pub(crate) prefer: Option<String>,
}

#[derive(Args, Debug, Clone)]
pub(crate) struct QuantArgs {
    #[arg(long)]
    pub(crate) strand: Option<rnaa_core::config::QuantStrand>,
    #[arg(long)]
    pub(crate) single_fragment_length: Option<f64>,
    #[arg(long)]
    pub(crate) single_fragment_sd: Option<f64>,
    #[arg(long)]
    pub(crate) engine: Option<String>,
    #[arg(long = "worker-budget")]
    pub(crate) worker_budget: Option<usize>,
    #[arg(long)]
    pub(crate) threads: Option<usize>,
    /// Exact-match cache allowance across mapping threads of each run, MiB.
    #[arg(long = "match-cache-mb")]
    pub(crate) match_cache_mb: Option<usize>,
    #[arg(long)]
    pub(crate) workers: Option<usize>,
    #[arg(long, default_value_t = false, conflicts_with = "no_preprocess")]
    pub(crate) preprocess: bool,
    #[arg(long = "fastq-retention")]
    pub(crate) fastq_retention: Option<String>,
    #[arg(long, default_value_t = false, conflicts_with = "preprocess")]
    pub(crate) no_preprocess: bool,
    #[arg(long, default_value_t = false, conflicts_with = "preprocess_bypass")]
    pub(crate) preprocess_strict: bool,
    #[arg(long, default_value_t = false, conflicts_with = "preprocess_strict")]
    pub(crate) preprocess_bypass: bool,
    #[arg(
        long = "preprocess-max-mb",
        help = "Maximum total decoded FASTQ input across preprocessing mates and lanes, MiB; 0 disables this additional cap"
    )]
    pub(crate) preprocess_max_mb: Option<usize>,
    #[arg(long = "preprocess-threads")]
    pub(crate) preprocess_threads: Option<usize>,
    #[arg(long)]
    pub(crate) cleanup: Option<String>,
    #[arg(long, value_parser = ["success", "always", "never"])]
    pub(crate) cleanup_on: Option<String>,
    #[arg(long = "trash-days")]
    pub(crate) trash_days: Option<u64>,
}

#[derive(Args, Debug)]
pub(crate) struct Deseq2Args {
    #[arg(long)]
    pub(crate) design: Option<String>,
    #[arg(long = "contrast", num_args = 3)]
    pub(crate) contrast: Vec<String>,
}

#[derive(Args, Debug)]
pub(crate) struct NormalizeArgs {
    #[arg(long)]
    pub(crate) design: Option<String>,
}

#[derive(Args, Debug)]
pub(crate) struct CorrArgs {
    #[arg(long)]
    pub(crate) matrix: Option<String>,
    #[arg(long)]
    pub(crate) adjust: Option<String>,
    #[arg(long)]
    pub(crate) model: Option<String>,
    #[arg(long)]
    pub(crate) method: Option<String>,
    #[arg(long)]
    pub(crate) geneset: Option<String>,
    #[arg(long = "out")]
    pub(crate) out_spec: Option<String>,
    #[arg(long = "module-gmt")]
    pub(crate) module_gmt: Option<PathBuf>,
    #[arg(long = "module-top", default_value_t = 500)]
    pub(crate) module_top: usize,
    #[arg(long = "module-min-size", default_value_t = 25)]
    pub(crate) module_min_size: usize,
    #[arg(long = "fgsea-permutations", default_value_t = 1000)]
    pub(crate) fgsea_permutations: usize,
    #[arg(long = "fgsea-seed", default_value_t = 42)]
    pub(crate) fgsea_seed: u64,
    #[arg(long = "fgsea-min-size", default_value_t = 10)]
    pub(crate) fgsea_min_size: usize,
    #[arg(long = "fgsea-max-size")]
    pub(crate) fgsea_max_size: Option<usize>,
}

#[derive(Args, Debug)]
pub(crate) struct QcArgs {
    #[command(subcommand)]
    pub(crate) command: QcCommand,
}

#[derive(Subcommand, Debug)]
pub(crate) enum QcCommand {
    Reevaluate(QcReevaluateArgs),
}

#[derive(Args, Debug)]
pub(crate) struct QcReevaluateArgs {}

#[derive(Args, Debug)]
pub(crate) struct RunArgs {
    /// Exact-match cache allowance across mapping threads of each run, MiB.
    #[arg(long = "match-cache-mb")]
    pub(crate) match_cache_mb: Option<usize>,
    #[arg(long, value_parser = ["none", "fastq", "sra", "all"])]
    pub(crate) cleanup: Option<String>,
    #[arg(long, value_parser = ["success", "always", "never"])]
    pub(crate) cleanup_on: Option<String>,
    #[arg(long)]
    pub(crate) trash_days: Option<u64>,
    #[arg(long, default_value_t = false)]
    pub(crate) no_download: bool,
    #[arg(long, default_value_t = false)]
    pub(crate) no_corr: bool,
    #[arg(long = "worker-budget")]
    pub(crate) worker_budget: Option<usize>,
    #[arg(long = "storage-limit")]
    pub(crate) storage_limit: Option<String>,
    #[arg(long = "fastq-retention")]
    pub(crate) fastq_retention: Option<String>,
    #[arg(long = "preprocess-threads")]
    pub(crate) preprocess_threads: Option<usize>,
}

#[derive(Args, Debug)]
pub(crate) struct ExportArgs {
    #[arg(long)]
    pub(crate) format: Option<String>,
    #[arg(long)]
    pub(crate) what: Option<String>,
    /// Refuse results whose recorded cohort differs from the current active cohort.
    #[arg(long)]
    pub(crate) require_current_cohort: bool,
    /// Also check current cohort, QC eligibility, abundance inputs and selected reference.
    /// This does not check model settings or metadata values.
    #[arg(long)]
    pub(crate) require_current_admission: bool,
    /// Also compare recorded model designs with current used sample metadata.
    /// Other scientific settings and inputs are not fully revalidated.
    #[arg(long)]
    pub(crate) require_current_design: bool,
    /// Reconstruct current scientific requests; unprovable inputs fail closed.
    #[arg(long)]
    pub(crate) require_current: bool,
}

#[derive(Args, Debug)]
pub(crate) struct SurveyArgs {
    #[arg(long)]
    pub(crate) stratify: Option<String>,
    /// Refuse results whose recorded cohort differs from the current active cohort.
    #[arg(long)]
    pub(crate) require_current_cohort: bool,
    /// Also check current cohort, QC eligibility, abundance inputs and selected reference.
    /// This does not check model settings or metadata values.
    #[arg(long)]
    pub(crate) require_current_admission: bool,
    /// Also compare recorded model designs with current used sample metadata.
    /// Other scientific settings and inputs are not fully revalidated.
    #[arg(long)]
    pub(crate) require_current_design: bool,
    /// Reconstruct current scientific requests; unprovable inputs fail closed.
    #[arg(long)]
    pub(crate) require_current: bool,
}

pub(crate) fn parse_merge_strategy(value: &str) -> Result<MetadataMergeStrategy> {
    match value {
        "override" => Ok(MetadataMergeStrategy::Override),
        "prefer_remote" => Ok(MetadataMergeStrategy::PreferRemote),
        "prefer_local" => Ok(MetadataMergeStrategy::PreferLocal),
        "union-with-prefix" | "union_with_prefix" => Ok(MetadataMergeStrategy::UnionWithPrefix),
        _ => bail!("unsupported merge strategy: {value}"),
    }
}

pub(crate) fn parse_download_preference(value: &str) -> Result<DownloadPreference> {
    match value {
        "fastq" | "ena_fastq" => Ok(DownloadPreference::Fastq),
        "sra" => Ok(DownloadPreference::Sra),
        _ => bail!("unsupported download preference: {value}"),
    }
}

pub(crate) fn parse_cleanup_mode(value: &str) -> Result<CleanupMode> {
    match value {
        "none" => Ok(CleanupMode::None),
        "fastq" => Ok(CleanupMode::Fastq),
        "sra" => Ok(CleanupMode::Sra),
        "all" => Ok(CleanupMode::All),
        _ => bail!("unsupported cleanup mode: {value}"),
    }
}

pub(crate) fn parse_cleanup_when(value: &str) -> Result<CleanupWhen> {
    match value {
        "success" => Ok(CleanupWhen::Success),
        "always" => Ok(CleanupWhen::Always),
        "never" => Ok(CleanupWhen::Never),
        _ => bail!("unsupported cleanup schedule: {value}"),
    }
}

pub(crate) fn parse_fastq_retention(value: &str) -> Result<FastqRetention> {
    match value {
        "both" => Ok(FastqRetention::Both),
        "original" => Ok(FastqRetention::Original),
        "trimmed" => Ok(FastqRetention::Trimmed),
        "none" => Ok(FastqRetention::None),
        _ => bail!("unsupported fastq retention mode: {value}"),
    }
}

pub(crate) fn fastq_retention_label(value: FastqRetention) -> &'static str {
    match value {
        FastqRetention::Both => "both",
        FastqRetention::Original => "original",
        FastqRetention::Trimmed => "trimmed",
        FastqRetention::None => "none",
    }
}

pub(crate) fn validate_where_predicate(value: &str) -> Result<()> {
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

pub(crate) fn parse_cli_contrasts(items: &[String]) -> Result<Vec<ContrastSpec>> {
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
