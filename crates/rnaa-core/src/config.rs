use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ProjectConfig {
    pub project: ProjectSection,
    #[serde(default)]
    pub resolver: ResolverConfig,
    #[serde(default)]
    pub download: DownloadConfig,
    #[serde(default)]
    pub refs: RefsConfig,
    #[serde(default)]
    pub storage: StorageConfig,
    #[serde(default)]
    pub quant: QuantConfig,
    #[serde(default)]
    pub deseq2: Deseq2Config,
    #[serde(default)]
    pub corr: CorrConfig,
}

impl ProjectConfig {
    pub fn load(path: &Path) -> Result<Self> {
        let content = fs::read_to_string(path)
            .with_context(|| format!("failed to read config at {}", path.display()))?;
        toml::from_str(&content)
            .with_context(|| format!("failed to parse TOML config {}", path.display()))
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        let content = toml::to_string_pretty(self).context("failed to serialize config")?;
        fs::write(path, content)
            .with_context(|| format!("failed to write config at {}", path.display()))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectSection {
    pub name: String,
    pub schema_version: u32,
}

impl Default for ProjectSection {
    fn default() -> Self {
        Self {
            name: "rnaa-project".to_string(),
            schema_version: 1,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResolverConfig {
    #[serde(default = "default_resolver_provider")]
    pub provider: String,
    #[serde(default = "default_resolver_fields")]
    pub fields: Vec<String>,
    #[serde(default)]
    pub merge_strategy: MetadataMergeStrategy,
}

impl Default for ResolverConfig {
    fn default() -> Self {
        Self {
            provider: default_resolver_provider(),
            fields: default_resolver_fields(),
            merge_strategy: MetadataMergeStrategy::PreferLocal,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DownloadConfig {
    #[serde(default = "default_download_concurrency")]
    pub concurrency: usize,
    #[serde(default = "default_download_retries")]
    pub retries: u32,
    #[serde(default = "default_backoff_seconds")]
    pub backoff_seconds: Vec<u64>,
    #[serde(default)]
    pub prefer: DownloadPreference,
    #[serde(default)]
    pub method: DownloadMethod,
}

impl Default for DownloadConfig {
    fn default() -> Self {
        Self {
            concurrency: default_download_concurrency(),
            retries: default_download_retries(),
            backoff_seconds: default_backoff_seconds(),
            prefer: DownloadPreference::Fastq,
            method: DownloadMethod::Reqwest,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefsConfig {
    #[serde(default = "default_organism")]
    pub organism: String,
    #[serde(default = "default_ensembl")]
    pub ensembl: String,
    #[serde(default)]
    pub custom_gtf: String,
    #[serde(default)]
    pub custom_cdna: String,
}

impl Default for RefsConfig {
    fn default() -> Self {
        Self {
            organism: default_organism(),
            ensembl: default_ensembl(),
            custom_gtf: String::new(),
            custom_cdna: String::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct StorageConfig {
    #[serde(default)]
    pub shared_root: String,
    #[serde(default)]
    pub max_active_size: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuantConfig {
    #[serde(default = "default_quant_engine")]
    pub engine: String,
    #[serde(default = "default_worker_budget")]
    pub worker_budget: usize,
    #[serde(default = "default_quant_threads")]
    pub threads: usize,
    #[serde(default = "default_quant_workers")]
    pub workers: usize,
    #[serde(default)]
    pub preprocess: bool,
    #[serde(default)]
    pub fastq_retention: FastqRetention,
    #[serde(default = "default_quant_preprocess_threads")]
    pub preprocess_threads: usize,
    #[serde(default)]
    pub qc_gate: PreprocessQcGate,
    #[serde(default = "default_quant_preprocess_strict")]
    pub preprocess_strict: bool,
    #[serde(default = "default_quant_preprocess_max_input_mb")]
    pub preprocess_max_input_mb: usize,
    #[serde(default)]
    pub cleanup: CleanupMode,
    #[serde(default = "default_trash_days")]
    pub trash_days: u64,
    #[serde(default)]
    pub cleanup_on: CleanupWhen,
}

impl Default for QuantConfig {
    fn default() -> Self {
        Self {
            engine: default_quant_engine(),
            worker_budget: default_worker_budget(),
            threads: default_quant_threads(),
            workers: default_quant_workers(),
            preprocess: false,
            fastq_retention: FastqRetention::Both,
            preprocess_threads: default_quant_preprocess_threads(),
            qc_gate: PreprocessQcGate::default(),
            preprocess_strict: default_quant_preprocess_strict(),
            preprocess_max_input_mb: default_quant_preprocess_max_input_mb(),
            cleanup: CleanupMode::None,
            trash_days: default_trash_days(),
            cleanup_on: CleanupWhen::Success,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Deseq2Config {
    #[serde(default = "default_design")]
    pub design: String,
    #[serde(default = "default_counts_from_abundance")]
    pub counts_from_abundance: String,
    #[serde(default = "default_transform")]
    pub transform: String,
}

impl Default for Deseq2Config {
    fn default() -> Self {
        Self {
            design: default_design(),
            counts_from_abundance: default_counts_from_abundance(),
            transform: default_transform(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorrConfig {
    #[serde(default = "default_adjust")]
    pub adjust: String,
    #[serde(default = "default_design")]
    pub model: String,
    #[serde(default = "default_corr_method")]
    pub method: String,
    #[serde(default = "default_geneset")]
    pub geneset: String,
    #[serde(default = "default_output")]
    pub output: String,
    #[serde(default = "default_corr_threads")]
    pub threads: usize,
}

impl Default for CorrConfig {
    fn default() -> Self {
        Self {
            adjust: default_adjust(),
            model: default_design(),
            method: default_corr_method(),
            geneset: default_geneset(),
            output: default_output(),
            threads: default_corr_threads(),
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum MetadataMergeStrategy {
    Override,
    PreferRemote,
    #[default]
    PreferLocal,
    #[serde(rename = "union-with-prefix")]
    UnionWithPrefix,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum DownloadPreference {
    #[default]
    Fastq,
    Sra,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum DownloadMethod {
    #[default]
    #[serde(alias = "curl", alias = "wget")]
    Reqwest,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum CleanupMode {
    #[default]
    None,
    Fastq,
    Sra,
    All,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum CleanupWhen {
    #[default]
    Success,
    Always,
    Never,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum FastqRetention {
    #[default]
    Both,
    Original,
    Trimmed,
    None,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PreprocessQcGate {
    #[serde(default)]
    pub min_pass_rate: Option<f64>,
    #[serde(default)]
    pub max_low_quality_rate: Option<f64>,
    #[serde(default)]
    pub max_low_complexity_rate: Option<f64>,
    #[serde(default)]
    pub max_too_many_n_rate: Option<f64>,
    #[serde(default)]
    pub max_too_short_rate: Option<f64>,
    #[serde(default)]
    pub max_duplication_rate: Option<f64>,
}

impl Default for PreprocessQcGate {
    fn default() -> Self {
        Self {
            // Conservative fail-fast defaults intended to block catastrophically bad inputs
            // without excluding runs for biologically variable duplication or complexity alone.
            min_pass_rate: Some(0.60),
            max_low_quality_rate: Some(0.30),
            max_low_complexity_rate: None,
            max_too_many_n_rate: Some(0.10),
            max_too_short_rate: Some(0.40),
            max_duplication_rate: None,
        }
    }
}

fn default_resolver_provider() -> String {
    "ena".to_string()
}

fn default_resolver_fields() -> Vec<String> {
    [
        "run_accession",
        "study_accession",
        "secondary_study_accession",
        "sample_accession",
        "secondary_sample_accession",
        "experiment_accession",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
        "submitted_ftp",
        "submitted_md5",
        "submitted_bytes",
    ]
    .into_iter()
    .map(ToString::to_string)
    .collect()
}

fn default_download_concurrency() -> usize {
    2
}

fn default_download_retries() -> u32 {
    6
}

fn default_backoff_seconds() -> Vec<u64> {
    vec![1, 2, 4, 8, 16, 32]
}

fn default_organism() -> String {
    "human".to_string()
}

fn default_ensembl() -> String {
    "latest".to_string()
}

fn default_quant_engine() -> String {
    "r-kallisto".to_string()
}

fn default_worker_budget() -> usize {
    0
}

fn default_quant_threads() -> usize {
    16
}

fn default_quant_workers() -> usize {
    0
}

fn default_quant_preprocess_strict() -> bool {
    true
}

fn default_quant_preprocess_threads() -> usize {
    16
}

fn default_quant_preprocess_max_input_mb() -> usize {
    1024
}

fn default_trash_days() -> u64 {
    7
}

fn default_design() -> String {
    "~ batch + condition".to_string()
}

fn default_counts_from_abundance() -> String {
    "lengthScaledTPM".to_string()
}

fn default_transform() -> String {
    "vst".to_string()
}

fn default_adjust() -> String {
    "residualize".to_string()
}

fn default_corr_method() -> String {
    "pearson".to_string()
}

fn default_geneset() -> String {
    "topvar:5000".to_string()
}

fn default_output() -> String {
    "edges:topk=50".to_string()
}

fn default_corr_threads() -> usize {
    0
}
