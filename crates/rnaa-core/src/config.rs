use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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
    pub resources: crate::resources::ResourceConfig,
    #[serde(default)]
    pub quant: QuantConfig,
    #[serde(default)]
    pub metadata: crate::metadata::MetadataConfig,
    #[serde(default)]
    pub deseq2: Deseq2Config,
    #[serde(default)]
    pub corr: CorrConfig,
}

impl ProjectConfig {
    pub fn validate(&self) -> Result<()> {
        if self.project.schema_version != 1 {
            anyhow::bail!(
                "unsupported project configuration schema_version {}; expected 1",
                self.project.schema_version
            );
        }
        self.quant.validate()?;
        self.resources.validate()?;
        self.corr.validate()?;
        self.metadata.validate()?;
        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self> {
        let content = fs::read_to_string(path)
            .with_context(|| format!("failed to read config at {}", path.display()))?;
        let config: Self = toml::from_str(&content)
            .with_context(|| format!("failed to parse TOML config {}", path.display()))?;
        config.validate()?;
        Ok(config)
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        self.validate()?;
        let content = toml::to_string_pretty(self).context("failed to serialize config")?;
        fs::write(path, content)
            .with_context(|| format!("failed to write config at {}", path.display()))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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
#[serde(deny_unknown_fields)]
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
    #[serde(default = "default_http_timeout")]
    pub connect_timeout_seconds: u64,
    #[serde(default = "default_http_timeout")]
    pub read_timeout_seconds: u64,
    #[serde(default = "default_transfer_timeout")]
    pub max_transfer_seconds: u64,
}

impl Default for DownloadConfig {
    fn default() -> Self {
        Self {
            concurrency: default_download_concurrency(),
            retries: default_download_retries(),
            backoff_seconds: default_backoff_seconds(),
            prefer: DownloadPreference::Fastq,
            method: DownloadMethod::Reqwest,
            connect_timeout_seconds: default_http_timeout(),
            read_timeout_seconds: default_http_timeout(),
            max_transfer_seconds: default_transfer_timeout(),
        }
    }
}

fn default_http_timeout() -> u64 {
    30
}
fn default_transfer_timeout() -> u64 {
    86_400
}

impl DownloadConfig {
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.concurrency == 0
            || self.connect_timeout_seconds == 0
            || self.read_timeout_seconds == 0
            || self.max_transfer_seconds == 0
        {
            anyhow::bail!(
                "download concurrency and connect/read/transfer timeouts must be positive"
            );
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RefsConfig {
    #[serde(default = "default_organism")]
    pub organism: String,
    #[serde(default = "default_ensembl")]
    pub ensembl: String,
    #[serde(default)]
    pub custom_gtf: String,
    #[serde(default)]
    pub custom_cdna: String,
    #[serde(default = "default_index_k")]
    pub index_k: usize,
    #[serde(default = "default_index_threads")]
    pub index_threads: usize,
    #[serde(default = "default_reference_line_mb")]
    pub max_source_line_mb: u64,
    #[serde(default)]
    pub max_decoded_mb: u64,
}

fn default_reference_line_mb() -> u64 {
    1
}

fn default_index_k() -> usize {
    31
}
fn default_index_threads() -> usize {
    1
}
impl RefsConfig {
    pub fn validate(&self) -> anyhow::Result<()> {
        self.max_source_line_mb
            .checked_mul(1024 * 1024)
            .filter(|n| *n > 0 && usize::try_from(*n).is_ok())
            .ok_or_else(|| {
                anyhow::anyhow!("refs.max_source_line_mb must be positive and addressable")
            })?;
        self.max_decoded_mb
            .checked_mul(1024 * 1024)
            .ok_or_else(|| anyhow::anyhow!("refs.max_decoded_mb overflows byte accounting"))?;
        if !(5..=31).contains(&self.index_k)
            || self.index_k.is_multiple_of(2)
            || self.index_threads == 0
        {
            anyhow::bail!(
                "reference index_k must be odd and between 5 and 31; index_threads must be positive"
            );
        }
        Ok(())
    }
}

impl Default for RefsConfig {
    fn default() -> Self {
        Self {
            index_k: default_index_k(),
            index_threads: default_index_threads(),
            max_source_line_mb: default_reference_line_mb(),
            max_decoded_mb: 0,
            organism: default_organism(),
            ensembl: default_ensembl(),
            custom_gtf: String::new(),
            custom_cdna: String::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct StorageConfig {
    #[serde(default)]
    pub shared_root: String,
    #[serde(default)]
    pub max_active_size: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct QuantConfig {
    #[serde(default = "default_quant_engine")]
    pub engine: String,
    #[serde(default)]
    pub strand: QuantStrand,
    #[serde(default = "default_single_fragment_length")]
    pub single_fragment_length: f64,
    #[serde(default = "default_single_fragment_sd")]
    pub single_fragment_sd: f64,
    #[serde(default = "default_worker_budget")]
    pub worker_budget: usize,
    #[serde(default = "default_quant_threads")]
    pub threads: usize,
    /// Aggregate exact-match lookup cache across mapping threads for one run, MiB.
    #[serde(default = "default_quant_match_cache_mb")]
    pub match_cache_mb: usize,
    #[serde(default = "default_fastq_max_line_mb")]
    pub fastq_max_line_mb: usize,
    #[serde(default)]
    pub fastq_max_decoded_mb: usize,
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
            strand: QuantStrand::Unstranded,
            single_fragment_length: default_single_fragment_length(),
            single_fragment_sd: default_single_fragment_sd(),
            worker_budget: default_worker_budget(),
            threads: default_quant_threads(),
            match_cache_mb: default_quant_match_cache_mb(),
            fastq_max_line_mb: default_fastq_max_line_mb(),
            fastq_max_decoded_mb: 0,
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
#[serde(deny_unknown_fields)]
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
#[serde(deny_unknown_fields)]
pub struct CorrConfig {
    #[serde(default = "default_corr_block_rows")]
    pub block_rows: usize,
    #[serde(default = "default_corr_memory_mb")]
    pub max_memory_mb: usize,
    #[serde(default = "default_corr_output_mb")]
    pub max_output_mb: usize,
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
    /// External matrix selected by the latest correlation invocation, or the current DE matrix.
    #[serde(default)]
    pub matrix: Option<PathBuf>,
    /// Options from the most recent correlation invocation; absent when modules were not requested.
    #[serde(default)]
    pub modules: Option<CorrModuleConfig>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CorrModuleConfig {
    pub gmt_path: PathBuf,
    pub top_genes: usize,
    pub min_module_size: usize,
    pub permutations: usize,
    pub seed: u64,
    pub pathway_min_size: usize,
    pub pathway_max_size: Option<usize>,
}

impl Default for CorrConfig {
    fn default() -> Self {
        Self {
            block_rows: default_corr_block_rows(),
            max_memory_mb: default_corr_memory_mb(),
            max_output_mb: default_corr_output_mb(),
            adjust: default_adjust(),
            model: default_design(),
            method: default_corr_method(),
            geneset: default_geneset(),
            output: default_output(),
            threads: default_corr_threads(),
            matrix: None,
            modules: None,
        }
    }
}

fn default_corr_block_rows() -> usize {
    64
}
fn default_corr_memory_mb() -> usize {
    1024
}
fn default_corr_output_mb() -> usize {
    65536
}
impl CorrConfig {
    fn validate(&self) -> Result<()> {
        if !(1..=1024).contains(&self.block_rows)
            || self.max_memory_mb == 0
            || self.max_output_mb == 0
        {
            anyhow::bail!(
                "corr.block_rows must be in 1..=1024 and correlation memory/output limits must be positive"
            );
        }
        Ok(())
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

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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
    "kallistors".to_string()
}

fn default_worker_budget() -> usize {
    0
}

fn default_quant_threads() -> usize {
    16
}

fn default_quant_match_cache_mb() -> usize {
    64
}
fn default_fastq_max_line_mb() -> usize {
    1
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

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QuantStrand {
    #[default]
    Unstranded,
    Forward,
    Reverse,
}

impl std::str::FromStr for QuantStrand {
    type Err = anyhow::Error;
    fn from_str(value: &str) -> Result<Self> {
        match value {
            "unstranded" => Ok(Self::Unstranded),
            "forward" => Ok(Self::Forward),
            "reverse" => Ok(Self::Reverse),
            _ => anyhow::bail!("strand must be unstranded, forward, or reverse"),
        }
    }
}
fn default_single_fragment_length() -> f64 {
    200.0
}
fn default_single_fragment_sd() -> f64 {
    20.0
}

impl QuantConfig {
    pub fn fastq_max_line_bytes(&self) -> Result<usize> {
        checked_input_mb(self.fastq_max_line_mb, "quant.fastq_max_line_mb", false)
    }
    pub fn fastq_max_decoded_bytes(&self) -> Result<Option<u64>> {
        optional_input_limit(self.fastq_max_decoded_mb, "quant.fastq_max_decoded_mb")
    }
    pub fn preprocess_max_decoded_bytes(&self) -> Result<Option<u64>> {
        let a = self.fastq_max_decoded_bytes()?;
        let b = optional_input_limit(
            self.preprocess_max_input_mb,
            "quant.preprocess_max_input_mb",
        )?;
        Ok(match (a, b) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (a, b) => a.or(b),
        })
    }
    pub fn match_cache_bytes(&self) -> Result<u64> {
        let bytes = self
            .match_cache_mb
            .checked_mul(1024 * 1024)
            .filter(|bytes| *bytes > 0 && *bytes <= isize::MAX as usize)
            .context(
                "quant.match_cache_mb must be positive and fit an addressable byte allowance",
            )?;
        Ok(u64::try_from(bytes)?)
    }

    pub fn validate(&self) -> Result<()> {
        self.match_cache_bytes()?;
        self.fastq_max_line_bytes()?;
        self.preprocess_max_decoded_bytes()?;
        if self.engine != "kallistors" {
            anyhow::bail!(
                "unsupported quant engine '{}'; expected kallistors",
                self.engine
            );
        }
        if self.threads == 0 || self.preprocess_threads == 0 {
            anyhow::bail!("quant and preprocess threads must be positive");
        }
        for (name, value) in [
            ("single_fragment_length", self.single_fragment_length),
            ("single_fragment_sd", self.single_fragment_sd),
        ] {
            if !value.is_finite() || value <= 0.0 {
                anyhow::bail!("{name} must be finite and positive");
            }
        }
        for threshold in [
            self.qc_gate.min_pass_rate,
            self.qc_gate.max_low_quality_rate,
            self.qc_gate.max_low_complexity_rate,
            self.qc_gate.max_too_many_n_rate,
            self.qc_gate.max_too_short_rate,
            self.qc_gate.max_duplication_rate,
        ]
        .into_iter()
        .flatten()
        {
            if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
                anyhow::bail!("QC thresholds must be finite rates between 0 and 1");
            }
        }
        Ok(())
    }
}

fn checked_input_mb(value: usize, name: &str, allow_zero: bool) -> Result<usize> {
    value
        .checked_mul(1024 * 1024)
        .filter(|bytes| (*bytes > 0 || allow_zero) && *bytes < isize::MAX as usize)
        .with_context(|| {
            format!(
                "{name} must fit an addressable byte limit{}",
                if allow_zero { "" } else { " and be positive" }
            )
        })
}
fn optional_input_limit(value: usize, name: &str) -> Result<Option<u64>> {
    if value == 0 {
        Ok(None)
    } else {
        Ok(Some(u64::try_from(checked_input_mb(value, name, true)?)?))
    }
}

#[cfg(test)]
mod quant_validation_tests {
    use super::*;
    #[test]
    fn decoded_fastq_limits_are_checked_and_legacy_configs_get_bounded_lines() {
        let mut config = QuantConfig::default();
        assert_eq!(config.fastq_max_line_bytes().unwrap(), 1024 * 1024);
        assert_eq!(config.fastq_max_decoded_bytes().unwrap(), None);
        config.fastq_max_decoded_mb = 2;
        config.preprocess_max_input_mb = 1;
        assert_eq!(
            config.preprocess_max_decoded_bytes().unwrap(),
            Some(1024 * 1024)
        );
        config.preprocess_max_input_mb = 0;
        assert_eq!(
            config.preprocess_max_decoded_bytes().unwrap(),
            Some(2 * 1024 * 1024)
        );
        for (line, decoded, preprocess) in [
            (0, 0, 0),
            (usize::MAX, 0, 0),
            (1, usize::MAX, 0),
            (1, 0, usize::MAX),
        ] {
            let mut bad = config.clone();
            bad.fastq_max_line_mb = line;
            bad.fastq_max_decoded_mb = decoded;
            bad.preprocess_max_input_mb = preprocess;
            assert!(bad.validate().is_err());
        }
        let mut saved = serde_json::to_value(QuantConfig::default()).unwrap();
        saved.as_object_mut().unwrap().remove("fastq_max_line_mb");
        saved
            .as_object_mut()
            .unwrap()
            .remove("fastq_max_decoded_mb");
        let loaded: QuantConfig = serde_json::from_value(saved).unwrap();
        assert_eq!(loaded.fastq_max_line_mb, 1);
        assert_eq!(loaded.fastq_max_decoded_mb, 0);
    }
    #[test]
    fn mapping_cache_defaults_and_invalid_allowances_are_checked_before_allocation() {
        let config = QuantConfig::default();
        assert_eq!(config.match_cache_bytes().unwrap(), 64 * 1024 * 1024);
        for mb in [0, usize::MAX] {
            let invalid = QuantConfig {
                match_cache_mb: mb,
                ..config.clone()
            };
            assert!(invalid.validate().is_err());
        }
        let mut legacy = serde_json::to_value(&config).unwrap();
        legacy.as_object_mut().unwrap().remove("match_cache_mb");
        let legacy: QuantConfig = serde_json::from_value(legacy).unwrap();
        assert_eq!(legacy.match_cache_mb, 64);
    }
    #[test]
    fn scientific_parameters_must_be_positive_and_finite() {
        let mut config = QuantConfig::default();
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            config.single_fragment_length = value;
            assert!(config.validate().is_err());
        }
        config.single_fragment_length = 250.0;
        config.single_fragment_sd = 30.0;
        assert!(config.validate().is_ok());
        config.qc_gate.min_pass_rate = Some(1.5);
        assert!(config.validate().is_err());
    }
    #[test]
    fn shipped_configuration_parses_and_typoed_sections_or_fields_fail() {
        let example = include_str!("../../../configs/rnaa.example.toml");
        let config: ProjectConfig = toml::from_str(example).unwrap();
        config.validate().unwrap();
        for bad in [
            example.replace("[metadata]", "[metdata]"),
            example.replace("minimum_replicates", "minumum_replicates"),
            example.replace("[deseq2]\ndesign", "[deseq2]\ndesgin"),
            example.replace("read_timeout_seconds", "read_timout_seconds"),
        ] {
            assert!(
                toml::from_str::<ProjectConfig>(&bad).is_err(),
                "configuration typo must not be silently ignored"
            );
        }
        let mut future = config;
        future.project.schema_version = 2;
        assert!(future.validate().is_err());
    }
}
