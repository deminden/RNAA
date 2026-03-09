use std::collections::BTreeMap;
use std::fmt::{Display, Formatter};
use std::path::PathBuf;
use std::str::FromStr;

use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};

use crate::config::MetadataMergeStrategy;
use crate::state::{ProjectStageState, RunState};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectRecord {
    pub project_id: String,
    pub root_dir: String,
    pub created_at: String,
    pub config_json: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum InputType {
    Run,
    Study,
    Project,
    Sample,
    SampleInclude,
    SampleExclude,
    FilterInclude,
    FilterExclude,
    Local,
}

impl InputType {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Run => "RUN",
            Self::Study => "STUDY",
            Self::Project => "PROJECT",
            Self::Sample => "SAMPLE",
            Self::SampleInclude => "SAMPLE_INCLUDE",
            Self::SampleExclude => "SAMPLE_EXCLUDE",
            Self::FilterInclude => "FILTER_INCLUDE",
            Self::FilterExclude => "FILTER_EXCLUDE",
            Self::Local => "LOCAL",
        }
    }

    pub fn from_accession(accession: &str) -> Result<Self> {
        let upper = accession.to_ascii_uppercase();
        if upper.starts_with("SRR") || upper.starts_with("ERR") || upper.starts_with("DRR") {
            return Ok(Self::Run);
        }
        if upper.starts_with("SRP") || upper.starts_with("ERP") || upper.starts_with("DRP") {
            return Ok(Self::Study);
        }
        if upper.starts_with("PRJ") {
            return Ok(Self::Project);
        }
        if upper.starts_with("SRS")
            || upper.starts_with("ERS")
            || upper.starts_with("DRS")
            || upper.starts_with("SAM")
        {
            return Ok(Self::Sample);
        }
        bail!("unsupported accession type: {accession}")
    }
}

impl Display for InputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for InputType {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "RUN" => Ok(Self::Run),
            "STUDY" => Ok(Self::Study),
            "PROJECT" => Ok(Self::Project),
            "SAMPLE" => Ok(Self::Sample),
            "SAMPLE_INCLUDE" => Ok(Self::SampleInclude),
            "SAMPLE_EXCLUDE" => Ok(Self::SampleExclude),
            "FILTER_INCLUDE" => Ok(Self::FilterInclude),
            "FILTER_EXCLUDE" => Ok(Self::FilterExclude),
            "LOCAL" => Ok(Self::Local),
            _ => bail!("unknown input type: {value}"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputRecord {
    pub project_id: String,
    pub input_id: String,
    pub input_type: InputType,
    pub added_at: String,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum LibraryLayout {
    Single,
    Paired,
    #[default]
    Unknown,
}

impl LibraryLayout {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Single => "SINGLE",
            Self::Paired => "PAIRED",
            Self::Unknown => "UNKNOWN",
        }
    }
}

impl Display for LibraryLayout {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for LibraryLayout {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value.to_ascii_uppercase().as_str() {
            "SINGLE" => Ok(Self::Single),
            "PAIRED" => Ok(Self::Paired),
            "" | "UNKNOWN" => Ok(Self::Unknown),
            _ => bail!("unknown library layout: {value}"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RemoteFileKind {
    Fastq,
    FastqR1,
    FastqR2,
    Sra,
    Submitted,
}

impl RemoteFileKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fastq => "FASTQ",
            Self::FastqR1 => "FASTQ_R1",
            Self::FastqR2 => "FASTQ_R2",
            Self::Sra => "SRA",
            Self::Submitted => "SUBMITTED",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RemoteFile {
    pub url: String,
    pub md5: Option<String>,
    pub bytes: Option<u64>,
    pub kind: RemoteFileKind,
    pub basename: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunRecord {
    pub project_id: String,
    pub study_accession: String,
    pub source_study_accession: Option<String>,
    pub run_accession: String,
    pub sample_accession: Option<String>,
    pub experiment_accession: Option<String>,
    pub library_layout: LibraryLayout,
    pub instrument_platform: Option<String>,
    pub instrument_model: Option<String>,
    pub remote_files: Vec<RemoteFile>,
    pub metadata: BTreeMap<String, String>,
    pub state: RunState,
    pub last_error: Option<String>,
    pub updated_at: String,
}

impl RunRecord {
    pub fn read1_remote(&self) -> Option<&RemoteFile> {
        self.remote_files
            .iter()
            .find(|file| file.kind == RemoteFileKind::FastqR1)
            .or_else(|| {
                self.remote_files.iter().find(|file| {
                    self.library_layout == LibraryLayout::Single
                        && file.kind == RemoteFileKind::Fastq
                })
            })
    }

    pub fn read2_remote(&self) -> Option<&RemoteFile> {
        self.remote_files
            .iter()
            .find(|file| file.kind == RemoteFileKind::FastqR2)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArtifactKind {
    FastqR1,
    FastqR2,
    FastqSingle,
    FastqTrimmedR1,
    FastqTrimmedR2,
    FastqTrimmedSingle,
    PreprocessReport,
    Sra,
    QuantDir,
    QuantAbundanceH5,
    QuantRunInfo,
    Counts,
    NormalizedCounts,
    Vst,
    VstRds,
    DeTable,
    CorrEdges,
    CorrDense,
    CorrModules,
    EnrichmentTable,
    AdjustedMatrix,
    ReferenceCdna,
    ReferenceGtf,
    ReferenceIndex,
    ReferenceTx2Gene,
    ReferenceGeneAnnotation,
    Samplesheet,
    Manifest,
    Trash,
}

impl ArtifactKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FastqR1 => "FASTQ_R1",
            Self::FastqR2 => "FASTQ_R2",
            Self::FastqSingle => "FASTQ_SINGLE",
            Self::FastqTrimmedR1 => "FASTQ_TRIMMED_R1",
            Self::FastqTrimmedR2 => "FASTQ_TRIMMED_R2",
            Self::FastqTrimmedSingle => "FASTQ_TRIMMED_SINGLE",
            Self::PreprocessReport => "PREPROCESS_REPORT",
            Self::Sra => "SRA",
            Self::QuantDir => "KALLISTO_DIR",
            Self::QuantAbundanceH5 => "ABUNDANCE_H5",
            Self::QuantRunInfo => "RUN_INFO_JSON",
            Self::Counts => "COUNTS",
            Self::NormalizedCounts => "NORM_COUNTS",
            Self::Vst => "VST",
            Self::VstRds => "VST_RDS",
            Self::DeTable => "DE_TABLE",
            Self::CorrEdges => "CORR_EDGES",
            Self::CorrDense => "CORR_DENSE",
            Self::CorrModules => "CORR_MODULES",
            Self::EnrichmentTable => "ENRICHMENT_TABLE",
            Self::AdjustedMatrix => "ADJUSTED_MATRIX",
            Self::ReferenceCdna => "REFERENCE_CDNA",
            Self::ReferenceGtf => "REFERENCE_GTF",
            Self::ReferenceIndex => "REFERENCE_INDEX",
            Self::ReferenceTx2Gene => "REFERENCE_TX2GENE",
            Self::ReferenceGeneAnnotation => "REFERENCE_GENE_ANNOTATION",
            Self::Samplesheet => "SAMPLESHEET",
            Self::Manifest => "MANIFEST",
            Self::Trash => "TRASH",
        }
    }
}

impl Display for ArtifactKind {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for ArtifactKind {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "FASTQ_R1" => Ok(Self::FastqR1),
            "FASTQ_R2" => Ok(Self::FastqR2),
            "FASTQ_SINGLE" => Ok(Self::FastqSingle),
            "FASTQ_TRIMMED_R1" => Ok(Self::FastqTrimmedR1),
            "FASTQ_TRIMMED_R2" => Ok(Self::FastqTrimmedR2),
            "FASTQ_TRIMMED_SINGLE" => Ok(Self::FastqTrimmedSingle),
            "PREPROCESS_REPORT" => Ok(Self::PreprocessReport),
            "SRA" => Ok(Self::Sra),
            "KALLISTO_DIR" => Ok(Self::QuantDir),
            "ABUNDANCE_H5" => Ok(Self::QuantAbundanceH5),
            "RUN_INFO_JSON" => Ok(Self::QuantRunInfo),
            "COUNTS" => Ok(Self::Counts),
            "NORM_COUNTS" => Ok(Self::NormalizedCounts),
            "VST" => Ok(Self::Vst),
            "VST_RDS" => Ok(Self::VstRds),
            "DE_TABLE" => Ok(Self::DeTable),
            "CORR_EDGES" => Ok(Self::CorrEdges),
            "CORR_DENSE" => Ok(Self::CorrDense),
            "CORR_MODULES" => Ok(Self::CorrModules),
            "ENRICHMENT_TABLE" => Ok(Self::EnrichmentTable),
            "ADJUSTED_MATRIX" => Ok(Self::AdjustedMatrix),
            "REFERENCE_CDNA" => Ok(Self::ReferenceCdna),
            "REFERENCE_GTF" => Ok(Self::ReferenceGtf),
            "REFERENCE_INDEX" => Ok(Self::ReferenceIndex),
            "REFERENCE_TX2GENE" => Ok(Self::ReferenceTx2Gene),
            "REFERENCE_GENE_ANNOTATION" => Ok(Self::ReferenceGeneAnnotation),
            "SAMPLESHEET" => Ok(Self::Samplesheet),
            "MANIFEST" => Ok(Self::Manifest),
            "TRASH" => Ok(Self::Trash),
            _ => bail!("unknown artifact kind: {value}"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArtifactRecord {
    pub project_id: String,
    pub run_accession: Option<String>,
    pub kind: ArtifactKind,
    pub path: String,
    #[serde(default)]
    pub blob_id: Option<String>,
    #[serde(default)]
    pub shared_path: Option<String>,
    pub checksum_type: String,
    pub checksum: String,
    pub bytes: u64,
    pub created_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SharedBlobRecord {
    pub blob_id: String,
    pub storage_path: String,
    pub checksum_type: String,
    pub checksum: String,
    pub bytes: u64,
    pub created_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetadataOverrideRecord {
    pub id: i64,
    pub project_id: String,
    pub source_path: String,
    pub stored_path: String,
    pub merge_strategy: MetadataMergeStrategy,
    pub added_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EventRecord {
    pub id: i64,
    pub project_id: String,
    pub ts: String,
    pub run_accession: Option<String>,
    pub stage: String,
    pub message: String,
    pub context: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectStageRecord {
    pub project_id: String,
    pub stage: String,
    pub state: ProjectStageState,
    pub last_error: Option<String>,
    pub updated_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContrastSpec {
    pub name: String,
    pub factor: String,
    pub level_a: String,
    pub level_b: String,
}

impl ContrastSpec {
    pub fn new(name: Option<String>, factor: String, level_a: String, level_b: String) -> Self {
        let auto_name = format!("{factor}_{level_a}_vs_{level_b}");
        Self {
            name: name.unwrap_or(auto_name),
            factor,
            level_a,
            level_b,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VerifiedFile {
    pub artifact_kind: ArtifactKind,
    pub path: PathBuf,
    pub checksum_type: String,
    pub checksum: String,
    pub bytes: u64,
    pub integrity_source: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceBundle {
    pub id: String,
    pub organism: String,
    pub ensembl_release: String,
    pub cdna_path: PathBuf,
    pub gtf_path: PathBuf,
    pub kallisto_index_path: PathBuf,
    pub tx2gene_path: PathBuf,
    pub gene_annotation_path: Option<PathBuf>,
    pub manifest_path: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuantArtifacts {
    pub run_accession: String,
    pub out_dir: PathBuf,
    pub abundance_h5: PathBuf,
    pub run_info_json: PathBuf,
    pub manifest_path: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PreprocessArtifacts {
    pub run_accession: String,
    pub out_dir: PathBuf,
    pub fastqs: Vec<VerifiedFile>,
    pub report_json: PathBuf,
    pub tool_name: String,
    pub tool_version: String,
    pub passed_reads: usize,
    pub failed_reads: usize,
    pub reused: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PreprocessQcRecord {
    pub project_id: String,
    pub run_accession: String,
    pub report_path: String,
    pub mode: String,
    pub threads: u64,
    pub total_reads_before: u64,
    pub passed_reads: u64,
    pub failed_reads: u64,
    pub pass_rate: f64,
    pub failed_rate: f64,
    pub low_quality_reads: u64,
    pub low_complexity_reads: u64,
    pub too_many_n_reads: u64,
    pub too_short_reads: u64,
    pub duplicated_reads: u64,
    pub duplication_rate: Option<f64>,
    pub adapter_trimmed_reads: u64,
    pub adapter_trimmed_bases: u64,
    pub gate_passed: bool,
    pub gate_reason: Option<String>,
    pub updated_at: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdjustedMatrix {
    pub matrix_path: PathBuf,
    pub genes: Vec<String>,
    pub samples: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum CorrelationMethod {
    Pearson,
    Spearman,
    Kendall,
    Bicor,
    Hellcor,
}

impl Display for CorrelationMethod {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let value = match self {
            Self::Pearson => "pearson",
            Self::Spearman => "spearman",
            Self::Kendall => "kendall",
            Self::Bicor => "bicor",
            Self::Hellcor => "hellcor",
        };
        f.write_str(value)
    }
}

impl FromStr for CorrelationMethod {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value.to_ascii_lowercase().as_str() {
            "pearson" => Ok(Self::Pearson),
            "spearman" => Ok(Self::Spearman),
            "kendall" => Ok(Self::Kendall),
            "bicor" => Ok(Self::Bicor),
            "hellcor" => Ok(Self::Hellcor),
            _ => bail!("unsupported correlation method: {value}"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum OutputMode {
    TopK { k: usize },
    Threshold { min_abs_r: f64 },
    Dense,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResolvedProject {
    pub runs: Vec<RunRecord>,
    pub metadata_columns: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SamplesheetRow {
    pub project_id: String,
    pub study_accession: String,
    pub run_accession: String,
    pub sample_accession: Option<String>,
    pub experiment_accession: Option<String>,
    pub library_layout: String,
    pub instrument_platform: Option<String>,
    pub instrument_model: Option<String>,
    pub read1_path: Option<String>,
    pub read2_path: Option<String>,
    pub condition: Option<String>,
    pub batch: Option<String>,
    pub metadata: BTreeMap<String, String>,
}

#[cfg(test)]
mod tests {
    use super::InputType;

    #[test]
    fn detects_sample_accessions() {
        assert_eq!(
            InputType::from_accession("SRS123").unwrap(),
            InputType::Sample
        );
        assert_eq!(
            InputType::from_accession("SAMN123").unwrap(),
            InputType::Sample
        );
    }
}
