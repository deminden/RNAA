use anyhow::Result;

use crate::config::ProjectConfig;
use crate::model::{
    AdjustedMatrix, ContrastSpec, CorrelationMethod, OutputMode, QuantArtifacts, ReferenceBundle,
    ResolvedProject, RunRecord, VerifiedFile,
};
use crate::paths::ProjectPaths;

pub trait Resolver: Send + Sync {
    fn resolve(&self, inputs: &[crate::model::InputRecord]) -> Result<ResolvedProject>;
}

pub trait Downloader: Send + Sync {
    fn ensure_downloaded(
        &self,
        run: &RunRecord,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<Vec<VerifiedFile>>;
}

pub trait Extractor: Send + Sync {
    fn ensure_fastq(&self, run: &RunRecord, paths: &ProjectPaths) -> Result<Vec<VerifiedFile>>;
}

pub trait ReferenceManager: Send + Sync {
    fn ensure_reference(
        &self,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<ReferenceBundle>;
}

pub trait Quantifier: Send + Sync {
    fn quantify(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        reference: &ReferenceBundle,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<QuantArtifacts>;
}

pub trait DifferentialExpression: Send + Sync {
    fn deseq2(
        &self,
        project_id: &str,
        reference: &ReferenceBundle,
        design: &str,
        contrasts: &[ContrastSpec],
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<Vec<std::path::PathBuf>>;
}

pub trait MatrixAdjuster: Send + Sync {
    fn adjust(
        &self,
        matrix_path: &std::path::Path,
        samplesheet_path: &std::path::Path,
        model: &str,
        paths: &ProjectPaths,
    ) -> Result<AdjustedMatrix>;
}

pub trait Correlator: Send + Sync {
    fn correlate(
        &self,
        matrix: &AdjustedMatrix,
        method: CorrelationMethod,
        output_mode: &OutputMode,
        paths: &ProjectPaths,
        project_id: &str,
    ) -> Result<Vec<std::path::PathBuf>>;
}

pub trait Cleaner: Send + Sync {
    fn cleanup(&self, run: &RunRecord, paths: &ProjectPaths, config: &ProjectConfig) -> Result<()>;
}
