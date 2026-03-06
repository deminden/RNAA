use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};
use fasterp::wasm::{
    WasmConfig, process_paired_end, process_single_end, version as fasterp_version,
};
use rnaa_core::config::ProjectConfig;
use rnaa_core::model::{ArtifactKind, PreprocessArtifacts, RunRecord, VerifiedFile};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::Preprocessor;
use rnaa_core::util::{compute_sha256, file_size};

#[derive(Debug, Clone, Default)]
pub struct FasterpInProcessPreprocessor;

impl Preprocessor for FasterpInProcessPreprocessor {
    fn preprocess(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<PreprocessArtifacts> {
        let read1 = fastqs
            .iter()
            .find(|item| {
                item.artifact_kind == ArtifactKind::FastqR1
                    || item.artifact_kind == ArtifactKind::FastqSingle
            })
            .ok_or_else(|| {
                anyhow!(
                    "no read1 FASTQ artifact found for run {}",
                    run.run_accession
                )
            })?;
        let read2 = fastqs
            .iter()
            .find(|item| item.artifact_kind == ArtifactKind::FastqR2);

        enforce_input_size_limit(read1, config)?;
        if let Some(read2) = read2 {
            enforce_input_size_limit(read2, config)?;
        }

        let out_dir = paths
            .quant_run_dir(&run.run_accession, &config.quant.engine)
            .join("fasterp");
        fs::create_dir_all(&out_dir)
            .with_context(|| format!("failed to create {}", out_dir.display()))?;
        let report_json = out_dir.join("fasterp.json");
        let out_read1 = if read2.is_some() {
            out_dir.join("trimmed.R1.fastq")
        } else {
            out_dir.join("trimmed.fastq")
        };
        let out_read2 = read2.map(|_| out_dir.join("trimmed.R2.fastq"));
        let tool_version = fasterp_version();

        if is_cached_output_ready(&out_read1, out_read2.as_deref(), &report_json) {
            let fastqs = build_preprocessed_fastqs(read1, &out_read1, out_read2.as_deref())?;
            return Ok(PreprocessArtifacts {
                run_accession: run.run_accession.clone(),
                out_dir,
                fastqs,
                report_json,
                tool_name: "fasterp".to_string(),
                tool_version,
                passed_reads: 0,
                failed_reads: 0,
                reused: true,
            });
        }

        let wasm_config = WasmConfig::new();
        let (passed_reads, failed_reads) = if let Some(read2) = read2 {
            let input1 = fs::read(&read1.path)
                .with_context(|| format!("failed to read {}", read1.path.display()))?;
            let input2 = fs::read(&read2.path)
                .with_context(|| format!("failed to read {}", read2.path.display()))?;
            let result = process_paired_end(&input1, &input2, &wasm_config).map_err(|err| {
                anyhow!(
                    "fasterp in-process paired preprocessing failed for {}: {err:?}",
                    run.run_accession
                )
            })?;

            fs::write(&out_read1, result.output1())
                .with_context(|| format!("failed to write {}", out_read1.display()))?;
            let out_read2_path = out_read2
                .as_ref()
                .ok_or_else(|| anyhow!("paired run missing read2 output path"))?;
            fs::write(out_read2_path, result.output2())
                .with_context(|| format!("failed to write {}", out_read2_path.display()))?;
            fs::write(&report_json, result.json_report())
                .with_context(|| format!("failed to write {}", report_json.display()))?;
            (result.passed_reads(), result.failed_reads())
        } else {
            let input1 = fs::read(&read1.path)
                .with_context(|| format!("failed to read {}", read1.path.display()))?;
            let result = process_single_end(&input1, &wasm_config).map_err(|err| {
                anyhow!(
                    "fasterp in-process single-end preprocessing failed for {}: {err:?}",
                    run.run_accession
                )
            })?;

            fs::write(&out_read1, result.output())
                .with_context(|| format!("failed to write {}", out_read1.display()))?;
            fs::write(&report_json, result.json_report())
                .with_context(|| format!("failed to write {}", report_json.display()))?;
            (result.passed_reads(), result.failed_reads())
        };

        if !is_cached_output_ready(&out_read1, out_read2.as_deref(), &report_json) {
            bail!(
                "fasterp output missing or empty for run {} after preprocessing",
                run.run_accession
            );
        }

        let fastqs = build_preprocessed_fastqs(read1, &out_read1, out_read2.as_deref())?;
        Ok(PreprocessArtifacts {
            run_accession: run.run_accession.clone(),
            out_dir,
            fastqs,
            report_json,
            tool_name: "fasterp".to_string(),
            tool_version,
            passed_reads,
            failed_reads,
            reused: false,
        })
    }
}

fn is_cached_output_ready(out_read1: &Path, out_read2: Option<&Path>, report_json: &Path) -> bool {
    if !is_nonempty(out_read1) || !is_nonempty(report_json) {
        return false;
    }
    out_read2.map(is_nonempty).unwrap_or(true)
}

fn is_nonempty(path: &Path) -> bool {
    path.exists() && file_size(path).unwrap_or_default() > 0
}

fn enforce_input_size_limit(input: &VerifiedFile, config: &ProjectConfig) -> Result<()> {
    let limit_mb = config.quant.preprocess_max_input_mb;
    if limit_mb == 0 {
        return Ok(());
    }
    let limit_bytes = (limit_mb as u64).saturating_mul(1024 * 1024);
    let bytes = file_size(&input.path)
        .with_context(|| format!("failed to stat {}", input.path.display()))?;
    if bytes > limit_bytes {
        bail!(
            "fasterp in-process input {} is {} bytes, exceeds preprocess_max_input_mb={} (set to 0 to disable limit)",
            input.path.display(),
            bytes,
            limit_mb
        );
    }
    Ok(())
}

fn build_preprocessed_fastqs(
    read1_input: &VerifiedFile,
    out_read1: &Path,
    out_read2: Option<&Path>,
) -> Result<Vec<VerifiedFile>> {
    let read1_kind = match read1_input.artifact_kind {
        ArtifactKind::FastqSingle => ArtifactKind::FastqTrimmedSingle,
        ArtifactKind::FastqR1 => ArtifactKind::FastqTrimmedR1,
        other => {
            bail!(
                "unsupported read1 kind for fasterp preprocessing: {}",
                other.as_str()
            );
        }
    };

    let mut outputs = vec![verified_output(read1_kind, out_read1)?];
    if let Some(read2_path) = out_read2 {
        outputs.push(verified_output(ArtifactKind::FastqTrimmedR2, read2_path)?);
    }
    Ok(outputs)
}

fn verified_output(kind: ArtifactKind, path: &Path) -> Result<VerifiedFile> {
    let bytes = file_size(path)?;
    if bytes == 0 {
        bail!("expected non-empty fasterp output {}", path.display());
    }
    Ok(VerifiedFile {
        artifact_kind: kind,
        path: PathBuf::from(path),
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(path)?,
        bytes,
        integrity_source: "fasterp_inprocess".to_string(),
    })
}
