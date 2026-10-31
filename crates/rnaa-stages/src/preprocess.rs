mod native;
mod resources;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};
use fasterp::{NativeProcessOptions, process_paths_with_options, version as fasterp_version};
use rnaa_core::config::ProjectConfig;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{ArtifactKind, LibraryLayout, PreprocessArtifacts, RunRecord, VerifiedFile};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::Preprocessor;
use rnaa_core::util::{compute_sha256, file_size, now_rfc3339, write_json_pretty};

#[derive(Debug, Clone, Default)]
pub struct FasterpInProcessPreprocessor;

pub fn request_digest(
    run: &RunRecord,
    inputs: &[VerifiedFile],
    config: &ProjectConfig,
) -> Result<String> {
    request_digest_from_inputs(run, &snapshot_inputs(inputs)?, config)
}

pub fn request_digest_from_inputs(
    run: &RunRecord,
    inputs: &[ManifestArtifact],
    config: &ProjectConfig,
) -> Result<String> {
    request_from_inputs(run, inputs, config)?.digest()
}

pub fn request_from_inputs(
    run: &RunRecord,
    inputs: &[ManifestArtifact],
    config: &ProjectConfig,
) -> Result<rnaa_core::fingerprint::StageRequest<serde_json::Value>> {
    use rnaa_core::fingerprint::{InputIdentity, StageRequest};
    Ok(StageRequest {
        schema_version: 1,
        stage: "preprocess".into(),
        algorithm_version: 3,
        engines: BTreeMap::from([
            ("fasterp".into(), fasterp_version()),
            (
                "dependency_lock".into(),
                crate::engines::DEPENDENCY_LOCK_SHA256.into(),
            ),
        ]),
        parameters: serde_json::json!({"layout":run.library_layout,
            "strict":config.quant.preprocess_strict, "max_input_mb":config.quant.preprocess_max_input_mb,
            "qc_gate":config.quant.qc_gate}),
        inputs: inputs.iter().map(InputIdentity::from).collect(),
        dependencies: BTreeMap::new(),
    })
}

fn snapshot_inputs(inputs: &[VerifiedFile]) -> Result<Vec<ManifestArtifact>> {
    let mut snapshots = inputs
        .iter()
        .map(|input| ManifestArtifact::snapshot(input.artifact_kind.as_str(), &input.path))
        .collect::<Result<Vec<_>>>()?;
    snapshots.sort_by(|a, b| a.path.cmp(&b.path));
    Ok(snapshots)
}

fn scientific_parameters(run: &RunRecord) -> serde_json::Value {
    serde_json::json!({"algorithm_version":3, "layout":run.library_layout, "engine_version":fasterp_version(),
        "dependency_lock":crate::engines::DEPENDENCY_LOCK_SHA256})
}

pub fn processing_identity_matches(run: &RunRecord, manifest: &StageManifest) -> bool {
    manifest.schema_version == 1
        && manifest.stage == "preprocess"
        && manifest.id == run.run_accession
        && manifest.parameters == scientific_parameters(run)
}

/// Verify processing identity/output witnesses without reopening pruned raw inputs.
/// Callers also verify the producer generation and current QC admission.
pub fn retained_processing_compatible(
    run: &RunRecord,
    paths: &ProjectPaths,
    config: &ProjectConfig,
) -> Result<bool> {
    let path = paths
        .preprocess_run_dir(&run.run_accession, &config.quant.engine)
        .join("rnaa_preprocess_manifest.json");
    let Some(manifest) = fs::read(path)
        .ok()
        .and_then(|bytes| serde_json::from_slice::<StageManifest>(&bytes).ok())
    else {
        return Ok(false);
    };
    Ok(manifest.schema_version == 1
        && manifest.stage == "preprocess"
        && manifest.id == run.run_accession
        && manifest.parameters == scientific_parameters(run)
        && !manifest.output_artifacts.is_empty()
        && manifest
            .output_artifacts
            .iter()
            .all(|file| file.matches_file().unwrap_or(false)))
}

fn output_paths(
    run: &RunRecord,
    out_dir: &Path,
) -> Result<(ArtifactKind, PathBuf, Option<PathBuf>)> {
    match run.library_layout {
        LibraryLayout::Single => Ok((
            ArtifactKind::FastqSingle,
            out_dir.join("trimmed.fastq"),
            None,
        )),
        LibraryLayout::Paired => Ok((
            ArtifactKind::FastqR1,
            out_dir.join("trimmed.R1.fastq"),
            Some(out_dir.join("trimmed.R2.fastq")),
        )),
        LibraryLayout::Unknown => bail!("preprocessing requires a known library layout"),
    }
}

/// A read-only lookup. A cache miss never creates directories, stages lanes or
/// invokes the engine against a previously committed output set.
pub fn reconcile_preprocess_artifacts(
    run: &RunRecord,
    fastqs: &[VerifiedFile],
    paths: &ProjectPaths,
    config: &ProjectConfig,
) -> Result<Option<PreprocessArtifacts>> {
    let out_dir = paths.preprocess_run_dir(&run.run_accession, &config.quant.engine);
    let manifest_path = out_dir.join("rnaa_preprocess_manifest.json");
    let Some(manifest) = fs::read(&manifest_path)
        .ok()
        .and_then(|bytes| serde_json::from_slice::<StageManifest>(&bytes).ok())
    else {
        return Ok(None);
    };
    if fastqs.is_empty() || fastqs.iter().any(|file| !file.path.is_file()) {
        return Ok(None);
    }
    if manifest.stage != "preprocess"
        || manifest.id != run.run_accession
        || manifest.parameters != scientific_parameters(run)
        || manifest.input_artifacts != snapshot_inputs(fastqs)?
        || !manifest.verify_local_artifacts()?
    {
        return Ok(None);
    }
    let report_json = out_dir.join("fasterp.json");
    let report: serde_json::Value = match fs::read(&report_json)
        .ok()
        .and_then(|raw| serde_json::from_slice(&raw).ok())
    {
        Some(report) => report,
        None => return Ok(None),
    };
    let (Some(passed), Some(failed)) = (
        report["passed_reads"].as_u64(),
        report["failed_reads"].as_u64(),
    ) else {
        return Ok(None);
    };
    let passed_reads = usize::try_from(passed)?;
    let failed_reads = usize::try_from(failed)?;
    let (kind, read1, read2) = output_paths(run, &out_dir)?;
    if !outputs_ready(&read1, read2.as_deref(), &report_json, passed_reads == 0) {
        return Ok(None);
    }
    let fastqs = build_preprocessed_fastqs(kind, &read1, read2.as_deref())?;
    if manifest.output_artifacts.len() != fastqs.len() + 1
        || fastqs
            .iter()
            .map(|file| file.path.as_path())
            .chain(std::iter::once(report_json.as_path()))
            .any(|path| {
                !manifest
                    .output_artifacts
                    .iter()
                    .any(|file| Path::new(&file.path) == path)
            })
    {
        return Ok(None);
    }
    Ok(Some(PreprocessArtifacts {
        run_accession: run.run_accession.clone(),
        out_dir,
        fastqs,
        report_json,
        tool_name: "fasterp".into(),
        tool_version: fasterp_version(),
        passed_reads,
        failed_reads,
        reused: true,
    }))
}

/// Called only for an explicitly permitted engine-error bypass in an unpublished
/// directory. Remove partial engine outputs so no stale trimmed reads become live.
pub fn bypass(
    run: &RunRecord,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    inputs: &[ManifestArtifact],
    error: &str,
) -> Result<PreprocessArtifacts> {
    if config.quant.preprocess_strict {
        bail!("preprocessing bypass is disabled");
    }
    verify_inputs(inputs)?;
    let out_dir = paths.preprocess_run_dir(&run.run_accession, &config.quant.engine);
    // The CLI supplies a fresh generation. Refuse to erase a canonical legacy path.
    let store = fs::canonicalize(&paths.generations_dir)?;
    if !out_dir.starts_with(&store) {
        bail!("bypass requires an unpublished generation directory");
    }
    let id = out_dir
        .strip_prefix(&store)?
        .components()
        .next()
        .context("generation has no ID")?;
    if store.join(id.as_os_str()).join("seal.json").exists() {
        bail!("bypass cannot replace a sealed generation");
    }
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir)?;
    }
    fs::create_dir_all(&out_dir)?;
    let report_json = out_dir.join("fasterp.json");
    write_json_pretty(
        &report_json,
        &serde_json::json!({"mode":"bypass", "engine_error":error,
        "passed_reads":0, "failed_reads":0, "threads":config.quant.preprocess_threads}),
    )?;
    let started_at = now_rfc3339();
    write_json_pretty(
        &out_dir.join("rnaa_preprocess_manifest.json"),
        &StageManifest {
            schema_version: 1,
            stage: "preprocess".into(),
            id: run.run_accession.clone(),
            started_at: started_at.clone(),
            finished_at: started_at,
            parameters: serde_json::json!({"bypass":true, "engine_error":error}),
            tool_versions: BTreeMap::from([("fasterp".into(), fasterp_version())]),
            input_artifacts: inputs.to_vec(),
            output_artifacts: vec![ManifestArtifact::snapshot(
                "PREPROCESS_REPORT",
                &report_json,
            )?],
            notes: vec!["Explicit engine-error bypass; quantification must use raw FASTQ.".into()],
        },
    )?;
    Ok(PreprocessArtifacts {
        run_accession: run.run_accession.clone(),
        out_dir,
        fastqs: vec![],
        report_json,
        tool_name: "fasterp".into(),
        tool_version: fasterp_version(),
        passed_reads: 0,
        failed_reads: 0,
        reused: false,
    })
}

fn verify_inputs(inputs: &[ManifestArtifact]) -> Result<()> {
    for input in inputs {
        if !input.matches_file()? {
            bail!(
                "preprocessing input changed during execution: {}",
                input.path
            );
        }
    }
    Ok(())
}

fn outputs_ready(read1: &Path, read2: Option<&Path>, report: &Path, allow_empty: bool) -> bool {
    let reads_ready =
        |path: &Path| path.is_file() && (allow_empty || file_size(path).unwrap_or_default() > 0);
    file_size(report).unwrap_or_default() > 0 && reads_ready(read1) && read2.is_none_or(reads_ready)
}

fn build_preprocessed_fastqs(
    read1_kind: ArtifactKind,
    out_read1: &Path,
    out_read2: Option<&Path>,
) -> Result<Vec<VerifiedFile>> {
    let kind = match read1_kind {
        ArtifactKind::FastqSingle => ArtifactKind::FastqTrimmedSingle,
        ArtifactKind::FastqR1 => ArtifactKind::FastqTrimmedR1,
        other => bail!(
            "unsupported read1 kind for fasterp preprocessing: {}",
            other.as_str()
        ),
    };
    let mut outputs = vec![verified_output(kind, out_read1)?];
    if let Some(read2) = out_read2 {
        outputs.push(verified_output(ArtifactKind::FastqTrimmedR2, read2)?);
    }
    Ok(outputs)
}

fn verified_output(kind: ArtifactKind, path: &Path) -> Result<VerifiedFile> {
    Ok(VerifiedFile {
        artifact_kind: kind,
        path: path.into(),
        checksum_type: "sha256".into(),
        checksum: compute_sha256(path)?,
        bytes: file_size(path)?,
        integrity_source: "fasterp_inprocess".into(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rnaa_core::model::LibraryLayout;
    use rnaa_core::state::RunState;

    #[test]
    fn cache_requires_unchanged_source_and_output_hashes() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let paths = ProjectPaths::new(temp.path());
        paths.ensure_layout()?;
        let source = temp.path().join("source.fastq.gz");
        let sequence = "ACGTTCGATCGGATCGATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCTAGCTGACTGATCGATCGTAGCTAGCTGATCGTAGCTAGCTAGCATCGA";
        let read = format!("@read\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()));
        // A single registered file may contain multiple concatenated gzip members.
        // Exercise the actual preprocessor reader rather than validation alone.
        let mut compressed = Vec::new();
        for _ in 0..2 {
            use std::io::Write;
            let mut member =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            member.write_all(read.as_bytes())?;
            compressed.extend(member.finish()?);
        }
        fs::write(&source, compressed)?;
        let run = RunRecord {
            project_id: "test".into(),
            study_accession: "study".into(),
            source_study_accession: None,
            run_accession: "RUN".into(),
            sample_accession: None,
            experiment_accession: None,
            library_layout: LibraryLayout::Single,
            instrument_platform: None,
            instrument_model: None,
            remote_files: vec![],
            metadata: BTreeMap::new(),
            state: RunState::Verified,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        let input = VerifiedFile {
            artifact_kind: ArtifactKind::FastqSingle,
            path: source.clone(),
            checksum_type: "sha256".into(),
            checksum: compute_sha256(&source)?,
            bytes: file_size(&source)?,
            integrity_source: "fixture".into(),
        };
        let mut config = ProjectConfig::default();
        config.quant.preprocess_threads = 1;
        let process = || {
            FasterpInProcessPreprocessor.preprocess(
                &run,
                std::slice::from_ref(&input),
                &paths,
                &config,
            )
        };
        let first = process()?;
        assert!(!first.reused);
        assert_eq!(first.passed_reads, 2);
        let reused = process()?;
        assert!(reused.reused);
        assert_eq!(reused.passed_reads, 2);
        let original = fs::read(&first.fastqs[0].path)?;
        fs::write(&source, read.repeat(3))?;
        assert!(
            reconcile_preprocess_artifacts(&run, std::slice::from_ref(&input), &paths, &config)?
                .is_none()
        );
        assert_eq!(
            fs::read(&first.fastqs[0].path)?,
            original,
            "cache misses must leave earlier outputs untouched"
        );
        let empty = ProjectPaths::new(temp.path().join("uncreated"));
        assert!(
            reconcile_preprocess_artifacts(&run, std::slice::from_ref(&input), &empty, &config)?
                .is_none()
        );
        assert!(
            !empty.root.exists(),
            "read-only reconciliation cannot stage or create outputs"
        );
        let changed = process()?;
        assert!(!changed.reused);
        assert_eq!(changed.passed_reads, 3);
        fs::write(&changed.report_json, b"invalid")?;
        assert!(!process()?.reused);
        fs::write(&changed.fastqs[0].path, b"corrupted")?;
        assert!(!process()?.reused);
        Ok(())
    }
}
