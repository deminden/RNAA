//! Owned retention of consumed FASTQs with durable removal evidence.
use crate::cli::fastq_retention_label;
use anyhow::Context;
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::config::{FastqRetention, ProjectConfig};
use rnaa_core::generations::OutputGroup;
use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::model::ArtifactKind;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::RunState;
use rnaa_core::util::compute_sha256;
use serde_json::json;
use std::fs;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};

fn retention_kinds(policy: FastqRetention) -> &'static [ArtifactKind] {
    match policy {
        FastqRetention::Both => &[],
        FastqRetention::Original => &[
            ArtifactKind::FastqTrimmedR1,
            ArtifactKind::FastqTrimmedR2,
            ArtifactKind::FastqTrimmedSingle,
        ],
        FastqRetention::Trimmed => &[
            ArtifactKind::FastqR1,
            ArtifactKind::FastqR2,
            ArtifactKind::FastqSingle,
        ],
        FastqRetention::None => &[
            ArtifactKind::FastqR1,
            ArtifactKind::FastqR2,
            ArtifactKind::FastqSingle,
            ArtifactKind::FastqTrimmedR1,
            ArtifactKind::FastqTrimmedR2,
            ArtifactKind::FastqTrimmedSingle,
        ],
    }
}

pub(crate) fn run_fastq_retention(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    run_accession: &str,
) -> Result<()> {
    // Explicit reversible FASTQ cleanup takes precedence over permanent retention.
    if matches!(
        config.quant.cleanup,
        rnaa_core::CleanupMode::Fastq | rnaa_core::CleanupMode::All
    ) && config.quant.cleanup_on != rnaa_core::CleanupWhen::Never
    {
        return Ok(());
    }
    let kinds = retention_kinds(config.quant.fastq_retention);
    if kinds.is_empty()
        || !db
            .list_artifacts_for_run(Some(run_accession))?
            .iter()
            .any(|artifact| kinds.contains(&artifact.kind))
    {
        return Ok(());
    }
    let Some(mut lease) = JobLease::try_acquire(
        db,
        "retention",
        Some(run_accession),
        &[
            RunState::QuantDone,
            RunState::DeRunning,
            RunState::DeDone,
            RunState::DeFailed,
            RunState::CorrRunning,
            RunState::CorrDone,
            RunState::CorrFailed,
            RunState::Cleaned,
        ],
    )?
    else {
        // A replacement producer owns the resource; reconciliation will retry the
        // policy against its committed result, never against this older snapshot.
        return Ok(());
    };
    let owned = lease.database();
    let result = (|| {
        let paths = paths.for_quant_run(&owned, run_accession)?;
        let manifest = paths
            .quant_run_dir(run_accession, &config.quant.engine)
            .join("rnaa_quant_manifest.json");
        let digest = rnaa_core::fingerprint::StageRequest {
            schema_version: 1,
            stage: "retention".into(),
            algorithm_version: 1,
            engines: Default::default(),
            parameters: json!({ "policy": fastq_retention_label(config.quant.fastq_retention) }),
            inputs: vec![],
            dependencies: std::collections::BTreeMap::from([(
                "quant_manifest".into(),
                compute_sha256(&manifest)?,
            )]),
        }
        .digest()?;
        owned.set_job_request(&digest)?;
        apply_fastq_retention(&owned, &paths, config, run_accession)
    })();
    match result {
        Ok(()) => lease.finish(JobOutcome::Succeeded, None),
        Err(error) => {
            lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))?;
            Err(error)
        }
    }
}

fn apply_fastq_retention(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    run_accession: &str,
) -> Result<()> {
    if db
        .generation_head(OutputGroup::Quant, run_accession)?
        .is_some()
        && !db.verify_generation(OutputGroup::Quant, run_accession)?
    {
        bail!("refusing retention against an invalid quantification generation");
    }
    let artifacts = db.list_artifacts_for_run(Some(run_accession))?;
    let remove_kinds = retention_kinds(config.quant.fastq_retention);

    let selected = artifacts
        .into_iter()
        .filter(|artifact| remove_kinds.contains(&artifact.kind))
        .collect::<Vec<_>>();
    if selected.is_empty() {
        return Ok(());
    }

    let quant_manifest = paths
        .quant_run_dir(run_accession, &config.quant.engine)
        .join("rnaa_quant_manifest.json");
    let generation = compute_sha256(&quant_manifest)?;
    let manifest: rnaa_core::StageManifest = serde_json::from_slice(&fs::read(&quant_manifest)?)?;
    let binding = db.generation_input_binding(OutputGroup::Quant, run_accession)?;
    let input_witnesses = binding
        .as_ref()
        .map(|binding| binding.inputs.as_slice())
        .unwrap_or(&manifest.input_artifacts);
    let mut deleted = Vec::new();
    for artifact in selected {
        let source = PathBuf::from(&artifact.path);
        let Some(snapshot) = input_witnesses
            .iter()
            .find(|item| item.path == artifact.path)
        else {
            bail!(
                "refusing retention of input absent from quant provenance: {}",
                artifact.path
            );
        };
        if !source.exists() {
            let authorized = db
                .quant_retention_evidence(run_accession)?
                .iter()
                .any(|record| {
                    record.generation_sha256 == generation && record.artifact == *snapshot
                });
            if !authorized {
                bail!(
                    "FASTQ disappeared without retention evidence: {}",
                    artifact.path
                );
            }
            db.record_fastq_retention(run_accession, &generation, snapshot)?;
            db.delete_artifact(&artifact.path)?;
            continue;
        }
        if !snapshot.matches_file()? {
            bail!(
                "FASTQ changed since quantification; retaining {}",
                artifact.path
            );
        }
        db.record_fastq_retention(run_accession, &generation, snapshot)?;
        fs::remove_file(&source)
            .with_context(|| format!("failed to permanently delete {}", source.display()))?;
        db.delete_artifact(&artifact.path)?;
        prune_empty_parent_dirs(&source, &paths.root)?;
        deleted.push(source.display().to_string());
    }
    db.append_event(
        "retention",
        Some(run_accession),
        "FASTQ retention policy applied",
        json!({
            "policy": fastq_retention_label(config.quant.fastq_retention),
            "deleted": deleted,
        }),
    )?;
    Ok(())
}

fn prune_empty_parent_dirs(path: &Path, stop_at: &Path) -> Result<()> {
    let mut current = path.parent();
    while let Some(dir) = current {
        if dir == stop_at || !dir.starts_with(stop_at) {
            break;
        }
        match fs::remove_dir(dir) {
            Ok(()) => {
                current = dir.parent();
            }
            Err(err) if err.kind() == ErrorKind::DirectoryNotEmpty => break,
            Err(err) if err.kind() == ErrorKind::NotFound => {
                current = dir.parent();
            }
            Err(err) => {
                return Err(err)
                    .with_context(|| format!("failed to prune empty dir {}", dir.display()));
            }
        }
    }
    Ok(())
}
