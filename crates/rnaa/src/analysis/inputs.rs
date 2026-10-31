//! Analysis input admission, metadata snapshots and design checks.
use crate::ingest::quant::quant_generations_valid;
use anyhow::Context;
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::generations::OutputGroup;
use rnaa_core::model::{ArtifactKind, RunRecord};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::RunState;
use rnaa_core::util::file_size;
use rnaa_stages::quant::external_quant_matches;
use serde_json::json;
use std::path::{Path, PathBuf};

pub(crate) fn validate_deseq2_quant_inputs(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::ReferenceBundle,
    runs: &[RunRecord],
) -> Result<Vec<RunRecord>> {
    let mut valid = Vec::new();
    for queued in runs {
        let accepted = db
            .with_idle_run(&queued.run_accession, |db| {
                let Some(run) = db.get_run(&queued.run_accession)? else {
                    return Ok(None);
                };
                if !rnaa_core::eligibility::quantified(run.state) {
                    return Ok(None);
                }
                validate_quant_run(db, paths, config, reference, &run)
            })?
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "run {} is busy; retry analysis after its worker finishes",
                    queued.run_accession
                )
            })?;
        if let Some(run) = accepted {
            valid.push(run);
        }
    }
    Ok(valid)
}
fn validate_quant_run(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    reference: &rnaa_core::ReferenceBundle,
    run: &RunRecord,
) -> Result<Option<RunRecord>> {
    // Report corruption affects admission, not the validity of completed counts.
    // The cohort snapshot below records the blocked decision independently.
    if !rnaa_core::eligibility::assess(db, config, &run.run_accession)?
        .descriptor
        .eligible()
    {
        return Ok(None);
    }

    let current = paths.for_quant_run(db, &run.run_accession)?;
    let dir = current.quant_run_dir(&run.run_accession, &config.quant.engine);
    if dir.join("rnaa_quant_import.json").exists() {
        if external_quant_matches(&run.run_accession, &dir, reference).unwrap_or(false) {
            return Ok(Some(run.clone()));
        } else {
            invalidate_quant_input(
                db,
                &run.run_accession,
                "external abundance declaration, reference or payload changed",
            )?;
        }
        return Ok(None);
    }
    let run_artifacts = db.list_artifacts_for_run(Some(&run.run_accession))?;
    let abundance_h5 = run_artifacts
        .iter()
        .find(|artifact| artifact.kind == ArtifactKind::QuantAbundanceH5)
        .map(|artifact| PathBuf::from(&artifact.path))
        .unwrap_or_else(|| {
            paths
                .quant_run_dir(&run.run_accession, &config.quant.engine)
                .join("abundance.h5")
        });
    let run_info_json = run_artifacts
        .iter()
        .find(|artifact| artifact.kind == ArtifactKind::QuantRunInfo)
        .map(|artifact| PathBuf::from(&artifact.path))
        .unwrap_or_else(|| {
            paths
                .quant_run_dir(&run.run_accession, &config.quant.engine)
                .join("run_info.json")
        });

    let missing = [
        (&abundance_h5, "abundance.h5"),
        (&run_info_json, "run_info.json"),
    ]
    .into_iter()
    .filter_map(|(path, label)| {
        if !path.exists() || file_size(path).unwrap_or_default() == 0 {
            Some(format!("{label} missing or empty at {}", path.display()))
        } else {
            None
        }
    })
    .collect::<Vec<_>>();
    let mut missing = missing;
    if !rnaa_stages::quant::source_binding_current(db, run, config)? {
        missing.push(
            "download source selection changed; restaging and quantification are required".into(),
        );
    }
    if !quant_generations_valid(db, &run.run_accession, false)? {
        missing.push("quantification or download generation failed verification".into());
    }
    if let Some(head) = db.generation_head(OutputGroup::Quant, &run.run_accession)?
        && (!abundance_h5.starts_with(&head.payload_path)
            || !run_info_json.starts_with(&head.payload_path))
    {
        missing.push("quantification records do not belong to the current generation".into());
    }

    if missing.is_empty()
        && crate::ingest::quant::reconcile_current_quant_artifacts(
            db, run, reference, &current, config,
        )?
        .is_none()
    {
        missing.push("quantification request does not match the current reference, scientific options or declared input set; requantify this run".into());
    }
    if missing.is_empty() {
        return Ok(Some(run.clone()));
    }

    let message = format!("quant outputs invalid for DE stage: {}", missing.join("; "));
    invalidate_quant_input(db, &run.run_accession, &message)?;
    Ok(None)
}

fn invalidate_quant_input(db: &Database, run: &str, message: &str) -> Result<()> {
    db.atomic(|db| {
        db.set_run_state(run, RunState::QuantFailed, Some(message))?;
        for stage in ["normalize", "deseq2", "corr"] {
            if db.get_project_stage_state(stage)?.is_some() {
                db.set_project_stage_state(
                    stage,
                    rnaa_core::ProjectStageState::Pending,
                    Some(message),
                )?;
            }
        }
        db.append_event(
            "deseq2",
            Some(run),
            "run excluded from DE due to invalid quant outputs",
            json!({"error":message}),
        )
    })
}

pub(crate) fn write_filtered_deseq2_samplesheet(
    paths: &ProjectPaths,
    runs: &[RunRecord],
    output_path: &Path,
) -> Result<rnaa_core::ManifestArtifact> {
    let source =
        rnaa_core::ManifestArtifact::snapshot("SAMPLESHEET_SOURCE", &paths.samplesheet_path())?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&source.path)?;
    let mut headers = reader
        .headers()?
        .iter()
        .map(str::to_string)
        .collect::<Vec<_>>();
    if headers
        .iter()
        .collect::<std::collections::BTreeSet<_>>()
        .len()
        != headers.len()
    {
        bail!("samplesheet contains duplicate column names");
    }
    let run_column = headers
        .iter()
        .position(|name| name == "run_accession")
        .context("samplesheet missing run_accession")?;
    let mut rows = std::collections::BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let id = row.get(run_column).unwrap_or_default().trim();
        if id.is_empty()
            || rows
                .insert(
                    id.to_string(),
                    row.iter().map(str::to_string).collect::<Vec<_>>(),
                )
                .is_some()
        {
            bail!("samplesheet contains duplicate or empty run ID '{id}'");
        }
    }
    let original_columns = headers.len();
    for field in ["biological_sample_id", "identity_source"] {
        if !headers.iter().any(|name| name == field) {
            headers.push(field.into());
        }
    }
    let bio_column = headers
        .iter()
        .position(|name| name == "biological_sample_id")
        .unwrap();
    let identity_column = headers
        .iter()
        .position(|name| name == "identity_source")
        .unwrap();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_path)?;
    writer.write_record(&headers)?;
    for run in runs {
        let mut row = rows.remove(&run.run_accession).with_context(|| {
            format!(
                "active run '{}' is missing from samplesheet; resolve or refresh its metadata",
                run.run_accession
            )
        })?;
        row.resize(headers.len(), String::new());
        let identity = rnaa_core::cohorts::BiologicalSample::for_run(run)?;
        if bio_column < original_columns
            && !row[bio_column].trim().is_empty()
            && row[bio_column].trim() != identity.id
        {
            bail!(
                "samplesheet biological identity for {} differs from its registered identity; update metadata and resolve again",
                run.run_accession
            );
        }
        row[bio_column] = identity.id;
        row[identity_column] = identity.source.as_str().into();
        writer.write_record(&row)?;
    }
    writer.flush()?;
    if !source.matches_file()? {
        bail!("samplesheet changed while capturing analysis metadata; retry analysis");
    }
    Ok(source)
}

pub(crate) fn bind_metadata<P>(
    request: &mut rnaa_stages::analysis_request::AnalysisRequest<P>,
    source: &rnaa_core::ManifestArtifact,
) -> Result<()> {
    request.inputs.push(source.clone());
    request.observe_input(source.clone());
    Ok(())
}

pub(crate) fn freeze_metadata_source<P>(
    request: &mut rnaa_stages::analysis_request::AnalysisRequest<P>,
    source: &rnaa_core::ManifestArtifact,
    root: &Path,
) -> Result<()> {
    let dir = root.join("metadata_source");
    std::fs::create_dir_all(&dir)?;
    let target = dir.join(
        Path::new(&source.path)
            .file_name()
            .context("metadata source has no file name")?,
    );
    std::fs::copy(&source.path, &target)?;
    request.relocate_input("SAMPLESHEET_SOURCE", &target)?;
    request.verify_inputs()
}

pub(crate) fn read_samplesheet_header(path: &Path) -> Result<std::collections::BTreeSet<String>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read samplesheet {}", path.display()))?;
    let headers = reader
        .headers()
        .with_context(|| format!("failed to parse samplesheet header {}", path.display()))?;
    Ok(headers.iter().map(ToString::to_string).collect())
}

pub(crate) fn validate_design_formula(
    header: &std::collections::BTreeSet<String>,
    design: &str,
    schema: &rnaa_core::metadata::MetadataConfig,
) -> Result<()> {
    rnaa_core::metadata::variables(design, &header.iter().cloned().collect::<Vec<_>>(), schema)?;
    Ok(())
}
