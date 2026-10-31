//! Admission and metadata for user-supplied correlation matrices.
use crate::analysis::inputs::write_filtered_deseq2_samplesheet;
use crate::quality;
use anyhow::{Result, bail};
use rnaa_core::{Database, ManifestArtifact, ProjectConfig, ProjectPaths, RunRecord};
use std::path::{Path, PathBuf};

pub(crate) struct CustomInputs {
    pub sheet: PathBuf,
    pub source: ManifestArtifact,
    pub runs: Vec<RunRecord>,
}

pub(crate) fn prepare(
    db: &Database,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    runs: &[RunRecord],
    columns: &[String],
    staging: &Path,
) -> Result<CustomInputs> {
    let mut runs = runs.to_vec();
    let column_set = columns
        .iter()
        .map(String::as_str)
        .collect::<std::collections::BTreeSet<_>>();
    let mut metadata_source =
        ManifestArtifact::snapshot("SAMPLESHEET_SOURCE", &paths.samplesheet_path())?;
    let mut custom_metadata = None;
    let active = db
        .list_active_runs()?
        .into_iter()
        .map(|run| run.run_accession)
        .collect::<std::collections::BTreeSet<_>>();
    let mut known = std::collections::BTreeMap::new();
    for run in db.list_runs()? {
        let sample = rnaa_core::cohorts::BiologicalSample::for_run(&run)?;
        let included = column_set.contains(run.run_accession.as_str())
            || column_set.contains(sample.id.as_str());
        if included && !active.contains(&run.run_accession) {
            bail!("custom matrix includes inactive run {}", run.run_accession);
        }
        if included && !quality::allowed(db, config, &run.run_accession)? {
            bail!(
                "custom matrix includes QC-ineligible run {}",
                run.run_accession
            );
        }
        if included {
            for alias in [&run.run_accession, &sample.id] {
                if column_set.contains(alias.as_str())
                    && let Some(previous) = known.insert(alias.clone(), sample.id.clone())
                    && previous != sample.id
                {
                    bail!(
                        "custom matrix identifier '{alias}' ambiguously names different biological samples"
                    );
                }
            }
        }
    }
    let mut independent = std::collections::BTreeSet::new();
    for column in columns {
        if let Some(sample) = known.get(column)
            && !independent.insert(sample)
        {
            bail!(
                "custom matrix contains multiple columns for biological sample '{sample}'; pool technical runs before correlation"
            );
        }
    }
    if columns
        .iter()
        .all(|column| known.get(column).is_some_and(|sample| sample == column))
    {
        runs.retain(|run| {
            rnaa_core::cohorts::BiologicalSample::for_run(run)
                .is_ok_and(|sample| column_set.contains(sample.id.as_str()))
        });
        let filtered = staging.join("samplesheet_quantified.tsv");
        metadata_source = write_filtered_deseq2_samplesheet(paths, &runs, &filtered)?;
        let plan = rnaa_core::sample_plan::SamplePlan::from_sheet_with_metadata(
            &filtered,
            if config.corr.adjust == "none" {
                "~ 1"
            } else {
                &config.corr.model
            },
            &[],
            &config.metadata,
        )?;
        let sheet = staging.join("sample_metadata.tsv");
        plan.write_sheet(&sheet)?;
        custom_metadata = Some(sheet);
    }
    Ok(CustomInputs {
        sheet: custom_metadata.unwrap_or_else(|| paths.samplesheet_path()),
        source: metadata_source,
        runs,
    })
}
