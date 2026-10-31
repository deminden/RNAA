//! Check recorded model designs against the current biological-sample metadata.
use anyhow::{Context, Result, bail};
use rnaa_core::{
    Database, ManifestArtifact, ProjectConfig, ProjectPaths, RunRecord,
    generations::{GenerationHead, OutputGroup},
    sample_plan::SamplePlan,
};
use rnaa_stages::{analysis_request::DeOptions, design::ScientificDesign};
use serde_json::{Value, json};
use std::{collections::BTreeSet, path::PathBuf};

enum Model {
    De(DeOptions),
    Corr {
        adjust: String,
        model: Option<String>,
        design: Option<ScientificDesign>,
        matrix: PathBuf,
        matrix_bytes: u64,
    },
}

pub(super) struct Guard {
    _scratch: tempfile::TempDir,
    filtered: PathBuf,
    sample_sheet: PathBuf,
    runs: Vec<RunRecord>,
    external_matrix_bytes: Option<u64>,
    source: ManifestArtifact,
    filtered_witness: ManifestArtifact,
    models: Vec<Model>,
    config_projection: Value,
    has_de: bool,
    has_corr: bool,
}

fn projection(config: &ProjectConfig, has_de: bool, has_corr: bool) -> Value {
    json!({
        "metadata": config.metadata,
        "de": has_de.then_some(&config.deseq2),
        "corr": has_corr.then(|| json!({"adjust":config.corr.adjust,"model":(config.corr.adjust != "none").then_some(&config.corr.model),"matrix":config.corr.matrix})),
    })
}

impl Guard {
    pub(super) fn external_matrix_bytes(&self) -> Result<u64> {
        self.external_matrix_bytes
            .context("external matrix was not admitted")
    }

    pub(super) fn filtered(&self) -> &std::path::Path {
        &self.filtered
    }

    pub(super) fn correlation_sheet(&self, config: &ProjectConfig) -> Result<&std::path::Path> {
        let formula = if config.corr.adjust == "none" {
            "~ 1"
        } else {
            &config.corr.model
        };
        let plan =
            SamplePlan::from_sheet_with_metadata(&self.filtered, formula, &[], &config.metadata)?;
        plan.write_sheet(&self.sample_sheet)?;
        Ok(&self.sample_sheet)
    }

    pub(super) fn custom_correlation_inputs(
        &self,
        paths: &ProjectPaths,
        db: &Database,
        config: &ProjectConfig,
        matrix: &std::path::Path,
    ) -> Result<crate::analysis::custom_matrix::CustomInputs> {
        let bytes = self.external_matrix_bytes()?;
        let columns = rnaa_formats::matrix::read_matrix_columns_exact(matrix, bytes)?;
        let staging = self._scratch.path().join("custom");
        std::fs::create_dir_all(&staging)?;
        crate::analysis::custom_matrix::prepare(db, paths, config, &self.runs, &columns, &staging)
    }
    pub(super) fn capture(
        paths: &ProjectPaths,
        db: &Database,
        config: &ProjectConfig,
        heads: &[GenerationHead],
        files: &[(GenerationHead, ManifestArtifact, bool)],
        eligible: &BTreeSet<String>,
        external_matrix_bytes: Option<u64>,
    ) -> Result<Self> {
        let runs = db
            .list_active_runs()?
            .into_iter()
            .filter(|run| eligible.contains(&run.run_accession))
            .collect::<Vec<_>>();
        if runs.is_empty() {
            bail!("no currently eligible quantified runs for design validation");
        }
        let scratch = tempfile::Builder::new()
            .prefix("rnaa-export-design-")
            .tempdir_in(&paths.support_dir)?;
        let filtered = scratch.path().join("filtered.tsv");
        let source =
            crate::analysis::inputs::write_filtered_deseq2_samplesheet(paths, &runs, &filtered)?;
        let filtered_witness = ManifestArtifact::snapshot("FILTERED_METADATA", &filtered)?;
        let mut models = Vec::new();
        let mut has_de = false;
        let mut has_corr = false;
        for head in heads {
            let request = db.scientific_request(&head.request_sha256)?.context("selected generation lacks a verified scientific request; cannot check current design")?;
            match head.output_group {
                OutputGroup::De => {
                    has_de = true;
                    models.push(Model::De(
                        serde_json::from_value(request.parameters)
                            .context("invalid recorded DE design")?,
                    ));
                }
                OutputGroup::Corr => {
                    has_corr = true;
                    let adjust = request.parameters["adjust"]
                        .as_str()
                        .context("missing recorded correlation adjustment")?
                        .to_owned();
                    let model = request.parameters["model"].as_str().map(str::to_owned);
                    let design: Option<ScientificDesign> = serde_json::from_value(
                        request
                            .parameters
                            .get("design_contract")
                            .context("missing recorded correlation design")?
                            .clone(),
                    )?;
                    let matrices = files
                        .iter()
                        .filter(|(owner, file, provenance)| {
                            owner.generation_id == head.generation_id
                                && !provenance
                                && file.kind == "ADJUSTED_MATRIX"
                                && std::path::Path::new(&file.path)
                                    .file_name()
                                    .is_some_and(|name| name == "adjusted_matrix.tsv")
                        })
                        .collect::<Vec<_>>();
                    if matrices.len() != 1 {
                        bail!("correlation design validation requires exactly one adjusted matrix");
                    }
                    models.push(Model::Corr {
                        adjust,
                        model,
                        design,
                        matrix: PathBuf::from(&matrices[0].1.path),
                        matrix_bytes: matrices[0].1.bytes,
                    });
                }
                _ => bail!("unsupported output group for design validation"),
            }
        }
        let config_projection = projection(config, has_de, has_corr);
        let sample_sheet = scratch.path().join("samples.tsv");
        let guard = Self {
            _scratch: scratch,
            filtered,
            sample_sheet,
            runs,
            external_matrix_bytes,
            source,
            filtered_witness,
            models,
            config_projection,
            has_de,
            has_corr,
        };
        guard.verify_current(paths, db, config)?;
        Ok(guard)
    }

    pub(super) fn verify_catalog(&self, db: &Database, config: &ProjectConfig) -> Result<()> {
        let persisted: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
        if projection(config, self.has_de, self.has_corr) != self.config_projection
            || projection(&persisted, self.has_de, self.has_corr) != self.config_projection
        {
            bail!("model configuration changed during export");
        }
        Ok(())
    }

    pub(super) fn verify_current(
        &self,
        paths: &ProjectPaths,
        db: &Database,
        config: &ProjectConfig,
    ) -> Result<()> {
        self.verify_catalog(db, config)?;
        let disk = ProjectConfig::load(&paths.config_path)?;
        if projection(&disk, self.has_de, self.has_corr) != self.config_projection {
            bail!("model configuration changed on disk during export");
        }
        if !self.source.matches_file()? || !self.filtered_witness.matches_file()? {
            bail!("sample metadata changed during export");
        }
        for model in &self.models {
            match model {
                Model::De(options) => {
                    if !options
                        .matches_metadata(&self.filtered, config)
                        .context("current DE design or metadata is invalid")?
                    {
                        bail!(
                            "selected DE generation does not match current design or used sample metadata"
                        );
                    }
                }
                Model::Corr {
                    adjust,
                    model,
                    design,
                    matrix,
                    matrix_bytes,
                } => {
                    if adjust != &config.corr.adjust {
                        bail!("selected correlation adjustment does not match current settings");
                    }
                    match adjust.as_str() {
                        "none" if model.is_none() && design.is_none() => {}
                        "residualize"
                            if model.as_deref() == Some(config.corr.model.as_str())
                                && design.is_some() =>
                        {
                            let custom = config
                                .corr
                                .matrix
                                .as_deref()
                                .map(|input| {
                                    self.custom_correlation_inputs(paths, db, config, input)
                                })
                                .transpose()?;
                            let sheet = if let Some(custom) = &custom {
                                custom.sheet.as_path()
                            } else {
                                self.correlation_sheet(config)?
                            };
                            let actual = rnaa_stages::corr::prepare_residualization_exact(
                                matrix,
                                *matrix_bytes,
                                sheet,
                                &config.corr.model,
                                &config.metadata,
                            )?
                            .contract
                            .scientific();
                            if let Some(custom) = custom
                                && !custom.source.matches_file()?
                            {
                                bail!("sample metadata changed during design validation");
                            }
                            if Some(&actual) != design.as_ref() {
                                bail!(
                                    "selected correlation generation does not match current design or used sample metadata"
                                );
                            }
                        }
                        _ => bail!(
                            "selected correlation generation does not match current model settings"
                        ),
                    }
                }
            }
        }
        if !self.source.matches_file()? || !self.filtered_witness.matches_file()? {
            bail!("sample metadata changed during design validation");
        }
        Ok(())
    }

    pub(super) fn proof(&self) -> Value {
        json!({"de_checked":self.has_de,"correlation_checked":self.has_corr,"samplesheet_sha256":self.source.checksum,"full_freshness":false})
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn external_header_is_bounded_and_validates_sample_ids() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let matrix = temp.path().join("matrix.tsv");
        std::fs::write(
            &matrix,
            b"gene_id\t\"sample one\"\tsample-two\nG1\tnot-a-number\t2\n",
        )?;
        let bytes = std::fs::metadata(&matrix)?.len();
        assert_eq!(
            rnaa_formats::matrix::read_matrix_columns_exact(&matrix, bytes)?,
            ["sample one", "sample-two"]
        );
        assert!(rnaa_formats::matrix::read_matrix_columns_exact(&matrix, bytes - 1).is_err());
        std::fs::write(&matrix, b"gene_id\tdup\tdup\nG1\t1\t2\n")?;
        assert!(
            rnaa_formats::matrix::read_matrix_columns_exact(
                &matrix,
                std::fs::metadata(&matrix)?.len()
            )
            .is_err()
        );
        Ok(())
    }

    #[test]
    fn correlation_design_projection_tracks_external_matrix_selection() {
        let mut config = ProjectConfig::default();
        let before = projection(&config, false, true);
        config.corr.matrix = Some("other.tsv".into());
        assert_ne!(projection(&config, false, true), before);
        assert_eq!(
            projection(&config, true, false),
            projection(&ProjectConfig::default(), true, false)
        );
    }
}
