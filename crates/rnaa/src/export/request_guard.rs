//! Rebuild supported current scientific requests instead of inferring freshness.
use super::{admission, design};
use anyhow::{Context, Result, bail};
use rnaa_core::{
    ContrastSpec, Database, ManifestArtifact, ProjectConfig, ProjectPaths, ReferenceBundle,
    RunRecord,
    generations::{GenerationHead, OutputGroup},
};
use serde_json::{Value, json};
use std::path::Path;

pub(super) struct Guard {
    heads: Vec<GenerationHead>,
    runs: Vec<RunRecord>,
    reference: ReferenceBundle,
    config_projection: Value,
    contrasts: Vec<Value>,
    checks_de_contrasts: bool,
    external_matrix: bool,
}

fn projection(config: &ProjectConfig) -> Value {
    json!({"quant":config.quant,"download_preference":config.download.prefer,"metadata":config.metadata,"de":config.deseq2,"corr":config.corr,"refs":config.refs})
}

fn contrast_values(contrasts: &[ContrastSpec]) -> Result<Vec<Value>> {
    let mut values = contrasts
        .iter()
        .map(serde_json::to_value)
        .collect::<std::result::Result<Vec<_>, _>>()?;
    values.sort_by(|left, right| left["name"].as_str().cmp(&right["name"].as_str()));
    Ok(values)
}

impl Guard {
    pub(super) fn capture(
        paths: &ProjectPaths,
        db: &Database,
        config: &ProjectConfig,
        heads: &[GenerationHead],
        admission: &admission::Guard,
        design: &design::Guard,
    ) -> Result<Self> {
        let selected = admission.reference();
        let reference_paths = selected.generation.paths(paths)?;
        let reference = ReferenceBundle {
            id: selected.reference_id.clone(),
            organism: config.refs.organism.clone(),
            ensembl_release: config.refs.ensembl.clone(),
            cdna_path: selected.cdna.path.clone().into(),
            gtf_path: selected.gtf.path.clone().into(),
            kallisto_index_path: selected.index.path.clone().into(),
            tx2gene_path: selected.mapping.path.clone().into(),
            gene_annotation_path: Some(selected.annotation.path.clone().into()),
            manifest_path: reference_paths.manifest_path("refs", &selected.reference_id),
        };
        let runs = db
            .list_active_runs()?
            .into_iter()
            .filter(|run| admission.eligible().contains(&run.run_accession))
            .collect::<Vec<_>>();
        if runs.is_empty() {
            bail!("no eligible quantified runs for current request validation");
        }
        let checks_de_contrasts = heads
            .iter()
            .any(|head| head.output_group == OutputGroup::De && head.stage == "deseq2");
        let contrasts = if checks_de_contrasts {
            contrast_values(&db.list_contrasts()?)?
        } else {
            Vec::new()
        };
        let guard = Self {
            heads: heads.to_vec(),
            runs,
            reference,
            config_projection: projection(config),
            contrasts,
            checks_de_contrasts,
            external_matrix: config.corr.matrix.is_some()
                && heads
                    .iter()
                    .any(|head| head.output_group == OutputGroup::Corr),
        };
        guard.verify_current(paths, db, config, admission, design)?;
        Ok(guard)
    }

    pub(super) fn verify_catalog(&self, db: &Database, config: &ProjectConfig) -> Result<()> {
        let persisted: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
        if projection(config) != self.config_projection
            || projection(&persisted) != self.config_projection
        {
            bail!("scientific configuration changed during current request validation");
        }
        if self.checks_de_contrasts && contrast_values(&db.list_contrasts()?)? != self.contrasts {
            bail!("contrast selection changed during current request validation");
        }
        for head in &self.heads {
            if db
                .generation_head(head.output_group, &head.subject)?
                .as_ref()
                != Some(head)
            {
                bail!("selected generation changed during current request validation");
            }
        }
        for run in &self.runs {
            rnaa_stages::quant::verify_quant_catalog(db, run, config)?;
        }
        Ok(())
    }

    fn validate_quant(
        &self,
        paths: &ProjectPaths,
        db: &Database,
        config: &ProjectConfig,
    ) -> Result<()> {
        for run in &self.runs {
            db.with_idle_run(&run.run_accession, |db| {
                rnaa_stages::quant::verify_quant_current(db, run, config)?;
                let current = paths.for_quant_run(db, &run.run_accession)?;
                let directory = current.quant_run_dir(&run.run_accession, &config.quant.engine);
                let valid = if directory.join("rnaa_quant_import.json").exists() {
                    rnaa_stages::quant::external_quant_matches(&run.run_accession, &directory, &self.reference)?
                } else {
                    for group in [OutputGroup::Quant, OutputGroup::Download, OutputGroup::Preprocess] {
                        if group == OutputGroup::Preprocess && !config.quant.preprocess { continue; }
                        if let Some(head) = db.generation_head(group, &run.run_accession)?
                            && !db.generation_lineage(&head)?.recorded {
                            bail!("{} generation for {} has unresolved input lineage", group.as_str(), run.run_accession);
                        }
                    }
                    db.generation_head(OutputGroup::Quant, &run.run_accession)?.is_some()
                        && crate::ingest::quant::quant_generations_valid(db, &run.run_accession, config.quant.preprocess)?
                        && rnaa_stages::quant::source_binding_current(db, run, config)?
                        && crate::ingest::quant::reconcile_current_quant_artifacts(db, run, &self.reference, &current, config)?.is_some()
                };
                if !valid { bail!("quantification for {} does not match current scientific inputs and settings", run.run_accession); }
                Ok(())
            })?.with_context(|| format!("run {} is busy; retry current request validation", run.run_accession))?;
        }
        Ok(())
    }

    pub(super) fn verify_current(
        &self,
        paths: &ProjectPaths,
        db: &Database,
        config: &ProjectConfig,
        admission: &admission::Guard,
        design: &design::Guard,
    ) -> Result<()> {
        self.verify_catalog(db, config)?;
        let disk = ProjectConfig::load(&paths.config_path)?;
        if projection(&disk) != self.config_projection {
            bail!("scientific configuration changed on disk during current request validation");
        }
        if !db
            .generation_lineage(&admission.reference().generation)?
            .recorded
        {
            bail!(
                "selected reference has unresolved input lineage; cannot establish current request"
            );
        }
        self.validate_quant(paths, db, config)?;
        let current_paths = paths.with_generation_heads(db)?;
        for head in &self.heads {
            if !db.generation_lineage(head)?.recorded {
                bail!(
                    "selected {} generation has unresolved input lineage; cannot establish current request",
                    head.output_group.as_str()
                );
            }
        }
        let de = self
            .heads
            .iter()
            .find(|head| head.output_group == OutputGroup::De);
        if let Some(head) = de {
            let saved = db
                .scientific_request(&head.request_sha256)?
                .context("DE scientific request missing")?;
            let contrasts: Vec<ContrastSpec> = serde_json::from_value(
                saved
                    .parameters
                    .get("contrasts")
                    .context("DE contrast selection missing")?
                    .clone(),
            )?;
            if head.stage == "deseq2" && contrast_values(&contrasts)? != self.contrasts {
                bail!("selected DE generation does not match current contrast selection");
            }
            if head.stage == "normalize" && !contrasts.is_empty() {
                bail!("normalization request unexpectedly contains contrasts");
            }
            let mut current = rnaa_stages::analysis_request::deseq_request(
                &current_paths,
                design.filtered(),
                &self.reference,
                &config.deseq2.design,
                &contrasts,
                config,
                &head.stage,
            )?;
            admission.bind(&mut current)?;
            current.verify_current_settings(db)?;
            current.verify_inputs()?;
            if current.request.digest()? != head.request_sha256 {
                bail!(
                    "selected DE generation does not match the reconstructed current scientific request"
                );
            }
        }
        for head in self
            .heads
            .iter()
            .filter(|head| head.output_group == OutputGroup::Corr)
        {
            let saved = db
                .scientific_request(&head.request_sha256)?
                .context("correlation scientific request missing")?;
            let recorded_modules = saved
                .parameters
                .get("modules")
                .context("correlation module settings missing")?;
            let modules = config
                .corr
                .modules
                .as_ref()
                .map(rnaa_stages::corr::ModuleGseaParams::from);
            if recorded_modules.is_null() != modules.is_none() {
                bail!(
                    "current module-analysis settings do not match the selected result; older results may need correlation to be rerun"
                );
            }
            let custom = config
                .corr
                .matrix
                .as_deref()
                .map(|matrix| design.custom_correlation_inputs(paths, db, config, matrix))
                .transpose()?;
            let (matrix_path, matrix_bytes) = if let Some(matrix) = &config.corr.matrix {
                (matrix.clone(), design.external_matrix_bytes()?)
            } else {
                let producer =
                    de.context("current request validation requires a selected DE matrix")?;
                let candidates = db
                    .generation_output_artifacts(producer)?
                    .into_iter()
                    .filter(|file| {
                        file.kind == "VST"
                            && Path::new(&file.path)
                                .file_name()
                                .is_some_and(|name| name == "vst.tsv")
                    })
                    .collect::<Vec<ManifestArtifact>>();
                if candidates.len() != 1 {
                    bail!("current request validation requires exactly one selected DE VST matrix");
                }
                let matrix = &candidates[0];
                let inputs = saved
                    .inputs
                    .iter()
                    .filter(|input| input.role == "MATRIX")
                    .collect::<Vec<_>>();
                if inputs.len() != 1
                    || inputs[0].sha256 != matrix.checksum
                    || inputs[0].bytes != matrix.bytes
                    || inputs[0].name != "vst.tsv"
                {
                    bail!("selected correlation input does not match the current DE matrix");
                }
                (matrix.path.clone().into(), matrix.bytes)
            };
            let sheet = if let Some(custom) = &custom {
                custom.sheet.as_path()
            } else {
                design.correlation_sheet(config)?
            };
            let mut current = rnaa_stages::analysis_request::correlation_request_exact(
                &matrix_path,
                sheet,
                self.reference.gene_annotation_path.as_deref(),
                config,
                modules.as_ref(),
                matrix_bytes,
            )?;
            if let Some(custom) = &custom {
                crate::analysis::inputs::bind_metadata(&mut current, &custom.source)?;
            }
            admission.bind(&mut current)?;
            current.verify_current_settings(db)?;
            current.verify_inputs()?;
            if current.request.digest()? != head.request_sha256 {
                bail!(
                    "selected correlation generation does not match the reconstructed current scientific request"
                );
            }
        }
        Ok(())
    }

    pub(super) fn proof(&self) -> Value {
        json!({"requests":self.heads.iter().map(|head| json!({"output_group":head.output_group.as_str(),"generation_id":head.generation_id,"request_sha256":head.request_sha256})).collect::<Vec<_>>(),"scope":"current supported DE and correlation requests","external_files_locked":false,"external_matrix":self.external_matrix.then(||json!({"bytes_match_request":true,"upstream_provenance_verified":false}))})
    }
}
