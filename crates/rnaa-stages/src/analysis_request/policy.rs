//! Configuration observations fence publication without putting execution paths
//! or thread allocations into scientific request identities.
use super::AnalysisRequest;
use anyhow::{Result, bail};
use rnaa_core::{Database, ProjectConfig, metadata::MetadataConfig};

#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub(super) struct QuantObservation {
    engine: String,
    strand: rnaa_core::config::QuantStrand,
    fragment_length: f64,
    fragment_sd: f64,
    preprocess: bool,
    source_preference: rnaa_core::config::DownloadPreference,
}
impl QuantObservation {
    fn capture(config: &ProjectConfig) -> Self {
        Self {
            engine: config.quant.engine.clone(),
            strand: config.quant.strand,
            fragment_length: config.quant.single_fragment_length,
            fragment_sd: config.quant.single_fragment_sd,
            preprocess: config.quant.preprocess,
            source_preference: config.download.prefer,
        }
    }
}

#[derive(Debug, Clone, serde::Serialize)]
#[serde(tag = "stage", rename_all = "snake_case")]
pub(super) enum ObservedConfig {
    De {
        quant: QuantObservation,
        schema: MetadataConfig,
        design: String,
        transform: String,
        counts: String,
    },
    Corr {
        quant: QuantObservation,
        schema: MetadataConfig,
        adjust: String,
        model: String,
        method: String,
        geneset: String,
        output: String,
        matrix: Option<std::path::PathBuf>,
        modules: Option<rnaa_core::config::CorrModuleConfig>,
    },
}
impl ObservedConfig {
    pub(super) fn de(config: &ProjectConfig, design: &str) -> Self {
        Self::De {
            quant: QuantObservation::capture(config),
            schema: config.metadata.clone(),
            design: design.into(),
            transform: config.deseq2.transform.clone(),
            counts: config.deseq2.counts_from_abundance.clone(),
        }
    }
    pub(super) fn corr(config: &ProjectConfig) -> Self {
        Self::Corr {
            quant: QuantObservation::capture(config),
            schema: config.metadata.clone(),
            adjust: config.corr.adjust.clone(),
            model: config.corr.model.clone(),
            method: config.corr.method.clone(),
            geneset: config.corr.geneset.clone(),
            output: config.corr.output.clone(),
            matrix: config.corr.matrix.clone(),
            modules: config.corr.modules.clone(),
        }
    }
    fn matches(&self, config: &ProjectConfig) -> bool {
        match self {
            Self::De {
                quant,
                schema,
                design,
                transform,
                counts,
            } => {
                *quant == QuantObservation::capture(config)
                    && *schema == config.metadata
                    && *design == config.deseq2.design
                    && *transform == config.deseq2.transform
                    && *counts == config.deseq2.counts_from_abundance
            }
            Self::Corr {
                quant,
                schema,
                adjust,
                model,
                method,
                geneset,
                output,
                matrix,
                modules,
            } => {
                *quant == QuantObservation::capture(config)
                    && *schema == config.metadata
                    && *adjust == config.corr.adjust
                    && *model == config.corr.model
                    && *method == config.corr.method
                    && *geneset == config.corr.geneset
                    && *output == config.corr.output
                    && *matrix == config.corr.matrix
                    && *modules == config.corr.modules
            }
        }
    }
}
impl<P> AnalysisRequest<P> {
    pub fn configuration_observation(&self) -> Result<serde_json::Value> {
        Ok(serde_json::to_value(&self.policy)?)
    }
    pub fn verify_settings_catalog(&self, db: &Database) -> Result<()> {
        if let Some(policy) = &self.policy {
            let config: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
            config.validate()?;
            if !policy.matches(&config) {
                bail!(
                    "model schema or scientific settings changed during analysis; retry with the current configuration"
                );
            }
        }
        Ok(())
    }
    pub fn verify_current_settings(&self, db: &Database) -> Result<()> {
        self.verify_settings_catalog(db)?;
        if let Some(policy) = &self.policy {
            let config = ProjectConfig::load(
                &std::path::Path::new(&db.project()?.root_dir).join("rnaa.toml"),
            )?;
            if !policy.matches(&config) {
                bail!(
                    "model schema or scientific settings changed on disk during analysis; retry with the current configuration"
                );
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
