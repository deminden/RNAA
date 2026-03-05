use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

#[derive(Debug, Clone)]
pub struct ProjectPaths {
    pub root: PathBuf,
    pub db_path: PathBuf,
    pub config_path: PathBuf,
    pub logs_dir: PathBuf,
    pub raw_dir: PathBuf,
    pub quant_dir: PathBuf,
    pub de_dir: PathBuf,
    pub corr_dir: PathBuf,
    pub refs_dir: PathBuf,
    pub manifests_dir: PathBuf,
    pub metadata_dir: PathBuf,
    pub overrides_dir: PathBuf,
    pub trash_dir: PathBuf,
    pub exports_dir: PathBuf,
}

impl ProjectPaths {
    pub fn new(root: impl AsRef<Path>) -> Self {
        let root = root.as_ref().to_path_buf();
        Self {
            db_path: root.join("state.sqlite"),
            config_path: root.join("rnaa.toml"),
            logs_dir: root.join("logs"),
            raw_dir: root.join("raw"),
            quant_dir: root.join("quant"),
            de_dir: root.join("de"),
            corr_dir: root.join("corr"),
            refs_dir: root.join("refs"),
            manifests_dir: root.join("manifests"),
            metadata_dir: root.join("metadata"),
            overrides_dir: root.join("metadata").join("overrides"),
            trash_dir: root.join("trash"),
            exports_dir: root.join("exports"),
            root,
        }
    }

    pub fn ensure_layout(&self) -> Result<()> {
        for dir in [
            &self.root,
            &self.logs_dir,
            &self.raw_dir,
            &self.quant_dir,
            &self.de_dir,
            &self.corr_dir,
            &self.refs_dir,
            &self.manifests_dir,
            &self.manifests_dir.join("resolve"),
            &self.manifests_dir.join("download"),
            &self.manifests_dir.join("refs"),
            &self.manifests_dir.join("quant"),
            &self.manifests_dir.join("deseq2"),
            &self.manifests_dir.join("corr"),
            &self.manifests_dir.join("cleanup"),
            &self.metadata_dir,
            &self.overrides_dir,
            &self.trash_dir,
            &self.exports_dir,
        ] {
            fs::create_dir_all(dir)
                .with_context(|| format!("failed to create {}", dir.display()))?;
        }
        Ok(())
    }

    pub fn samplesheet_path(&self) -> PathBuf {
        self.metadata_dir.join("samplesheet.tsv")
    }

    pub fn samplesheet_schema_path(&self) -> PathBuf {
        self.metadata_dir.join("samplesheet.schema.json")
    }

    pub fn column_map_path(&self) -> PathBuf {
        self.metadata_dir.join("column_map.tsv")
    }

    pub fn raw_run_dir(&self, run_accession: &str) -> PathBuf {
        self.raw_dir.join(run_accession)
    }

    pub fn quant_run_dir(&self, run_accession: &str, engine: &str) -> PathBuf {
        self.quant_dir.join(run_accession).join(engine)
    }

    pub fn de_project_dir(&self, project_id: &str) -> PathBuf {
        self.de_dir.join(project_id)
    }

    pub fn corr_project_dir(&self, project_id: &str) -> PathBuf {
        self.corr_dir.join(project_id)
    }

    pub fn reference_dir(&self, reference_id: &str) -> PathBuf {
        self.refs_dir.join(reference_id)
    }

    pub fn manifest_path(&self, stage: &str, id: &str) -> PathBuf {
        self.manifests_dir.join(stage).join(format!("{id}.json"))
    }

    pub fn trash_run_dir(&self, run_accession: &str) -> PathBuf {
        self.trash_dir.join(run_accession)
    }
}
