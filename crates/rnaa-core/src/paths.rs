use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

#[derive(Debug, Clone)]
pub struct ProjectPaths {
    generation_pin: Option<std::sync::Arc<crate::generations::GenerationStorePin>>,
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
    pub support_dir: PathBuf,
    pub trash_dir: PathBuf,
    pub exports_dir: PathBuf,
    pub generations_dir: PathBuf,
    pub transfers_dir: PathBuf,
    manifest_overrides: std::collections::BTreeMap<String, PathBuf>,
    quant_overrides: std::collections::BTreeMap<String, PathBuf>,
    preprocess_overrides: std::collections::BTreeMap<String, PathBuf>,
    reference_overrides: std::collections::BTreeMap<String, PathBuf>,
    raw_overrides: std::collections::BTreeMap<String, PathBuf>,
    run_manifest_overrides: std::collections::BTreeMap<(String, String), PathBuf>,
}

impl ProjectPaths {
    pub fn new(root: impl AsRef<Path>) -> Self {
        let root = root.as_ref().to_path_buf();
        Self {
            generation_pin: None,
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
            support_dir: root.join(".rnaa_support"),
            trash_dir: root.join("trash"),
            exports_dir: root.join("exports"),
            generations_dir: root.join(".rnaa_generations"),
            transfers_dir: root.join(".rnaa_transfers"),
            manifest_overrides: std::collections::BTreeMap::new(),
            quant_overrides: std::collections::BTreeMap::new(),
            preprocess_overrides: std::collections::BTreeMap::new(),
            reference_overrides: std::collections::BTreeMap::new(),
            raw_overrides: std::collections::BTreeMap::new(),
            run_manifest_overrides: std::collections::BTreeMap::new(),
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
            &self.support_dir,
            &self.support_dir.join("r"),
            &self.trash_dir,
            &self.exports_dir,
            &self.generations_dir,
            &self.transfers_dir,
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
        self.raw_overrides
            .get(run_accession)
            .cloned()
            .unwrap_or_else(|| self.raw_dir.join(run_accession))
    }

    pub fn quant_run_dir(&self, run_accession: &str, engine: &str) -> PathBuf {
        self.quant_overrides
            .get(run_accession)
            .cloned()
            .unwrap_or_else(|| self.quant_dir.join(run_accession))
            .join(engine)
    }

    pub fn preprocess_run_dir(&self, run_accession: &str, engine: &str) -> PathBuf {
        self.preprocess_overrides
            .get(run_accession)
            .cloned()
            .unwrap_or_else(|| {
                self.quant_dir
                    .join(run_accession)
                    .join(engine)
                    .join("fasterp")
            })
    }

    pub fn de_project_dir(&self, project_id: &str) -> PathBuf {
        self.de_dir.join(project_id)
    }

    pub fn corr_project_dir(&self, project_id: &str) -> PathBuf {
        self.corr_dir.join(project_id)
    }

    pub fn reference_dir(&self, reference_id: &str) -> PathBuf {
        self.reference_overrides
            .get(reference_id)
            .cloned()
            .unwrap_or_else(|| self.refs_dir.join(reference_id))
    }

    pub fn manifest_path(&self, stage: &str, id: &str) -> PathBuf {
        if let Some(path) = self.run_manifest_overrides.get(&(stage.into(), id.into())) {
            return path.clone();
        }
        self.manifest_overrides
            .get(stage)
            .cloned()
            .unwrap_or_else(|| self.manifests_dir.join(stage))
            .join(format!("{id}.json"))
    }

    pub fn with_generation_heads(&self, db: &crate::Database) -> Result<Self> {
        let mut paths = Self::new(&self.root).pin_generation_store()?;
        for head in db.generation_heads()? {
            paths = paths.for_generation(&head)?;
        }
        Ok(paths)
    }

    /// Pin one run for an ingest worker without loading every other run's catalog.
    /// Analysis callers resolve the complete project snapshot instead.
    pub fn for_quant_run(&self, db: &crate::Database, run: &str) -> Result<Self> {
        let mut paths = Self::new(&self.root).pin_generation_store()?;
        if let Some(head) = db.generation_head(crate::generations::OutputGroup::Download, run)? {
            paths = paths.for_generation(&head)?;
        }
        if let Some(head) = db.generation_head(crate::generations::OutputGroup::Quant, run)? {
            paths = paths.for_generation(&head)?;
        }
        if let Some(head) = db.generation_head(crate::generations::OutputGroup::Preprocess, run)? {
            paths = paths.for_generation(&head)?;
        }
        Ok(paths)
    }

    pub fn for_reference(&self, db: &crate::Database, reference: &str) -> Result<Self> {
        let paths = Self::new(&self.root).pin_generation_store()?;
        Ok(
            match db.generation_head(crate::generations::OutputGroup::Refs, reference)? {
                Some(head) => paths.for_generation(&head)?,
                None => paths,
            },
        )
    }

    pub fn for_download_transfer(&self, run: &str, request: &str) -> Self {
        let mut paths = Self::new(&self.root);
        let transfer = self.transfers_dir.join("download").join(request).join(run);
        paths.raw_overrides.insert(run.into(), transfer.join("raw"));
        paths.run_manifest_overrides.insert(
            ("download".into(), run.into()),
            transfer.join("manifest.json"),
        );
        paths
    }

    fn pin_generation_store(mut self) -> Result<Self> {
        if self.generation_pin.is_none() {
            self.generation_pin = Some(crate::generations::GenerationStorePin::acquire(
                &self.db_path,
                "generation paths",
            )?);
        }
        Ok(self)
    }
    pub(crate) fn for_generation(&self, head: &crate::generations::GenerationHead) -> Result<Self> {
        let mut paths = self.clone().pin_generation_store()?;
        let connection = rusqlite::Connection::open_with_flags(
            &self.db_path,
            rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
        )?;
        let retired: bool = connection.query_row(
            "SELECT EXISTS(SELECT 1 FROM generation_retirements WHERE generation_id=?1)",
            [&head.generation_id],
            |row| row.get(0),
        )?;
        if retired {
            anyhow::bail!("generation is retired; restore it before opening its paths");
        }
        match head.output_group {
            crate::generations::OutputGroup::Refs => {
                paths.reference_overrides.insert(
                    head.subject.clone(),
                    head.payload_path.join("refs").join(&head.subject),
                );
                paths.run_manifest_overrides.insert(
                    ("refs".into(), head.subject.clone()),
                    head.payload_path
                        .join("manifests/refs")
                        .join(format!("{}.json", head.subject)),
                );
            }
            crate::generations::OutputGroup::Download => {
                paths.raw_overrides.insert(
                    head.subject.clone(),
                    head.payload_path.join("raw").join(&head.subject),
                );
                paths.run_manifest_overrides.insert(
                    ("download".into(), head.subject.clone()),
                    head.payload_path
                        .join("manifests/download")
                        .join(format!("{}.json", head.subject)),
                );
            }
            crate::generations::OutputGroup::Preprocess => {
                paths.preprocess_overrides.insert(
                    head.subject.clone(),
                    head.payload_path
                        .join("preprocess")
                        .join(&head.subject)
                        .join("fasterp"),
                );
            }
            crate::generations::OutputGroup::Quant => {
                paths.quant_overrides.insert(
                    head.subject.clone(),
                    head.payload_path.join("quant").join(&head.subject),
                );
            }
            crate::generations::OutputGroup::De => paths.de_dir = head.payload_path.join("de"),
            crate::generations::OutputGroup::Corr => {
                paths.corr_dir = head.payload_path.join("corr");
                paths
                    .manifest_overrides
                    .insert("corr".into(), head.payload_path.join("manifests/corr"));
            }
        }
        Ok(paths)
    }

    pub fn support_r_script_path(&self, script_name: &str) -> PathBuf {
        self.support_dir.join("r").join(script_name)
    }

    pub fn trash_run_dir(&self, run_accession: &str) -> PathBuf {
        self.trash_dir.join(run_accession)
    }
}
