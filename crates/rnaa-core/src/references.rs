//! Selected reference identity and stage-specific admission, independent of paths.
use crate::{
    Database, ManifestArtifact, ProjectConfig, ReferenceBundle,
    fingerprint::{InputIdentity, StageRequest},
    generations::{GenerationHead, OutputGroup},
    util::now_rfc3339,
};
use anyhow::{Context, Result, bail};
use rusqlite::{Connection, OptionalExtension, params};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, path::Path};

/// Opaque compare-and-swap identity; deliberately distinct from scientific keys.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceChoiceToken(String);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceScope {
    Index,
    Analysis,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferencePolicy {
    organism: String,
    release: String,
    cdna: String,
    gtf: String,
    k: usize,
}
impl ReferencePolicy {
    pub fn from_config(config: &ProjectConfig) -> Self {
        Self {
            organism: config.refs.organism.clone(),
            release: config.refs.ensembl.clone(),
            cdna: config.refs.custom_cdna.clone(),
            gtf: config.refs.custom_gtf.clone(),
            k: config.refs.index_k,
        }
    }
    fn matches(&self, other: &Self, scope: ReferenceScope) -> bool {
        self.organism == other.organism
            && self.release == other.release
            && self.cdna == other.cdna
            && self.k == other.k
            && (scope == ReferenceScope::Index || self.gtf == other.gtf)
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceSelection {
    pub schema_version: u32,
    pub reference_id: String,
    pub generation: GenerationHead,
    pub policy: ReferencePolicy,
    pub cdna: ManifestArtifact,
    pub gtf: ManifestArtifact,
    pub index: ManifestArtifact,
    pub mapping: ManifestArtifact,
    pub annotation: ManifestArtifact,
}
impl ReferenceSelection {
    pub fn capture(
        bundle: &ReferenceBundle,
        generation: GenerationHead,
        config: &ProjectConfig,
    ) -> Result<Self> {
        Ok(Self {
            schema_version: 1,
            reference_id: bundle.id.clone(),
            generation,
            policy: ReferencePolicy::from_config(config),
            cdna: ManifestArtifact::snapshot("REFERENCE_CDNA", &bundle.cdna_path)?,
            gtf: ManifestArtifact::snapshot("REFERENCE_GTF", &bundle.gtf_path)?,
            index: ManifestArtifact::snapshot("REFERENCE_INDEX", &bundle.kallisto_index_path)?,
            mapping: ManifestArtifact::snapshot("REFERENCE_TX2GENE", &bundle.tx2gene_path)?,
            annotation: ManifestArtifact::snapshot(
                "REFERENCE_GENE_ANNOTATION",
                bundle
                    .gene_annotation_path
                    .as_ref()
                    .context("reference annotation missing")?,
            )?,
        })
    }
    pub fn fingerprint(&self, scope: ReferenceScope) -> Result<String> {
        let files = match scope {
            ReferenceScope::Index => vec![&self.index],
            ReferenceScope::Analysis => vec![&self.index, &self.mapping, &self.annotation],
        };
        StageRequest {
            schema_version: 1,
            stage: "reference_binding".to_string(),
            algorithm_version: 1,
            engines: BTreeMap::new(),
            parameters: match scope {
                ReferenceScope::Index => "index",
                ReferenceScope::Analysis => "analysis",
            },
            inputs: files
                .into_iter()
                .map(|file| InputIdentity {
                    role: file.kind.clone(),
                    name: file.kind.clone(),
                    sha256: file.checksum.clone(),
                    bytes: file.bytes,
                })
                .collect(),
            dependencies: BTreeMap::new(),
        }
        .digest()
    }
    /// Choice identity includes provenance so an equivalent but older producer
    /// cannot replace a later explicit selection. Scientific cache keys use scopes.
    pub fn selection_token(&self) -> Result<ReferenceChoiceToken> {
        StageRequest {
            schema_version: 1,
            stage: "reference_choice".to_string(),
            algorithm_version: 1,
            engines: BTreeMap::new(),
            parameters: self,
            inputs: Vec::new(),
            dependencies: BTreeMap::new(),
        }
        .digest()
        .map(ReferenceChoiceToken)
    }
    pub fn verify_catalog(
        &self,
        db: &Database,
        config: &ProjectConfig,
        scope: ReferenceScope,
    ) -> Result<()> {
        let persisted: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
        if !self
            .policy
            .matches(&ReferencePolicy::from_config(config), scope)
            || !self
                .policy
                .matches(&ReferencePolicy::from_config(&persisted), scope)
        {
            bail!(
                "reference settings changed during computation; retry with the current reference"
            );
        }
        let selected = db
            .selected_reference()?
            .context("selected reference disappeared during computation")?;
        if self.fingerprint(scope)? != selected.fingerprint(scope)? {
            bail!(
                "selected reference changed during computation; retry with the current reference"
            );
        }
        Ok(())
    }
    pub fn verify_current(
        &self,
        db: &Database,
        config: &ProjectConfig,
        scope: ReferenceScope,
    ) -> Result<()> {
        self.verify_catalog(db, config, scope)?;
        let disk = ProjectConfig::load(&Path::new(&db.project()?.root_dir).join("rnaa.toml"))?;
        if !self
            .policy
            .matches(&ReferencePolicy::from_config(&disk), scope)
        {
            bail!("reference settings changed on disk during computation");
        }
        if !db.verify_generation_head(&self.generation)? {
            bail!("reference generation failed integrity verification");
        }
        for (path, file) in [
            (&disk.refs.custom_cdna, &self.cdna),
            (&disk.refs.custom_gtf, &self.gtf),
        ] {
            if !path.is_empty()
                && (scope == ReferenceScope::Analysis || file.kind == "REFERENCE_CDNA")
                && crate::util::compute_sha256(Path::new(path))? != file.checksum
            {
                bail!("reference source changed during computation: {path}");
            }
        }
        Ok(())
    }
}
fn decode(
    connection: &Connection,
    project: &str,
    current: bool,
) -> Result<Option<ReferenceSelection>> {
    let row:Option<(String,String,String)>=connection.query_row("SELECT generation_id,fingerprint,descriptor_json FROM reference_selection WHERE project_id=?1",[project],|row|Ok((row.get(0)?,row.get(1)?,row.get(2)?))).optional()?;
    row.map(|(generation, fingerprint, json)| {
        let selection: ReferenceSelection = serde_json::from_str(&json)?;
        if selection.schema_version != 1
            || selection.generation.generation_id != generation
            || selection.fingerprint(ReferenceScope::Analysis)? != fingerprint
        {
            bail!("selected reference descriptor is corrupt");
        }
        validate_catalog(connection, project, &selection, current)?;
        Ok(selection)
    })
    .transpose()
}
fn validate_catalog(
    connection: &Connection,
    project: &str,
    selection: &ReferenceSelection,
    current: bool,
) -> Result<()> {
    if selection.schema_version != 1
        || selection.generation.output_group != OutputGroup::Refs
        || selection.generation.subject != selection.reference_id
    {
        bail!("selected reference belongs to another output group or subject");
    }
    let row:Option<(String,String)>=connection.query_row("SELECT g.inventory_json,g.artifacts_json FROM artifact_generations g
        WHERE g.project_id=?1 AND g.output_group='refs' AND g.subject=?2 AND g.generation_id=?3 AND g.state='COMMITTED'
        AND (?4=0 OR EXISTS(SELECT 1 FROM generation_heads h WHERE h.project_id=g.project_id AND h.output_group=g.output_group AND h.subject=g.subject AND h.generation_id=g.generation_id))",
        params![project,selection.reference_id,selection.generation.generation_id,current],|row|Ok((row.get(0)?,row.get(1)?))).optional()?;
    let (seal, records) = row.context("selected reference has no matching committed head")?;
    let seal: serde_json::Value = serde_json::from_str(&seal)?;
    if seal["schema_version"] != 1 || seal["head"] != serde_json::to_value(&selection.generation)? {
        bail!("selected reference does not match its generation receipt");
    }
    let inventory: Vec<ManifestArtifact> = serde_json::from_value(seal["files"].clone())?;
    let records: Vec<crate::ArtifactRecord> = serde_json::from_str(&records)?;
    for (file, kind) in [
        (&selection.cdna, "REFERENCE_CDNA"),
        (&selection.gtf, "REFERENCE_GTF"),
        (&selection.index, "REFERENCE_INDEX"),
        (&selection.mapping, "REFERENCE_TX2GENE"),
        (&selection.annotation, "REFERENCE_GENE_ANNOTATION"),
    ] {
        if file.kind != kind
            || file.checksum_type != "sha256"
            || !Path::new(&file.path).starts_with(&selection.generation.payload_path)
            || !inventory.iter().any(|item| {
                item.path == file.path
                    && item.checksum_type == file.checksum_type
                    && item.checksum == file.checksum
                    && item.bytes == file.bytes
            })
            || !records.iter().any(|record| {
                record.kind.as_str() == file.kind
                    && record.path == file.path
                    && record.checksum == file.checksum
                    && record.bytes == file.bytes
                    && record.checksum_type == file.checksum_type
            })
        {
            bail!(
                "selected reference artifact does not match its committed registration: {}",
                file.path
            );
        }
    }
    Ok(())
}
impl Database {
    pub fn selected_reference(&self) -> Result<Option<ReferenceSelection>> {
        self.read(|connection| decode(connection, self.project_id(), true))
    }
    /// Called in reference publication's transaction, or as an atomic cache selection.
    /// The expected choice token prevents a slower producer from replacing a newer choice.
    pub fn select_reference(
        &self,
        selection: &ReferenceSelection,
        expected: Option<&ReferenceChoiceToken>,
    ) -> Result<bool> {
        self.write(|connection| {
            let retired:bool=connection.query_row("SELECT EXISTS(SELECT 1 FROM generation_retirements WHERE generation_id=?1)",[&selection.generation.generation_id],|row|row.get(0))?;
            if retired {bail!("restore the retired reference before selecting it");}
            validate_catalog(connection, self.project_id(), selection, true)?;
            let previous = decode(connection, self.project_id(), false)?;
            let fingerprint = selection.fingerprint(ReferenceScope::Analysis)?;
            let old = previous.as_ref().map(|value| value.fingerprint(ReferenceScope::Analysis)).transpose()?;
            let semantic_change = old.as_deref() != Some(fingerprint.as_str());
            let previous_token=previous.as_ref().map(ReferenceSelection::selection_token).transpose()?;
            if previous_token.as_ref() != expected && previous.as_ref() != Some(selection) {
                bail!("reference selection changed while preparing a bundle; retry");
            }
            let config: String = connection.query_row(
                "SELECT config_json FROM projects WHERE project_id = ?1",
                [self.project_id()], |row| row.get(0),
            )?;
            let config: ProjectConfig = serde_json::from_str(&config)?;
            config.refs.validate()?;
            if selection.policy != ReferencePolicy::from_config(&config) {
                bail!("reference settings changed while preparing a bundle; retry");
            }
            if previous.as_ref() == Some(selection) { return Ok(false); }
            connection.execute(
                "INSERT INTO reference_selection(project_id,generation_id,fingerprint,descriptor_json,updated_at)
                 VALUES(?1,?2,?3,?4,?5) ON CONFLICT(project_id) DO UPDATE SET
                 generation_id=excluded.generation_id,fingerprint=excluded.fingerprint,
                 descriptor_json=excluded.descriptor_json,updated_at=excluded.updated_at",
                params![self.project_id(), selection.generation.generation_id, fingerprint,
                    serde_json::to_string(selection)?, now_rfc3339()],
            )?;
            if semantic_change {
                connection.execute(
                    "UPDATE project_stages SET state='PENDING',last_error='selected reference changed',updated_at=?1
                     WHERE project_id=?2 AND stage IN ('normalize','deseq2','corr')",
                    params![now_rfc3339(), self.project_id()],
                )?;
            }
            let context = serde_json::json!({
                "reference_id":selection.reference_id,
                "generation_id":selection.generation.generation_id,
                "fingerprint":fingerprint,
                "semantic_change":semantic_change,
            });
            connection.execute(
                "INSERT INTO events(project_id,ts,stage,message,context_json)
                 VALUES(?1,?2,'refs','reference selected',?3)",
                params![self.project_id(), now_rfc3339(), context.to_string()],
            )?;
            Ok(true)
        })
    }
}

#[cfg(test)]
mod tests;
