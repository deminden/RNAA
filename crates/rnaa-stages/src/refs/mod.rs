//! Verified reference preparation and immutable bundle publication.
mod annotation;
mod build;
mod observation;
mod planning;
mod reuse;
mod source;
use observation::ReferenceObservation;
use reuse::{cached_reference, reuse_reference};
#[cfg(test)]
mod tests;
use annotation::{build_annotation_tables, validate_index_mapping};
use anyhow::{Context, Result, anyhow, bail};
use reqwest::blocking::Client;
use rnaa_core::config::ProjectConfig;
use rnaa_core::fingerprint::{InputIdentity, StageRequest};
use rnaa_core::generations::{Generation, OutputGroup};
use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::ReferenceBundle;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::references::{ReferenceChoiceToken, ReferenceSelection};
use rnaa_core::traits::ReferenceManager;
use rnaa_core::util::{now_rfc3339, write_json_pretty};
use rnaa_core::{ArtifactKind, ArtifactRecord, Database};
use serde_json::json;
use source::{ReferenceSpec, acquire_source, resolve_reference_spec, verify_originals};
use std::collections::BTreeMap;
use std::fs::{self};
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct EnsemblReferenceManager {
    client: Client,
}

struct Preparation<'a> {
    spec: &'a ReferenceSpec,
    expected: Option<&'a ReferenceChoiceToken>,
    observation: ReferenceObservation,
}

impl Default for EnsemblReferenceManager {
    fn default() -> Self {
        Self::new().expect("failed to construct EnsemblReferenceManager")
    }
}

impl EnsemblReferenceManager {
    pub fn new() -> Result<Self> {
        Ok(Self {
            client: Client::builder()
                .timeout(std::time::Duration::from_secs(90))
                .build()
                .context("failed to construct HTTP client")?,
        })
    }
}

impl ReferenceManager for EnsemblReferenceManager {
    fn ensure_reference(
        &self,
        paths: &ProjectPaths,
        parent: &Database,
        config: &ProjectConfig,
    ) -> Result<ReferenceBundle> {
        config.refs.validate()?;
        config.download.validate()?;
        parent.check_job()?;
        let observation = ReferenceObservation::capture(parent, config)?;
        let expected = parent
            .selected_reference()?
            .map(|selection| selection.selection_token())
            .transpose()?;
        let spec = resolve_reference_spec(&self.client, config)?;
        let preparation = Preparation {
            spec: &spec,
            expected: expected.as_ref(),
            observation,
        };
        loop {
            parent.check_job()?;
            preparation.observation.verify_current(parent)?;
            // A no-change read needs no new attempt. Changing the project choice
            // or its observations uses the owned atomic reuse path.
            if let Some(bundle) = cached_reference(paths, parent, config, &spec)?
                && parent.selected_reference()?.is_some_and(|selected| {
                    selected.policy == rnaa_core::references::ReferencePolicy::from_config(config)
                        && bundle
                            .cdna_path
                            .starts_with(&selected.generation.payload_path)
                })
            {
                parent.check_job()?;
                select_bundle(parent, &bundle, config, &preparation)?;
                return Ok(bundle);
            }
            let Some(mut lease) = JobLease::try_acquire_reference(parent, &spec.reference_id)?
            else {
                std::thread::sleep(std::time::Duration::from_millis(100));
                continue;
            };
            let db = lease.database();
            // A competing producer may have committed between lookup and claim.
            let result = (|| {
                if let Some(bundle) = reuse_reference(paths, &db, parent, config, &preparation)? {
                    return Ok(bundle);
                }
                self.build_reference(paths, &db, parent, config, &preparation, lease.attempt_id())
            })();
            return match result {
                Ok(bundle) => {
                    lease.finish(JobOutcome::Succeeded, None)?;
                    Ok(bundle)
                }
                Err(error) => {
                    lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))?;
                    Err(error)
                }
            };
        }
    }
}

fn select_bundle(
    db: &Database,
    bundle: &ReferenceBundle,
    config: &ProjectConfig,
    preparation: &Preparation<'_>,
) -> Result<()> {
    let head = db
        .generation_head(OutputGroup::Refs, &bundle.id)?
        .context("reference head disappeared")?;
    if !bundle.cdna_path.starts_with(&head.payload_path) {
        bail!("reference head changed during selection; retry");
    }
    let selection = ReferenceSelection::capture(bundle, head, config)?;
    preparation.observation.verify_current(db)?;
    db.select_reference(&selection, preparation.expected)?;
    Ok(())
}

#[derive(Debug, Clone, serde::Serialize)]
struct IndexSettings {
    k: usize,
    g: usize,
    make_unique: bool,
    ec_max_size: i32,
    aa: bool,
}
impl IndexSettings {
    fn from_config(config: &ProjectConfig) -> Self {
        Self {
            k: config.refs.index_k,
            g: kallistors::index::builder::default_minimizer_len(config.refs.index_k),
            make_unique: false,
            ec_max_size: -1,
            aa: false,
        }
    }
    fn options(&self, threads: usize) -> kallistors::index::IndexBuildOptions {
        kallistors::index::IndexBuildOptions {
            k: self.k,
            g: Some(self.g),
            threads,
            make_unique: self.make_unique,
            ec_max_size: self.ec_max_size,
            aa: self.aa,
        }
    }
}
#[derive(serde::Serialize)]
struct ReferenceSettings<'a> {
    source_decode_contract: &'static str,
    provider: &'a str,
    organism: &'a str,
    release: &'a str,
    id_policy: &'a str,
    annotation_schema: u32,
    index: IndexSettings,
}
fn reference_request<'a>(
    spec: &'a ReferenceSpec,
    config: &ProjectConfig,
    inputs: &[ManifestArtifact],
) -> StageRequest<ReferenceSettings<'a>> {
    StageRequest {
        schema_version: 1,
        stage: "refs".into(),
        algorithm_version: 2,
        engines: BTreeMap::from([
            ("kallistors".into(), crate::engines::KALLISTORS.into()),
            (
                "dependency_lock".into(),
                crate::engines::DEPENDENCY_LOCK_SHA256.into(),
            ),
        ]),
        parameters: ReferenceSettings {
            source_decode_contract: "complete_multimember_fasta_gtf_v1",
            provider: &spec.provider,
            organism: &spec.organism,
            release: &spec.release,
            id_policy: "numeric_version_preserve_par_y_v1",
            annotation_schema: 1,
            index: IndexSettings::from_config(config),
        },
        inputs: inputs
            .iter()
            .map(|input| {
                let mut identity = InputIdentity::from(input);
                identity.name = input.kind.clone();
                identity
            })
            .collect(),
        dependencies: [
            ("transcriptome", &spec.transcript_source),
            ("annotations", &spec.gtf_source),
        ]
        .into_iter()
        .filter_map(|(role, source)| {
            source
                .md5
                .as_ref()
                .map(|md5| (format!("provider_md5:{role}"), md5.clone()))
        })
        .collect(),
    }
}
fn bundle(paths: &ProjectPaths, spec: &ReferenceSpec) -> ReferenceBundle {
    let dir = paths.reference_dir(&spec.reference_id);
    ReferenceBundle {
        id: spec.reference_id.clone(),
        organism: spec.organism.clone(),
        ensembl_release: spec.release.clone(),
        cdna_path: dir.join(if spec.transcript_source.gzip {
            "transcriptome.fa.gz"
        } else {
            "transcriptome.fa"
        }),
        gtf_path: dir.join(if spec.gtf_source.gzip {
            "annotations.gtf.gz"
        } else {
            "annotations.gtf"
        }),
        kallisto_index_path: dir.join("kallisto.index"),
        tx2gene_path: dir.join("tx2gene.tsv"),
        gene_annotation_path: Some(dir.join("gene_annotation.tsv")),
        manifest_path: paths.manifest_path("refs", &spec.reference_id),
    }
}
fn bundle_outputs(bundle: &ReferenceBundle) -> Vec<(ArtifactKind, PathBuf)> {
    vec![
        (ArtifactKind::ReferenceCdna, bundle.cdna_path.clone()),
        (ArtifactKind::ReferenceGtf, bundle.gtf_path.clone()),
        (
            ArtifactKind::ReferenceIndex,
            bundle.kallisto_index_path.clone(),
        ),
        (ArtifactKind::ReferenceTx2Gene, bundle.tx2gene_path.clone()),
        (
            ArtifactKind::ReferenceGeneAnnotation,
            bundle.gene_annotation_path.clone().unwrap(),
        ),
        (ArtifactKind::Manifest, bundle.manifest_path.clone()),
    ]
}
fn source_snapshots(bundle: &ReferenceBundle) -> Result<Vec<ManifestArtifact>> {
    Ok(vec![
        ManifestArtifact::snapshot("REFERENCE_CDNA", &bundle.cdna_path)?,
        ManifestArtifact::snapshot("REFERENCE_GTF", &bundle.gtf_path)?,
    ])
}
