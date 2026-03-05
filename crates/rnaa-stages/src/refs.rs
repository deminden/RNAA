use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result, anyhow, bail};
use flate2::read::MultiGzDecoder;
use regex::Regex;
use reqwest::blocking::Client;
use reqwest::header::USER_AGENT;
use rnaa_core::config::ProjectConfig;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::ReferenceBundle;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::ReferenceManager;
use rnaa_core::util::{
    command_version, compute_sha256, copy_if_needed, file_size, now_rfc3339, write_json_pretty,
};
use serde_json::json;

#[derive(Debug, Clone)]
pub struct EnsemblReferenceManager {
    client: Client,
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
        config: &ProjectConfig,
    ) -> Result<ReferenceBundle> {
        let started_at = now_rfc3339();
        let spec = resolve_reference_spec(&self.client, config)?;
        let reference_dir = paths.reference_dir(&spec.reference_id);
        fs::create_dir_all(&reference_dir)
            .with_context(|| format!("failed to create {}", reference_dir.display()))?;

        let cdna_path = reference_dir.join(
            spec.cdna_source
                .basename
                .clone()
                .unwrap_or_else(|| "transcriptome.cdna.fa.gz".to_string()),
        );
        let gtf_path = reference_dir.join(
            spec.gtf_source
                .basename
                .clone()
                .unwrap_or_else(|| "annotations.gtf.gz".to_string()),
        );
        ensure_source_materialized(&spec.cdna_source, &cdna_path)?;
        ensure_source_materialized(&spec.gtf_source, &gtf_path)?;

        let tx2gene_path = reference_dir.join("tx2gene.tsv");
        let gene_annotation_path = reference_dir.join("gene_annotation.tsv");
        if !tx2gene_path.exists() || !gene_annotation_path.exists() {
            build_annotation_tables(&gtf_path, &tx2gene_path, &gene_annotation_path)?;
        }

        let kallisto_index_path = reference_dir.join("kallisto.index");
        if !kallisto_index_path.exists() || file_size(&kallisto_index_path).unwrap_or_default() == 0
        {
            let status = Command::new("kallisto")
                .arg("index")
                .arg("-i")
                .arg(&kallisto_index_path)
                .arg(&cdna_path)
                .status()
                .context("failed to spawn kallisto")?;
            if !status.success() {
                bail!("kallisto index failed with status {status}");
            }
        }

        let manifest_path = paths.manifest_path("refs", &spec.reference_id);
        let manifest = StageManifest {
            schema_version: 1,
            stage: "refs".to_string(),
            id: spec.reference_id.clone(),
            started_at,
            finished_at: now_rfc3339(),
            parameters: json!({
                "organism": spec.organism,
                "ensembl": spec.release,
                "cdna_source": spec.cdna_source.describe(),
                "gtf_source": spec.gtf_source.describe(),
            }),
            tool_versions: BTreeMap::from_iter(
                [(
                    "kallisto".to_string(),
                    command_version("kallisto", &["version"]),
                )]
                .into_iter()
                .filter_map(|(name, version)| version.map(|version| (name, version))),
            ),
            input_artifacts: spec
                .inputs
                .iter()
                .map(|source| ManifestArtifact {
                    kind: source.kind.to_string(),
                    path: source.describe(),
                    checksum_type: "none".to_string(),
                    checksum: String::new(),
                    bytes: 0,
                })
                .collect(),
            output_artifacts: vec![
                manifest_artifact("REFERENCE_CDNA", &cdna_path)?,
                manifest_artifact("REFERENCE_GTF", &gtf_path)?,
                manifest_artifact("REFERENCE_TX2GENE", &tx2gene_path)?,
                manifest_artifact("REFERENCE_GENE_ANNOTATION", &gene_annotation_path)?,
                manifest_artifact("REFERENCE_INDEX", &kallisto_index_path)?,
            ],
            notes: vec![
                "Reference caching is deterministic by organism/release or custom source fingerprint."
                    .to_string(),
                "tx2gene.tsv and gene_annotation.tsv are derived from the GTF and reused by DE/correlation stages."
                    .to_string(),
            ],
        };
        write_json_pretty(&manifest_path, &manifest)?;

        Ok(ReferenceBundle {
            id: spec.reference_id,
            organism: spec.organism,
            ensembl_release: spec.release,
            cdna_path,
            gtf_path,
            kallisto_index_path,
            tx2gene_path,
            gene_annotation_path: Some(gene_annotation_path),
            manifest_path,
        })
    }
}

#[derive(Debug, Clone)]
struct ReferenceSpec {
    reference_id: String,
    organism: String,
    release: String,
    cdna_source: SourceSpec,
    gtf_source: SourceSpec,
    inputs: Vec<SourceSpec>,
}

#[derive(Debug, Clone)]
enum SourceKind {
    RemoteUrl,
    LocalFile,
}

impl std::fmt::Display for SourceKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RemoteUrl => f.write_str("remote_url"),
            Self::LocalFile => f.write_str("local_file"),
        }
    }
}

#[derive(Debug, Clone)]
struct SourceSpec {
    kind: SourceKind,
    location: String,
    basename: Option<String>,
}

impl SourceSpec {
    fn describe(&self) -> String {
        self.location.clone()
    }
}

#[derive(Debug, Clone, Copy)]
struct OrganismPreset {
    key: &'static str,
    species_dir: &'static str,
    prefix: &'static str,
}

const ORGANISM_PRESETS: &[OrganismPreset] = &[
    OrganismPreset {
        key: "human",
        species_dir: "homo_sapiens",
        prefix: "Homo_sapiens",
    },
    OrganismPreset {
        key: "mouse",
        species_dir: "mus_musculus",
        prefix: "Mus_musculus",
    },
];

fn resolve_reference_spec(client: &Client, config: &ProjectConfig) -> Result<ReferenceSpec> {
    if !config.refs.custom_cdna.is_empty() || !config.refs.custom_gtf.is_empty() {
        if config.refs.custom_cdna.is_empty() || config.refs.custom_gtf.is_empty() {
            bail!("custom references require both `custom_cdna` and `custom_gtf`");
        }
        let cdna = PathBuf::from(&config.refs.custom_cdna);
        let gtf = PathBuf::from(&config.refs.custom_gtf);
        if !cdna.exists() || !gtf.exists() {
            bail!("custom reference files do not exist");
        }
        let cdna_sha = compute_sha256(&cdna)?;
        let gtf_sha = compute_sha256(&gtf)?;
        let reference_id = format!(
            "custom_{}_{}",
            &cdna_sha[..12.min(cdna_sha.len())],
            &gtf_sha[..12.min(gtf_sha.len())]
        );
        let cdna_source = SourceSpec {
            kind: SourceKind::LocalFile,
            location: cdna.display().to_string(),
            basename: cdna
                .file_name()
                .map(|value| value.to_string_lossy().to_string()),
        };
        let gtf_source = SourceSpec {
            kind: SourceKind::LocalFile,
            location: gtf.display().to_string(),
            basename: gtf
                .file_name()
                .map(|value| value.to_string_lossy().to_string()),
        };
        return Ok(ReferenceSpec {
            reference_id,
            organism: config.refs.organism.clone(),
            release: "custom".to_string(),
            cdna_source: cdna_source.clone(),
            gtf_source: gtf_source.clone(),
            inputs: vec![cdna_source, gtf_source],
        });
    }

    let preset = ORGANISM_PRESETS
        .iter()
        .find(|preset| preset.key == config.refs.organism)
        .copied()
        .ok_or_else(|| anyhow!("unsupported organism preset {}", config.refs.organism))?;

    let latest = config.refs.ensembl == "latest";
    let cdna_dir = if latest {
        format!(
            "https://ftp.ensembl.org/pub/current_fasta/{}/cdna/",
            preset.species_dir
        )
    } else {
        format!(
            "https://ftp.ensembl.org/pub/release-{}/fasta/{}/cdna/",
            config.refs.ensembl, preset.species_dir
        )
    };
    let gtf_dir = if latest {
        format!(
            "https://ftp.ensembl.org/pub/current_gtf/{}/",
            preset.species_dir
        )
    } else {
        format!(
            "https://ftp.ensembl.org/pub/release-{}/gtf/{}/",
            config.refs.ensembl, preset.species_dir
        )
    };

    let cdna_listing = fetch_listing(client, &cdna_dir)?;
    let gtf_listing = fetch_listing(client, &gtf_dir)?;
    let cdna_file = find_cdna_file(&cdna_listing, preset.prefix)?;
    let gtf_file = find_gtf_file(&gtf_listing, preset.prefix)?;
    let release =
        extract_release_from_gtf(&gtf_file).unwrap_or_else(|| config.refs.ensembl.clone());
    let reference_id = format!("{}_ensembl_{}", preset.key, release);

    let cdna_source = SourceSpec {
        kind: SourceKind::RemoteUrl,
        location: format!("{cdna_dir}{cdna_file}"),
        basename: Some(cdna_file),
    };
    let gtf_source = SourceSpec {
        kind: SourceKind::RemoteUrl,
        location: format!("{gtf_dir}{gtf_file}"),
        basename: Some(gtf_file),
    };

    Ok(ReferenceSpec {
        reference_id,
        organism: preset.key.to_string(),
        release,
        cdna_source: cdna_source.clone(),
        gtf_source: gtf_source.clone(),
        inputs: vec![cdna_source, gtf_source],
    })
}

fn fetch_listing(client: &Client, url: &str) -> Result<String> {
    client
        .get(url)
        .header(USER_AGENT, "RNAA/0.1.0")
        .send()
        .and_then(reqwest::blocking::Response::error_for_status)
        .with_context(|| format!("failed to fetch {url}"))?
        .text()
        .with_context(|| format!("failed to read listing from {url}"))
}

fn find_cdna_file(listing: &str, prefix: &str) -> Result<String> {
    let regex = Regex::new(&format!(r#"href="({prefix}\.[^"]*cdna\.all\.fa\.gz)""#))?;
    regex
        .captures(listing)
        .and_then(|captures| captures.get(1))
        .map(|capture| capture.as_str().to_string())
        .ok_or_else(|| anyhow!("could not find cdna.all file in Ensembl listing"))
}

fn find_gtf_file(listing: &str, prefix: &str) -> Result<String> {
    let regex = Regex::new(&format!(r#"href="({prefix}\.[^"]*\.gtf\.gz)""#))?;
    let candidate = regex
        .captures_iter(listing)
        .filter_map(|captures| captures.get(1).map(|capture| capture.as_str().to_string()))
        .find(|name| {
            !name.contains(".chr.")
                && !name.contains("abinitio")
                && !name.contains("chr_patch_hapl_scaff")
        });
    candidate.ok_or_else(|| anyhow!("could not find canonical GTF file in Ensembl listing"))
}

fn extract_release_from_gtf(file_name: &str) -> Option<String> {
    Regex::new(r"\.(\d+)\.gtf\.gz$")
        .ok()
        .and_then(|regex| regex.captures(file_name))
        .and_then(|captures| captures.get(1).map(|value| value.as_str().to_string()))
}

fn ensure_source_materialized(source: &SourceSpec, dest: &Path) -> Result<()> {
    if dest.exists() && file_size(dest).unwrap_or_default() > 0 {
        return Ok(());
    }
    match source.kind {
        SourceKind::LocalFile => {
            copy_if_needed(Path::new(&source.location), dest)?;
        }
        SourceKind::RemoteUrl => {
            download_with_curl(&source.location, dest)?;
        }
    }
    Ok(())
}

fn download_with_curl(url: &str, dest: &Path) -> Result<()> {
    if let Some(parent) = dest.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let part = dest.with_extension(format!(
        "{}.part",
        dest.extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or_default()
    ));
    let status = Command::new("curl")
        .arg("--fail")
        .arg("--location")
        .arg("--continue-at")
        .arg("-")
        .arg("--retry")
        .arg("0")
        .arg("--output")
        .arg(&part)
        .arg(url)
        .status()
        .context("failed to spawn curl for reference download")?;
    if !status.success() {
        bail!("curl download failed for {url} with status {status}");
    }
    fs::rename(&part, dest)
        .with_context(|| format!("failed to rename {} to {}", part.display(), dest.display()))?;
    Ok(())
}

fn build_annotation_tables(
    gtf_path: &Path,
    tx2gene_path: &Path,
    gene_annotation_path: &Path,
) -> Result<()> {
    let file =
        File::open(gtf_path).with_context(|| format!("failed to open {}", gtf_path.display()))?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let transcript_re = Regex::new(r#"transcript_id "([^"]+)""#)?;
    let gene_re = Regex::new(r#"gene_id "([^"]+)""#)?;
    let gene_name_re = Regex::new(r#"gene_name "([^"]+)""#)?;
    let biotype_re = Regex::new(r#"(gene_biotype|gene_type) "([^"]+)""#)?;

    let mut tx2gene = BTreeSet::new();
    let mut gene_annotation: BTreeMap<String, (String, String)> = BTreeMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let _seqname = fields.next();
        let _source = fields.next();
        let feature = fields.next().unwrap_or_default();
        let _start = fields.next();
        let _end = fields.next();
        let _score = fields.next();
        let _strand = fields.next();
        let _frame = fields.next();
        let attributes = fields.next().unwrap_or_default();

        let Some(gene_id) = gene_re
            .captures(attributes)
            .and_then(|captures| captures.get(1).map(|value| value.as_str().to_string()))
        else {
            continue;
        };
        if let Some(transcript_id) = transcript_re
            .captures(attributes)
            .and_then(|captures| captures.get(1).map(|value| value.as_str().to_string()))
        {
            tx2gene.insert((transcript_id, gene_id.clone()));
        }

        if feature == "gene" || !gene_annotation.contains_key(&gene_id) {
            let gene_name = gene_name_re
                .captures(attributes)
                .and_then(|captures| captures.get(1).map(|value| value.as_str().to_string()))
                .unwrap_or_default();
            let biotype = biotype_re
                .captures(attributes)
                .and_then(|captures| captures.get(2).map(|value| value.as_str().to_string()))
                .unwrap_or_default();
            gene_annotation.insert(gene_id, (gene_name, biotype));
        }
    }

    let mut tx_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(tx2gene_path)
        .with_context(|| format!("failed to create {}", tx2gene_path.display()))?;
    tx_writer.write_record(["transcript_id", "gene_id"])?;
    for (transcript_id, gene_id) in tx2gene {
        tx_writer.write_record([transcript_id, gene_id])?;
    }
    tx_writer.flush()?;

    let mut gene_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(gene_annotation_path)
        .with_context(|| format!("failed to create {}", gene_annotation_path.display()))?;
    gene_writer.write_record(["gene_id", "gene_name", "gene_biotype"])?;
    for (gene_id, (gene_name, biotype)) in gene_annotation {
        gene_writer.write_record([gene_id, gene_name, biotype])?;
    }
    gene_writer.flush()?;

    Ok(())
}

fn manifest_artifact(kind: &str, path: &Path) -> Result<ManifestArtifact> {
    Ok(ManifestArtifact {
        kind: kind.to_string(),
        path: path.display().to_string(),
        checksum_type: "sha256".to_string(),
        checksum: compute_sha256(path)?,
        bytes: file_size(path)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_release_from_gtf_name() {
        assert_eq!(
            extract_release_from_gtf("Homo_sapiens.GRCh38.115.gtf.gz").as_deref(),
            Some("115")
        );
    }

    #[test]
    fn finds_cdna_and_gtf_from_listings() {
        let cdna_listing = r#"<a href="Homo_sapiens.GRCh38.cdna.all.fa.gz">"#;
        let gtf_listing = r#"<a href="Homo_sapiens.GRCh38.115.gtf.gz">"#;
        assert_eq!(
            find_cdna_file(cdna_listing, "Homo_sapiens").unwrap(),
            "Homo_sapiens.GRCh38.cdna.all.fa.gz"
        );
        assert_eq!(
            find_gtf_file(gtf_listing, "Homo_sapiens").unwrap(),
            "Homo_sapiens.GRCh38.115.gtf.gz"
        );
    }
}
