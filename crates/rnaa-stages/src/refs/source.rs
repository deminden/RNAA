//! Source resolution, provider checksums, durable transfers and source snapshots.
use anyhow::{Context, Result, anyhow, bail};
use regex::Regex;
use reqwest::blocking::Client;
use reqwest::header::USER_AGENT;
use rnaa_core::config::ProjectConfig;
use rnaa_core::manifest::ManifestArtifact;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::util::compute_sha256;
use serde_json::json;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
#[cfg(test)]
mod tests;

pub(super) fn verify_originals(spec: &ReferenceSpec) -> Result<()> {
    for source in [&spec.transcript_source, &spec.gtf_source] {
        if let Some(expected) = &source.sha256
            && compute_sha256(Path::new(&source.location))? != *expected
        {
            bail!(
                "reference source changed during preparation: {}",
                source.location
            );
        }
    }
    Ok(())
}
pub(super) fn transfer_path(
    paths: &ProjectPaths,
    source: &SourceSpec,
    subject: &str,
) -> Result<PathBuf> {
    use sha2::{Digest, Sha256};
    let key = hex::encode(Sha256::digest(serde_json::to_vec(&(
        &source.location,
        &source.md5,
        subject,
    ))?));
    Ok(paths.transfers_dir.join("refs").join(key).join("source"))
}
pub(super) struct MaterializedSource {
    pub(super) input: ManifestArtifact,
    pub(super) http: Option<serde_json::Value>,
    transfer: Option<PathBuf>,
    pub(super) execution: Option<crate::http_transfer::TransferExecution>,
}
impl MaterializedSource {
    pub(super) fn copy_to(&self, dest: &Path, check: &impl Fn() -> Result<()>) -> Result<()> {
        rnaa_core::util::copy_bounded(Path::new(&self.input.path), dest, self.input.bytes, check)?;
        if compute_sha256(dest)? != self.input.checksum || !self.input.matches_file()? {
            bail!(
                "reference source changed while copying: {}",
                self.input.path
            );
        }
        check()
    }
    pub(super) fn discard_transfer(&self) -> Result<()> {
        if let Some(path) = &self.transfer {
            fs::remove_file(path)?;
            crate::http_transfer::remove_checkpoint(path)?;
        }
        Ok(())
    }
}
pub(super) fn acquire_source(
    source: &SourceSpec,
    kind: rnaa_core::ArtifactKind,
    paths: &ProjectPaths,
    config: &ProjectConfig,
    subject: &str,
    admission: &rnaa_core::resources::ExecutionResources,
    check: &impl Fn() -> Result<()>,
) -> Result<MaterializedSource> {
    check()?;
    let mut verification = Some(admission.reserve(
        "source_verify",
        rnaa_core::resources::ResourceRequest {
            cpu_threads: 1,
            memory_bytes: 16 * 1024 * 1024,
            ..Default::default()
        },
        check,
    )?);
    let (origin, http, execution) = match source.kind {
        SourceKind::LocalFile => (PathBuf::from(&source.location), None, None),
        SourceKind::RemoteUrl => {
            let path = transfer_path(paths, source, subject)?;
            let verified = path.is_file()
                && source.md5.as_ref().is_some_and(|md5| {
                    rnaa_core::util::compute_md5(&path).is_ok_and(|actual| actual == *md5)
                });
            let execution = if !verified {
                drop(verification.take());
                let execution = crate::http_transfer::download_admitted(
                    &source.location,
                    &path,
                    source.md5.as_deref(),
                    None,
                    &config.download,
                    Some(admission),
                    check,
                )?;
                verification = Some(admission.reserve(
                    "source_verify",
                    rnaa_core::resources::ResourceRequest {
                        cpu_threads: 1,
                        memory_bytes: 16 * 1024 * 1024,
                        ..Default::default()
                    },
                    check,
                )?);
                execution
            } else {
                None
            };
            let receipt = crate::http_transfer::receipt(&path).ok();
            (path, receipt, execution)
        }
    };
    let _verification = verification;
    let input = ManifestArtifact::snapshot(kind.as_str(), &origin)?;
    if let Some(expected) = &source.sha256
        && input.checksum != *expected
    {
        bail!(
            "reference source changed before copying: {}",
            source.location
        );
    }
    check()?;
    let transfer = matches!(source.kind, SourceKind::RemoteUrl).then_some(origin);
    Ok(MaterializedSource {
        input,
        http,
        transfer,
        execution,
    })
}
#[derive(Debug, Clone)]
pub(super) struct ReferenceSpec {
    pub(super) reference_id: String,
    pub(super) provider: String,
    pub(super) organism: String,
    pub(super) release: String,
    pub(super) transcript_source: SourceSpec,
    pub(super) gtf_source: SourceSpec,
    pub(super) checksum_catalog: Option<serde_json::Value>,
}

#[derive(Debug, Clone, serde::Serialize)]
pub(super) enum SourceKind {
    RemoteUrl,
    LocalFile,
}

#[derive(Debug, Clone, serde::Serialize)]
pub(super) struct SourceSpec {
    pub(super) kind: SourceKind,
    pub(super) location: String,
    pub(super) gzip: bool,
    pub(super) sha256: Option<String>,
    pub(super) md5: Option<String>,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct OrganismPreset {
    pub(super) key: &'static str,
    pub(super) gencode_dir: &'static str,
    pub(super) gencode_prefix: &'static str,
}

pub(super) const ORGANISM_PRESETS: &[OrganismPreset] = &[
    OrganismPreset {
        key: "human",
        gencode_dir: "Gencode_human",
        gencode_prefix: "gencode.v",
    },
    OrganismPreset {
        key: "mouse",
        gencode_dir: "Gencode_mouse",
        gencode_prefix: "gencode.vM",
    },
];

pub(super) const PRESET_REFERENCE_LAYOUT_VERSION: &str = "txall_v2";

pub(super) fn resolve_reference_spec(
    client: &Client,
    config: &ProjectConfig,
) -> Result<ReferenceSpec> {
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
            gzip: is_gzip(&cdna)?,
            sha256: Some(cdna_sha),
            md5: None,
        };
        let gtf_source = SourceSpec {
            kind: SourceKind::LocalFile,
            location: gtf.display().to_string(),
            gzip: is_gzip(&gtf)?,
            sha256: Some(gtf_sha),
            md5: None,
        };
        return Ok(ReferenceSpec {
            reference_id,
            provider: "custom".to_string(),
            organism: config.refs.organism.clone(),
            release: "custom".to_string(),
            transcript_source: cdna_source.clone(),
            gtf_source: gtf_source.clone(),
            checksum_catalog: None,
        });
    }

    let preset = ORGANISM_PRESETS
        .iter()
        .find(|preset| preset.key == config.refs.organism)
        .copied()
        .ok_or_else(|| anyhow!("unsupported organism preset {}", config.refs.organism))?;

    let requested_release = normalize_gencode_release(&config.refs.ensembl, preset)?;
    let folder = requested_release
        .as_ref()
        .map(|release| format!("release_{}", release.trim_start_matches('v')))
        .unwrap_or_else(|| "latest_release".into());
    let base = format!(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/{}/{folder}/",
        preset.gencode_dir
    );
    let listing = fetch_listing(client, &base)?;
    let (transcript_file, gtf_file, release) = matching_release_files(
        &listing,
        preset.gencode_prefix,
        requested_release.as_deref(),
    )?;
    // Resolve latest to a concrete release before fetching either source/checksum.
    let base = format!(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/{}/release_{}/",
        preset.gencode_dir,
        release.trim_start_matches('v')
    );
    let catalog_url = format!("{base}MD5SUMS");
    let catalog = fetch_listing(client, &catalog_url)?;
    let transcript_source = SourceSpec {
        kind: SourceKind::RemoteUrl,
        location: format!("{base}{transcript_file}"),
        gzip: true,
        sha256: None,
        md5: Some(provider_md5(&catalog, &transcript_file)?),
    };
    let gtf_source = SourceSpec {
        kind: SourceKind::RemoteUrl,
        location: format!("{base}{gtf_file}"),
        gzip: true,
        sha256: None,
        md5: Some(provider_md5(&catalog, &gtf_file)?),
    };
    Ok(ReferenceSpec {
        reference_id: format!(
            "{}_gencode_{}_{}",
            preset.key,
            sanitize_release_id(&release),
            PRESET_REFERENCE_LAYOUT_VERSION
        ),
        provider: "gencode".into(),
        organism: preset.key.into(),
        release,
        transcript_source,
        gtf_source,
        checksum_catalog: Some(json!({"url":catalog_url,"contents":catalog})),
    })
}
pub(super) fn matching_release_files(
    listing: &str,
    prefix: &str,
    requested: Option<&str>,
) -> Result<(String, String, String)> {
    let regex = Regex::new(&format!(r#"href="(({prefix}\d+)\.transcripts\.fa\.gz)""#))?;
    let mut files = regex
        .captures_iter(listing)
        .filter_map(|capture| capture.get(1).map(|value| value.as_str().to_string()))
        .filter_map(|name| extract_gencode_release(&name).map(|release| (name, release)))
        .filter(|(_, release)| {
            requested.is_none_or(|value| {
                value.trim_start_matches('v') == release.trim_start_matches('v')
            })
        })
        .collect::<Vec<_>>();
    files.sort_by_key(|(_, release)| {
        release
            .trim_start_matches('v')
            .trim_start_matches('M')
            .parse::<u32>()
            .unwrap_or(0)
    });
    let (transcript, release) = files
        .pop()
        .context("reference listing does not contain the requested transcript release")?;
    let gtf = format!(
        "{prefix}{}.annotation.gtf.gz",
        release.trim_start_matches('v').trim_start_matches('M')
    );
    if !listing.contains(&format!("href=\"{gtf}\"")) {
        bail!("reference FASTA and GTF do not belong to the same release");
    }
    Ok((transcript, gtf, release))
}
pub(super) fn provider_md5(catalog: &str, name: &str) -> Result<String> {
    let mut result = None;
    for line in catalog.lines() {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() == 2 && fields[1].trim_start_matches('*').trim_start_matches("./") == name {
            if result.is_some()
                || fields[0].len() != 32
                || !fields[0].bytes().all(|value| value.is_ascii_hexdigit())
            {
                bail!("invalid or duplicate provider checksum for {name}");
            }
            result = Some(fields[0].to_ascii_lowercase());
        }
    }
    result.with_context(|| format!("provider checksum missing for {name}"))
}
pub(super) fn is_gzip(path: &Path) -> Result<bool> {
    use std::io::Read;
    let mut input = File::open(path)?;
    let mut magic = [0; 2];
    Ok(input.read(&mut magic)? == 2 && magic == [0x1f, 0x8b])
}

pub(super) fn fetch_listing(client: &Client, url: &str) -> Result<String> {
    use std::io::Read;
    let response = client
        .get(url)
        .header(USER_AGENT, "RNAA/0.1.0")
        .send()
        .and_then(reqwest::blocking::Response::error_for_status)
        .with_context(|| format!("failed to fetch {url}"))?;
    let mut text = String::new();
    response.take(1024 * 1024 + 1).read_to_string(&mut text)?;
    if text.len() > 1024 * 1024 {
        bail!("reference provider listing/checksum catalog exceeds 1 MiB");
    }
    Ok(text)
}

pub(super) fn extract_gencode_release(file_name: &str) -> Option<String> {
    Regex::new(r"(vM?\d+)")
        .ok()
        .and_then(|regex| regex.captures(file_name))
        .and_then(|captures| captures.get(1).map(|value| value.as_str().to_string()))
}

pub(super) fn normalize_gencode_release(
    release: &str,
    preset: OrganismPreset,
) -> Result<Option<String>> {
    let release = release.trim();
    if release.eq_ignore_ascii_case("latest") {
        return Ok(None);
    }
    let number = release.trim_start_matches('v');
    let number = if preset.key == "mouse" {
        number.trim_start_matches('M')
    } else {
        number
    };
    if number.is_empty()
        || !number.bytes().all(|value| value.is_ascii_digit())
        || number.parse::<u32>().is_err()
        || number == "0"
    {
        bail!("invalid GENCODE reference release {release}");
    }
    Ok(Some(if preset.key == "mouse" {
        format!("M{number}")
    } else {
        format!("v{number}")
    }))
}
pub(super) fn sanitize_release_id(release: &str) -> String {
    release.replace('.', "_")
}
