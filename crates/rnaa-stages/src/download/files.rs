//! Private verified transfer caches and bounded receipt publication.
use super::*;

pub(super) fn verify_existing_or_download(
    remote: &RemoteFile,
    target: &Path,
    part: &Path,
    config: &rnaa_core::config::DownloadConfig,
    admission: Option<&rnaa_core::resources::ExecutionResources>,
    check: &impl Fn() -> Result<()>,
) -> Result<VerifiedFile> {
    let cap = resources::receipt_cap(remote, target)?;
    let verification = || resources::verify_request(remote, target);
    let guard = resources::reserve(admission, "verify", verification()?, check)?;
    let receipt_path = target.with_added_extension("rnaa-verified.json");
    if let Some(existing) = verify_local_file(remote, target)? {
        let receipt = resources::read_json::<serde_json::Value>(&receipt_path, cap).ok();
        if let Some(receipt) = receipt
            && receipt["source"] == serde_json::to_value(remote)?
            && receipt["artifact"]["checksum"] == existing.checksum
            && receipt["artifact"]["bytes"] == existing.bytes
            && receipt["http"]["prefix_sha256"] == existing.checksum
        {
            // A private cache can survive a failed multi-file attempt. Local
            // hashes alone do not establish that an unchecksummed URL is current.
            if remote.md5.is_some()
                || crate::http_transfer::validate_saved_version(
                    &remote.url,
                    &existing.checksum,
                    existing.bytes,
                    remote.bytes,
                    &receipt["http"],
                    config,
                    check,
                )?
                .is_some()
            {
                return Ok(existing);
            }
        }
    }
    drop(guard);
    let execution = crate::http_transfer::download_admitted(
        &remote.url,
        part,
        remote.md5.as_deref(),
        remote.bytes,
        config,
        admission,
        check,
    )?;
    let _guard = resources::reserve(admission, "verify", verification()?, check)?;
    let Some(mut verified) = verify_local_file(remote, part)? else {
        bail!("downloaded file failed verification");
    };
    fs::rename(part, target).with_context(|| {
        format!(
            "failed to publish private download cache {}",
            target.display()
        )
    })?;
    verified.path = target.to_path_buf();
    resources::write_json(
        &receipt_path,
        &json!({"source":remote,
            "artifact":ManifestArtifact::snapshot(verified.artifact_kind.as_str(), target)?,
            "http":crate::http_transfer::receipt(part)?,
            "resource_execution":execution,
        }),
        cap,
    )?;
    crate::http_transfer::remove_checkpoint(part)?;
    Ok(verified)
}

pub(super) fn verify_local_file(remote: &RemoteFile, path: &Path) -> Result<Option<VerifiedFile>> {
    if !path.exists() {
        return Ok(None);
    }
    let bytes = file_size(path)?;
    if let Some(expected_bytes) = remote.bytes
        && bytes != expected_bytes
    {
        return Ok(None);
    }

    let integrity_source = if let Some(expected_md5) = &remote.md5 {
        if !compute_md5(path)?.eq_ignore_ascii_case(expected_md5) {
            return Ok(None);
        }
        "authoritative_md5+local_sha256"
    } else if remote.bytes.is_some() {
        "remote_size+local_sha256"
    } else {
        "local_sha256_only"
    };
    let checksum = compute_sha256(path)?;
    Ok(Some(VerifiedFile {
        artifact_kind: match remote.kind {
            RemoteFileKind::FastqR1 => ArtifactKind::FastqR1,
            RemoteFileKind::FastqR2 => ArtifactKind::FastqR2,
            RemoteFileKind::Fastq => ArtifactKind::FastqSingle,
            RemoteFileKind::Sra | RemoteFileKind::Submitted => ArtifactKind::Sra,
        },
        path: path.to_path_buf(),
        checksum_type: "sha256".into(),
        checksum,
        bytes,
        integrity_source: integrity_source.into(),
    }))
}
