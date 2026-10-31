//! Bounded, validator-aware HTTP transfers. Resume metadata is a durable prefix
//! checkpoint; bytes beyond it are never trusted after interruption.
mod admission;
mod version;
pub(crate) use admission::TransferExecution;
use admission::{ReadmitLength, TransferBudget};
use anyhow::{Context, Result, bail};
use reqwest::{
    StatusCode,
    blocking::Client,
    header::{
        ACCEPT_ENCODING, CONTENT_ENCODING, CONTENT_RANGE, ETAG, IF_RANGE, LAST_MODIFIED, RANGE,
    },
};
use rnaa_core::{
    config::DownloadConfig,
    util::{compute_md5, write_json_pretty},
};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    time::{Duration, Instant},
};
pub(crate) use version::validate_saved_version;

const CHECKPOINT_BYTES: u64 = 8 * 1024 * 1024;
const CHECKPOINT_METADATA_LIMIT: u64 = 1024 * 1024;

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Checkpoint {
    schema_version: u32,
    url: String,
    expected_md5: Option<String>,
    expected_bytes: Option<u64>,
    effective_url: String,
    etag: Option<String>,
    last_modified: Option<String>,
    total_bytes: Option<u64>,
    committed_bytes: u64,
    prefix_sha256: String,
}

fn sidecar(part: &Path) -> PathBuf {
    part.with_added_extension("rnaa-resume.json")
}

fn digest(hash: &Sha256) -> String {
    hex::encode(hash.clone().finalize())
}

#[cfg(test)]
fn hash_prefix(path: &Path, bytes: u64) -> Result<Sha256> {
    hash_prefix_checked(path, bytes, &|| Ok(()))
}

fn hash_prefix_checked(path: &Path, bytes: u64, check: &impl Fn() -> Result<()>) -> Result<Sha256> {
    let mut input = File::open(path)?.take(bytes);
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 65536];
    let mut total = 0_u64;
    let mut next_check = 0_u64;
    loop {
        if total >= next_check {
            check()?;
            next_check = total.saturating_add(CHECKPOINT_BYTES);
        }
        let count = input.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
        total += count as u64;
    }
    if total != bytes {
        bail!("resume prefix is truncated");
    }
    Ok(hasher)
}

fn load_prefix(
    url: &str,
    part: &Path,
    md5: Option<&str>,
    bytes: Option<u64>,
    check: &impl Fn() -> Result<()>,
) -> Result<(Option<Checkpoint>, Sha256)> {
    let saved = checkpoint_bytes(part)
        .ok()
        .and_then(|raw| serde_json::from_slice::<Checkpoint>(&raw).ok());
    if let Some(saved) = saved
        && saved.schema_version == 1
        && saved.url == url
        && saved.expected_md5.as_deref() == md5
        && saved.expected_bytes == bytes
        && (md5.is_some() || saved.validator().is_some())
        && let Ok(hash) = hash_prefix_checked(part, saved.committed_bytes, check)
        && digest(&hash) == saved.prefix_sha256
        && bytes.is_none_or(|total| saved.committed_bytes <= total)
        && saved
            .total_bytes
            .is_none_or(|total| saved.committed_bytes <= total)
    {
        check()?;
        // A process may die after flushing data but before publishing its next
        // checkpoint. Discard that uncommitted tail before asking for a range.
        OpenOptions::new()
            .write(true)
            .open(part)?
            .set_len(saved.committed_bytes)?;
        return Ok((Some(saved), hash));
    }
    // A cancellation during prefix hashing must not become a destructive reset.
    check()?;
    if part.exists() {
        let file = OpenOptions::new().write(true).open(part)?;
        file.set_len(0)?;
        file.sync_all()?;
    }
    Ok((None, Sha256::new()))
}

impl Checkpoint {
    fn validator(&self) -> Option<&str> {
        self.etag
            .as_deref()
            .filter(|etag| etag.starts_with('"') && etag.ends_with('"'))
            .or_else(|| {
                self.expected_md5
                    .as_ref()
                    .and(self.last_modified.as_deref())
            })
    }
}

#[derive(Debug, PartialEq, Eq)]
struct ContentRange {
    start: u64,
    end: u64,
    total: u64,
}

fn content_range(raw: &str) -> Result<ContentRange> {
    let (span, total) = raw
        .strip_prefix("bytes ")
        .context("invalid Content-Range unit")?
        .split_once('/')
        .context("invalid Content-Range")?;
    let (start, end) = span.split_once('-').context("invalid Content-Range span")?;
    let (start, end, total): (u64, u64, u64) = (start.parse()?, end.parse()?, total.parse()?);
    if start > end || end >= total {
        bail!("invalid Content-Range bounds");
    }
    Ok(ContentRange { start, end, total })
}

fn checkpoint(
    file: &File,
    path: &Path,
    saved: &mut Checkpoint,
    hash: &Sha256,
    bytes: u64,
) -> Result<()> {
    file.sync_all()?;
    saved.committed_bytes = bytes;
    saved.prefix_sha256 = digest(hash);
    if serde_json::to_vec_pretty(saved)?.len() as u64 > CHECKPOINT_METADATA_LIMIT {
        bail!("HTTP checkpoint metadata exceeds its bounded allowance");
    }
    write_json_pretty(&sidecar(path), saved)
}

fn checkpoint_bytes(part: &Path) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    File::open(sidecar(part))?
        .take(CHECKPOINT_METADATA_LIMIT + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 > CHECKPOINT_METADATA_LIMIT {
        bail!("HTTP checkpoint metadata exceeds its bounded allowance");
    }
    Ok(bytes)
}

pub(crate) fn remove_checkpoint(part: &Path) -> Result<()> {
    match fs::remove_file(sidecar(part)) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error.into()),
    }
}

pub(crate) fn receipt(part: &Path) -> Result<serde_json::Value> {
    let checkpoint: Checkpoint = serde_json::from_slice(&checkpoint_bytes(part)?)?;
    Ok(serde_json::to_value(checkpoint)?)
}

#[cfg(test)]
pub(crate) fn download(
    url: &str,
    part: &Path,
    md5: Option<&str>,
    expected_bytes: Option<u64>,
    config: &DownloadConfig,
    check: &impl Fn() -> Result<()>,
) -> Result<()> {
    download_admitted(url, part, md5, expected_bytes, config, None, check).map(|_| ())
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn download_admitted(
    url: &str,
    part: &Path,
    md5: Option<&str>,
    expected_bytes: Option<u64>,
    config: &DownloadConfig,
    admission: Option<&rnaa_core::resources::ExecutionResources>,
    check: &impl Fn() -> Result<()>,
) -> Result<Option<TransferExecution>> {
    config.validate()?;
    let deadline = Instant::now()
        .checked_add(Duration::from_secs(config.max_transfer_seconds))
        .context("transfer deadline is too large")?;
    let mut attempt = 0_u32;
    let mut learned_total = None;
    let fence = || {
        check()?;
        if Instant::now() >= deadline {
            bail!("HTTP transfer exceeded max_transfer_seconds");
        }
        Ok(())
    };
    loop {
        fence()?;
        let scan = admission
            .map(|c| -> Result<_> {
                c.reserve("transfer_scan", admission::request(part, 0, true)?, &fence)
            })
            .transpose()?;
        if let Some(parent) = part.parent() {
            fs::create_dir_all(parent)?;
        }
        let prefix = load_prefix(url, part, md5, expected_bytes, &fence)?;
        let start = prefix.0.as_ref().map_or(0, |p| p.committed_bytes);
        let total = expected_bytes
            .or(learned_total)
            .or_else(|| prefix.0.as_ref().and_then(|p| p.total_bytes));
        drop(scan);
        let mut budget = TransferBudget::reserve(admission, part, start, total, &fence)?;
        match request(
            url,
            part,
            md5,
            expected_bytes,
            config,
            deadline,
            &fence,
            prefix,
            &mut budget,
        ) {
            Ok(true) => return Ok(budget.execution()),
            Ok(false) => continue, // A valid shorter range made progress.
            Err(error) => {
                // The request and its response buffers have unwound. Release
                // all phase claims before readmission or network backoff.
                drop(budget);
                fence()?;
                if let Some(ReadmitLength(total)) = error.downcast_ref::<ReadmitLength>() {
                    learned_total = Some(*total);
                    continue;
                }
                if error.is::<rnaa_core::resources::ResourceAdmissionError>() {
                    return Err(error);
                }
                if attempt >= config.retries {
                    return Err(error).context("HTTP transfer exhausted retries");
                }
                let seconds = config
                    .backoff_seconds
                    .get(attempt as usize)
                    .copied()
                    .unwrap_or(30);
                attempt += 1;
                tracing::warn!(url, attempt, "HTTP transfer retry: {error:#}");
                let until = Instant::now()
                    .checked_add(Duration::from_secs(seconds))
                    .context("retry backoff is too large")?
                    .min(deadline);
                while Instant::now() < until {
                    check()?;
                    std::thread::sleep(
                        Duration::from_millis(100)
                            .min(until.saturating_duration_since(Instant::now())),
                    );
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn request(
    url: &str,
    part: &Path,
    md5: Option<&str>,
    expected_bytes: Option<u64>,
    config: &DownloadConfig,
    deadline: Instant,
    check: &impl Fn() -> Result<()>,
    prefix: (Option<Checkpoint>, Sha256),
    budget: &mut TransferBudget<'_>,
) -> Result<bool> {
    let (previous, mut hash) = prefix;
    let offset = previous.as_ref().map_or(0, |saved| saved.committed_bytes);
    let remaining = deadline.saturating_duration_since(Instant::now());
    if remaining.is_zero() {
        bail!("HTTP transfer deadline reached");
    }
    let client = Client::builder()
        .connect_timeout(Duration::from_secs(config.connect_timeout_seconds).min(remaining))
        .timeout(Duration::from_secs(config.read_timeout_seconds).min(remaining))
        .no_gzip()
        .no_brotli()
        .no_deflate()
        .no_zstd()
        .build()?;
    let mut request = client.get(url).header(ACCEPT_ENCODING, "identity");
    if offset > 0 {
        request = request.header(RANGE, format!("bytes={offset}-"));
        if let Some(validator) = previous.as_ref().and_then(Checkpoint::validator) {
            request = request.header(IF_RANGE, validator);
        }
    }
    check()?;
    let mut response = request
        .send()
        .with_context(|| format!("failed to request {url}"))?;
    if response
        .headers()
        .get(CONTENT_ENCODING)
        .is_some_and(|value| value != "identity")
    {
        bail!("HTTP download requires an unencoded representation for byte verification");
    }
    let etag = response
        .headers()
        .get(ETAG)
        .map(|value| value.to_str().map(str::to_owned))
        .transpose()?;
    let modified = response
        .headers()
        .get(LAST_MODIFIED)
        .map(|value| value.to_str().map(str::to_owned))
        .transpose()?;
    if response.status() == StatusCode::RANGE_NOT_SATISFIABLE && offset > 0 {
        let total: Option<u64> = response
            .headers()
            .get(CONTENT_RANGE)
            .and_then(|value| value.to_str().ok())
            .and_then(|value| value.strip_prefix("bytes */"))
            .and_then(|value| value.parse().ok());
        let unchanged = previous.as_ref().is_some_and(|saved| {
            saved.etag.is_some() && saved.etag == etag
                || saved.last_modified.is_some() && saved.last_modified == modified
        });
        if total == Some(offset)
            && expected_bytes.is_none_or(|expected| expected == offset)
            && previous
                .as_ref()
                .is_some_and(|saved| saved.total_bytes == total)
            && (unchanged || md5.is_some())
        {
            if let Some(md5) = md5
                && !compute_md5(part)?.eq_ignore_ascii_case(md5)
            {
                reset(part)?;
                bail!("completed range has the wrong MD5");
            }
            budget.finish(offset)?;
            return Ok(true);
        }
        reset(part)?;
        bail!("HTTP 416 did not prove that the partial file is complete");
    }
    if !matches!(
        response.status(),
        StatusCode::OK | StatusCode::PARTIAL_CONTENT
    ) {
        bail!(
            "HTTP download failed for {url} with status {}",
            response.status()
        );
    }
    let partial = response.status() == StatusCode::PARTIAL_CONTENT;
    let (start, total, expected_body) = if partial {
        let raw = response
            .headers()
            .get(CONTENT_RANGE)
            .context("206 response lacks Content-Range")?
            .to_str()?;
        let range = content_range(raw)?;
        if range.start != offset {
            bail!("HTTP range starts at {}, expected {offset}", range.start);
        }
        if let Some(saved) = &previous
            && (saved.total_bytes.is_some_and(|total| total != range.total)
                || saved.etag.is_some() && saved.etag != etag
                || saved.etag.is_none()
                    && saved.last_modified.is_some()
                    && saved.last_modified != modified
                || saved.effective_url != response.url().as_str())
        {
            reset(part)?;
            bail!("HTTP entity changed during resume");
        }
        let length = range.end - range.start + 1;
        if response
            .content_length()
            .is_some_and(|value| value != length)
        {
            bail!("Content-Length disagrees with Content-Range");
        }
        (range.start, Some(range.total), Some(length))
    } else {
        (0, response.content_length(), response.content_length())
    };
    if expected_bytes.is_some() && total.is_some() && expected_bytes != total {
        bail!("remote size disagrees with the expected file size");
    }
    if start == 0 {
        hash = Sha256::new();
        // A full response replaces the private prefix. Reclaim its actual
        // storage before reserving the replacement's complete future writes.
        if offset > 0 {
            reset(part)?;
        }
    }
    budget.begin(start, total.or(expected_bytes))?;
    let mut file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(start == 0)
        .append(start > 0)
        .open(part)?;
    let mut saved = Checkpoint {
        schema_version: 1,
        url: url.into(),
        expected_md5: md5.map(str::to_owned),
        expected_bytes,
        effective_url: response.url().to_string(),
        etag,
        last_modified: modified,
        total_bytes: total.or(expected_bytes),
        committed_bytes: start,
        prefix_sha256: digest(&hash),
    };
    checkpoint(&file, part, &mut saved, &hash, start)?;
    let mut buffer = [0_u8; 65536];
    let mut received = 0_u64;
    let mut written = start;
    let mut last_checkpoint = Instant::now();
    loop {
        check()?;
        if Instant::now() >= deadline {
            checkpoint(&file, part, &mut saved, &hash, written)?;
            bail!("HTTP transfer deadline reached");
        }
        let count = match response.read(&mut buffer) {
            Ok(count) => count,
            Err(error) => {
                checkpoint(&file, part, &mut saved, &hash, written)?;
                return Err(error).context("HTTP body was interrupted or stalled");
            }
        };
        if count == 0 {
            break;
        }
        received = received
            .checked_add(count as u64)
            .context("HTTP byte count overflow")?;
        let next_written = written
            .checked_add(count as u64)
            .context("HTTP byte count overflow")?;
        if expected_body.is_some_and(|total| received > total)
            || saved.total_bytes.is_some_and(|total| next_written > total)
        {
            reset(part)?;
            bail!("HTTP body exceeds its declared bounds");
        }
        if let Err(error) = budget.before_write(written, next_written) {
            checkpoint(&file, part, &mut saved, &hash, written)?;
            return Err(error);
        }
        file.write_all(&buffer[..count])?;
        hash.update(&buffer[..count]);
        written = next_written;
        if written - saved.committed_bytes >= CHECKPOINT_BYTES
            || last_checkpoint.elapsed() >= Duration::from_secs(1)
        {
            checkpoint(&file, part, &mut saved, &hash, written)?;
            budget.progress(written)?;
            last_checkpoint = Instant::now();
        }
    }
    checkpoint(&file, part, &mut saved, &hash, written)?;
    if expected_body.is_some_and(|total| received != total) {
        bail!("HTTP body was truncated");
    }
    if let Some(total) = saved.total_bytes {
        if written < total {
            return Ok(false);
        }
        if written != total {
            reset(part)?;
            bail!("HTTP size verification failed");
        }
    }
    if let Some(expected) = md5
        && !compute_md5(part)?.eq_ignore_ascii_case(expected)
    {
        reset(part)?;
        bail!("downloaded MD5 differs from the expected digest");
    }
    budget.finish(written)?;
    Ok(true)
}

fn reset(part: &Path) -> Result<()> {
    if part.exists() {
        File::create(part)?.sync_all()?;
    }
    remove_checkpoint(part)
}

#[cfg(test)]
mod tests;
