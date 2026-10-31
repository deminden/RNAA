use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use chrono::Utc;
use md5::{Digest as Md5Digest, Md5};
use serde::Serialize;
use sha2::Sha256;
use walkdir::WalkDir;

pub fn now_rfc3339() -> String {
    Utc::now().to_rfc3339()
}

pub fn write_json_pretty<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let content = serde_json::to_string_pretty(value).context("failed to serialize json")?;
    publish_file(path, |staged| {
        std::fs::write(staged, content)
            .with_context(|| format!("failed to write {}", path.display()))
    })
}

pub fn compute_md5(path: &Path) -> Result<String> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Md5::new();
    let mut buffer = [0_u8; 1024 * 1024];
    loop {
        let read = reader
            .read(&mut buffer)
            .with_context(|| format!("failed to read {}", path.display()))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(hex::encode(hasher.finalize()))
}

pub fn compute_sha256(path: &Path) -> Result<String> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 1024 * 1024];
    loop {
        let read = reader
            .read(&mut buffer)
            .with_context(|| format!("failed to read {}", path.display()))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(hex::encode(hasher.finalize()))
}

pub fn file_size(path: &Path) -> Result<u64> {
    Ok(std::fs::metadata(path)
        .with_context(|| format!("failed to stat {}", path.display()))?
        .len())
}

pub fn dir_size(path: &Path) -> u64 {
    WalkDir::new(path)
        .into_iter()
        .filter_map(Result::ok)
        .filter_map(|entry| entry.metadata().ok())
        .filter(|meta| meta.is_file())
        .map(|meta| meta.len())
        .sum()
}

pub fn parse_size_bytes(raw: &str) -> Result<u64> {
    let value = raw.trim();
    if value.is_empty() {
        anyhow::bail!("size string is empty");
    }

    let split_idx = value
        .find(|ch: char| !ch.is_ascii_digit() && ch != '.')
        .unwrap_or(value.len());
    let number = value[..split_idx].trim();
    let suffix = value[split_idx..].trim().to_ascii_lowercase();
    let scalar = number
        .parse::<f64>()
        .with_context(|| format!("invalid size value '{}'", value))?;
    if !scalar.is_finite() || scalar < 0.0 {
        anyhow::bail!("invalid size value '{}'", value);
    }

    let multiplier = match suffix.as_str() {
        "" | "b" => 1_f64,
        "k" | "kb" => 1_000_f64,
        "m" | "mb" => 1_000_000_f64,
        "g" | "gb" => 1_000_000_000_f64,
        "t" | "tb" => 1_000_000_000_000_f64,
        "kib" => 1024_f64,
        "mib" => 1024_f64 * 1024_f64,
        "gib" => 1024_f64 * 1024_f64 * 1024_f64,
        "tib" => 1024_f64 * 1024_f64 * 1024_f64 * 1024_f64,
        _ => anyhow::bail!("unsupported size suffix '{}'", suffix),
    };
    let bytes = scalar * multiplier;
    if bytes > u64::MAX as f64 {
        anyhow::bail!("size '{}' exceeds supported range", value);
    }
    Ok(bytes.round() as u64)
}

pub fn sanitize_basename(raw: &str) -> String {
    raw.chars()
        .map(|ch| match ch {
            'A'..='Z' | 'a'..='z' | '0'..='9' | '.' | '_' | '-' => ch,
            _ => '_',
        })
        .collect()
}

pub fn relative_or_absolute(path: &Path, root: &Path) -> String {
    path.strip_prefix(root)
        .map(|value| value.display().to_string())
        .unwrap_or_else(|_| path.display().to_string())
}

pub fn normalize_url(url: &str) -> String {
    if url.starts_with("http://") || url.starts_with("https://") {
        url.to_string()
    } else if url.starts_with("ftp://") {
        url.replacen("ftp://", "https://", 1)
    } else {
        format!("https://{}", url.trim_start_matches('/'))
    }
}

pub fn map_to_json_strings(map: &BTreeMap<String, String>) -> serde_json::Value {
    let object = map
        .iter()
        .map(|(key, value)| (key.clone(), serde_json::Value::String(value.clone())))
        .collect();
    serde_json::Value::Object(object)
}

pub fn copy_if_needed(src: &Path, dest: &Path) -> Result<PathBuf> {
    if src == dest {
        return Ok(dest.to_path_buf());
    }
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    std::fs::copy(src, dest)
        .with_context(|| format!("failed to copy {} to {}", src.display(), dest.display()))?;
    Ok(dest.to_path_buf())
}

pub fn stable_map_from_pairs(
    pairs: impl IntoIterator<Item = (String, String)>,
) -> BTreeMap<String, String> {
    pairs.into_iter().collect()
}

/// Normalize provider FASTA identifiers while preserving meaningful suffixes
/// such as GENCODE's _PAR_Y. Only numeric version components are removed.
pub fn normalize_feature_id(raw: &str) -> String {
    let id = raw.trim().split('|').next().unwrap_or_default();
    let (base, suffix) = id
        .strip_suffix("_PAR_Y")
        .map(|base| (base, "_PAR_Y"))
        .unwrap_or((id, ""));
    let normalized = match base.rsplit_once('.') {
        Some((prefix, version))
            if !version.is_empty() && version.bytes().all(|b| b.is_ascii_digit()) =>
        {
            prefix
        }
        _ => base,
    };
    format!("{normalized}{suffix}")
}

/// Write a replacement on a fresh inode and publish it only after success.
/// Existing hardlinks continue to reference the previous immutable contents.
pub fn publish_file(path: &Path, writer: impl FnOnce(&Path) -> Result<()>) -> Result<()> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let staged = tempfile::NamedTempFile::new_in(parent)
        .with_context(|| format!("failed to stage {}", path.display()))?;
    writer(staged.path())?;
    staged.as_file().sync_all()?;
    staged
        .persist(path)
        .with_context(|| format!("failed to publish {}", path.display()))?;
    #[cfg(unix)]
    File::open(parent)?.sync_all()?;
    Ok(())
}

pub fn copy_bounded(
    source: &Path,
    target: &Path,
    expected: u64,
    check: &impl Fn() -> Result<()>,
) -> Result<()> {
    use std::io::{Read, Write};
    let mut input = File::open(source)?.take(expected);
    let mut output = File::create(target)?;
    let mut buffer = [0_u8; 65536];
    let mut written = 0_u64;
    let mut next_check = 0_u64;
    loop {
        if written >= next_check {
            check()?;
            next_check = written.saturating_add(8 * 1024 * 1024);
        }
        let count = input.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        output.write_all(&buffer[..count])?;
        written += count as u64;
    }
    if written != expected || input.into_inner().read(&mut [0_u8; 1])? != 0 {
        anyhow::bail!("source changed beyond its publication byte promise");
    }
    check()
}

#[cfg(test)]
mod publication_tests {
    use super::*;

    #[test]
    fn replacement_preserves_preexisting_shared_hardlinks() -> Result<()> {
        let temp = tempfile::TempDir::new()?;
        let output = temp.path().join("output");
        let shared = temp.path().join("shared");
        std::fs::write(&output, b"original")?;
        std::fs::hard_link(&output, &shared)?;
        publish_file(&output, |staged| {
            std::fs::write(staged, b"replacement")?;
            Ok(())
        })?;
        assert_eq!(std::fs::read(shared)?, b"original");
        assert_eq!(std::fs::read(output)?, b"replacement");
        Ok(())
    }

    #[test]
    fn unsuccessful_writer_keeps_previous_generation() -> Result<()> {
        let temp = tempfile::TempDir::new()?;
        let output = temp.path().join("output");
        std::fs::write(&output, b"original")?;
        assert!(
            publish_file(&output, |staged| {
                std::fs::write(staged, b"partial")?;
                anyhow::bail!("interrupted")
            })
            .is_err()
        );
        assert_eq!(std::fs::read(output)?, b"original");
        assert_eq!(std::fs::read_dir(temp.path())?.count(), 1);
        Ok(())
    }

    #[test]
    fn identifiers_preserve_nonversion_suffixes() {
        assert_eq!(normalize_feature_id("ENST1.2|GENE|label|"), "ENST1");
        assert_eq!(normalize_feature_id("ENST1.2_PAR_Y"), "ENST1_PAR_Y");
        assert_eq!(normalize_feature_id("custom.isoform"), "custom.isoform");
    }
}
