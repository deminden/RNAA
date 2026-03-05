use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

use anyhow::{Context, Result, anyhow, bail};
use rnaa_core::config::{DownloadMethod, DownloadPreference, ProjectConfig};
use rnaa_core::db::Database;
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{
    ArtifactKind, ArtifactRecord, RemoteFile, RemoteFileKind, RunRecord, VerifiedFile,
};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::RunState;
use rnaa_core::traits::Downloader;
use rnaa_core::util::{
    compute_md5, compute_sha256, file_size, now_rfc3339, sanitize_basename, write_json_pretty,
};
use serde_json::json;
use tracing::{error, info, warn};

#[derive(Debug, Clone, Default)]
pub struct ShellDownloader;

#[derive(Debug, Clone, PartialEq, Eq)]
struct DownloaderCommand {
    binary: String,
    args: Vec<String>,
}

impl ShellDownloader {
    pub fn run_worker(
        &self,
        db: Database,
        paths: ProjectPaths,
        config: ProjectConfig,
        concurrency: usize,
        forever: bool,
        prefer: Option<DownloadPreference>,
    ) -> Result<()> {
        let concurrency = concurrency.max(1);
        loop {
            let runs = db.list_runs_in_states(&[
                RunState::Resolved,
                RunState::DownloadFailed,
                RunState::VerifyFailed,
                RunState::Downloaded,
            ])?;
            if runs.is_empty() {
                if forever {
                    thread::sleep(Duration::from_secs(30));
                    continue;
                }
                break;
            }

            let queue = Arc::new(Mutex::new(runs));
            let mut handles = Vec::new();
            for _ in 0..concurrency {
                let queue = Arc::clone(&queue);
                let db = db.clone();
                let paths = paths.clone();
                let mut config = config.clone();
                if let Some(prefer) = prefer {
                    config.download.prefer = prefer;
                }
                let downloader = self.clone();
                handles.push(thread::spawn(move || {
                    loop {
                        let run = {
                            let mut queue = queue.lock().expect("queue poisoned");
                            queue.pop()
                        };
                        let Some(run) = run else { break };
                        if let Err(err) = downloader.process_run(&db, &paths, &config, &run) {
                            error!("download failed for {}: {err:#}", run.run_accession);
                            let _ = db.set_run_state(
                                &run.run_accession,
                                RunState::DownloadFailed,
                                Some(&format!("{err:#}")),
                            );
                            let _ = db.append_event(
                                "download",
                                Some(&run.run_accession),
                                "download failed",
                                json!({"error": format!("{err:#}")}),
                            );
                        }
                    }
                }));
            }

            for handle in handles {
                handle
                    .join()
                    .map_err(|_| anyhow!("download thread panicked"))?;
            }

            if !forever {
                break;
            }
        }
        Ok(())
    }

    fn process_run(
        &self,
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        run: &RunRecord,
    ) -> Result<Vec<VerifiedFile>> {
        db.set_run_state(&run.run_accession, RunState::Downloading, None)?;
        db.append_event(
            "download",
            Some(&run.run_accession),
            "download started",
            json!({"prefer": config.download.prefer, "method": config.download.method}),
        )?;

        let verified = self.ensure_downloaded(run, paths, config)?;
        for item in &verified {
            let artifact = ArtifactRecord {
                project_id: db.project_id().to_string(),
                run_accession: Some(run.run_accession.clone()),
                kind: item.artifact_kind,
                path: item.path.display().to_string(),
                checksum_type: item.checksum_type.clone(),
                checksum: item.checksum.clone(),
                bytes: item.bytes,
                created_at: now_rfc3339(),
            };
            db.record_artifact(&artifact)?;
        }

        db.set_run_state(&run.run_accession, RunState::Verified, None)?;
        let total_bytes = verified.iter().map(|item| item.bytes).sum::<u64>();
        db.append_event(
            "download",
            Some(&run.run_accession),
            "download verified",
            json!({"files": verified.len(), "bytes": total_bytes}),
        )?;
        Ok(verified)
    }

    fn select_remote_files<'a>(
        &self,
        run: &'a RunRecord,
        preference: DownloadPreference,
    ) -> Vec<&'a RemoteFile> {
        let fastqs = run
            .remote_files
            .iter()
            .filter(|file| {
                matches!(
                    file.kind,
                    RemoteFileKind::Fastq | RemoteFileKind::FastqR1 | RemoteFileKind::FastqR2
                )
            })
            .collect::<Vec<_>>();
        let sras = run
            .remote_files
            .iter()
            .filter(|file| file.kind == RemoteFileKind::Sra)
            .collect::<Vec<_>>();

        match preference {
            DownloadPreference::Fastq if !fastqs.is_empty() => fastqs,
            DownloadPreference::Sra if !sras.is_empty() => sras,
            DownloadPreference::Fastq if !sras.is_empty() => sras,
            DownloadPreference::Sra if !fastqs.is_empty() => fastqs,
            _ => run.remote_files.iter().collect(),
        }
    }

    fn local_target_path(
        &self,
        paths: &ProjectPaths,
        run: &RunRecord,
        remote: &RemoteFile,
    ) -> PathBuf {
        let basename = sanitize_basename(&remote.basename);
        paths.raw_run_dir(&run.run_accession).join(basename)
    }

    fn download_to_part(&self, method: DownloadMethod, url: &str, part_path: &Path) -> Result<()> {
        let command = build_download_command(method, url, part_path);
        let status = Command::new(&command.binary)
            .args(&command.args)
            .status()
            .with_context(|| format!("failed to spawn {}", command.binary))?;
        if !status.success() {
            bail!("downloader exited with status {status}");
        }
        Ok(())
    }
}

impl Downloader for ShellDownloader {
    fn ensure_downloaded(
        &self,
        run: &RunRecord,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<Vec<VerifiedFile>> {
        let remote_files = self.select_remote_files(run, config.download.prefer);
        if remote_files.is_empty() {
            bail!("no downloadable files resolved for {}", run.run_accession);
        }

        let run_dir = paths.raw_run_dir(&run.run_accession);
        fs::create_dir_all(&run_dir)
            .with_context(|| format!("failed to create {}", run_dir.display()))?;
        let started_at = now_rfc3339();
        let mut outputs = Vec::new();
        let mut manifest_outputs = Vec::new();

        for remote in &remote_files {
            let target_path = self.local_target_path(paths, run, remote);
            let part_path = target_path.with_extension(format!(
                "{}.part",
                target_path
                    .extension()
                    .and_then(|ext| ext.to_str())
                    .unwrap_or_default()
            ));
            let mut last_error = None;

            for attempt in 0..=config.download.retries {
                match verify_existing_or_download(
                    self,
                    config.download.method,
                    remote,
                    &target_path,
                    &part_path,
                ) {
                    Ok(verified) => {
                        manifest_outputs.push(ManifestArtifact {
                            kind: verified.artifact_kind.as_str().to_string(),
                            path: verified.path.display().to_string(),
                            checksum_type: verified.checksum_type.clone(),
                            checksum: verified.checksum.clone(),
                            bytes: verified.bytes,
                        });
                        outputs.push(verified);
                        last_error = None;
                        break;
                    }
                    Err(err) => {
                        last_error = Some(format!("{err:#}"));
                        warn!(
                            "download attempt {} failed for {}: {err:#}",
                            attempt + 1,
                            remote.url
                        );
                        if attempt < config.download.retries {
                            let backoff = config
                                .download
                                .backoff_seconds
                                .get(attempt as usize)
                                .copied()
                                .unwrap_or(30);
                            thread::sleep(Duration::from_secs(backoff));
                        }
                    }
                }
            }

            if let Some(error) = last_error {
                return Err(anyhow!(
                    "failed to download {} after retries: {error}",
                    remote.url
                ));
            }
        }

        let manifest = StageManifest {
            schema_version: 1,
            stage: "download".to_string(),
            id: run.run_accession.clone(),
            started_at,
            finished_at: now_rfc3339(),
            parameters: json!({
                "preference": config.download.prefer,
                "method": config.download.method,
                "retries": config.download.retries,
                "backoff_seconds": config.download.backoff_seconds,
            }),
            tool_versions: BTreeMap::from_iter(
                [
                    ("curl".to_string(), rnaa_core::util::command_version("curl", &["--version"])),
                    ("wget".to_string(), rnaa_core::util::command_version("wget", &["--version"])),
                ]
                .into_iter()
                .filter_map(|(name, version)| version.map(|version| (name, version))),
            ),
            input_artifacts: remote_files
                .iter()
                .map(|remote| ManifestArtifact {
                    kind: remote.kind.as_str().to_string(),
                    path: remote.url.clone(),
                    checksum_type: if remote.md5.is_some() {
                        "md5".to_string()
                    } else {
                        "none".to_string()
                    },
                    checksum: remote.md5.clone().unwrap_or_default(),
                    bytes: remote.bytes.unwrap_or_default(),
                })
                .collect(),
            output_artifacts: manifest_outputs,
            notes: vec![
                "Downloads use a .part file and atomic rename.".to_string(),
                "If authoritative md5 is unavailable, RNAA records the local sha256 plus any remote size check."
                    .to_string(),
            ],
        };
        write_json_pretty(
            &paths.manifest_path("download", &run.run_accession),
            &manifest,
        )?;

        Ok(outputs)
    }
}

fn build_download_command(
    method: DownloadMethod,
    url: &str,
    part_path: &Path,
) -> DownloaderCommand {
    match method {
        DownloadMethod::Curl => DownloaderCommand {
            binary: "curl".to_string(),
            args: vec![
                "--fail".to_string(),
                "--location".to_string(),
                "--continue-at".to_string(),
                "-".to_string(),
                "--retry".to_string(),
                "0".to_string(),
                "--output".to_string(),
                part_path.display().to_string(),
                url.to_string(),
            ],
        },
        DownloadMethod::Wget => DownloaderCommand {
            binary: "wget".to_string(),
            args: vec![
                "-c".to_string(),
                "-O".to_string(),
                part_path.display().to_string(),
                url.to_string(),
            ],
        },
    }
}

fn verify_existing_or_download(
    downloader: &ShellDownloader,
    method: DownloadMethod,
    remote: &RemoteFile,
    target_path: &Path,
    part_path: &Path,
) -> Result<VerifiedFile> {
    if let Some(existing) = verify_local_file(remote, target_path)? {
        return Ok(existing);
    }

    if part_path.exists() {
        info!("resuming partial download {}", part_path.display());
    }
    downloader.download_to_part(method, &remote.url, part_path)?;
    fs::rename(part_path, target_path).with_context(|| {
        format!(
            "failed to rename {} to {}",
            part_path.display(),
            target_path.display()
        )
    })?;
    match verify_local_file(remote, target_path)? {
        Some(verified) => Ok(verified),
        None => bail!(
            "downloaded file {} did not pass verification",
            target_path.display()
        ),
    }
}

fn verify_local_file(remote: &RemoteFile, path: &Path) -> Result<Option<VerifiedFile>> {
    if !path.exists() {
        return Ok(None);
    }
    let bytes = file_size(path)?;
    if let Some(expected_bytes) = remote.bytes
        && bytes != expected_bytes
    {
        return Ok(None);
    }

    let (checksum_type, checksum, integrity_source) = if let Some(expected_md5) = &remote.md5 {
        let actual_md5 = compute_md5(path)?;
        if actual_md5 != *expected_md5 {
            return Ok(None);
        }
        (
            "md5".to_string(),
            actual_md5,
            "authoritative_md5".to_string(),
        )
    } else {
        let sha256 = compute_sha256(path)?;
        let integrity_source = if remote.bytes.is_some() {
            "remote_size+local_sha256"
        } else {
            "local_sha256_only"
        };
        ("sha256".to_string(), sha256, integrity_source.to_string())
    };

    Ok(Some(VerifiedFile {
        artifact_kind: match remote.kind {
            RemoteFileKind::FastqR1 => ArtifactKind::FastqR1,
            RemoteFileKind::FastqR2 => ArtifactKind::FastqR2,
            RemoteFileKind::Fastq => ArtifactKind::FastqSingle,
            RemoteFileKind::Sra | RemoteFileKind::Submitted => ArtifactKind::Sra,
        },
        path: path.to_path_buf(),
        checksum_type,
        checksum,
        bytes,
        integrity_source,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rnaa_core::LibraryLayout;

    #[test]
    fn remote_selection_prefers_fastq() {
        let run = RunRecord {
            project_id: "p".to_string(),
            study_accession: "SRP1".to_string(),
            source_study_accession: None,
            run_accession: "SRR1".to_string(),
            sample_accession: None,
            experiment_accession: None,
            library_layout: LibraryLayout::Paired,
            instrument_platform: None,
            instrument_model: None,
            remote_files: vec![
                RemoteFile {
                    url: "https://example.com/a_1.fastq.gz".to_string(),
                    md5: None,
                    bytes: None,
                    kind: RemoteFileKind::FastqR1,
                    basename: "a_1.fastq.gz".to_string(),
                },
                RemoteFile {
                    url: "https://example.com/a_2.fastq.gz".to_string(),
                    md5: None,
                    bytes: None,
                    kind: RemoteFileKind::FastqR2,
                    basename: "a_2.fastq.gz".to_string(),
                },
                RemoteFile {
                    url: "https://example.com/a.sra".to_string(),
                    md5: None,
                    bytes: None,
                    kind: RemoteFileKind::Sra,
                    basename: "a.sra".to_string(),
                },
            ],
            metadata: BTreeMap::new(),
            state: RunState::Resolved,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        let selected = ShellDownloader.select_remote_files(&run, DownloadPreference::Fastq);
        assert_eq!(selected.len(), 2);
        assert_eq!(selected[0].kind, RemoteFileKind::FastqR1);
    }

    #[test]
    fn build_curl_command_template() {
        let part_path = Path::new("/tmp/SRR1.fastq.gz.part");
        let command = build_download_command(
            DownloadMethod::Curl,
            "https://example.com/a.fastq.gz",
            part_path,
        );
        assert_eq!(command.binary, "curl");
        assert_eq!(
            command.args,
            vec![
                "--fail",
                "--location",
                "--continue-at",
                "-",
                "--retry",
                "0",
                "--output",
                "/tmp/SRR1.fastq.gz.part",
                "https://example.com/a.fastq.gz",
            ]
        );
    }
}
