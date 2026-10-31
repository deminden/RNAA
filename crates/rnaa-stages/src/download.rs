use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

use anyhow::{Context, Result, anyhow, bail};
use rnaa_core::config::{DownloadPreference, ProjectConfig};
use rnaa_core::db::Database;
use rnaa_core::generations::{Generation, OutputGroup};
use rnaa_core::jobs::{JobLease, JobOutcome};
use rnaa_core::manifest::{ManifestArtifact, StageManifest};
use rnaa_core::model::{
    ArtifactKind, ArtifactRecord, RemoteFile, RemoteFileKind, RunRecord, VerifiedFile,
};
use rnaa_core::paths::ProjectPaths;
use rnaa_core::state::RunState;
use rnaa_core::traits::Downloader;
use rnaa_core::util::{
    compute_md5, compute_sha256, dir_size, file_size, now_rfc3339, parse_size_bytes, publish_file,
    sanitize_basename,
};
use serde_json::json;
use tracing::{error, info, warn};
mod files;
mod publication;
mod resources;
mod reuse;
use files::{verify_existing_or_download, verify_local_file};
use reuse::DownloadObservation;

#[derive(Debug, Clone, Default)]
pub struct ReqwestDownloader;

impl ReqwestDownloader {
    pub fn fastq_recovery_available(
        &self,
        run: &RunRecord,
        preference: DownloadPreference,
    ) -> bool {
        let files = self.select_remote_files(run, preference);
        let singles = files
            .iter()
            .filter(|file| file.kind == RemoteFileKind::Fastq)
            .count();
        let first = files
            .iter()
            .filter(|file| file.kind == RemoteFileKind::FastqR1)
            .count();
        let second = files
            .iter()
            .filter(|file| file.kind == RemoteFileKind::FastqR2)
            .count();
        match run.library_layout {
            rnaa_core::LibraryLayout::Single => singles > 0 && singles == files.len(),
            rnaa_core::LibraryLayout::Paired => {
                first > 0 && first == second && first + second == files.len()
            }
            _ => false,
        }
    }

    pub fn source_request_digest(
        &self,
        run: &RunRecord,
        preference: DownloadPreference,
    ) -> Result<String> {
        download_request(run, preference)?.digest()
    }

    /// Recover one requested run through the normal owned, resumable publisher.
    pub fn recover_inputs(
        &self,
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        run: &str,
    ) -> Result<()> {
        let selected = db.get_run(run)?.context("recovery run disappeared")?;
        if !self.fastq_recovery_available(&selected, config.download.prefer) {
            bail!(
                "no complete downloadable FASTQ group for {run}; restore/register original reads or resolve a supported FASTQ source"
            );
        }
        self.process_run(db, paths, config, &selected)?;
        while db
            .get_run(run)?
            .is_some_and(|run| run.state == RunState::Downloading)
        {
            db.recover_jobs()?;
            thread::sleep(Duration::from_millis(200));
        }
        if !db.get_run(run)?.is_some_and(|run| {
            rnaa_core::eligibility::quantified(run.state)
                || matches!(
                    run.state,
                    RunState::Verified
                        | RunState::PreprocessDone
                        | RunState::QuantRunning
                        | RunState::QuantDone
                )
        }) {
            bail!(
                "input recovery for {run} did not complete; inspect its current attempt/error and retry"
            );
        }
        Ok(())
    }

    pub fn run_worker(
        &self,
        db: Database,
        paths: ProjectPaths,
        mut config: ProjectConfig,
        concurrency: usize,
        forever: bool,
        prefer: Option<DownloadPreference>,
    ) -> Result<()> {
        if let Some(preference) = prefer {
            config.download.prefer = preference;
        }
        config.download.concurrency = concurrency;
        config.download.validate()?;
        if concurrency == 0 {
            bail!("download concurrency must be positive");
        }
        let mut deep_check = true;
        let concurrency = concurrency.max(1);
        let storage_limit_bytes = configured_storage_limit_bytes(&config)?;
        loop {
            self.reconcile_completed(&db, &config, deep_check)?;
            deep_check = false;
            db.recover_jobs()?;
            db.recover_generations()?;
            let runs = db.list_active_runs_in_states(&[
                RunState::Resolved,
                RunState::DownloadFailed,
                RunState::VerifyFailed,
                RunState::Downloaded,
            ])?;
            if runs.is_empty() {
                if !db
                    .list_active_runs_in_states(&[RunState::Downloading])?
                    .is_empty()
                {
                    thread::sleep(Duration::from_millis(200));
                    continue;
                }
                if forever {
                    thread::sleep(Duration::from_secs(30));
                    continue;
                }
                break;
            }

            let queue = Arc::new(Mutex::new(runs));
            let failures = Arc::new(AtomicUsize::new(0));
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
                let failures = Arc::clone(&failures);
                handles.push(thread::spawn(move || {
                    loop {
                        // An empty queue is finished even when the final download
                        // crosses the active-storage threshold. Only queued work
                        // waits for completed-count cleanup to free space.
                        let run = {
                            let mut queue = queue.lock().expect("queue poisoned");
                            queue.pop()
                        };
                        let Some(run) = run else { break };
                        if let Some(limit_bytes) = storage_limit_bytes {
                            while project_active_storage_bytes(&paths) >= limit_bytes {
                                info!(
                                    "storage limit reached ({} >= {}), pausing downloads",
                                    project_active_storage_bytes(&paths),
                                    limit_bytes
                                );
                                thread::sleep(Duration::from_secs(30));
                            }
                        }
                        if let Err(err) = downloader.process_run(&db, &paths, &config, &run) {
                            error!("download failed for {}: {err:#}", run.run_accession);
                            failures.fetch_add(1, Ordering::Relaxed);
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
                let failed = failures.load(Ordering::Relaxed);
                if failed > 0 {
                    bail!("download failed for {failed} requested run(s)");
                }
                // A second finite downloader may be finishing other claims. Do
                // not announce completion while it still owns requested work.
                while !db
                    .list_active_runs_in_states(&[RunState::Downloading])?
                    .is_empty()
                {
                    db.recover_jobs()?;
                    thread::sleep(Duration::from_millis(200));
                }
                let failed = db
                    .list_active_runs_in_states(&[
                        RunState::DownloadFailed,
                        RunState::VerifyFailed,
                    ])?
                    .len();
                if failed > 0 {
                    bail!("download failed for {failed} requested run(s)");
                }
                break;
            }
            thread::sleep(Duration::from_secs(30));
        }
        Ok(())
    }

    fn process_run(
        &self,
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        run: &RunRecord,
    ) -> Result<()> {
        if !db
            .list_active_runs()?
            .iter()
            .any(|current| current.run_accession == run.run_accession)
        {
            return Ok(());
        }
        let Some(mut lease) = JobLease::try_acquire(
            db,
            "download",
            Some(&run.run_accession),
            &[
                RunState::Resolved,
                RunState::DownloadFailed,
                RunState::VerifyFailed,
                RunState::Downloaded,
            ],
        )?
        else {
            return Ok(());
        };
        let owned = lease.database();
        match self.process_run_owned(&owned, paths, config, run, lease.attempt_id()) {
            Ok(_) => lease.finish(JobOutcome::Succeeded, None),
            Err(error) => {
                let message = format!("{error:#}");
                // This handle cannot overwrite a newer owner, including during
                // error handling. Cancellation may intentionally reject these.
                let _ = owned.atomic(|db| {
                    db.set_run_state(&run.run_accession, RunState::DownloadFailed, Some(&message))?;
                    db.append_event(
                        "download",
                        Some(&run.run_accession),
                        "download failed",
                        json!({"error":message}),
                    )?;
                    db.complete_owned_job(JobOutcome::Failed, Some(&message))
                });
                lease.finish(JobOutcome::Failed, Some(&message))?;
                Err(error)
            }
        }
    }

    fn process_run_owned(
        &self,
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        queued: &RunRecord,
        attempt_id: &str,
    ) -> Result<Vec<VerifiedFile>> {
        let run = db
            .get_run(&queued.run_accession)?
            .context("claimed download run disappeared")?;
        let request = download_request(&run, config.download.prefer)?;
        let digest = request.digest()?;
        db.set_job_request(&digest)?;
        let observation = DownloadObservation::capture(db)?;
        if let Some(files) = self.reuse_download(db, paths, config, &run, &request, &observation)? {
            return Ok(files);
        }
        let admission = rnaa_core::resources::ExecutionResources::open(
            &config.resources,
            rnaa_core::resources::ReservationOwner {
                project_root: fs::canonicalize(&paths.root)?,
                attempt_id: attempt_id.into(),
                stage: "download".into(),
                run_accession: Some(run.run_accession.clone()),
            },
        )?;
        db.append_event(
            "download",
            Some(&run.run_accession),
            "download started",
            json!({"prefer":config.download.prefer, "method":config.download.method}),
        )?;
        let transfer = paths.for_download_transfer(&run.run_accession, &digest);
        let check = || db.check_job();
        let cached = self.download_files(&run, &transfer, config, Some(&admission), &check)?;
        publication::Publication {
            db,
            paths,
            config,
            run: &run,
            request: &request,
            observation: &observation,
            admission: &admission,
        }
        .publish(&transfer, cached)
    }

    fn reconcile_completed(&self, db: &Database, config: &ProjectConfig, deep: bool) -> Result<()> {
        db.with_idle_analysis(|db| {
            for queued in db.list_active_runs()? {
                if matches!(queued.state, RunState::Resolved | RunState::Downloading | RunState::DownloadFailed | RunState::VerifyFailed | RunState::Downloaded) { continue; }
                db.with_idle_run(&queued.run_accession, |db| {
                    let Some(run) = db.get_run(&queued.run_accession)? else { return Ok(()) };
                    let Some(head) = db.generation_head(OutputGroup::Download, &run.run_accession)? else { return Ok(()) };
                    if head.request_sha256 == download_request(&run, config.download.prefer)?.digest()?
                        && (!deep || db.verify_generation(OutputGroup::Download, &run.run_accession)?) { return Ok(()); }
                    let reason = "download source selection or verified generation changed; redownload required";
                    db.atomic(|db| {
                        db.set_run_state(&run.run_accession, RunState::Resolved, Some(reason))?;
                        for stage in ["normalize", "deseq2", "corr"] {
                            if db.get_project_stage_state(stage)?.is_some() { db.set_project_stage_state(stage, rnaa_core::ProjectStageState::Pending, Some(reason))?; }
                        }
                        db.append_event("download", Some(&run.run_accession), "completed download invalidated", json!({"reason":reason}))
                    })
                })?;
            }
            Ok(())
        })?;
        Ok(())
    }

    fn download_files(
        &self,
        run: &RunRecord,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        admission: Option<&rnaa_core::resources::ExecutionResources>,
        check: &impl Fn() -> Result<()>,
    ) -> Result<Vec<VerifiedFile>> {
        let remotes = self.select_remote_files(run, config.download.prefer);
        if remotes.is_empty() {
            bail!("no downloadable files resolved for {}", run.run_accession);
        }
        let mut names = std::collections::BTreeSet::new();
        for remote in &remotes {
            if matches!(remote.basename.as_str(), "" | "." | "..")
                || !names.insert(sanitize_basename(&remote.basename))
            {
                bail!("download basenames are empty or collide after sanitizing");
            }
            if let Some(md5) = &remote.md5
                && (md5.len() != 32 || !md5.bytes().all(|byte| byte.is_ascii_hexdigit()))
            {
                bail!("invalid authoritative MD5");
            }
        }
        let singles = remotes
            .iter()
            .filter(|file| file.kind == RemoteFileKind::Fastq)
            .count();
        let first = remotes
            .iter()
            .filter(|file| file.kind == RemoteFileKind::FastqR1)
            .count();
        let second = remotes
            .iter()
            .filter(|file| file.kind == RemoteFileKind::FastqR2)
            .count();
        if singles + first + second > 0
            && !matches!(run.library_layout,
            rnaa_core::LibraryLayout::Single if singles == remotes.len())
            && !(run.library_layout == rnaa_core::LibraryLayout::Paired
                && singles == 0
                && first > 0
                && first == second
                && first + second == remotes.len())
        {
            bail!(
                "download selection is not a complete FASTQ group for {} layout",
                run.library_layout
            );
        }
        let mut verified = Vec::new();
        for remote in remotes {
            check()?;
            let target = self.local_target_path(paths, run, remote);
            let part = target.with_added_extension("part");
            verified.push(verify_existing_or_download(
                remote,
                &target,
                &part,
                &config.download,
                admission,
                check,
            )?);
        }
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
        // Private per-file directories keep payload basenames from colliding with
        // another file's .part/checkpoint/receipt names.
        paths
            .raw_run_dir(&run.run_accession)
            .join(&basename)
            .join(basename)
    }
}

fn configured_storage_limit_bytes(config: &ProjectConfig) -> Result<Option<u64>> {
    let raw = config.storage.max_active_size.trim();
    if raw.is_empty() {
        return Ok(None);
    }
    Ok(Some(parse_size_bytes(raw)?))
}

fn project_active_storage_bytes(paths: &ProjectPaths) -> u64 {
    [
        "raw",
        "quant",
        "de",
        "corr",
        "refs",
        "metadata",
        ".rnaa_support",
        ".rnaa_generations",
        ".rnaa_transfers",
    ]
    .into_iter()
    .map(|directory| dir_size(&paths.root.join(directory)))
    .sum()
}

impl Downloader for ReqwestDownloader {
    fn ensure_downloaded(
        &self,
        run: &RunRecord,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<Vec<VerifiedFile>> {
        self.download_files(run, paths, config, None, &|| Ok(()))
    }
}

fn assert_download_request(
    db: &Database,
    run: &str,
    preference: DownloadPreference,
    expected: &str,
) -> Result<()> {
    let latest = db.get_run(run)?.context("download run disappeared")?;
    if download_request(&latest, preference)?.digest()? != expected {
        bail!("download source metadata changed during transfer");
    }
    Ok(())
}

fn download_request(
    run: &RunRecord,
    preference: DownloadPreference,
) -> Result<rnaa_core::fingerprint::StageRequest<serde_json::Value>> {
    let mut files = ReqwestDownloader.select_remote_files(run, preference);
    files.sort_by(|a, b| a.basename.cmp(&b.basename).then(a.url.cmp(&b.url)));
    Ok(rnaa_core::fingerprint::StageRequest {
        schema_version: 1,
        stage: "download".into(),
        algorithm_version: 2,
        engines: BTreeMap::new(),
        parameters: json!({"layout":run.library_layout, "preference":preference, "files":files}),
        inputs: vec![],
        dependencies: BTreeMap::new(),
    })
}

fn artifact_record(db: &Database, run: &str, file: &VerifiedFile) -> ArtifactRecord {
    ArtifactRecord {
        project_id: db.project_id().into(),
        run_accession: Some(run.into()),
        kind: file.artifact_kind,
        path: file.path.display().to_string(),
        blob_id: None,
        shared_path: None,
        checksum_type: file.checksum_type.clone(),
        checksum: file.checksum.clone(),
        bytes: file.bytes,
        created_at: now_rfc3339(),
    }
}

fn verified_snapshot(snapshot: &ManifestArtifact) -> Result<VerifiedFile> {
    Ok(VerifiedFile {
        artifact_kind: snapshot.kind.parse()?,
        path: snapshot.path.clone().into(),
        checksum_type: snapshot.checksum_type.clone(),
        checksum: snapshot.checksum.clone(),
        bytes: snapshot.bytes,
        integrity_source: "verified_generation_sha256".into(),
    })
}

fn complete_download(db: &Database, run: &str, files: &[VerifiedFile], reused: bool) -> Result<()> {
    db.set_run_state(run, RunState::Verified, None)?;
    db.append_event(
        "download",
        Some(run),
        "download verified",
        json!({"files":files.len(),
        "bytes":files.iter().map(|file| file.bytes).sum::<u64>(), "reused":reused}),
    )
}

#[cfg(test)]
mod tests;
