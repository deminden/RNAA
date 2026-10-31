//! Complete-copy admission, provenance and atomic generation publication.
use super::*;
use rnaa_core::resources::ExecutionResources;

pub(super) struct Publication<'a> {
    pub db: &'a Database,
    pub paths: &'a ProjectPaths,
    pub config: &'a ProjectConfig,
    pub run: &'a RunRecord,
    pub request: &'a rnaa_core::fingerprint::StageRequest<serde_json::Value>,
    pub observation: &'a DownloadObservation,
    pub admission: &'a ExecutionResources,
}
impl Publication<'_> {
    pub fn publish(
        self,
        transfer: &ProjectPaths,
        cached: Vec<VerifiedFile>,
    ) -> Result<Vec<VerifiedFile>> {
        let Self {
            db,
            paths,
            config,
            run,
            request,
            observation,
            admission,
        } = self;
        let digest = request.digest()?;
        let check = || db.check_job();
        if cached.iter().all(|file| {
            matches!(
                file.artifact_kind,
                ArtifactKind::FastqR1 | ArtifactKind::FastqR2 | ArtifactKind::FastqSingle
            )
        }) {
            let scan =
                crate::fastq::scan_request(crate::fastq::FastqLimits::try_from(&config.quant)?)?;
            let _guard = admission.reserve("fastq_scan", scan, &check)?;
            crate::fastq::FastqPlan::inspect(
                run.library_layout,
                &cached,
                crate::fastq::FastqLimits::try_from(&config.quant)?,
                &check,
            )?;
        }
        let publication_plan = resources::PublicationPlan::inspect(&cached, run)?;
        let publication = admission.reserve(
            "publication",
            publication_plan.request(&paths.generations_dir, &paths.root)?,
            &check,
        )?;
        let mut generation = Generation::begin(db, OutputGroup::Download, "download")?;
        let output_paths = generation.paths(paths)?;
        let mut verified = Vec::new();
        let source_versions = cached
            .iter()
            .map(|file| {
                let receipt_path = file.path.with_added_extension("rnaa-verified.json");
                let receipt: serde_json::Value =
                    resources::read_json(&receipt_path, publication_plan.metadata_bytes)?;
                Ok(json!({
                    "http":receipt["http"],
                    "resource_execution":receipt["resource_execution"]
                }))
            })
            .collect::<Result<Vec<serde_json::Value>>>()?;
        fs::create_dir_all(output_paths.raw_run_dir(&run.run_accession))?;
        for source in cached {
            check()?;
            let target = output_paths.raw_run_dir(&run.run_accession).join(
                source
                    .path
                    .file_name()
                    .context("download has no basename")?,
            );
            publish_file(&target, |temporary| {
                resources::copy_bounded(&source.path, temporary, source.bytes, &check)
            })?;
            let snapshot = ManifestArtifact::snapshot(source.artifact_kind.as_str(), &target)?;
            if snapshot.checksum != source.checksum || snapshot.bytes != source.bytes {
                bail!("download changed while preparing publication");
            }
            verified.push(VerifiedFile {
                path: target,
                ..source
            });
        }
        let manifest_path = output_paths.manifest_path("download", &run.run_accession);
        let manifest = StageManifest {
            schema_version:1, stage:"download".into(), id:run.run_accession.clone(), started_at:now_rfc3339(), finished_at:now_rfc3339(),
            parameters:json!({"request":request, "request_sha256":digest, "method":config.download.method,
                "http_sources":source_versions.iter().map(|v| &v["http"]).collect::<Vec<_>>(),
                "resource_execution":{"publication_estimate":publication_plan,"publication_reservation":publication.snapshot(),
                    "transfers":source_versions.iter().map(|v| &v["resource_execution"]).collect::<Vec<_>>()}}),
            tool_versions:BTreeMap::from([("native_http".into(), "reqwest-blocking".into())]),
            input_artifacts:ReqwestDownloader.select_remote_files(run, config.download.prefer).iter().map(|remote| ManifestArtifact {
                kind:remote.kind.as_str().into(), path:remote.url.clone(), checksum_type:if remote.md5.is_some() {"md5"} else {"none"}.into(),
                checksum:remote.md5.clone().unwrap_or_default(), bytes:remote.bytes.unwrap_or_default(),
            }).collect(),
            output_artifacts:verified.iter().map(|file| ManifestArtifact::snapshot(file.artifact_kind.as_str(), &file.path)).collect::<Result<_>>()?,
            notes:vec!["Complete run published atomically; HTTP prefix checkpoints remain outside the generation.".into(),
                "Authoritative MD5 is checked when supplied. Local SHA256 alone does not establish the publisher's content identity.".into()],
        };
        resources::write_json(&manifest_path, &manifest, publication_plan.metadata_bytes)?;
        let mut records = verified
            .iter()
            .map(|file| artifact_record(db, &run.run_accession, file))
            .collect::<Vec<_>>();
        records.push(artifact_record(
            db,
            &run.run_accession,
            &verified_snapshot(&ManifestArtifact::snapshot("MANIFEST", &manifest_path)?)?,
        ));
        let required = records
            .iter()
            .map(|record| PathBuf::from(&record.path))
            .collect::<Vec<_>>();
        generation.seal_request(&required, request, &[])?;
        observation.verify_current(db)?;
        assert_download_request(db, &run.run_accession, config.download.prefer, &digest)?;
        generation.commit(&records, |db| {
            observation.verify_catalog(db)?;
            assert_download_request(db, &run.run_accession, config.download.prefer, &digest)?;
            complete_download(db, &run.run_accession, &verified, false)
        })?;
        // This is private transfer cache, with independent inodes and no published
        // readers. The run writer lock remains held until process_run returns.
        if let Some(directory) = transfer.raw_run_dir(&run.run_accession).parent()
            && let Err(error) = fs::remove_dir_all(directory)
        {
            warn!("could not remove completed transfer cache: {error}");
        }
        Ok(verified)
    }
}
