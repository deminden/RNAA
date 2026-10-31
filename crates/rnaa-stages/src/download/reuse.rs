//! Scoped historical download selection and source authority observations.
use super::*;

pub(super) struct DownloadObservation {
    catalog: DownloadPreference,
    disk: Option<DownloadPreference>,
}
impl DownloadObservation {
    pub(super) fn capture(db: &Database) -> Result<Self> {
        let catalog: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
        Ok(Self {
            catalog: catalog.download.prefer,
            disk: Self::disk_preference(db)?,
        })
    }
    fn disk_preference(db: &Database) -> Result<Option<DownloadPreference>> {
        let path = ProjectPaths::new(
            db.path()
                .parent()
                .context("download database has no project root")?,
        )
        .config_path;
        if path.exists() {
            Ok(Some(ProjectConfig::load(&path)?.download.prefer))
        } else {
            Ok(None)
        }
    }
    pub(super) fn verify_catalog(&self, db: &Database) -> Result<()> {
        let current: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
        if current.download.prefer != self.catalog {
            bail!("download preference changed during the attempt; retry");
        }
        Ok(())
    }
    pub(super) fn verify_current(&self, db: &Database) -> Result<()> {
        self.verify_catalog(db)?;
        if Self::disk_preference(db)? != self.disk {
            bail!("download preference changed on disk during the attempt; retry");
        }
        Ok(())
    }
}

impl ReqwestDownloader {
    pub(super) fn reuse_download(
        &self,
        db: &Database,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        run: &RunRecord,
        request: &rnaa_core::fingerprint::StageRequest<serde_json::Value>,
        observation: &DownloadObservation,
    ) -> Result<Option<Vec<VerifiedFile>>> {
        let digest = request.digest()?;
        let selected = db.generation_head(OutputGroup::Download, &run.run_accession)?;
        for head in db.generation_candidates(
            OutputGroup::Download,
            "download",
            &run.run_accession,
            &digest,
        )? {
            let Some(prepared) = db.prepare_generation_reuse(&head)? else {
                continue;
            };
            let pinned = head.paths(paths)?;
            let Some(manifest) = fs::read(pinned.manifest_path("download", &run.run_accession))
                .ok()
                .and_then(|bytes| serde_json::from_slice::<StageManifest>(&bytes).ok())
            else {
                continue;
            };
            let remotes = self.select_remote_files(run, config.download.prefer);
            if manifest.schema_version != 1
                || manifest.stage != "download"
                || manifest.id != run.run_accession
                || manifest.parameters["request_sha256"] != digest
                || manifest.parameters["request"] != serde_json::to_value(request)?
                || manifest.output_artifacts.is_empty()
                || manifest.output_artifacts.len() != remotes.len()
            {
                continue;
            }
            let mut files = Vec::new();
            let mut validations = Vec::new();
            let mut used = std::collections::BTreeSet::new();
            for remote in remotes {
                let basename = sanitize_basename(&remote.basename);
                let Some(snapshot) = manifest.output_artifacts.iter().find(|file| {
                    Path::new(&file.path)
                        .file_name()
                        .is_some_and(|name| name == basename.as_str())
                }) else {
                    break;
                };
                if !used.insert(&snapshot.path)
                    || !Path::new(&snapshot.path)
                        .starts_with(pinned.raw_run_dir(&run.run_accession))
                    || snapshot.checksum_type != "sha256"
                {
                    break;
                }
                let Some(mut file) = verify_local_file(remote, Path::new(&snapshot.path))? else {
                    break;
                };
                if file.artifact_kind.as_str() != snapshot.kind
                    || file.checksum != snapshot.checksum
                    || file.bytes != snapshot.bytes
                {
                    break;
                }
                if remote.md5.is_none() && selected.as_ref() != Some(&head) {
                    let Some(receipt) =
                        manifest.parameters["http_sources"]
                            .as_array()
                            .and_then(|versions| {
                                versions.iter().find(|value| {
                                    value["url"] == remote.url
                                        && value["prefix_sha256"] == snapshot.checksum
                                        && value["committed_bytes"] == snapshot.bytes
                                })
                            })
                    else {
                        break;
                    };
                    let Some(validation) = crate::http_transfer::validate_saved_version(
                        &remote.url,
                        &snapshot.checksum,
                        snapshot.bytes,
                        remote.bytes,
                        receipt,
                        &config.download,
                        &|| db.check_job(),
                    )?
                    else {
                        break;
                    };
                    validations.push(validation);
                    file.integrity_source = "current_strong_etag+sealed_sha256".into();
                } else if remote.md5.is_some() {
                    validations.push(json!({"kind":"publisher_md5","url":remote.url,"md5":remote.md5,"sha256":snapshot.checksum,"bytes":snapshot.bytes}));
                } else {
                    validations.push(json!({"kind":"selected_snapshot","generation_id":head.generation_id,"url":remote.url,
                        "sha256":snapshot.checksum,"bytes":snapshot.bytes}));
                }
                files.push(file);
            }
            if files.len() != manifest.output_artifacts.len() {
                continue;
            }
            observation.verify_current(db)?;
            assert_download_request(db, &run.run_accession, config.download.prefer, &digest)?;
            db.reattach_generation(prepared,|db| {
                observation.verify_catalog(db)?;
                assert_download_request(db,&run.run_accession,config.download.prefer,&digest)?;
                complete_download(db,&run.run_accession,&files,true)?;
                db.append_event("download",Some(&run.run_accession),"download selected from verified history",json!({
                    "generation_id":head.generation_id,"request_sha256":digest,"source_validations":validations,
                }))
            })?;
            return Ok(Some(files));
        }
        Ok(None)
    }
}
