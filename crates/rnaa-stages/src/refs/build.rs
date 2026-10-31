//! Acquired-source planning, admitted native work and fenced publication.
use super::*;

impl EnsemblReferenceManager {
    pub(super) fn build_reference(
        &self,
        paths: &ProjectPaths,
        db: &Database,
        parent: &Database,
        config: &ProjectConfig,
        preparation: &Preparation<'_>,
        attempt_id: &str,
    ) -> Result<ReferenceBundle> {
        let spec = preparation.spec;
        let expected = preparation.expected;
        preparation.observation.verify_current(db)?;
        let started = now_rfc3339();
        let check = || {
            parent.check_job()?;
            db.check_job()
        };
        let admission = rnaa_core::resources::ExecutionResources::open(
            &config.resources,
            rnaa_core::resources::ReservationOwner {
                project_root: fs::canonicalize(&paths.root)?,
                attempt_id: attempt_id.into(),
                stage: "refs".into(),
                run_accession: None,
            },
        )?;
        let cdna = acquire_source(
            &spec.transcript_source,
            ArtifactKind::ReferenceCdna,
            paths,
            config,
            &spec.reference_id,
            &admission,
            &check,
        )?;
        let gtf = acquire_source(
            &spec.gtf_source,
            ArtifactKind::ReferenceGtf,
            paths,
            config,
            &spec.reference_id,
            &admission,
            &check,
        )?;
        check()?;
        verify_originals(spec)?;
        let setup = admission.reserve(
            "planning_setup",
            rnaa_core::resources::ResourceRequest {
                cpu_threads: 1,
                memory_bytes: 16 * 1024 * 1024,
                disks: vec![
                    rnaa_core::resources::DiskRequest {
                        path: paths.generations_dir.clone(),
                        bytes: 4 * 1024 * 1024,
                    },
                    rnaa_core::resources::DiskRequest {
                        path: paths.root.clone(),
                        bytes: 4 * 1024 * 1024,
                    },
                ],
            },
            &check,
        )?;
        let mut generation = Generation::begin(db, OutputGroup::Refs, "refs")?;
        let pinned = generation.paths(paths)?;
        let result = bundle(&pinned, spec);
        let workspace = generation.head().payload_path.join(".rnaa_reference_work");
        fs::create_dir_all(&workspace)?;
        drop(setup);
        let plan = planning::Plan::inspect(
            &cdna.input,
            &gtf.input,
            &config.refs,
            &admission,
            &workspace,
            &check,
        )?;
        let build_guard = admission.reserve(
            "build",
            plan.request(&paths.generations_dir, &paths.root, &workspace)?,
            &check,
        )?;
        fs::create_dir_all(pinned.reference_dir(&spec.reference_id))?;
        cdna.copy_to(&result.cdna_path, &check)?;
        gtf.copy_to(&result.gtf_path, &check)?;
        let inputs = source_snapshots(&result)?;
        let request = reference_request(spec, config, &inputs);
        let digest = request.digest()?;
        db.set_job_request(&digest)?;
        build_annotation_tables(
            &result.gtf_path,
            &result.tx2gene_path,
            result.gene_annotation_path.as_ref().unwrap(),
            &plan.gtf,
            &check,
        )?;
        check()?;
        let settings = IndexSettings::from_config(config);
        let staging = if plan.staging_bytes > 0 {
            Some(
                tempfile::Builder::new()
                    .prefix("rnaa-reference-decoded-")
                    .tempdir_in(&workspace)?,
            )
        } else {
            None
        };
        let fasta_path = if let Some(directory) = &staging {
            let path = directory.path().join("transcriptome.fa");
            plan.stage_fasta(&result.cdna_path, &path, &check)?;
            path
        } else {
            result.cdna_path.clone()
        };
        let options = settings.options(plan.threads);
        let build_report = std::thread::scope(|scope| -> Result<_> {
            std::thread::Builder::new()
                .spawn_scoped(scope, || {
                    kallistors::index::build_index_with_report(
                        &result.kallisto_index_path,
                        std::slice::from_ref(&fasta_path),
                        options,
                    )
                })?
                .join()
                .map_err(|_| anyhow!("reference index thread panicked"))?
                .map_err(|err| anyhow!("kallistors index build failed: {err}"))
        })?;
        check()?;
        plan.validate_outputs(&result, &build_report)?;
        crate::quant::QuantIndexPlan::inspect(&result.kallisto_index_path, &check)?;
        let mapping = validate_index_mapping(&result)?;
        for input in &inputs {
            if !input.matches_file()? {
                bail!("reference source snapshot changed during index construction");
            }
        }
        let outputs = bundle_outputs(&result);
        let manifest=StageManifest {schema_version:1,stage:"refs".into(),id:spec.reference_id.clone(),started_at:started,finished_at:now_rfc3339(),
            parameters:json!({"request":request,"request_sha256":digest,"transcript_source":spec.transcript_source,
                "gtf_source":spec.gtf_source,"checksum_catalog":spec.checksum_catalog,"mapping":mapping,"index_threads":plan.threads,"requested_index_threads":config.refs.index_threads,
                "http_sources":[cdna.http,gtf.http],
                "resource_execution":{"estimate":plan,"build_reservation":build_guard.snapshot(),"source_transfers":[cdna.execution,gtf.execution],
                    "build_report":{"transcripts":build_report.transcripts,"unitigs":build_report.unitigs,"kmers":build_report.kmers,"minimizers":build_report.minimizers,"edge_scratch_bytes":build_report.edge_scratch_bytes}}}),
            tool_versions:BTreeMap::from([("kallistors".into(),crate::engines::KALLISTORS.into())]),input_artifacts:inputs,
            output_artifacts:outputs.iter().filter(|(kind,_)|*kind!=ArtifactKind::Manifest).map(|(kind,path)|ManifestArtifact::snapshot(kind.as_str(),path)).collect::<Result<_>>()?,
            notes:vec!["Sources, annotation tables and index publish as one immutable reference generation.".into()]};
        if serde_json::to_vec_pretty(&manifest)?.len() as u64 > plan.metadata_bytes {
            bail!("reference manifest exceeds its planned metadata write allowance");
        }
        write_json_pretty(&result.manifest_path, &manifest)?;
        let records = outputs
            .iter()
            .map(|(kind, path)| {
                let artifact = ManifestArtifact::snapshot(kind.as_str(), path)?;
                Ok(ArtifactRecord {
                    project_id: db.project_id().into(),
                    run_accession: None,
                    kind: *kind,
                    path: artifact.path,
                    blob_id: None,
                    shared_path: None,
                    checksum_type: artifact.checksum_type,
                    checksum: artifact.checksum,
                    bytes: artifact.bytes,
                    created_at: now_rfc3339(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        verify_originals(spec)?;
        drop(staging);
        fs::remove_dir(&workspace)?;
        generation.seal_request(
            &outputs
                .iter()
                .map(|(_, path)| path.clone())
                .collect::<Vec<_>>(),
            &request,
            &[cdna.input.clone(), gtf.input.clone()],
        )?;
        let selection = ReferenceSelection::capture(&result, generation.head().clone(), config)?;
        preparation.observation.verify_current(db)?;
        let generation_id = generation.head().generation_id.clone();
        generation.commit(&records,|db| {
            preparation.observation.verify_catalog(db)?;
            parent.check_job()?;
            db.select_reference(&selection,expected)?;
            db.append_event("refs",None,"reference bundle prepared",json!({"reference_id":result.id,"generation_id":generation_id,"request_sha256":digest}))
        })?;
        // Retain failed transfers for retry; successful generations own separate copies.
        for source in [cdna, gtf] {
            if let Err(error) = source.discard_transfer() {
                tracing::warn!("reference transfer cache cleanup failed: {error:#}");
            }
        }
        Ok(result)
    }
}
