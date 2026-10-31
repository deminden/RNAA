//! Native preprocessing execution with optional owned phase admission.
use super::*;

impl Preprocessor for FasterpInProcessPreprocessor {
    fn preprocess(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<PreprocessArtifacts> {
        self.preprocess_with_resources(run, fastqs, paths, config, None, &|| Ok(()))
    }
}

impl FasterpInProcessPreprocessor {
    /// Production callers provide an owned host context; direct library callers
    /// may retain the compatible unadmitted API through the Preprocessor trait.
    #[allow(clippy::too_many_arguments)]
    pub fn preprocess_with_resources(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        paths: &ProjectPaths,
        config: &ProjectConfig,
        admission: Option<&rnaa_core::resources::ExecutionResources>,
        check: &impl Fn() -> Result<()>,
    ) -> Result<PreprocessArtifacts> {
        if let Some(cached) = reconcile_preprocess_artifacts(run, fastqs, paths, config)? {
            return Ok(cached);
        }
        let started_at = now_rfc3339();
        // Capture before decoding/concatenation, so staging cannot hide an input
        // mutation behind a fingerprint of bytes the engine did not consume.
        let mut limits = crate::fastq::FastqLimits::try_from(&config.quant)?;
        limits.max_decoded_bytes = config.quant.preprocess_max_decoded_bytes()?;
        let scanning = admission
            .map(|a| {
                a.reserve(
                    "scan",
                    resources::preflight(limits)
                        .map_err(rnaa_core::resources::ResourceAdmissionError::new)?,
                    check,
                )
            })
            .transpose()?;
        check()?;
        let input_artifacts = snapshot_inputs(fastqs)?;
        let plan = crate::fastq::FastqPlan::inspect(run.library_layout, fastqs, limits, check)
            .with_context(|| {
                format!(
                    "decoded preprocessing admission failed (preprocess_max_input_mb={})",
                    config.quant.preprocess_max_input_mb
                )
            })?;
        // Never retain a small scan claim while waiting to grow into the
        // engine claim: concurrent workers could otherwise deadlock each other.
        drop(scanning);
        let threads = admission.map_or(config.quant.preprocess_threads, |a| {
            config
                .quant
                .preprocess_threads
                .min(usize::try_from(a.limits().cpu_threads).unwrap_or(usize::MAX))
        });
        let mut native_options = NativeProcessOptions::for_threads(threads);
        native_options.max_record_bytes = usize::try_from(plan.summary().max_record_bytes)?
            .checked_add(4)
            .context("native record allowance overflows")?;
        native_options.min_record_bytes =
            usize::try_from(plan.summary().min_record_bytes.saturating_sub(4))?;
        native_options.max_input_records = Some(plan.summary().fragments);
        native_options.input_byte_limits = Some([
            plan.summary().r1_normalized_bytes,
            plan.summary().r2_normalized_bytes,
        ]);
        let estimate = resources::Estimate::for_plan(&plan, native_options)
            .map_err(rnaa_core::resources::ResourceAdmissionError::new)?;
        let staging_parent = std::env::temp_dir();
        let out_dir = paths.preprocess_run_dir(&run.run_accession, &config.quant.engine);
        let processing = admission
            .map(|a| {
                a.reserve(
                    "engine",
                    estimate
                        .request(&staging_parent, &out_dir)
                        .map_err(rnaa_core::resources::ResourceAdmissionError::new)?,
                    check,
                )
            })
            .transpose()?;
        check()?;
        let prepared = plan.prepare_in(&staging_parent, check)?;
        fs::create_dir_all(&out_dir)
            .with_context(|| format!("failed to create {}", out_dir.display()))?;
        let report_json = out_dir.join("fasterp.json");
        let (kind, read1, read2) = output_paths(run, &out_dir)?;
        let result = process_paths_with_options(
            &prepared.read1.display().to_string(),
            &read1.display().to_string(),
            prepared
                .read2
                .as_ref()
                .map(|path| path.display().to_string())
                .as_deref(),
            read2
                .as_ref()
                .map(|path| path.display().to_string())
                .as_deref(),
            native_options,
        )
        .map_err(|err| {
            anyhow!(
                "fasterp native preprocessing failed for {}: {err:#}",
                run.run_accession
            )
        })?;
        check()?;
        let mut report: serde_json::Value = serde_json::from_str(&result.json_report)?;
        report["execution"]["resource_estimate"] = serde_json::to_value(&estimate)?;
        if let Some(guard) = &processing {
            report["execution"]["resource_reservation"] = serde_json::to_value(guard.snapshot())?;
        }
        let mut metadata_remaining = estimate.metadata_bytes;
        resources::write_metadata(&report_json, &report, &mut metadata_remaining)?;
        if !outputs_ready(
            &read1,
            read2.as_deref(),
            &report_json,
            result.passed_reads == 0,
        ) {
            bail!(
                "fasterp outputs missing or inconsistent with passing reads for {}",
                run.run_accession
            );
        }
        let fastqs = build_preprocessed_fastqs(kind, &read1, read2.as_deref())?;
        verify_inputs(&input_artifacts)?;
        plan.verify_sources(check)?;
        let output_artifacts = fastqs
            .iter()
            .map(|file| file.path.as_path())
            .chain(std::iter::once(report_json.as_path()))
            .map(|path| ManifestArtifact::snapshot("PREPROCESS_OUTPUT", path))
            .collect::<Result<Vec<_>>>()?;
        resources::write_metadata(
            &out_dir.join("rnaa_preprocess_manifest.json"),
            &StageManifest {
                schema_version: 1,
                stage: "preprocess".into(),
                id: run.run_accession.clone(),
                started_at,
                finished_at: now_rfc3339(),
                parameters: scientific_parameters(run),
                tool_versions: BTreeMap::from([("fasterp".into(), fasterp_version())]),
                input_artifacts,
                output_artifacts,
                notes: vec![format!(
                    "All lanes validated and consumed; decoded preflight: {}",
                    serde_json::to_string(plan.summary())?
                )],
            },
            &mut metadata_remaining,
        )?;
        Ok(PreprocessArtifacts {
            run_accession: run.run_accession.clone(),
            out_dir,
            fastqs,
            report_json,
            tool_name: "fasterp".into(),
            tool_version: fasterp_version(),
            passed_reads: result.passed_reads,
            failed_reads: result.failed_reads,
            reused: false,
        })
    }
}
