//! Native quantification, bounded input preflight and ownership checkpoints.
use super::*;

impl Quantifier for QuantIndexGroup {
    fn quantify(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        reference: &ReferenceBundle,
        paths: &ProjectPaths,
        config: &ProjectConfig,
    ) -> Result<QuantArtifacts> {
        self.quantify_checked(run, fastqs, reference, paths, config, &|| Ok(()))
    }
}

impl QuantIndexGroup {
    pub fn quantify_checked(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        reference: &ReferenceBundle,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        check: &impl Fn() -> Result<()>,
    ) -> Result<QuantArtifacts> {
        self.quantify_with_resources(run, fastqs, reference, paths, config, None, check)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn quantify_with_resources(
        &self,
        run: &RunRecord,
        fastqs: &[VerifiedFile],
        reference: &ReferenceBundle,
        paths: &ProjectPaths,
        config: &ProjectConfig,
        admission: Option<&rnaa_core::resources::ExecutionResources>,
        check: &impl Fn() -> Result<()>,
    ) -> Result<QuantArtifacts> {
        let mut effective = config.clone();
        if let Some(context) = admission {
            effective.quant.threads = effective
                .quant
                .threads
                .max(1)
                .min(usize::try_from(context.limits().cpu_threads).unwrap_or(usize::MAX));
        }
        let config = &effective;
        config.quant.validate()?;
        let match_cache =
            kallistors::pseudoalign::MatchCacheSize::from_bytes(config.quant.match_cache_bytes()?)
                .map_err(|error| anyhow::anyhow!("invalid quant match cache: {error}"))?;
        let cache_allocation = match_cache
            .allocation(config.quant.threads)
            .map_err(|error| anyhow::anyhow!("invalid quant match cache allocation: {error}"))?;
        let fastq_limits = crate::fastq::FastqLimits::try_from(&config.quant)?;
        let scanning =
            self.acquire_scan(admission, resources::scan_request(fastq_limits)?, check)?;
        let index_plan = self.preflight(reference, check)?;
        index_plan.check_memory_limit(
            config.resources.limits()?.memory_bytes,
            config.quant.threads,
        )?;
        let plan =
            crate::fastq::FastqPlan::inspect(run.library_layout, fastqs, fastq_limits, check)?;
        let estimate = resources::QuantResources::estimate(
            &index_plan,
            &plan,
            config.quant.threads,
            cache_allocation.total_bytes() as u64,
        )
        .map_err(rnaa_core::resources::ResourceAdmissionError::new)?;
        drop(scanning);
        let staging_parent = std::env::temp_dir();
        let out_dir = paths.quant_run_dir(&run.run_accession, &config.quant.engine);
        let index = self.acquire(
            reference,
            &estimate,
            estimate.run_request(&staging_parent, &out_dir),
            admission,
            check,
        )?;
        check()?;
        let prepared = plan.prepare_in(&staging_parent, check)?;
        let read2 = prepared.read2.as_ref();
        fs::create_dir_all(&out_dir)
            .with_context(|| format!("failed to create {}", out_dir.display()))?;

        let started_at = now_rfc3339();
        let start_time = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs()
            .to_string();
        check()?;

        let strand = match config.quant.strand {
            rnaa_core::config::QuantStrand::Unstranded => {
                kallistors::pseudoalign::Strand::Unstranded
            }
            rnaa_core::config::QuantStrand::Forward => kallistors::pseudoalign::Strand::Forward,
            rnaa_core::config::QuantStrand::Reverse => kallistors::pseudoalign::Strand::Reverse,
        };
        let options = kallistors::pseudoalign::PseudoalignOptions {
            match_cache_size: match_cache,
            ..Default::default()
        };
        let (
            ec_list,
            counts,
            reads_processed,
            reads_aligned,
            unique,
            bias,
            fragment_length_hist,
            fragment_length_stats,
        ) = std::thread::scope(|scope| -> Result<_> {
            std::thread::Builder::new()
                .spawn_scoped(scope, || -> Result<_> {
                    let mut reader1 = input::open(&prepared.read1, &prepared.summary, 0)?;
                    let ec_counts = if let Some(read2) = read2 {
                        let mut reader2 = input::open(read2, &prepared.summary, 1)?;
                        kallistors::pseudoalign::pseudoalign_paired_bifrost_with_options_threaded(
                            &index,
                            &mut reader1,
                            &mut reader2,
                            strand,
                            options,
                            config.quant.threads.max(1),
                        )
                    } else {
                        kallistors::pseudoalign::pseudoalign_single_end_bifrost_with_options_threaded(
                            &index,
                            &mut reader1,
                            strand,
                            None,
                            options,
                            config.quant.threads.max(1),
                        )
                    };
                    let ec_counts = ec_counts.map_err(|err| {
                        anyhow::anyhow!("kallistors pseudoalignment failed: {err}")
                    })?;
                    let unique = kallistors::pseudoalign::unique_pseudoaligned_reads(&ec_counts);
                    Ok((
                        ec_counts.ec_list,
                        ec_counts.counts,
                        ec_counts.reads_processed,
                        ec_counts.reads_aligned,
                        unique,
                        ec_counts.bias,
                        ec_counts.fragment_length_hist,
                        ec_counts.fragment_length_stats,
                    ))
                })?
                .join()
                .map_err(|_| anyhow::anyhow!("quant mapping thread panicked"))?
        })?;
        let ec_input = kallistors::quant::EcCountsInput {
            ec_list: kallistors::ec::EcList { classes: ec_list },
            counts,
        };
        let (mean_fragment_length, fragment_length_sd) = if read2.is_some() {
            // kallistors estimates paired effective lengths from the observed fragment histogram.
            let mean = fragment_length_stats
                .as_ref()
                .and_then(|stats| stats.mean())
                .unwrap_or(200.0);
            let sd = fragment_length_stats
                .as_ref()
                .and_then(|stats| stats.sd())
                .unwrap_or(20.0);
            (mean, sd)
        } else {
            // Preserve kallisto's conservative single-end defaults used by RNAA's previous stage.
            (
                config.quant.single_fragment_length,
                config.quant.single_fragment_sd,
            )
        };
        check()?;
        let quant_options = kallistors::quant::QuantOptions {
            mean_fragment_length,
            fragment_length_sd,
            ..kallistors::quant::QuantOptions::default()
        };
        let result = kallistors::quant::em_quantify_with_fragment_histogram_threaded(
            &ec_input,
            &index.transcript_lengths,
            None,
            bias.as_ref(),
            fragment_length_hist.as_deref(),
            quant_options,
            config.quant.threads.max(1),
        )
        .map_err(|err| anyhow::anyhow!("kallistors quantification failed: {err}"))?;

        check()?;
        let abundance_tsv = out_dir.join("abundance.tsv");
        let abundance_h5 = out_dir.join("abundance.h5");
        let run_info_json = out_dir.join("run_info.json");
        let manifest_path = out_dir.join("rnaa_quant_manifest.json");
        let target_names = index
            .transcript_names
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        kallistors::quant::write_abundance_tsv_from_parts(
            &abundance_tsv,
            &target_names,
            &index.transcript_lengths,
            &result,
        )
        .map_err(|err| anyhow::anyhow!("failed to write kallistors abundance table: {err}"))?;
        let call = format!("rnaa quant {}", run.run_accession);
        rnaa_core::util::publish_file(&abundance_h5, |staged| {
            kallistors::quant::write_abundance_h5(
                staged,
                kallistors::quant::H5QuantOutput {
                    target_names: &target_names,
                    target_lengths: &index.transcript_lengths,
                    result: &result,
                    bootstrap_est_counts: &[],
                    reads_processed,
                    fragment_length_hist: fragment_length_hist.as_deref(),
                    bias_observed: bias.as_ref().map(|value| value.counts.as_slice()),
                    index_version: index.index_version,
                    kallisto_version: &format!("kallistors {}", crate::engines::KALLISTORS),
                    call: &call,
                    start_time: &start_time,
                },
            )
            .map_err(|err| anyhow::anyhow!("failed to write kallistors abundance HDF5: {err}"))
        })?;
        kallistors::quant::RunInfo {
            n_targets: index.transcript_names.len(),
            n_bootstraps: 0,
            n_processed: reads_processed,
            n_pseudoaligned: reads_aligned,
            n_unique: unique,
            p_pseudoaligned: if reads_processed == 0 {
                0.0
            } else {
                reads_aligned as f64 * 100.0 / reads_processed as f64
            },
            p_unique: if reads_processed == 0 {
                0.0
            } else {
                unique as f64 * 100.0 / reads_processed as f64
            },
            kallisto_version: format!("kallistors {}", crate::engines::KALLISTORS),
            index_version: index.index_version,
            kmer_length: index.k as u32,
            start_time,
            call,
        }
        .write_json(&run_info_json)
        .map_err(|err| anyhow::anyhow!("failed to write kallistors run info: {err}"))?;

        let upstream_inputs = if matches!(
            prepared.read1_kind,
            rnaa_core::ArtifactKind::FastqTrimmedR1 | rnaa_core::ArtifactKind::FastqTrimmedSingle
        ) {
            fs::read(
                paths
                    .preprocess_run_dir(&run.run_accession, &config.quant.engine)
                    .join("rnaa_preprocess_manifest.json"),
            )
            .ok()
            .and_then(|raw| serde_json::from_slice::<StageManifest>(&raw).ok())
            .map(|manifest| manifest.input_artifacts)
            .unwrap_or_default()
        } else {
            Vec::new()
        };
        let tool_versions = BTreeMap::from([
            (
                "kallistors".to_string(),
                crate::engines::KALLISTORS.to_string(),
            ),
            ("rnaa".to_string(), env!("CARGO_PKG_VERSION").to_string()),
        ]);
        self.verify(reference)?;
        plan.verify_sources(check)?;
        let index_execution = self.stats();
        let manifest = StageManifest {
            schema_version: 1,
            stage: "quant".to_string(),
            id: run.run_accession.clone(),
            started_at,
            finished_at: now_rfc3339(),
            parameters: json!({
                "engine": "kallistors",
                "engine_version": crate::engines::KALLISTORS,
                "algorithm_version": 3,
                "preprocess_requested": config.quant.preprocess,
                "strand": config.quant.strand,
                "threads": config.quant.threads.max(1),
                "match_cache": {
                    "scope": "per_run_across_mapping_threads",
                    "requested_bytes": cache_allocation.requested_bytes(),
                    "allocated_bytes": cache_allocation.total_bytes(),
                    "mapping_workers": cache_allocation.workers(),
                },
                "layout": if read2.is_some() { "paired" } else { "single" },
                "single_fragment_length": if read2.is_none() { Some(config.quant.single_fragment_length) } else { None },
                "single_fragment_sd": if read2.is_none() { Some(config.quant.single_fragment_sd) } else { None },
                "index": reference.kallisto_index_path.display().to_string(),
                "index_sha256": index_execution.index_sha256,
                "index_execution": index_execution,
                "resource_execution": {"estimate":estimate,"reservations":index.snapshots()},
                "fastq_execution": prepared.summary,
                "reference_manifest": reference.manifest_path.display().to_string()
            }),
            tool_versions,
            input_artifacts: plan.inputs().iter().cloned().chain(upstream_inputs).collect(),
            output_artifacts: [&abundance_h5, &abundance_tsv, &run_info_json].into_iter().map(|path| {
                Ok(ManifestArtifact {
                    kind: path.file_name().unwrap_or_default().to_string_lossy().to_ascii_uppercase(),
                    path: path.display().to_string(),
                    checksum_type: "sha256".to_string(),
                    checksum: compute_sha256(path)?,
                    bytes: file_size(path).unwrap_or_default(),
                })
            }).collect::<Result<Vec<_>>>()?,
            notes: vec!["Quantification performed in-process by kallistors; abundance.tsv is retained as the tximport-compatible transcript table.".to_string()],
        };
        write_json_pretty(&manifest_path, &manifest)?;

        Ok(QuantArtifacts {
            run_accession: run.run_accession.clone(),
            out_dir,
            abundance_h5,
            run_info_json,
            manifest_path,
        })
    }
}
