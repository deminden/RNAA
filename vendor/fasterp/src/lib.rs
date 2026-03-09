//! fasterp - A high-performance FASTQ preprocessor
//!
//! This library provides FASTQ quality control and adapter trimming.
//!
//! # Features
//! - `native` (default): Full functionality with threading, file I/O, S3, HTTP
//! - `wasm`: WebAssembly support with buffer-based processing
#![allow(
    dead_code,
    unused_imports,
    unused_variables,
    unused_unsafe,
    private_interfaces
)]

#[cfg(feature = "python")]
use pyo3::prelude::*;

// Core modules - always available
pub mod adapter;
pub mod bloom;
pub mod dedup;
pub mod kmer;
pub mod overlap;
pub mod processor;
pub mod simd;
pub mod stats;
pub mod trimming;
pub mod umi;
pub mod util;

// Native-only modules (require threading, file I/O, etc.)
#[cfg(feature = "native")]
pub mod html;
#[cfg(feature = "native")]
pub mod io;
#[cfg(feature = "native")]
pub mod pipeline;
#[cfg(feature = "native")]
pub mod split;

// WASM module
#[cfg(feature = "wasm")]
pub mod wasm;

pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

#[cfg(feature = "native")]
#[derive(Debug, Clone)]
pub struct NativeProcessResult {
    pub passed_reads: usize,
    pub failed_reads: usize,
    pub json_report: String,
}

#[cfg(feature = "native")]
pub fn process_paths(
    input: &str,
    output: &str,
    input2: Option<&str>,
    output2: Option<&str>,
    threads: usize,
) -> anyhow::Result<NativeProcessResult> {
    use crate::adapter::AdapterConfig;
    use crate::io::open_input;
    use crate::pipeline::{
        Batch, PairedBatch, PairedWorkerResult, WorkerResult, merger_thread,
        paired_merger_thread, paired_producer_thread_with_paths, paired_worker_thread,
        producer_thread, worker_thread,
    };
    use crate::processor::{process_fastq_stream, process_paired_fastq_stream};
    use crate::split::{SplitConfig, SplitWriter};
    use crate::trimming::TrimmingConfig;
    use crossbeam_channel::bounded;
    use serde_json::json;
    use std::thread;

    let threads = threads.max(1);
    let compression = 6;
    let parallel_compression = false;
    let split_config = SplitConfig::default();
    let adapter_config = AdapterConfig::new();
    let trimming_config = TrimmingConfig {
        enable_trim_front: false,
        enable_trim_tail: true,
        cut_mean_quality: 20,
        cut_window_size: 4,
        trim_front_bases: 0,
        trim_tail_bases: 0,
        max_len: 0,
        enable_poly_g: true,
        enable_poly_x: false,
        poly_min_len: 10,
        adapter_config: adapter_config.clone(),
    };

    if let (Some(input2), Some(output2)) = (input2, output2) {
        let mut trimming_config_r2 = trimming_config.clone();
        trimming_config_r2.adapter_config.adapter_seq = adapter_config.adapter_seq_r2.clone();
        let acc = if threads == 1 {
            let reader1 = open_input(input)?;
            let reader2 = open_input(input2)?;
            let mut writer1 =
                SplitWriter::new(output, split_config.clone(), compression, parallel_compression)?;
            let mut writer2 =
                SplitWriter::new(output2, split_config.clone(), compression, parallel_compression)?;
            let acc = process_paired_fastq_stream(
                reader1,
                reader2,
                &mut writer1,
                &mut writer2,
                None::<&mut Vec<u8>>,
                false,
                false,
                15,
                5,
                15,
                40,
                0,
                false,
                30,
                &trimming_config,
                &trimming_config_r2,
                None,
                None,
                None,
                None,
            )?;
            writer1.finish()?;
            writer2.finish()?;
            acc
        } else {
            let backlog = threads + 1;
            let batch_bytes = 32 * 1024 * 1024;
            let (batch_tx, batch_rx) = bounded::<Option<PairedBatch>>(backlog);
            let (result_tx, result_rx) = bounded::<Option<PairedWorkerResult>>(backlog);
            let producer = thread::spawn({
                let input1 = input.to_string();
                let input2 = input2.to_string();
                move || {
                    paired_producer_thread_with_paths(
                        input1, input2, batch_bytes, batch_tx, threads,
                    )
                }
            });
            let mut workers = Vec::new();
            for _ in 0..threads {
                let batch_rx = batch_rx.clone();
                let result_tx = result_tx.clone();
                let trim1 = trimming_config.clone();
                let trim2 = trimming_config_r2.clone();
                workers.push(thread::spawn(move || {
                    paired_worker_thread(
                        batch_rx, result_tx, 15, 5, 15, 40, 0, false, 30, false, trim1, trim2,
                        None, None, None, None,
                    );
                }));
            }
            drop(batch_rx);
            drop(result_tx);
            let merger = thread::spawn({
                let output1 = output.to_string();
                let output2 = output2.to_string();
                let split_config = split_config.clone();
                move || {
                    paired_merger_thread(
                        result_rx,
                        output1,
                        output2,
                        threads,
                        Some(compression),
                        split_config,
                        parallel_compression,
                    )
                }
            });
            producer.join().expect("paired producer thread panicked")?;
            for worker in workers {
                worker.join().expect("paired worker thread panicked");
            }
            merger.join().expect("paired merger thread panicked")?
        };
        let passed_reads = acc.after_r1.total_reads;
        let failed_reads = acc.before_r1.total_reads.saturating_sub(passed_reads);
        return Ok(NativeProcessResult {
            passed_reads,
            failed_reads,
            json_report: json!({
                "mode": "paired",
                "threads": threads,
                "passed_reads": passed_reads,
                "failed_reads": failed_reads,
                "duplicated_reads": acc.duplicated,
                "summary": {
                    "before_filtering": { "total_reads": acc.before_r1.total_reads },
                    "after_filtering": { "total_reads": passed_reads },
                },
                "filtering_result": {
                    "passed_filter_reads": passed_reads,
                    "low_quality_reads": acc.low_quality * 2,
                    "low_complexity_reads": acc.low_complexity * 2,
                    "too_many_N_reads": acc.too_many_n * 2,
                    "too_short_reads": acc.too_short * 2,
                    "too_long_reads": 0,
                },
                "duplication": {
                    "rate": if acc.before_r1.total_reads == 0 {
                        0.0
                    } else {
                        acc.duplicated as f64 / acc.before_r1.total_reads as f64
                    }
                },
                "adapter_cutting": {
                    "adapter_trimmed_reads": acc.adapter_trimmed_reads,
                    "adapter_trimmed_bases": acc.adapter_trimmed_bases,
                },
                "adapter_trimmed_reads": acc.adapter_trimmed_reads,
                "adapter_trimmed_bases": acc.adapter_trimmed_bases,
            })
            .to_string(),
        });
    }

    let acc = if threads == 1 {
        let reader = open_input(input)?;
        let mut writer =
            SplitWriter::new(output, split_config.clone(), compression, parallel_compression)?;
        let acc = process_fastq_stream(
            reader,
            &mut writer,
            15,
            5,
            15,
            40,
            0,
            false,
            30,
            &trimming_config,
        )?;
        writer.finish()?;
        acc
    } else {
        let backlog = threads + 1;
        let batch_bytes = 32 * 1024 * 1024;
        let (batch_tx, batch_rx) = bounded::<Option<Batch>>(backlog);
        let (result_tx, result_rx) = bounded::<Option<WorkerResult>>(backlog);
        let producer = thread::spawn({
            let input = input.to_string();
            move || producer_thread(input, batch_bytes, batch_tx)
        });
        let mut workers = Vec::new();
        for _ in 0..threads {
            let batch_rx = batch_rx.clone();
            let result_tx = result_tx.clone();
            let trim = trimming_config.clone();
            workers.push(thread::spawn(move || {
                worker_thread(batch_rx, result_tx, 15, 5, 15, 40, 0, false, 30, false, trim);
            }));
        }
        drop(batch_rx);
        drop(result_tx);
        let merger = thread::spawn({
            let output = output.to_string();
            let split_config = split_config.clone();
            move || {
                merger_thread(
                    result_rx,
                    output,
                    threads,
                    Some(compression),
                    split_config,
                    parallel_compression,
                )
            }
        });
        producer.join().expect("single-end producer thread panicked")?;
        for worker in workers {
            worker.join().expect("single-end worker thread panicked");
        }
        merger.join().expect("single-end merger thread panicked")?
    };
    let passed_reads = acc.after.total_reads;
    let failed_reads = acc.before.total_reads.saturating_sub(passed_reads);
    Ok(NativeProcessResult {
        passed_reads,
        failed_reads,
        json_report: json!({
            "mode": "single",
            "threads": threads,
            "passed_reads": passed_reads,
            "failed_reads": failed_reads,
            "duplicated_reads": 0,
            "summary": {
                "before_filtering": { "total_reads": acc.before.total_reads },
                "after_filtering": { "total_reads": passed_reads },
            },
            "filtering_result": {
                "passed_filter_reads": passed_reads,
                "low_quality_reads": acc.low_quality,
                "low_complexity_reads": acc.low_complexity,
                "too_many_N_reads": acc.too_many_n,
                "too_short_reads": acc.too_short,
                "too_long_reads": 0,
            },
            "duplication": { "rate": 0.0 },
            "adapter_cutting": {
                "adapter_trimmed_reads": acc.adapter_trimmed_reads,
                "adapter_trimmed_bases": acc.adapter_trimmed_bases,
            },
            "adapter_trimmed_reads": acc.adapter_trimmed_reads,
            "adapter_trimmed_bases": acc.adapter_trimmed_bases,
        })
        .to_string(),
    })
}

#[cfg(feature = "python")]
#[pyclass]
#[derive(Clone)]
pub struct ProcessResult {
    #[pyo3(get)]
    pub passed_reads: usize,
    #[pyo3(get)]
    pub failed_reads: usize,
    #[pyo3(get)]
    pub total_bases_before: usize,
    #[pyo3(get)]
    pub total_bases_after: usize,
    #[pyo3(get)]
    pub q20_rate: f64,
    #[pyo3(get)]
    pub q30_rate: f64,
    #[pyo3(get)]
    pub gc_content: f64,
    #[pyo3(get)]
    pub duplication_rate: f64,
    #[pyo3(get)]
    pub adapter_trimmed_reads: usize,
    #[pyo3(get)]
    pub adapter_trimmed_bases: usize,
    #[pyo3(get)]
    pub json_report: String,
}

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (
    input,
    output,
    input2=None,
    output2=None,
    json=None,
    threads=1,
    min_length=15,
    max_length=0,
    average_qual=0,
    n_base_limit=5,
    qualified_quality_phred=15,
    unqualified_percent_limit=40,
    low_complexity_filter=false,
    complexity_threshold=30,
    trim_front=0,
    trim_tail=0,
    cut_front=false,
    cut_tail=true,
    cut_mean_quality=20,
    cut_window_size=4,
    trim_poly_g=true,
    trim_poly_x=false,
    poly_g_min_len=10,
    adapter_sequence=None,
    adapter_sequence_r2=None,
    disable_adapter_trimming=false,
    compression_level=4
))]
#[allow(clippy::too_many_arguments)]
pub fn process(
    input: String,
    output: String,
    input2: Option<String>,
    output2: Option<String>,
    json: Option<String>,
    threads: usize,
    min_length: usize,
    max_length: usize,
    average_qual: u8,
    n_base_limit: usize,
    qualified_quality_phred: u8,
    unqualified_percent_limit: usize,
    low_complexity_filter: bool,
    complexity_threshold: usize,
    trim_front: usize,
    trim_tail: usize,
    cut_front: bool,
    cut_tail: bool,
    cut_mean_quality: u8,
    cut_window_size: usize,
    trim_poly_g: bool,
    trim_poly_x: bool,
    poly_g_min_len: usize,
    adapter_sequence: Option<String>,
    adapter_sequence_r2: Option<String>,
    disable_adapter_trimming: bool,
    compression_level: u32,
) -> PyResult<ProcessResult> {
    use crate::adapter::AdapterConfig;
    use crate::io::open_output;
    use crate::processor::{process_fastq_stream, process_paired_fastq_stream};
    use crate::stats::{
        AdapterCuttingStats, DetailedReadStats, DuplicationStats, FasterpReport, FilteringResult,
        Summary,
    };
    use crate::trimming::TrimmingConfig;

    // Create adapter configuration
    let mut adapter_config = AdapterConfig::new();
    if disable_adapter_trimming {
        adapter_config.detect_adapter_for_pe = false;
    } else {
        if let Some(ref seq) = adapter_sequence {
            adapter_config.adapter_seq = Some(seq.as_bytes().to_vec());
        }
        if let Some(ref seq) = adapter_sequence_r2 {
            adapter_config.adapter_seq_r2 = Some(seq.as_bytes().to_vec());
        }
    }

    // Create trimming configuration
    let trimming_config = TrimmingConfig {
        enable_trim_front: cut_front && cut_mean_quality > 0,
        enable_trim_tail: cut_tail && cut_mean_quality > 0,
        cut_mean_quality,
        cut_window_size,
        trim_front_bases: trim_front,
        trim_tail_bases: trim_tail,
        max_len: max_length,
        enable_poly_g: trim_poly_g,
        enable_poly_x: trim_poly_x,
        poly_min_len: poly_g_min_len,
        adapter_config: adapter_config.clone(),
    };

    // Check if paired-end
    let is_paired = input2.is_some() && output2.is_some();

    if is_paired {
        // Paired-end processing
        let input2_path = input2.unwrap();
        let output2_path = output2.unwrap();

        // Open input files
        let reader1 = crate::io::open_input(&input).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open input: {}", e))
        })?;
        let reader2 = crate::io::open_input(&input2_path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open input2: {}", e))
        })?;

        // Open output files
        let mut writer1 = open_output(&output, Some(compression_level), false).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create output: {}", e))
        })?;
        let mut writer2 = open_output(&output2_path, Some(compression_level), false).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create output2: {}", e))
        })?;

        // Create trimming config for R2
        let mut trimming_config_r2 = trimming_config.clone();
        trimming_config_r2.adapter_config.adapter_seq = adapter_config.adapter_seq_r2.clone();

        // Dummy writer for merged output (not used)
        let mut dummy_writer: Vec<u8> = Vec::new();

        // Process paired-end
        let pe_acc = process_paired_fastq_stream(
            reader1,
            reader2,
            &mut writer1,
            &mut writer2,
            None::<&mut Vec<u8>>, // merged_writer
            false, // merge_enabled
            false, // include_unmerged
            min_length,
            n_base_limit,
            qualified_quality_phred,
            unqualified_percent_limit,
            average_qual,
            low_complexity_filter,
            complexity_threshold,
            &trimming_config,
            &trimming_config_r2,
            None, // overlap_config
            None, // umi_config
            None, // dedup_config
            None, // dup_detector
        )
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Processing failed: {}", e))
        })?;

        writer1.finish().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to finish output: {}", e))
        })?;
        writer2.finish().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to finish output2: {}", e))
        })?;

        // Build report
        let before_stats_r1 = pe_acc.before_r1.to_read_stats();
        let after_stats_r1 = pe_acc.after_r1.to_read_stats();
        let before_stats_r2 = pe_acc.before_r2.to_read_stats();
        let after_stats_r2 = pe_acc.after_r2.to_read_stats();

        let kmer_map_r1 = pe_acc.kmer_table_to_map_r1();
        let kmer_map_r2 = pe_acc.kmer_table_to_map_r2();
        let dup_rate_r1 = stats::calculate_duplication_rate(&kmer_map_r1);
        let dup_rate_r2 = stats::calculate_duplication_rate(&kmer_map_r2);
        let duplication_rate = (dup_rate_r1 + dup_rate_r2) / 2.0;

        let total_before = before_stats_r1.total_reads + before_stats_r2.total_reads;
        let total_after = after_stats_r1.total_reads + after_stats_r2.total_reads;
        let failed = total_before - total_after;

        let report = FasterpReport {
            summary: Summary {
                fastp_version: env!("CARGO_PKG_VERSION").to_string(),
                sequencing: format!("paired end ({} cycles + {} cycles)", pe_acc.max_cycle_r1, pe_acc.max_cycle_r2),
                before_filtering: before_stats_r1.clone(),
                after_filtering: after_stats_r1.clone(),
                read1_before_filtering: Some(before_stats_r1.clone()),
                read2_before_filtering: Some(before_stats_r2.clone()),
                read1_after_filtering: Some(after_stats_r1.clone()),
                read2_after_filtering: Some(after_stats_r2.clone()),
            },
            filtering_result: FilteringResult {
                passed_filter_reads: pe_acc.after_r1.total_reads,
                low_quality_reads: pe_acc.low_quality,
                low_complexity_reads: pe_acc.low_complexity,
                too_many_n_reads: pe_acc.too_many_n,
                too_short_reads: pe_acc.too_short,
                too_long_reads: 0,
            },
            read1_before_filtering: DetailedReadStats {
                total_reads: before_stats_r1.total_reads,
                total_bases: before_stats_r1.total_bases,
                q20_bases: before_stats_r1.q20_bases,
                q30_bases: before_stats_r1.q30_bases,
                quality_curves: pe_acc.pos_r1.to_quality_curves(),
                content_curves: pe_acc.pos_r1.to_content_curves(),
                qual_hist: pe_acc.pos_r1.to_qual_hist(),
                kmer_count: kmer_map_r1,
            },
            read2_before_filtering: Some(DetailedReadStats {
                total_reads: before_stats_r2.total_reads,
                total_bases: before_stats_r2.total_bases,
                q20_bases: before_stats_r2.q20_bases,
                q30_bases: before_stats_r2.q30_bases,
                quality_curves: pe_acc.pos_r2.to_quality_curves(),
                content_curves: pe_acc.pos_r2.to_content_curves(),
                qual_hist: pe_acc.pos_r2.to_qual_hist(),
                kmer_count: kmer_map_r2,
            }),
            read1_after_filtering: Some(DetailedReadStats {
                total_reads: after_stats_r1.total_reads,
                total_bases: after_stats_r1.total_bases,
                q20_bases: after_stats_r1.q20_bases,
                q30_bases: after_stats_r1.q30_bases,
                quality_curves: pe_acc.pos_r1_after.to_quality_curves(),
                content_curves: pe_acc.pos_r1_after.to_content_curves(),
                qual_hist: pe_acc.pos_r1_after.to_qual_hist(),
                kmer_count: pe_acc.kmer_table_to_map_r1_after(),
            }),
            read2_after_filtering: Some(DetailedReadStats {
                total_reads: after_stats_r2.total_reads,
                total_bases: after_stats_r2.total_bases,
                q20_bases: after_stats_r2.q20_bases,
                q30_bases: after_stats_r2.q30_bases,
                quality_curves: pe_acc.pos_r2_after.to_quality_curves(),
                content_curves: pe_acc.pos_r2_after.to_content_curves(),
                qual_hist: pe_acc.pos_r2_after.to_qual_hist(),
                kmer_count: pe_acc.kmer_table_to_map_r2_after(),
            }),
            duplication: Some(DuplicationStats { rate: duplication_rate }),
            adapter_cutting: if pe_acc.adapter_trimmed_reads > 0 {
                Some(AdapterCuttingStats {
                    adapter_trimmed_reads: pe_acc.adapter_trimmed_reads,
                    adapter_trimmed_bases: pe_acc.adapter_trimmed_bases,
                })
            } else {
                None
            },
            insert_size: None,
        };

        // Write JSON report if requested
        let json_report = serde_json::to_string_pretty(&report).unwrap_or_default();
        if let Some(json_path) = json {
            std::fs::write(&json_path, &json_report).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to write JSON: {}", e))
            })?;
        }

        Ok(ProcessResult {
            passed_reads: total_after,
            failed_reads: failed,
            total_bases_before: before_stats_r1.total_bases + before_stats_r2.total_bases,
            total_bases_after: after_stats_r1.total_bases + after_stats_r2.total_bases,
            q20_rate: after_stats_r1.q20_rate,
            q30_rate: after_stats_r1.q30_rate,
            gc_content: after_stats_r1.gc_content,
            duplication_rate,
            adapter_trimmed_reads: pe_acc.adapter_trimmed_reads,
            adapter_trimmed_bases: pe_acc.adapter_trimmed_bases,
            json_report,
        })
    } else {
        // Single-end processing
        let reader = crate::io::open_input(&input).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open input: {}", e))
        })?;

        // Open output file
        let mut writer = open_output(&output, Some(compression_level), false).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create output: {}", e))
        })?;

        // Process single-end
        let acc = process_fastq_stream(
            reader,
            &mut writer,
            min_length,
            n_base_limit,
            qualified_quality_phred,
            unqualified_percent_limit,
            average_qual,
            low_complexity_filter,
            complexity_threshold,
            &trimming_config,
        )
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Processing failed: {}", e))
        })?;

        writer.finish().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to finish output: {}", e))
        })?;

        // Build report
        let before_stats = acc.before.to_read_stats();
        let after_stats = acc.after.to_read_stats();
        let quality_curves = acc.pos.to_quality_curves();
        let content_curves = acc.pos.to_content_curves();
        let qual_hist = acc.pos.to_qual_hist();
        let quality_curves_after = acc.pos_after.to_quality_curves();
        let content_curves_after = acc.pos_after.to_content_curves();
        let qual_hist_after = acc.pos_after.to_qual_hist();
        let kmer_map = acc.kmer_table_to_map();
        let kmer_map_after = acc.kmer_table_after_to_map();
        let duplication_rate = stats::calculate_duplication_rate(&kmer_map);

        let report = FasterpReport {
            summary: Summary {
                fastp_version: env!("CARGO_PKG_VERSION").to_string(),
                sequencing: format!("single end ({} cycles)", acc.max_cycle),
                before_filtering: before_stats.clone(),
                after_filtering: after_stats.clone(),
                read1_before_filtering: None,
                read2_before_filtering: None,
                read1_after_filtering: None,
                read2_after_filtering: None,
            },
            filtering_result: FilteringResult {
                passed_filter_reads: acc.after.total_reads,
                low_quality_reads: acc.low_quality,
                low_complexity_reads: acc.low_complexity,
                too_many_n_reads: acc.too_many_n,
                too_short_reads: acc.too_short,
                too_long_reads: 0,
            },
            read1_before_filtering: DetailedReadStats {
                total_reads: before_stats.total_reads,
                total_bases: before_stats.total_bases,
                q20_bases: before_stats.q20_bases,
                q30_bases: before_stats.q30_bases,
                quality_curves,
                content_curves,
                qual_hist,
                kmer_count: kmer_map,
            },
            read2_before_filtering: None,
            read1_after_filtering: Some(DetailedReadStats {
                total_reads: after_stats.total_reads,
                total_bases: after_stats.total_bases,
                q20_bases: after_stats.q20_bases,
                q30_bases: after_stats.q30_bases,
                quality_curves: quality_curves_after,
                content_curves: content_curves_after,
                qual_hist: qual_hist_after,
                kmer_count: kmer_map_after,
            }),
            read2_after_filtering: None,
            duplication: Some(DuplicationStats { rate: duplication_rate }),
            adapter_cutting: if acc.adapter_trimmed_reads > 0 {
                Some(AdapterCuttingStats {
                    adapter_trimmed_reads: acc.adapter_trimmed_reads,
                    adapter_trimmed_bases: acc.adapter_trimmed_bases,
                })
            } else {
                None
            },
            insert_size: None,
        };

        // Write JSON report if requested
        let json_report = serde_json::to_string_pretty(&report).unwrap_or_default();
        if let Some(json_path) = json {
            std::fs::write(&json_path, &json_report).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to write JSON: {}", e))
            })?;
        }

        Ok(ProcessResult {
            passed_reads: acc.after.total_reads,
            failed_reads: acc.before.total_reads - acc.after.total_reads,
            total_bases_before: before_stats.total_bases,
            total_bases_after: after_stats.total_bases,
            q20_rate: after_stats.q20_rate,
            q30_rate: after_stats.q30_rate,
            gc_content: after_stats.gc_content,
            duplication_rate,
            adapter_trimmed_reads: acc.adapter_trimmed_reads,
            adapter_trimmed_bases: acc.adapter_trimmed_bases,
            json_report,
        })
    }
}

#[cfg(feature = "python")]
#[pymodule]
fn fasterp(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ProcessResult>()?;
    m.add_function(wrap_pyfunction!(process, m)?)?;
    Ok(())
}
