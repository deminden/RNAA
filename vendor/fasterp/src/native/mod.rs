//! Bounded native FASTQ processing and execution observations.
mod batches;
mod dispatch;
#[cfg(test)]
mod tests;
mod window;
pub(crate) use window::BatchPermit;

#[derive(Debug, Clone, Copy, serde::Serialize)]
pub struct NativeProcessOptions {
    pub threads: usize,
    pub batch_bytes: usize,
    pub max_in_flight_batches: usize,
    pub max_record_bytes: usize,
    /// Optional total decoded byte limits for R1 and R2. Extra input is an error.
    pub input_byte_limits: Option<[u64; 2]>,
    pub min_record_bytes: usize,
    pub max_input_records: Option<u64>,
}
impl NativeProcessOptions {
    pub fn for_threads(threads: usize) -> Self {
        let threads = threads.max(1);
        Self {
            threads,
            batch_bytes: 32 * 1024 * 1024,
            max_in_flight_batches: threads.saturating_add(1),
            max_record_bytes: 4 * 1024 * 1024,
            input_byte_limits: None,
            min_record_bytes: 0,
            max_input_records: None,
        }
    }
    fn validate(self) -> anyhow::Result<()> {
        if self.threads == 0
            || self.batch_bytes == 0
            || self.max_in_flight_batches == 0
            || self.max_record_bytes < 4
            || self.min_record_bytes > self.max_record_bytes
        {
            anyhow::bail!(
                "native preprocessing threads, batch/window sizes and record limit must be positive"
            );
        }
        if self
            .batch_bytes
            .checked_add(self.max_record_bytes)
            .and_then(|bytes| bytes.checked_add(1))
            .and_then(|bytes| bytes.checked_mul(self.max_in_flight_batches))
            .is_none_or(|bytes| bytes > isize::MAX as usize)
        {
            anyhow::bail!("native preprocessing window exceeds addressable storage");
        }
        Ok(())
    }
}

pub fn process_paths(
    input: &str,
    output: &str,
    input2: Option<&str>,
    output2: Option<&str>,
    threads: usize,
) -> anyhow::Result<NativeProcessResult> {
    process_paths_with_options(
        input,
        output,
        input2,
        output2,
        NativeProcessOptions::for_threads(threads),
    )
}
pub fn process_paths_with_options(
    input: &str,
    output: &str,
    input2: Option<&str>,
    output2: Option<&str>,
    options: NativeProcessOptions,
) -> anyhow::Result<NativeProcessResult> {
    options.validate()?;
    if input2.is_some() != output2.is_some() {
        anyhow::bail!("paired native processing requires both input2 and output2");
    }
    let window = window::BatchWindow::new(options.max_in_flight_batches);
    let mut result = run_paths(input, output, input2, output2, options, &window)?;
    let stats = window.stats()?;
    if stats.active != 0 {
        anyhow::bail!("native preprocessing retained batch permits after its workers ended");
    }
    let mut report: serde_json::Value = serde_json::from_str(&result.json_report)?;
    report["execution"] = serde_json::json!({"options":options,"batch_window":stats});
    result.json_report = serde_json::to_string(&report)?;
    Ok(result)
}

#[cfg(feature = "native")]
#[derive(Debug, Clone)]
pub struct NativeProcessResult {
    pub passed_reads: usize,
    pub failed_reads: usize,
    pub json_report: String,
}

#[cfg(feature = "native")]
fn run_paths(
    input: &str,
    output: &str,
    input2: Option<&str>,
    output2: Option<&str>,
    options: NativeProcessOptions,
    window: &std::sync::Arc<window::BatchWindow>,
) -> anyhow::Result<NativeProcessResult> {
    use crate::processor::{process_fastq_stream, process_paired_fastq_stream};
    use crate::split::{SplitConfig, SplitWriter};
    use serde_json::json;

    let threads = options.threads;
    let compression = 6;
    let parallel_compression = false;
    let split_config = SplitConfig::default();
    let trimming_config = default_trimming();

    if let (Some(input2), Some(output2)) = (input2, output2) {
        let mut trimming_config_r2 = trimming_config.clone();
        trimming_config_r2.adapter_config.adapter_seq =
            trimming_config.adapter_config.adapter_seq_r2.clone();
        let acc = if threads == 1 {
            let reader1 = batches::input(input, options, 0)?;
            let reader2 = batches::input(input2, options, 1)?;
            let mut writer1 = SplitWriter::new(
                output,
                split_config.clone(),
                compression,
                parallel_compression,
            )?;
            let mut writer2 = SplitWriter::new(
                output2,
                split_config.clone(),
                compression,
                parallel_compression,
            )?;
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
            dispatch::paired(
                input,
                input2,
                output,
                output2,
                options,
                trimming_config,
                trimming_config_r2,
                window,
            )?
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
        let reader = batches::input(input, options, 0)?;
        let mut writer = SplitWriter::new(
            output,
            split_config.clone(),
            compression,
            parallel_compression,
        )?;
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
        dispatch::single(input, output, options, trimming_config, window)?
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

fn default_trimming() -> crate::trimming::TrimmingConfig {
    crate::trimming::TrimmingConfig {
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
        adapter_config: crate::adapter::AdapterConfig::new(),
    }
}
