//! Structural payload estimates and capped metadata writes for native defaults.
use crate::fastq::{FastqLimits, FastqPlan};
use anyhow::{Context, Result};
use fasterp::NativeProcessOptions;
use rnaa_core::resources::{DiskRequest, ResourceAdmissionError, ResourceRequest};
use std::path::Path;
const MIB: u64 = 1024 * 1024;

#[derive(Debug, serde::Serialize)]
pub(super) struct Estimate {
    pub memory_bytes: u64,
    pub staging_bytes: u64,
    pub read_output_bytes: u64,
    pub metadata_bytes: u64,
    pub native_threads: usize,
    pub model: &'static str,
}
pub(super) fn preflight(limits: FastqLimits) -> Result<ResourceRequest> {
    crate::fastq::scan_request(limits)
}
impl Estimate {
    pub fn for_plan(plan: &FastqPlan, options: NativeProcessOptions) -> Result<Self> {
        let s = plan.summary();
        let mates = if s.paired { 2_u128 } else { 1 };
        let window = if options.threads > 1 {
            options.max_in_flight_batches as u128
        } else {
            1
        };
        let record = options.max_record_bytes as u128;
        let total = s.normalized_bytes()? as u128;
        let batch = (options.batch_bytes as u128 + record).min(total);
        // CRLF may remove four bytes per record. Use the smaller normalized
        // denominator conservatively; no payload estimate uses compressed size.
        let min_record = s.min_record_bytes.saturating_sub(4).max(1) as u128;
        let records = (batch / min_record + 1).min(s.fragments as u128);
        let cycles = (options.max_record_bytes as u128).max(150);
        // Vec spare capacities, original/trimmed buffers, record offsets and
        // result pieces; before/after position arrays plus merger; JSON trees
        // and serialization; decoder/writer buffers, stacks and fixed overhead.
        let payload = mates * window * (batch * 8 + records * 256);
        let positions = cycles * mates * 320 * (window + 1);
        let descriptors = serde_json::to_vec(plan.inputs())?.len() as u128;
        let metadata = 8 * MIB as u128 + cycles * mates * 4096 + descriptors * 4;
        let memory = payload
            + positions
            + metadata * 4
            + 64 * MIB as u128
            + (options.threads as u128 + 2) * 2 * MIB as u128;
        Ok(Self {
            memory_bytes: u64::try_from(memory)
                .context("preprocessing memory estimate overflows")?,
            staging_bytes: s.staging_bytes()?,
            read_output_bytes: u64::try_from(total)?,
            metadata_bytes: u64::try_from(metadata)
                .context("preprocessing metadata estimate overflows")?,
            native_threads: options.threads,
            model: "fasterp-native-defaults-v1",
        })
    }
    pub fn request(&self, staging: &Path, output: &Path) -> Result<ResourceRequest> {
        let mut disks = Vec::new();
        if self.staging_bytes > 0 {
            disks.push(DiskRequest {
                path: staging.into(),
                bytes: self.staging_bytes,
            });
        }
        disks.push(DiskRequest {
            path: output.into(),
            bytes: self
                .read_output_bytes
                .checked_add(self.metadata_bytes)
                .context("preprocessing output estimate overflows")?,
        });
        Ok(ResourceRequest {
            cpu_threads: self.native_threads as u64,
            memory_bytes: self.memory_bytes,
            disks,
        })
    }
}
pub(super) fn write_metadata(
    path: &Path,
    value: &impl serde::Serialize,
    remaining: &mut u64,
) -> Result<()> {
    let content = serde_json::to_vec_pretty(value)?;
    let bytes = content.len() as u64;
    if bytes > *remaining {
        return Err(ResourceAdmissionError::new(anyhow::anyhow!(
            "preprocessing metadata exceeds its reserved byte promise"
        ))
        .into());
    }
    rnaa_core::util::publish_file(path, |staged| Ok(std::fs::write(staged, &content)?))?;
    *remaining -= bytes;
    Ok(())
}

#[cfg(test)]
mod tests;
