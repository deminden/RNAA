//! Decoded request bounds for one mapping/EM invocation plus its shared graph.
use super::QuantIndexPlan;
use crate::fastq::{FastqLimits, FastqPlan};
use anyhow::{Context, Result, bail};
use rnaa_core::resources::{DiskRequest, ResourceRequest};
use std::path::Path;
const MIB: u128 = 1024 * 1024;

#[derive(Clone, Debug, serde::Serialize)]
pub(super) struct QuantResources {
    pub model: &'static str,
    pub index_loading_bytes: u64,
    pub index_loading_reservation_bytes: u64,
    pub loading_setup_bytes: u64,
    pub index_resident_bytes: u64,
    pub run_memory_bytes: u64,
    pub match_cache_bytes: u64,
    pub runtime_cache_bytes: u64,
    pub ec_classes_upper: u64,
    pub ec_members_upper: u64,
    pub staging_bytes: u64,
    pub output_bytes: u64,
    pub threads: usize,
}
pub(super) fn scan_request(limits: FastqLimits) -> Result<ResourceRequest> {
    crate::fastq::scan_request(limits)
}
fn subset_bound(targets: u64, card: u64, fragments: u64) -> u64 {
    let (mut choose, mut sum) = (1_u128, 0_u128);
    for k in 1..=targets.min(card) {
        choose = choose * u128::from(targets - k + 1) / u128::from(k);
        sum += choose;
        if sum >= u128::from(fragments) {
            return fragments;
        }
    }
    sum as u64
}
impl QuantResources {
    pub fn estimate(
        index: &QuantIndexPlan,
        fastq: &FastqPlan,
        threads: usize,
        cache_bytes: u64,
    ) -> Result<Self> {
        let s = fastq.summary();
        if s.fragments > u64::from(u32::MAX) {
            bail!("quant run exceeds the engine's 32-bit EC counter limit");
        }
        let workers = threads.max(1) as u128;
        let targets = u128::from(index.transcripts);
        let card = if index.has_shade_names {
            index.transcripts
        } else {
            index.max_ec_cardinality.min(index.transcripts)
        };
        let mut classes = subset_bound(index.transcripts, card, s.fragments);
        if !index.has_shade_names {
            classes = classes.min(
                index
                    .transcripts
                    .saturating_add(index.non_single_subset_classes_upper),
            );
        }
        let mut members = u128::from(classes) * u128::from(card);
        if !index.has_shade_names {
            members = members.min(
                u128::from(index.transcripts) + u128::from(index.non_single_subset_members_upper),
            );
        }
        let mates = if s.paired { 2_u128 } else { 1 };
        let batch = if s.paired { 6000_u128 } else { 4096 };
        let sequence = s.max_sequence_bytes as u128;
        if sequence * batch > u128::from(u32::MAX) {
            bail!(
                "quant FASTQ batch exceeds the engine's 32-bit sequence span; shorter read batches require engine support"
            );
        }
        let slots = if threads > 1 { workers + 5 } else { 2 };
        let batch_records = u128::from(s.fragments).min(batch);
        let per_batch = mates * ((batch * 300).max(batch_records * sequence * 2) + batch * 32);
        let per_job_batches = if s.paired && threads > 1 {
            // The paired mapper groups parser batches into an 8 MiB sampling
            // job. Short reads can retain several preallocated parser buffers.
            let minimum = mates * (s.min_sequence_bytes.max(1) as u128 + 1);
            let records = u128::from(s.fragments).min((8 * MIB).div_ceil(minimum));
            records.div_ceil(batch) + 1
        } else {
            1
        };
        let packed = slots * per_job_batches * per_batch;
        // Each input contributes to at most one worker's key list. Pool keys and
        // merged ECs also coexist; EM keeps sparse edge/weight/transpose arrays.
        let worker_classes = (u128::from(classes) * workers).min(u128::from(s.fragments));
        let ec = (worker_classes + u128::from(classes) * 3) * 256 + members * 160;
        // Small budgets restore the engine's default 8 MiB TLS table on scope
        // exit. Dedicated mapper-thread teardown is part of the guarded lifetime.
        let runtime_cache = u128::from(cache_bytes).max(workers * 8 * MIB);
        // The nominal aggregate cache is 70 MiB; replacement entries can coexist
        // transiently on up to 64 workers before the engine enforces its budget.
        let pairs = if card >= 128 { 72 * MIB } else { 0 };
        let scratch = sequence * workers * (u128::from(card) * 16 + 128);
        let inventory = serde_json::to_vec(fastq.inputs())?.len() as u128 * 8;
        let memory = 64 * MIB
            + runtime_cache
            + workers * 32 * MIB
            + packed
            + ec
            + targets * 1024
            + scratch
            + pairs;
        let memory = memory + inventory;
        let loading_setup = 64 * MIB + inventory;
        let index_loading = u128::from(index.load_memory_bytes(threads)?);
        let index_resident = u128::from(index.resident_estimate_bytes);
        // The first run cannot allocate mapping state until loading has ended.
        // Admit the larger phase plus its retained setup, not the sum of both
        // peaks. The run's unused allowance supports the load temporarily.
        let loading_claim = (index_loading + loading_setup).max(index_resident + memory) - memory;
        let output = 64 * MIB
            + u128::from(index.transcript_name_bytes) * 12
            + targets * 4096
            + (serde_json::to_vec(fastq.inputs())?.len() as u128) * 8;
        Ok(Self {
            model: "kallistors-mapping-em-v1",
            index_loading_bytes: index.load_memory_bytes(threads)?,
            index_loading_reservation_bytes: u64::try_from(loading_claim)
                .context("quant loading reservation overflows")?,
            loading_setup_bytes: u64::try_from(loading_setup)?,
            index_resident_bytes: index.resident_estimate_bytes,
            run_memory_bytes: u64::try_from(memory)
                .context("quant run memory estimate overflows")?,
            match_cache_bytes: cache_bytes,
            runtime_cache_bytes: u64::try_from(runtime_cache)?,
            ec_classes_upper: classes,
            ec_members_upper: u64::try_from(members)
                .context("quant EC member estimate overflows")?,
            staging_bytes: s.staging_bytes()?,
            output_bytes: u64::try_from(output).context("quant output estimate overflows")?,
            threads: threads.max(1),
        })
    }
    pub fn index_request(&self, loading: bool) -> ResourceRequest {
        ResourceRequest {
            memory_bytes: if loading {
                self.index_loading_reservation_bytes
            } else {
                self.index_resident_bytes
            },
            ..Default::default()
        }
    }
    pub fn run_request(&self, staging: &Path, output: &Path) -> ResourceRequest {
        let mut disks = vec![DiskRequest {
            path: output.into(),
            bytes: self.output_bytes,
        }];
        if self.staging_bytes > 0 {
            disks.push(DiskRequest {
                path: staging.into(),
                bytes: self.staging_bytes,
            });
        }
        ResourceRequest {
            cpu_threads: self.threads as u64,
            memory_bytes: self.run_memory_bytes,
            disks,
        }
    }
}

#[cfg(test)]
mod tests;
