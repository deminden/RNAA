//! One lazy, immutable counting index for a compatible group of quant workers.
use super::QuantIndexPlan;
mod admission;
pub(super) use admission::ResidentIndex;
use anyhow::{Result, bail};
use kallistors::pseudoalign::BifrostIndex;
use rnaa_core::{ManifestArtifact, ReferenceBundle};
use std::{
    path::Path,
    sync::{
        Arc, Mutex, OnceLock,
        atomic::{AtomicUsize, Ordering},
    },
};

static NEXT_GROUP: AtomicUsize = AtomicUsize::new(1);

/// Execution-only observations; these do not participate in scientific reuse.
#[derive(Debug, Clone, serde::Serialize)]
pub struct QuantIndexStats {
    pub group_id: String,
    pub loader: &'static str,
    pub engine_version: &'static str,
    pub index_sha256: String,
    pub index_bytes: u64,
    pub load_attempts: usize,
    pub successful_loads: usize,
    pub layout_plan: Option<QuantIndexPlan>,
    pub evictions: usize,
}

/// A group owns exactly one frozen index identity. Construction does not load it.
/// All mapping/EM state remains local to each quantification call; only the
/// engine's read-only graph, transcript metadata and lookup tables are shared.
pub struct QuantIndexGroup {
    input: ManifestArtifact,
    id: String,
    index: Mutex<Option<Arc<ResidentIndex>>>,
    attempts: AtomicUsize,
    loads: AtomicUsize,
    plan: OnceLock<QuantIndexPlan>,
    evictions: AtomicUsize,
}

impl QuantIndexGroup {
    pub fn new(input: ManifestArtifact) -> Result<Self> {
        if input.kind != "REFERENCE_INDEX"
            || input.path.is_empty()
            || input.checksum_type != "sha256"
            || input.checksum.len() != 64
            || !input.checksum.bytes().all(|byte| byte.is_ascii_hexdigit())
            || input.bytes == 0
        {
            bail!("quant worker group requires a frozen SHA-256 reference index descriptor");
        }
        Ok(Self {
            input,
            id: format!(
                "{}-{}",
                std::process::id(),
                NEXT_GROUP.fetch_add(1, Ordering::Relaxed)
            ),
            index: Mutex::new(None),
            attempts: AtomicUsize::new(0),
            loads: AtomicUsize::new(0),
            plan: OnceLock::new(),
            evictions: AtomicUsize::new(0),
        })
    }

    pub fn stats(&self) -> QuantIndexStats {
        QuantIndexStats {
            group_id: self.id.clone(),
            loader: "kallistors-counting-v1",
            engine_version: crate::engines::KALLISTORS,
            index_sha256: self.input.checksum.clone(),
            index_bytes: self.input.bytes,
            load_attempts: self.attempts.load(Ordering::Relaxed),
            successful_loads: self.loads.load(Ordering::Relaxed),
            layout_plan: self.plan.get().cloned(),
            evictions: self.evictions.load(Ordering::Relaxed),
        }
    }

    pub fn preflight(
        &self,
        reference: &ReferenceBundle,
        check: &impl Fn() -> Result<()>,
    ) -> Result<QuantIndexPlan> {
        check()?;
        self.verify(reference)?;
        if let Some(plan) = self.plan.get() {
            check()?;
            return Ok(plan.clone());
        }
        let plan = QuantIndexPlan::inspect(&reference.kallisto_index_path, check)
            .map_err(|error| anyhow::anyhow!("quant index preflight failed: {error:#}"))?;
        self.verify(reference)?;
        check()?;
        let _ = self.plan.set(plan);
        Ok(self.plan.get().expect("installed layout plan").clone())
    }

    pub(super) fn verify(&self, reference: &ReferenceBundle) -> Result<()> {
        if reference.kallisto_index_path != Path::new(&self.input.path)
            || !self.input.matches_file()?
        {
            bail!(
                "quant worker group's frozen reference index changed; retry with the current reference"
            );
        }
        Ok(())
    }

    pub(super) fn load(
        &self,
        reference: &ReferenceBundle,
        threads: usize,
    ) -> Result<Arc<ResidentIndex>> {
        self.verify(reference)?;
        let mut slot = self
            .index
            .lock()
            .map_err(|_| anyhow::anyhow!("quant worker group's index loader panicked"))?;
        if let Some(index) = slot.as_ref() {
            if index.reservation.is_some() {
                bail!("an admitted quant group requires owned run admission");
            }
            return Ok(Arc::clone(index));
        }
        let index = Arc::new(ResidentIndex {
            index: self.build(reference, threads)?,
            reservation: None,
            pool_root: None,
        });
        *slot = Some(Arc::clone(&index));
        self.loads.fetch_add(1, Ordering::Relaxed);
        Ok(index)
    }
    fn build(&self, reference: &ReferenceBundle, threads: usize) -> Result<BifrostIndex> {
        self.attempts.fetch_add(1, Ordering::Relaxed);
        if self.plan.get().is_none() {
            self.preflight(reference, &|| Ok(()))?;
        }
        let index = kallistors::pseudoalign::build_bifrost_index_quant_threaded(
            &reference.kallisto_index_path,
            false,
            threads.max(1),
        )
        .map_err(|error| anyhow::anyhow!("failed to load kallistors index: {error}"))?;
        // Loading a file does not establish that it is still the frozen input.
        self.verify(reference)?;
        Ok(index)
    }
}

#[cfg(test)]
mod tests;
