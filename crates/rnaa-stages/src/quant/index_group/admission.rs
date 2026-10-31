//! Resident graph ownership, atomic first-run admission and idle eviction.
use super::*;
use crate::quant::resources::QuantResources;
use anyhow::Context;
use rnaa_core::resources::{
    Admission, AdmissionWait, BundleAdmission, ExecutionResources, ReservationSnapshot,
    ResourceAdmissionError, ResourceRequest, ResourceReservation,
};
use std::{
    ops::Deref,
    path::PathBuf,
    sync::TryLockError,
    time::{Duration, Instant},
};

pub(in crate::quant) struct ResidentIndex {
    pub index: BifrostIndex,
    // Field order: free the graph before releasing its memory promise.
    pub reservation: Option<ResourceReservation>,
    pub pool_root: Option<PathBuf>,
}
impl Deref for ResidentIndex {
    type Target = BifrostIndex;
    fn deref(&self) -> &BifrostIndex {
        &self.index
    }
}

pub(in crate::quant) struct RunLease {
    pub index: Arc<ResidentIndex>,
    pub reservation: Option<ResourceReservation>,
}
impl Deref for RunLease {
    type Target = BifrostIndex;
    fn deref(&self) -> &BifrostIndex {
        &self.index.index
    }
}
impl RunLease {
    pub fn snapshots(&self) -> Vec<ReservationSnapshot> {
        self.index
            .reservation
            .iter()
            .chain(self.reservation.iter())
            .map(|g| g.snapshot().clone())
            .collect()
    }
}

fn wait(
    context: &ExecutionResources,
    started: Instant,
    reason: Option<AdmissionWait>,
) -> Result<()> {
    if started.elapsed() >= context.timeout() {
        return Err(ResourceAdmissionError::new(anyhow::anyhow!(
            "quant resource admission timed out waiting for {}",
            reason.map_or("index loader".to_string(), |w| w.resource.to_string())
        ))
        .into());
    }
    std::thread::sleep(Duration::from_millis(50));
    Ok(())
}
fn announce(phase: &str, reason: &AdmissionWait, announced: &mut bool) {
    if !*announced {
        tracing::info!(
            phase,
            resource = reason.resource,
            requested = reason.requested,
            available = reason.available,
            live_reservations = reason.live_reservations,
            "quantification waiting for shared resources"
        );
        *announced = true;
    }
}
impl QuantIndexGroup {
    /// Safe only when the cache owns the last graph reference. Used before
    /// waiting for run capacity and when preprocessing must feed later waves.
    pub fn release_idle(&self) -> Result<bool> {
        let mut slot = match self.index.try_lock() {
            Ok(slot) => slot,
            Err(TryLockError::WouldBlock) => return Ok(false),
            Err(TryLockError::Poisoned(_)) => bail!("quant index loader panicked"),
        };
        if slot
            .as_ref()
            .is_some_and(|index| Arc::strong_count(index) == 1)
        {
            drop(slot.take());
            self.evictions.fetch_add(1, Ordering::Relaxed);
            return Ok(true);
        }
        Ok(false)
    }
    pub(in crate::quant) fn acquire_scan(
        &self,
        context: Option<&ExecutionResources>,
        request: ResourceRequest,
        check: &impl Fn() -> Result<()>,
    ) -> Result<Option<ResourceReservation>> {
        let Some(context) = context else {
            return Ok(None);
        };
        let started = Instant::now();
        let mut announced = false;
        loop {
            check()?;
            match context.try_reserve("scan", request.clone())? {
                Admission::Granted(guard) => {
                    check()?;
                    return Ok(Some(guard));
                }
                Admission::Wait(reason) => {
                    announce("scan", &reason, &mut announced);
                    if reason.resource == "memory" && context.pool().snapshot()?.cpu_threads == 0 {
                        self.release_idle()?;
                    }
                    wait(context, started, Some(reason))?;
                }
            }
        }
    }
    pub(in crate::quant) fn acquire(
        &self,
        reference: &ReferenceBundle,
        estimate: &QuantResources,
        run: ResourceRequest,
        context: Option<&ExecutionResources>,
        check: &impl Fn() -> Result<()>,
    ) -> Result<RunLease> {
        let Some(context) = context else {
            return Ok(RunLease {
                index: self.load(reference, estimate.threads)?,
                reservation: None,
            });
        };
        let started = Instant::now();
        let mut announced = false;
        loop {
            check()?;
            let mut slot = match self.index.try_lock() {
                Ok(slot) => slot,
                Err(TryLockError::WouldBlock) => {
                    wait(context, started, None)?;
                    continue;
                }
                Err(TryLockError::Poisoned(_)) => bail!("quant index loader panicked"),
            };
            if let Some(index) = slot.as_ref() {
                if index.pool_root.as_deref() != Some(context.pool().root()) {
                    bail!(
                        "quant group cannot mix admitted/unadmitted use or different resource pools"
                    );
                }
                let resident = index
                    .reservation
                    .as_ref()
                    .context("resident index lacks its resource owner")?
                    .snapshot()
                    .requested
                    .memory_bytes;
                if resident
                    .checked_add(run.memory_bytes)
                    .is_none_or(|n| n > context.limits().memory_bytes)
                {
                    return Err(ResourceAdmissionError::new(anyhow::anyhow!(
                        "resident index plus quant run exceeds its own memory capacity"
                    ))
                    .into());
                }
                match context.try_reserve("run", run.clone())? {
                    Admission::Granted(guard) => {
                        check()?;
                        self.verify(reference)?;
                        return Ok(RunLease {
                            index: Arc::clone(index),
                            reservation: Some(guard),
                        });
                    }
                    Admission::Wait(reason) => {
                        announce("run", &reason, &mut announced);
                        if reason.resource == "memory"
                            && Arc::strong_count(index) == 1
                            && context.pool().snapshot()?.cpu_threads == 0
                        {
                            drop(slot.take());
                            self.evictions.fetch_add(1, Ordering::Relaxed);
                            drop(slot);
                            // Retry cold with no resident allocation held. The
                            // first run and index must then fit atomically.
                            continue;
                        }
                        drop(slot);
                        wait(context, started, Some(reason))?;
                        continue;
                    }
                }
            }
            match context.try_reserve_bundle(vec![
                ("index", estimate.index_request(true)),
                ("run", run.clone()),
            ])? {
                BundleAdmission::Wait(reason) => {
                    announce("index_and_run", &reason, &mut announced);
                    drop(slot);
                    wait(context, started, Some(reason))?;
                }
                BundleAdmission::Granted(mut guards) => {
                    let run_guard = guards.pop().expect("run claim");
                    let mut index_guard = guards.pop().expect("index claim");
                    check()?;
                    let index = self.build(reference, estimate.threads)?;
                    check()?;
                    if index_guard
                        .try_resize(estimate.index_request(false))?
                        .is_some()
                    {
                        bail!("cannot reduce index loading reservation to resident ownership");
                    }
                    let index = Arc::new(ResidentIndex {
                        index,
                        reservation: Some(index_guard),
                        pool_root: Some(context.pool().root().into()),
                    });
                    *slot = Some(Arc::clone(&index));
                    self.loads.fetch_add(1, Ordering::Relaxed);
                    return Ok(RunLease {
                        index,
                        reservation: Some(run_guard),
                    });
                }
            }
        }
    }
}
