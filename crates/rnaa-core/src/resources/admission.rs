use super::{
    Admission, AdmissionWait, PoolLimits, PoolSnapshot, ReservationOwner, ResourcePool,
    ResourceRequest, ResourceReservation, platform,
};
use anyhow::{Context, Result, bail};
use std::{
    path::Path,
    time::{Duration, Instant},
};

pub(super) fn denied(
    usage: &PoolSnapshot,
    limits: PoolLimits,
    request: &ResourceRequest,
    probe: &impl Fn(&Path) -> Result<u64>,
) -> Result<Option<AdmissionWait>> {
    let limits = usage
        .effective_limits
        .map_or(limits, |active| limits.intersect(active));
    for (resource, requested, used, capacity) in [
        (
            "CPU",
            request.cpu_threads,
            usage.cpu_threads,
            limits.cpu_threads,
        ),
        (
            "memory",
            request.memory_bytes,
            usage.memory_bytes,
            limits.memory_bytes,
        ),
    ] {
        if used > capacity {
            return Ok(Some(AdmissionWait {
                resource,
                requested: used,
                available: capacity,
                live_reservations: usage.reservations.len(),
            }));
        }
        let available = capacity - used;
        if requested > available {
            return Ok(Some(AdmissionWait {
                resource,
                requested,
                available,
                live_reservations: usage.reservations.len(),
            }));
        }
    }
    // Paths can share physical backing even when their filesystem identifiers
    // differ (APFS volumes, subvolumes and overlays). Sum every pending promise
    // instead of treating such paths as independent copies of free space.
    let requested = request.disk_bytes()?;
    for disk in &request.disks {
        let free = probe(&disk.path)?;
        let available = free
            .saturating_sub(limits.disk_headroom_bytes)
            .saturating_sub(usage.pending_disk_bytes);
        if requested > available {
            return Ok(Some(AdmissionWait {
                resource: "disk",
                requested,
                available,
                live_reservations: usage.reservations.len(),
            }));
        }
    }
    Ok(None)
}

impl ResourcePool {
    pub fn try_reserve(
        &self,
        limits: PoolLimits,
        owner: ReservationOwner,
        request: ResourceRequest,
    ) -> Result<Admission> {
        self.reserve_with_probe(limits, owner, request, &platform::available_disk_bytes)
    }

    pub(super) fn reserve_with_probe(
        &self,
        limits: PoolLimits,
        owner: ReservationOwner,
        request: ResourceRequest,
        probe: &impl Fn(&Path) -> Result<u64>,
    ) -> Result<Admission> {
        match self.reserve_bundle_with_probe(
            limits,
            vec![super::ReservationClaim { owner, request }],
            probe,
        )? {
            super::BundleAdmission::Granted(mut guards) => Ok(Admission::Granted(
                guards.pop().expect("one admitted claim"),
            )),
            super::BundleAdmission::Wait(wait) => Ok(Admission::Wait(wait)),
        }
    }

    /// Wait without holding a catalog transaction. Callers supply their own
    /// cancellation/ownership fence, which is checked before every admission.
    pub fn reserve(
        &self,
        limits: PoolLimits,
        owner: ReservationOwner,
        request: ResourceRequest,
        timeout: Duration,
        mut check: impl FnMut() -> Result<()>,
    ) -> Result<ResourceReservation> {
        let started = Instant::now();
        let mut announced = false;
        loop {
            check()?;
            match self.try_reserve(limits, owner.clone(), request.clone())? {
                Admission::Granted(reservation) => {
                    check()?;
                    return Ok(reservation);
                }
                Admission::Wait(wait) if started.elapsed() >= timeout => bail!(
                    "resource admission timed out waiting for {}: requested {}, available {}, live reservations {}",
                    wait.resource,
                    wait.requested,
                    wait.available,
                    wait.live_reservations
                ),
                Admission::Wait(wait) => {
                    if !announced {
                        tracing::info!(stage=%owner.stage, attempt=%owner.attempt_id,
                            resource=wait.resource, requested=wait.requested, available=wait.available,
                            live_reservations=wait.live_reservations,
                            "waiting for shared resource capacity");
                        announced = true;
                    }
                    std::thread::sleep(Duration::from_millis(50));
                }
            }
        }
    }
}

impl ResourceReservation {
    pub fn try_resize(&mut self, request: ResourceRequest) -> Result<Option<AdmissionWait>> {
        request.validate(self.snapshot.limits)?;
        let request = self.pool.normalize_request(request)?;
        if self.owner_file.is_none() {
            bail!("resource reservation has been released");
        }
        let mut catalog = self.pool.catalog()?;
        let transaction = catalog.connection.transaction()?;
        let (mut live, reclaimed) = self.pool.live_rows(&transaction)?;
        let position = live
            .iter()
            .position(|row| row.reservation_id == self.snapshot.reservation_id)
            .context("resource reservation lost its catalog record")?;
        if live.remove(position) != *self.snapshot {
            bail!("resource reservation changed outside its owner");
        }
        let usage = PoolSnapshot::from_rows(live)?;
        if let Some(wait) = denied(
            &usage,
            self.snapshot.limits,
            &request,
            &platform::available_disk_bytes,
        )? {
            transaction.commit()?;
            self.pool.clean_reclaimed(reclaimed);
            return Ok(Some(wait));
        }
        let mut snapshot = self.snapshot.clone();
        snapshot.requested = request;
        transaction.execute(
            "UPDATE reservations SET descriptor_json=?1 WHERE reservation_id=?2",
            rusqlite::params![serde_json::to_string(&snapshot)?, snapshot.reservation_id],
        )?;
        transaction.commit()?;
        self.pool.clean_reclaimed(reclaimed);
        self.snapshot = snapshot;
        Ok(None)
    }
}
