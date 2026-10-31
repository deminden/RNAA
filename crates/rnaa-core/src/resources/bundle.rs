//! Atomic admission for allocations whose ownership lifetimes will diverge.
use super::{
    BundleAdmission, PoolLimits, PoolSnapshot, ReservationClaim, ReservationSnapshot, ResourcePool,
    ResourceRequest, ResourceReservation, admission::denied, platform,
};
use anyhow::{Context, Result, bail};
use std::{
    fs::{File, OpenOptions},
    path::Path,
    time::{Duration, Instant},
};

// Deliberately not a ResourceReservation: rollback must close its file without
// attempting another catalog transaction while this admission owns the lock.
struct Pending {
    snapshot: Box<ReservationSnapshot>,
    file: File,
}
impl ResourcePool {
    pub fn try_reserve_bundle(
        &self,
        limits: PoolLimits,
        claims: Vec<ReservationClaim>,
    ) -> Result<BundleAdmission> {
        self.reserve_bundle_with_probe(limits, claims, &platform::available_disk_bytes)
    }
    pub(super) fn reserve_bundle_with_probe(
        &self,
        limits: PoolLimits,
        claims: Vec<ReservationClaim>,
        probe: &impl Fn(&Path) -> Result<u64>,
    ) -> Result<BundleAdmission> {
        if claims.is_empty() {
            bail!("resource admission bundle must not be empty");
        }
        let mut total = ResourceRequest::default();
        let mut normalized = Vec::with_capacity(claims.len());
        for claim in claims {
            claim.owner.validate()?;
            claim.request.validate(limits)?;
            let request = self.normalize_request(claim.request)?;
            total.cpu_threads = total
                .cpu_threads
                .checked_add(request.cpu_threads)
                .context("bundle CPU total overflows")?;
            total.memory_bytes = total
                .memory_bytes
                .checked_add(request.memory_bytes)
                .context("bundle memory total overflows")?;
            total.disks.extend(request.disks.iter().cloned());
            normalized.push(ReservationClaim {
                owner: claim.owner,
                request,
            });
        }
        // A bundle that cannot fit by itself must not wait for capacity that
        // no other owner can release. Admission checks the sum, never a prefix.
        total.validate(limits)?;
        let mut catalog = self.catalog()?;
        let transaction = catalog.connection.transaction()?;
        let (live, reclaimed) = self.live_rows(&transaction)?;
        let usage = PoolSnapshot::from_rows(live)?;
        if let Some(wait) = denied(&usage, limits, &total, probe)? {
            transaction.commit()?;
            self.clean_reclaimed(reclaimed);
            return Ok(BundleAdmission::Wait(wait));
        }
        let mut pending = Vec::with_capacity(normalized.len());
        for claim in normalized {
            let id = uuid::Uuid::new_v4().to_string();
            let path = self.owner_path(&id)?;
            let mut options = OpenOptions::new();
            options.read(true).write(true).create_new(true);
            #[cfg(unix)]
            {
                use std::os::unix::fs::OpenOptionsExt;
                options.mode(0o600);
            }
            let file = options.open(&path)?;
            file.try_lock()
                .context("cannot hold a bundled resource owner lock")?;
            file.sync_all()?;
            let snapshot = Box::new(ReservationSnapshot {
                reservation_id: id,
                owner: claim.owner,
                limits,
                requested: claim.request,
                created_at: crate::util::now_rfc3339(),
                owner_file_identity: super::files::identity(&file.metadata()?)?,
            });
            transaction.execute(
                "INSERT INTO reservations(reservation_id,descriptor_json) VALUES(?1,?2)",
                rusqlite::params![snapshot.reservation_id, serde_json::to_string(&snapshot)?],
            )?;
            pending.push(Pending { snapshot, file });
        }
        File::open(self.root.join("owners"))?.sync_all()?;
        transaction.commit()?;
        self.clean_reclaimed(reclaimed);
        // Only committed rows become independently releasing guards. A failed
        // transaction leaves unlocked empty files for explicit orphan recovery.
        Ok(BundleAdmission::Granted(
            pending
                .into_iter()
                .map(|p| ResourceReservation {
                    pool: self.clone(),
                    snapshot: p.snapshot,
                    owner_file: Some(p.file),
                })
                .collect(),
        ))
    }
    pub fn reserve_bundle(
        &self,
        limits: PoolLimits,
        claims: Vec<ReservationClaim>,
        timeout: Duration,
        mut check: impl FnMut() -> Result<()>,
    ) -> Result<Vec<ResourceReservation>> {
        let started = Instant::now();
        loop {
            check()?;
            match self.try_reserve_bundle(limits, claims.clone())? {
                BundleAdmission::Granted(guards) => {
                    check()?;
                    return Ok(guards);
                }
                BundleAdmission::Wait(wait) if started.elapsed() >= timeout => bail!(
                    "resource bundle admission timed out waiting for {}: requested {}, available {}, live reservations {}",
                    wait.resource,
                    wait.requested,
                    wait.available,
                    wait.live_reservations
                ),
                BundleAdmission::Wait(_) => std::thread::sleep(Duration::from_millis(50)),
            }
        }
    }
}

#[cfg(test)]
mod tests;
