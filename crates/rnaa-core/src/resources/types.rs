use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Cooperative capacities. A new admission must honor both its own limits and
/// every live reservation's limits; these are not operating-system RSS caps.
pub struct PoolLimits {
    pub cpu_threads: u64,
    pub memory_bytes: u64,
    pub disk_headroom_bytes: u64,
}
impl PoolLimits {
    pub(super) fn validate(self) -> Result<()> {
        if self.cpu_threads == 0 || self.memory_bytes == 0 {
            bail!("resource pool CPU and memory capacities must be positive");
        }
        Ok(())
    }
    pub(super) fn intersect(self, other: Self) -> Self {
        Self {
            cpu_threads: self.cpu_threads.min(other.cpu_threads),
            memory_bytes: self.memory_bytes.min(other.memory_bytes),
            disk_headroom_bytes: self.disk_headroom_bytes.max(other.disk_headroom_bytes),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiskRequest {
    pub path: PathBuf,
    pub bytes: u64,
}
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// CPU and memory held for the lifetime of work, plus additional bytes still to
/// be written. Callers must bound their work and reserve before allocation/growth.
pub struct ResourceRequest {
    pub cpu_threads: u64,
    pub memory_bytes: u64,
    pub disks: Vec<DiskRequest>,
}
impl ResourceRequest {
    pub fn disk_bytes(&self) -> Result<u64> {
        self.disks.iter().try_fold(0_u64, |sum, disk| {
            sum.checked_add(disk.bytes)
                .ok_or_else(|| anyhow::anyhow!("disk reservation total overflows"))
        })
    }
    pub(super) fn validate(&self, limits: PoolLimits) -> Result<()> {
        limits.validate()?;
        if self.cpu_threads > limits.cpu_threads || self.memory_bytes > limits.memory_bytes {
            bail!(
                "resource request exceeds its own CPU or memory capacity; reduce the request or increase its declared limit"
            );
        }
        if self.cpu_threads == 0 && self.memory_bytes == 0 && self.disk_bytes()? == 0 {
            bail!("resource request must reserve CPU, memory or disk");
        }
        if self
            .disks
            .iter()
            .any(|disk| disk.path.as_os_str().is_empty() || disk.bytes == 0)
        {
            bail!("disk reservation paths and byte allowances must be nonempty");
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReservationOwner {
    pub project_root: PathBuf,
    pub attempt_id: String,
    pub stage: String,
    pub run_accession: Option<String>,
}

/// All claims are admitted together, but each returned owner can release or
/// resize its allocation independently after the shared preparation finishes.
#[derive(Debug, Clone)]
pub struct ReservationClaim {
    pub owner: ReservationOwner,
    pub request: ResourceRequest,
}
pub enum BundleAdmission {
    Granted(Vec<super::ResourceReservation>),
    Wait(AdmissionWait),
}
impl ReservationOwner {
    pub(super) fn validate(&self) -> Result<()> {
        if !self.project_root.is_absolute() || self.attempt_id.is_empty() || self.stage.is_empty() {
            bail!(
                "resource reservation requires an absolute project root, attempt identity and stage"
            );
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReservationSnapshot {
    pub reservation_id: String,
    pub owner: ReservationOwner,
    pub limits: PoolLimits,
    pub requested: ResourceRequest,
    pub created_at: String,
    pub owner_file_identity: OwnerFileIdentity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OwnerFileIdentity {
    pub device: u64,
    pub inode: u64,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize)]
pub struct ResourceRecovery {
    pub abandoned_reservations: usize,
    pub orphan_lock_files: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct PoolSnapshot {
    pub reservations: Vec<ReservationSnapshot>,
    pub cpu_threads: u64,
    pub memory_bytes: u64,
    pub pending_disk_bytes: u64,
    pub effective_limits: Option<PoolLimits>,
}
impl PoolSnapshot {
    pub(super) fn from_rows(reservations: Vec<ReservationSnapshot>) -> Result<Self> {
        let mut value = Self {
            reservations: vec![],
            cpu_threads: 0,
            memory_bytes: 0,
            pending_disk_bytes: 0,
            effective_limits: None,
        };
        for row in &reservations {
            row.owner.validate()?;
            row.requested.validate(row.limits)?;
            value.cpu_threads = value
                .cpu_threads
                .checked_add(row.requested.cpu_threads)
                .ok_or_else(|| anyhow::anyhow!("CPU reservation total overflows"))?;
            value.memory_bytes = value
                .memory_bytes
                .checked_add(row.requested.memory_bytes)
                .ok_or_else(|| anyhow::anyhow!("memory reservation total overflows"))?;
            value.pending_disk_bytes = value
                .pending_disk_bytes
                .checked_add(row.requested.disk_bytes()?)
                .ok_or_else(|| anyhow::anyhow!("disk reservation total overflows"))?;
            value.effective_limits = Some(
                value
                    .effective_limits
                    .map_or(row.limits, |limits| limits.intersect(row.limits)),
            );
        }
        if value.effective_limits.is_some_and(|limits| {
            value.cpu_threads > limits.cpu_threads || value.memory_bytes > limits.memory_bytes
        }) {
            bail!("resource catalog contains reservations beyond its live capacities");
        }
        value.reservations = reservations;
        Ok(value)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct AdmissionWait {
    pub resource: &'static str,
    pub requested: u64,
    pub available: u64,
    pub live_reservations: usize,
}
pub enum Admission {
    Granted(super::ResourceReservation),
    Wait(AdmissionWait),
}
