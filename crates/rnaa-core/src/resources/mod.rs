//! Cross-process admission accounting, held by kernel locks rather than timeouts.
mod admission;
mod bundle;
mod files;
mod platform;
mod policy;
mod reservation;
mod store;
mod types;
pub use policy::{ExecutionResources, ResourceAdmissionError, ResourceConfig};
pub use reservation::ResourceReservation;
pub use store::ResourcePool;
pub use types::{
    Admission, AdmissionWait, BundleAdmission, DiskRequest, OwnerFileIdentity, PoolLimits,
    PoolSnapshot, ReservationClaim, ReservationOwner, ReservationSnapshot, ResourceRecovery,
    ResourceRequest,
};

/// Available bytes on the filesystem containing the nearest existing ancestor.
pub fn available_disk_bytes(path: &std::path::Path) -> anyhow::Result<u64> {
    platform::available_disk_bytes(&files::existing_path(path)?)
}

#[cfg(test)]
mod tests;
