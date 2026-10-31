use super::{ReservationSnapshot, ResourcePool};
use anyhow::{Result, bail};
use std::fs::File;

pub struct ResourceReservation {
    pub(super) pool: ResourcePool,
    pub(super) snapshot: Box<ReservationSnapshot>,
    pub(super) owner_file: Option<File>,
}
impl ResourceReservation {
    pub fn snapshot(&self) -> &ReservationSnapshot {
        &self.snapshot
    }
    pub fn release(mut self) -> Result<()> {
        self.release_inner()
    }

    fn release_inner(&mut self) -> Result<()> {
        if self.owner_file.is_none() {
            return Ok(());
        }
        let catalog = self.pool.catalog()?;
        let deleted = catalog.connection.execute(
            "DELETE FROM reservations WHERE reservation_id=?1 AND descriptor_json=?2",
            rusqlite::params![
                self.snapshot.reservation_id,
                serde_json::to_string(&self.snapshot)?
            ],
        )?;
        if deleted != 1 {
            bail!("resource reservation changed outside its owner; release is held");
        }
        drop(self.owner_file.take());
        let path = self.pool.owner_path(&self.snapshot.reservation_id)?;
        super::files::remove_owner(&path, self.snapshot.owner_file_identity)?;
        File::open(path.parent().unwrap())?.sync_all()?;
        Ok(())
    }
}
impl Drop for ResourceReservation {
    fn drop(&mut self) {
        if let Err(error) = self.release_inner() {
            tracing::warn!(reservation=%self.snapshot.reservation_id,%error,"resource reservation release will require lock-based recovery");
        }
    }
}
