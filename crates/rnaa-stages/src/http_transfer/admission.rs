//! Future-write windows and lifetime claims for a bounded HTTP request.
use anyhow::{Context, Result};
use rnaa_core::resources::{
    DiskRequest, ExecutionResources, ReservationSnapshot, ResourceAdmissionError, ResourceRequest,
    ResourceReservation,
};
use std::path::{Path, PathBuf};

const MIB: u64 = 1024 * 1024;
pub(super) const METADATA_BYTES: u64 = 4 * MIB;
const WINDOW_BYTES: u64 = 8 * MIB;

#[derive(Debug)]
pub(super) struct ReadmitLength(pub u64);
impl std::fmt::Display for ReadmitLength {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "HTTP size discovery requires readmission for {} total bytes",
            self.0
        )
    }
}
impl std::error::Error for ReadmitLength {}

#[derive(Debug, Clone, serde::Serialize)]
pub(crate) struct TransferExecution {
    model: &'static str,
    max_pending_disk_bytes: u64,
    reservation: ReservationSnapshot,
}

pub(super) fn request(path: &Path, future: u64, scanning: bool) -> Result<ResourceRequest> {
    Ok(ResourceRequest {
        cpu_threads: 1,
        memory_bytes: if scanning { 16 * MIB } else { 64 * MIB },
        disks: vec![DiskRequest {
            path: path.into(),
            bytes: future
                .checked_add(METADATA_BYTES)
                .context("HTTP disk allowance overflows")?,
        }],
    })
}

pub(super) struct TransferBudget<'a> {
    context: Option<&'a ExecutionResources>,
    guard: Option<ResourceReservation>,
    path: PathBuf,
    end: u64,
    total: Option<u64>,
    maximum_pending: u64,
}
impl<'a> TransferBudget<'a> {
    pub fn reserve(
        context: Option<&'a ExecutionResources>,
        path: &Path,
        start: u64,
        total: Option<u64>,
        check: &impl Fn() -> Result<()>,
    ) -> Result<Self> {
        let future = Self::future(start, total)?;
        let claim = request(path, future, false)?;
        let guard = context
            .map(|c| c.reserve("transfer", claim, check))
            .transpose()?;
        Ok(Self {
            context,
            guard,
            path: path.into(),
            end: start
                .checked_add(future)
                .context("HTTP window end overflows")?,
            total,
            maximum_pending: future + METADATA_BYTES,
        })
    }
    fn future(written: u64, total: Option<u64>) -> Result<u64> {
        total.map_or(Ok(WINDOW_BYTES), |n| {
            n.checked_sub(written)
                .context("HTTP prefix exceeds total size")
        })
    }
    fn resize(&mut self, written: u64, total: Option<u64>) -> Result<bool> {
        let future = Self::future(written, total)?;
        let claim = request(&self.path, future, false)?;
        if let Some(guard) = &mut self.guard
            && guard
                .try_resize(claim)
                .map_err(ResourceAdmissionError::new)?
                .is_some()
        {
            return Ok(false);
        }
        self.maximum_pending = self.maximum_pending.max(future + METADATA_BYTES);
        self.end = written
            .checked_add(future)
            .context("HTTP window end overflows")?;
        self.total = total;
        Ok(true)
    }
    pub fn begin(&mut self, start: u64, total: Option<u64>) -> Result<()> {
        if !self.resize(start, total)? {
            // Unwind the response and guard before waiting for the discovered
            // complete size. A second request then runs under the full claim.
            if let Some(total) = total {
                return Err(ReadmitLength(total).into());
            }
            return Err(ResourceAdmissionError::new(anyhow::anyhow!(
                "HTTP initial write window lacks shared disk capacity"
            ))
            .into());
        }
        Ok(())
    }
    pub fn before_write(&mut self, written: u64, next: u64) -> Result<()> {
        if self.context.is_some() && next > self.end && !self.resize(written, self.total)? {
            // Unknown-length growth cannot wait while retaining a response and
            // CPU/memory ownership. The caller saves its proven prefix and exits.
            return Err(ResourceAdmissionError::new(anyhow::anyhow!("HTTP write window cannot grow within shared disk capacity; saved prefix can be retried")).into());
        }
        Ok(())
    }
    pub fn progress(&mut self, written: u64) -> Result<()> {
        if self.total.is_some() && !self.resize(written, self.total)? {
            return Err(ResourceAdmissionError::new(anyhow::anyhow!(
                "HTTP remaining-write claim no longer fits the live policy"
            ))
            .into());
        }
        Ok(())
    }
    pub fn finish(&mut self, written: u64) -> Result<()> {
        self.total = Some(written);
        self.progress(written)
    }
    pub fn execution(&self) -> Option<TransferExecution> {
        self.guard.as_ref().map(|g| TransferExecution {
            model: "http-streaming-window-v1",
            max_pending_disk_bytes: self.maximum_pending,
            reservation: g.snapshot().clone(),
        })
    }
}
