//! Owned phase estimates; resource settings never enter source request identity.
use anyhow::{Context, Result, bail};
use rnaa_core::{
    VerifiedFile,
    resources::{DiskRequest, ExecutionResources, ResourceRequest, ResourceReservation},
};
use serde::{Serialize, de::DeserializeOwned};
use std::{fs::File, io::Read, path::Path};
const MIB: u64 = 1024 * 1024;

pub(super) fn receipt_cap(remote: &rnaa_core::RemoteFile, target: &Path) -> Result<u64> {
    let descriptor = serde_json::to_vec(remote)?.len() as u128 + target.as_os_str().len() as u128;
    u64::try_from(4 * u128::from(MIB) + descriptor * 16)
        .context("download receipt allowance overflows")
}
pub(super) fn verify_request(
    remote: &rnaa_core::RemoteFile,
    target: &Path,
) -> Result<ResourceRequest> {
    let metadata = receipt_cap(remote, target)?;
    Ok(ResourceRequest {
        cpu_threads: 1,
        memory_bytes: metadata
            .checked_mul(4)
            .and_then(|n| n.checked_add(16 * MIB))
            .context("download verification memory overflows")?,
        disks: vec![DiskRequest {
            path: target.into(),
            bytes: metadata
                .checked_mul(2)
                .context("download receipt disk allowance overflows")?,
        }],
    })
}
pub(super) fn reserve(
    context: Option<&ExecutionResources>,
    phase: &str,
    request: ResourceRequest,
    check: &impl Fn() -> Result<()>,
) -> Result<Option<ResourceReservation>> {
    context
        .map(|c| c.reserve(phase, request, check))
        .transpose()
}

#[derive(Debug, Serialize)]
pub(super) struct PublicationPlan {
    pub model: &'static str,
    pub copy_bytes: u64,
    pub metadata_bytes: u64,
    pub memory_bytes: u64,
}
impl PublicationPlan {
    pub fn inspect(files: &[VerifiedFile], run: &rnaa_core::RunRecord) -> Result<Self> {
        let copy = files.iter().try_fold(0_u64, |n, f| {
            n.checked_add(f.bytes)
                .context("download publication bytes overflow")
        })?;
        let descriptors = serde_json::to_vec(run)?.len() as u128
            + files
                .iter()
                .map(|f| f.path.as_os_str().len() as u128 + f.checksum.len() as u128 + 256)
                .sum::<u128>();
        let metadata =
            u64::try_from(64 * 1024 + files.len() as u128 * 4 * u128::from(MIB) + descriptors * 16)
                .context("download publication metadata overflows")?;
        let memory = metadata
            .checked_mul(4)
            .and_then(|n| n.checked_add(64 * MIB))
            .context("download publication memory overflows")?;
        Ok(Self {
            model: "download-publication-v1",
            copy_bytes: copy,
            metadata_bytes: metadata,
            memory_bytes: memory,
        })
    }
    pub fn request(&self, generations: &Path, catalog: &Path) -> Result<ResourceRequest> {
        Ok(ResourceRequest {
            cpu_threads: 1,
            memory_bytes: self.memory_bytes,
            disks: vec![
                DiskRequest {
                    path: generations.into(),
                    bytes: self
                        .copy_bytes
                        .checked_add(self.metadata_bytes)
                        .context("download publication disk allowance overflows")?,
                },
                DiskRequest {
                    path: catalog.into(),
                    bytes: self.metadata_bytes,
                },
            ],
        })
    }
}
pub(super) fn read_json<T: DeserializeOwned>(path: &Path, cap: u64) -> Result<T> {
    let mut raw = Vec::new();
    File::open(path)?
        .take(cap.checked_add(1).context("receipt cap overflows")?)
        .read_to_end(&mut raw)?;
    if raw.len() as u64 > cap {
        bail!("download receipt exceeds its metadata allowance");
    }
    Ok(serde_json::from_slice(&raw)?)
}
pub(super) fn write_json(path: &Path, value: &impl Serialize, cap: u64) -> Result<()> {
    if serde_json::to_vec_pretty(value)?.len() as u64 > cap {
        bail!("download metadata exceeds its write allowance");
    }
    rnaa_core::util::write_json_pretty(path, value)
}

pub(super) use rnaa_core::util::copy_bounded;
