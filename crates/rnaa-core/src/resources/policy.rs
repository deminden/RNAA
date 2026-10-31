//! Common host policy and owned phase admission; independent of scientific keys.
use super::{PoolLimits, ReservationOwner, ResourcePool, ResourceRequest, ResourceReservation};
use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::{
    path::{Path, PathBuf},
    time::Duration,
};

const MIB: u64 = 1024 * 1024;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ResourceConfig {
    /// Empty selects a common per-user host pool. Overrides must be absolute.
    pub pool_root: String,
    /// Zero uses the host's available parallelism.
    pub cpu_threads: u64,
    /// Zero uses three quarters of detected physical/container capacity.
    pub memory_mb: u64,
    pub disk_headroom_mb: u64,
    pub wait_seconds: u64,
}
impl Default for ResourceConfig {
    fn default() -> Self {
        Self {
            pool_root: String::new(),
            cpu_threads: 0,
            memory_mb: 0,
            disk_headroom_mb: 1024,
            wait_seconds: 300,
        }
    }
}
impl ResourceConfig {
    pub fn validate(&self) -> Result<()> {
        if !self.pool_root.is_empty() && !Path::new(&self.pool_root).is_absolute() {
            bail!("resources.pool_root must be absolute or empty for the shared host default");
        }
        self.memory_mb
            .checked_mul(MIB)
            .context("resources.memory_mb overflows byte accounting")?;
        self.disk_headroom_mb
            .checked_mul(MIB)
            .context("resources.disk_headroom_mb overflows byte accounting")?;
        if self.wait_seconds == 0 {
            bail!("resources.wait_seconds must be positive");
        }
        Ok(())
    }
    fn pool_path(&self) -> Result<PathBuf> {
        if !self.pool_root.is_empty() {
            return Ok(PathBuf::from(&self.pool_root));
        }
        let home = std::env::var_os("HOME")
            .map(PathBuf::from)
            .filter(|p| p.is_absolute())
            .context(
                "cannot select a common resource pool: configure absolute resources.pool_root",
            )?;
        Ok(home.join(".cache/rnaa/resources/v1"))
    }
    pub fn limits(&self) -> Result<PoolLimits> {
        self.validate()?;
        let cpu_threads = if self.cpu_threads == 0 {
            std::thread::available_parallelism()?.get() as u64
        } else {
            self.cpu_threads
        };
        let memory_bytes = if self.memory_mb == 0 {
            let capacity = super::platform::memory_capacity_bytes()?;
            capacity - capacity / 4
        } else {
            self.memory_mb * MIB
        };
        let limits = PoolLimits {
            cpu_threads,
            memory_bytes,
            disk_headroom_bytes: self.disk_headroom_mb * MIB,
        };
        limits.validate()?;
        Ok(limits)
    }
}

/// Admission failures must not become an engine bypass that escapes accounting.
#[derive(Debug, thiserror::Error)]
#[error("resource admission failed: {0:#}")]
pub struct ResourceAdmissionError(#[source] anyhow::Error);
impl ResourceAdmissionError {
    pub fn new(error: anyhow::Error) -> Self {
        Self(error)
    }
}

#[derive(Clone, Debug)]
pub struct ExecutionResources {
    pool: ResourcePool,
    limits: PoolLimits,
    owner: ReservationOwner,
    timeout: Duration,
}
impl ExecutionResources {
    pub fn open(config: &ResourceConfig, owner: ReservationOwner) -> Result<Self> {
        (|| {
            config.validate()?;
            owner.validate()?;
            let limits = config.limits()?;
            let pool = ResourcePool::open(&config.pool_path()?)?;
            Ok(Self {
                pool,
                limits,
                owner,
                timeout: Duration::from_secs(config.wait_seconds),
            })
        })()
        .map_err(|e| ResourceAdmissionError(e).into())
    }
    pub fn limits(&self) -> PoolLimits {
        self.limits
    }
    pub fn pool(&self) -> &ResourcePool {
        &self.pool
    }
    pub fn timeout(&self) -> Duration {
        self.timeout
    }
    fn phase_owner(&self, phase: &str) -> ReservationOwner {
        let mut owner = self.owner.clone();
        owner.stage = format!("{}.{}", owner.stage, phase);
        owner
    }
    pub fn try_reserve(&self, phase: &str, request: ResourceRequest) -> Result<super::Admission> {
        self.pool
            .try_reserve(self.limits, self.phase_owner(phase), request)
            .map_err(|e| ResourceAdmissionError(e).into())
    }
    pub fn try_reserve_bundle(
        &self,
        phases: Vec<(&str, ResourceRequest)>,
    ) -> Result<super::BundleAdmission> {
        let claims = phases
            .into_iter()
            .map(|(phase, request)| super::ReservationClaim {
                owner: self.phase_owner(phase),
                request,
            })
            .collect();
        self.pool
            .try_reserve_bundle(self.limits, claims)
            .map_err(|e| ResourceAdmissionError(e).into())
    }
    pub fn reserve(
        &self,
        phase: &str,
        request: ResourceRequest,
        check: &impl Fn() -> Result<()>,
    ) -> Result<ResourceReservation> {
        self.pool
            .reserve(
                self.limits,
                self.phase_owner(phase),
                request,
                self.timeout,
                check,
            )
            .map_err(|e| ResourceAdmissionError(e).into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn legacy_policy_defaults_and_overflow_are_explicit() -> Result<()> {
        let config: crate::ProjectConfig = serde_json::from_value(serde_json::json!({
            "project":{"name":"legacy","schema_version":1}
        }))?;
        assert_eq!(config.resources.wait_seconds, 300);
        assert!(config.resources.limits()?.memory_bytes > 0);
        let mut policy = ResourceConfig {
            pool_root: "relative".into(),
            ..Default::default()
        };
        assert!(policy.validate().is_err());
        policy.pool_root.clear();
        policy.memory_mb = u64::MAX;
        assert!(policy.validate().is_err());
        policy.memory_mb = 1;
        policy.disk_headroom_mb = u64::MAX;
        assert!(policy.validate().is_err());
        policy.disk_headroom_mb = 0;
        policy.wait_seconds = 0;
        assert!(policy.validate().is_err());
        Ok(())
    }
    #[test]
    fn owned_phase_admission_preserves_fences_and_typed_errors() -> Result<()> {
        let root = tempfile::tempdir()?;
        let resources = ExecutionResources::open(
            &ResourceConfig {
                pool_root: root.path().join("pool").display().to_string(),
                cpu_threads: 1,
                memory_mb: 1,
                disk_headroom_mb: 0,
                wait_seconds: 1,
            },
            ReservationOwner {
                project_root: root.path().into(),
                attempt_id: "attempt".into(),
                stage: "preprocess".into(),
                run_accession: Some("RUN1".into()),
            },
        )?;
        let request = ResourceRequest {
            cpu_threads: 1,
            memory_bytes: 512,
            ..Default::default()
        };
        let guard = resources.reserve("scan", request.clone(), &|| Ok(()))?;
        let snapshot = resources.pool().snapshot()?;
        assert_eq!(snapshot.reservations[0].owner.stage, "preprocess.scan");
        let error = resources.reserve("engine", request.clone(), &|| anyhow::bail!("cancelled"));
        assert!(
            error
                .err()
                .unwrap()
                .downcast_ref::<ResourceAdmissionError>()
                .is_some()
        );
        drop(guard);
        assert!(resources.pool().snapshot()?.reservations.is_empty());
        let error = resources.reserve(
            "engine",
            ResourceRequest {
                memory_bytes: 2 * MIB,
                ..request
            },
            &|| Ok(()),
        );
        assert!(
            error
                .err()
                .unwrap()
                .downcast_ref::<ResourceAdmissionError>()
                .is_some()
        );
        Ok(())
    }
}
