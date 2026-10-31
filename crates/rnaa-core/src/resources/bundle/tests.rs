use super::super::{Admission, ReservationOwner};
use super::*;
use std::{
    fs,
    sync::{Arc, Barrier},
};
fn limits() -> PoolLimits {
    PoolLimits {
        cpu_threads: 2,
        memory_bytes: 100,
        disk_headroom_bytes: 0,
    }
}
fn claim(root: &Path, label: &str, cpu: u64, memory: u64) -> ReservationClaim {
    ReservationClaim {
        owner: ReservationOwner {
            project_root: root.into(),
            attempt_id: label.into(),
            stage: label.into(),
            run_accession: None,
        },
        request: ResourceRequest {
            cpu_threads: cpu,
            memory_bytes: memory,
            disks: Vec::new(),
        },
    }
}
fn granted(value: BundleAdmission) -> Vec<ResourceReservation> {
    match value {
        BundleAdmission::Granted(guards) => guards,
        BundleAdmission::Wait(_) => panic!("unexpected bundle wait"),
    }
}

#[test]
fn bundle_admission_is_atomic_and_owners_release_independently() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let mut held = granted(pool.try_reserve_bundle(
        limits(),
        vec![
            claim(root.path(), "index", 0, 40),
            claim(root.path(), "run", 1, 50),
        ],
    )?);
    assert_eq!(pool.snapshot()?.memory_bytes, 90);
    assert_eq!(pool.snapshot()?.reservations.len(), 2);
    let rejected = pool.try_reserve_bundle(
        limits(),
        vec![
            claim(root.path(), "other-index", 0, 5),
            claim(root.path(), "other-run", 1, 10),
        ],
    )?;
    assert!(matches!(rejected, BundleAdmission::Wait(_)));
    assert_eq!(pool.snapshot()?.reservations.len(), 2);
    drop(held.pop());
    assert_eq!(pool.snapshot()?.memory_bytes, 40);
    assert_eq!(pool.snapshot()?.cpu_threads, 0);
    assert_eq!(pool.snapshot()?.reservations[0].owner.stage, "index");
    assert!(
        held[0]
            .try_resize(ResourceRequest {
                memory_bytes: 20,
                ..Default::default()
            })?
            .is_none()
    );
    assert_eq!(pool.snapshot()?.memory_bytes, 20);
    drop(held);
    assert!(pool.snapshot()?.reservations.is_empty());
    assert!(
        pool.try_reserve_bundle(
            limits(),
            vec![
                claim(root.path(), "a", 0, 60),
                claim(root.path(), "b", 1, 60)
            ]
        )
        .is_err()
    );
    assert!(pool.try_reserve_bundle(limits(), vec![]).is_err());
    Ok(())
}

#[test]
fn failed_bundle_insert_rolls_back_every_owner_and_recovers_orphans() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    pool.catalog()?.connection.execute_batch(
        "CREATE TRIGGER reject_bundle BEFORE INSERT ON reservations
        WHEN json_extract(NEW.descriptor_json,'$.owner.stage')='reject'
        BEGIN SELECT RAISE(ABORT,'injected second owner failure'); END;",
    )?;
    assert!(
        pool.try_reserve_bundle(
            limits(),
            vec![
                claim(root.path(), "index", 0, 40),
                claim(root.path(), "reject", 1, 50)
            ]
        )
        .is_err()
    );
    assert!(pool.snapshot()?.reservations.is_empty());
    assert_eq!(fs::read_dir(root.path().join("owners"))?.count(), 2);
    assert_eq!(pool.recover()?.orphan_lock_files, 2);
    assert_eq!(fs::read_dir(root.path().join("owners"))?.count(), 0);
    Ok(())
}

#[test]
fn disk_promises_and_post_admission_cancellation_cover_entire_bundle() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let mut a = claim(root.path(), "a", 0, 20);
    a.request.disks.push(super::super::DiskRequest {
        path: root.path().into(),
        bytes: 60,
    });
    let mut b = claim(root.path(), "b", 1, 20);
    b.request.disks.push(super::super::DiskRequest {
        path: root.path().join("future"),
        bytes: 50,
    });
    assert!(matches!(
        pool.reserve_bundle_with_probe(limits(), vec![a, b], &|_| Ok(100))?,
        BundleAdmission::Wait(_)
    ));
    assert!(pool.snapshot()?.reservations.is_empty());
    let mut calls = 0;
    let result = pool.reserve_bundle(
        limits(),
        vec![
            claim(root.path(), "index", 0, 40),
            claim(root.path(), "run", 1, 50),
        ],
        Duration::from_secs(1),
        || {
            calls += 1;
            if calls == 2 {
                bail!("cancelled after admission");
            }
            Ok(())
        },
    );
    assert!(result.is_err());
    assert!(pool.snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn competing_bundles_never_expose_partial_capacity() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let barrier = Arc::new(Barrier::new(8));
    std::thread::scope(|scope| -> Result<()> {
        let handles = (0..8)
            .map(|i| {
                let barrier = Arc::clone(&barrier);
                let pool = pool.clone();
                let root = root.path();
                scope.spawn(move || -> Result<()> {
                    barrier.wait();
                    let held = pool.reserve_bundle(
                        limits(),
                        vec![
                            claim(root, &format!("i{i}"), 0, 40),
                            claim(root, &format!("r{i}"), 1, 50),
                        ],
                        Duration::from_secs(10),
                        || Ok(()),
                    )?;
                    let snapshot = pool.snapshot()?;
                    assert_eq!(snapshot.reservations.len(), 2);
                    assert_eq!(snapshot.memory_bytes, 90);
                    assert!(matches!(
                        pool.try_reserve(
                            limits(),
                            claim(root, "probe", 1, 20).owner,
                            ResourceRequest {
                                cpu_threads: 1,
                                memory_bytes: 20,
                                disks: Vec::new()
                            }
                        )?,
                        Admission::Wait(_)
                    ));
                    drop(held);
                    Ok(())
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            handle.join().unwrap()?;
        }
        Ok(())
    })?;
    assert!(pool.snapshot()?.reservations.is_empty());
    Ok(())
}
