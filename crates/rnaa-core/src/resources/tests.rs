use super::*;
use anyhow::{Result, bail};
use std::{
    fs,
    path::Path,
    process::{Child, Command, Stdio},
    sync::{
        Arc, Barrier,
        atomic::{AtomicUsize, Ordering},
    },
    time::{Duration, Instant},
};

fn limits(cpu: u64, memory: u64) -> PoolLimits {
    PoolLimits {
        cpu_threads: cpu,
        memory_bytes: memory,
        disk_headroom_bytes: 10,
    }
}
fn owner(root: &Path, label: &str) -> ReservationOwner {
    ReservationOwner {
        project_root: root.to_path_buf(),
        attempt_id: label.into(),
        stage: "quant".into(),
        run_accession: None,
    }
}
fn request(cpu: u64, memory: u64) -> ResourceRequest {
    ResourceRequest {
        cpu_threads: cpu,
        memory_bytes: memory,
        disks: vec![],
    }
}
fn granted(value: Admission) -> ResourceReservation {
    match value {
        Admission::Granted(value) => value,
        Admission::Wait(wait) => panic!("unexpected admission wait: {wait:?}"),
    }
}
fn waiting(value: Admission, resource: &str) -> AdmissionWait {
    match value {
        Admission::Granted(_) => panic!("unexpected admission"),
        Admission::Wait(wait) => {
            assert_eq!(wait.resource, resource);
            wait
        }
    }
}

#[test]
fn cpu_memory_and_disk_are_admitted_atomically_across_pool_handles() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let second = ResourcePool::open(root.path())?;
    let first =
        granted(pool.try_reserve(limits(4, 100), owner(root.path(), "first"), request(2, 80))?);
    let wait = waiting(
        second.try_reserve(limits(4, 100), owner(root.path(), "memory"), request(1, 30))?,
        "memory",
    );
    assert_eq!(wait.available, 20);
    assert_eq!(pool.snapshot()?.reservations.len(), 1);
    let other =
        granted(second.try_reserve(limits(4, 100), owner(root.path(), "other"), request(2, 20))?);
    waiting(
        pool.try_reserve(limits(4, 100), owner(root.path(), "cpu"), request(1, 1))?,
        "CPU",
    );
    drop(other);
    drop(first);
    assert!(pool.snapshot()?.reservations.is_empty());
    assert_eq!(pool.recover()?.orphan_lock_files, 0);
    let disk = |bytes| ResourceRequest {
        disks: vec![DiskRequest {
            path: root.path().join("not-created/output"),
            bytes,
        }],
        ..request(0, 1)
    };
    let held = granted(pool.reserve_with_probe(
        limits(4, 100),
        owner(root.path(), "disk-a"),
        disk(60),
        &|_| Ok(100),
    )?);
    let mut next = disk(40);
    next.disks[0].path = root.path().join("another/path");
    let wait = waiting(
        second.reserve_with_probe(limits(4, 100), owner(root.path(), "disk-b"), next, &|_| {
            Ok(100)
        })?,
        "disk",
    );
    assert_eq!(wait.available, 30);
    assert_eq!(pool.snapshot()?.pending_disk_bytes, 60);
    drop(held);
    assert!(super::available_disk_bytes(&root.path().join("future/data"))? > 0);
    Ok(())
}

#[test]
fn restrictive_live_limits_cannot_be_relaxed_by_another_project() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let first =
        granted(pool.try_reserve(limits(2, 100), owner(root.path(), "small"), request(1, 50))?);
    waiting(
        pool.try_reserve(limits(8, 1000), owner(root.path(), "large"), request(2, 50))?,
        "CPU",
    );
    let second = granted(pool.try_reserve(
        limits(8, 1000),
        owner(root.path(), "remaining"),
        request(1, 50),
    )?);
    assert_eq!(pool.snapshot()?.effective_limits, Some(limits(2, 100)));
    drop(first);
    drop(second);
    let large = granted(pool.try_reserve(
        limits(8, 1000),
        owner(root.path(), "large"),
        request(8, 1000),
    )?);
    assert_eq!(pool.snapshot()?.cpu_threads, 8);
    assert!(
        pool.try_reserve(
            limits(2, 100),
            owner(root.path(), "impossible"),
            request(3, 1)
        )
        .is_err()
    );
    drop(large);
    Ok(())
}

#[test]
fn zero_requested_dimensions_cannot_retroactively_reduce_capacity_below_live_usage() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let held =
        granted(pool.try_reserve(limits(4, 100), owner(root.path(), "held"), request(4, 80))?);
    waiting(
        pool.try_reserve(limits(2, 100), owner(root.path(), "no-cpu"), request(0, 1))?,
        "CPU",
    );
    drop(held);
    let held =
        granted(pool.try_reserve(limits(4, 100), owner(root.path(), "held"), request(1, 80))?);
    waiting(
        pool.try_reserve(
            limits(4, 50),
            owner(root.path(), "no-memory"),
            request(1, 0),
        )?,
        "memory",
    );
    assert_eq!(pool.snapshot()?.effective_limits, Some(limits(4, 100)));
    drop(held);
    Ok(())
}

#[test]
fn resize_preserves_old_reservation_on_failure_and_releases_unused_capacity() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let mut a =
        granted(pool.try_reserve(limits(4, 100), owner(root.path(), "a"), request(1, 70))?);
    let b = granted(pool.try_reserve(limits(4, 100), owner(root.path(), "b"), request(1, 20))?);
    let original = a.snapshot().clone();
    assert_eq!(a.try_resize(request(2, 90))?.unwrap().resource, "memory");
    assert_eq!(a.snapshot(), &original);
    assert_eq!(pool.snapshot()?.memory_bytes, 90);
    assert!(a.try_resize(request(1, 20))?.is_none());
    let c = granted(pool.try_reserve(limits(4, 100), owner(root.path(), "c"), request(2, 60))?);
    assert_eq!(pool.snapshot()?.memory_bytes, 100);
    c.release()?;
    b.release()?;
    a.release()?;
    assert!(pool.snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn competing_threads_cannot_overbook_memory_even_when_cpu_is_available() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let barrier = Arc::new(Barrier::new(8));
    let active = AtomicUsize::new(0);
    let peak = AtomicUsize::new(0);
    std::thread::scope(|scope| -> Result<()> {
        let workers = (0..8)
            .map(|worker| {
                let pool = &pool;
                let barrier = &barrier;
                let active = &active;
                let peak = &peak;
                let path = root.path();
                scope.spawn(move || -> Result<()> {
                    barrier.wait();
                    for iteration in 0..3 {
                        let held = pool.reserve(
                            limits(8, 128),
                            owner(path, &format!("{worker}-{iteration}")),
                            request(1, 64),
                            Duration::from_secs(10),
                            || Ok(()),
                        )?;
                        let now = active.fetch_add(1, Ordering::SeqCst) + 1;
                        peak.fetch_max(now, Ordering::SeqCst);
                        assert!(now <= 2);
                        let snapshot = pool.snapshot()?;
                        assert!(snapshot.cpu_threads <= 2 && snapshot.memory_bytes <= 128);
                        std::thread::sleep(Duration::from_millis(10));
                        active.fetch_sub(1, Ordering::SeqCst);
                        drop(held);
                    }
                    Ok(())
                })
            })
            .collect::<Vec<_>>();
        for worker in workers {
            worker.join().unwrap()?;
        }
        Ok(())
    })?;
    assert!(peak.load(Ordering::SeqCst) > 0);
    assert!(pool.snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn waiting_checks_cancellation_and_timeout_without_publishing_a_reservation() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let _held =
        granted(pool.try_reserve(limits(1, 100), owner(root.path(), "held"), request(1, 1))?);
    let mut checks = 0;
    let cancelled = pool.reserve(
        limits(1, 100),
        owner(root.path(), "cancelled"),
        request(1, 1),
        Duration::from_secs(5),
        || {
            checks += 1;
            if checks == 2 {
                bail!("cancelled by owner");
            }
            Ok(())
        },
    );
    assert!(
        cancelled
            .err()
            .unwrap()
            .to_string()
            .contains("cancelled by owner")
    );
    assert_eq!(checks, 2);
    assert!(
        pool.reserve(
            limits(1, 100),
            owner(root.path(), "timeout"),
            request(1, 1),
            Duration::ZERO,
            || Ok(())
        )
        .err()
        .unwrap()
        .to_string()
        .contains("timed out")
    );
    assert_eq!(pool.snapshot()?.reservations.len(), 1);
    Ok(())
}

struct Worker(Child);
impl Drop for Worker {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

#[test]
fn reservation_child() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_RESERVATION_TEST_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let pool = ResourcePool::open(root)?;
    if std::env::var_os("RNAA_RESERVATION_CONTENTION").is_some() {
        let label = std::process::id().to_string();
        fs::write(root.join(format!("started-{label}")), [])?;
        while !root.join("start-gate").exists() {
            std::thread::sleep(Duration::from_millis(10));
        }
        let _held = if std::env::var("RNAA_RESERVATION_BUNDLE").as_deref() == Ok("1") {
            pool.reserve_bundle(
                limits(8, 128),
                vec![
                    ReservationClaim {
                        owner: owner(root, &label),
                        request: request(0, 24),
                    },
                    ReservationClaim {
                        owner: owner(root, &label),
                        request: request(1, 40),
                    },
                ],
                Duration::from_secs(10),
                || Ok(()),
            )?
        } else {
            vec![pool.reserve(
                limits(8, 128),
                owner(root, &label),
                request(1, 64),
                Duration::from_secs(10),
                || Ok(()),
            )?]
        };
        fs::write(root.join(format!("claimed-{label}")), [])?;
        while !root.join(format!("release-{label}")).exists() {
            std::thread::sleep(Duration::from_millis(10));
        }
        return Ok(());
    }
    let held = granted(pool.try_reserve(limits(2, 100), owner(root, "child"), request(2, 80))?);
    crate::util::publish_file(&root.join("ready"), |staged| {
        fs::write(staged, &held.snapshot().reservation_id)?;
        Ok(())
    })?;
    loop {
        std::thread::sleep(Duration::from_millis(50));
    }
}

#[test]
fn competing_processes_share_one_memory_budget_through_multiple_waves() -> Result<()> {
    process_contention(false)
}
#[test]
fn competing_bundled_processes_share_one_budget_without_partial_admission() -> Result<()> {
    process_contention(true)
}
fn process_contention(bundled: bool) -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let mut workers = (0..6)
        .map(|_| {
            Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "resources::tests::reservation_child",
                    "--nocapture",
                ])
                .env("RNAA_RESERVATION_TEST_ROOT", root.path())
                .env("RNAA_RESERVATION_CONTENTION", "1")
                .env("RNAA_RESERVATION_BUNDLE", if bundled { "1" } else { "0" })
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .spawn()
                .map(Worker)
        })
        .collect::<std::io::Result<Vec<_>>>()?;
    let deadline = Instant::now() + Duration::from_secs(15);
    while fs::read_dir(root.path())?
        .filter_map(|entry| entry.ok())
        .filter(|entry| entry.file_name().to_string_lossy().starts_with("started-"))
        .count()
        != 6
    {
        assert!(Instant::now() < deadline, "workers did not start");
        std::thread::sleep(Duration::from_millis(10));
    }
    fs::write(root.path().join("start-gate"), [])?;
    let mut released = std::collections::BTreeSet::new();
    while released.len() < 6 {
        assert!(Instant::now() < deadline, "workers did not pass admission");
        let snapshot = pool.snapshot()?;
        assert!(
            snapshot.reservations.len() <= if bundled { 4 } else { 2 }
                && snapshot.memory_bytes <= 128
        );
        let mut fresh = std::collections::BTreeMap::new();
        for row in snapshot.reservations {
            if !released.contains(&row.owner.attempt_id) {
                *fresh.entry(row.owner.attempt_id).or_insert(0) += 1;
            }
        }
        assert!(
            fresh
                .values()
                .all(|&rows| rows == if bundled { 2 } else { 1 }),
            "new owner exposed only part of its bundle"
        );
        if fresh.len() == 2 || fresh.len() + released.len() == 6 {
            for (id, _) in fresh {
                released.insert(id.clone());
                fs::write(root.path().join(format!("release-{id}")), [])?;
            }
        }
        std::thread::sleep(Duration::from_millis(10));
    }
    for worker in &mut workers {
        assert!(worker.0.wait()?.success());
    }
    assert!(pool.snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn a_live_process_keeps_its_reservation_and_killing_it_allows_recovery() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let mut child = Worker(
        Command::new(std::env::current_exe()?)
            .args([
                "--exact",
                "resources::tests::reservation_child",
                "--nocapture",
            ])
            .env("RNAA_RESERVATION_TEST_ROOT", root.path())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()?,
    );
    let deadline = Instant::now() + Duration::from_secs(10);
    while !root.path().join("ready").exists() {
        assert!(
            child.0.try_wait()?.is_none(),
            "reservation child exited before admission"
        );
        assert!(
            Instant::now() < deadline,
            "reservation child never claimed capacity"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
    let old = fs::read_to_string(root.path().join("ready"))?;
    assert!(
        !old.is_empty(),
        "reservation child published an empty identity"
    );
    for _ in 0..3 {
        waiting(
            pool.try_reserve(limits(2, 100), owner(root.path(), "waiting"), request(1, 1))?,
            "CPU",
        );
        std::thread::sleep(Duration::from_millis(30));
    }
    assert_eq!(pool.snapshot()?.reservations[0].reservation_id, old);
    child.0.kill()?;
    child.0.wait()?;
    let replacement = granted(pool.try_reserve(
        limits(2, 100),
        owner(root.path(), "replacement"),
        request(2, 100),
    )?);
    assert_ne!(replacement.snapshot().reservation_id, old);
    assert_eq!(pool.snapshot()?.reservations.len(), 1);
    assert!(!pool.owner_path(&old)?.exists());
    Ok(())
}

#[test]
fn catalog_failures_roll_back_and_unknown_or_missing_owner_proof_is_held() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    {
        let catalog = pool.catalog()?;
        catalog.connection.execute_batch("CREATE TRIGGER reject_admission BEFORE INSERT ON reservations BEGIN SELECT RAISE(ABORT,'injected reservation failure'); END;")?;
    }
    assert!(
        pool.try_reserve(limits(2, 100), owner(root.path(), "failed"), request(1, 1))
            .is_err()
    );
    assert!(pool.snapshot()?.reservations.is_empty());
    assert_eq!(pool.recover()?.orphan_lock_files, 1);
    {
        let catalog = pool.catalog()?;
        catalog
            .connection
            .execute_batch("DROP TRIGGER reject_admission")?;
    }
    let held =
        granted(pool.try_reserve(limits(2, 100), owner(root.path(), "held"), request(1, 1))?);
    let id = held.snapshot().reservation_id.clone();
    fs::remove_file(pool.owner_path(&id)?)?;
    assert!(
        pool.snapshot()
            .err()
            .unwrap()
            .to_string()
            .contains("cannot establish ownership")
    );
    drop(held);
    // Explicit release can retire its own unchanged row even if the empty lock
    // pathname was removed. No other claim may infer that absence means death.
    assert!(pool.snapshot()?.reservations.is_empty());
    {
        let catalog = pool.catalog()?;
        catalog
            .connection
            .execute("UPDATE pool_schema SET version=2", [])?;
    }
    assert!(ResourcePool::open(root.path()).is_err());
    Ok(())
}

#[test]
fn replacing_a_live_owner_lock_cannot_reclaim_its_capacity() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let held =
        granted(pool.try_reserve(limits(1, 100), owner(root.path(), "live"), request(1, 100))?);
    let path = pool.owner_path(&held.snapshot().reservation_id)?;
    fs::remove_file(&path)?;
    fs::write(&path, [])?;
    assert!(
        pool.snapshot()
            .err()
            .unwrap()
            .to_string()
            .contains("cannot establish ownership")
    );
    assert!(
        pool.try_reserve(limits(1, 100), owner(root.path(), "other"), request(1, 100))
            .is_err()
    );
    assert_eq!(held.snapshot().requested.memory_bytes, 100);
    drop(held);
    assert!(
        path.exists(),
        "release must leave a substituted owner file alone"
    );
    assert!(pool.snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn cancellation_after_admission_releases_the_newly_claimed_capacity() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let mut checks = 0;
    let result = pool.reserve(
        limits(1, 100),
        owner(root.path(), "cancel"),
        request(1, 100),
        Duration::from_secs(1),
        || {
            checks += 1;
            if checks == 2 {
                bail!("cancelled while entering admission");
            }
            Ok(())
        },
    );
    assert!(
        result
            .err()
            .unwrap()
            .to_string()
            .contains("cancelled while entering admission")
    );
    assert!(pool.snapshot()?.reservations.is_empty());
    Ok(())
}

#[test]
fn disk_growth_is_reserved_before_writes_and_shrinking_releases_future_bytes() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let disk = |bytes| ResourceRequest {
        disks: vec![DiskRequest {
            path: root.path().to_path_buf(),
            bytes,
        }],
        ..request(0, 1)
    };
    let mut held =
        granted(pool.try_reserve(limits(1, 100), owner(root.path(), "disk"), disk(60))?);
    assert_eq!(held.try_resize(disk(u64::MAX))?.unwrap().resource, "disk");
    assert_eq!(held.snapshot().requested.disk_bytes()?, 60);
    assert!(held.try_resize(disk(16))?.is_none());
    let other = granted(pool.reserve_with_probe(
        limits(1, 100),
        owner(root.path(), "other"),
        disk(32),
        &|_| Ok(100),
    )?);
    assert_eq!(pool.snapshot()?.pending_disk_bytes, 48);
    drop(other);
    drop(held);
    Ok(())
}

#[test]
fn orphan_recovery_holds_live_locks_unknown_names_and_nonempty_files() -> Result<()> {
    let root = tempfile::tempdir()?;
    let pool = ResourcePool::open(root.path())?;
    let held =
        granted(pool.try_reserve(limits(1, 100), owner(root.path(), "held"), request(1, 10))?);
    let orphan = pool.owner_path(&uuid::Uuid::new_v4().to_string())?;
    let locked = fs::File::create(&orphan)?;
    locked.try_lock()?;
    let nonempty = pool.owner_path(&uuid::Uuid::new_v4().to_string())?;
    fs::write(&nonempty, b"retained diagnostics")?;
    let unknown = pool.root().join("owners/unknown");
    fs::write(&unknown, [])?;
    assert_eq!(pool.recover()?.orphan_lock_files, 0);
    drop(locked);
    assert_eq!(pool.recover()?.orphan_lock_files, 1);
    assert!(nonempty.is_file() && unknown.is_file());
    assert_eq!(pool.snapshot()?.reservations.len(), 1);
    drop(held);
    Ok(())
}
