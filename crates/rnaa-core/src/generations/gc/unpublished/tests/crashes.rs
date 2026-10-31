use super::*;
use std::process::{Child, Command, Stdio};
struct KillOnDrop(Child);
impl Drop for KillOnDrop {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

fn collect_when_idle(db: &Database, opts: GenerationGcOptions) -> Result<UnpublishedCleanupReport> {
    let start = std::time::Instant::now();
    loop {
        let report = db.unpublished_generation_gc(opts)?;
        if !report.busy {
            return Ok(report);
        }
        assert!(
            start.elapsed() < std::time::Duration::from_secs(5),
            "unpublished cleanup remained busy"
        );
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
}
#[test]
fn crash_worker() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_UNPUBLISHED_CRASH_ROOT") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let db = Database::open(root)?;
    let phase = std::env::var("RNAA_UNPUBLISHED_CRASH_PHASE")?;
    let id = std::env::var("RNAA_UNPUBLISHED_CRASH_ID")?;
    let _worker;
    let _generation;
    let _maintenance;
    let ready_id;
    if phase == "building" {
        let worker = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
        let generation = Generation::begin(&worker.database(), OutputGroup::De, "deseq2")?;
        fs::write(
            generation.head().payload_path.join("partial.tsv"),
            "unfinished live worker",
        )?;
        ready_id = generation.head().generation_id.clone();
        _worker = Some(worker);
        _generation = Some(generation);
        _maintenance = None;
    } else {
        _worker = None;
        _generation = None;
        _maintenance = Some(crate::generations::pins::maintenance_lock(&db)?.unwrap());
        let mut db = db.clone();
        db.generation_maintenance = true;
        ready_id = id.clone();
        if phase.starts_with("trash") {
            let mut row = GenerationRow::unpublished(&db)?[&id].clone();
            let original_state = row.state.clone();
            row.state = "ABANDONED".into();
            let now = now_ms()?;
            let plan = Descriptor {
                schema_version: 1,
                operation_id: uuid::Uuid::new_v4().to_string(),
                generation: row,
                original_state,
                tree: tree::observe(&tree::container(&db, &id, false, false)?, false)?,
                trashed_at_ms: now,
                purge_after_ms: now,
            };
            db.mark_unpublished(&[plan])?;
            if phase == "trash_moved" {
                tree::sync_move(
                    &tree::container(&db, &id, false, false)?,
                    &tree::container(&db, &id, true, true)?,
                )?;
            }
        } else if phase.starts_with("restore") {
            let row = &Receipt::all(&db)?[&id];
            row.transition(&db, "RESTORING")?;
            if phase == "restore_moved" {
                tree::sync_move(
                    &tree::container(&db, &id, true, false)?,
                    &tree::container(&db, &id, false, true)?,
                )?;
            }
        } else if phase == "purge_partial" {
            let row = &Receipt::all(&db)?[&id];
            row.transition(&db, "PURGING")?;
            fs::remove_file(tree::container(&db, &id, true, false)?.join("payload/partial.tsv"))?;
        } else {
            bail!("unknown unpublished crash phase");
        }
    }
    crate::util::publish_file(&root.join("unpublished-ready"), |staged| {
        fs::write(staged, &ready_id)?;
        Ok(())
    })?;
    loop {
        std::thread::sleep(std::time::Duration::from_secs(1));
    }
}
#[test]
fn actual_killed_builders_and_maintenance_resume_only_durable_owned_intents() -> Result<()> {
    for phase in [
        "building",
        "trash_marked",
        "trash_moved",
        "restore_marked",
        "restore_moved",
        "purge_partial",
    ] {
        let (root, db) = setup()?;
        let initial = if phase == "building" {
            None
        } else {
            Some(abandoned(&db, true)?)
        };
        if phase.starts_with("restore") || phase == "purge_partial" {
            let mut opts = options(true);
            opts.trash_days = 0;
            let report = collect_when_idle(&db, opts)?;
            assert!(report.errors.is_empty(), "{phase}: {:?}", report.errors);
            assert!(
                Receipt::all(&db)?.contains_key(&initial.as_ref().unwrap().generation_id),
                "{phase}: cleanup did not create its receipt"
            );
        }
        let mut child = KillOnDrop(
            Command::new(std::env::current_exe()?)
                .args([
                    "--exact",
                    "generations::gc::unpublished::tests::crashes::crash_worker",
                    "--nocapture",
                ])
                .env("RNAA_UNPUBLISHED_CRASH_ROOT", root.path())
                .env("RNAA_UNPUBLISHED_CRASH_PHASE", phase)
                .env(
                    "RNAA_UNPUBLISHED_CRASH_ID",
                    initial
                        .as_ref()
                        .map_or("", |head| head.generation_id.as_str()),
                )
                .stdout(Stdio::null())
                .stderr(Stdio::piped())
                .spawn()?,
        );
        let started = std::time::Instant::now();
        while !root.path().join("unpublished-ready").exists() {
            assert!(child.0.try_wait()?.is_none(), "worker exited at {phase}");
            assert!(
                started.elapsed() < std::time::Duration::from_secs(15),
                "worker not ready at {phase}"
            );
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        let id = fs::read_to_string(root.path().join("unpublished-ready"))?;
        uuid::Uuid::parse_str(&id)?;
        if let Some(initial) = &initial {
            assert_eq!(
                id, initial.generation_id,
                "wrong ready generation at {phase}"
            );
        }
        assert!(db.unpublished_generation_gc(options(false))?.busy);
        child.0.kill()?;
        child.0.wait()?;
        if phase == "building" {
            let report = collect_when_idle(&db, options(false))?;
            assert_eq!(
                report
                    .entries
                    .iter()
                    .find(|row| row.generation_id == id)
                    .unwrap()
                    .action,
                "held_invalid"
            );
            db.write(|conn| {
                conn.execute(
                    "UPDATE job_attempts SET lease_expires_ms=0 WHERE state='RUNNING'",
                    [],
                )?;
                Ok(())
            })?;
            db.recover_jobs()?;
            db.recover_generations()?;
        }
        if phase.starts_with("restore") {
            db.restore_unpublished_generation(&id)?;
            assert!(Receipt::all(&db)?.is_empty());
            assert!(
                tree::container(&db, &id, false, false)?
                    .join("payload/partial.tsv")
                    .is_file()
            );
            assert_eq!(GenerationRow::unpublished(&db)?[&id].state, "ABANDONED");
        } else {
            let report = collect_when_idle(&db, options(true))?;
            assert!(report.errors.is_empty(), "{phase}: {:?}", report.errors);
            let result = report
                .entries
                .iter()
                .find(|entry| entry.generation_id == id)
                .unwrap_or_else(|| {
                    panic!(
                        "{phase}: missing {id} from report {}",
                        serde_json::to_string(&report).unwrap()
                    )
                });
            assert_eq!(
                result.state,
                if phase == "purge_partial" {
                    "PURGED"
                } else {
                    "TRASHED"
                }
            );
        }
        assert!(db.generation_heads()?.is_empty());
    }
    Ok(())
}
