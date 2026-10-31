use super::*;

fn publish(db: &Database, value: &str) -> Result<GenerationHead> {
    let (mut lease, mut generation, records) = preparing(db, value)?;
    generation.seal(
        &records
            .iter()
            .map(|r| PathBuf::from(&r.path))
            .collect::<Vec<_>>(),
        &"a".repeat(64),
    )?;
    let head = generation.head().clone();
    generation.commit(&records, |_| Ok(()))?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(head)
}

fn registry(db: &Database) -> Result<serde_json::Value> {
    Ok(serde_json::to_value(db.list_artifacts()?)?)
}

#[test]
fn historical_selection_restores_registry_without_rewriting_producer_or_payload() -> Result<()> {
    let (_temp, db) = setup()?;
    let first = publish(&db, "first")?;
    let witness = first.payload_path.parent().unwrap().join("seal.json");
    let original = fs::read(&witness)?;
    let original_registry = registry(&db)?;
    let second = publish(&db, "second")?;
    let candidates =
        db.generation_candidates(OutputGroup::De, "deseq2", db.project_id(), &"a".repeat(64))?;
    assert_eq!(candidates[0], second);
    assert!(candidates.contains(&first));
    for (group, stage, subject, request) in [
        (OutputGroup::Corr, "deseq2", db.project_id(), "a".repeat(64)),
        (
            OutputGroup::De,
            "normalize",
            db.project_id(),
            "a".repeat(64),
        ),
        (OutputGroup::De, "deseq2", "another subject", "a".repeat(64)),
        (OutputGroup::De, "deseq2", db.project_id(), "b".repeat(64)),
    ] {
        assert!(
            db.generation_candidates(group, stage, subject, &request)?
                .is_empty()
        );
    }
    let pinned_second = second.clone();
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&first.request_sha256)?;
    let verified = owned.prepare_generation_reuse(&first)?.unwrap();
    assert_eq!(verified.head(), &first);
    assert_eq!(verified.artifact_count(), 2);
    owned.reattach_generation(verified, |db| {
        db.set_project_stage_state("deseq2", ProjectStageState::Done, None)
    })?;
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .unwrap(),
        first
    );
    assert_eq!(registry(&db)?, original_registry);
    assert_eq!(fs::read(witness)?, original);
    assert!(db.verify_generation_head(&pinned_second)?);
    let event = db
        .list_events(100)?
        .into_iter()
        .find(|event| event.message == "selected verified committed generation")
        .unwrap();
    assert_eq!(event.context["producer_attempt_id"], first.attempt_id);
    assert_ne!(event.context["reuse_attempt_id"], first.attempt_id);
    assert_eq!(
        event.context["previous_generation_id"],
        second.generation_id
    );
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(())
}

#[test]
fn callback_failure_rolls_back_head_registry_events_state_and_attempt_success() -> Result<()> {
    let (_temp, db) = setup()?;
    let first = publish(&db, "first")?;
    let second = publish(&db, "second")?;
    let before = registry(&db)?;
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&first.request_sha256)?;
    owned.set_project_stage_state("deseq2", ProjectStageState::Running, None)?;
    let events = db.list_events(100)?.len();
    let verified = owned.prepare_generation_reuse(&first)?.unwrap();
    let result: Result<()> = owned.reattach_generation(verified, |db| {
        db.set_project_stage_state("deseq2", ProjectStageState::Done, None)?;
        db.append_event("deseq2", None, "must roll back", serde_json::json!({}))?;
        bail!("injected failure after reattachment")
    });
    assert!(result.is_err());
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .unwrap(),
        second
    );
    assert_eq!(registry(&db)?, before);
    assert_eq!(db.list_events(100)?.len(), events);
    owned.check_job()?;
    assert_eq!(
        db.get_project_stage_state("deseq2")?.unwrap().state,
        ProjectStageState::Running
    );
    lease.finish(JobOutcome::Failed, None)?;
    Ok(())
}

#[test]
fn stale_selection_and_changed_catalog_cannot_reattach_a_verified_candidate() -> Result<()> {
    let (_temp, db) = setup()?;
    let first = publish(&db, "first")?;
    let verified = db.prepare_generation_reuse(&first)?.unwrap();
    let second = publish(&db, "second")?;
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&first.request_sha256)?;
    assert!(owned.reattach_generation(verified, |_| Ok(())).is_err());
    let verified = owned.prepare_generation_reuse(&first)?.unwrap();
    db.write(|conn| {
        conn.execute(
            "UPDATE artifact_generations SET artifacts_json='[]' WHERE generation_id=?1",
            [&first.generation_id],
        )?;
        Ok(())
    })?;
    assert!(owned.reattach_generation(verified, |_| Ok(())).is_err());
    assert!(owned.prepare_generation_reuse(&first)?.is_none());
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .unwrap(),
        second
    );
    lease.finish(JobOutcome::Failed, None)?;
    Ok(())
}

#[test]
fn corrupt_payloads_and_changed_payloads_after_preparation_are_rejected() -> Result<()> {
    let (_temp, db) = setup()?;
    let first = publish(&db, "first")?;
    let verified = db.prepare_generation_reuse(&first)?.unwrap();
    let second = publish(&db, "second")?;
    fs::write(first.payload_path.join("counts.tsv"), "corruption")?;
    assert!(db.prepare_generation_reuse(&first)?.is_none());
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&first.request_sha256)?;
    let error = owned.reattach_generation(verified, |_| Ok(())).unwrap_err();
    assert!(error.to_string().contains("payload changed"));
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .unwrap(),
        second
    );
    lease.finish(JobOutcome::Failed, None)?;
    Ok(())
}

#[test]
fn wrong_attempt_request_stage_project_and_cancelled_owner_cannot_select_history() -> Result<()> {
    let (_temp, db) = setup()?;
    let first = publish(&db, "first")?;
    assert!(
        db.reattach_generation(db.prepare_generation_reuse(&first)?.unwrap(), |_| Ok(()))
            .is_err()
    );
    for stage in ["deseq2", "normalize"] {
        let mut lease = JobLease::try_acquire(&db, stage, None, &[])?.unwrap();
        let owned = lease.database();
        owned.set_job_request(&if stage == "deseq2" {
            "b".repeat(64)
        } else {
            first.request_sha256.clone()
        })?;
        assert!(
            owned
                .reattach_generation(owned.prepare_generation_reuse(&first)?.unwrap(), |_| Ok(()))
                .is_err()
        );
        lease.finish(JobOutcome::Failed, None)?;
    }
    let (_other_temp, other) = setup()?;
    let verified = db.prepare_generation_reuse(&first)?.unwrap();
    assert!(other.reattach_generation(verified, |_| Ok(())).is_err());
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&first.request_sha256)?;
    let verified = owned.prepare_generation_reuse(&first)?.unwrap();
    let attempts = db.list_job_attempts(1)?;
    db.cancel_job(&attempts[0].attempt_id)?;
    assert!(owned.reattach_generation(verified, |_| Ok(())).is_err());
    lease.finish(JobOutcome::Interrupted, None)?;
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .unwrap(),
        first
    );
    Ok(())
}

#[test]
fn selection_without_a_current_head_retires_legacy_registration_and_keeps_its_file() -> Result<()> {
    let (temp, db) = setup()?;
    let first = publish(&db, "first")?;
    let original = registry(&db)?;
    let legacy = temp.path().join("de/legacy.tsv");
    fs::create_dir_all(legacy.parent().unwrap())?;
    fs::write(&legacy, "old legacy payload")?;
    let mut record = db.list_artifacts()?.remove(0);
    record.path = legacy.to_string_lossy().into_owned();
    record.checksum = compute_sha256(&legacy)?;
    record.bytes = fs::metadata(&legacy)?.len();
    db.record_artifact(&record)?;
    db.write(|conn| {
        conn.execute(
            "DELETE FROM generation_heads WHERE project_id=?1 AND output_group='de'",
            [db.project_id()],
        )?;
        Ok(())
    })?;
    let mut lease = JobLease::try_acquire(&db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&first.request_sha256)?;
    let verified = owned.prepare_generation_reuse(&first)?.unwrap();
    owned.reattach_generation(verified, |_| Ok(()))?;
    assert_eq!(registry(&db)?, original);
    assert_eq!(fs::read_to_string(legacy)?, "old legacy payload");
    assert_eq!(
        db.generation_head(OutputGroup::De, db.project_id())?
            .unwrap(),
        first
    );
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(())
}
