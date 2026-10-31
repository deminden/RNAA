//! QC evidence evaluation, frozen cohort admission and reevaluation commands.
use crate::cli::{QcArgs, QcCommand, QcReevaluateArgs};
use crate::ingest::inputs::recover_interrupted_work;
use anyhow::{Context, Result, bail};
use rnaa_core::ProjectPaths;
use rnaa_core::{
    ArtifactKind, Database, ProjectConfig, RunRecord, RunState,
    eligibility::{self, EligibilityChange, EligibilityState, QcPolicy},
};

#[derive(Default)]
pub struct Summary {
    pub evaluated: u64,
    pub eligible: u64,
    pub rejected: u64,
    pub required: u64,
    pub invalid: u64,
    pub busy: u64,
    pub changed: u64,
}

pub fn refresh(db: &Database, config: &ProjectConfig, run: &str) -> Result<EligibilityChange> {
    let assessment = eligibility::assess(db, config, run)?;
    db.atomic(|db| {
        let change = db.record_eligibility(run, &assessment)?;
        if change.changed
            && let Some(qc) = &assessment.qc
        {
            db.record_preprocess_qc(qc)?;
        }
        Ok(change)
    })
}
pub fn refresh_all(db: &Database, config: &ProjectConfig) -> Result<Summary> {
    let mut summary = Summary::default();
    for run in db.list_active_runs()? {
        let Some(change) = db.with_idle_run(&run.run_accession, |db| {
            refresh(db, config, &run.run_accession)
        })?
        else {
            summary.busy += 1;
            continue;
        };
        summary.evaluated += 1;
        summary.changed += u64::from(change.changed);
        match change.current.descriptor.state {
            EligibilityState::Eligible => summary.eligible += 1,
            EligibilityState::Rejected => summary.rejected += 1,
            EligibilityState::PreprocessRequired => summary.required += 1,
            EligibilityState::InvalidEvidence => summary.invalid += 1,
        }
    }
    Ok(summary)
}
pub fn allowed(db: &Database, config: &ProjectConfig, run: &str) -> Result<bool> {
    if !config.quant.preprocess {
        return Ok(true);
    }
    Ok(db.eligibility(run)?.is_some_and(|value| {
        value.descriptor.policy == QcPolicy::from_config(config) && value.descriptor.eligible()
    }))
}
pub fn raw_available(db: &Database, run: &str) -> Result<bool> {
    Ok(db.list_artifacts_for_run(Some(run))?.iter().any(|file| {
        matches!(
            file.kind,
            ArtifactKind::FastqR1 | ArtifactKind::FastqR2 | ArtifactKind::FastqSingle
        ) && std::path::Path::new(&file.path).is_file()
    }))
}
pub fn needs_repair(db: &Database, config: &ProjectConfig, run: &RunRecord) -> Result<bool> {
    if !config.quant.preprocess || !raw_available(db, &run.run_accession)? {
        return Ok(false);
    }
    Ok(db.eligibility(&run.run_accession)?.is_some_and(|value| {
        value.descriptor.policy == QcPolicy::from_config(config)
            && value.descriptor.needs_preprocessing()
    }))
}
pub fn blocked_count(db: &Database, config: &ProjectConfig) -> Result<usize> {
    if !config.quant.preprocess {
        return Ok(0);
    }
    let mut blocked = 0;
    for run in db.list_active_runs()? {
        if matches!(
            run.state,
            RunState::Verified
                | RunState::PreprocessDone
                | RunState::PreprocessFailed
                | RunState::QuantFailed
        ) || eligibility::quantified(run.state)
        {
            blocked += usize::from(!allowed(db, config, &run.run_accession)?);
        }
    }
    Ok(blocked)
}

/// Freeze both admitted and excluded quantified runs. A newly quantified run or
/// an eligibility change invalidates publication even if it was not fitted.
pub struct Snapshot {
    cohort: rnaa_core::cohorts::CohortSnapshot,
    policy: QcPolicy,
    decisions: std::collections::BTreeMap<String, String>,
    abundances: std::collections::BTreeMap<String, Option<rnaa_core::fingerprint::InputIdentity>>,
    quant_heads: std::collections::BTreeMap<String, Option<String>>,
    pub eligible: std::collections::BTreeSet<String>,
    reports: Vec<rnaa_core::ManifestArtifact>,
}
impl Snapshot {
    pub fn capture(db: &Database, config: &ProjectConfig) -> Result<Self> {
        Self::capture_impl(db, config, true)
    }

    /// Observe existing eligibility without changing QC decisions or stage state.
    /// Stale or absent evidence is rejected rather than refreshed implicitly.
    pub(crate) fn observe(db: &Database, config: &ProjectConfig) -> Result<Self> {
        Self::capture_impl(db, config, false)
    }

    fn capture_impl(db: &Database, config: &ProjectConfig, refresh_evidence: bool) -> Result<Self> {
        let mut snapshot = Self {
            cohort: rnaa_core::cohorts::CohortSnapshot::capture(db)?,
            policy: QcPolicy::from_config(config),
            decisions: Default::default(),
            abundances: Default::default(),
            quant_heads: Default::default(),
            eligible: Default::default(),
            reports: Vec::new(),
        };
        for run in db
            .list_active_runs()?
            .into_iter()
            .filter(|run| eligibility::quantified(run.state))
        {
            let decision = db
                .with_idle_run(&run.run_accession, |db| {
                    if refresh_evidence {
                        Ok(refresh(db, config, &run.run_accession)?.current)
                    } else {
                        let stored = db.eligibility(&run.run_accession)?.ok_or_else(|| {
                            anyhow::anyhow!(
                                "QC eligibility is missing for {}; rerun analysis or reevaluate QC",
                                run.run_accession
                            )
                        })?;
                        let current =
                            eligibility::assess(db, config, &run.run_accession)?.descriptor;
                        if stored.descriptor != current
                            || stored.fingerprint != current.fingerprint()?
                        {
                            anyhow::bail!(
                                "QC eligibility is stale for {}; rerun analysis or reevaluate QC",
                                run.run_accession
                            );
                        }
                        Ok(stored)
                    }
                })?
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "run {} is busy; retry analysis after its worker finishes",
                        run.run_accession
                    )
                })?;
            if decision.descriptor.eligible() {
                snapshot.eligible.insert(run.run_accession.clone());
            }
            if let Some(evidence) = &decision.descriptor.evidence
                && evidence.report.matches_file().unwrap_or(false)
            {
                let mut report = evidence.report.clone();
                report.kind = format!("QC_REPORT:{}", run.run_accession);
                snapshot.reports.push(report);
            }
            snapshot.quant_heads.insert(
                run.run_accession.clone(),
                db.generation_head(
                    rnaa_core::generations::OutputGroup::Quant,
                    &run.run_accession,
                )?
                .map(|head| head.generation_id),
            );
            snapshot.abundances.insert(
                run.run_accession.clone(),
                current_abundance(db, config, &run.run_accession)?,
            );
            snapshot
                .decisions
                .insert(run.run_accession, decision.fingerprint);
        }
        snapshot.verify_catalog(db, config)?;
        Ok(snapshot)
    }
    pub fn digest(&self) -> Result<String> {
        rnaa_core::fingerprint::StageRequest {
            schema_version: 1,
            stage: "eligibility_cohort".to_string(),
            algorithm_version: 2,
            engines: Default::default(),
            parameters: serde_json::json!({"policy":self.policy,"decisions":self.decisions,"abundances":self.abundances,"active_cohort":self.cohort.digest()?}),
            inputs: Vec::new(),
            dependencies: Default::default(),
        }
        .digest()
    }
    pub fn bind<P>(
        &self,
        request: &mut rnaa_stages::analysis_request::AnalysisRequest<P>,
    ) -> Result<()> {
        request
            .request
            .dependencies
            .insert("eligibility".into(), self.digest()?);
        request
            .request
            .dependencies
            .insert("active_cohort".into(), self.cohort.digest()?);
        for report in &self.reports {
            request
                .request
                .inputs
                .push(rnaa_core::fingerprint::InputIdentity::from(report));
            request.inputs.push(report.clone());
        }
        Ok(())
    }
    /// Only catalog reads belong in the final publication transaction. Large
    /// payload hashes are checked immediately before it while heartbeats can renew.
    pub fn verify_catalog(&self, db: &Database, config: &ProjectConfig) -> Result<()> {
        self.cohort.verify_catalog(db)?;
        let persisted: ProjectConfig = serde_json::from_str(&db.project()?.config_json)?;
        if self.policy != QcPolicy::from_config(config)
            || self.policy != QcPolicy::from_config(&persisted)
        {
            anyhow::bail!("QC policy changed during analysis");
        }
        let mut decisions = std::collections::BTreeMap::new();
        let mut heads = std::collections::BTreeMap::new();
        for run in db
            .list_active_runs()?
            .into_iter()
            .filter(|run| eligibility::quantified(run.state))
        {
            let decision = db
                .eligibility(&run.run_accession)?
                .ok_or_else(|| anyhow::anyhow!("QC eligibility disappeared during analysis"))?;
            heads.insert(
                run.run_accession.clone(),
                db.generation_head(
                    rnaa_core::generations::OutputGroup::Quant,
                    &run.run_accession,
                )?
                .map(|head| head.generation_id),
            );
            decisions.insert(run.run_accession, decision.fingerprint);
        }
        if decisions != self.decisions || heads != self.quant_heads {
            anyhow::bail!(
                "quantified cohort, count generation or QC eligibility changed during analysis"
            );
        }
        Ok(())
    }
    pub fn verify_current(&self, db: &Database, config: &ProjectConfig) -> Result<()> {
        self.verify_catalog(db, config)?;
        let project = db.project()?;
        let persisted: ProjectConfig = serde_json::from_str(&project.config_json)?;
        let disk = ProjectConfig::load(&std::path::Path::new(&project.root_dir).join("rnaa.toml"))?;
        self.cohort.verify_config(config)?;
        self.cohort.verify_config(&disk)?;
        if QcPolicy::from_config(&persisted) != self.policy
            || QcPolicy::from_config(&disk) != self.policy
        {
            anyhow::bail!("QC policy changed during analysis; retry with the current policy");
        }
        let mut decisions = std::collections::BTreeMap::new();
        let mut abundances = std::collections::BTreeMap::new();
        for run in db
            .list_active_runs()?
            .into_iter()
            .filter(|run| eligibility::quantified(run.state))
        {
            let stored = db
                .eligibility(&run.run_accession)?
                .ok_or_else(|| anyhow::anyhow!("QC eligibility disappeared during analysis"))?;
            let current = eligibility::assess(db, config, &run.run_accession)?.descriptor;
            if current != stored.descriptor {
                anyhow::bail!(
                    "QC evidence changed during analysis for {}; reevaluate QC and retry",
                    run.run_accession
                );
            }
            abundances.insert(
                run.run_accession.clone(),
                current_abundance(db, config, &run.run_accession)?,
            );
            decisions.insert(run.run_accession, stored.fingerprint);
        }
        if self.policy != QcPolicy::from_config(config)
            || decisions != self.decisions
            || abundances != self.abundances
        {
            anyhow::bail!(
                "quantified cohort or QC eligibility changed during analysis; retry with the current cohort"
            );
        }
        Ok(())
    }
}

fn current_abundance(
    db: &Database,
    config: &ProjectConfig,
    run: &str,
) -> Result<Option<rnaa_core::fingerprint::InputIdentity>> {
    let paths = rnaa_core::ProjectPaths::new(db.project()?.root_dir).for_quant_run(db, run)?;
    let path = paths
        .quant_run_dir(run, &config.quant.engine)
        .join("abundance.tsv");
    if !path.exists() {
        return Ok(None);
    }
    Ok(Some(rnaa_core::fingerprint::InputIdentity::from(
        &rnaa_core::ManifestArtifact::snapshot(&format!("ABUNDANCE:{run}"), &path)?,
    )))
}

pub(crate) fn cmd_qc(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    args: QcArgs,
) -> Result<()> {
    match args.command {
        QcCommand::Reevaluate(args) => cmd_qc_reevaluate(paths, db, config, args),
    }
}

fn cmd_qc_reevaluate(
    _paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    _args: QcReevaluateArgs,
) -> Result<()> {
    recover_interrupted_work(db)?;
    db.with_idle_analysis(|db| reevaluate_qc_idle(db, config))?
        .context("QC reevaluation requires the current analysis writer to finish")
}

fn reevaluate_qc_idle(db: &Database, config: &ProjectConfig) -> Result<()> {
    let summary = refresh_all(db, config)?;
    println!("qc_reevaluated\t{}", summary.evaluated);
    println!("qc_passed\t{}", summary.eligible);
    println!(
        "qc_failed\t{}",
        summary.rejected + summary.required + summary.invalid
    );
    println!("qc_state_changes\t0");
    println!("qc_gate_outcome_changed\t{}", summary.changed);
    println!("qc_missing_reports\t{}", summary.required);
    println!("qc_invalid_evidence\t{}", summary.invalid);
    println!("qc_busy_runs\t{}", summary.busy);
    println!(
        "project_outputs_invalidated\t{}",
        if summary.changed > 0 { "yes" } else { "no" }
    );
    if summary.invalid > 0 {
        bail!(
            "QC evidence failed verification for {} run(s)",
            summary.invalid
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rnaa_core::{
        ArtifactRecord, LibraryLayout,
        util::{compute_sha256, now_rfc3339},
    };
    #[test]
    fn frozen_admission_detects_new_members_evidence_and_policy_changes_and_busy_runs() -> Result<()>
    {
        let root = tempfile::tempdir()?;
        let mut config = ProjectConfig::default();
        config.quant.preprocess = true;
        let db = Database::create(root.path(), &config)?;
        config.save(&root.path().join("rnaa.toml"))?;
        let run = RunRecord {
            project_id: db.project_id().into(),
            study_accession: "synthetic".into(),
            source_study_accession: None,
            run_accession: "RUN1".into(),
            sample_accession: None,
            experiment_accession: None,
            library_layout: LibraryLayout::Single,
            instrument_platform: None,
            instrument_model: None,
            remote_files: Vec::new(),
            metadata: Default::default(),
            state: RunState::QuantDone,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        db.upsert_run(&run)?;
        let report = root.path().join("qc.json");
        let original = r#"{"mode":"single","passed_reads":8,"failed_reads":2}"#;
        std::fs::write(&report, original)?;
        db.record_artifact(&ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some("RUN1".into()),
            kind: ArtifactKind::PreprocessReport,
            path: report.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".into(),
            checksum: compute_sha256(&report)?,
            bytes: original.len() as u64,
            created_at: now_rfc3339(),
        })?;
        assert!(Snapshot::observe(&db, &config).is_err());
        assert!(db.eligibility("RUN1")?.is_none());
        let frozen = Snapshot::capture(&db, &config)?;
        assert_eq!(Snapshot::observe(&db, &config)?.digest()?, frozen.digest()?);
        frozen.verify_current(&db, &config)?;
        let abundance = rnaa_core::ProjectPaths::new(root.path())
            .quant_run_dir("RUN1", "kallistors")
            .join("abundance.tsv");
        std::fs::create_dir_all(abundance.parent().unwrap())?;
        std::fs::write(&abundance, "replacement counts")?;
        assert!(
            frozen.verify_current(&db, &config).is_err(),
            "a new count result must not publish an analysis of the previous counts"
        );
        std::fs::remove_file(abundance)?;
        db.set_project_stage_state("deseq2", rnaa_core::ProjectStageState::Running, None)?;
        frozen.verify_current(&db, &config)?;
        db.set_project_stage_state("deseq2", rnaa_core::ProjectStageState::Done, None)?;
        std::fs::write(&report, "corrupt")?;
        assert!(frozen.verify_current(&db, &config).is_err());
        let prior = serde_json::to_value(db.eligibility("RUN1")?)?;
        assert!(Snapshot::observe(&db, &config).is_err());
        assert_eq!(serde_json::to_value(db.eligibility("RUN1")?)?, prior);
        std::fs::write(&report, original)?;
        frozen.verify_current(&db, &config)?;
        config.quant.qc_gate.min_pass_rate = Some(0.95);
        config.save(&root.path().join("rnaa.toml"))?;
        assert!(frozen.verify_current(&db, &config).is_err());
        config.quant.qc_gate.min_pass_rate = Some(0.6);
        config.save(&root.path().join("rnaa.toml"))?;
        frozen.verify_current(&db, &config)?;
        let mut lease = rnaa_core::jobs::JobLease::try_acquire(
            &db,
            "retention",
            Some("RUN1"),
            &[RunState::QuantDone],
        )?
        .unwrap();
        assert!(Snapshot::capture(&db, &config).is_err());
        assert!(Snapshot::observe(&db, &config).is_err());
        lease.finish(rnaa_core::jobs::JobOutcome::Succeeded, None)?;
        drop(lease);
        let mut added = run;
        added.run_accession = "RUN2".into();
        db.upsert_run(&added)?;
        assert!(frozen.verify_current(&db, &config).is_err());
        Ok(())
    }
}
