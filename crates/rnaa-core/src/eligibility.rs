//! QC eligibility is a versioned decision about verified evidence. It never
//! rewrites ingest/quantification progress or published scientific outputs.
use crate::{
    ArtifactKind, Database, ManifestArtifact, PreprocessQcRecord, ProjectConfig, RunState,
    fingerprint::{InputIdentity, StageRequest},
    generations::OutputGroup,
    util::now_rfc3339,
};
use anyhow::{Context, Result, bail};
use rusqlite::{Connection, OptionalExtension, params};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, fs, path::Path};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct QcPolicy {
    pub enabled: bool,
    pub strict: Option<bool>,
    pub gate: Option<crate::PreprocessQcGate>,
}
impl QcPolicy {
    pub fn from_config(config: &ProjectConfig) -> Self {
        Self {
            enabled: config.quant.preprocess,
            strict: config
                .quant
                .preprocess
                .then_some(config.quant.preprocess_strict),
            gate: config
                .quant
                .preprocess
                .then(|| config.quant.qc_gate.clone()),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum EligibilityState {
    Eligible,
    Rejected,
    PreprocessRequired,
    InvalidEvidence,
}
impl EligibilityState {
    pub fn eligible(self) -> bool {
        self == Self::Eligible
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Evidence {
    pub report: ManifestArtifact,
    pub producer_generation: Option<String>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Eligibility {
    pub policy: QcPolicy,
    pub state: EligibilityState,
    pub reason: Option<String>,
    pub evidence: Option<Evidence>,
}
impl Eligibility {
    pub fn fingerprint(&self) -> Result<String> {
        let inputs = self
            .evidence
            .iter()
            .filter(|item| item.report.checksum_type == "sha256")
            .map(|item| InputIdentity {
                role: "QC_REPORT".into(),
                name: "preprocess_report".into(),
                sha256: item.report.checksum.clone(),
                bytes: item.report.bytes,
            })
            .collect();
        StageRequest { schema_version:1, stage:"eligibility".into(), algorithm_version:1,
            engines:BTreeMap::new(), parameters:serde_json::json!({"policy":self.policy, "state":self.state, "reason":self.reason}),
            inputs, dependencies:BTreeMap::new() }.digest()
    }
    pub fn eligible(&self) -> bool {
        self.state.eligible()
    }
    pub fn needs_preprocessing(&self) -> bool {
        matches!(
            self.state,
            EligibilityState::PreprocessRequired | EligibilityState::InvalidEvidence
        )
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StoredEligibility {
    pub decision_id: String,
    pub project_id: String,
    pub run_accession: String,
    pub fingerprint: String,
    pub descriptor: Eligibility,
    pub created_at: String,
}
pub struct EligibilityChange {
    pub previous: Option<StoredEligibility>,
    pub current: StoredEligibility,
    pub changed: bool,
}
pub struct Assessment {
    pub descriptor: Eligibility,
    pub qc: Option<PreprocessQcRecord>,
}

impl Assessment {
    pub fn from_qc(config: &ProjectConfig, evidence: Evidence, mut qc: PreprocessQcRecord) -> Self {
        let state = if qc.mode == "bypass" {
            qc.gate_passed = !config.quant.preprocess_strict;
            qc.gate_reason = config
                .quant
                .preprocess_strict
                .then(|| "preprocessing bypass is disabled by current policy".into());
            if qc.gate_passed {
                EligibilityState::Eligible
            } else {
                EligibilityState::PreprocessRequired
            }
        } else if qc.gate_passed {
            EligibilityState::Eligible
        } else {
            EligibilityState::Rejected
        };
        Self {
            descriptor: Eligibility {
                policy: QcPolicy::from_config(config),
                state,
                reason: qc.gate_reason.clone(),
                evidence: Some(evidence),
            },
            qc: Some(qc),
        }
    }

    fn blocked(
        policy: QcPolicy,
        evidence: Option<Evidence>,
        state: EligibilityState,
        reason: impl Into<String>,
    ) -> Self {
        Self {
            descriptor: Eligibility {
                policy,
                evidence,
                state,
                reason: Some(reason.into()),
            },
            qc: None,
        }
    }
}

/// Verify the persisted report and producer, including intentional FASTQ pruning.
/// Missing evidence is not inferred from an earlier progress state or a QC row.
pub fn assess(db: &Database, config: &ProjectConfig, run: &str) -> Result<Assessment> {
    config.quant.validate()?;
    let policy = QcPolicy::from_config(config);
    if !policy.enabled {
        return Ok(Assessment {
            descriptor: Eligibility {
                policy,
                state: EligibilityState::Eligible,
                reason: None,
                evidence: None,
            },
            qc: None,
        });
    }
    let artifacts = db.list_artifacts_for_run(Some(run))?;
    let reports = artifacts
        .iter()
        .filter(|item| item.kind == ArtifactKind::PreprocessReport)
        .collect::<Vec<_>>();
    if reports.is_empty() {
        return Ok(Assessment::blocked(
            policy,
            None,
            EligibilityState::PreprocessRequired,
            "no verified preprocessing report is available",
        ));
    }
    if reports.len() != 1 {
        return Ok(Assessment::blocked(
            policy,
            None,
            EligibilityState::InvalidEvidence,
            "more than one active preprocessing report",
        ));
    }
    let report = &reports[0];
    let producer = db.generation_head(OutputGroup::Preprocess, run)?;
    let evidence = Evidence {
        report: ManifestArtifact {
            kind: report.kind.as_str().into(),
            path: report.path.clone(),
            checksum_type: report.checksum_type.clone(),
            checksum: report.checksum.clone(),
            bytes: report.bytes,
        },
        producer_generation: producer.as_ref().map(|head| head.generation_id.clone()),
    };
    if let Some(producer) = producer
        && (!Path::new(&report.path).starts_with(&producer.payload_path)
            || !db.verify_generation(OutputGroup::Preprocess, run)?)
    {
        return Ok(Assessment::blocked(
            policy,
            Some(evidence),
            EligibilityState::InvalidEvidence,
            "preprocessing generation failed integrity verification",
        ));
    }
    if !evidence.report.matches_file().unwrap_or(false) {
        return Ok(Assessment::blocked(
            policy,
            Some(evidence),
            EligibilityState::InvalidEvidence,
            "preprocessing report is missing or changed since publication",
        ));
    }
    if evidence.producer_generation.is_some() {
        let sources = db
            .generation_inputs(OutputGroup::Preprocess, run)?
            .unwrap_or_default();
        let retained = db.verified_retained_inputs(run)?;
        let active = artifacts
            .iter()
            .filter(|item| {
                matches!(
                    item.kind,
                    ArtifactKind::FastqR1 | ArtifactKind::FastqR2 | ArtifactKind::FastqSingle
                )
            })
            .collect::<Vec<_>>();
        let source_changed = sources.is_empty()
            || active.iter().any(|item| {
                !sources
                    .iter()
                    .any(|old| old.path == item.path && old.kind == item.kind.as_str())
            })
            || sources.iter().any(|source| {
                let registered = active
                    .iter()
                    .any(|item| item.path == source.path && item.kind.as_str() == source.kind);
                let removed = retained.contains(source);
                (!registered && !removed)
                    || if Path::new(&source.path).exists() {
                        !source.matches_file().unwrap_or(false)
                    } else {
                        !removed
                    }
            });
        if source_changed {
            return Ok(Assessment::blocked(
                policy,
                Some(evidence),
                EligibilityState::PreprocessRequired,
                "preprocessing source set changed or disappeared without retention evidence",
            ));
        }
    }
    let qc = match parse_report(
        db.project_id(),
        run,
        Path::new(&report.path),
        &config.quant.qc_gate,
    ) {
        Ok(qc) => qc,
        Err(_) => {
            return Ok(Assessment::blocked(
                policy,
                Some(evidence),
                EligibilityState::InvalidEvidence,
                "preprocessing report is malformed",
            ));
        }
    };
    if !evidence.report.matches_file()? {
        bail!("preprocessing report changed during eligibility evaluation");
    }
    Ok(Assessment::from_qc(config, evidence, qc))
}

fn read_current(
    connection: &Connection,
    project: &str,
    run: &str,
) -> Result<Option<StoredEligibility>> {
    let row = connection.query_row(
        "SELECT d.decision_id,d.fingerprint,d.descriptor_json,d.created_at,d.project_id,d.run_accession
        FROM eligibility_heads h LEFT JOIN eligibility_decisions d ON d.decision_id=h.decision_id
        WHERE h.project_id=?1 AND h.run_accession=?2",
        params![project,run], |row| Ok((row.get::<_,Option<String>>(0)?,row.get::<_,Option<String>>(1)?,
            row.get::<_,Option<String>>(2)?,row.get::<_,Option<String>>(3)?,row.get::<_,Option<String>>(4)?,row.get::<_,Option<String>>(5)?))
    ).optional()?;
    let row = row
        .map(
            |(id, fingerprint, json, created_at, owner, subject)| -> Result<_> {
                if owner.as_deref() != Some(project) || subject.as_deref() != Some(run) {
                    bail!("eligibility head points to a missing or unrelated decision");
                }
                Ok((
                    id.ok_or_else(|| anyhow::anyhow!("missing eligibility decision ID"))?,
                    fingerprint
                        .ok_or_else(|| anyhow::anyhow!("missing eligibility fingerprint"))?,
                    json.ok_or_else(|| anyhow::anyhow!("missing eligibility descriptor"))?,
                    created_at
                        .ok_or_else(|| anyhow::anyhow!("missing eligibility creation time"))?,
                ))
            },
        )
        .transpose()?;
    row.map(|(decision_id, fingerprint, json, created_at)| {
        let descriptor: Eligibility = serde_json::from_str(&json)?;
        if descriptor.fingerprint()? != fingerprint {
            bail!("stored eligibility fingerprint is corrupt");
        }
        Ok(StoredEligibility {
            decision_id,
            fingerprint,
            descriptor,
            created_at,
            project_id: project.into(),
            run_accession: run.into(),
        })
    })
    .transpose()
}

impl Database {
    pub fn eligibility(&self, run: &str) -> Result<Option<StoredEligibility>> {
        self.read(|connection| read_current(connection, self.project_id(), run))
    }
    pub fn eligibility_history(&self, run: &str) -> Result<Vec<StoredEligibility>> {
        self.read(|connection| {
            let mut query=connection.prepare("SELECT decision_id,fingerprint,descriptor_json,created_at FROM eligibility_decisions
                WHERE project_id=?1 AND run_accession=?2 ORDER BY created_at,decision_id")?;
            query.query_map(params![self.project_id(),run],|row|Ok((row.get::<_,String>(0)?,row.get::<_,String>(1)?,row.get::<_,String>(2)?,row.get::<_,String>(3)?)))?
                .map(|row|{let (decision_id,fingerprint,json,created_at)=row?; let descriptor:Eligibility=serde_json::from_str(&json)?;
                    if descriptor.fingerprint()? != fingerprint {bail!("stored eligibility fingerprint is corrupt");}
                    Ok(StoredEligibility {decision_id,fingerprint,descriptor,created_at,project_id:self.project_id().into(),run_accession:run.into()})}).collect()
        })
    }
    pub fn record_eligibility(
        &self,
        run: &str,
        assessment: &Assessment,
    ) -> Result<EligibilityChange> {
        self.write(|connection| {
            if let Some(id)=assessment.descriptor.evidence.as_ref().and_then(|evidence|evidence.producer_generation.as_deref()) {
                let retired:bool=connection.query_row("SELECT EXISTS(SELECT 1 FROM generation_retirements WHERE generation_id=?1)",[id],|row|row.get(0))?;
                if retired {bail!("restore retired QC evidence before selecting it");}
            }
            let previous=read_current(connection,self.project_id(),run)?;
            if let Some(previous)=&previous && previous.descriptor == assessment.descriptor {
                return Ok(EligibilityChange {previous:Some(previous.clone()),current:previous.clone(),changed:false});
            }
            let current=StoredEligibility {decision_id:uuid::Uuid::new_v4().to_string(),project_id:self.project_id().into(),run_accession:run.into(),
                fingerprint:assessment.descriptor.fingerprint()?,descriptor:assessment.descriptor.clone(),created_at:now_rfc3339()};
            connection.execute("INSERT INTO eligibility_decisions(decision_id,project_id,run_accession,fingerprint,descriptor_json,created_at)
                VALUES(?1,?2,?3,?4,?5,?6)",params![current.decision_id,self.project_id(),run,current.fingerprint,serde_json::to_string(&current.descriptor)?,current.created_at])?;
            connection.execute("INSERT INTO eligibility_heads(project_id,run_accession,decision_id) VALUES(?1,?2,?3)
                ON CONFLICT(project_id,run_accession) DO UPDATE SET decision_id=excluded.decision_id",params![self.project_id(),run,current.decision_id])?;
            // Eligibility changes never demote or restore scientific progress.
            connection.execute("UPDATE project_stages SET state='PENDING', last_error='QC eligibility changed',updated_at=?1
                WHERE project_id=?2 AND stage IN ('normalize','deseq2','corr')",params![now_rfc3339(),self.project_id()])?;
            connection.execute("INSERT INTO events(project_id,ts,stage,run_accession,message,context_json) VALUES(?1,?2,'qc',?3,'QC eligibility updated',?4)",
                params![self.project_id(),now_rfc3339(),run,serde_json::json!({"decision_id":current.decision_id,"fingerprint":current.fingerprint,
                    "eligible":current.descriptor.eligible(),"decision":current.descriptor.state,"reason":current.descriptor.reason,
                    "previous_decision_id":previous.as_ref().map(|old|&old.decision_id)}).to_string()])?;
            Ok(EligibilityChange { previous,current,changed:true })
        })
    }
}

pub fn quantified(state: RunState) -> bool {
    matches!(
        state,
        RunState::QuantDone
            | RunState::DeRunning
            | RunState::DeDone
            | RunState::DeFailed
            | RunState::CorrRunning
            | RunState::CorrDone
            | RunState::CorrFailed
            | RunState::Cleaned
    )
}
pub fn parse_report(
    project_id: &str,
    run_accession: &str,
    report_path: &Path,
    gate: &crate::PreprocessQcGate,
) -> Result<crate::PreprocessQcRecord> {
    let report_text = fs::read_to_string(report_path).with_context(|| {
        format!(
            "failed to read preprocess QC report {}",
            report_path.display()
        )
    })?;
    let report: serde_json::Value = serde_json::from_str(&report_text).with_context(|| {
        format!(
            "failed to parse preprocess QC report {}",
            report_path.display()
        )
    })?;

    let passed_reads = json_u64(&report, &["filtering_result", "passed_filter_reads"])
        .or_else(|| json_u64(&report, &["passed_reads"]))
        .unwrap_or(0);
    let failed_reads = json_u64(&report, &["failed_reads"]).unwrap_or_else(|| {
        json_u64(&report, &["filtering_result", "low_quality_reads"]).unwrap_or(0)
            + json_u64(&report, &["filtering_result", "low_complexity_reads"]).unwrap_or(0)
            + json_u64(&report, &["filtering_result", "too_many_N_reads"]).unwrap_or(0)
            + json_u64(&report, &["filtering_result", "too_short_reads"]).unwrap_or(0)
    });
    let total_reads_before = json_u64(&report, &["summary", "before_filtering", "total_reads"])
        .unwrap_or(passed_reads.saturating_add(failed_reads));
    let total_reads_before = total_reads_before.max(passed_reads.saturating_add(failed_reads));
    let pass_rate = ratio(passed_reads, total_reads_before);
    let failed_rate = ratio(failed_reads, total_reads_before);
    let low_quality_reads =
        json_u64(&report, &["filtering_result", "low_quality_reads"]).unwrap_or(0);
    let low_complexity_reads =
        json_u64(&report, &["filtering_result", "low_complexity_reads"]).unwrap_or(0);
    let too_many_n_reads =
        json_u64(&report, &["filtering_result", "too_many_N_reads"]).unwrap_or(0);
    let too_short_reads = json_u64(&report, &["filtering_result", "too_short_reads"]).unwrap_or(0);
    let duplicated_reads = json_u64(&report, &["duplicated_reads"]).unwrap_or(0);
    let duplication_rate = json_f64(&report, &["duplication", "rate"]);
    let adapter_trimmed_reads = json_u64(&report, &["adapter_cutting", "adapter_trimmed_reads"])
        .unwrap_or_else(|| json_u64(&report, &["adapter_trimmed_reads"]).unwrap_or(0));
    let adapter_trimmed_bases = json_u64(&report, &["adapter_cutting", "adapter_trimmed_bases"])
        .unwrap_or_else(|| json_u64(&report, &["adapter_trimmed_bases"]).unwrap_or(0));
    let mode = json_str(&report, &["mode"])
        .map(ToString::to_string)
        .unwrap_or_else(|| "unknown".to_string());
    let threads = json_u64(&report, &["threads"]).unwrap_or(1);

    let mut failures = Vec::new();
    if passed_reads == 0 {
        failures.push("no reads passed preprocessing".into());
    }
    if let Some(min_pass_rate) = gate.min_pass_rate
        && pass_rate < min_pass_rate
    {
        failures.push(format!("pass_rate {:.4} < {:.4}", pass_rate, min_pass_rate));
    }
    if let Some(max_low_quality_rate) = gate.max_low_quality_rate
        && ratio(low_quality_reads, total_reads_before) > max_low_quality_rate
    {
        failures.push(format!(
            "low_quality_rate {:.4} > {:.4}",
            ratio(low_quality_reads, total_reads_before),
            max_low_quality_rate
        ));
    }
    if let Some(max_low_complexity_rate) = gate.max_low_complexity_rate
        && ratio(low_complexity_reads, total_reads_before) > max_low_complexity_rate
    {
        failures.push(format!(
            "low_complexity_rate {:.4} > {:.4}",
            ratio(low_complexity_reads, total_reads_before),
            max_low_complexity_rate
        ));
    }
    if let Some(max_too_many_n_rate) = gate.max_too_many_n_rate
        && ratio(too_many_n_reads, total_reads_before) > max_too_many_n_rate
    {
        failures.push(format!(
            "too_many_n_rate {:.4} > {:.4}",
            ratio(too_many_n_reads, total_reads_before),
            max_too_many_n_rate
        ));
    }
    if let Some(max_too_short_rate) = gate.max_too_short_rate
        && ratio(too_short_reads, total_reads_before) > max_too_short_rate
    {
        failures.push(format!(
            "too_short_rate {:.4} > {:.4}",
            ratio(too_short_reads, total_reads_before),
            max_too_short_rate
        ));
    }
    if let Some(max_duplication_rate) = gate.max_duplication_rate
        && duplication_rate.unwrap_or(0.0) > max_duplication_rate
    {
        failures.push(format!(
            "duplication_rate {:.4} > {:.4}",
            duplication_rate.unwrap_or(0.0),
            max_duplication_rate
        ));
    }

    Ok(crate::PreprocessQcRecord {
        project_id: project_id.to_string(),
        run_accession: run_accession.to_string(),
        report_path: report_path.display().to_string(),
        mode,
        threads,
        total_reads_before,
        passed_reads,
        failed_reads,
        pass_rate,
        failed_rate,
        low_quality_reads,
        low_complexity_reads,
        too_many_n_reads,
        too_short_reads,
        duplicated_reads,
        duplication_rate,
        adapter_trimmed_reads,
        adapter_trimmed_bases,
        gate_passed: failures.is_empty(),
        gate_reason: if failures.is_empty() {
            None
        } else {
            Some(failures.join("; "))
        },
        updated_at: now_rfc3339(),
    })
}

fn json_u64(value: &serde_json::Value, path: &[&str]) -> Option<u64> {
    let mut current = value;
    for key in path {
        current = current.get(*key)?;
    }
    current.as_u64()
}

fn json_f64(value: &serde_json::Value, path: &[&str]) -> Option<f64> {
    let mut current = value;
    for key in path {
        current = current.get(*key)?;
    }
    current.as_f64()
}

fn json_str<'a>(value: &'a serde_json::Value, path: &[&str]) -> Option<&'a str> {
    let mut current = value;
    for key in path {
        current = current.get(*key)?;
    }
    current.as_str()
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ArtifactRecord, LibraryLayout, RunRecord,
        util::{compute_sha256, file_size},
    };
    fn setup() -> Result<(tempfile::TempDir, Database, ProjectConfig)> {
        let root = tempfile::tempdir()?;
        let mut config = ProjectConfig::default();
        config.quant.preprocess = true;
        let db = Database::create(root.path(), &config)?;
        db.upsert_run(&RunRecord {
            project_id: db.project_id().into(),
            study_accession: "study".into(),
            source_study_accession: None,
            run_accession: "RUN1".into(),
            sample_accession: None,
            experiment_accession: None,
            library_layout: LibraryLayout::Single,
            instrument_platform: None,
            instrument_model: None,
            remote_files: vec![],
            metadata: BTreeMap::new(),
            state: RunState::QuantDone,
            last_error: None,
            updated_at: now_rfc3339(),
        })?;
        let report = root.path().join("qc.json");
        fs::write(
            &report,
            r#"{"mode":"single","passed_reads":8,"failed_reads":2}"#,
        )?;
        db.record_artifact(&ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: Some("RUN1".into()),
            kind: ArtifactKind::PreprocessReport,
            path: report.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".into(),
            checksum: compute_sha256(&report)?,
            bytes: file_size(&report)?,
            created_at: now_rfc3339(),
        })?;
        Ok((root, db, config))
    }
    #[test]
    fn decisions_are_independent_idempotent_and_rolled_back_with_their_head_and_event() -> Result<()>
    {
        let (_root, db, mut config) = setup()?;
        let first = db.record_eligibility("RUN1", &assess(&db, &config, "RUN1")?)?;
        assert!(first.current.descriptor.eligible());
        assert!(
            !db.record_eligibility("RUN1", &assess(&db, &config, "RUN1")?)?
                .changed
        );
        config.quant.qc_gate.min_pass_rate = Some(0.95);
        let rejected = assess(&db, &config, "RUN1")?;
        let events = db.list_events(100)?.len();
        let failed: Result<()> = db.atomic(|db| {
            db.record_eligibility("RUN1", &rejected)?;
            bail!("injected rollback")
        });
        assert!(failed.is_err());
        assert_eq!(db.eligibility_history("RUN1")?.len(), 1);
        assert_eq!(db.list_events(100)?.len(), events);
        assert_eq!(
            db.eligibility("RUN1")?.unwrap().decision_id,
            first.current.decision_id
        );
        assert!(db.record_eligibility("RUN1", &rejected)?.changed);
        assert_eq!(db.get_run("RUN1")?.unwrap().state, RunState::QuantDone);
        assert_eq!(db.eligibility_history("RUN1")?.len(), 2);
        let mut relocated = rejected.descriptor.clone();
        relocated.evidence.as_mut().unwrap().report.path = "/elsewhere/qc.json".into();
        relocated.evidence.as_mut().unwrap().producer_generation = Some("other-publication".into());
        assert_eq!(relocated.fingerprint()?, rejected.descriptor.fingerprint()?);
        Ok(())
    }
    #[test]
    fn a_head_cannot_silently_refer_to_another_runs_decision() -> Result<()> {
        let (_root, db, config) = setup()?;
        let first = db.record_eligibility("RUN1", &assess(&db, &config, "RUN1")?)?;
        db.write(|connection| {
            connection.execute(
                "UPDATE eligibility_decisions SET run_accession='RUN2' WHERE decision_id=?1",
                [&first.current.decision_id],
            )?;
            Ok(())
        })?;
        assert!(db.eligibility("RUN1").is_err());
        Ok(())
    }
    #[test]
    fn strict_bypass_requires_processing_but_zero_passing_reads_are_always_rejected() -> Result<()>
    {
        let (_root, db, mut config) = setup()?;
        let assessment = assess(&db, &config, "RUN1")?;
        let evidence = assessment.descriptor.evidence.unwrap();
        let mut qc = assessment.qc.unwrap();
        qc.mode = "bypass".into();
        config.quant.preprocess_strict = true;
        assert_eq!(
            Assessment::from_qc(&config, evidence.clone(), qc.clone())
                .descriptor
                .state,
            EligibilityState::PreprocessRequired
        );
        config.quant.preprocess_strict = false;
        assert!(
            Assessment::from_qc(&config, evidence.clone(), qc)
                .descriptor
                .eligible()
        );
        let report = Path::new(&evidence.report.path);
        fs::write(
            report,
            r#"{"mode":"single","passed_reads":0,"failed_reads":1}"#,
        )?;
        let gate = crate::PreprocessQcGate {
            min_pass_rate: None,
            max_low_complexity_rate: None,
            max_low_quality_rate: None,
            max_too_many_n_rate: None,
            max_too_short_rate: None,
            max_duplication_rate: None,
        };
        let qc = parse_report(db.project_id(), "RUN1", report, &gate)?;
        assert!(!qc.gate_passed);
        Ok(())
    }
}
