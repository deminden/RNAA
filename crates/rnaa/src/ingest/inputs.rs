//! Registered FASTQ lookup and recovery of interrupted ingest work.
use anyhow::Result;
use rnaa_core::Database;
use rnaa_core::model::{ArtifactKind, VerifiedFile};
use rnaa_core::state::RunState;
use serde_json::json;
use std::path::PathBuf;

pub(crate) fn reconcile_verified_download_state(db: &Database) -> Result<()> {
    let runs = db.list_active_runs()?;
    let mut reconciled = 0_u64;
    for queued in runs {
        db.with_idle_run(&queued.run_accession, |db| {
            let Some(run) = db.get_run(&queued.run_accession)? else {
                return Ok(());
            };
            if matches!(
                run.state,
                RunState::Verified
                    | RunState::PreprocessRunning
                    | RunState::PreprocessDone
                    | RunState::PreprocessFailed
                    | RunState::QuantRunning
                    | RunState::QuantDone
                    | RunState::QuantFailed
                    | RunState::DeRunning
                    | RunState::DeDone
                    | RunState::DeFailed
                    | RunState::CorrRunning
                    | RunState::CorrDone
                    | RunState::CorrFailed
                    | RunState::Cleaned
            ) {
                return Ok(());
            }
            let artifacts = db.list_artifacts_for_run(Some(&run.run_accession))?;
            let has_verified_raw = artifacts.iter().any(|artifact| {
                matches!(
                    artifact.kind,
                    ArtifactKind::FastqR1
                        | ArtifactKind::FastqR2
                        | ArtifactKind::FastqSingle
                        | ArtifactKind::Sra
                )
            });
            if has_verified_raw {
                db.set_run_state(&run.run_accession, RunState::Verified, None)?;
                reconciled += 1;
            }
            Ok(())
        })?;
    }
    if reconciled > 0 {
        db.append_event(
            "download",
            None,
            "download state reconciled from existing artifacts",
            json!({ "runs_reconciled": reconciled }),
        )?;
    }
    Ok(())
}

pub(crate) fn recover_interrupted_work(db: &Database) -> Result<()> {
    if !db.recover_input_cleanup()? {
        anyhow::bail!(
            "pending input cleanup is in use; retry storage files gc --apply after readers finish"
        );
    }
    db.recover_jobs()?;
    db.recover_generations()?;
    Ok(())
}

pub(crate) fn has_pending_download_runs(db: &Database) -> Result<bool> {
    Ok(!db
        .list_active_runs_in_states(&[
            RunState::Resolved,
            RunState::Downloading,
            RunState::Downloaded,
        ])?
        .is_empty())
}

pub(crate) fn load_quant_input_fastqs(
    db: &Database,
    run_accession: &str,
    preprocess_enabled: bool,
) -> Result<Vec<VerifiedFile>> {
    if preprocess_enabled {
        let trimmed = load_artifacts_by_kinds(
            db,
            run_accession,
            &[
                ArtifactKind::FastqTrimmedR1,
                ArtifactKind::FastqTrimmedR2,
                ArtifactKind::FastqTrimmedSingle,
            ],
        )?;
        if !trimmed.is_empty() {
            return Ok(trimmed);
        }
    }
    load_raw_fastq_artifacts(db, run_accession)
}

pub(crate) fn load_trimmed_fastq_artifacts(
    db: &Database,
    run_accession: &str,
) -> Result<Vec<VerifiedFile>> {
    load_artifacts_by_kinds(
        db,
        run_accession,
        &[
            ArtifactKind::FastqTrimmedR1,
            ArtifactKind::FastqTrimmedR2,
            ArtifactKind::FastqTrimmedSingle,
        ],
    )
}

pub(crate) fn load_raw_fastq_artifacts(
    db: &Database,
    run_accession: &str,
) -> Result<Vec<VerifiedFile>> {
    load_artifacts_by_kinds(
        db,
        run_accession,
        &[
            ArtifactKind::FastqR1,
            ArtifactKind::FastqR2,
            ArtifactKind::FastqSingle,
        ],
    )
}

fn load_artifacts_by_kinds(
    db: &Database,
    run_accession: &str,
    kinds: &[ArtifactKind],
) -> Result<Vec<VerifiedFile>> {
    let artifacts = db.list_artifacts_for_run(Some(run_accession))?;
    let mut verified = artifacts
        .into_iter()
        .filter(|artifact| kinds.contains(&artifact.kind))
        .map(|artifact| VerifiedFile {
            artifact_kind: artifact.kind,
            path: PathBuf::from(artifact.path),
            checksum_type: artifact.checksum_type,
            checksum: artifact.checksum,
            bytes: artifact.bytes,
            integrity_source: "artifact_db".to_string(),
        })
        .collect::<Vec<_>>();
    verified.sort_by(|left, right| left.path.cmp(&right.path));
    Ok(verified)
}

#[cfg(test)]
mod tests {
    use super::reconcile_verified_download_state;
    use rnaa_core::config::ProjectConfig;
    use rnaa_core::eligibility::parse_report as parse_preprocess_qc_report;
    use rnaa_core::model::{ArtifactKind, ArtifactRecord, LibraryLayout, RunRecord};
    use rnaa_core::state::RunState;
    use rnaa_core::util::{compute_sha256, file_size, now_rfc3339};
    use rnaa_core::{Database, PreprocessQcGate, ProjectPaths};
    use std::collections::BTreeMap;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn parses_normalized_preprocess_qc_and_computes_rates() {
        let temp = TempDir::new().expect("failed to create tempdir");
        let report = temp.path().join("fasterp.json");
        fs::write(
            &report,
            r#"{
              "mode":"single",
              "threads":16,
              "summary":{"before_filtering":{"total_reads":100},"after_filtering":{"total_reads":92}},
              "filtering_result":{
                "passed_filter_reads":92,
                "low_quality_reads":3,
                "low_complexity_reads":2,
                "too_many_N_reads":1,
                "too_short_reads":2,
                "too_long_reads":0
              },
              "duplicated_reads":7,
              "duplication":{"rate":0.07},
              "adapter_cutting":{"adapter_trimmed_reads":40,"adapter_trimmed_bases":120}
            }"#,
        )
        .expect("failed to write report");

        let qc = parse_preprocess_qc_report("P1", "RUN1", &report, &PreprocessQcGate::default())
            .expect("failed to parse qc");
        assert_eq!(qc.threads, 16);
        assert_eq!(qc.total_reads_before, 100);
        assert_eq!(qc.passed_reads, 92);
        assert_eq!(qc.failed_reads, 8);
        assert_eq!(qc.low_quality_reads, 3);
        assert_eq!(qc.too_many_n_reads, 1);
        assert_eq!(qc.adapter_trimmed_reads, 40);
        assert_eq!(qc.adapter_trimmed_bases, 120);
        assert!(qc.gate_passed);
        assert_eq!(qc.pass_rate, 0.92);
        assert_eq!(qc.failed_rate, 0.08);
    }

    #[test]
    fn qc_gate_rejects_run_when_thresholds_fail() {
        let temp = TempDir::new().expect("failed to create tempdir");
        let report = temp.path().join("fasterp.json");
        fs::write(
            &report,
            r#"{
              "mode":"single",
              "threads":8,
              "summary":{"before_filtering":{"total_reads":10},"after_filtering":{"total_reads":7}},
              "filtering_result":{
                "passed_filter_reads":7,
                "low_quality_reads":2,
                "low_complexity_reads":1,
                "too_many_N_reads":0,
                "too_short_reads":0,
                "too_long_reads":0
              },
              "duplicated_reads":0,
              "duplication":{"rate":0.0},
              "adapter_cutting":{"adapter_trimmed_reads":2,"adapter_trimmed_bases":8}
            }"#,
        )
        .expect("failed to write report");

        let qc = parse_preprocess_qc_report(
            "P1",
            "RUN1",
            &report,
            &PreprocessQcGate {
                min_pass_rate: Some(0.8),
                max_low_quality_rate: Some(0.1),
                ..PreprocessQcGate::default()
            },
        )
        .expect("failed to parse qc");
        assert!(!qc.gate_passed);
        let reason = qc.gate_reason.expect("expected gate reason");
        assert!(reason.contains("pass_rate"));
        assert!(reason.contains("low_quality_rate"));
    }

    #[test]
    fn default_qc_gate_uses_conservative_fail_fast_thresholds() {
        let gate = PreprocessQcGate::default();
        assert_eq!(gate.min_pass_rate, Some(0.60));
        assert_eq!(gate.max_low_quality_rate, Some(0.30));
        assert_eq!(gate.max_too_many_n_rate, Some(0.10));
        assert_eq!(gate.max_too_short_rate, Some(0.40));
        assert_eq!(gate.max_low_complexity_rate, None);
        assert_eq!(gate.max_duplication_rate, None);
    }

    #[test]
    fn reconciles_existing_raw_artifacts_back_to_verified() {
        let temp = TempDir::new().expect("failed to create tempdir");
        let root = temp.path();
        let config = ProjectConfig::default();
        let db = Database::create(root, &config).expect("failed to create database");
        let paths = ProjectPaths::new(root);
        let project_id = db.project().expect("failed to load project").project_id;
        let run = RunRecord {
            project_id: project_id.clone(),
            study_accession: "SRPTEST".to_string(),
            source_study_accession: None,
            run_accession: "RUN1".to_string(),
            sample_accession: Some("SAMP1".to_string()),
            experiment_accession: Some("EXP1".to_string()),
            library_layout: LibraryLayout::Single,
            instrument_platform: Some("ILLUMINA".to_string()),
            instrument_model: Some("TEST".to_string()),
            remote_files: Vec::new(),
            metadata: BTreeMap::new(),
            state: RunState::Resolved,
            last_error: None,
            updated_at: now_rfc3339(),
        };
        db.upsert_run(&run).expect("failed to upsert run");

        let raw_dir = paths.raw_run_dir("RUN1");
        fs::create_dir_all(&raw_dir).expect("failed to create raw dir");
        let fastq = raw_dir.join("RUN1.fastq");
        fs::write(&fastq, b"@r1\nACGT\n+\n!!!!\n").expect("failed to write fastq");
        db.record_artifact(&ArtifactRecord {
            project_id,
            run_accession: Some("RUN1".to_string()),
            kind: ArtifactKind::FastqSingle,
            path: fastq.display().to_string(),
            blob_id: None,
            shared_path: None,
            checksum_type: "sha256".to_string(),
            checksum: compute_sha256(&fastq).expect("failed to hash fastq"),
            bytes: file_size(&fastq).expect("failed to stat fastq"),
            created_at: now_rfc3339(),
        })
        .expect("failed to record artifact");

        reconcile_verified_download_state(&db).expect("failed to reconcile download state");

        let run = db
            .list_runs()
            .expect("failed to list runs")
            .into_iter()
            .find(|item| item.run_accession == "RUN1")
            .expect("missing run");
        assert_eq!(run.state, RunState::Verified);
    }
}
