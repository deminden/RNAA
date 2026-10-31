//! Normalized preprocessing QC catalog access.
use super::Database;
use crate::model::PreprocessQcRecord;
use anyhow::{Context, Result};
use rusqlite::{OptionalExtension, Row, params};

impl Database {
    pub fn record_preprocess_qc(&self, qc: &PreprocessQcRecord) -> Result<()> {
        self.write(|conn| {
            conn.execute(
                "INSERT INTO preprocess_qc (
                   project_id, run_accession, report_path, mode, threads, total_reads_before,
                   passed_reads, failed_reads, pass_rate, failed_rate, low_quality_reads,
                   low_complexity_reads, too_many_n_reads, too_short_reads, duplicated_reads,
                   duplication_rate, adapter_trimmed_reads, adapter_trimmed_bases, gate_passed,
                   gate_reason, updated_at
                 ) VALUES (
                   ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16,
                   ?17, ?18, ?19, ?20, ?21
                 )
                 ON CONFLICT(project_id, run_accession) DO UPDATE SET
                   report_path = excluded.report_path,
                   mode = excluded.mode,
                   threads = excluded.threads,
                   total_reads_before = excluded.total_reads_before,
                   passed_reads = excluded.passed_reads,
                   failed_reads = excluded.failed_reads,
                   pass_rate = excluded.pass_rate,
                   failed_rate = excluded.failed_rate,
                   low_quality_reads = excluded.low_quality_reads,
                   low_complexity_reads = excluded.low_complexity_reads,
                   too_many_n_reads = excluded.too_many_n_reads,
                   too_short_reads = excluded.too_short_reads,
                   duplicated_reads = excluded.duplicated_reads,
                   duplication_rate = excluded.duplication_rate,
                   adapter_trimmed_reads = excluded.adapter_trimmed_reads,
                   adapter_trimmed_bases = excluded.adapter_trimmed_bases,
                   gate_passed = excluded.gate_passed,
                   gate_reason = excluded.gate_reason,
                   updated_at = excluded.updated_at",
                params![
                    qc.project_id,
                    qc.run_accession,
                    qc.report_path,
                    qc.mode,
                    qc.threads as i64,
                    qc.total_reads_before as i64,
                    qc.passed_reads as i64,
                    qc.failed_reads as i64,
                    qc.pass_rate,
                    qc.failed_rate,
                    qc.low_quality_reads as i64,
                    qc.low_complexity_reads as i64,
                    qc.too_many_n_reads as i64,
                    qc.too_short_reads as i64,
                    qc.duplicated_reads as i64,
                    qc.duplication_rate,
                    qc.adapter_trimmed_reads as i64,
                    qc.adapter_trimmed_bases as i64,
                    if qc.gate_passed { 1_i64 } else { 0_i64 },
                    qc.gate_reason,
                    qc.updated_at,
                ],
            )
            .with_context(|| format!("failed to record preprocess QC for {}", qc.run_accession))?;
            Ok(())
        })
    }

    pub fn get_preprocess_qc(&self, run_accession: &str) -> Result<Option<PreprocessQcRecord>> {
        self.read(|conn| {
            conn.query_row(
                "SELECT project_id, run_accession, report_path, mode, threads, total_reads_before,
                        passed_reads, failed_reads, pass_rate, failed_rate, low_quality_reads,
                        low_complexity_reads, too_many_n_reads, too_short_reads, duplicated_reads,
                        duplication_rate, adapter_trimmed_reads, adapter_trimmed_bases, gate_passed,
                        gate_reason, updated_at
                 FROM preprocess_qc
                 WHERE project_id = ?1 AND run_accession = ?2",
                params![self.project_id, run_accession],
                parse_preprocess_qc,
            )
            .optional()
            .context("failed to load preprocess qc")
        })
    }

    pub fn list_preprocess_qc(&self) -> Result<Vec<PreprocessQcRecord>> {
        self.read(|conn| {
            let mut stmt = conn.prepare(
                "SELECT project_id, run_accession, report_path, mode, threads, total_reads_before,
                        passed_reads, failed_reads, pass_rate, failed_rate, low_quality_reads,
                        low_complexity_reads, too_many_n_reads, too_short_reads, duplicated_reads,
                        duplication_rate, adapter_trimmed_reads, adapter_trimmed_bases, gate_passed,
                        gate_reason, updated_at
                 FROM preprocess_qc
                 WHERE project_id = ?1
                 ORDER BY run_accession",
            )?;
            let rows = stmt.query_map(params![self.project_id], parse_preprocess_qc)?;
            let mut records = Vec::new();
            for row in rows {
                records.push(row?);
            }
            Ok(records)
        })
    }
}

fn parse_preprocess_qc(row: &Row<'_>) -> rusqlite::Result<PreprocessQcRecord> {
    Ok(PreprocessQcRecord {
        project_id: row.get(0)?,
        run_accession: row.get(1)?,
        report_path: row.get(2)?,
        mode: row.get(3)?,
        threads: row.get::<_, i64>(4)? as u64,
        total_reads_before: row.get::<_, i64>(5)? as u64,
        passed_reads: row.get::<_, i64>(6)? as u64,
        failed_reads: row.get::<_, i64>(7)? as u64,
        pass_rate: row.get(8)?,
        failed_rate: row.get(9)?,
        low_quality_reads: row.get::<_, i64>(10)? as u64,
        low_complexity_reads: row.get::<_, i64>(11)? as u64,
        too_many_n_reads: row.get::<_, i64>(12)? as u64,
        too_short_reads: row.get::<_, i64>(13)? as u64,
        duplicated_reads: row.get::<_, i64>(14)? as u64,
        duplication_rate: row.get(15)?,
        adapter_trimmed_reads: row.get::<_, i64>(16)? as u64,
        adapter_trimmed_bases: row.get::<_, i64>(17)? as u64,
        gate_passed: row.get::<_, i64>(18)? != 0,
        gate_reason: row.get(19)?,
        updated_at: row.get(20)?,
    })
}
