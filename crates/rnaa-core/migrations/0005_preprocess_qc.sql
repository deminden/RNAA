CREATE TABLE IF NOT EXISTS preprocess_qc (
  project_id TEXT NOT NULL,
  run_accession TEXT NOT NULL,
  report_path TEXT NOT NULL,
  mode TEXT NOT NULL,
  threads INTEGER NOT NULL,
  total_reads_before INTEGER NOT NULL,
  passed_reads INTEGER NOT NULL,
  failed_reads INTEGER NOT NULL,
  pass_rate REAL NOT NULL,
  failed_rate REAL NOT NULL,
  low_quality_reads INTEGER NOT NULL,
  low_complexity_reads INTEGER NOT NULL,
  too_many_n_reads INTEGER NOT NULL,
  too_short_reads INTEGER NOT NULL,
  duplicated_reads INTEGER NOT NULL,
  duplication_rate REAL,
  adapter_trimmed_reads INTEGER NOT NULL,
  adapter_trimmed_bases INTEGER NOT NULL,
  gate_passed INTEGER NOT NULL,
  gate_reason TEXT,
  updated_at TEXT NOT NULL,
  PRIMARY KEY (project_id, run_accession)
);

CREATE INDEX IF NOT EXISTS idx_preprocess_qc_gate
  ON preprocess_qc(gate_passed);
