CREATE TABLE IF NOT EXISTS run_metadata (
  project_id TEXT NOT NULL,
  run_accession TEXT NOT NULL,
  field_name TEXT NOT NULL,
  field_value TEXT NOT NULL,
  updated_at TEXT NOT NULL,
  PRIMARY KEY (project_id, run_accession, field_name)
);

CREATE INDEX IF NOT EXISTS idx_run_metadata_run
  ON run_metadata(run_accession);

CREATE INDEX IF NOT EXISTS idx_run_metadata_field
  ON run_metadata(field_name);
