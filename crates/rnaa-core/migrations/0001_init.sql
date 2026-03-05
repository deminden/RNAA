CREATE TABLE IF NOT EXISTS schema_migrations (
  version TEXT PRIMARY KEY,
  applied_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS projects (
  project_id TEXT PRIMARY KEY,
  root_dir TEXT NOT NULL,
  created_at TEXT NOT NULL,
  config_json TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS inputs (
  project_id TEXT NOT NULL,
  input_id TEXT NOT NULL,
  input_type TEXT NOT NULL,
  added_at TEXT NOT NULL,
  PRIMARY KEY (project_id, input_id)
);

CREATE TABLE IF NOT EXISTS metadata_overrides (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  project_id TEXT NOT NULL,
  source_path TEXT NOT NULL,
  stored_path TEXT NOT NULL,
  merge_strategy TEXT NOT NULL,
  added_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS runs (
  project_id TEXT NOT NULL,
  study_accession TEXT NOT NULL,
  source_study_accession TEXT,
  run_accession TEXT PRIMARY KEY,
  sample_accession TEXT,
  experiment_accession TEXT,
  library_layout TEXT NOT NULL,
  instrument_platform TEXT,
  instrument_model TEXT,
  urls_json TEXT NOT NULL,
  md5_json TEXT NOT NULL,
  bytes_json TEXT NOT NULL,
  remote_files_json TEXT NOT NULL,
  metadata_json TEXT NOT NULL,
  state TEXT NOT NULL,
  last_error TEXT,
  updated_at TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_runs_state ON runs(state);
CREATE INDEX IF NOT EXISTS idx_runs_study ON runs(study_accession);

CREATE TABLE IF NOT EXISTS artifacts (
  project_id TEXT NOT NULL,
  run_accession TEXT,
  kind TEXT NOT NULL,
  path TEXT NOT NULL,
  checksum_type TEXT NOT NULL,
  checksum TEXT NOT NULL,
  bytes INTEGER NOT NULL,
  created_at TEXT NOT NULL,
  PRIMARY KEY (project_id, path)
);

CREATE INDEX IF NOT EXISTS idx_artifacts_run ON artifacts(run_accession);
CREATE INDEX IF NOT EXISTS idx_artifacts_kind ON artifacts(kind);

CREATE TABLE IF NOT EXISTS contrasts (
  project_id TEXT NOT NULL,
  name TEXT NOT NULL,
  factor TEXT NOT NULL,
  level_a TEXT NOT NULL,
  level_b TEXT NOT NULL,
  PRIMARY KEY (project_id, name)
);

CREATE TABLE IF NOT EXISTS events (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  project_id TEXT NOT NULL,
  ts TEXT NOT NULL,
  run_accession TEXT,
  stage TEXT NOT NULL,
  message TEXT NOT NULL,
  context_json TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_events_stage ON events(stage);
CREATE INDEX IF NOT EXISTS idx_events_run ON events(run_accession);
