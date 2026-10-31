CREATE TABLE retained_inputs (
  project_id TEXT NOT NULL,
  run_accession TEXT NOT NULL,
  generation_sha256 TEXT NOT NULL,
  path TEXT NOT NULL,
  artifact_json TEXT NOT NULL,
  recorded_at TEXT NOT NULL,
  PRIMARY KEY (project_id, run_accession, generation_sha256, path)
);
