CREATE TABLE artifact_generations (
  generation_id TEXT PRIMARY KEY,
  project_id TEXT NOT NULL,
  output_group TEXT NOT NULL,
  stage TEXT NOT NULL,
  subject TEXT NOT NULL,
  attempt_id TEXT NOT NULL UNIQUE REFERENCES job_attempts(attempt_id),
  payload_path TEXT NOT NULL UNIQUE,
  request_sha256 TEXT,
  state TEXT NOT NULL CHECK(state IN ('BUILDING','SEALED','COMMITTED','ABANDONED')),
  inventory_json TEXT,
  artifacts_json TEXT,
  created_at TEXT NOT NULL,
  committed_at TEXT
);
CREATE INDEX generation_attempt_state ON artifact_generations(state, attempt_id);
CREATE TABLE generation_heads (
  project_id TEXT NOT NULL,
  output_group TEXT NOT NULL,
  subject TEXT NOT NULL,
  generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
  PRIMARY KEY(project_id, output_group, subject)
);
