CREATE TABLE job_resources (
  project_id TEXT NOT NULL,
  resource TEXT NOT NULL,
  fencing_token INTEGER NOT NULL DEFAULT 0,
  PRIMARY KEY (project_id, resource)
);

CREATE TABLE job_attempts (
  attempt_id TEXT PRIMARY KEY,
  project_id TEXT NOT NULL,
  resource TEXT NOT NULL,
  stage TEXT NOT NULL,
  run_accession TEXT,
  request_sha256 TEXT,
  owner TEXT NOT NULL,
  fencing_token INTEGER NOT NULL,
  state TEXT NOT NULL CHECK (state IN ('RUNNING', 'SUCCEEDED', 'FAILED', 'CANCELLED', 'INTERRUPTED')),
  cancel_requested INTEGER NOT NULL DEFAULT 0 CHECK (cancel_requested IN (0, 1)),
  started_at TEXT NOT NULL,
  heartbeat_ms INTEGER NOT NULL,
  lease_expires_ms INTEGER NOT NULL,
  finished_at TEXT,
  last_error TEXT,
  UNIQUE (project_id, resource, fencing_token)
);
CREATE UNIQUE INDEX one_live_attempt_per_resource
  ON job_attempts(project_id, resource) WHERE state = 'RUNNING';
CREATE INDEX jobs_by_state ON job_attempts(project_id, state, lease_expires_ms);
