CREATE TABLE IF NOT EXISTS project_stages (
  project_id TEXT NOT NULL,
  stage TEXT NOT NULL,
  state TEXT NOT NULL,
  last_error TEXT,
  updated_at TEXT NOT NULL,
  PRIMARY KEY (project_id, stage)
);

CREATE INDEX IF NOT EXISTS idx_project_stages_state ON project_stages(state);
