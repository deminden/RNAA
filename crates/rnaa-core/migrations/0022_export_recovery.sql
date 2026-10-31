CREATE TABLE export_publications (
  attempt_id TEXT PRIMARY KEY REFERENCES job_attempts(attempt_id),
  project_id TEXT NOT NULL,
  stage TEXT NOT NULL CHECK (stage IN ('export', 'survey')),
  staged_name TEXT NOT NULL,
  device INTEGER NOT NULL,
  inode INTEGER NOT NULL,
  state TEXT NOT NULL CHECK (state IN ('STAGED', 'PUBLISHED', 'RECOVERED', 'ORPHANED', 'MISSING')),
  UNIQUE (project_id, staged_name)
);
CREATE INDEX export_publications_recovery ON export_publications(project_id, state);
