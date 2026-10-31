CREATE TABLE file_trash (
  trash_id TEXT PRIMARY KEY,
  project_id TEXT NOT NULL,
  run_accession TEXT NOT NULL,
  canonical_source TEXT NOT NULL,
  operation_id TEXT NOT NULL,
  state TEXT NOT NULL CHECK(state IN ('MOVING','TRASHED','RESTORING','RESTORED','PURGING','PURGED')),
  descriptor_json TEXT NOT NULL,
  descriptor_sha256 TEXT NOT NULL,
  created_at_ms INTEGER NOT NULL,
  purge_after_ms INTEGER NOT NULL
);
CREATE UNIQUE INDEX file_trash_active_source ON file_trash(project_id,canonical_source)
  WHERE state IN ('MOVING','TRASHED','RESTORING','PURGING');
CREATE INDEX file_trash_operation ON file_trash(project_id,operation_id,state);

CREATE INDEX file_trash_run ON file_trash(project_id,run_accession,created_at_ms);
CREATE INDEX file_trash_source ON file_trash(project_id,canonical_source,created_at_ms);
