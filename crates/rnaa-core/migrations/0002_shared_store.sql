ALTER TABLE artifacts ADD COLUMN blob_id TEXT;
ALTER TABLE artifacts ADD COLUMN shared_path TEXT;

CREATE INDEX IF NOT EXISTS idx_artifacts_blob_id ON artifacts(blob_id);

CREATE TABLE IF NOT EXISTS shared_blobs (
  blob_id TEXT PRIMARY KEY,
  storage_path TEXT NOT NULL,
  checksum_type TEXT NOT NULL,
  checksum TEXT NOT NULL,
  bytes INTEGER NOT NULL,
  created_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS shared_blob_refs (
  blob_id TEXT NOT NULL,
  project_id TEXT NOT NULL,
  artifact_path TEXT NOT NULL,
  referenced_at TEXT NOT NULL,
  PRIMARY KEY (blob_id, project_id, artifact_path),
  FOREIGN KEY (blob_id) REFERENCES shared_blobs(blob_id)
);

CREATE INDEX IF NOT EXISTS idx_shared_blob_refs_project ON shared_blob_refs(project_id);
