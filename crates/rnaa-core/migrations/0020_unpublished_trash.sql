-- Unpublished data never becomes a scientific result through cleanup or restore.
CREATE TABLE unpublished_generation_trash (
  generation_id TEXT PRIMARY KEY REFERENCES artifact_generations(generation_id),
  project_id TEXT NOT NULL,
  state TEXT NOT NULL CHECK(state IN ('TRASHING','TRASHED','RESTORING','PURGING','PURGED')),
  descriptor_json TEXT NOT NULL,
  descriptor_sha256 TEXT NOT NULL,
  trashed_at_ms INTEGER NOT NULL,
  purge_after_ms INTEGER NOT NULL
);
