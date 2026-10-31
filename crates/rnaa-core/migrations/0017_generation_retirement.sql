CREATE TABLE generation_store_pins (
  pin_id TEXT PRIMARY KEY,
  owner TEXT NOT NULL,
  created_at TEXT NOT NULL
);
CREATE TABLE generation_keeps (
  generation_id TEXT PRIMARY KEY REFERENCES artifact_generations(generation_id),
  label TEXT NOT NULL,
  created_at TEXT NOT NULL
);
CREATE TABLE generation_retirements (
  generation_id TEXT PRIMARY KEY REFERENCES artifact_generations(generation_id),
  project_id TEXT NOT NULL,
  state TEXT NOT NULL CHECK(state IN ('TRASHING','TRASHED','RESTORING','PURGING','PURGED')),
  inventory_sha256 TEXT NOT NULL,
  artifacts_sha256 TEXT NOT NULL,
  payload_inventory_json TEXT NOT NULL,
  seal_sha256 TEXT NOT NULL,
  operation_id TEXT NOT NULL,
  trashed_at_ms INTEGER NOT NULL,
  purge_after_ms INTEGER NOT NULL
);
