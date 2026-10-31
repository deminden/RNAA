CREATE TABLE generation_lineage_imports (
  generation_id TEXT PRIMARY KEY REFERENCES artifact_generations(generation_id),
  project_id TEXT NOT NULL,
  import_id TEXT NOT NULL UNIQUE,
  descriptor_json TEXT NOT NULL,
  descriptor_sha256 TEXT NOT NULL,
  created_at TEXT NOT NULL
);
CREATE TABLE generation_import_dependencies (
  consumer_generation_id TEXT NOT NULL REFERENCES generation_lineage_imports(generation_id),
  producer_generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
  producer_seal_sha256 TEXT NOT NULL,
  role TEXT NOT NULL,
  path TEXT NOT NULL,
  artifact_json TEXT NOT NULL,
  PRIMARY KEY (consumer_generation_id, role, path),
  CHECK(consumer_generation_id <> producer_generation_id)
);
CREATE INDEX generation_import_dependencies_producer
  ON generation_import_dependencies(producer_generation_id);
