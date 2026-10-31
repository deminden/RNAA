CREATE TABLE generation_files (
  generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
  path TEXT NOT NULL,
  sha256 TEXT NOT NULL,
  bytes INTEGER NOT NULL,
  PRIMARY KEY (generation_id, path)
);
CREATE INDEX generation_files_identity ON generation_files(path, sha256, bytes);
INSERT INTO generation_files(generation_id, path, sha256, bytes)
SELECT g.generation_id, json_extract(f.value, '$.path'),
       json_extract(f.value, '$.checksum'), json_extract(f.value, '$.bytes')
FROM artifact_generations g, json_each(g.inventory_json, '$.files') f
WHERE g.state = 'COMMITTED' AND json_valid(g.inventory_json);

CREATE TABLE generation_dependencies (
  consumer_generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
  producer_generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
  producer_seal_sha256 TEXT NOT NULL,
  role TEXT NOT NULL,
  path TEXT NOT NULL,
  artifact_json TEXT NOT NULL,
  PRIMARY KEY (consumer_generation_id, role, path),
  CHECK(consumer_generation_id <> producer_generation_id)
);
CREATE INDEX generation_dependencies_producer ON generation_dependencies(producer_generation_id);
