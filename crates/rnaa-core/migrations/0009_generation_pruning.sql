CREATE TABLE generation_pruned_files (
    project_id TEXT NOT NULL,
    run_accession TEXT NOT NULL,
    producer_generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
    consumer_generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
    quant_manifest_sha256 TEXT NOT NULL,
    path TEXT NOT NULL,
    artifact_json TEXT NOT NULL,
    recorded_at TEXT NOT NULL,
    PRIMARY KEY (producer_generation_id, path, consumer_generation_id)
);
CREATE INDEX generation_pruning_producer ON generation_pruned_files(project_id, producer_generation_id, path);
