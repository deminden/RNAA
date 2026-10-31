CREATE TABLE reference_selection (
    project_id TEXT PRIMARY KEY,
    generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
    fingerprint TEXT NOT NULL,
    descriptor_json TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
