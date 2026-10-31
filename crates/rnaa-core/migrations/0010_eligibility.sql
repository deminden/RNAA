CREATE TABLE eligibility_decisions (
    decision_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    run_accession TEXT NOT NULL,
    fingerprint TEXT NOT NULL,
    descriptor_json TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE INDEX eligibility_history ON eligibility_decisions(project_id, run_accession, created_at);
CREATE TABLE eligibility_heads (
    project_id TEXT NOT NULL,
    run_accession TEXT NOT NULL,
    decision_id TEXT NOT NULL REFERENCES eligibility_decisions(decision_id),
    PRIMARY KEY(project_id, run_accession)
);
