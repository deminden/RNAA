-- Resolutions are immutable. Replacing the head never deletes runs or results.
-- One accession can be both a query and a selector; role must be part of its key.
ALTER TABLE inputs RENAME TO inputs_legacy;
CREATE TABLE inputs (
  project_id TEXT NOT NULL,
  input_id TEXT NOT NULL,
  input_type TEXT NOT NULL,
  added_at TEXT NOT NULL,
  PRIMARY KEY (project_id, input_id, input_type)
);
INSERT INTO inputs SELECT project_id, input_id, input_type, added_at FROM inputs_legacy;
DROP TABLE inputs_legacy;

CREATE TABLE cohorts (
  cohort_id TEXT PRIMARY KEY,
  project_id TEXT NOT NULL REFERENCES projects(project_id),
  descriptor_json TEXT NOT NULL,
  created_at TEXT NOT NULL
);

CREATE TABLE cohort_heads (
  project_id TEXT PRIMARY KEY REFERENCES projects(project_id),
  cohort_id TEXT NOT NULL REFERENCES cohorts(cohort_id)
);
