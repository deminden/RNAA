CREATE TABLE metadata_publications (
  publication_id TEXT PRIMARY KEY,
  project_id TEXT NOT NULL REFERENCES projects(project_id),
  descriptor_json TEXT NOT NULL,
  created_at TEXT NOT NULL
);
CREATE TABLE metadata_heads (
  project_id TEXT PRIMARY KEY REFERENCES projects(project_id),
  publication_id TEXT NOT NULL REFERENCES metadata_publications(publication_id),
  installed_id TEXT REFERENCES metadata_publications(publication_id)
);
