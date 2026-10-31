CREATE TABLE IF NOT EXISTS legacy_imports (
    import_id TEXT PRIMARY KEY,
    descriptor_json TEXT NOT NULL,
    descriptor_sha256 TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS legacy_blob_imports (
    sha256 TEXT PRIMARY KEY REFERENCES blobs(sha256),
    import_id TEXT NOT NULL REFERENCES legacy_imports(import_id),
    released_hold INTEGER NOT NULL CHECK(released_hold IN (0,1))
);
CREATE TABLE IF NOT EXISTS legacy_projects (
    project_db TEXT NOT NULL,
    project_id TEXT NOT NULL,
    first_import_id TEXT NOT NULL REFERENCES legacy_imports(import_id),
    PRIMARY KEY(project_db,project_id)
);
PRAGMA user_version=2;
