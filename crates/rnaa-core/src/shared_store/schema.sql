CREATE TABLE IF NOT EXISTS blobs (
    sha256 TEXT PRIMARY KEY,
    bytes INTEGER NOT NULL CHECK(bytes >= 0),
    state TEXT NOT NULL CHECK(state IN ('LIVE','TRASHING','TRASHED','RESTORING','DELETING','DELETED')),
    legacy_hold INTEGER NOT NULL CHECK(legacy_hold IN (0,1)),
    trashed_at_ms INTEGER,
    purge_after_ms INTEGER,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS owners (
    sha256 TEXT NOT NULL REFERENCES blobs(sha256),
    project_db TEXT NOT NULL,
    project_id TEXT NOT NULL,
    artifact_path TEXT NOT NULL,
    state TEXT NOT NULL CHECK(state IN ('PREPARED','COMMITTED')),
    generation_id TEXT,
    created_at TEXT NOT NULL,
    PRIMARY KEY(sha256,project_db,project_id,artifact_path)
);
CREATE TABLE IF NOT EXISTS pins (
    pin_id TEXT PRIMARY KEY,
    sha256 TEXT NOT NULL REFERENCES blobs(sha256),
    owner TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS pins_blob ON pins(sha256);
PRAGMA user_version=1;
