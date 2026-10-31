-- Preserve legacy cohort labels before restoring run-only progress semantics.
CREATE TABLE legacy_run_analysis_states (
  run_accession TEXT PRIMARY KEY,
  project_id TEXT NOT NULL,
  analysis_state TEXT NOT NULL,
  analysis_error TEXT,
  source_updated_at TEXT NOT NULL,
  projected_at TEXT NOT NULL
);
INSERT INTO legacy_run_analysis_states
SELECT run_accession, project_id, state, last_error, updated_at,
       strftime('%Y-%m-%dT%H:%M:%fZ', 'now')
FROM runs WHERE state IN ('DE_RUNNING','DE_DONE','DE_FAILED','CORR_RUNNING','CORR_DONE','CORR_FAILED');
UPDATE runs SET state='QUANT_DONE', last_error=NULL,
                updated_at=strftime('%Y-%m-%dT%H:%M:%fZ', 'now')
WHERE run_accession IN (SELECT run_accession FROM legacy_run_analysis_states);
CREATE TRIGGER legacy_analysis_state_immutable_update BEFORE UPDATE ON legacy_run_analysis_states
BEGIN SELECT RAISE(ABORT, 'legacy analysis-state receipts are immutable'); END;
CREATE TRIGGER legacy_analysis_state_immutable_delete BEFORE DELETE ON legacy_run_analysis_states
BEGIN SELECT RAISE(ABORT, 'legacy analysis-state receipts are immutable'); END;
CREATE TRIGGER run_progress_reject_analysis_insert BEFORE INSERT ON runs
WHEN NEW.state IN ('DE_RUNNING','DE_DONE','DE_FAILED','CORR_RUNNING','CORR_DONE','CORR_FAILED')
BEGIN SELECT RAISE(ABORT, 'cohort analysis states cannot be stored as run progress'); END;
CREATE TRIGGER run_progress_reject_analysis_update BEFORE UPDATE OF state ON runs
WHEN NEW.state IN ('DE_RUNNING','DE_DONE','DE_FAILED','CORR_RUNNING','CORR_DONE','CORR_FAILED')
BEGIN SELECT RAISE(ABORT, 'cohort analysis states cannot be stored as run progress'); END;
