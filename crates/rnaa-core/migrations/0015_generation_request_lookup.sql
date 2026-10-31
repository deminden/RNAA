-- Requests are scoped by project, output contract, stage and subject. Selecting
-- an old result never changes its immutable producer attempt or commit time.
CREATE INDEX generation_request_lookup ON artifact_generations
  (project_id, output_group, stage, subject, request_sha256, state, committed_at);
