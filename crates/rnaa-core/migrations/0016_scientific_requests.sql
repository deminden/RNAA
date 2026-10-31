CREATE TABLE scientific_requests (
  project_id TEXT NOT NULL,
  request_sha256 TEXT NOT NULL,
  stage TEXT NOT NULL,
  policy_sha256 TEXT NOT NULL,
  request_json TEXT NOT NULL,
  PRIMARY KEY (project_id, request_sha256)
);
CREATE INDEX scientific_request_policy ON scientific_requests(project_id,stage,policy_sha256);

CREATE TABLE generation_input_bindings (
  binding_id TEXT PRIMARY KEY,
  project_id TEXT NOT NULL,
  generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
  attempt_id TEXT NOT NULL REFERENCES job_attempts(attempt_id),
  binding_sha256 TEXT NOT NULL,
  binding_json TEXT NOT NULL
);
CREATE INDEX generation_binding_consumer ON generation_input_bindings(project_id,generation_id);
CREATE TABLE generation_binding_heads (
  project_id TEXT NOT NULL,
  output_group TEXT NOT NULL,
  subject TEXT NOT NULL,
  binding_id TEXT NOT NULL REFERENCES generation_input_bindings(binding_id),
  PRIMARY KEY(project_id,output_group,subject)
);
CREATE TABLE generation_binding_pruning (
  producer_generation_id TEXT NOT NULL REFERENCES artifact_generations(generation_id),
  path TEXT NOT NULL,
  binding_id TEXT NOT NULL REFERENCES generation_input_bindings(binding_id),
  binding_sha256 TEXT NOT NULL,
  quant_manifest_sha256 TEXT NOT NULL,
  artifact_json TEXT NOT NULL,
  PRIMARY KEY(producer_generation_id,path,binding_id)
);
