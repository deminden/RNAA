use super::super::*;
use crate::{
    ProjectConfig,
    generations::{Generation, GenerationGcOptions, GenerationHead},
    jobs::{JobLease, JobOutcome},
};

pub(super) fn setup() -> Result<(tempfile::TempDir, Database)> {
    let root = tempfile::tempdir()?;
    let db = Database::create(root.path(), &ProjectConfig::default())?;
    Ok((root, db))
}
pub(super) fn publish(
    db: &Database,
    value: &str,
    inputs: &[ManifestArtifact],
) -> Result<GenerationHead> {
    publish_with_provenance(db, value, inputs, &[])
}
pub(super) fn publish_with_provenance(
    db: &Database,
    value: &str,
    inputs: &[ManifestArtifact],
    provenance: &[ManifestArtifact],
) -> Result<GenerationHead> {
    let mut lease = JobLease::try_acquire(db, "deseq2", None, &[])?.unwrap();
    let owned = lease.database();
    owned.set_job_request(&"a".repeat(64))?;
    let mut generation = Generation::begin(&owned, OutputGroup::De, "deseq2")?;
    let data = generation.head().payload_path.join("counts.tsv");
    fs::write(&data, value)?;
    let file = ManifestArtifact::snapshot("COUNTS", &data)?;
    let manifest = generation.head().payload_path.join("rnaa_manifest.json");
    write_json_pretty(
        &manifest,
        &StageManifest {
            schema_version: 1,
            stage: "deseq2".into(),
            id: db.project_id().into(),
            started_at: now_rfc3339(),
            finished_at: now_rfc3339(),
            parameters: serde_json::json!({"input_lineage":provenance}),
            tool_versions: BTreeMap::new(),
            input_artifacts: inputs.to_vec(),
            output_artifacts: vec![file],
            notes: Vec::new(),
        },
    )?;
    let records = [
        (ArtifactKind::Counts, data.clone()),
        (ArtifactKind::Manifest, manifest.clone()),
    ]
    .into_iter()
    .map(|(kind, path)| -> Result<ArtifactRecord> {
        let file = ManifestArtifact::snapshot(kind.as_str(), &path)?;
        Ok(ArtifactRecord {
            project_id: db.project_id().into(),
            run_accession: None,
            kind,
            path: file.path,
            blob_id: None,
            shared_path: None,
            checksum_type: file.checksum_type,
            checksum: file.checksum,
            bytes: file.bytes,
            created_at: now_rfc3339(),
        })
    })
    .collect::<Result<Vec<_>>>()?;
    generation.seal_with_inputs(&[data, manifest], &"a".repeat(64), inputs)?;
    let head = generation.head().clone();
    generation.commit(&records, |_| Ok(()))?;
    lease.finish(JobOutcome::Succeeded, None)?;
    Ok(head)
}
pub(super) fn input(head: &GenerationHead) -> Result<ManifestArtifact> {
    ManifestArtifact::snapshot("COUNTS", &head.payload_path.join("counts.tsv"))
}
pub(super) fn downgrade(
    db: &Database,
    head: &GenerationHead,
    omit_inputs: bool,
) -> Result<Vec<u8>> {
    let mut seal = db.committed_seal(&head.generation_id)?.unwrap();
    seal.lineage_version = 0;
    seal.dependencies.clear();
    if omit_inputs {
        seal.inputs.clear();
    }
    db.write(|conn| {
        conn.execute(
            "DELETE FROM generation_dependencies WHERE consumer_generation_id=?1",
            [&head.generation_id],
        )?;
        conn.execute(
            "UPDATE artifact_generations SET inventory_json=?1 WHERE generation_id=?2",
            params![serde_json::to_string(&seal)?, head.generation_id],
        )?;
        Ok(())
    })?;
    let path = head.payload_path.parent().unwrap().join("seal.json");
    write_json_pretty(&path, &seal)?;
    Ok(fs::read(path)?)
}
pub(super) fn gc(apply: bool) -> GenerationGcOptions {
    GenerationGcOptions {
        apply,
        purge_expired: false,
        trash_days: 7,
    }
}
