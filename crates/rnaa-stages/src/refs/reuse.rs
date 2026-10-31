//! Current and historical reference validation/selection share one bundle contract.
use super::*;
use std::collections::BTreeSet;

pub(super) fn cached_reference(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    spec: &ReferenceSpec,
) -> Result<Option<ReferenceBundle>> {
    let Some(head) = db.generation_head(OutputGroup::Refs, &spec.reference_id)? else {
        return Ok(None);
    };
    validate_bundle(paths, db, config, spec, &head)
}

fn validate_bundle(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    spec: &ReferenceSpec,
    head: &rnaa_core::generations::GenerationHead,
) -> Result<Option<ReferenceBundle>> {
    if !db.verify_generation_head(head)? {
        return Ok(None);
    }
    let pinned = head.paths(paths)?;
    let result = bundle(&pinned, spec);
    let Ok(manifest) = serde_json::from_slice::<StageManifest>(&fs::read(&result.manifest_path)?)
    else {
        return Ok(None);
    };
    let inputs = source_snapshots(&result)?;
    if [&spec.transcript_source, &spec.gtf_source]
        .iter()
        .zip(&inputs)
        .any(|(source, input)| {
            source
                .sha256
                .as_ref()
                .is_some_and(|expected| *expected != input.checksum)
        })
    {
        return Ok(None);
    }
    let request = reference_request(spec, config, &inputs);
    let digest = request.digest()?;
    let expected = bundle_outputs(&result)
        .into_iter()
        .filter(|(kind, _)| *kind != ArtifactKind::Manifest)
        .map(|(_, path)| path)
        .collect::<BTreeSet<_>>();
    let actual = manifest
        .output_artifacts
        .iter()
        .map(|item| PathBuf::from(&item.path))
        .collect::<BTreeSet<_>>();
    if head.request_sha256 != digest
        || manifest.stage != "refs"
        || manifest.id != spec.reference_id
        || manifest.schema_version != 1
        || manifest.parameters["request_sha256"] != digest
        || manifest.parameters["request"] != serde_json::to_value(&request)?
        || manifest.input_artifacts != inputs
        || actual != expected
        || manifest.output_artifacts.len() != expected.len()
        || !manifest.verify_local_artifacts()?
    {
        return Ok(None);
    }
    verify_originals(spec)?;
    Ok(Some(result))
}

pub(super) fn reuse_reference(
    paths: &ProjectPaths,
    db: &Database,
    parent: &Database,
    config: &ProjectConfig,
    preparation: &Preparation<'_>,
) -> Result<Option<ReferenceBundle>> {
    let spec = preparation.spec;
    let expected = preparation.expected;
    // Backfill only full request objects whose digest is already sealed. Source
    // observations and old payloads are never rewritten to enable reuse.
    for head in db.committed_generation_history(OutputGroup::Refs, &spec.reference_id)? {
        if db.scientific_request(&head.request_sha256)?.is_some() {
            continue;
        }
        let result = bundle(&head.paths(paths)?, spec);
        let Some(manifest) = fs::read(&result.manifest_path)
            .ok()
            .and_then(|bytes| serde_json::from_slice::<StageManifest>(&bytes).ok())
        else {
            continue;
        };
        let value = &manifest.parameters["request"];
        let Ok(request) = serde_json::from_value::<StageRequest<serde_json::Value>>(value.clone())
        else {
            continue;
        };
        if request.stage == "refs" && request.digest()? == head.request_sha256 {
            db.index_generation_request(&head, value)?;
        }
    }
    let policy = reference_request(spec, config, &[]).digest()?;
    for head in
        db.generation_policy_candidates(OutputGroup::Refs, "refs", &spec.reference_id, &policy)?
    {
        parent.check_job()?;
        let Some(result) = validate_bundle(paths, db, config, spec, &head)? else {
            continue;
        };
        let Some(verified) = db.prepare_generation_reuse(&head)? else {
            continue;
        };
        let selection = ReferenceSelection::capture(&result, head.clone(), config)?;
        verify_originals(spec)?;
        preparation.observation.verify_current(db)?;
        db.set_job_request(&head.request_sha256)?;
        db.reattach_generation(verified,|db| {
            preparation.observation.verify_catalog(db)?;
            parent.check_job()?;
            db.select_reference(&selection,expected)?;
            db.append_event("refs",None,"reference bundle selected from verified history",json!({"reference_id":result.id,"generation_id":head.generation_id,"request_sha256":head.request_sha256,
                "source_observations":{"transcriptome":spec.transcript_source,"annotations":spec.gtf_source}}))
        })?;
        return Ok(Some(result));
    }
    Ok(None)
}
