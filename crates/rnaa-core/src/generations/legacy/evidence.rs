//! Only immutable, checksummed declarations supply historical input identities.
use super::*;

fn strong(input: &ManifestArtifact) -> bool {
    input.checksum_type == "sha256"
        && input.checksum.len() == 64
        && input
            .checksum
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn sort_inputs(values: &mut [ManifestArtifact]) {
    values.sort_by(|a, b| (&a.kind, &a.path).cmp(&(&b.kind, &b.path)));
}
pub(super) struct InputSet {
    pub inputs: Vec<ManifestArtifact>,
    pub ancestry: Vec<ManifestArtifact>,
}
fn add_input(
    values: &mut BTreeMap<(String, String), ManifestArtifact>,
    input: ManifestArtifact,
) -> Result<()> {
    if !strong(&input)
        || !Path::new(&input.path).is_absolute()
        || Path::new(&input.path)
            .components()
            .any(|part| matches!(part, std::path::Component::ParentDir))
    {
        bail!("legacy input lacks an absolute path and exact SHA-256 identity");
    }
    if let Some(previous) = values.insert((input.kind.clone(), input.path.clone()), input.clone())
        && previous != input
    {
        bail!("legacy input declarations disagree");
    }
    Ok(())
}
impl Database {
    pub(super) fn legacy_manifest_evidence(&self, seal: &Seal) -> Result<Vec<Evidence>> {
        let (_, registry) = self.legacy_catalog_documents(&seal.head.generation_id)?;
        let records: Vec<ArtifactRecord> = serde_json::from_str(&registry)?;
        let mut evidence = Vec::new();
        for record in records
            .iter()
            .filter(|record| record.kind == ArtifactKind::Manifest)
        {
            let file = seal
                .files
                .iter()
                .find(|file| {
                    file.path == record.path
                        && file.checksum == record.checksum
                        && file.bytes == record.bytes
                })
                .context("legacy manifest registration differs from its sealed payload")?;
            let json = fs::read_to_string(&file.path)?;
            let item = Evidence {
                file: file.clone(),
                json,
            };
            self.verify_legacy_manifest(seal, &item)?;
            self.legacy_descriptor_inputs(seal, &item)?;
            evidence.push(item);
        }
        evidence.sort_by(|a, b| a.file.path.cmp(&b.file.path));
        Ok(evidence)
    }
    pub(super) fn verify_legacy_manifest(
        &self,
        seal: &Seal,
        evidence: &Evidence,
    ) -> Result<Option<StageManifest>> {
        if !seal.files.contains(&evidence.file)
            || evidence.file.checksum_type != "sha256"
            || digest(evidence.json.as_bytes()) != evidence.file.checksum
            || evidence.json.len() as u64 != evidence.file.bytes
        {
            bail!("legacy manifest evidence differs from its sealed bytes");
        }
        let document: serde_json::Value = serde_json::from_str(&evidence.json)?;
        if document.get("input_artifacts").is_none() {
            return Ok(None);
        }
        let manifest: StageManifest = serde_json::from_str(&evidence.json)?;
        if manifest.schema_version != 1
            || manifest.id != seal.head.subject
            || !(manifest.stage == seal.head.stage
                || seal.head.output_group == OutputGroup::De && manifest.stage == "deseq2"
                || seal.head.output_group == OutputGroup::Corr && manifest.stage == "module_gsea")
            || manifest.output_artifacts.is_empty()
        {
            bail!("legacy manifest stage, subject or output declaration is unsupported");
        }
        if let Some(request) = manifest.parameters.get("request_sha256")
            && request.as_str() != Some(&seal.head.request_sha256)
        {
            bail!("legacy manifest request differs from its generation");
        }
        for output in &manifest.output_artifacts {
            if !seal.files.iter().any(|file| {
                file.path == output.path
                    && file.bytes == output.bytes
                    && (output.checksum_type == "none"
                        || output.checksum_type == "sha256" && output.checksum == file.checksum)
            }) {
                bail!("legacy stage output is not bound by the generation seal");
            }
        }
        Ok(Some(manifest))
    }
    fn legacy_descriptor_inputs(
        &self,
        seal: &Seal,
        evidence: &Evidence,
    ) -> Result<Option<Vec<ManifestArtifact>>> {
        let value: serde_json::Value = serde_json::from_str(&evidence.json)?;
        let Some(inputs) = value.get("input_witnesses") else {
            return Ok(None);
        };
        let request: crate::fingerprint::StageRequest<serde_json::Value> =
            serde_json::from_value(value["request"].clone())?;
        if value["schema_version"] != 1
            || value["request_sha256"].as_str() != Some(&seal.head.request_sha256)
            || request.schema_version != 1
            || request.stage != seal.head.stage
            || request.digest()? != seal.head.request_sha256
        {
            bail!("legacy request descriptor is not bound to its generation request");
        }
        let inputs: Vec<ManifestArtifact> = serde_json::from_value(inputs.clone())?;
        if inputs.is_empty() {
            bail!("legacy analysis descriptor has no complete physical input witnesses");
        }
        Ok(Some(inputs))
    }
    pub(super) fn legacy_inputs(
        &self,
        seal: &Seal,
        evidence: &[Evidence],
        observed: Option<&[ManifestArtifact]>,
    ) -> Result<InputSet> {
        let mut manifests = Vec::new();
        let mut descriptor: Option<Vec<ManifestArtifact>> = None;
        for item in evidence {
            if let Some(inputs) = self.legacy_descriptor_inputs(seal, item)? {
                if descriptor.as_ref().is_some_and(|prior| prior != &inputs) {
                    bail!("legacy physical input descriptors disagree");
                }
                descriptor = Some(inputs);
            }
            if let Some(manifest) = self.verify_legacy_manifest(seal, item)? {
                manifests.push(manifest);
            }
        }
        let strong_sources = seal
            .inputs
            .iter()
            .chain(&seal.files)
            .chain(descriptor.iter().flatten())
            .chain(
                manifests
                    .iter()
                    .flat_map(|manifest| &manifest.input_artifacts),
            )
            .filter(|input| strong(input))
            .collect::<Vec<_>>();
        let mut values = BTreeMap::new();
        for input in seal.inputs.iter().chain(descriptor.iter().flatten()) {
            add_input(&mut values, input.clone())?;
        }
        let mut primary: Option<Vec<ManifestArtifact>> = None;
        for manifest in &manifests {
            if seal.head.output_group == OutputGroup::Refs {
                for input in self.legacy_reference_origins(seal, manifest)? {
                    add_input(&mut values, input)?;
                }
            }
            let mut normalized = Vec::new();
            for input in &manifest.input_artifacts {
                if seal.head.output_group == OutputGroup::Download {
                    let remote = ["http://", "https://", "ftp://", "s3://"]
                        .iter()
                        .any(|scheme| input.path.starts_with(scheme));
                    let checksum = strong(input)
                        || input.checksum_type == "none" && input.checksum.is_empty()
                        || input.checksum_type == "md5"
                            && input.checksum.len() == 32
                            && input.checksum.bytes().all(|byte| byte.is_ascii_hexdigit());
                    if !remote || !checksum {
                        bail!(
                            "legacy download origin is not a supported remote-source declaration"
                        );
                    }
                    continue;
                }
                let mut input = input.clone();
                if !strong(&input) {
                    let source = strong_sources.iter().find(|source| source.path == input.path && source.bytes == input.bytes)
                        .context("legacy manifest input has no sealed SHA-256 witness; recomputation is required")?;
                    input.checksum_type = source.checksum_type.clone();
                    input.checksum = source.checksum.clone();
                }
                add_input(&mut values, input.clone())?;
                normalized.push(input);
            }
            if manifest.stage == seal.head.stage
                || seal.head.output_group == OutputGroup::De && manifest.stage == "deseq2"
            {
                if primary.as_ref().is_some_and(|prior| prior != &normalized) {
                    bail!("legacy primary stage manifests disagree on their consumed inputs");
                }
                primary = Some(normalized);
            }
            if let Some(lineage) = manifest.parameters.get("input_lineage") {
                for input in serde_json::from_value::<Vec<ManifestArtifact>>(lineage.clone())? {
                    add_input(&mut values, input)?;
                }
            }
            if seal.head.output_group == OutputGroup::Quant
                && let (Some(path), Some(sha)) = (
                    manifest.parameters.get("index").and_then(|v| v.as_str()),
                    manifest
                        .parameters
                        .get("index_sha256")
                        .and_then(|v| v.as_str()),
                )
            {
                let known = strong_sources
                    .iter()
                    .copied()
                    .chain(observed.unwrap_or_default())
                    .find(|input| {
                        input.kind == "REFERENCE_INDEX"
                            && input.path == path
                            && input.checksum == sha
                    })
                    .cloned();
                let input = if let Some(input) = known {
                    input
                } else {
                    self.legacy_index_input(path, sha)?
                };
                add_input(&mut values, input)?;
            }
        }
        let inputs = if !seal.inputs.is_empty() {
            seal.inputs.clone()
        } else if let Some(inputs) = descriptor.clone() {
            inputs
        } else {
            primary.context(
                "legacy seal has no complete input inventory or supported primary manifest",
            )?
        };
        if seal.head.output_group == OutputGroup::Corr
            && seal.inputs.is_empty()
            && descriptor.is_none()
            && !manifests.iter().any(|manifest| {
                manifest.stage == "corr"
                    && manifest.input_artifacts.iter().any(|input| {
                        strong(input)
                            && !Path::new(&input.path).starts_with(&seal.head.payload_path)
                    })
            })
        {
            bail!(
                "legacy correlation lacks a complete checksummed source inventory for adjustment and modules"
            );
        }
        if inputs.is_empty() && seal.head.output_group != OutputGroup::Download {
            bail!("legacy non-download stage has no verified input declarations");
        }
        let mut ancestry = values.into_values().collect::<Vec<_>>();
        sort_inputs(&mut ancestry);
        if seal.head.output_group == OutputGroup::Quant
            && !ancestry.iter().any(|input| input.kind == "REFERENCE_INDEX")
        {
            bail!("legacy counts lack a checksummed reference index witness");
        }
        Ok(InputSet { inputs, ancestry })
    }
    fn legacy_reference_origins(
        &self,
        seal: &Seal,
        manifest: &StageManifest,
    ) -> Result<Vec<ManifestArtifact>> {
        let mut origins = Vec::new();
        for (field, role) in [
            ("transcript_source", "REFERENCE_CDNA"),
            ("gtf_source", "REFERENCE_GTF"),
        ] {
            let Some(source) = manifest.parameters.get(field) else {
                if seal.inputs.is_empty()
                    && manifest.input_artifacts.iter().any(|input| {
                        input.kind == role
                            && Path::new(&input.path).starts_with(&seal.head.payload_path)
                    })
                {
                    bail!(
                        "legacy materialized reference source is missing its original source descriptor"
                    );
                }
                continue;
            };
            match source.get("kind").and_then(|value| value.as_str()) {
                Some("RemoteUrl") => {
                    let location = source
                        .get("location")
                        .and_then(|value| value.as_str())
                        .context("legacy remote reference origin missing")?;
                    if !["http://", "https://", "ftp://", "s3://"]
                        .iter()
                        .any(|scheme| location.starts_with(scheme))
                    {
                        bail!("unsupported legacy remote reference source");
                    }
                }
                Some("LocalFile") => {
                    let path = source
                        .get("location")
                        .and_then(|value| value.as_str())
                        .context("legacy reference origin path missing")?;
                    let sha = source
                        .get("sha256")
                        .and_then(|value| value.as_str())
                        .context("legacy local reference origin lacks a checksum")?;
                    let materialized = manifest
                        .input_artifacts
                        .iter()
                        .find(|input| input.kind == role && input.checksum == sha && strong(input))
                        .context(
                            "legacy local reference origin differs from its sealed source copy",
                        )?;
                    origins.push(ManifestArtifact {
                        kind: role.into(),
                        path: path.into(),
                        checksum_type: "sha256".into(),
                        checksum: sha.into(),
                        bytes: materialized.bytes,
                    });
                }
                _ => bail!("unsupported legacy reference origin descriptor"),
            }
        }
        Ok(origins)
    }
    fn legacy_index_input(&self, path: &str, sha: &str) -> Result<ManifestArtifact> {
        let lengths=self.read(|conn| {let mut query=conn.prepare("SELECT DISTINCT f.bytes FROM generation_files f JOIN artifact_generations g ON g.generation_id=f.generation_id WHERE g.project_id=?1 AND g.state='COMMITTED' AND f.path=?2 AND f.sha256=?3")?;Ok(query.query_map(params![self.project_id(),path,sha],|row|row.get::<_,i64>(0))?.collect::<std::result::Result<Vec<_>,_>>()?)})?;
        let bytes = if let [length] = lengths.as_slice() {
            u64::try_from(*length)?
        } else if lengths.is_empty() {
            if compute_sha256(Path::new(path))? != sha {
                bail!("legacy external index differs from its sealed checksum");
            }
            fs::metadata(path)?.len()
        } else {
            bail!("legacy reference index has ambiguous sizes");
        };
        Ok(ManifestArtifact {
            kind: "REFERENCE_INDEX".into(),
            path: path.into(),
            checksum_type: "sha256".into(),
            checksum: sha.into(),
            bytes,
        })
    }
}
