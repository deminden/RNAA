//! Explicit external abundance inputs. A declaration records supplied data and
//! reference binding; it never claims that RNAA quantified the underlying reads.
use anyhow::{Context, Result, bail};
use rnaa_core::{
    ManifestArtifact, ReferenceBundle, StageManifest,
    fingerprint::{InputIdentity, StageRequest},
    util::{compute_sha256, now_rfc3339, write_json_pretty},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
};
const FILE: &str = "rnaa_quant_import.json";
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ExternalOptions {
    run: String,
    provenance: String,
    format: String,
}
fn request(
    run: &str,
    provenance: &str,
    inputs: &[ManifestArtifact],
) -> StageRequest<ExternalOptions> {
    StageRequest {
        schema_version: 1,
        stage: "quant_import".into(),
        algorithm_version: 1,
        engines: BTreeMap::new(),
        parameters: ExternalOptions {
            run: run.into(),
            provenance: provenance.into(),
            format: "kallisto_abundance_tsv_v1".into(),
        },
        inputs: inputs.iter().map(InputIdentity::from).collect(),
        dependencies: BTreeMap::new(),
    }
}
/// Atomically declare external kallisto-compatible TSV data and its reference.
/// The caller supplies provenance explicitly. This does not create a native quant
/// generation, assert FASTQ lineage, or turn external data into a quant cache hit.
pub fn declare_external_quant(
    run: &str,
    dir: &Path,
    reference: &ReferenceBundle,
    provenance: &str,
) -> Result<()> {
    if run.trim().is_empty() || provenance.trim().is_empty() {
        bail!("external quantification requires a run ID and explicit provenance");
    }
    if dir.join("rnaa_quant_manifest.json").exists() {
        bail!(
            "native quantification provenance already exists; declare external data in a separate directory"
        );
    }
    let index = kallistors::index::Index::load(&reference.kallisto_index_path)
        .map_err(|error| anyhow::anyhow!("invalid reference index: {error}"))?;
    let expected = index
        .transcripts
        .into_iter()
        .map(|target| (target.name, target.length))
        .collect::<BTreeMap<_, _>>();
    let abundance = dir.join("abundance.tsv");
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&abundance)?;
    let header = reader.headers()?.clone();
    if header.len() != 5 || header.iter().collect::<BTreeSet<_>>().len() != 5 {
        bail!("external abundance requires five unique kallisto column names");
    }
    let column = |name: &str| {
        header
            .iter()
            .position(|value| value == name)
            .with_context(|| format!("external abundance missing {name}"))
    };
    let target = column("target_id")?;
    let length = column("length")?;
    let effective = column("eff_length")?;
    let counts = column("est_counts")?;
    let tpm = column("tpm")?;
    let mut seen = BTreeSet::new();
    for row in reader.records() {
        let row = row?;
        let id = row.get(target).unwrap_or_default();
        let length = row.get(length).unwrap_or_default().parse::<u32>()?;
        if !seen.insert(id.to_string()) || expected.get(id) != Some(&length) {
            bail!("external abundance target or length does not match reference: {id}");
        }
        for col in [effective, counts, tpm] {
            let value = row.get(col).unwrap_or_default().parse::<f64>()?;
            if !value.is_finite() || value < 0.0 {
                bail!("external abundance contains invalid quantitative values");
            }
        }
    }
    if seen.len() != expected.len() {
        bail!("external abundance transcript set does not match reference");
    }
    let inputs = vec![
        ManifestArtifact::snapshot("REFERENCE_INDEX", &reference.kallisto_index_path)?,
        ManifestArtifact::snapshot("EXTERNAL_ABUNDANCE", &abundance)?,
    ];
    let request = request(run, provenance, &inputs);
    let manifest=StageManifest {schema_version:1,stage:"quant_import".into(),id:run.into(),started_at:now_rfc3339(),finished_at:now_rfc3339(),
        parameters:serde_json::json!({"request_sha256":request.digest()?,"request":request,"reference_index_sha256":compute_sha256(&reference.kallisto_index_path)?}),
        tool_versions:BTreeMap::new(),input_artifacts:inputs,output_artifacts:vec![ManifestArtifact::snapshot("EXTERNAL_ABUNDANCE",&abundance)?],
        notes:vec!["Externally supplied transcript abundances; FASTQ/quantifier provenance is supplied by the caller, not reconstructed by RNAA.".into()]};
    write_json_pretty(&dir.join(FILE), &manifest)
}
pub fn external_quant_matches(run: &str, dir: &Path, reference: &ReferenceBundle) -> Result<bool> {
    if dir.join("rnaa_quant_manifest.json").exists() {
        return Ok(false);
    }
    let Some(manifest) = fs::read(dir.join(FILE))
        .ok()
        .and_then(|bytes| serde_json::from_slice::<StageManifest>(&bytes).ok())
    else {
        return Ok(false);
    };
    let Ok(saved) = serde_json::from_value::<StageRequest<ExternalOptions>>(
        manifest.parameters["request"].clone(),
    ) else {
        return Ok(false);
    };
    if manifest.schema_version != 1
        || manifest.stage != "quant_import"
        || manifest.id != run
        || saved.parameters.provenance.trim().is_empty()
        || manifest.input_artifacts.len() != 2
        || manifest.output_artifacts.len() != 1
    {
        return Ok(false);
    }
    let inputs = vec![
        ManifestArtifact::snapshot("REFERENCE_INDEX", &reference.kallisto_index_path)?,
        ManifestArtifact::snapshot("EXTERNAL_ABUNDANCE", &dir.join("abundance.tsv"))?,
    ];
    let expected = request(run, &saved.parameters.provenance, &inputs);
    Ok(
        manifest.parameters["request"] == serde_json::to_value(&saved)?
            && saved.digest()?
                == request(run, &saved.parameters.provenance, &manifest.input_artifacts)
                    .digest()?
            && saved.digest()? == expected.digest()?
            && manifest.parameters["request_sha256"] == saved.digest()?
            && manifest.parameters["reference_index_sha256"] == inputs[0].checksum
            && manifest.output_artifacts[0] == inputs[1]
            && manifest.verify_local_artifacts()?,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use rnaa_core::{Database, ProjectConfig, ProjectPaths, traits::ReferenceManager};
    #[test]
    fn declarations_bind_payload_and_reference_without_claiming_native_quantification() -> Result<()>
    {
        let root = tempfile::tempdir()?;
        let paths = ProjectPaths::new(root.path());
        paths.ensure_layout()?;
        let cdna = root.path().join("source.fa");
        let gtf = root.path().join("source.gtf");
        fs::write(&cdna, format!(">tx1\n{}\n", "ACGT".repeat(25)))?;
        fs::write(
            &gtf,
            "chr1\ts\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\";\n",
        )?;
        let mut config = ProjectConfig::default();
        config.refs.custom_cdna = cdna.display().to_string();
        config.refs.custom_gtf = gtf.display().to_string();
        config.save(&paths.config_path)?;
        let db = Database::create(root.path(), &config)?;
        let reference =
            crate::refs::EnsemblReferenceManager::new()?.ensure_reference(&paths, &db, &config)?;
        let dir = root.path().join("provided");
        fs::create_dir_all(&dir)?;
        let abundance = dir.join("abundance.tsv");
        let data = "target_id\tlength\teff_length\test_counts\ttpm\ntx1\t100\t50\t10\t1000000\n";
        fs::write(&abundance, data)?;
        assert!(declare_external_quant("RUN1", &dir, &reference, "").is_err());
        declare_external_quant("RUN1", &dir, &reference, "explicit external fixture")?;
        assert!(external_quant_matches("RUN1", &dir, &reference)?);
        assert!(!external_quant_matches("OTHER", &dir, &reference)?);
        let mut forged: serde_json::Value = serde_json::from_slice(&fs::read(dir.join(FILE))?)?;
        forged["parameters"]["request"]["parameters"]["provenance"] =
            serde_json::json!("changed declaration without digest");
        write_json_pretty(&dir.join(FILE), &forged)?;
        assert!(!external_quant_matches("RUN1", &dir, &reference)?);
        declare_external_quant("RUN1", &dir, &reference, "explicit external fixture")?;
        let original = fs::read(dir.join(FILE))?;
        fs::write(&abundance, data.replace("\t10\t", "\t11\t"))?;
        assert!(!external_quant_matches("RUN1", &dir, &reference)?);
        declare_external_quant(
            "RUN1",
            &dir,
            &reference,
            "explicit revised external fixture",
        )?;
        assert_ne!(fs::read(dir.join(FILE))?, original);
        config.refs.index_k = 29;
        config.save(&paths.config_path)?;
        db.set_project_config(&config)?;
        let changed =
            crate::refs::EnsemblReferenceManager::new()?.ensure_reference(&paths, &db, &config)?;
        assert!(!external_quant_matches("RUN1", &dir, &changed)?);
        fs::write(&abundance, data.replace("tx1\t", "unknown\t"))?;
        assert!(declare_external_quant("RUN1", &dir, &changed, "wrong targets").is_err());
        fs::write(&abundance, data.replace("\t10\t", "\tNaN\t"))?;
        assert!(declare_external_quant("RUN1", &dir, &changed, "invalid values").is_err());
        Ok(())
    }
}
