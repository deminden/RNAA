use anyhow::{Context, Result};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};
pub fn reference(name: &str) -> Result<PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../validation/scientific")
        .join(name);
    let metadata: serde_json::Value =
        serde_json::from_slice(&fs::read(root.join("provenance.json"))?)?;
    assert_eq!(metadata["schema_version"], 1);
    if let Some(generator) = metadata["generator"]["path"].as_str() {
        assert_eq!(
            Some(
                rnaa_core::util::compute_sha256(
                    &root
                        .parent()
                        .context("reference parent missing")?
                        .join(generator)
                )?
                .as_str()
            ),
            metadata["generator"]["sha256"].as_str(),
            "reference generator integrity"
        );
    }
    if let Some(inputs) = metadata["generator"]["inputs"].as_object() {
        for (path, expected) in inputs {
            let input = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../..")
                .join(path);
            assert_eq!(
                Some(rnaa_core::util::compute_sha256(&input)?.as_str()),
                expected["sha256"].as_str(),
                "oracle generator input integrity: {path}"
            );
        }
    }
    for (file, proof) in metadata["files"]
        .as_object()
        .context("reference files missing")?
    {
        let observed = rnaa_core::util::compute_sha256(&root.join(file))?;
        assert_eq!(
            Some(observed.as_str()),
            proof["sha256"].as_str(),
            "reference file integrity: {file}"
        );
    }
    Ok(root)
}
pub fn read_table(path: &Path) -> Result<Vec<BTreeMap<String, String>>> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let header = reader.headers()?.clone();
    reader
        .records()
        .map(|row| {
            let row = row?;
            Ok(header
                .iter()
                .zip(&row)
                .map(|(key, value)| (key.into(), value.into()))
                .collect())
        })
        .collect()
}
pub fn assert_close(actual: f64, expected: f64, absolute: f64, relative: f64, label: &str) {
    let limit = absolute + relative * expected.abs();
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= limit,
        "{label}: actual={actual:.17e}, expected={expected:.17e}, error={:.3e}, allowed={limit:.3e}",
        (actual - expected).abs()
    );
}
