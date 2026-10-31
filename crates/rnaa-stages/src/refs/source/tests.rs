use super::*;
#[test]
fn reference_subjects_own_independent_transfer_caches_and_cleanup_cannot_remove_another_input()
-> Result<()> {
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    let spec = SourceSpec {
        kind: SourceKind::RemoteUrl,
        location: "http://origin/shared.fa".into(),
        gzip: false,
        sha256: None,
        md5: None,
    };
    let mut sources = Vec::new();
    for subject in ["reference-A", "reference-B"] {
        let path = transfer_path(&paths, &spec, subject)?;
        fs::create_dir_all(path.parent().unwrap())?;
        fs::write(&path, b"identical source")?;
        sources.push(MaterializedSource {
            input: ManifestArtifact::snapshot("REFERENCE_CDNA", &path)?,
            http: None,
            transfer: Some(path),
            execution: None,
        });
    }
    assert_ne!(sources[0].input.path, sources[1].input.path);
    sources[0].discard_transfer()?;
    assert!(sources[1].input.matches_file()?);
    let copied = root.path().join("generation-source.fa");
    sources[1].copy_to(&copied, &|| Ok(()))?;
    assert_eq!(fs::read(copied)?, b"identical source");
    sources[1].discard_transfer()?;
    Ok(())
}
