use super::annotation::normalize_stable_id;
use super::source::{
    ORGANISM_PRESETS, extract_gencode_release, matching_release_files, normalize_gencode_release,
    provider_md5,
};
use super::*;
use rnaa_core::util::compute_sha256;
use std::path::Path;

fn build_annotation_tables(source: &Path, tx: &Path, genes: &Path) -> Result<()> {
    let limits = super::planning::Limits::try_from(&rnaa_core::config::RefsConfig::default())?;
    let plan = super::planning::scan_gtf(source, limits, &|| Ok(()))?;
    super::annotation::build_annotation_tables(source, tx, genes, &plan, &|| Ok(()))
}

#[test]
fn admitted_complete_gzip_reference_publishes_all_members_under_one_build_owner() -> Result<()> {
    use std::io::Write;
    let (root, paths, db, mut config) = local_project()?;
    let sequence = "ACGT".repeat(150);
    let mut encoded = Vec::new();
    for (name, text) in [
        ("tx1", sequence.clone()),
        ("tx2", sequence.replace('T', "U")),
    ] {
        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), Default::default());
        gzip.write_all(format!(">{name}\n{text}\n").as_bytes())?;
        encoded.extend(gzip.finish()?);
    }
    let fasta = root.path().join("full.fa.gz");
    fs::write(&fasta, encoded)?;
    config.refs.custom_cdna = fasta.display().to_string();
    fs::write(
        &config.refs.custom_gtf,
        "chr1\ttest\ttranscript\t1\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\";\nchr1\ttest\ttranscript\t1\t600\t.\t+\t.\tgene_id \"G2\"; transcript_id \"tx2\";\n",
    )?;
    config.resources.pool_root = root.path().join("pool").display().to_string();
    config.resources.memory_mb = 512;
    config.resources.cpu_threads = 2;
    config.resources.disk_headroom_mb = 0;
    config.refs.index_threads = 4;
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    let result = EnsemblReferenceManager::new()?.ensure_reference(&paths, &db, &config)?;
    let index =
        kallistors::index::Index::load(&result.kallisto_index_path).map_err(|e| anyhow!("{e}"))?;
    assert_eq!(
        index
            .transcripts
            .iter()
            .map(|t| t.name.as_str())
            .collect::<Vec<_>>(),
        vec!["tx1", "tx2"]
    );
    assert_eq!(
        fs::read_to_string(&result.tx2gene_path)?,
        "transcript_id\tgene_id\ntx1\tG1\ntx2\tG2\n"
    );
    let manifest: StageManifest = serde_json::from_slice(&fs::read(&result.manifest_path)?)?;
    assert_eq!(manifest.parameters["request"]["algorithm_version"], 2);
    assert_eq!(manifest.parameters["mapping"]["mapped_transcripts"], 2);
    assert_eq!(manifest.parameters["index_threads"], 2);
    assert_eq!(manifest.parameters["requested_index_threads"], 4);
    assert_eq!(
        manifest.parameters["resource_execution"]["build_reservation"]["owner"]["stage"],
        "refs.build"
    );
    assert!(
        db.verify_generation_head(&db.generation_head(OutputGroup::Refs, &result.id)?.unwrap())?
    );
    assert!(
        rnaa_core::resources::ResourcePool::open(Path::new(&config.resources.pool_root))?
            .snapshot()?
            .reservations
            .is_empty()
    );
    Ok(())
}

#[test]
fn build_capacity_rejection_keeps_prior_bundle_and_historical_reuse_precedes_pool_open()
-> Result<()> {
    let (root, paths, db, mut config) = local_project()?;
    config.resources.pool_root = root.path().join("pool").display().to_string();
    config.resources.memory_mb = 256;
    config.resources.cpu_threads = 1;
    config.resources.disk_headroom_mb = 0;
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    let manager = EnsemblReferenceManager::new()?;
    let first = manager.ensure_reference(&paths, &db, &config)?;
    let head = db.generation_head(OutputGroup::Refs, &first.id)?.unwrap();
    let before = fs::read(&first.manifest_path)?;
    config.refs.index_k = 29;
    config.resources.memory_mb = 100;
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    let error = manager.ensure_reference(&paths, &db, &config).unwrap_err();
    assert!(
        format!("{error:#}").contains("resource admission failed"),
        "{error:#}"
    );
    assert_eq!(
        db.generation_head(OutputGroup::Refs, &first.id)?.unwrap(),
        head
    );
    assert_eq!(fs::read(&first.manifest_path)?, before);
    assert!(
        rnaa_core::resources::ResourcePool::open(Path::new(&config.resources.pool_root))?
            .snapshot()?
            .reservations
            .is_empty()
    );
    config.resources.memory_mb = 256;
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    let second = manager.ensure_reference(&paths, &db, &config)?;
    assert_ne!(first.manifest_path, second.manifest_path);
    let unavailable = root.path().join("not-a-directory");
    fs::write(&unavailable, b"held")?;
    config.resources.pool_root = unavailable.display().to_string();
    config.resources.memory_mb = 1;
    config.refs.index_k = 31;
    config.refs.max_decoded_mb = 1;
    config.save(&paths.config_path)?;
    db.set_project_config(&config)?;
    let restored = manager.ensure_reference(&paths, &db, &config)?;
    assert_eq!(first.manifest_path, restored.manifest_path);
    assert_eq!(
        db.generation_head(OutputGroup::Refs, &first.id)?.unwrap(),
        head
    );
    assert_eq!(fs::read(&first.manifest_path)?, before);
    Ok(())
}

#[test]
fn extracts_gencode_release_name() {
    assert_eq!(
        extract_gencode_release("gencode.v49.annotation.gtf.gz").as_deref(),
        Some("v49")
    );
}

#[test]
fn finds_gencode_transcript_and_gtf_from_listings() {
    let listing = r#"<a href="gencode.v49.transcripts.fa.gz">tx</a><a href="gencode.v49.annotation.gtf.gz">gtf</a>"#;
    let (tx, gtf, release) = matching_release_files(listing, "gencode.v", None).unwrap();
    assert_eq!(tx, "gencode.v49.transcripts.fa.gz");
    assert_eq!(gtf, "gencode.v49.annotation.gtf.gz");
    assert_eq!(release, "v49");
}
#[test]
fn normalizes_requested_gencode_release() {
    assert_eq!(
        normalize_gencode_release("49", ORGANISM_PRESETS[0]).unwrap(),
        Some("v49".to_string())
    );
    assert_eq!(
        normalize_gencode_release("M38", ORGANISM_PRESETS[1]).unwrap(),
        Some("M38".to_string())
    );
}

#[test]
fn strips_versions_from_stable_ids() {
    assert_eq!(normalize_stable_id("ENST00000335137.4"), "ENST00000335137");
    assert_eq!(normalize_stable_id("ENSG00000186092.7"), "ENSG00000186092");
    assert_eq!(normalize_stable_id("ENSG00000186092"), "ENSG00000186092");
}
fn local_project() -> Result<(tempfile::TempDir, ProjectPaths, Database, ProjectConfig)> {
    let root = tempfile::tempdir()?;
    let paths = ProjectPaths::new(root.path());
    paths.ensure_layout()?;
    let cdna = root.path().join("input.fa");
    let gtf = root.path().join("input.gtf");
    fs::write(
        &cdna,
        format!(
            ">tx1.2|G1.1|label|\n{}\n",
            "ACGTTCGATCGGATCGATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCTAGCTGACTGATCGATCGTAGCTAGCTGATCGTAGCTAGCTAGCATCGA"
        ),
    )?;
    fs::write(
        &gtf,
        "chr1\tsynthetic\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1.1\"; transcript_id \"tx1.2\"; gene_name \"Gene1\"; gene_type \"protein_coding\";\n",
    )?;
    let mut config = ProjectConfig::default();
    config.refs.custom_cdna = cdna.display().to_string();
    config.refs.custom_gtf = gtf.display().to_string();
    let db = Database::create(root.path(), &config)?;
    Ok((root, paths, db, config))
}
#[test]
fn reference_reuse_verifies_all_payloads_and_binds_sources_and_scientific_settings() -> Result<()> {
    let (_root, paths, db, mut config) = local_project()?;
    let manager = EnsemblReferenceManager::new()?;
    // Legacy loadable files and annotations have no proof of their inputs.
    let spec = resolve_reference_spec(&manager.client, &config)?;
    fs::create_dir_all(paths.reference_dir(&spec.reference_id))?;
    let legacy = paths.reference_dir(&spec.reference_id).join("tx2gene.tsv");
    fs::write(&legacy, "obsolete\n")?;
    let first = manager.ensure_reference(&paths, &db, &config)?;
    let head = db.generation_head(OutputGroup::Refs, &first.id)?.unwrap();
    assert!(db.verify_generation_head(&head)?);
    assert_eq!(db.list_artifacts()?.len(), 6);
    assert_eq!(fs::read_to_string(&legacy)?, "obsolete\n");
    let manifest = fs::read(&first.manifest_path)?;
    let attempts = db.list_job_attempts(100)?.len();
    assert_eq!(
        fs::read_to_string(&first.tx2gene_path)?,
        "transcript_id\tgene_id\ntx1\tG1\n"
    );
    config.refs.index_threads = 2;
    let reused = manager.ensure_reference(&paths, &db, &config)?;
    assert_eq!(reused.manifest_path, first.manifest_path);
    assert_eq!(fs::read(&first.manifest_path)?, manifest);
    assert_eq!(db.list_job_attempts(100)?.len(), attempts);
    let old_index = fs::read(&first.kallisto_index_path)?;
    config.refs.index_k = 29;
    db.set_project_config(&config)?;
    let second = manager.ensure_reference(&paths, &db, &config)?;
    assert_ne!(second.manifest_path, first.manifest_path);
    assert_eq!(second.id, first.id);
    assert_eq!(fs::read(&first.kallisto_index_path)?, old_index);
    assert!(db.verify_generation_head(&head)?);
    assert_eq!(
        kallistors::index::Index::load(&second.kallisto_index_path)
            .unwrap()
            .k,
        29
    );
    // A syntactically valid index built for another option is still invalid here.
    fs::copy(&first.kallisto_index_path, &second.kallisto_index_path)?;
    assert!(kallistors::index::Index::load(&second.kallisto_index_path).is_ok());
    let repaired = manager.ensure_reference(&paths, &db, &config)?;
    assert_ne!(repaired.manifest_path, second.manifest_path);
    assert_eq!(
        kallistors::index::Index::load(&repaired.kallisto_index_path)
            .unwrap()
            .k,
        29
    );
    fs::write(&repaired.tx2gene_path, "transcript_id\tgene_id\ntx1\tBAD\n")?;
    let annotations = manager.ensure_reference(&paths, &db, &config)?;
    assert_ne!(annotations.manifest_path, repaired.manifest_path);
    assert_eq!(
        fs::read_to_string(&annotations.tx2gene_path)?,
        "transcript_id\tgene_id\ntx1\tG1\n"
    );
    let source = Path::new(&config.refs.custom_cdna);
    let text = fs::read_to_string(source)?;
    fs::write(source, text.replacen("ACGTTCGA", "TCGTTCGA", 1))?;
    let changed = manager.ensure_reference(&paths, &db, &config)?;
    assert_ne!(changed.id, first.id);
    assert_ne!(
        compute_sha256(&changed.cdna_path)?,
        compute_sha256(&first.cdna_path)?
    );
    assert!(
        db.list_job_attempts(100)?
            .iter()
            .all(|job| job.stage == "refs"
                && job.state == "SUCCEEDED"
                && job.request_sha256.is_some())
    );
    Ok(())
}
#[test]
fn changing_sources_and_cancelled_parents_cannot_publish_a_reference_replacement() -> Result<()> {
    let (_root, paths, db, config) = local_project()?;
    let manager = EnsemblReferenceManager::new()?;
    let first = manager.ensure_reference(&paths, &db, &config)?;
    let head = db.generation_head(OutputGroup::Refs, &first.id)?.unwrap();
    let manifest = fs::read(&first.manifest_path)?;
    let spec = resolve_reference_spec(&manager.client, &config)?;
    let mut lease = JobLease::try_acquire_reference(&db, &spec.reference_id)?.unwrap();
    let expected = db.selected_reference()?.unwrap().selection_token()?;
    let preparation = Preparation {
        spec: &spec,
        expected: Some(&expected),
        observation: ReferenceObservation::capture(&db, &config)?,
    };
    fs::write(&config.refs.custom_gtf, "changed after resolution")?;
    let error = manager
        .build_reference(
            &paths,
            &lease.database(),
            &db,
            &config,
            &preparation,
            lease.attempt_id(),
        )
        .unwrap_err();
    assert!(format!("{error:#}").contains("source changed"));
    lease.finish(JobOutcome::Failed, Some("source changed"))?;
    drop(lease);
    assert_eq!(
        db.generation_head(OutputGroup::Refs, &first.id)?.unwrap(),
        head
    );
    assert_eq!(fs::read(&first.manifest_path)?, manifest);
    let mut parent = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    let mut child =
        JobLease::try_acquire_reference(&parent.database(), &spec.reference_id)?.unwrap();
    assert!(db.cancel_job(parent.attempt_id())?);
    let error = manager
        .build_reference(
            &paths,
            &child.database(),
            &parent.database(),
            &config,
            &preparation,
            child.attempt_id(),
        )
        .unwrap_err();
    assert!(format!("{error:#}").contains("cancelled"));
    child.finish(JobOutcome::Failed, Some("parent cancelled"))?;
    parent.finish(JobOutcome::Failed, Some("cancelled"))?;
    assert!(db.verify_generation_head(&head)?);
    assert_eq!(
        db.generation_head(OutputGroup::Refs, &first.id)?.unwrap(),
        head
    );
    assert_eq!(db.list_artifacts()?.len(), 6);
    Ok(())
}
#[test]
fn provider_pairing_checksums_and_release_paths_cannot_mix_versions() -> Result<()> {
    let listing = r#"<a href="gencode.v48.transcripts.fa.gz">old</a><a href="gencode.v49.transcripts.fa.gz">new</a><a href="gencode.v49.annotation.gtf.gz">annotation</a>"#;
    let (fasta, gtf, release) = matching_release_files(listing, "gencode.v", None)?;
    assert_eq!(release, "v49");
    assert_eq!(fasta, "gencode.v49.transcripts.fa.gz");
    assert_eq!(gtf, "gencode.v49.annotation.gtf.gz");
    assert!(matching_release_files(listing, "gencode.v", Some("v48")).is_err());
    let mixed = r#"href="gencode.vM38.transcripts.fa.gz" href="gencode.vM37.annotation.gtf.gz""#;
    assert!(matching_release_files(mixed, "gencode.vM", None).is_err());
    assert!(normalize_gencode_release("../../bad", ORGANISM_PRESETS[0]).is_err());
    let checksum = "a".repeat(32);
    let line = format!("{checksum}  {fasta}\n");
    assert_eq!(provider_md5(&line, &fasta)?, checksum);
    assert!(provider_md5(&line.repeat(2), &fasta).is_err());
    assert!(provider_md5("bad  x.gz", "x.gz").is_err());
    assert!(provider_md5(&line, "missing.gz").is_err());
    Ok(())
}
#[test]
fn annotation_derivation_rejects_malformed_and_conflicting_normalized_mappings() -> Result<()> {
    let (root, paths, _db, config) = local_project()?;
    let tx = root.path().join("table.tsv");
    let genes = root.path().join("genes.tsv");
    let source = Path::new(&config.refs.custom_gtf);
    fs::write(source, "literal\\ttabs\\tand no records")?;
    assert!(build_annotation_tables(source, &tx, &genes).is_err());
    fs::write(
        source,
        "chr1\ts\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1.1\";\nchr1\ts\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G2\"; transcript_id \"tx1.2\";\n",
    )?;
    assert!(build_annotation_tables(source, &tx, &genes).is_err());
    assert!(paths.generations_dir.is_dir());
    Ok(())
}
