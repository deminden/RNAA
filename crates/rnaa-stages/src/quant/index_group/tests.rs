use super::*;
use crate::quant::{RustKallistoQuantifier, quant_request_digest, snapshot_quant_inputs};
use rnaa_core::{
    LibraryLayout, ProjectConfig, ProjectPaths, RunRecord, RunState,
    model::{ArtifactKind, VerifiedFile},
    traits::Quantifier,
};
use std::{collections::BTreeMap, fs, sync::Barrier};
mod admission;

fn reference(root: &Path) -> Result<(ReferenceBundle, String)> {
    let mut state = 71_u64;
    let sequence = (0..1200)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[((state >> 32) & 3) as usize] as char
        })
        .collect::<String>();
    let fasta = root.join("transcripts.fa");
    fs::write(
        &fasta,
        format!(">tx1\n{sequence}\n>tx2\n{}\n", "ACGT".repeat(300)),
    )?;
    let index = root.join("reference.index");
    kallistors::index::build_index(
        &index,
        std::slice::from_ref(&fasta),
        kallistors::index::IndexBuildOptions {
            k: 21,
            ..Default::default()
        },
    )
    .map_err(|error| anyhow::anyhow!("{error}"))?;
    Ok((
        ReferenceBundle {
            id: "fixture".into(),
            organism: "synthetic".into(),
            ensembl_release: "none".into(),
            cdna_path: fasta,
            gtf_path: root.join("unused.gtf"),
            kallisto_index_path: index,
            tx2gene_path: root.join("unused-mapping.tsv"),
            gene_annotation_path: None,
            manifest_path: root.join("unused-manifest.json"),
        },
        sequence,
    ))
}
fn group(reference: &ReferenceBundle) -> Result<QuantIndexGroup> {
    QuantIndexGroup::new(ManifestArtifact::snapshot(
        "REFERENCE_INDEX",
        &reference.kallisto_index_path,
    )?)
}

#[test]
fn competing_workers_load_once_and_the_group_releases_its_index() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, _) = reference(root.path())?;
    let group = Arc::new(group(&reference)?);
    assert_eq!(group.stats().load_attempts, 0);
    assert_eq!(group.stats().successful_loads, 0);
    let barrier = Barrier::new(8);
    let handles = std::thread::scope(|scope| {
        let workers = (0..8)
            .map(|worker| {
                let group = &group;
                let barrier = &barrier;
                let reference = &reference;
                scope.spawn(move || {
                    barrier.wait();
                    group.load(reference, worker % 2 + 1)
                })
            })
            .collect::<Vec<_>>();
        workers
            .into_iter()
            .map(|worker| worker.join().unwrap())
            .collect::<Result<Vec<_>>>()
    })?;
    for index in &handles[1..] {
        assert!(Arc::ptr_eq(&handles[0], index));
    }
    assert_eq!(group.stats().load_attempts, 1);
    assert_eq!(group.stats().successful_loads, 1);
    let weak = Arc::downgrade(&handles[0]);
    drop(handles);
    assert!(weak.upgrade().is_some());
    drop(group);
    assert!(weak.upgrade().is_none());
    Ok(())
}

#[test]
fn cached_index_cannot_hide_same_size_corruption_or_a_different_reference() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, _) = reference(root.path())?;
    let group = group(&reference)?;
    let original = fs::read(&reference.kallisto_index_path)?;
    let loaded = group.load(&reference, 1)?;
    let mut corrupt = original.clone();
    *corrupt.last_mut().unwrap() ^= 1;
    fs::write(&reference.kallisto_index_path, corrupt)?;
    assert!(
        group
            .load(&reference, 1)
            .err()
            .unwrap()
            .to_string()
            .contains("frozen reference index changed")
    );
    fs::write(&reference.kallisto_index_path, original)?;
    assert!(Arc::ptr_eq(&loaded, &group.load(&reference, 2)?));
    let mut moved = reference.clone();
    moved.kallisto_index_path = root.path().join("another.index");
    fs::copy(&reference.kallisto_index_path, &moved.kallisto_index_path)?;
    assert!(group.load(&moved, 1).is_err());
    assert_eq!(group.stats().successful_loads, 1);
    Ok(())
}

#[test]
fn failed_loading_does_not_install_an_index_or_poison_later_requests() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, _) = reference(root.path())?;
    let good = fs::read(&reference.kallisto_index_path)?;
    fs::write(&reference.kallisto_index_path, b"invalid index")?;
    let failed = group(&reference)?;
    for attempt in 1..=2 {
        assert!(failed.load(&reference, 1).is_err());
        assert_eq!(failed.stats().load_attempts, attempt);
        assert_eq!(failed.stats().successful_loads, 0);
    }
    fs::write(&reference.kallisto_index_path, good)?;
    assert!(failed.load(&reference, 1).is_err());
    let replacement = group(&reference)?;
    assert!(replacement.load(&reference, 1).is_ok());
    assert_eq!(replacement.stats().successful_loads, 1);
    Ok(())
}

#[test]
fn group_construction_requires_a_valid_frozen_descriptor_but_does_not_load() -> Result<()> {
    let root = tempfile::tempdir()?;
    let path = root.path().join("not-an-index");
    fs::write(&path, b"not an index")?;
    let input = ManifestArtifact::snapshot("REFERENCE_INDEX", &path)?;
    let mut bad = input.clone();
    bad.kind = "OTHER".into();
    assert!(QuantIndexGroup::new(bad).is_err());
    let mut bad = input.clone();
    bad.checksum_type = "md5".into();
    assert!(QuantIndexGroup::new(bad).is_err());
    let mut bad = input.clone();
    bad.checksum = "x".repeat(64);
    assert!(QuantIndexGroup::new(bad).is_err());
    let mut bad = input.clone();
    bad.bytes = 0;
    assert!(QuantIndexGroup::new(bad).is_err());
    let group = QuantIndexGroup::new(input)?;
    assert_eq!(group.stats().load_attempts, 0);
    assert!(group.index.lock().unwrap().is_none());
    Ok(())
}

fn inputs(root: &Path, sequence: &str, paired: bool) -> Result<Vec<VerifiedFile>> {
    let mut first = String::new();
    let mut second = String::new();
    for read in 0..24 {
        let start = read * 11;
        let a = &sequence[start..start + 80];
        let b = sequence[start + 120..start + 200]
            .bytes()
            .rev()
            .map(|byte| match byte {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' => 'A',
                _ => unreachable!(),
            })
            .collect::<String>();
        first.push_str(&format!("@read{read}/1\n{a}\n+\n{}\n", "I".repeat(80)));
        second.push_str(&format!("@read{read}/2\n{b}\n+\n{}\n", "I".repeat(80)));
    }
    let files = if paired {
        vec![
            ("reads_R1.fastq", ArtifactKind::FastqR1, first),
            ("reads_R2.fastq", ArtifactKind::FastqR2, second),
        ]
    } else {
        vec![("reads.fastq", ArtifactKind::FastqSingle, first)]
    };
    files
        .into_iter()
        .map(|(name, artifact_kind, text)| {
            let path = root.join(name);
            fs::write(&path, text)?;
            let input = ManifestArtifact::snapshot(artifact_kind.as_str(), &path)?;
            Ok(VerifiedFile {
                artifact_kind,
                path,
                checksum_type: input.checksum_type,
                checksum: input.checksum,
                bytes: input.bytes,
                integrity_source: "fixture".into(),
            })
        })
        .collect()
}

#[test]
fn memory_rejection_and_owned_cancellation_precede_index_loading() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, sequence) = reference(root.path())?;
    let group = group(&reference)?;
    let files = inputs(root.path(), &sequence, false)?;
    let run = RunRecord {
        project_id: "test".into(),
        study_accession: "test".into(),
        source_study_accession: None,
        run_accession: "RUN".into(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Single,
        instrument_platform: None,
        instrument_model: None,
        remote_files: Vec::new(),
        metadata: Default::default(),
        state: RunState::Verified,
        last_error: None,
        updated_at: "fixture".into(),
    };
    let paths = ProjectPaths::new(root.path().join("outputs"));
    let mut config = ProjectConfig::default();
    config.resources.memory_mb = 1;
    config.quant.match_cache_mb = 1;
    let error = group
        .quantify(&run, &files, &reference, &paths, &config)
        .err()
        .unwrap();
    assert!(
        error
            .downcast_ref::<rnaa_core::resources::ResourceAdmissionError>()
            .is_some()
    );
    assert_eq!(group.stats().load_attempts, 0);
    assert_eq!(group.stats().successful_loads, 0);
    assert!(!paths.root.exists());
    config.resources.memory_mb = 512;
    let error = group
        .quantify_checked(&run, &files, &reference, &paths, &config, &|| {
            anyhow::bail!("cancelled owner")
        })
        .err()
        .unwrap();
    assert!(format!("{error:#}").contains("cancelled owner"));
    assert_eq!(group.stats().load_attempts, 0);
    assert!(!paths.root.exists());
    group.quantify(&run, &files, &reference, &paths, &config)?;
    assert_eq!(group.stats().successful_loads, 1);
    assert_eq!(group.stats().layout_plan.unwrap().transcripts, 2);
    Ok(())
}

#[test]
fn shared_mapping_matches_independent_loads_for_layout_strand_and_thread_options() -> Result<()> {
    let root = tempfile::tempdir()?;
    let (reference, sequence) = reference(root.path())?;
    let group = group(&reference)?;
    let single = inputs(root.path(), &sequence, false)?;
    let paired = inputs(root.path(), &sequence, true)?;
    let barrier = Barrier::new(6);
    std::thread::scope(|scope| -> Result<()> {
        let workers = (0..6)
            .map(|worker| {
                let group = &group;
                let reference = &reference;
                let barrier = &barrier;
                let paired_layout = worker >= 3;
                let inputs = if paired_layout { &paired } else { &single };
                let root = root.path();
                scope.spawn(move || -> Result<()> {
                    let mut config = ProjectConfig::default();
                    config.quant.threads = worker % 2 + 1;
                    config.quant.match_cache_mb = 1;
                    config.quant.strand = [
                        rnaa_core::config::QuantStrand::Unstranded,
                        rnaa_core::config::QuantStrand::Forward,
                        rnaa_core::config::QuantStrand::Reverse,
                    ][worker % 3];
                    let run = RunRecord {
                        project_id: "fixture".into(),
                        study_accession: "fixture".into(),
                        source_study_accession: None,
                        run_accession: format!("RUN{worker}"),
                        sample_accession: None,
                        experiment_accession: None,
                        library_layout: if paired_layout {
                            LibraryLayout::Paired
                        } else {
                            LibraryLayout::Single
                        },
                        instrument_platform: None,
                        instrument_model: None,
                        remote_files: vec![],
                        metadata: BTreeMap::new(),
                        state: RunState::Verified,
                        last_error: None,
                        updated_at: "fixture".into(),
                    };
                    barrier.wait();
                    let shared = group.quantify(
                        &run,
                        inputs,
                        reference,
                        &ProjectPaths::new(root.join(format!("shared-{worker}"))),
                        &config,
                    )?;
                    let mut baseline = config.clone();
                    baseline.quant.match_cache_mb = if worker == 0 { 512 } else { 8 };
                    let separate = RustKallistoQuantifier.quantify(
                        &run,
                        inputs,
                        reference,
                        &ProjectPaths::new(root.join(format!("separate-{worker}"))),
                        &baseline,
                    )?;
                    assert_eq!(
                        fs::read(shared.out_dir.join("abundance.tsv"))?,
                        fs::read(separate.out_dir.join("abundance.tsv"))?
                    );
                    let mut a: serde_json::Value =
                        serde_json::from_slice(&fs::read(&shared.run_info_json)?)?;
                    let mut b: serde_json::Value =
                        serde_json::from_slice(&fs::read(&separate.run_info_json)?)?;
                    a.as_object_mut().unwrap().remove("start_time");
                    b.as_object_mut().unwrap().remove("start_time");
                    assert_eq!(a, b);
                    if worker == 0 {
                        assert_eq!(a["n_pseudoaligned"], 24);
                    }
                    let shared_manifest: serde_json::Value =
                        serde_json::from_slice(&fs::read(shared.manifest_path)?)?;
                    assert_eq!(
                        shared_manifest["parameters"]["index_execution"]["group_id"],
                        group.stats().group_id
                    );
                    assert!(
                        shared_manifest["parameters"]["match_cache"]["allocated_bytes"]
                            .as_u64()
                            .unwrap()
                            <= 1024 * 1024
                    );
                    let frozen = snapshot_quant_inputs(inputs)?;
                    let key = quant_request_digest(&run, reference, &config, &frozen)?;
                    config.quant.threads += 1;
                    config.quant.match_cache_mb = 32;
                    assert_eq!(
                        quant_request_digest(&run, reference, &config, &frozen)?,
                        key
                    );
                    Ok(())
                })
            })
            .collect::<Vec<_>>();
        for worker in workers {
            worker.join().unwrap()?;
        }
        Ok(())
    })?;
    assert_eq!(group.stats().load_attempts, 1);
    assert_eq!(group.stats().successful_loads, 1);
    Ok(())
}
