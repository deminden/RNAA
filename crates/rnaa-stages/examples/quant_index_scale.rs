//! Compare retained independent indexes with one shared worker-group index.
//! Build the reference in a separate process so build RSS cannot mask load RSS.
use anyhow::{Context, Result, bail};
use rnaa_core::{
    LibraryLayout, ManifestArtifact, ProjectConfig, ProjectPaths, ReferenceBundle, RunRecord,
    RunState,
    model::{ArtifactKind, VerifiedFile},
    traits::Quantifier,
};
use rnaa_stages::quant::QuantIndexGroup;
use std::{
    collections::{BTreeMap, BTreeSet},
    env, fs,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    sync::{Arc, Barrier},
    time::Instant,
};

fn build(root: &Path, transcripts: usize, length: usize) -> Result<()> {
    if transcripts == 0 || length < 250 {
        bail!("require transcripts > 0 and length >= 250");
    }
    fs::create_dir_all(root)?;
    let fasta = root.join("transcripts.fa");
    let mut out = BufWriter::new(fs::File::create(&fasta)?);
    let mut random = 17_u64;
    for transcript in 0..transcripts {
        writeln!(out, ">tx{transcript}")?;
        for _ in 0..length {
            random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
            out.write_all(&[b"ACGT"[((random >> 32) & 3) as usize]])?;
        }
        writeln!(out)?;
    }
    out.flush()?;
    kallistors::index::build_index(&root.join("reference.index"), &[fasta], Default::default())
        .map_err(|error| anyhow::anyhow!("index build failed: {error}"))?;
    println!(
        "{}",
        serde_json::json!({"transcripts":transcripts,"length":length})
    );
    Ok(())
}

fn measure(root: &Path, workers: usize, mode: &str, cache_mb: usize) -> Result<()> {
    if workers == 0 || workers > 64 || !matches!(mode, "shared" | "independent") {
        bail!("require mode shared|independent and workers in 1..=64");
    }
    let fasta = root.join("transcripts.fa");
    let sequence = BufReader::new(fs::File::open(&fasta)?)
        .lines()
        .nth(1)
        .context("missing sequence")??;
    let reference = ReferenceBundle {
        id: "scale".into(),
        organism: "synthetic".into(),
        ensembl_release: "none".into(),
        cdna_path: fasta,
        gtf_path: root.join("unused.gtf"),
        kallisto_index_path: root.join("reference.index"),
        tx2gene_path: root.join("unused.tsv"),
        gene_annotation_path: None,
        manifest_path: root.join("unused.json"),
    };
    let input = ManifestArtifact::snapshot("REFERENCE_INDEX", &reference.kallisto_index_path)?;
    let output = tempfile::tempdir()?;
    let raw = output.path().join("reads.fastq");
    let mut reads = String::new();
    for read in 0..64 {
        reads.push_str(&format!(
            "@read{read}\n{}\n+\n{}\n",
            &sequence[..100],
            "I".repeat(100)
        ));
    }
    fs::write(&raw, reads)?;
    let fastq = ManifestArtifact::snapshot("FASTQ_SINGLE", &raw)?;
    let files = [VerifiedFile {
        artifact_kind: ArtifactKind::FastqSingle,
        path: raw,
        checksum_type: fastq.checksum_type,
        checksum: fastq.checksum,
        bytes: fastq.bytes,
        integrity_source: "synthetic".into(),
    }];
    let groups = if mode == "shared" {
        let shared = Arc::new(QuantIndexGroup::new(input.clone())?);
        vec![shared; workers]
    } else {
        (0..workers)
            .map(|_| QuantIndexGroup::new(input.clone()).map(Arc::new))
            .collect::<Result<Vec<_>>>()?
    };
    let barrier = Barrier::new(workers);
    let started = Instant::now();
    let hashes = std::thread::scope(|scope| {
        let handles = groups
            .iter()
            .enumerate()
            .map(|(worker, group)| {
                let reference = &reference;
                let files = &files;
                let barrier = &barrier;
                let out = output.path().join(format!("worker-{worker}"));
                scope.spawn(move || -> Result<String> {
                    let mut config = ProjectConfig::default();
                    config.quant.threads = 1;
                    config.quant.match_cache_mb = cache_mb;
                    let run = RunRecord {
                        project_id: "scale".into(),
                        study_accession: "scale".into(),
                        source_study_accession: None,
                        run_accession: format!("RUN{worker}"),
                        sample_accession: None,
                        experiment_accession: None,
                        library_layout: LibraryLayout::Single,
                        instrument_platform: None,
                        instrument_model: None,
                        remote_files: vec![],
                        metadata: BTreeMap::new(),
                        state: RunState::Verified,
                        last_error: None,
                        updated_at: "scale".into(),
                    };
                    barrier.wait();
                    let result =
                        group.quantify(&run, files, reference, &ProjectPaths::new(out), &config)?;
                    let info: serde_json::Value =
                        serde_json::from_slice(&fs::read(result.run_info_json)?)?;
                    if info["n_pseudoaligned"] != 64 {
                        bail!("scale fixture did not map all reads: {info}");
                    }
                    rnaa_core::util::compute_sha256(&result.out_dir.join("abundance.tsv"))
                })
            })
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .map(|handle| {
                handle
                    .join()
                    .map_err(|_| anyhow::anyhow!("scale worker panicked"))?
            })
            .collect::<Result<BTreeSet<_>>>()
    })?;
    if hashes.len() != 1 {
        bail!("worker output hashes differ");
    }
    let unique = groups
        .iter()
        .map(|group| (group.stats().group_id, group.stats()))
        .collect::<BTreeMap<_, _>>();
    println!(
        "{}",
        serde_json::json!({
            "mode":mode,"workers":workers,"index_serialized_bytes":input.bytes,
            "match_cache_mb_per_run":cache_mb,
            "index_sha256":input.checksum,"loaded_index_groups":unique.len(),
            "successful_index_loads":unique.values().map(|stats|stats.successful_loads).sum::<usize>(),
            "load_and_quant_wall_seconds":started.elapsed().as_secs_f64(),
            "abundance_sha256":hashes.first(),"pseudoaligned_reads_per_worker":64,
        })
    );
    Ok(())
}
fn main() -> Result<()> {
    let args = env::args().skip(1).collect::<Vec<_>>();
    match args.as_slice() {
        [mode, root, transcripts, length] if mode == "build" => {
            build(Path::new(root), transcripts.parse()?, length.parse()?)
        }
        [mode, root, workers, cache] => {
            measure(Path::new(root), workers.parse()?, mode, cache.parse()?)
        }
        [mode, root, workers] => measure(
            Path::new(root),
            workers.parse()?,
            mode,
            ProjectConfig::default().quant.match_cache_mb,
        ),
        _ => bail!(
            "usage: quant_index_scale build DIR TRANSCRIPTS LENGTH | shared|independent DIR WORKERS [CACHE_MB]"
        ),
    }
}
