use super::*;
use rnaa_core::resources::{ReservationOwner, ResourceConfig};
use std::io::Write;
fn gzip(text: &str) -> Result<Vec<u8>> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), Default::default());
    encoder.write_all(text.as_bytes())?;
    Ok(encoder.finish()?)
}

#[test]
fn ambiguity_and_poly_a_normalization_stay_within_the_native_graph_envelope() -> Result<()> {
    let root = tempfile::tempdir()?;
    let fasta = root.path().join("ambiguous.fa");
    let gtf = root.path().join("genes.gtf");
    let text = format!(
        ">tx1\n{}N{}NN{}\n",
        "acgu".repeat(30),
        "TGCA".repeat(30),
        "A".repeat(50)
    );
    std::fs::write(&fasta, text)?;
    std::fs::write(
        &gtf,
        "chr1\ttest\ttranscript\t1\t300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\";\n",
    )?;
    let context = ExecutionResources::open(
        &ResourceConfig {
            pool_root: root.path().join("pool").display().to_string(),
            cpu_threads: 2,
            memory_mb: 512,
            disk_headroom_mb: 0,
            wait_seconds: 5,
        },
        ReservationOwner {
            project_root: root.path().into(),
            attempt_id: "ambiguity".into(),
            stage: "refs".into(),
            run_accession: None,
        },
    )?;
    let plan = Plan::inspect(
        &ManifestArtifact::snapshot("REFERENCE_CDNA", &fasta)?,
        &ManifestArtifact::snapshot("REFERENCE_GTF", &gtf)?,
        &RefsConfig::default(),
        &context,
        root.path(),
        &|| Ok(()),
    )?;
    assert!(plan.fasta.uncertain_kmers > 0);
    for threads in [1, 2] {
        let report = kallistors::index::build_index_with_report(
            &root.path().join(format!("index-{threads}")),
            std::slice::from_ref(&fasta),
            kallistors::index::IndexBuildOptions {
                threads,
                ..Default::default()
            },
        )
        .map_err(|e| anyhow::anyhow!("{e}"))?;
        assert!(report.kmers as u64 <= plan.fasta.unique_kmers_upper);
        assert!(report.unitigs as u64 <= plan.fasta.unique_kmers_upper);
    }
    Ok(())
}
#[test]
fn full_gzip_sources_unique_cardinality_staging_and_native_shape_agree() -> Result<()> {
    let root = tempfile::tempdir()?;
    let fasta = root.path().join("transcripts.fa.gz");
    let gtf = root.path().join("genes.gtf.gz");
    let sequence = "ACGT".repeat(150);
    let first = format!(">tx1\r\n{sequence}\r\n");
    let second = format!(">tx2\n{}\n", sequence.replace('T', "U"));
    let mut bytes = gzip(&first)?;
    bytes.extend(gzip(&second)?);
    std::fs::write(&fasta, bytes)?;
    let annotation = "chr1\ttest\ttranscript\t1\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx1\";\nchr1\ttest\ttranscript\t1\t600\t.\t+\t.\tgene_id \"G2\"; transcript_id \"tx2\";\n";
    std::fs::write(&gtf, gzip(annotation)?)?;
    let context = ExecutionResources::open(
        &ResourceConfig {
            pool_root: root.path().join("pool").display().to_string(),
            cpu_threads: 2,
            memory_mb: 512,
            disk_headroom_mb: 0,
            wait_seconds: 5,
        },
        ReservationOwner {
            project_root: root.path().into(),
            attempt_id: "plan".into(),
            stage: "refs".into(),
            run_accession: None,
        },
    )?;
    let config = RefsConfig {
        index_threads: 4,
        ..Default::default()
    };
    let plan = Plan::inspect(
        &ManifestArtifact::snapshot("REFERENCE_CDNA", &fasta)?,
        &ManifestArtifact::snapshot("REFERENCE_GTF", &gtf)?,
        &config,
        &context,
        root.path(),
        &|| Ok(()),
    )?;
    assert_eq!(plan.fasta.transcripts, 2);
    assert_eq!(plan.fasta.bases, 1200);
    assert!(plan.fasta.clean_unique_kmers < 10);
    assert_eq!(plan.fasta.uncertain_kmers, 0);
    assert_eq!(plan.threads, 2);
    assert_eq!(plan.gtf.rows, 2);
    assert_eq!(plan.staging_bytes, (first.len() + second.len()) as u64);
    let plain = root.path().join("decoded.fa");
    plan.stage_fasta(&fasta, &plain, &|| Ok(()))?;
    assert_eq!(std::fs::read_to_string(&plain)?, first + &second);
    let index = root.path().join("reference.index");
    let report = kallistors::index::build_index_with_report(
        &index,
        &[plain],
        super::super::IndexSettings::from_config(&rnaa_core::ProjectConfig::default()).options(2),
    )
    .map_err(|e| anyhow::anyhow!("{e}"))?;
    assert_eq!(report.transcripts, 2);
    assert_eq!(report.kmers as u64, plan.fasta.clean_unique_kmers);
    assert!(std::fs::metadata(index)?.len() <= plan.index_bytes);
    assert!(context.pool().snapshot()?.reservations.is_empty());
    Ok(())
}
#[test]
fn bounded_decoded_lines_groups_and_cancellation_fail_before_growth() -> Result<()> {
    let root = tempfile::tempdir()?;
    let fasta = root.path().join("fasta");
    std::fs::write(&fasta, ">tx\nACGTACGTACGT\n")?;
    let check = || Ok(());
    let mut lines = Lines::new(
        &fasta,
        Limits {
            line: 8,
            decoded: None,
        },
        &check,
    )?;
    let mut line = Vec::new();
    assert!(lines.next(&mut line)?);
    assert!(lines.next(&mut line).is_err());
    let mut lines = Lines::new(
        &fasta,
        Limits {
            line: 64,
            decoded: Some(4),
        },
        &check,
    )?;
    assert!(lines.next(&mut line)?);
    assert!(lines.next(&mut line).is_err());
    let cancel = || bail!("cancelled");
    let mut lines = Lines::new(
        &fasta,
        Limits {
            line: 64,
            decoded: None,
        },
        &cancel,
    )?;
    assert!(lines.next(&mut line).is_err());
    Ok(())
}
