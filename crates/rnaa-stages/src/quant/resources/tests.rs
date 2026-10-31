use super::*;
use rnaa_core::{ArtifactKind, LibraryLayout, VerifiedFile};
use std::fs;

#[test]
fn decoded_request_accounts_for_short_read_jobs_tls_and_sequential_loading() -> Result<()> {
    let root = tempfile::tempdir()?;
    let seq = "ACGT".repeat(150);
    let fasta = root.path().join("transcripts.fa");
    fs::write(&fasta, format!(">base\n{seq}\n>other\n{seq}\n"))?;
    let path = root.path().join("reference.index");
    kallistors::index::build_index(
        &path,
        &[fasta],
        kallistors::index::IndexBuildOptions {
            k: 21,
            ..Default::default()
        },
    )
    .map_err(|e| anyhow::anyhow!("{e}"))?;
    let index = QuantIndexPlan::inspect(&path, &|| Ok(()))?;
    assert_eq!(index.max_ec_cardinality, 2);
    assert!(!index.has_shade_names);
    let mut files = Vec::new();
    for (mate, kind) in [(1, ArtifactKind::FastqR1), (2, ArtifactKind::FastqR2)] {
        let raw = root.path().join(format!("reads_R{mate}.fastq"));
        let record = format!(
            "@read/{mate}\n{}\n+\n{}\n",
            "ACGT".repeat(20),
            "I".repeat(80)
        );
        fs::write(&raw, record.repeat(100))?;
        files.push(VerifiedFile {
            artifact_kind: kind,
            path: raw,
            bytes: 0,
            checksum_type: "none".into(),
            checksum: String::new(),
            integrity_source: "fixture".into(),
        });
    }
    let plan = FastqPlan::inspect(
        LibraryLayout::Paired,
        &files,
        FastqLimits::default(),
        &|| Ok(()),
    )?;
    assert_eq!(plan.summary().min_sequence_bytes, 80);
    let serial = QuantResources::estimate(&index, &plan, 1, 1024 * 1024)?;
    let parallel = QuantResources::estimate(&index, &plan, 4, 1024 * 1024)?;
    assert_eq!(serial.ec_classes_upper, 3); // {base}, {other}, {base,other}
    assert!((4..=6).contains(&serial.ec_members_upper));
    assert_eq!(parallel.runtime_cache_bytes, 4 * 8 * 1024 * 1024);
    assert!(parallel.run_memory_bytes > serial.run_memory_bytes);
    assert_eq!(
        serial.index_loading_reservation_bytes + serial.run_memory_bytes,
        (serial.index_loading_bytes + serial.loading_setup_bytes)
            .max(serial.index_resident_bytes + serial.run_memory_bytes)
    );
    assert!(serial.index_loading_reservation_bytes >= serial.index_resident_bytes);
    assert!(serial.output_bytes > index.transcript_name_bytes);
    Ok(())
}

#[test]
fn shade_name_detection_survives_buffer_boundaries_and_uses_union_bounds() -> Result<()> {
    let root = tempfile::tempdir()?;
    let fasta = root.path().join("shade.fa");
    let name = format!("{}_shade_1", "b".repeat(65534));
    fs::write(
        &fasta,
        format!(
            ">base\n{}\n>{name}\n{}\n",
            "ACGT".repeat(150),
            "TGCA".repeat(150)
        ),
    )?;
    let index = root.path().join("shade.index");
    kallistors::index::build_index(
        &index,
        &[fasta],
        kallistors::index::IndexBuildOptions {
            k: 21,
            ..Default::default()
        },
    )
    .map_err(|e| anyhow::anyhow!("{e}"))?;
    let layout = QuantIndexPlan::inspect(&index, &|| Ok(()))?;
    assert!(layout.has_shade_names);
    assert!(layout.max_transcript_name_bytes > 65536);
    let raw = root.path().join("reads.fastq");
    fs::write(
        &raw,
        format!("@r\n{}\n+\n{}\n", "ACGT".repeat(20), "I".repeat(80)),
    )?;
    let plan = FastqPlan::inspect(
        LibraryLayout::Single,
        &[VerifiedFile {
            artifact_kind: ArtifactKind::FastqSingle,
            path: raw,
            bytes: 0,
            checksum_type: "none".into(),
            checksum: String::new(),
            integrity_source: "fixture".into(),
        }],
        FastqLimits::default(),
        &|| Ok(()),
    )?;
    let estimate = QuantResources::estimate(&layout, &plan, 1, 1024 * 1024)?;
    assert_eq!(estimate.ec_members_upper, 2);
    Ok(())
}
