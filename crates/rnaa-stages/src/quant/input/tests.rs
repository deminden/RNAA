use super::*;
use crate::fastq::{FastqLimits, FastqPlan};
use anyhow::Context;
use kallistors::io::{PackedSeqBatchSource, ReadSource};
use rnaa_core::{ArtifactKind, LibraryLayout, VerifiedFile};
use std::fs;
fn text(records: usize, length: usize, header: usize) -> String {
    (0..records)
        .map(|i| {
            format!(
                "@read{i}{}\n{}\n+\n{}\n",
                "x".repeat(header),
                "A".repeat(length),
                "I".repeat(length)
            )
        })
        .collect()
}
fn plan(path: &Path, contents: &str) -> Result<FastqPlan> {
    fs::write(path, contents)?;
    FastqPlan::inspect(
        LibraryLayout::Single,
        &[VerifiedFile {
            artifact_kind: ArtifactKind::FastqSingle,
            path: path.into(),
            bytes: contents.len() as u64,
            checksum_type: "none".into(),
            checksum: String::new(),
            integrity_source: "fixture".into(),
        }],
        FastqLimits::default(),
        &|| Ok(()),
    )
}
#[test]
fn packed_and_record_paths_enforce_preflight_after_source_shape_changes() -> Result<()> {
    let root = tempfile::tempdir()?;
    let path = root.path().join("reads.fastq");
    let original = text(3, 100, 1000);
    let plan = plan(&path, &original)?;
    for (changed, reason) in [
        (text(4, 100, 0), "record count"),
        (text(3, 101, 1000), "sequence length"),
        (text(3, 99, 1000), "minimum sequence length"),
        (text(3, 100, 3000), "line bytes"),
        (
            text(3, 100, 1000).replace("\n+\n", &format!("\n+{}\n", "p".repeat(1000))),
            "record bytes",
        ),
    ] {
        fs::write(&path, changed)?;
        let mut packed = open(&path, plan.summary(), 0)?;
        let error = packed
            .next_packed_seq_batch(4096)
            .err()
            .context("shape change unexpectedly accepted")?;
        assert!(format!("{error:#}").contains(reason), "{error:#}");
        let mut records = open(&path, plan.summary(), 0)?;
        let mut error = None;
        while let Some(record) = records.next_record() {
            if let Err(e) = record {
                error = Some(e);
                break;
            }
        }
        assert!(error.is_some());
    }
    fs::write(&path, &original)?;
    let mut reader = open(&path, plan.summary(), 0)?;
    assert_eq!(reader.next_packed_seq_batch(4096)?.unwrap().len(), 3);
    assert!(reader.next_packed_seq_batch(4096)?.is_none());
    Ok(())
}

#[test]
fn total_decoded_growth_and_concatenated_gzip_members_are_bounded() -> Result<()> {
    use std::io::Write;
    let root = tempfile::tempdir()?;
    let path = root.path().join("reads.fastq");
    let original = text(2, 100, 0);
    let plan = plan(&path, &original)?;
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(original.as_bytes())?;
    let mut compressed = encoder.finish()?;
    let first = compressed.clone();
    fs::write(&path, &compressed)?;
    let mut reader = open(&path, plan.summary(), 0)?;
    assert_eq!(reader.next_packed_seq_batch(4096)?.unwrap().len(), 2);
    assert!(reader.next_packed_seq_batch(4096)?.is_none());
    compressed.extend(first);
    fs::write(&path, compressed)?;
    let mut reader = open(&path, plan.summary(), 0)?;
    let error = reader
        .next_packed_seq_batch(4096)
        .err()
        .context("gzip growth accepted")?;
    assert!(format!("{error:#}").contains("decoded bytes"));
    Ok(())
}
