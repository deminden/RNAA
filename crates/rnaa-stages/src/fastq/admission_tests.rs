use super::*;
use anyhow::Result;
use flate2::{Compression, write::GzEncoder};
use std::{cell::Cell, fs, io::Write, path::Path};

fn input(root: &Path, name: &str, kind: ArtifactKind, text: &[u8]) -> VerifiedFile {
    let path = root.join(name);
    fs::write(&path, text).unwrap();
    VerifiedFile {
        artifact_kind: kind,
        path,
        checksum_type: "none".into(),
        checksum: String::new(),
        bytes: text.len() as u64,
        integrity_source: "fixture".into(),
    }
}
fn gzip(text: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(text).unwrap();
    encoder.finish().unwrap()
}
fn reads(name: &str, records: usize, length: usize) -> String {
    let seq = "ACGT".repeat(length.div_ceil(4));
    let seq = &seq[..length];
    (0..records)
        .map(|read| format!("@{name}{read}\n{seq}\n+\n{}\n", "I".repeat(length)))
        .collect()
}

#[test]
fn preflight_counts_all_members_lanes_and_exact_normalized_staging_before_writes() -> Result<()> {
    let root = tempfile::tempdir()?;
    let staging = root.path().join("staging");
    fs::create_dir(&staging)?;
    let first = "@read0/1\r\nACGT\r\n+\r\nIIII\r\n";
    let second = "@read0/2\nTGCA\n+\nIIII";
    let mut members = gzip(first.as_bytes());
    members.extend(gzip(first.replace("read0", "read1").as_bytes()));
    let files = vec![
        input(
            root.path(),
            "lane0_R1.fq.gz",
            ArtifactKind::FastqR1,
            &members,
        ),
        input(
            root.path(),
            "lane0_R2.fq",
            ArtifactKind::FastqR2,
            second.as_bytes(),
        ),
    ];
    assert!(
        FastqPlan::inspect(
            LibraryLayout::Paired,
            &files,
            FastqLimits::default(),
            &|| Ok(())
        )
        .is_err()
    );
    let mut files = files;
    let right = format!("{second}\n{}", second.replace("read0", "read1"));
    fs::write(&files[1].path, &right)?;
    files[1].bytes = right.len() as u64;
    let plan = FastqPlan::inspect(
        LibraryLayout::Paired,
        &files,
        FastqLimits::default(),
        &|| Ok(()),
    )?;
    assert_eq!(plan.summary().fragments, 2);
    assert_eq!(plan.summary().read_records, 4);
    assert_eq!(
        plan.summary().decoded_bytes,
        (first.len() * 2 + right.len()) as u64
    );
    assert!(plan.summary().needs_staging);
    assert_eq!(fs::read_dir(&staging)?.count(), 0);
    let prepared = plan.prepare_in(&staging, &|| Ok(()))?;
    let a = fs::read(&prepared.read1)?;
    let b = fs::read(prepared.read2.as_ref().unwrap())?;
    assert_eq!(plan.summary().staging_bytes()?, (a.len() + b.len()) as u64);
    assert!(!a.contains(&b'\r') && !b.contains(&b'\r'));
    assert_eq!(a.iter().filter(|&&byte| byte == b'\n').count(), 8);
    drop(prepared);
    assert_eq!(fs::read_dir(&staging)?.count(), 0);
    Ok(())
}

#[test]
fn decoded_limit_uses_the_whole_group_instead_of_each_compressed_file() -> Result<()> {
    let root = tempfile::tempdir()?;
    let text = reads("read", 4000, 100);
    let files = [
        input(
            root.path(),
            "lane_a.gz",
            ArtifactKind::FastqSingle,
            &gzip(text.as_bytes()),
        ),
        input(
            root.path(),
            "lane_b.gz",
            ArtifactKind::FastqSingle,
            &gzip(text.as_bytes()),
        ),
    ];
    assert!(files.iter().all(|file| file.bytes < 1024 * 1024));
    assert!(text.len() < 1024 * 1024 && text.len() * 2 > 1024 * 1024);
    let error = FastqPlan::inspect(
        LibraryLayout::Single,
        &files,
        FastqLimits {
            max_decoded_bytes: Some(1024 * 1024),
            ..Default::default()
        },
        &|| Ok(()),
    )
    .unwrap_err();
    assert!(format!("{error:#}").contains("decoded FASTQ group exceeds"));
    let plan = FastqPlan::inspect(
        LibraryLayout::Single,
        &files,
        FastqLimits::default(),
        &|| Ok(()),
    )?;
    assert_eq!(plan.summary().decoded_bytes, (text.len() * 2) as u64);
    assert_eq!(plan.summary().lanes, 2);
    Ok(())
}

#[test]
fn line_and_staging_limits_reject_before_temp_payloads_exist() -> Result<()> {
    let root = tempfile::tempdir()?;
    let long = input(
        root.path(),
        "long.gz",
        ArtifactKind::FastqSingle,
        &gzip(reads("read", 1, 5000).as_bytes()),
    );
    let error = FastqPlan::inspect(
        LibraryLayout::Single,
        &[long],
        FastqLimits {
            max_line_bytes: 128,
            ..Default::default()
        },
        &|| Ok(()),
    )
    .unwrap_err();
    assert!(format!("{error:#}").contains("max_line_bytes=128"));
    let text = reads("read", 2, 40);
    let files = [
        input(
            root.path(),
            "a.fastq",
            ArtifactKind::FastqSingle,
            text.as_bytes(),
        ),
        input(
            root.path(),
            "b.fastq",
            ArtifactKind::FastqSingle,
            text.as_bytes(),
        ),
    ];
    assert!(
        FastqPlan::inspect(
            LibraryLayout::Single,
            &files,
            FastqLimits {
                max_staging_bytes: Some(text.len() as u64),
                ..Default::default()
            },
            &|| Ok(())
        )
        .is_err()
    );
    Ok(())
}

#[test]
fn mutation_before_or_during_materialization_cannot_outgrow_the_preflight_promise() -> Result<()> {
    let root = tempfile::tempdir()?;
    let staging = root.path().join("staging");
    fs::create_dir(&staging)?;
    let text = reads("read", 2, 40);
    let files = [
        input(
            root.path(),
            "a.fastq",
            ArtifactKind::FastqSingle,
            text.as_bytes(),
        ),
        input(
            root.path(),
            "b.fastq",
            ArtifactKind::FastqSingle,
            text.as_bytes(),
        ),
    ];
    let plan = FastqPlan::inspect(
        LibraryLayout::Single,
        &files,
        FastqLimits::default(),
        &|| Ok(()),
    )?;
    fs::write(&files[0].path, text.replace("ACGT", "TGCA"))?;
    assert!(plan.prepare_in(&staging, &|| Ok(())).is_err());
    assert_eq!(fs::read_dir(&staging)?.count(), 0);
    fs::write(&files[0].path, &text)?;
    let checks = Cell::new(0);
    let error = plan
        .prepare_in(&staging, &|| {
            checks.set(checks.get() + 1);
            if checks.get() == 3 {
                let mut file = fs::OpenOptions::new().append(true).open(&files[0].path)?;
                file.write_all(text.repeat(20).as_bytes())?;
            }
            Ok(())
        })
        .err()
        .unwrap();
    assert!(format!("{error:#}").contains("decoded FASTQ group exceeds"));
    assert_eq!(fs::read_dir(&staging)?.count(), 0);
    Ok(())
}

#[test]
fn cancellation_and_gzip_corruption_leave_no_staged_files() -> Result<()> {
    let root = tempfile::tempdir()?;
    let bytes = gzip(reads("read", 5000, 40).as_bytes());
    let file = input(root.path(), "reads.gz", ArtifactKind::FastqSingle, &bytes);
    let checks = Cell::new(0);
    let result = FastqPlan::inspect(
        LibraryLayout::Single,
        std::slice::from_ref(&file),
        FastqLimits::default(),
        &|| {
            checks.set(checks.get() + 1);
            if checks.get() == 3 {
                anyhow::bail!("cancelled during validation");
            }
            Ok(())
        },
    );
    assert!(format!("{:#}", result.unwrap_err()).contains("cancelled during validation"));
    let mut corrupt = bytes;
    *corrupt.last_mut().unwrap() ^= 1;
    fs::write(&file.path, corrupt)?;
    assert!(
        FastqPlan::inspect(
            LibraryLayout::Single,
            &[file],
            FastqLimits::default(),
            &|| Ok(())
        )
        .is_err()
    );
    Ok(())
}
