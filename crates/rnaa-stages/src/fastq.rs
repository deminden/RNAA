//! Bounded decoded FASTQ admission and normalized staging for complete groups.
#[cfg(test)]
mod admission_tests;
mod layout;
mod plan;
mod reader;
use anyhow::{Context, Result, bail};
pub use plan::FastqPlan;
pub(crate) use reader::open as open_reader;
use rnaa_core::{
    LibraryLayout,
    config::QuantConfig,
    model::{ArtifactKind, VerifiedFile},
};
use std::path::PathBuf;

pub(crate) fn scan_request(limits: FastqLimits) -> Result<rnaa_core::resources::ResourceRequest> {
    let memory_bytes = (limits.max_line_bytes as u64)
        .checked_mul(16)
        .and_then(|n| n.checked_add(16 * 1024 * 1024))
        .context("FASTQ preflight memory estimate overflows")?;
    Ok(rnaa_core::resources::ResourceRequest {
        cpu_threads: 1,
        memory_bytes,
        disks: Vec::new(),
    })
}

#[derive(Debug, Clone, Copy)]
pub struct FastqLimits {
    pub max_line_bytes: usize,
    pub max_decoded_bytes: Option<u64>,
    pub max_staging_bytes: Option<u64>,
}
impl Default for FastqLimits {
    fn default() -> Self {
        Self {
            max_line_bytes: 1024 * 1024,
            max_decoded_bytes: None,
            max_staging_bytes: None,
        }
    }
}
impl TryFrom<&QuantConfig> for FastqLimits {
    type Error = anyhow::Error;
    fn try_from(config: &QuantConfig) -> Result<Self> {
        Ok(Self {
            max_line_bytes: config.fastq_max_line_bytes()?,
            max_decoded_bytes: config.fastq_max_decoded_bytes()?,
            max_staging_bytes: None,
        })
    }
}
impl FastqLimits {
    fn validate(self) -> Result<()> {
        if self.max_line_bytes == 0 || self.max_line_bytes >= isize::MAX as usize {
            bail!("FASTQ line limit must be positive and addressable");
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize)]
pub struct FastqSummary {
    pub lanes: usize,
    pub paired: bool,
    pub fragments: u64,
    pub read_records: u64,
    pub decoded_bytes: u64,
    pub r1_normalized_bytes: u64,
    pub r2_normalized_bytes: u64,
    pub max_line_bytes: usize,
    pub max_record_bytes: u64,
    pub min_record_bytes: u64,
    pub max_sequence_bytes: usize,
    pub min_sequence_bytes: usize,
    pub needs_staging: bool,
}
impl FastqSummary {
    pub fn normalized_bytes(&self) -> Result<u64> {
        self.r1_normalized_bytes
            .checked_add(self.r2_normalized_bytes)
            .context("normalized FASTQ byte count overflows")
    }
    pub fn staging_bytes(&self) -> Result<u64> {
        if self.needs_staging {
            self.normalized_bytes()
        } else {
            Ok(0)
        }
    }
}
pub struct PreparedFastqs {
    pub read1: PathBuf,
    pub read2: Option<PathBuf>,
    pub read1_kind: ArtifactKind,
    pub summary: FastqSummary,
    _staging: Option<tempfile::TempDir>,
}

pub fn prepare_fastqs(layout: LibraryLayout, files: &[VerifiedFile]) -> Result<PreparedFastqs> {
    FastqPlan::inspect(layout, files, FastqLimits::default(), &|| Ok(()))?
        .prepare_in(&std::env::temp_dir(), &|| Ok(()))
}
pub fn validate_fastqs(
    layout: LibraryLayout,
    files: &[VerifiedFile],
    check: &impl Fn() -> Result<()>,
) -> Result<()> {
    FastqPlan::inspect(layout, files, FastqLimits::default(), check).map(|_| ())
}

#[cfg(test)]
use std::{fs, io::Write, path::Path};
#[cfg(test)]
mod tests {
    use super::*;
    fn file(root: &Path, name: &str, kind: ArtifactKind, records: &str) -> VerifiedFile {
        let path = root.join(name);
        fs::write(&path, records).unwrap();
        VerifiedFile {
            artifact_kind: kind,
            path,
            checksum_type: "none".into(),
            checksum: String::new(),
            bytes: records.len() as u64,
            integrity_source: "test".into(),
        }
    }
    #[test]
    fn consumes_all_single_end_lanes_with_missing_final_newline() {
        let root = tempfile::tempdir().unwrap();
        let a = file(
            root.path(),
            "a.fastq",
            ArtifactKind::FastqSingle,
            "@a\nACGT\n+\nIIII",
        );
        let b = file(
            root.path(),
            "b.fastq",
            ArtifactKind::FastqSingle,
            "@b\nTGCA\n+\nIIII\n",
        );
        let prepared = prepare_fastqs(LibraryLayout::Single, &[b, a]).unwrap();
        assert_eq!(
            fs::read_to_string(&prepared.read1).unwrap(),
            "@a\nACGT\n+\nIIII\n@b\nTGCA\n+\nIIII\n"
        );
    }
    #[test]
    fn matches_paired_lanes_and_consumes_all_mates() {
        let root = tempfile::tempdir().unwrap();
        let mut files = Vec::new();
        for lane in [2, 1] {
            for (mate, kind) in [(1, ArtifactKind::FastqR1), (2, ArtifactKind::FastqR2)] {
                files.push(file(
                    root.path(),
                    &format!("L{lane}_R{mate}.fastq"),
                    kind,
                    &format!("@read{lane}/{mate}\nACGT\n+\nIIII\n"),
                ));
            }
        }
        let prepared = prepare_fastqs(LibraryLayout::Paired, &files).unwrap();
        for path in [&prepared.read1, prepared.read2.as_ref().unwrap()] {
            let data = fs::read_to_string(path).unwrap();
            assert!(data.find("@read1").unwrap() < data.find("@read2").unwrap());
        }
    }
    #[test]
    fn rejects_missing_mates_duplicate_inputs_and_wrong_layout() {
        let root = tempfile::tempdir().unwrap();
        let a = file(
            root.path(),
            "R1.fastq",
            ArtifactKind::FastqR1,
            "@a/1\nACGT\n+\nIIII\n",
        );
        assert!(prepare_fastqs(LibraryLayout::Paired, std::slice::from_ref(&a)).is_err());
        assert!(prepare_fastqs(LibraryLayout::Single, std::slice::from_ref(&a)).is_err());
        assert!(prepare_fastqs(LibraryLayout::Single, &[a.clone(), a]).is_err());
    }
    #[test]
    fn rejects_unequal_counts_mismatched_ids_and_truncated_reads() {
        let root = tempfile::tempdir().unwrap();
        let a = file(
            root.path(),
            "R1.fastq",
            ArtifactKind::FastqR1,
            "@a/1\nACGT\n+\nIIII\n",
        );
        for content in [
            "@b/2\nACGT\n+\nIIII\n",
            "@a/2\nACGT\n+\nIIII\n@b/2\nACGT\n+\nIIII\n",
            "@a/2\nACGT\n+\n",
        ] {
            let b = file(root.path(), "R2.fastq", ArtifactKind::FastqR2, content);
            assert!(prepare_fastqs(LibraryLayout::Paired, &[a.clone(), b]).is_err());
        }
    }
    #[test]
    fn reads_concatenated_gzip_members() {
        use flate2::{Compression, write::GzEncoder};
        let root = tempfile::tempdir().unwrap();
        let mut a = file(root.path(), "reads.gz", ArtifactKind::FastqSingle, "");
        let mut bytes = Vec::new();
        for id in ["a", "b"] {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder
                .write_all(format!("@{id}\nACGT\n+\nIIII\n").as_bytes())
                .unwrap();
            bytes.extend(encoder.finish().unwrap());
        }
        fs::write(&a.path, bytes).unwrap();
        a.bytes = fs::metadata(&a.path).unwrap().len();
        assert!(prepare_fastqs(LibraryLayout::Single, &[a]).is_ok());
    }
}
