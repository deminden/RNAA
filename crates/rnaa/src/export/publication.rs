//! Bounded table writing and no-replacement publication.
use anyhow::{Context, Result, bail};
use rnaa_core::ManifestArtifact;
use sha2::{Digest, Sha256};
use std::{
    fs::{File, OpenOptions},
    io::{Read, Write},
    path::Path,
};

pub(super) fn publish(
    db: &rnaa_core::Database,
    source: &Path,
    destination: &Path,
    complete: impl FnOnce(&rnaa_core::Database) -> Result<()>,
) -> Result<()> {
    let mut renamed = false;
    let result = db.atomic(|transaction| {
        transaction.check_job()?;
        rename_exclusive(source, destination)?;
        renamed = true;
        complete(transaction)
    });
    if let Err(error) = result {
        if renamed {
            rename_exclusive(destination, source).with_context(|| {
                format!(
                    "export finalization failed ({error:#}); complete output may remain at {}",
                    destination.display()
                )
            })?;
        }
        return Err(error);
    }
    Ok(())
}

struct Budget<W> {
    inner: W,
    remaining: u64,
}

struct VerifiedInput<R> {
    inner: R,
    hash: Sha256,
    bytes: u64,
}

impl<R: Read> Read for VerifiedInput<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let count = self.inner.read(buffer)?;
        self.hash.update(&buffer[..count]);
        self.bytes = self
            .bytes
            .checked_add(count as u64)
            .ok_or_else(|| std::io::Error::other("export input byte count overflow"))?;
        Ok(count)
    }
}

impl<R> VerifiedInput<R> {
    fn verify(self, source: &ManifestArtifact) -> Result<()> {
        if self.bytes != source.bytes || hex::encode(self.hash.finalize()) != source.checksum {
            bail!("export source bytes changed while copying: {}", source.path);
        }
        Ok(())
    }
}
impl<W: Write> Write for Budget<W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() as u64 > self.remaining {
            return Err(std::io::Error::other(
                "export output exceeded its reserved byte limit",
            ));
        }
        let count = self.inner.write(bytes)?;
        self.remaining -= count as u64;
        Ok(count)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

pub(super) fn write_json(
    path: &Path,
    value: &serde_json::Value,
    limit: u64,
    check: &impl Fn() -> Result<()>,
) -> Result<()> {
    struct CheckedWriter<'a, W, F> {
        inner: W,
        check: &'a F,
        unchecked: usize,
    }
    impl<W: Write, F: Fn() -> Result<()>> Write for CheckedWriter<'_, W, F> {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if self.unchecked >= 1024 * 1024 {
                (self.check)().map_err(|error| std::io::Error::other(format!("{error:#}")))?;
                self.unchecked = 0;
            }
            let count = self.inner.write(&bytes[..bytes.len().min(64 * 1024)])?;
            self.unchecked = self.unchecked.saturating_add(count);
            Ok(count)
        }
        fn flush(&mut self) -> std::io::Result<()> {
            self.inner.flush()
        }
    }
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    let mut output = CheckedWriter {
        inner: Budget {
            inner: file,
            remaining: limit,
        },
        check,
        unchecked: 1024 * 1024,
    };
    serde_json::to_writer_pretty(&mut output, value)?;
    output.write_all(b"\n")?;
    output.flush()?;
    output.inner.inner.sync_all()?;
    Ok(())
}

pub(super) fn write(
    source: &ManifestArtifact,
    target: &Path,
    delimiter: u8,
    provenance: bool,
    limit: u64,
    check: &impl Fn() -> Result<()>,
) -> Result<(Option<u64>, Option<u64>)> {
    if !source.matches_file()? {
        bail!("export source changed: {}", source.path);
    }
    let input = VerifiedInput {
        inner: File::open(&source.path)?.take(source.bytes.saturating_add(1)),
        hash: Sha256::new(),
        bytes: 0,
    };
    let output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(target)?;
    let mut output = Budget {
        inner: output,
        remaining: limit,
    };
    let (rows, columns) = if provenance {
        let mut input = input;
        let mut buffer = [0_u8; 65536];
        let mut unchecked = 8 * 1024 * 1024_u64;
        loop {
            if unchecked >= 8 * 1024 * 1024 {
                check()?;
                unchecked = 0;
            }
            let count = input.read(&mut buffer)?;
            if count == 0 {
                break;
            }
            output.write_all(&buffer[..count])?;
            unchecked += count as u64;
        }
        input.verify(source)?;
        (None, None)
    } else {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(input);
        let mut writer = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .has_headers(false)
            .from_writer(output);
        let mut records = 0_u64;
        let mut unchecked = 0_u64;
        let mut columns = None;
        for row in reader.byte_records() {
            if records.is_multiple_of(1024) || unchecked >= 8 * 1024 * 1024 {
                check()?;
                unchecked = 0;
            }
            let row = row?;
            unchecked = unchecked
                .saturating_add(row.as_slice().len() as u64)
                .saturating_add(row.len() as u64);
            columns.get_or_insert(row.len() as u64);
            writer.write_byte_record(&row)?;
            records = records
                .checked_add(1)
                .context("export row count overflow")?;
        }
        reader.into_inner().verify(source)?;
        writer.flush()?;
        output = writer.into_inner().map_err(|error| error.into_error())?;
        if records == 0 {
            bail!("cannot export an empty table: {}", source.path);
        }
        (Some(records - 1), columns)
    };
    output.flush()?;
    output.inner.sync_all()?;
    check()?;
    if !source.matches_file()? {
        bail!("export source changed while copying: {}", source.path);
    }
    Ok((rows, columns))
}

#[cfg(any(target_os = "macos", target_os = "linux"))]
pub(super) fn rename_exclusive(source: &Path, destination: &Path) -> Result<()> {
    use std::{ffi::CString, os::unix::ffi::OsStrExt};
    let source = CString::new(source.as_os_str().as_bytes())?;
    let destination = CString::new(destination.as_os_str().as_bytes())?;
    // SAFETY: both paths are valid NUL-terminated strings; the operating system
    // performs the no-replacement rename on the same filesystem.
    #[cfg(target_os = "macos")]
    let result =
        unsafe { libc::renamex_np(source.as_ptr(), destination.as_ptr(), libc::RENAME_EXCL) };
    // SAFETY: same path contract as above; AT_FDCWD resolves relative paths.
    #[cfg(target_os = "linux")]
    let result = unsafe {
        libc::renameat2(
            libc::AT_FDCWD,
            source.as_ptr(),
            libc::AT_FDCWD,
            destination.as_ptr(),
            libc::RENAME_NOREPLACE,
        )
    };
    if result != 0 {
        return Err(std::io::Error::last_os_error())
            .context("cannot publish export without replacement");
    }
    Ok(())
}

#[cfg(not(any(target_os = "macos", target_os = "linux")))]
pub(super) fn rename_exclusive(_source: &Path, _destination: &Path) -> Result<()> {
    bail!("atomic export publication is supported on macOS and Linux")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn json_manifest_cannot_exceed_its_reserved_budget() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let value = serde_json::json!({"metadata": "x".repeat(1024)});
        let complete = temp.path().join("complete.json");
        write_json(&complete, &value, 4096, &|| Ok(()))?;
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&fs::read(&complete)?)?,
            value
        );

        let limited = temp.path().join("limited.json");
        let error = write_json(&limited, &value, 64, &|| Ok(())).unwrap_err();
        assert!(format!("{error:#}").contains("reserved byte limit"));
        assert!(fs::metadata(&limited)?.len() <= 64);
        assert!(write_json(&complete, &serde_json::json!({}), 4096, &|| Ok(())).is_err());
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&fs::read(&complete)?)?,
            value
        );
        Ok(())
    }

    #[test]
    fn json_manifest_checks_cancellation_while_writing() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let value = serde_json::json!({"metadata": "x".repeat(3 * 1024 * 1024)});
        let checks = std::cell::Cell::new(0);
        let error = write_json(
            &temp.path().join("cancelled.json"),
            &value,
            4 * 1024 * 1024,
            &|| {
                checks.set(checks.get() + 1);
                if checks.get() > 1 {
                    bail!("cancelled");
                }
                Ok(())
            },
        )
        .unwrap_err();
        assert!(format!("{error:#}").contains("cancelled"));
        assert_eq!(checks.get(), 2);
        Ok(())
    }

    #[test]
    fn conversion_preserves_quoted_fields_and_enforces_output_budget() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let input = temp.path().join("input.tsv");
        fs::write(&input, b"id\tvalue\nrow\t\"comma, and \"\"quote\"\"\"\n")?;
        let source = ManifestArtifact::snapshot("TEST", &input)?;
        let target = temp.path().join("out.csv");
        assert_eq!(
            write(&source, &target, b',', false, 1024, &|| Ok(()))?,
            (Some(1), Some(2))
        );
        let mut reader = csv::Reader::from_path(target)?;
        let row = reader.records().next().context("missing data row")??;
        assert_eq!(&row[1], "comma, and \"quote\"");
        assert!(
            write(
                &source,
                &temp.path().join("limited.csv"),
                b',',
                false,
                3,
                &|| Ok(())
            )
            .is_err()
        );
        Ok(())
    }

    #[test]
    fn cancellation_is_checked_during_table_conversion() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let input = temp.path().join("input.tsv");
        let mut data = String::from("id\tvalue\n");
        for index in 0..3000 {
            data.push_str(&format!("{index}\tvalue\n"));
        }
        fs::write(&input, data.as_bytes())?;
        let source = ManifestArtifact::snapshot("TEST", &input)?;
        let checks = std::cell::Cell::new(0);
        let result = write(
            &source,
            &temp.path().join("staged.csv"),
            b',',
            false,
            1024 * 1024,
            &|| {
                checks.set(checks.get() + 1);
                if checks.get() == 2 {
                    bail!("cancelled");
                }
                Ok(())
            },
        );
        assert!(format!("{:#}", result.unwrap_err()).contains("cancelled"));
        assert_eq!(checks.get(), 2);
        assert_eq!(fs::read(&input)?, data.as_bytes());
        Ok(())
    }

    #[test]
    fn conversion_rejects_temporary_source_changes_even_after_restoration() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let input = temp.path().join("input.tsv");
        let mut original = String::from("id\tvalue\n");
        for index in 0..3000 {
            original.push_str(&format!("{index}\toriginal\n"));
        }
        let modified = original.replace("original", "modified");
        fs::write(&input, &original)?;
        let source = ManifestArtifact::snapshot("TEST", &input)?;
        let checks = std::cell::Cell::new(0);
        let result = write(
            &source,
            &temp.path().join("staged.csv"),
            b',',
            false,
            1024 * 1024,
            &|| {
                checks.set(checks.get() + 1);
                match checks.get() {
                    1 => fs::write(&input, &modified)?,
                    2 => fs::write(&input, &original)?,
                    _ => {}
                }
                Ok(())
            },
        );
        assert!(format!("{:#}", result.unwrap_err()).contains("source bytes changed"));
        assert!(source.matches_file()?);
        Ok(())
    }

    #[cfg(any(target_os = "macos", target_os = "linux"))]
    #[test]
    fn exclusive_publication_never_replaces_an_existing_directory() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let source = temp.path().join("source");
        let destination = temp.path().join("destination");
        fs::create_dir(&source)?;
        fs::create_dir(&destination)?;
        fs::write(source.join("new"), b"new")?;
        fs::write(destination.join("old"), b"old")?;
        assert!(rename_exclusive(&source, &destination).is_err());
        assert_eq!(fs::read(source.join("new"))?, b"new");
        assert_eq!(fs::read(destination.join("old"))?, b"old");
        Ok(())
    }

    #[cfg(any(target_os = "macos", target_os = "linux"))]
    #[test]
    fn failed_database_finalization_restores_the_staged_directory() -> Result<()> {
        use rnaa_core::jobs::{JobLease, JobOutcome};
        let temp = tempfile::tempdir()?;
        let db = rnaa_core::Database::create(temp.path(), &rnaa_core::ProjectConfig::default())?;
        let mut lease =
            JobLease::try_acquire(&db, "export", None, &[])?.context("missing lease")?;
        let source = temp.path().join("staged");
        let destination = temp.path().join("published");
        fs::create_dir(&source)?;
        fs::write(source.join("result"), b"content")?;
        let result = publish(&lease.database(), &source, &destination, |_| {
            bail!("injected transaction failure")
        });
        assert!(result.is_err());
        assert!(!destination.exists());
        assert_eq!(fs::read(source.join("result"))?, b"content");
        lease.finish(JobOutcome::Failed, Some("injected"))?;
        Ok(())
    }

    #[cfg(any(target_os = "macos", target_os = "linux"))]
    #[test]
    fn process_exit_after_rename_rolls_back_job_completion() -> Result<()> {
        use rnaa_core::jobs::{JobLease, JobOutcome};
        const CHILD_ROOT: &str = "RNAA_TEST_PUBLICATION_CRASH_ROOT";
        if let Some(root) = std::env::var_os(CHILD_ROOT) {
            let root = std::path::PathBuf::from(root);
            let db = rnaa_core::Database::open(&root)?;
            let lease =
                JobLease::try_acquire(&db, "export", None, &[])?.context("missing lease")?;
            let source = root.join("staged");
            fs::create_dir(&source)?;
            fs::write(source.join("result"), b"content")?;
            publish(
                &lease.database(),
                &source,
                &root.join("published"),
                |transaction| {
                    transaction.complete_owned_job(JobOutcome::Succeeded, None)?;
                    std::process::exit(73);
                },
            )?;
            unreachable!();
        }
        let temp = tempfile::tempdir()?;
        let db = rnaa_core::Database::create(temp.path(), &rnaa_core::ProjectConfig::default())?;
        let status = std::process::Command::new(std::env::current_exe()?)
            .arg("--exact")
            .arg("export::publication::tests::process_exit_after_rename_rolls_back_job_completion")
            .env(CHILD_ROOT, temp.path())
            .status()?;
        assert_eq!(status.code(), Some(73));
        assert!(!temp.path().join("staged").exists());
        assert_eq!(fs::read(temp.path().join("published/result"))?, b"content");
        let attempts = db.list_job_attempts(10)?;
        assert_eq!(attempts.len(), 1);
        assert_eq!(attempts[0].stage, "export");
        assert_eq!(attempts[0].state, "RUNNING");
        Ok(())
    }
}
