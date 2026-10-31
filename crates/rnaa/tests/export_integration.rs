mod support;
use rnaa_core::{ProjectConfig, util::compute_sha256};
use std::{
    fs,
    path::{Path, PathBuf},
};
use support::native::{command, setup, success};

struct Worker(std::process::Child);
impl Drop for Worker {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

fn exported(root: &Path, args: &[&str]) -> PathBuf {
    let result = command(root, args);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    PathBuf::from(
        String::from_utf8(result.stdout)
            .unwrap()
            .lines()
            .find_map(|line| line.strip_prefix("export\t"))
            .unwrap(),
    )
}

fn rows(path: &Path, delimiter: u8) -> Vec<csv::ByteRecord> {
    csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(false)
        .from_path(path)
        .unwrap()
        .byte_records()
        .map(|row| row.unwrap())
        .collect()
}

#[path = "export_integration/current.rs"]
mod current;
#[path = "export_integration/custom_matrix.rs"]
mod custom_matrix;
#[path = "export_integration/publication.rs"]
mod publication;
