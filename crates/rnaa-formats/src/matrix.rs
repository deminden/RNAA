use std::fs::File;
use std::io::{Cursor, Read};
use std::path::Path;

use anyhow::{Context, Result, bail};
use csv::{ReaderBuilder, WriterBuilder};
use flate2::read::GzDecoder;
use tar::Archive;

#[derive(Debug, Clone)]
pub struct MatrixData {
    pub row_ids: Vec<String>,
    pub col_ids: Vec<String>,
    pub values: Vec<Vec<f64>>,
}

impl MatrixData {
    pub fn validate(&self) -> Result<()> {
        use std::collections::BTreeSet;
        if self.row_ids.is_empty()
            || self.col_ids.is_empty()
            || self.row_ids.len() != self.values.len()
        {
            bail!("matrix must have nonempty rows/columns and one value row per gene ID");
        }
        for (label, ids) in [("gene", &self.row_ids), ("sample", &self.col_ids)] {
            let unique = ids.iter().collect::<BTreeSet<_>>();
            if unique.len() != ids.len() || ids.iter().any(|id| id.trim().is_empty()) {
                bail!("matrix contains duplicate or empty {label} IDs");
            }
        }
        for row in &self.values {
            if row.len() != self.col_ids.len() || row.iter().any(|value| !value.is_finite()) {
                bail!("matrix values must be finite and match the sample dimensions");
            }
        }
        Ok(())
    }
}

pub fn read_matrix_tsv(path: &Path) -> Result<MatrixData> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    read_matrix_from_reader(file)
}

/// Read only sample identifiers, rejecting files that differ from the admitted size.
pub fn read_matrix_columns_exact(path: &Path, admitted_bytes: u64) -> Result<Vec<String>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let metadata = file.metadata()?;
    if !metadata.is_file() || metadata.len() != admitted_bytes {
        bail!("matrix size changed after resource admission");
    }
    let limit = admitted_bytes
        .checked_add(1)
        .context("matrix size overflow")?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.take(limit));
    let columns = reader
        .headers()?
        .iter()
        .skip(1)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    if reader.position().byte() > admitted_bytes
        || reader.get_ref().get_ref().metadata()?.len() != admitted_bytes
    {
        bail!("matrix size changed after resource admission");
    }
    let unique = columns.iter().collect::<std::collections::BTreeSet<_>>();
    if columns.is_empty()
        || unique.len() != columns.len()
        || columns.iter().any(|id| id.trim().is_empty())
    {
        bail!("matrix contains duplicate or empty sample IDs");
    }
    Ok(columns)
}

pub fn read_matrix_maybe_gzip(path: &Path) -> Result<MatrixData> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    if path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        read_matrix_from_reader(GzDecoder::new(file))
    } else {
        read_matrix_from_reader(file)
    }
}

pub fn read_first_matrix_from_tar_gz(path: &Path) -> Result<MatrixData> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let gz = GzDecoder::new(file);
    let mut archive = Archive::new(gz);
    for entry in archive.entries()? {
        let mut entry = entry?;
        if entry.header().entry_type().is_file() {
            let mut buf = Vec::new();
            entry.read_to_end(&mut buf)?;
            return read_matrix_from_reader(Cursor::new(buf));
        }
    }
    bail!("no matrix file found in {}", path.display())
}

pub fn write_matrix_tsv(path: &Path, matrix: &MatrixData) -> Result<()> {
    matrix.validate()?;
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut header = vec!["gene_id".to_string()];
    header.extend(matrix.col_ids.iter().cloned());
    writer.write_record(header)?;

    for (row_id, values) in matrix.row_ids.iter().zip(&matrix.values) {
        if values.len() != matrix.col_ids.len() {
            bail!(
                "row {} has {} values but header has {} columns",
                row_id,
                values.len(),
                matrix.col_ids.len()
            );
        }
        let mut record = vec![row_id.clone()];
        record.extend(values.iter().map(|value| value.to_string()));
        writer.write_record(record)?;
    }
    writer.flush()?;
    Ok(())
}

fn read_matrix_from_reader<R: Read>(reader: R) -> Result<MatrixData> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(reader);
    let headers = reader.headers()?.clone();
    let col_ids = headers.iter().skip(1).map(ToString::to_string).collect();

    let mut row_ids = Vec::new();
    let mut values = Vec::new();
    for record in reader.records() {
        let record = record?;
        if record.is_empty() {
            continue;
        }
        row_ids.push(record.get(0).unwrap_or_default().to_string());
        let row_values = record
            .iter()
            .skip(1)
            .map(|field| {
                field
                    .parse::<f64>()
                    .with_context(|| format!("failed to parse float '{field}'"))
            })
            .collect::<Result<Vec<_>>>()?;
        values.push(row_values);
    }
    let matrix = MatrixData {
        row_ids,
        col_ids,
        values,
    };
    matrix.validate()?;
    Ok(matrix)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_nonfinite_values_duplicate_ids_and_wrong_dimensions() {
        for input in [
            "gene_id\tS\nG\tNaN\n",
            "gene_id\tS\tS\nG\t1\t2\n",
            "gene_id\tS\nG\t1\nG\t2\n",
            "gene_id\tS\nG\t1\t2\n",
        ] {
            assert!(read_matrix_from_reader(input.as_bytes()).is_err());
        }
        let matrix = MatrixData {
            row_ids: vec!["G".into()],
            col_ids: vec!["S".into()],
            values: vec![],
        };
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("matrix.tsv");
        assert!(write_matrix_tsv(&path, &matrix).is_err());
        assert!(!path.exists());
    }
}
