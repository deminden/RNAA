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

pub fn read_matrix_tsv(path: &Path) -> Result<MatrixData> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    read_matrix_from_reader(file)
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
    Ok(MatrixData {
        row_ids,
        col_ids,
        values,
    })
}
