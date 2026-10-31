//! Conservative payload estimates checked before loading the numerical matrix.
use anyhow::{Result, bail};
use rnaa_core::model::{CorrelationMethod, OutputMode};
use serde::Serialize;
use std::path::Path;

#[derive(Clone, Copy, Debug, Serialize)]
pub struct CorrelationLimits {
    pub block_rows: usize,
    pub max_memory_mb: usize,
    pub max_output_mb: usize,
}
impl Default for CorrelationLimits {
    fn default() -> Self {
        Self::from(&rnaa_core::config::CorrConfig::default())
    }
}
impl From<&rnaa_core::config::CorrConfig> for CorrelationLimits {
    fn from(config: &rnaa_core::config::CorrConfig) -> Self {
        Self {
            block_rows: config.block_rows,
            max_memory_mb: config.max_memory_mb,
            max_output_mb: config.max_output_mb,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct CorrelationEstimate {
    pub genes: usize,
    pub samples: usize,
    pub input_bytes: u64,
    pub working_bytes: u128,
    pub output_bytes: u128,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct AdjustmentEstimate {
    pub genes: usize,
    pub samples: usize,
    pub input_bytes: u64,
    pub gene_id_bytes: u64,
    pub memory_bytes: u64,
    pub output_bytes: u64,
}

/// Estimate row-wise adjustment and its complete adjusted matrix output.
pub fn estimate_adjustment(
    path: &Path,
    sheet: &Path,
    residualize: bool,
) -> Result<AdjustmentEstimate> {
    let shape = inspect_shape(path)?;
    let n = shape.genes as u128;
    let s = shape.samples as u128;
    let add = |a: u128, b: u128| {
        a.checked_add(b)
            .ok_or_else(|| anyhow::anyhow!("adjustment estimate overflow"))
    };
    let mul = |a: u128, b: u128| {
        a.checked_mul(b)
            .ok_or_else(|| anyhow::anyhow!("adjustment estimate overflow"))
    };
    let input_bytes = std::fs::metadata(path)?.len();
    let mut memory = 64_u128 * 1024 * 1024;
    for term in [
        mul(shape.id_bytes, 8)?,
        mul(n, 512)?,
        mul(shape.max_record_bytes, 4)?,
        mul(shape.header_bytes, 4)?,
        mul(s, 1024)?,
        mul(u128::from(input_bytes), 4)?,
    ] {
        memory = add(memory, term)?;
    }
    if residualize {
        let metadata_bytes = std::fs::metadata(sheet)?.len() as u128;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(sheet)?;
        let header = reader.headers()?.clone();
        let header_bytes = header.as_slice().len() as u128;
        let header_fields = header.len() as u128;
        let mut rows = 0_u128;
        for record in reader.records() {
            record?;
            rows = add(rows, 1)?;
        }
        for term in [
            mul(metadata_bytes, 32)?,
            mul(rows, add(mul(header_bytes, 2)?, mul(header_fields, 128)?)?)?,
            mul(s, mul(s, 256)?)?,
        ] {
            memory = add(memory, term)?;
        }
    }
    let mut output = mul(mul(n, s)?, 330)?;
    for term in [
        mul(shape.id_bytes, 2)?,
        mul(shape.header_bytes, 2)?,
        mul(n, 4)?,
        mul(s, 4)?,
        1024,
    ] {
        output = add(output, term)?;
    }
    Ok(AdjustmentEstimate {
        genes: shape.genes,
        samples: shape.samples,
        input_bytes,
        gene_id_bytes: u64::try_from(shape.id_bytes).map_err(|_| {
            anyhow::anyhow!("adjustment gene identifiers exceed resource accounting")
        })?,
        memory_bytes: u64::try_from(memory).map_err(|_| {
            anyhow::anyhow!("adjustment memory estimate exceeds resource accounting")
        })?,
        output_bytes: u64::try_from(output).map_err(|_| {
            anyhow::anyhow!("adjustment output estimate exceeds resource accounting")
        })?,
    })
}
impl CorrelationLimits {
    pub fn preflight(
        &self,
        path: &Path,
        method: CorrelationMethod,
        mode: &OutputMode,
        threads: usize,
    ) -> Result<CorrelationEstimate> {
        if self.block_rows == 0
            || self.block_rows > 1024
            || self.max_memory_mb == 0
            || self.max_output_mb == 0
        {
            bail!("correlation requires block_rows in 1..=1024 and positive memory/output limits");
        }
        match mode {
            OutputMode::TopK { k: 0 } => bail!("top-k correlation requires k greater than zero"),
            OutputMode::Threshold { min_abs_r }
                if !min_abs_r.is_finite() || !(0.0..=1.0).contains(min_abs_r) =>
            {
                bail!("correlation threshold must be finite and within [0, 1]")
            }
            _ => {}
        }
        let Shape {
            genes,
            samples,
            max_id,
            id_bytes,
            header_bytes,
            max_record_bytes,
        } = inspect_shape(path)?;
        let n = genes as u128;
        let s = samples as u128;
        let block = self.block_rows.min(genes) as u128;
        let k = match mode {
            OutputMode::TopK { k } => (*k).min(genes - 1) as u128,
            _ => 0,
        };
        // Include input/flatten overlap, IDs, containers, tile, row preparation and
        // growable stripe/heap capacities. HellCor has (20+1)^2 scratch per sample.
        let engine_scratch = if method == CorrelationMethod::Hellcor {
            (block * 2 * 32 * s + threads.max(1) as u128 * (448 * s + 1024)) * 8
        } else {
            (block * 2 * s + threads.max(1) as u128 * s * 32) * 8
        };
        let stripe = if matches!(mode, OutputMode::TopK { .. }) {
            block * k * 32
        } else {
            block * n * 32
        };
        let working_bytes = n * s * 32
            + n * 256
            + id_bytes * 3
            + header_bytes * 3
            + max_record_bytes * 4
            + block * block * 16
            + engine_scratch
            + stripe
            + s * 128;
        let edges = match mode {
            OutputMode::Dense => n * n,
            OutputMode::TopK { .. } => n * k,
            OutputMode::Threshold { .. } => n * (n - 1) / 2,
        };
        // f64 Display may emit a subnormal using over 320 decimal places. IDs may
        // need TSV quoting. Use a conservative upper bound, not an average value.
        let output_bytes = edges * (330 + max_id as u128 * 4 + 8) + id_bytes * 2 + n * 8 + 1024;
        if working_bytes > self.max_memory_mb as u128 * 1048576 {
            bail!(
                "correlation working estimate {working_bytes} bytes exceeds max_memory_mb={}; reduce block_rows/threads/genes or raise the limit",
                self.max_memory_mb
            );
        }
        if output_bytes > self.max_output_mb as u128 * 1048576 {
            bail!(
                "correlation output estimate {output_bytes} bytes exceeds max_output_mb={}; select fewer edges/genes or raise the limit",
                self.max_output_mb
            );
        }
        Ok(CorrelationEstimate {
            genes,
            samples,
            input_bytes: std::fs::metadata(path)?.len(),
            working_bytes,
            output_bytes,
        })
    }
}

#[derive(Clone, Copy)]
pub(super) struct Shape {
    pub genes: usize,
    pub samples: usize,
    pub max_id: usize,
    pub id_bytes: u128,
    pub header_bytes: u128,
    pub max_record_bytes: u128,
}
pub(super) fn inspect_shape(path: &Path) -> Result<Shape> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let samples = reader.headers()?.len().saturating_sub(1);
    let header_bytes = reader.headers()?.as_slice().len() as u128;
    let mut genes = 0usize;
    let mut max_id = 0usize;
    let mut id_bytes = 0u128;
    let mut max_record_bytes = header_bytes;
    for record in reader.records() {
        let record = record?;
        max_record_bytes = max_record_bytes.max(record.as_slice().len() as u128);
        let id = record.get(0).unwrap_or_default().len();
        max_id = max_id.max(id);
        id_bytes += id as u128;
        genes += 1;
    }
    if genes < 2 || samples < 2 {
        bail!("correlation requires at least two genes and two samples");
    }
    Ok(Shape {
        genes,
        samples,
        max_id,
        id_bytes,
        header_bytes,
        max_record_bytes,
    })
}
