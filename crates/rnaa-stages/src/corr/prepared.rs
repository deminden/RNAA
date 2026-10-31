//! Reusable row preparation and bounded calls into the pinned correlation engine.
use mincorr::{hellcor, kendall, rank_data};
use ndarray::{Array1, Array2, s};
use rayon::prelude::*;
use rnaa_core::model::CorrelationMethod;
use std::ops::Range;

pub(super) struct Prepared {
    values: Array2<f64>,
    method: CorrelationMethod,
    pub constant: Vec<bool>,
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_unstable_by(f64::total_cmp);
    let mid = values.len() / 2;
    if values.len().is_multiple_of(2) {
        values[mid - 1] * 0.5 + values[mid] * 0.5
    } else {
        values[mid]
    }
}

fn unit_length(row: &mut [f64]) {
    let norm = row.iter().map(|v| v * v).sum::<f64>().sqrt();
    if norm > 0.0 && norm.is_finite() {
        for value in row {
            *value /= norm;
        }
    } else {
        row.fill(f64::NAN);
    }
}

fn prepare_dot(row: &mut [f64], method: CorrelationMethod) {
    if method == CorrelationMethod::Spearman {
        // Rank original values: scaling first can merge adjacent floating values.
        let ranks = rank_data(&Array1::from_vec(row.to_vec()));
        row.copy_from_slice(ranks.as_slice().expect("contiguous ranks"));
    }
    // Scaling first avoids overflow while preserving these affine-invariant methods.
    let scale = row.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
    if scale > 0.0 {
        for value in row.iter_mut() {
            *value /= scale;
        }
    }
    if method == CorrelationMethod::Bicor {
        // WGCNA defaults: maxPOutliers=1, cosine=FALSE, individual Pearson fallback.
        // With maxPOutliers=1 the extra outlier-quantile rescaling is the identity.
        let center = median(&mut row.to_vec());
        let mad = median(
            &mut row
                .iter()
                .map(|value| (value - center).abs())
                .collect::<Vec<_>>(),
        );
        if mad > 0.0 {
            for value in row.iter_mut() {
                let deviation = *value - center;
                let u = (deviation / (9.0 * mad)).abs().min(1.0);
                *value = deviation * (1.0 - u * u).powi(2);
            }
            unit_length(row);
            return;
        }
    }
    let center = row.iter().sum::<f64>() / row.len() as f64;
    for value in row.iter_mut() {
        *value -= center;
    }
    unit_length(row);
}

impl Prepared {
    pub fn correlations_with(&self, row: &[f64]) -> Vec<f64> {
        assert_eq!(row.len(), self.values.ncols());
        let mut prepared = row.to_vec();
        prepare_dot(&mut prepared, self.method);
        let prepared = Array1::from_vec(prepared);
        self.values
            .rows()
            .into_iter()
            .enumerate()
            .map(|(i, row)| {
                if self.constant[i] {
                    f64::NAN
                } else {
                    canonical(row.dot(&prepared))
                }
            })
            .collect()
    }
    pub fn new(rows: Vec<Vec<f64>>, method: CorrelationMethod) -> Self {
        let shape = (rows.len(), rows[0].len());
        let constant = rows
            .iter()
            .map(|row| row.iter().all(|v| *v == row[0]))
            .collect();
        let mut values = Array2::from_shape_vec(shape, rows.into_iter().flatten().collect())
            .expect("validated matrix shape");
        if matches!(
            method,
            CorrelationMethod::Pearson | CorrelationMethod::Spearman | CorrelationMethod::Bicor
        ) {
            for mut row in values.rows_mut() {
                prepare_dot(row.as_slice_mut().expect("contiguous row"), method);
            }
        }
        Self {
            values,
            method,
            constant,
        }
    }

    pub fn block(&self, rows: Range<usize>, cols: Range<usize>) -> Array2<f64> {
        let left = self.values.slice(s![rows.clone(), ..]);
        let right = self.values.slice(s![cols.clone(), ..]);
        let mut result = match self.method {
            CorrelationMethod::Kendall => kendall::correlation_cross_matrix(&left, &right),
            CorrelationMethod::Hellcor => hellcor::correlation_cross_matrix(&left, &right),
            _ => {
                let mut result = Array2::zeros((rows.len(), cols.len()));
                result
                    .as_slice_mut()
                    .expect("contiguous tile")
                    .par_chunks_mut(cols.len())
                    .enumerate()
                    .for_each(|(i, output)| {
                        for (j, value) in output.iter_mut().enumerate() {
                            // A fixed-length vector dot keeps rounding independent of
                            // tile dimensions, preserving edge ties across allocations.
                            *value = left.row(i).dot(&right.row(j));
                        }
                    });
                result
            }
        };
        for ((i, j), value) in result.indexed_iter_mut() {
            let (gene, target) = (rows.start + i, cols.start + j);
            if self.constant[gene] || self.constant[target] || !value.is_finite() {
                *value = f64::NAN;
            } else if gene == target {
                *value = 1.0;
            } else {
                *value = canonical(*value);
            }
        }
        result
    }
}

fn canonical(value: f64) -> f64 {
    if !value.is_finite() {
        return f64::NAN;
    }
    let value = value.clamp(-1.0, 1.0);
    if 1.0 - value.abs() <= 8.0 * f64::EPSILON {
        value.signum()
    } else {
        value
    }
}
