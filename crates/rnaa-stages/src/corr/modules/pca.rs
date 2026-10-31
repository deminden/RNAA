//! Standardized first-component scores with reproducible sign and tied-space policy.
use anyhow::{Result, bail};
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct ModuleScore {
    pub values: Vec<f64>,
    pub variance_explained: f64,
    pub leading_multiplicity: usize,
    pub orientation: &'static str,
}

/// Rows are module genes, columns are samples. Every gene must be variable.
/// Scores have mean zero and sample variance one. A tied leading eigenspace uses
/// the first nonzero projection of a canonical sample axis before sign alignment.
pub fn standardized_module_score(rows: &[Vec<f64>]) -> Result<ModuleScore> {
    let samples = rows.first().map_or(0, Vec::len);
    if rows.len() < 2 || samples < 3 {
        bail!("module PCA requires at least two genes and three samples");
    }
    let mut gram = DMatrix::<f64>::zeros(samples, samples);
    let mut average = DVector::<f64>::zeros(samples);
    for row in rows {
        if row.len() != samples || row.iter().any(|value| !value.is_finite()) {
            bail!("module PCA rows must have equal length and finite values");
        }
        let scale = row.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            bail!("module PCA cannot include a constant gene");
        }
        let mut centered = DVector::from_iterator(samples, row.iter().map(|value| value / scale));
        let mean = centered.iter().sum::<f64>() / samples as f64;
        centered.add_scalar_mut(-mean);
        let norm = centered.norm();
        if norm == 0.0 || !norm.is_finite() {
            bail!("module PCA cannot include a constant or numerically degenerate gene");
        }
        centered *= ((samples - 1) as f64).sqrt() / norm;
        average += &centered;
        gram += &centered * centered.transpose();
    }
    let trace = gram.trace();
    let eig = SymmetricEigen::try_new(gram, 1e-14, 20_000)
        .ok_or_else(|| anyhow::anyhow!("module PCA eigendecomposition did not converge"))?;
    let largest = eig
        .eigenvalues
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    if !largest.is_finite() || largest <= 0.0 {
        bail!("module PCA has no positive finite component");
    }
    let tolerance = 64.0 * f64::EPSILON * samples as f64 * largest;
    let leading = eig
        .eigenvalues
        .iter()
        .enumerate()
        .filter_map(|(i, value)| (largest - *value <= tolerance).then_some(i))
        .collect::<Vec<_>>();
    let mut score = if leading.len() == 1 {
        eig.eigenvectors.column(leading[0]).into_owned()
    } else {
        let mut selected = None;
        for axis in 0..samples {
            let mut projection = DVector::zeros(samples);
            for &column in &leading {
                projection += eig.eigenvectors.column(column) * eig.eigenvectors[(axis, column)];
            }
            if projection.norm() > 64.0 * f64::EPSILON * samples as f64 {
                selected = Some(projection.normalize());
                break;
            }
        }
        selected
            .ok_or_else(|| anyhow::anyhow!("module PCA leading eigenspace has no stable axis"))?
    };
    let alignment = score.dot(&average);
    let orientation = if alignment.abs() > 64.0 * f64::EPSILON * average.norm() * score.norm() {
        if alignment < 0.0 {
            score *= -1.0;
        }
        "along_standardized_average"
    } else {
        let largest_coordinate = score
            .iter()
            .map(|value| value.abs())
            .fold(0.0_f64, f64::max);
        let coordinate_tolerance = largest_coordinate * 64.0 * f64::EPSILON * samples as f64;
        let anchor = score
            .iter()
            .position(|value| value.abs() >= largest_coordinate - coordinate_tolerance)
            .unwrap();
        if score[anchor] < 0.0 {
            score *= -1.0;
        }
        "first_maximal_sample_coordinate_positive"
    };
    let mean = score.iter().sum::<f64>() / samples as f64;
    score.add_scalar_mut(-mean);
    score *= ((samples - 1) as f64).sqrt() / score.norm();
    Ok(ModuleScore {
        values: score.iter().copied().collect(),
        variance_explained: (largest / trace).clamp(0.0, 1.0),
        leading_multiplicity: leading.len(),
        orientation,
    })
}
