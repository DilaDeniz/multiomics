use anyhow::Result;
use ndarray::{Array1, Array2};

/// Result of dimensionality reduction.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PcaResult {
    /// Projected coordinates, shape (n_samples, n_components).
    pub points: Vec<[f64; 2]>,
    /// Fraction of variance explained by each component.
    pub explained_variance_ratio: Vec<f64>,
}

/// Project an N×M feature matrix to 2D using PCA.
///
/// Performs zero-mean unit-variance standardization before projection.
/// Uses `linfa-reduction` for the actual PCA computation.
///
/// Returns the 2D projection and explained variance ratios, or a fallback
/// identity projection when the matrix is too small or linfa fails.
pub fn run_pca(features: &Array2<f64>, n_components: usize) -> Result<PcaResult> {
    let n_samples = features.nrows();
    let n_features = features.ncols();

    if n_samples < 2 || n_features < 2 {
        // Not enough data for PCA — return identity projection
        let points = (0..n_samples)
            .map(|i| [i as f64, 0.0])
            .collect();
        return Ok(PcaResult {
            points,
            explained_variance_ratio: vec![1.0, 0.0],
        });
    }

    // Standardize: subtract mean, divide by std
    let standardized = standardize(features);

    // Try linfa PCA
    match try_linfa_pca(&standardized, n_components) {
        Ok(result) => Ok(result),
        Err(e) => {
            log::warn!("linfa PCA failed ({}), using manual SVD approximation", e);
            Ok(manual_pca(&standardized, n_components))
        }
    }
}

fn standardize(data: &Array2<f64>) -> Array2<f64> {
    let mut out = data.to_owned();
    let n = data.nrows() as f64;
    for mut col in out.columns_mut() {
        let mean = col.sum() / n;
        let std = col.iter().map(|x| (x - mean).powi(2)).sum::<f64>().sqrt() / n.sqrt();
        if std > 1e-12 {
            col.iter_mut().for_each(|x| *x = (*x - mean) / std);
        } else {
            col.iter_mut().for_each(|x| *x -= mean);
        }
    }
    out
}

fn try_linfa_pca(data: &Array2<f64>, n_components: usize) -> Result<PcaResult> {
    use linfa::prelude::*;
    use linfa_reduction::Pca;

    let n = data.nrows();
    let targets = Array1::<usize>::zeros(n);
    let dataset = linfa::Dataset::new(data.to_owned(), targets);

    let components = n_components.min(data.nrows() - 1).min(data.ncols());
    if components == 0 {
        anyhow::bail!("Cannot run PCA with 0 components");
    }

    let model = Pca::params(components)
        .fit(&dataset)
        .map_err(|e| anyhow::anyhow!("PCA fit failed: {:?}", e))?;

    // In linfa 0.7, predict returns Array2<F> directly
    let projected: ndarray::Array2<f64> = model.predict(&dataset);

    let points: Vec<[f64; 2]> = (0..n)
        .map(|i| {
            let x = projected.get((i, 0)).copied().unwrap_or(0.0);
            let y = projected.get((i, 1)).copied().unwrap_or(0.0);
            [x, y]
        })
        .collect();

    let explained: Vec<f64> = model
        .explained_variance_ratio()
        .iter()
        .copied()
        .collect();

    Ok(PcaResult {
        points,
        explained_variance_ratio: explained,
    })
}

/// Simple 2-component PCA via covariance matrix eigenvectors (power iteration).
/// Used as fallback when linfa is unavailable or fails.
fn manual_pca(data: &Array2<f64>, _n_components: usize) -> PcaResult {
    let n = data.nrows();
    let d = data.ncols();

    // Covariance matrix C = (1/n) * X^T * X
    let mut cov = Array2::<f64>::zeros((d, d));
    for i in 0..d {
        for j in 0..d {
            let col_i = data.column(i);
            let col_j = data.column(j);
            cov[[i, j]] = col_i.iter().zip(col_j.iter()).map(|(a, b)| a * b).sum::<f64>()
                / n as f64;
        }
    }

    // Power iteration for first two eigenvectors
    let pc1 = power_iteration(&cov, 100);
    let pc2 = deflate_and_iterate(&cov, &pc1, 100);

    // Project data onto the two PCs
    let points: Vec<[f64; 2]> = (0..n)
        .map(|i| {
            let row = data.row(i);
            let x: f64 = row.iter().zip(pc1.iter()).map(|(a, b)| a * b).sum();
            let y: f64 = row.iter().zip(pc2.iter()).map(|(a, b)| a * b).sum();
            [x, y]
        })
        .collect();

    PcaResult {
        points,
        explained_variance_ratio: vec![0.0, 0.0],
    }
}

fn power_iteration(cov: &Array2<f64>, iterations: usize) -> Vec<f64> {
    let d = cov.nrows();
    let mut v: Vec<f64> = (0..d).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
    for _ in 0..iterations {
        let mut w = vec![0.0f64; d];
        for i in 0..d {
            for j in 0..d {
                w[i] += cov[[i, j]] * v[j];
            }
        }
        let norm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-12);
        v = w.iter().map(|x| x / norm).collect();
    }
    v
}

fn deflate_and_iterate(cov: &Array2<f64>, pc1: &[f64], iterations: usize) -> Vec<f64> {
    let d = cov.nrows();
    // Deflate: C' = C - λ₁ * v₁ * v₁ᵀ
    let eigenvalue: f64 = {
        let mut lam = 0.0f64;
        for i in 0..d {
            let mut row_sum = 0.0f64;
            for j in 0..d {
                row_sum += cov[[i, j]] * pc1[j];
            }
            lam += pc1[i] * row_sum;
        }
        lam
    };

    let mut deflated = cov.to_owned();
    for i in 0..d {
        for j in 0..d {
            deflated[[i, j]] -= eigenvalue * pc1[i] * pc1[j];
        }
    }

    power_iteration(&deflated, iterations)
}
