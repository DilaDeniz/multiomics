//! Highly variable gene (HVG) selection.
//!
//! Implements the Seurat v3 method (Stuart et al. 2019, Cell 177:1888):
//! bin genes by mean expression, standardize variance within bins,
//! and rank by standardized dispersion.

use ndarray::Array2;

const N_BINS: usize = 20;

/// Select the top `n_top` highly variable genes using the Seurat v3 method.
///
/// Returns the gene indices (column indices into `norm_matrix`) sorted by
/// descending standardized variance.
///
/// `norm_matrix`: `[n_cells × n_genes]`.
pub fn select_hvg(norm_matrix: &Array2<f32>, n_top: usize) -> Vec<usize> {
    let (means, variances) = gene_stats(norm_matrix);
    let n_genes = means.len();
    if n_genes == 0 {
        return Vec::new();
    }

    // Bin genes by log10(mean + 1)
    let log_means: Vec<f64> = means.iter().map(|&m| (m as f64 + 1.0).log10()).collect();
    let min_lm = log_means.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_lm = log_means.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let bin_width = if (max_lm - min_lm).abs() < 1e-12 {
        1.0
    } else {
        (max_lm - min_lm) / N_BINS as f64
    };

    let bin_of: Vec<usize> = log_means
        .iter()
        .map(|&lm| {
            let b = ((lm - min_lm) / bin_width) as usize;
            b.min(N_BINS - 1)
        })
        .collect();

    // For each bin, collect log(variance) values; compute mean and std
    let log_vars: Vec<f64> = variances
        .iter()
        .map(|&v| (v as f64 + 1e-10).ln())
        .collect();

    let mut bin_log_vars: Vec<Vec<f64>> = vec![Vec::new(); N_BINS];
    for (g, &b) in bin_of.iter().enumerate() {
        bin_log_vars[b].push(log_vars[g]);
    }

    let bin_mean: Vec<f64> = bin_log_vars
        .iter()
        .map(|vs| {
            if vs.is_empty() {
                0.0
            } else {
                vs.iter().sum::<f64>() / vs.len() as f64
            }
        })
        .collect();

    let bin_std: Vec<f64> = bin_log_vars
        .iter()
        .zip(bin_mean.iter())
        .map(|(vs, &m)| {
            if vs.len() < 2 {
                1.0
            } else {
                let var = vs.iter().map(|&v| (v - m).powi(2)).sum::<f64>() / (vs.len() - 1) as f64;
                var.sqrt().max(1e-10)
            }
        })
        .collect();

    // Compute standardized variance per gene, clip at 10
    let mut scores: Vec<(usize, f64)> = (0..n_genes)
        .map(|g| {
            let b = bin_of[g];
            let z = (log_vars[g] - bin_mean[b]) / bin_std[b];
            (g, z.min(10.0))
        })
        .collect();

    scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    scores.truncate(n_top);
    scores.into_iter().map(|(g, _)| g).collect()
}

/// Compute per-gene (column) means and variances using Welford's online algorithm.
///
/// `matrix`: `[n_cells × n_genes]`.
pub fn gene_stats(matrix: &Array2<f32>) -> (Vec<f32>, Vec<f32>) {
    let (n_cells, n_genes) = matrix.dim();
    let mut means = vec![0.0f64; n_genes];
    let mut m2 = vec![0.0f64; n_genes];

    for i in 0..n_cells {
        for g in 0..n_genes {
            let x = matrix[[i, g]] as f64;
            let delta = x - means[g];
            means[g] += delta / (i + 1) as f64;
            let delta2 = x - means[g];
            m2[g] += delta * delta2;
        }
    }

    let variances: Vec<f32> = m2
        .iter()
        .map(|&s| if n_cells > 1 { (s / (n_cells - 1) as f64) as f32 } else { 0.0 })
        .collect();

    let means_f32: Vec<f32> = means.iter().map(|&m| m as f32).collect();
    (means_f32, variances)
}
