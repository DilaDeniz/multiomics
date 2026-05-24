//! Multi-Omics Factor Analysis Plus (MOFA+) — variational Bayes implementation.
//!
//! MOFA+ learns K latent factors that explain co-variation across M modalities
//! observed in the same N samples. Each factor can be active in any subset of
//! modalities (automatic relevance determination), revealing shared and
//! modality-specific biological signals.
//!
//! # Algorithm
//! Variational EM with automatic relevance determination (ARD) sparsity priors
//! on the factor loadings. The generative model is:
//!
//! ```text
//! X_m  ≈  Z · W_m^T  +  ε_m        (X_m: n × d_m, Z: n × K, W_m: d_m × K)
//! w_mk  ~  N(0, 1/α_mk)             (ARD prior — drives inactive factors to 0)
//! ε_m   ~  N(0, τ_m^{-1} I)         (modality-specific noise precision)
//! ```
//!
//! The VB updates iterate until ELBO convergence (relative change < ε).
//!
//! # References
//! * Argelaguet R et al. (2020) MOFA+: a statistical framework for
//!   comprehensive integration of multi-modal single-cell data.
//!   Genome Biol. 21:111. <https://doi.org/10.1186/s13059-020-02015-1>
//! * Argelaguet R et al. (2018) Multi-Omics Factor Analysis — a framework for
//!   unsupervised integration of multi-omics data sets. Mol. Syst. Biol. 14:e8124.

use ndarray::{Array1, Array2};

/// Configuration for the MOFA+ algorithm.
#[derive(Debug, Clone)]
pub struct MofaConfig {
    /// Number of latent factors to learn.
    pub n_factors: usize,
    /// Maximum VB iterations.
    pub max_iter: usize,
    /// Convergence threshold on relative ELBO change.
    pub tol: f64,
    /// Minimum ARD precision — factors with α > ard_threshold are considered inactive.
    pub ard_threshold: f64,
    /// Random seed for initialisation.
    pub seed: u64,
}

impl Default for MofaConfig {
    fn default() -> Self {
        Self {
            n_factors: 10,
            max_iter: 1000,
            tol: 1e-6,
            ard_threshold: 1e3,
            seed: 0xDEAD_BEEF_C0FF_EE42,
        }
    }
}

/// Factor scores and per-modality loadings returned by MOFA+.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MofaResult {
    /// Factor scores for each sample: shape `(n_samples, n_factors)`.
    pub scores: Vec<Vec<f64>>,
    /// Per-modality factor loadings: `loadings[m]` has shape `(d_m, n_factors)`.
    pub loadings: Vec<Vec<Vec<f64>>>,
    /// ARD precision per (modality, factor): large values → factor inactive in modality.
    pub ard: Vec<Vec<f64>>,
    /// R² explained by each factor in each modality: `r2[m][k]`.
    pub r2_per_modality_factor: Vec<Vec<f64>>,
    /// Active factors per modality (ard < ard_threshold).
    pub active_factors: Vec<Vec<bool>>,
    /// Number of VB iterations run.
    pub n_iter: usize,
    /// Final ELBO value.
    pub elbo: f64,
}

/// Run MOFA+ on a set of modality matrices.
///
/// # Arguments
/// - `views`: list of `(name, data)` pairs. Each `data` is an `n_samples × d_m`
///   matrix. All views must have the same number of rows (samples).
/// - `cfg`: algorithm configuration.
///
/// # Returns
/// `MofaResult` containing factor scores, per-modality loadings, and
/// R² explained per factor per modality.
pub fn run_mofa(views: &[(&str, &Array2<f64>)], cfg: &MofaConfig) -> anyhow::Result<MofaResult> {
    let n_views = views.len();
    if n_views == 0 {
        anyhow::bail!("run_mofa: no views provided");
    }

    let n_samples = views[0].1.nrows();
    for (name, mat) in views {
        if mat.nrows() != n_samples {
            anyhow::bail!(
                "run_mofa: view '{}' has {} rows but expected {}",
                name,
                mat.nrows(),
                n_samples
            );
        }
    }

    let k = cfg.n_factors.min(n_samples - 1).max(1);
    let dims: Vec<usize> = views.iter().map(|(_, m)| m.ncols()).collect();

    // Standardize each view to zero mean, unit variance (per feature).
    let views_std: Vec<Array2<f64>> = views.iter().map(|(_, m)| standardize_cols(m)).collect();

    // Initialise factor scores Z with seeded pseudo-random normals.
    let mut z = init_z(n_samples, k, cfg.seed);

    // Initialise loadings W_m (d_m × k) for each view.
    let mut w: Vec<Array2<f64>> = dims
        .iter()
        .enumerate()
        .map(|(m, &d)| init_w(d, k, cfg.seed ^ (m as u64 * 0x9E3779B97F4A7C15)))
        .collect();

    // ARD precisions α[m][k]: one per (view, factor).
    let mut alpha: Vec<Array1<f64>> = dims.iter().map(|&d| Array1::ones(k) * (d as f64)).collect();

    // Noise precisions τ[m]: one per view (initialise to 1).
    let mut tau: Vec<f64> = vec![1.0; n_views];

    let mut prev_elbo = f64::NEG_INFINITY;

    let mut n_iter = 0usize;
    for iter in 0..cfg.max_iter {
        n_iter = iter + 1;

        // ── VB E-step: update factor scores Z ─────────────────────────────
        // Posterior: q(Z) = N(μ_Z, Σ_Z)
        // Σ_Z = (I + Σ_m τ_m * W_m^T W_m)^{-1}   (k × k)
        // μ_Z = Σ_Z * Σ_m τ_m * X_m * W_m          (n × k)
        let mut sigma_z_inv = Array2::<f64>::eye(k); // precision (k×k)
        let mut rhs = Array2::<f64>::zeros((n_samples, k));

        for m in 0..n_views {
            let wm = &w[m]; // d_m × k
            let tm = tau[m];
            // WtW (k×k)
            let wtw = wm.t().dot(wm);
            sigma_z_inv = sigma_z_inv + wtw * tm;
            // X_m * W_m  (n × k)
            let xw = views_std[m].dot(wm);
            rhs = rhs + xw * tm;
        }

        // Solve Σ_Z = sigma_z_inv^{-1} by Cholesky (k is small, ≤ 50).
        let sigma_z = chol_inv_k(&sigma_z_inv, k);
        // μ_Z = rhs * Σ_Z^T  (n × k)
        z = rhs.dot(&sigma_z.t());

        // ── VB M-step: update loadings W_m ─────────────────────────────────
        // Σ_W_m = (diag(α_m) + τ_m * Z^T Z)^{-1}   (k × k)
        // μ_W_m = X_m^T * Z * Σ_W_m * τ_m           (d_m × k)
        let ztz = z.t().dot(&z); // k × k
        for m in 0..n_views {
            let tm = tau[m];
            let mut sigma_w_inv = Array2::<f64>::zeros((k, k));
            for (fac, &a) in alpha[m].iter().enumerate() {
                sigma_w_inv[[fac, fac]] = a;
            }
            sigma_w_inv = sigma_w_inv + ztz.clone() * tm;
            let sigma_w = chol_inv_k(&sigma_w_inv, k);
            // X_m^T Z (d_m × k)
            let xtz = views_std[m].t().dot(&z);
            w[m] = xtz.dot(&sigma_w) * tm;
        }

        // ── Update ARD precisions α_mk ──────────────────────────────────────
        // α_mk = d_m / Σ_d w_dm^2   (column 2-norm of W_m)
        for m in 0..n_views {
            let d = dims[m] as f64;
            for (fac, a) in alpha[m].iter_mut().enumerate() {
                let col_norm_sq: f64 = w[m].column(fac).iter().map(|&x| x * x).sum();
                *a = d / (col_norm_sq + 1e-10);
            }
        }

        // ── Update noise precisions τ_m ──────────────────────────────────────
        // τ_m = (n * d_m) / ||X_m - Z * W_m^T||^2_F
        for m in 0..n_views {
            let xhat = z.dot(&w[m].t()); // n × d_m
            let diff = &views_std[m] - &xhat;
            let ss: f64 = diff.iter().map(|&x| x * x).sum();
            let n_d = (n_samples * dims[m]) as f64;
            tau[m] = n_d / (ss + 1e-10);
        }

        // ── ELBO (simplified reconstruction term) ──────────────────────────
        let elbo = compute_elbo(&views_std, &z, &w, &alpha, &tau, n_samples, &dims);

        let rel_change = if prev_elbo.is_finite() {
            (elbo - prev_elbo).abs() / (prev_elbo.abs() + 1.0)
        } else {
            f64::INFINITY
        };

        prev_elbo = elbo;

        if rel_change < cfg.tol {
            log::debug!("MOFA+ converged at iteration {iter}, ELBO={elbo:.4}");
            break;
        }
    }

    // ── Compute R² per (modality, factor) ─────────────────────────────────
    // R²_mk = 1 - ||X_m - z_k * w_mk^T||^2 / ||X_m||^2
    let r2_per_modality_factor = compute_r2_per_factor(&views_std, &z, &w, &dims, k);

    // ── Determine active factors ───────────────────────────────────────────
    let active_factors: Vec<Vec<bool>> = (0..n_views)
        .map(|m| (0..k).map(|ki| alpha[m][ki] < cfg.ard_threshold).collect())
        .collect();

    // Convert ndarray structures to plain vecs for serialisation.
    let scores: Vec<Vec<f64>> = (0..n_samples)
        .map(|i| z.row(i).iter().copied().collect())
        .collect();

    let loadings: Vec<Vec<Vec<f64>>> = w
        .iter()
        .map(|wm| {
            (0..wm.nrows())
                .map(|d| wm.row(d).iter().copied().collect())
                .collect()
        })
        .collect();

    let ard: Vec<Vec<f64>> = alpha.iter().map(|a| a.iter().copied().collect()).collect();

    Ok(MofaResult {
        scores,
        loadings,
        ard,
        r2_per_modality_factor,
        active_factors,
        n_iter,
        elbo: prev_elbo,
    })
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Standardize columns of a matrix to zero mean, unit variance.
fn standardize_cols(mat: &Array2<f64>) -> Array2<f64> {
    let mut out = mat.to_owned();
    let n = mat.nrows() as f64;
    for mut col in out.columns_mut() {
        let mean = col.iter().sum::<f64>() / n;
        let var = col.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
        let std = var.sqrt();
        if std > 1e-12 {
            col.iter_mut().for_each(|x| *x = (*x - mean) / std);
        } else {
            col.iter_mut().for_each(|x| *x -= mean);
        }
    }
    out
}

/// Initialise a matrix with small pseudo-random values (LCG, no external deps).
fn init_z(n: usize, k: usize, seed: u64) -> Array2<f64> {
    let mut state = seed;
    let mut data = vec![0.0f64; n * k];
    for x in &mut data {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        // Map high bits to [-0.1, 0.1]
        *x = ((state >> 33) as f64 / (u32::MAX as f64) - 0.5) * 0.2;
    }
    Array2::from_shape_vec((n, k), data).unwrap()
}

fn init_w(d: usize, k: usize, seed: u64) -> Array2<f64> {
    init_z(d, k, seed)
}

/// Invert a small k×k positive-definite matrix via Cholesky decomposition.
/// Falls back to regularised inverse if Cholesky fails.
fn chol_inv_k(a: &Array2<f64>, k: usize) -> Array2<f64> {
    // Attempt Cholesky: L L^T = A, then solve L L^T X = I column by column.
    let l = cholesky(a, k);
    if l.is_none() {
        // Regularise and retry once.
        let mut reg = a.to_owned();
        for i in 0..k {
            reg[[i, i]] += 1e-6;
        }
        let l2 = cholesky(&reg, k);
        return l2.map(|l| chol_solve_identity(&l, k)).unwrap_or_else(|| {
            // Last resort: return identity (very ill-conditioned system).
            Array2::eye(k)
        });
    }
    chol_solve_identity(&l.unwrap(), k)
}

/// Lower-triangular Cholesky factor of a k×k PD matrix.
fn cholesky(a: &Array2<f64>, k: usize) -> Option<Array2<f64>> {
    let mut l = Array2::<f64>::zeros((k, k));
    for i in 0..k {
        for j in 0..=i {
            let mut s = a[[i, j]];
            for kk in 0..j {
                s -= l[[i, kk]] * l[[j, kk]];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[[i, j]] = s.sqrt();
            } else {
                l[[i, j]] = s / l[[j, j]];
            }
        }
    }
    Some(l)
}

/// Solve L L^T X = I for X, given lower triangular L.
fn chol_solve_identity(l: &Array2<f64>, k: usize) -> Array2<f64> {
    let mut inv = Array2::<f64>::zeros((k, k));
    for col in 0..k {
        // Forward substitution: L y = e_col
        let mut y = vec![0.0f64; k];
        for i in 0..k {
            let mut s = if i == col { 1.0 } else { 0.0 };
            for j in 0..i {
                s -= l[[i, j]] * y[j];
            }
            y[i] = s / l[[i, i]];
        }
        // Back substitution: L^T x = y
        let mut x = vec![0.0f64; k];
        for i in (0..k).rev() {
            let mut s = y[i];
            for j in (i + 1)..k {
                s -= l[[j, i]] * x[j];
            }
            x[i] = s / l[[i, i]];
        }
        for i in 0..k {
            inv[[i, col]] = x[i];
        }
    }
    inv
}

/// Simplified ELBO: sum of log-likelihood terms (reconstruction + ARD penalty).
fn compute_elbo(
    views: &[Array2<f64>],
    z: &Array2<f64>,
    w: &[Array2<f64>],
    alpha: &[Array1<f64>],
    tau: &[f64],
    n: usize,
    dims: &[usize],
) -> f64 {
    let mut elbo = 0.0;
    for m in 0..views.len() {
        let xhat = z.dot(&w[m].t());
        let diff = &views[m] - &xhat;
        let ss: f64 = diff.iter().map(|&x| x * x).sum();
        let d = dims[m] as f64;
        let tau_m = tau[m];
        // Expected log-likelihood: 0.5 * n*d*ln(τ) - 0.5*τ*||X - ZW^T||^2
        elbo += 0.5 * (n as f64) * d * tau_m.ln() - 0.5 * tau_m * ss;
        // ARD penalty: -0.5 * Σ_k α_mk * ||w_mk||^2
        for (fac, &a) in alpha[m].iter().enumerate() {
            let col_sq: f64 = w[m].column(fac).iter().map(|&x| x * x).sum();
            elbo -= 0.5 * a * col_sq;
        }
    }
    elbo
}

/// R² explained by each individual factor k in modality m.
/// Removes factor k from the reconstruction and measures the drop in SS.
fn compute_r2_per_factor(
    views: &[Array2<f64>],
    z: &Array2<f64>,
    w: &[Array2<f64>],
    dims: &[usize],
    k: usize,
) -> Vec<Vec<f64>> {
    views
        .iter()
        .enumerate()
        .map(|(m, xm)| {
            let total_ss: f64 = xm.iter().map(|&x| x * x).sum();
            if total_ss < 1e-12 {
                return vec![0.0; k];
            }
            let xhat_full = z.dot(&w[m].t());
            let resid_full: f64 = (xm - &xhat_full).iter().map(|&x| x * x).sum();
            let _ = dims[m];

            (0..k)
                .map(|ki| {
                    // Reconstruction without factor ki
                    let mut xhat_no_k = xhat_full.clone();
                    let z_col = z.column(ki);
                    let w_col = w[m].column(ki);
                    // Subtract contribution of factor ki: z_k ⊗ w_mk
                    for i in 0..z.nrows() {
                        for d in 0..w[m].nrows() {
                            xhat_no_k[[i, d]] -= z_col[i] * w_col[d];
                        }
                    }
                    let resid_no_k: f64 = (xm - &xhat_no_k).iter().map(|&x| x * x).sum();
                    // R²_k = (SS_no_k - SS_full) / SS_total
                    ((resid_no_k - resid_full) / total_ss).max(0.0)
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn mofa_basic_convergence() {
        // 5 samples, 2 views of 4 and 3 features, 2 factors
        let x1 = array![
            [1.0, 0.5, -1.0, 0.2],
            [2.0, 1.0, -2.0, 0.4],
            [3.0, 1.5, -3.0, 0.6],
            [-1.0, -0.5, 1.0, -0.2],
            [-2.0, -1.0, 2.0, -0.4],
        ];
        let x2 = array![
            [0.1, 0.9, -0.5],
            [0.2, 1.8, -1.0],
            [0.3, 2.7, -1.5],
            [-0.1, -0.9, 0.5],
            [-0.2, -1.8, 1.0],
        ];

        let cfg = MofaConfig {
            n_factors: 2,
            max_iter: 200,
            tol: 1e-5,
            ..Default::default()
        };
        let views = vec![("genomics", &x1), ("epigenomics", &x2)];
        let result = run_mofa(&views, &cfg).unwrap();

        assert_eq!(result.scores.len(), 5, "5 samples");
        assert_eq!(result.scores[0].len(), 2, "2 factors");
        assert_eq!(result.loadings.len(), 2, "2 modalities");
        assert!(result.elbo.is_finite(), "ELBO must be finite");
    }

    #[test]
    fn mofa_r2_nonnegative() {
        let x1 = array![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [-1.0, -2.0, -3.0],
        ];
        let cfg = MofaConfig {
            n_factors: 2,
            max_iter: 100,
            ..Default::default()
        };
        let views = vec![("rna", &x1)];
        let result = run_mofa(&views, &cfg).unwrap();
        for &r2 in &result.r2_per_modality_factor[0] {
            assert!(r2 >= 0.0, "R² must be non-negative, got {r2}");
        }
    }

    #[test]
    fn mofa_mismatched_rows_errors() {
        let x1 = array![[1.0, 2.0], [3.0, 4.0]];
        let x2 = array![[1.0], [2.0], [3.0]];
        let cfg = MofaConfig::default();
        let views = vec![("a", &x1), ("b", &x2)];
        assert!(run_mofa(&views, &cfg).is_err());
    }
}
