//! Dispersion estimation (Love et al. 2014, "Dispersion estimation" section):
//! gene-wise maximum-likelihood estimates, a parametric mean–dispersion trend,
//! an empirical-Bayes prior on log dispersion and maximum a posteriori (MAP)
//! shrinkage towards the trend. Constants follow the DESeq2 defaults.

use rayon::prelude::*;

use biomics_core::statistics::trigamma;

use super::glm::{fit_beta, GlmSettings};
use super::linalg::OlsOperators;
use super::nbinom::{grid_search, DispersionProblem, LineSearch};

/// Lower bound on any dispersion estimate.
pub(crate) const MIN_DISP: f64 = 1e-8;
/// Floor applied to fitted means before dispersion estimation.
const MIN_MU: f64 = 0.5;

/// Upper bound on dispersion estimates: max(10, number of samples).
pub(crate) fn max_disp(n_samples: usize) -> f64 {
    (n_samples as f64).max(10.0)
}

fn line_search() -> LineSearch {
    LineSearch {
        kappa_0: 1.0,
        tol: 1e-6,
        max_iter: 100,
        min_log_alpha: (MIN_DISP / 10.0).ln(),
    }
}

/// Number of distinct rows of the design matrix.
pub(crate) fn n_design_groups(x: &[Vec<f64>]) -> usize {
    let mut rows: Vec<Vec<u64>> = x
        .iter()
        .map(|r| r.iter().map(|v| v.to_bits()).collect())
        .collect();
    rows.sort();
    rows.dedup();
    rows.len()
}

/// Gene-wise estimates and the fitted means used to obtain them.
pub(crate) struct GeneWise {
    pub disp: Vec<f64>,
    /// Fitted means (floored at 0.5), reused for the MAP step.
    pub mu: Vec<Vec<f64>>,
}

/// Per-gene inputs for dispersion estimation (all-zero genes excluded).
pub(crate) struct GeneData<'a> {
    pub counts: &'a [f64],
    pub norm: &'a [f64],
    pub base_mean: f64,
    pub base_var: f64,
}

/// Gene-wise dispersion estimates by Cox-Reid adjusted maximum likelihood.
pub(crate) fn gene_wise(
    genes: &[GeneData<'_>],
    sf: &[f64],
    x: &[Vec<f64>],
    ols: &OlsOperators,
) -> GeneWise {
    let m = sf.len();
    let p = x[0].len();
    let max_disp = max_disp(m);
    // Means from a linear model on normalised counts suffice when every
    // coefficient corresponds to one group of identical design rows.
    let linear_mu = n_design_groups(x) == p;
    let xim = sf.iter().map(|s| 1.0 / s).sum::<f64>() / m as f64;
    let search = line_search();
    let glm = GlmSettings::default();

    let per_gene: Vec<(f64, Vec<f64>)> = genes
        .par_iter()
        .map(|g| {
            // Starting value: the smaller of a rough moments estimate around the
            // linear-model fit and a plain method-of-moments estimate.
            let lm = ols.fitted(g.norm);
            let rough = (g
                .norm
                .iter()
                .zip(&lm)
                .map(|(&yn, &fit)| {
                    let mu = fit.max(1.0);
                    ((yn - mu).powi(2) - mu) / (mu * mu)
                })
                .sum::<f64>()
                / (m - p) as f64)
                .max(0.0);
            let moments = (g.base_var - xim * g.base_mean) / (g.base_mean * g.base_mean);
            let alpha_init = rough.min(moments).clamp(MIN_DISP, max_disp);

            let mu: Vec<f64> = if linear_mu {
                lm.iter()
                    .zip(sf)
                    .map(|(f, s)| (f * s).max(MIN_MU))
                    .collect()
            } else {
                let log_norm: Vec<f64> = g.norm.iter().map(|v| (v + 0.1).ln()).collect();
                let init = ols.coefficients(&log_norm);
                fit_beta(g.counts, sf, x, alpha_init, &init, &glm)
                    .mu
                    .into_iter()
                    .map(|v| v.max(MIN_MU))
                    .collect()
            };

            let problem = DispersionProblem {
                y: g.counts,
                mu: &mu,
                x,
                prior_mean: 0.0,
                prior_var: None,
            };
            let fit = search.run(&problem, alpha_init.ln());
            let mut disp = fit.log_alpha.exp().min(max_disp);
            if fit.last_lp < fit.initial_lp + fit.initial_lp.abs() / 1e6 {
                disp = alpha_init;
            }
            let converged = fit.iterations < search.max_iter && fit.iterations != 1;
            if !converged && disp > MIN_DISP * 10.0 {
                disp = grid_search(&problem, MIN_DISP.ln(), max_disp.ln()).exp();
            }
            (disp.clamp(MIN_DISP, max_disp), mu)
        })
        .collect();

    let (disp, mu) = per_gene.into_iter().unzip();
    GeneWise { disp, mu }
}

/// How the mean–dispersion trend was obtained.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrendKind {
    /// α(μ) = a₀ + a₁/μ fitted by a Gamma-family GLM (DESeq2 `fitType = "parametric"`).
    Parametric,
    /// Constant trend: trimmed mean of gene-wise estimates (DESeq2 `fitType = "mean"`).
    /// Used when the parametric fit fails; R/DESeq2 would use a local regression.
    Mean,
}

/// Fit α(μ) = a₀ + a₁/μ by iterating a Gamma(identity) GLM on the gene-wise
/// estimates, dropping points whose ratio to the current fit is outside
/// (1e-4, 15). Returns `None` if the fit fails or does not converge.
pub(crate) fn parametric_trend(means: &[f64], disps: &[f64]) -> Option<(f64, f64)> {
    let mut coefs = (0.1, 1.0);
    for _ in 0..=10 {
        let good: Vec<usize> = (0..means.len())
            .filter(|&i| {
                let r = disps[i] / (coefs.0 + coefs.1 / means[i]);
                r > 1e-4 && r < 15.0
            })
            .collect();
        let y: Vec<f64> = good.iter().map(|&i| disps[i]).collect();
        let xinv: Vec<f64> = good.iter().map(|&i| 1.0 / means[i]).collect();
        let (new, converged) = gamma_identity_glm(&y, &xinv, coefs)?;
        if !(new.0 > 0.0 && new.1 > 0.0) {
            return None;
        }
        let change = (new.0 / coefs.0).ln().powi(2) + (new.1 / coefs.1).ln().powi(2);
        coefs = new;
        if change < 1e-6 && converged {
            return Some(coefs);
        }
    }
    None
}

/// Gamma-family GLM with identity link, y ~ a + b·x, fitted by IRLS from
/// `start` with step halving whenever a fitted mean leaves (0, ∞).
fn gamma_identity_glm(y: &[f64], x: &[f64], start: (f64, f64)) -> Option<((f64, f64), bool)> {
    if y.len() < 2 {
        return None;
    }
    let deviance = |a: f64, b: f64| -> Option<f64> {
        let mut dev = 0.0;
        for (&yi, &xi) in y.iter().zip(x) {
            let mu = a + b * xi;
            if mu <= 0.0 || !mu.is_finite() {
                return None;
            }
            dev += -2.0 * ((yi / mu).ln() - (yi - mu) / mu);
        }
        Some(dev)
    };
    let mut coef_old = start;
    let mut dev_old = deviance(start.0, start.1)?;
    for _ in 0..25 {
        // Working weights 1/μ², working response y (identity link).
        let (mut s0, mut s1, mut s2, mut t0, mut t1) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for (&yi, &xi) in y.iter().zip(x) {
            let mu = coef_old.0 + coef_old.1 * xi;
            let w = 1.0 / (mu * mu);
            s0 += w;
            s1 += w * xi;
            s2 += w * xi * xi;
            t0 += w * yi;
            t1 += w * xi * yi;
        }
        let det = s0 * s2 - s1 * s1;
        if det.abs() < 1e-300 {
            return None;
        }
        let mut coef = ((s2 * t0 - s1 * t1) / det, (s0 * t1 - s1 * t0) / det);
        let mut dev = deviance(coef.0, coef.1);
        let mut halvings = 0;
        while dev.is_none() {
            halvings += 1;
            if halvings > 25 {
                return None;
            }
            coef = ((coef.0 + coef_old.0) / 2.0, (coef.1 + coef_old.1) / 2.0);
            dev = deviance(coef.0, coef.1);
        }
        let dev = dev?;
        if (dev - dev_old).abs() / (dev.abs() + 0.1) < 1e-8 {
            return Some((coef, true));
        }
        dev_old = dev;
        coef_old = coef;
    }
    Some((coef_old, false))
}

/// R's `mean(x, trim = trim)`.
pub(crate) fn trimmed_mean(values: &[f64], trim: f64) -> f64 {
    let n = values.len();
    if n == 0 {
        return f64::NAN;
    }
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if trim >= 0.5 {
        return median_of_sorted(&v);
    }
    let lo = (n as f64 * trim).floor() as usize; // zero-based start
    let hi = n - lo;
    v[lo..hi].iter().sum::<f64>() / (hi - lo) as f64
}

pub(crate) fn median_of_sorted(v: &[f64]) -> f64 {
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Median absolute deviation with R's consistency constant 1.4826.
fn mad(values: &[f64]) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = median_of_sorted(&v);
    let mut dev: Vec<f64> = v.iter().map(|x| (x - med).abs()).collect();
    dev.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    1.4826 * median_of_sorted(&dev)
}

/// Squared MAD of log residuals from the trend over genes whose gene-wise
/// estimate is at least 100 × the minimum dispersion.
pub(crate) fn var_log_disp_ests(disp_gene: &[f64], disp_fit: &[f64]) -> Option<f64> {
    let resid: Vec<f64> = disp_gene
        .iter()
        .zip(disp_fit)
        .filter(|(&g, _)| g >= MIN_DISP * 100.0)
        .map(|(&g, &f)| g.ln() - f.ln())
        .collect();
    if resid.is_empty() {
        None
    } else {
        Some(mad(&resid).powi(2))
    }
}

/// Variance of the normal prior on log dispersion: the observed spread of log
/// residuals minus the expected sampling variance trigamma((m − p)/2), floored
/// at 0.25.
///
/// For m − p ≤ 3, R/DESeq2 instead matches the residual distribution by Monte
/// Carlo simulation; the closed form is used here for every m > p.
pub(crate) fn prior_variance(var_log: f64, m: usize, p: usize) -> f64 {
    if m > p {
        (var_log - trigamma((m - p) as f64 / 2.0)).max(0.25)
    } else {
        var_log
    }
}

/// MAP dispersion estimates and the dispersion-outlier flags.
pub(crate) struct MapDispersions {
    pub dispersion: Vec<f64>,
    pub outlier: Vec<bool>,
}

/// Shrink gene-wise estimates towards the trend with a normal prior on log α.
/// Genes whose gene-wise estimate lies more than two prior standard deviations
/// above the trend keep their gene-wise estimate (dispersion outliers).
#[allow(clippy::too_many_arguments)]
pub(crate) fn map_dispersions(
    genes: &[GeneData<'_>],
    mu: &[Vec<f64>],
    x: &[Vec<f64>],
    disp_gene: &[f64],
    disp_fit: &[f64],
    prior_var: f64,
    var_log: f64,
    n_samples: usize,
) -> MapDispersions {
    let max_disp = max_disp(n_samples);
    let search = line_search();
    let outlier_sd = 2.0 * var_log.sqrt();

    let (dispersion, outlier) = (0..genes.len())
        .into_par_iter()
        .map(|i| {
            let (dg, df) = (disp_gene[i], disp_fit[i]);
            let init = if dg > 0.1 * df { dg } else { df };
            let problem = DispersionProblem {
                y: genes[i].counts,
                mu: &mu[i],
                x,
                prior_mean: df.ln(),
                prior_var: Some(prior_var),
            };
            let fit = search.run(&problem, init.ln());
            let log_alpha = if fit.iterations < search.max_iter {
                fit.log_alpha
            } else {
                grid_search(&problem, MIN_DISP.ln(), max_disp.ln())
            };
            let map = log_alpha.exp().clamp(MIN_DISP, max_disp);
            let is_outlier = dg.ln() > df.ln() + outlier_sd;
            (if is_outlier { dg } else { map }, is_outlier)
        })
        .unzip();
    MapDispersions {
        dispersion,
        outlier,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trimmed_mean_matches_r() {
        // R: mean(c(5, 1, 9), trim = 1/3) == 5 (the median)
        assert_eq!(trimmed_mean(&[5.0, 1.0, 9.0], 1.0 / 3.0), 5.0);
        // R: mean(1:10, trim = 0.1) == 5.5
        let v: Vec<f64> = (1..=10).map(f64::from).collect();
        assert!((trimmed_mean(&v, 0.1) - 5.5).abs() < 1e-12);
    }

    #[test]
    fn parametric_trend_recovers_known_curve() {
        let means: Vec<f64> = (1..200).map(|i| i as f64 * 3.0).collect();
        let disps: Vec<f64> = means.iter().map(|m| 0.05 + 4.0 / m).collect();
        let (a0, a1) = parametric_trend(&means, &disps).expect("fit");
        assert!(
            (a0 - 0.05).abs() < 1e-6 && (a1 - 4.0).abs() < 1e-4,
            "{a0} {a1}"
        );
    }
}
