//! Negative-binomial GLM fitting with a log link and fixed dispersion.
//!
//! Coefficients are fitted by iteratively reweighted least squares with a small
//! ridge penalty (DESeq2's "wide prior", λ = 1e-6 on the log₂ scale), stopping
//! when the relative change in deviance falls below `tol`. Genes that fail to
//! converge, diverge (|β| > 30 on the natural-log scale) or yield non-positive
//! variances are refitted by a bounded Newton method on the same objective.

use super::linalg::{chol_inverse, chol_solve, cholesky, mat_mul, xtwx, xtwz, Matrix};
use super::nbinom::nb_log_pmf;

/// Ridge penalty on the log₂ scale used for every coefficient (DESeq2 default).
pub(crate) const RIDGE_LOG2: f64 = 1e-6;
/// |β| bound: divergence threshold for IRLS (natural-log scale) and box
/// constraint for the Newton fallback (log₂ scale), as in DESeq2.
const LARGE: f64 = 30.0;

/// GLM fit for one gene. Coefficients and variances are on the natural-log scale.
pub(crate) struct GlmFit {
    pub beta: Vec<f64>,
    /// Sandwich covariance (XᵀWX + Λ)⁻¹ XᵀWX (XᵀWX + Λ)⁻¹.
    pub sigma: Matrix,
    /// Diagonal of the hat matrix W½X(XᵀWX + Λ)⁻¹XᵀW½.
    pub hat: Vec<f64>,
    /// Fitted means s_j·exp(x_jᵀβ), not floored at `min_mu`.
    pub mu: Vec<f64>,
    pub converged: bool,
}

/// Settings shared by every gene's fit.
pub(crate) struct GlmSettings {
    pub tol: f64,
    pub max_iter: usize,
    pub min_mu: f64,
}

impl Default for GlmSettings {
    fn default() -> Self {
        Self {
            tol: 1e-8,
            max_iter: 100,
            min_mu: 0.5,
        }
    }
}

fn linear_predictor(x: &[Vec<f64>], beta: &[f64]) -> Vec<f64> {
    x.iter()
        .map(|row| row.iter().zip(beta).map(|(a, b)| a * b).sum())
        .collect()
}

fn means(x: &[Vec<f64>], beta: &[f64], sf: &[f64]) -> Vec<f64> {
    linear_predictor(x, beta)
        .iter()
        .zip(sf)
        .map(|(eta, s)| s * eta.exp())
        .collect()
}

fn deviance(y: &[f64], mu: &[f64], alpha: f64) -> f64 {
    let size = 1.0 / alpha;
    -2.0 * y
        .iter()
        .zip(mu)
        .map(|(&k, &m)| nb_log_pmf(k, size, m))
        .sum::<f64>()
}

/// Covariance, hat diagonal and unfloored means at coefficients `beta`.
#[allow(clippy::too_many_arguments)]
fn finish(
    y_len: usize,
    x: &[Vec<f64>],
    sf: &[f64],
    alpha: f64,
    beta: Vec<f64>,
    ridge: &[f64],
    min_mu: f64,
    converged: bool,
) -> GlmFit {
    let mu = means(x, &beta, sf);
    let w: Vec<f64> = mu
        .iter()
        .map(|&m| {
            let m = m.max(min_mu);
            m / (1.0 + alpha * m)
        })
        .collect();
    let info = xtwx(x, &w);
    let mut penalised = info.clone();
    for (i, r) in ridge.iter().enumerate() {
        penalised[i][i] += r;
    }
    let p = beta.len();
    let (sigma, hat) = match cholesky(&penalised) {
        Some(l) => {
            let inv = chol_inverse(&l);
            let sigma = mat_mul(&mat_mul(&inv, &info), &inv);
            let hat = (0..y_len)
                .map(|j| {
                    let xj = &x[j];
                    let q: f64 = (0..p)
                        .map(|a| (0..p).map(|b| xj[a] * inv[a][b] * xj[b]).sum::<f64>())
                        .sum();
                    w[j] * q
                })
                .collect();
            (sigma, hat)
        }
        None => (vec![vec![f64::NAN; p]; p], vec![f64::NAN; y_len]),
    };
    GlmFit {
        beta,
        sigma,
        hat,
        mu,
        converged,
    }
}

/// Closed-form fit for an intercept-only design, as DESeq2 does: the
/// coefficient is the log of the mean normalised count (no IRLS).
pub(crate) fn fit_intercept_only(sf: &[f64], mean_norm: f64, alpha: f64) -> GlmFit {
    let mu: Vec<f64> = sf.iter().map(|s| s * mean_norm).collect();
    let w: Vec<f64> = mu.iter().map(|&m| 1.0 / (1.0 / m + alpha)).collect();
    let info: f64 = w.iter().sum();
    GlmFit {
        beta: vec![mean_norm.ln()],
        sigma: vec![vec![1.0 / info]],
        hat: w.iter().map(|wj| wj / info).collect(),
        mu,
        converged: true,
    }
}

/// Fit the NB GLM for one gene by ridge-penalised IRLS.
pub(crate) fn fit_beta(
    y: &[f64],
    sf: &[f64],
    x: &[Vec<f64>],
    alpha: f64,
    beta_init: &[f64],
    settings: &GlmSettings,
) -> GlmFit {
    let p = beta_init.len();
    let ridge = vec![RIDGE_LOG2 / std::f64::consts::LN_2.powi(2); p];
    let floor =
        |mu: Vec<f64>| -> Vec<f64> { mu.into_iter().map(|m| m.max(settings.min_mu)).collect() };

    let mut beta = beta_init.to_vec();
    let mut mu = floor(means(x, &beta, sf));
    let mut dev_old = 0.0;
    let mut converged = false;

    for t in 0..settings.max_iter {
        let w: Vec<f64> = mu.iter().map(|&m| m / (1.0 + alpha * m)).collect();
        let z: Vec<f64> = y
            .iter()
            .zip(&mu)
            .zip(sf)
            .map(|((&k, &m), &s)| (m / s).ln() + (k - m) / m)
            .collect();
        let mut a = xtwx(x, &w);
        for (i, r) in ridge.iter().enumerate() {
            a[i][i] += r;
        }
        let Some(l) = cholesky(&a) else { break };
        beta = chol_solve(&l, &xtwz(x, &w, &z));
        // Divergence check on the natural-log scale, as in DESeq2's IRLS.
        if beta.iter().any(|b| !b.is_finite() || b.abs() > LARGE) {
            break;
        }
        mu = floor(means(x, &beta, sf));
        let dev = deviance(y, &mu, alpha);
        let conv = (dev - dev_old).abs() / (dev.abs() + 0.1);
        if conv.is_nan() {
            break;
        }
        if t > 0 && conv < settings.tol {
            converged = true;
            break;
        }
        dev_old = dev;
    }

    let fit = finish(
        y.len(),
        x,
        sf,
        alpha,
        beta,
        &ridge,
        settings.min_mu,
        converged,
    );
    let stable = fit.beta.iter().all(|b| b.is_finite());
    let var_positive = (0..p).all(|i| fit.sigma[i][i] > 0.0);
    if fit.converged && stable && var_positive {
        return fit;
    }
    let start = if stable
        && fit
            .beta
            .iter()
            .all(|b| (b / std::f64::consts::LN_2).abs() < LARGE)
    {
        fit.beta.clone()
    } else {
        beta_init.to_vec()
    };
    bounded_newton(y, sf, x, alpha, &start, &ridge, settings)
}

/// Maximise the ridge-penalised NB log-likelihood with |β| ≤ 30 (log₂ scale)
/// by projected Fisher scoring with step halving.
fn bounded_newton(
    y: &[f64],
    sf: &[f64],
    x: &[Vec<f64>],
    alpha: f64,
    start: &[f64],
    ridge: &[f64],
    settings: &GlmSettings,
) -> GlmFit {
    let bound = LARGE * std::f64::consts::LN_2;
    let size = 1.0 / alpha;
    let objective = |beta: &[f64]| -> f64 {
        let mu = means(x, beta, sf);
        let ll: f64 = y
            .iter()
            .zip(&mu)
            .map(|(&k, &m)| nb_log_pmf(k, size, m))
            .sum();
        let pen: f64 = beta.iter().zip(ridge).map(|(b, r)| 0.5 * r * b * b).sum();
        ll - pen
    };
    let mut beta: Vec<f64> = start.iter().map(|b| b.clamp(-bound, bound)).collect();
    let mut f = objective(&beta);
    let mut converged = false;
    for _ in 0..500 {
        let mu = means(x, &beta, sf);
        // Score and Fisher information of the penalised log-likelihood.
        let resid: Vec<f64> = y
            .iter()
            .zip(&mu)
            .map(|(&k, &m)| (k - m) / (1.0 + alpha * m))
            .collect();
        let ones = vec![1.0; y.len()];
        let mut grad = xtwz(x, &ones, &resid);
        for ((g, b), r) in grad.iter_mut().zip(&beta).zip(ridge) {
            *g -= r * b;
        }
        let w: Vec<f64> = mu.iter().map(|&m| m / (1.0 + alpha * m)).collect();
        let mut info = xtwx(x, &w);
        for (i, r) in ridge.iter().enumerate() {
            info[i][i] += r;
        }
        let Some(l) = cholesky(&info) else { break };
        let step = chol_solve(&l, &grad);
        let mut t = 1.0;
        let mut improved = false;
        while t > 1e-10 {
            let cand: Vec<f64> = beta
                .iter()
                .zip(&step)
                .map(|(b, s)| (b + t * s).clamp(-bound, bound))
                .collect();
            let fc = objective(&cand);
            if fc >= f {
                let gain = fc - f;
                beta = cand;
                f = fc;
                improved = gain > 1e-12 * f.abs().max(1.0);
                break;
            }
            t /= 2.0;
        }
        if !improved {
            converged = true;
            break;
        }
    }
    finish(
        y.len(),
        x,
        sf,
        alpha,
        beta,
        ridge,
        settings.min_mu,
        converged,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_group_fit_recovers_group_means() {
        // With negligible dispersion the MLE equals log(group mean / size factor).
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let y = [10.0, 14.0, 40.0, 48.0];
        let sf = [1.0; 4];
        let fit = fit_beta(&y, &sf, &x, 1e-6, &[2.0, 1.0], &GlmSettings::default());
        assert!(fit.converged);
        assert!((fit.beta[0] - 12f64.ln()).abs() < 1e-5);
        assert!((fit.beta[1] - (44f64 / 12.0).ln()).abs() < 1e-5);
        assert!(fit.hat.iter().all(|h| (0.0..=1.0).contains(h)));
    }

    #[test]
    fn all_zero_group_is_bounded() {
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let fit = fit_beta(
            &[0.0, 0.0, 20.0, 25.0],
            &[1.0; 4],
            &x,
            0.1,
            &[-2.3, 5.4],
            &GlmSettings::default(),
        );
        assert!(fit.beta.iter().all(|b| b.is_finite()));
        assert!((fit.beta[1] / std::f64::consts::LN_2).abs() <= LARGE + 1e-9);
    }
}
