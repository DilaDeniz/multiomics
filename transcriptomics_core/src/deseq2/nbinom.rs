//! Negative-binomial likelihood and dispersion optimisation.
//!
//! Implements the dispersion estimation of Love, Huber & Anders (2014),
//! "Moderated estimation of fold change and dispersion for RNA-seq data with
//! DESeq2", Genome Biology 15:550: the Cox-Reid adjusted profile likelihood of
//! log α, optionally with a normal prior on log α, maximised by a backtracking
//! (Armijo) line search along the gradient with a grid-search fallback.
//! Tuning constants follow the DESeq2 defaults so results can be compared
//! one-to-one with R/DESeq2.

use biomics_core::statistics::{digamma, ln_gamma};

use super::linalg::{chol_inverse, chol_logdet, cholesky, trace_of_product, xtwx};

/// Below this count, ln Γ(y + r) − ln Γ(r) and ψ(y + r) − ψ(r) are evaluated
/// as exact finite sums, which avoids cancellation when r = 1/α is huge.
const EXACT_SUM_LIMIT: f64 = 64.0;

/// Σ_{k=0}^{y−1} ln((r + k) / (r + μ)) = ln Γ(y + r) − ln Γ(r) − y ln(r + μ).
fn lgamma_ratio_minus_log(y: f64, r: f64, mu: f64) -> f64 {
    if y == 0.0 {
        return 0.0;
    }
    if y < EXACT_SUM_LIMIT && y.fract() == 0.0 {
        let denom = r + mu;
        (0..y as u32)
            .map(|k| ((k as f64 - mu) / denom).ln_1p())
            .sum()
    } else {
        ln_gamma(y + r) - ln_gamma(r) - y * (r + mu).ln()
    }
}

/// ψ(r) − ψ(y + r).
fn digamma_diff(y: f64, r: f64) -> f64 {
    if y == 0.0 {
        return 0.0;
    }
    if y < EXACT_SUM_LIMIT && y.fract() == 0.0 {
        -(0..y as u32).map(|k| 1.0 / (r + k as f64)).sum::<f64>()
    } else {
        digamma(r) - digamma(y + r)
    }
}

/// log P(Y = y) for Y ~ NB(mean μ, size r) (variance μ + μ²/r).
pub(crate) fn nb_log_pmf(y: f64, size: f64, mu: f64) -> f64 {
    if mu <= 0.0 {
        return if y == 0.0 { 0.0 } else { f64::NEG_INFINITY };
    }
    let size_term = -size * (mu / size).ln_1p();
    if y == 0.0 {
        return size_term;
    }
    lgamma_ratio_minus_log(y, size, mu) + y * mu.ln() + size_term - ln_gamma(y + 1.0)
}

/// Inputs shared by every evaluation of one gene's dispersion posterior.
pub(crate) struct DispersionProblem<'a> {
    pub y: &'a [f64],
    pub mu: &'a [f64],
    pub x: &'a [Vec<f64>],
    /// Prior mean of log α (ignored unless `prior_var` is set).
    pub prior_mean: f64,
    /// Prior variance of log α; `None` for maximum likelihood.
    pub prior_var: Option<f64>,
}

impl DispersionProblem<'_> {
    /// Cox-Reid adjusted log posterior of log α (data-only constants dropped).
    pub(crate) fn log_posterior(&self, log_alpha: f64) -> f64 {
        let alpha = log_alpha.exp();
        let r = 1.0 / alpha;
        let w: Vec<f64> = self.mu.iter().map(|&m| 1.0 / (1.0 / m + alpha)).collect();
        let cr = match cholesky(&xtwx(self.x, &w)) {
            Some(l) => -0.5 * chol_logdet(&l),
            None => return f64::NAN,
        };
        let ll: f64 = self
            .y
            .iter()
            .zip(self.mu)
            .map(|(&y, &m)| lgamma_ratio_minus_log(y, r, m) - r * (m * alpha).ln_1p())
            .sum();
        let prior = self
            .prior_var
            .map_or(0.0, |v| -0.5 * (log_alpha - self.prior_mean).powi(2) / v);
        ll + cr + prior
    }

    /// Derivative of [`Self::log_posterior`] with respect to log α.
    pub(crate) fn dlog_posterior(&self, log_alpha: f64) -> f64 {
        let alpha = log_alpha.exp();
        let r = 1.0 / alpha;
        let w: Vec<f64> = self.mu.iter().map(|&m| 1.0 / (1.0 / m + alpha)).collect();
        let dw: Vec<f64> = w.iter().map(|&wj| -wj * wj).collect();
        let cr = match cholesky(&xtwx(self.x, &w)) {
            // d/dα [−½ ln det B] = −½ tr(B⁻¹ dB/dα)
            Some(l) => -0.5 * trace_of_product(&chol_inverse(&l), &xtwx(self.x, &dw)),
            None => return f64::NAN,
        };
        let ll = r
            * r
            * self
                .y
                .iter()
                .zip(self.mu)
                .map(|(&y, &m)| {
                    let ma = m * alpha;
                    digamma_diff(y, r) + ma.ln_1p() - ma / (1.0 + ma) + y / (m + r)
                })
                .sum::<f64>();
        let prior = self
            .prior_var
            .map_or(0.0, |v| -(log_alpha - self.prior_mean) / v);
        (ll + cr) * alpha + prior
    }
}

/// Result of the line search for one gene.
pub(crate) struct LineSearchFit {
    pub log_alpha: f64,
    pub iterations: usize,
    pub initial_lp: f64,
    pub last_lp: f64,
}

/// Line-search settings (DESeq2 defaults: κ₀ = 1, tol = 1e-6, maxit = 100,
/// min log α = ln(minDisp / 10)).
pub(crate) struct LineSearch {
    pub kappa_0: f64,
    pub tol: f64,
    pub max_iter: usize,
    pub min_log_alpha: f64,
}

/// Armijo sufficient-increase constant.
const ARMIJO_EPSILON: f64 = 1.0e-4;

impl LineSearch {
    /// Maximise the log posterior over log α by gradient ascent with an Armijo
    /// backtracking line search, starting from `log_alpha`.
    pub(crate) fn run(&self, problem: &DispersionProblem<'_>, log_alpha: f64) -> LineSearchFit {
        let mut a = log_alpha;
        let mut lp = problem.log_posterior(a);
        let mut dlp = problem.dlog_posterior(a);
        let initial_lp = lp;
        let mut kappa = self.kappa_0;
        let mut iterations = 0;
        let mut accepted = 0usize;

        for _ in 0..self.max_iter {
            iterations += 1;
            // Keep log α inside [−30, 10] where ln Γ stays numerically stable.
            let proposal = a + kappa * dlp;
            if proposal < -30.0 {
                kappa = (-30.0 - a) / dlp;
            }
            if proposal > 10.0 {
                kappa = (10.0 - a) / dlp;
            }
            let theta = -problem.log_posterior(a + kappa * dlp);
            let theta_hat = -lp - kappa * ARMIJO_EPSILON * dlp * dlp;
            if theta <= theta_hat {
                accepted += 1;
                a += kappa * dlp;
                let lp_new = problem.log_posterior(a);
                let change = lp_new - lp;
                lp = lp_new;
                if change < self.tol || a < self.min_log_alpha {
                    break;
                }
                dlp = problem.dlog_posterior(a);
                kappa = (kappa * 1.1).min(self.kappa_0);
                if accepted.is_multiple_of(5) {
                    kappa /= 2.0;
                }
            } else {
                kappa /= 2.0;
            }
        }

        LineSearchFit {
            log_alpha: a,
            iterations,
            initial_lp,
            last_lp: lp,
        }
    }
}

/// Grid search for log α: a 20-point coarse grid over [`min_log_alpha`,
/// `max_log_alpha`], then a 20-point fine grid of ± one coarse step around the
/// best coarse point.
pub(crate) fn grid_search(
    problem: &DispersionProblem<'_>,
    min_log_alpha: f64,
    max_log_alpha: f64,
) -> f64 {
    const N: usize = 20;
    let linspace = |lo: f64, hi: f64| -> Vec<f64> {
        (0..N)
            .map(|i| lo + (hi - lo) * i as f64 / (N - 1) as f64)
            .collect()
    };
    let argmax = |grid: &[f64]| -> f64 {
        let mut best = grid[0];
        let mut best_lp = f64::NEG_INFINITY;
        for &a in grid {
            let lp = problem.log_posterior(a);
            if lp > best_lp {
                best_lp = lp;
                best = a;
            }
        }
        best
    };
    let coarse = linspace(min_log_alpha, max_log_alpha);
    let delta = coarse[1] - coarse[0];
    let a_hat = argmax(&coarse);
    argmax(&linspace(a_hat - delta, a_hat + delta))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nb_log_pmf_matches_reference_values() {
        // R 4.2: dnbinom(c(0, 3, 50), size = 2, mu = 10, log = TRUE)
        let r = [
            -3.583_518_938_456_11,
            -2.744_189_247_718_084,
            -8.767_771_145_429_514,
        ];
        for (y, expected) in [0.0, 3.0, 50.0].iter().zip(r) {
            let got = nb_log_pmf(*y, 2.0, 10.0);
            assert!((got - expected).abs() < 1e-12, "y={y}: {got} vs {expected}");
        }
        // Near-Poisson regime (size 1e10, mu 4, y 2): the exact value exceeds
        // dpois(2, 4, log = TRUE) by ~1e-10. (R 4.2's dnbinom is off by 4e-8 here.)
        let pois = 2.0 * 4f64.ln() - 4.0 - 2f64.ln();
        let got = nb_log_pmf(2.0, 1e10, 4.0);
        assert!((got - pois).abs() < 1e-9, "{got} vs {pois}");
    }

    #[test]
    fn dlog_posterior_matches_finite_difference() {
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let y = [3.0, 9.0, 40.0, 25.0];
        let mu = [6.0, 6.0, 32.5, 32.5];
        for prior_var in [None, Some(0.5)] {
            let p = DispersionProblem {
                y: &y,
                mu: &mu,
                x: &x,
                prior_mean: -1.0,
                prior_var,
            };
            for a in [-4.0, -1.5, 0.3] {
                let h = 1e-5;
                let fd = (p.log_posterior(a + h) - p.log_posterior(a - h)) / (2.0 * h);
                let an = p.dlog_posterior(a);
                assert!((fd - an).abs() < 1e-6 * an.abs().max(1.0), "{fd} vs {an}");
            }
        }
    }
}
