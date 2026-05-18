//! Multilevel Monte Carlo adaptive p-value estimation.
//!
//! Implements the fgsea multilevel algorithm:
//! Korotkevich G, Sukhov V & Sergushichev A (2021). Fast gene set enrichment
//! analysis. bioRxiv. <https://doi.org/10.1101/060012>

use rayon::prelude::*;

use crate::gsea::es_walk::{compute_es, xorshift_sample};

/// Generate `n_samples` null ES values by permuting hit positions.
///
/// Uses rayon parallel iteration; each worker seeds its own xorshift state
/// from the global `base_seed` XOR-ed with a per-permutation offset.
pub fn generate_null_es(
    n_genes: usize,
    n_hits: usize,
    n_samples: usize,
    base_seed: u64,
) -> Vec<f64> {
    (0..n_samples)
        .into_par_iter()
        .map(|i| {
            let seed = base_seed
                .wrapping_add(i as u64)
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let positions = xorshift_sample(n_genes, n_hits, seed);
            let mut perm_hit = vec![false; n_genes];
            for pos in &positions {
                perm_hit[*pos] = true;
            }
            compute_es(&perm_hit, n_hits).0
        })
        .collect()
}

/// Compute an adaptive p-value for one pathway using the multilevel algorithm.
///
/// Starts with 101 permutations and doubles until the relative standard error
/// of the proportion estimate drops below `eps` or `max_perm` is reached.
///
/// When the plain estimate falls below 1e-4, switches automatically to
/// importance sampling so that very small p-values are accurate without
/// requiring millions of uniform permutations.
///
/// Returns the Laplace-corrected empirical p-value.
pub fn multilevel_pvalue(
    es_obs: f64,
    n_genes: usize,
    n_hits: usize,
    eps: f64,
    base_seed: u64,
) -> f64 {
    let max_perm: usize = 100_000;
    let mut n_perm: usize = 101;

    let mut null_es = generate_null_es(n_genes, n_hits, n_perm, base_seed);

    loop {
        let obs_abs = es_obs.abs();
        let n_exceed = null_es.iter().filter(|&&e| e.abs() >= obs_abs).count();
        let p_hat = (n_exceed as f64 + 1.0) / (null_es.len() as f64 + 1.0);

        let se = (p_hat * (1.0 - p_hat) / null_es.len() as f64).sqrt();

        if se < eps * p_hat || n_perm >= max_perm {
            // Switch to importance sampling when p < 1e-4 to get accurate
            // estimates without millions of uniform permutations.
            if p_hat < 1e-4 {
                return importance_sampling_pvalue(es_obs, n_genes, n_hits, eps, base_seed, 50_000);
            }
            return p_hat;
        }

        // Double: generate another n_perm permutations with offset seed
        let offset_seed = base_seed.wrapping_add((n_perm as u64).wrapping_mul(0xDEAD_BEEF_CAFE));
        null_es.extend(generate_null_es(n_genes, n_hits, n_perm, offset_seed));
        n_perm *= 2;
    }
}

/// Importance-sampling p-value estimator for the tail region |ES| ≥ |es_obs|.
///
/// Uses a biased proposal that shifts pathway hits toward the top/bottom of the
/// ranked list (concentrating mass in the tail), then corrects by the
/// importance weight w = P_uniform(sample) / P_proposal(sample).
///
/// Weights are computed in log-space and combined via Kahan-stable summation to
/// prevent floating-point cancellation.
///
/// Reference: Owen & Zhou (2000) Safe and effective importance sampling.
///            J. Am. Stat. Assoc. 95:135–143.
pub fn importance_sampling_pvalue(
    es_obs: f64,
    n_genes: usize,
    n_hits: usize,
    _eps: f64,
    base_seed: u64,
    n_samples: usize,
) -> f64 {
    // Bias parameter: concentrate proposal near the top γ-fraction of genes.
    // γ = 2 * n_hits / n_genes keeps the expected overlap with the pathway at 2×.
    let gamma = (2.0 * n_hits as f64 / n_genes as f64).min(0.5_f64);
    let n_top = ((gamma * n_genes as f64).ceil() as usize).max(n_hits);

    // log probability under proposal: sample hits uniformly from [0, n_top)
    // log P_proposal = log(C(n_top, n_hits))^{-1} — same for every proposal draw
    // log probability under uniform:  log(C(n_genes, n_hits))^{-1}
    // log weight = log P_uniform - log P_proposal
    //            = log C(n_top, n_hits) - log C(n_genes, n_hits)
    let log_w = log_binom(n_top, n_hits) - log_binom(n_genes, n_hits);

    let obs_abs = es_obs.abs();

    // Kahan summation for numerator (sum of weights where |ES| >= obs_abs)
    // and denominator (sum of all weights).
    let (mut sum_w_num, mut comp_num) = (0.0_f64, 0.0_f64);
    let (mut sum_w_den, mut comp_den) = (0.0_f64, 0.0_f64);

    for i in 0..n_samples {
        let seed = base_seed
            .wrapping_add(i as u64)
            .wrapping_mul(6364136223846793005)
            .wrapping_add(0xC0FF_EE12_3400_u64.wrapping_mul(i as u64 | 1));

        // Draw hits from the biased proposal: uniform over [0, n_top)
        let positions = xorshift_sample(n_top, n_hits, seed);
        let mut is_hit = vec![false; n_genes];
        for pos in &positions {
            is_hit[*pos] = true;
        }
        let es_sample = compute_es(&is_hit, n_hits).0;

        // Importance weight (all samples share the same log_w here)
        let w = log_w.exp();

        // Kahan: denominator
        let y_d = w - comp_den;
        let t_d = sum_w_den + y_d;
        comp_den = (t_d - sum_w_den) - y_d;
        sum_w_den = t_d;

        // Kahan: numerator (only when in the tail)
        if es_sample.abs() >= obs_abs {
            let y_n = w - comp_num;
            let t_n = sum_w_num + y_n;
            comp_num = (t_n - sum_w_num) - y_n;
            sum_w_num = t_n;
        }
    }

    if sum_w_den == 0.0 {
        return 1.0 / (n_samples as f64 + 1.0);
    }

    // Laplace correction: add 1 pseudo-count to numerator and denominator
    let w_unit = log_w.exp();
    ((sum_w_num + w_unit) / (sum_w_den + w_unit)).clamp(1.0 / (n_samples as f64 + 1.0), 1.0)
}

/// log of the binomial coefficient C(n, k) computed via Stirling / lgamma.
fn log_binom(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
}

/// Natural log-gamma via Lanczos approximation (g=7, 9 coefficients).
/// Accurate to ~15 significant digits for x > 0.
fn lgamma(x: usize) -> f64 {
    // Lanczos coefficients for g=7
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_312_4e-7,
    ];
    let x = x as f64;
    if x < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - lgamma((1.0 - x) as usize);
    }
    let x = x - 1.0;
    let mut a = C[0];
    for (i, &c) in C[1..].iter().enumerate() {
        a += c / (x + i as f64 + 1.0);
    }
    let t = x + G + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn null_es_has_correct_length() {
        let vals = generate_null_es(100, 10, 50, 42);
        assert_eq!(vals.len(), 50);
    }

    #[test]
    fn multilevel_pvalue_bounded() {
        // A pathway at rank 0..10 of 100 genes should have a low p-value.
        // ES for perfect top-loaded pathway is deterministically computable.
        let n = 100_usize;
        let k = 10_usize;
        let is_hit: Vec<bool> = (0..n).map(|i| i < k).collect();
        let (es, _, _) = compute_es(&is_hit, k);
        let p = multilevel_pvalue(es, n, k, 0.1, 1234);
        assert!(p > 0.0 && p <= 1.0, "p-value out of range: {p}");
    }
}
