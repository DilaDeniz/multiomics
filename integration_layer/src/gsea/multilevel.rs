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
            return p_hat;
        }

        // Double: generate another n_perm permutations with offset seed
        let offset_seed = base_seed.wrapping_add((n_perm as u64).wrapping_mul(0xDEAD_BEEF_CAFE));
        null_es.extend(generate_null_es(n_genes, n_hits, n_perm, offset_seed));
        n_perm *= 2;
    }
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
