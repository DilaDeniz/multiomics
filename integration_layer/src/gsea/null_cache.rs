//! Null ES distribution cache: one sorted distribution per pathway size.
//!
//! Pathways of the same size `k` share the same null ES distribution under
//! label-permutation testing.  Computing the distribution once and reusing it
//! for all pathways of that size is the primary speed-up in fgsea.
//!
//! Reference:
//! Korotkevich G, Sukhov V & Sergushichev A (2021). Fast gene set enrichment
//! analysis. bioRxiv. <https://doi.org/10.1101/060012>

use ahash::AHashMap;

use crate::gsea::multilevel::generate_null_es;

/// Cache of sorted absolute null-ES vectors keyed by pathway size.
pub struct NullDistCache {
    /// key = pathway size k, value = sorted |ES_null| values (ascending).
    cache: AHashMap<usize, Vec<f64>>,
}

impl NullDistCache {
    /// Create an empty cache.
    pub fn new() -> Self {
        Self {
            cache: AHashMap::new(),
        }
    }

    /// Return (or compute and cache) the sorted null |ES| distribution for
    /// pathway size `k` in a ranked list of `n_genes` genes.
    ///
    /// `n_samples` controls how many permutations are generated on first call.
    pub fn get_or_compute(
        &mut self,
        k: usize,
        n_genes: usize,
        n_samples: usize,
        base_seed: u64,
    ) -> &[f64] {
        self.cache.entry(k).or_insert_with(|| {
            let mut vals: Vec<f64> = generate_null_es(n_genes, k, n_samples, base_seed)
                .into_iter()
                .map(f64::abs)
                .collect();
            vals.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            vals
        })
    }

    /// Compute p-value for `es_obs` by binary-searching the cached distribution.
    ///
    /// Returns a Laplace-corrected p-value so that p is never exactly 0.
    pub fn pvalue_from_cache(&self, k: usize, es_obs: f64) -> Option<f64> {
        let sorted = self.cache.get(&k)?;
        let obs_abs = es_obs.abs();
        // partition_point gives first index where sorted[i] >= obs_abs
        let n_below = sorted.partition_point(|&x| x < obs_abs);
        let n_exceed = sorted.len() - n_below;
        let p = (n_exceed as f64 + 1.0) / (sorted.len() as f64 + 1.0);
        Some(p)
    }

    /// Compute NES by normalising `es` against the mean of the cached null |ES|.
    pub fn normalize_es(&self, k: usize, es: f64) -> f64 {
        let Some(sorted) = self.cache.get(&k) else {
            return es;
        };
        if sorted.is_empty() {
            return es;
        }
        let mean_null: f64 = sorted.iter().sum::<f64>() / sorted.len() as f64;
        if mean_null < 1e-9 {
            es
        } else {
            es / mean_null
        }
    }
}

impl Default for NullDistCache {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cache_computes_and_returns_sorted_distribution() {
        let mut cache = NullDistCache::new();
        let dist = cache.get_or_compute(10, 100, 200, 42);
        assert_eq!(dist.len(), 200);
        // Verify ascending sort
        for w in dist.windows(2) {
            assert!(w[0] <= w[1], "distribution not sorted: {} > {}", w[0], w[1]);
        }
    }

    #[test]
    fn pvalue_from_cache_is_in_range() {
        let mut cache = NullDistCache::new();
        cache.get_or_compute(10, 100, 500, 99);
        // A very large ES should give near-minimum p
        let p = cache.pvalue_from_cache(10, 1e9).unwrap();
        assert!(p > 0.0 && p <= 1.0, "p out of [0,1]: {p}");
    }
}
