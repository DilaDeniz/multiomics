//! GSEA pre-ranked enrichment — fgsea multilevel Monte Carlo.
//!
//! Implements the adaptive multilevel sampling algorithm for accurate p-values
//! with far fewer permutations than the classic fixed-1000 approach.
//!
//! # References
//! * Korotkevich G, Sukhov V & Sergushichev A (2021). Fast gene set enrichment
//!   analysis. bioRxiv. <https://doi.org/10.1101/060012>
//! * Subramanian A, et al. (2005). Gene set enrichment analysis: A
//!   knowledge-based approach for interpreting genome-wide expression profiles.
//!   PNAS 102(43):15545–15550. <https://doi.org/10.1073/pnas.0506580102>

pub mod es_walk;
pub mod multilevel;
pub mod null_cache;
pub mod types;

pub use types::GseaResult;

use ahash::AHashSet;
use biomics_core::statistics::benjamini_hochberg;

use es_walk::compute_es;
use multilevel::multilevel_pvalue;
use null_cache::NullDistCache;

/// Relative precision target for the multilevel adaptive p-value estimator.
const EPS: f64 = 0.1;

/// Run GSEA pre-ranked enrichment analysis using the fgsea multilevel algorithm.
///
/// # Arguments
/// - `ranked_genes`: genes sorted **descending** by enrichment metric
///   (e.g. sign(log2FC) × −log10(p)).
/// - `pathways`: list of `(id, name, genes)` tuples.
/// - `min_size`: minimum number of pathway genes present in the ranked list.
/// - `max_size`: maximum number of pathway genes present in the ranked list.
/// - `n_perm`: used as the null-cache sample count (default ≥ 1 is enforced).
///   Adaptive multilevel sampling is then used per pathway; the cache provides
///   size-matched NES normalisation and a secondary p-value bound.
///
/// # Returns
/// A `Vec<GseaResult>` sorted by NES descending.
/// Benjamini-Hochberg FDR correction is applied across all returned results.
pub fn gsea_preranked(
    ranked_genes: &[(String, f64)],
    pathways: &[(&str, &str, &[&str])],
    min_size: usize,
    max_size: usize,
    n_perm: usize,
) -> Vec<GseaResult> {
    if ranked_genes.is_empty() || pathways.is_empty() {
        return Vec::new();
    }

    let cache_size = n_perm.max(1);
    let n_genes = ranked_genes.len();
    let base_seed: u64 = 0xBEEF_CAFE_1234_5678;

    // Normalize ranked gene names to uppercase once — avoids per-pathway allocation.
    let ranked_upper_owned: Vec<String> =
        ranked_genes.iter().map(|(g, _)| g.to_uppercase()).collect();
    let ranked_upper_refs: Vec<&str> = ranked_upper_owned.iter().map(|s| s.as_str()).collect();

    // Build size-matched null distribution cache (one distribution per unique k).
    let mut null_cache = NullDistCache::new();

    // Determine unique pathway sizes that pass the filter, then pre-warm cache.
    let unique_sizes: Vec<usize> = {
        let mut sizes: Vec<usize> = pathways
            .iter()
            .map(|&(_, _, pgenes)| {
                let set: AHashSet<&str> = pgenes.iter().copied().collect();
                // Compare against pre-uppercased ranked list.
                let set_upper: AHashSet<String> = set.iter().map(|g| g.to_uppercase()).collect();
                ranked_upper_refs
                    .iter()
                    .filter(|&&g| set_upper.contains(g))
                    .count()
            })
            .filter(|&k| k >= min_size && k <= max_size)
            .collect();
        sizes.sort_unstable();
        sizes.dedup();
        sizes
    };

    for k in &unique_sizes {
        null_cache.get_or_compute(*k, n_genes, cache_size, base_seed);
    }

    // Compute results for each eligible pathway.
    let mut results: Vec<GseaResult> = pathways
        .iter()
        .enumerate()
        .filter_map(|(pw_idx, &(pid, pname, pgenes))| {
            // Pathway genes uppercased once per pathway.
            let pathway_set: AHashSet<String> = pgenes.iter().map(|g| g.to_uppercase()).collect();

            let n_hits = ranked_upper_refs
                .iter()
                .filter(|&&g| pathway_set.contains(g))
                .count();

            if n_hits < min_size || n_hits > max_size {
                return None;
            }

            let is_hit: Vec<bool> = ranked_upper_refs
                .iter()
                .map(|&g| pathway_set.contains(g))
                .collect();

            let (es, running_sum, peak_idx) = compute_es(&is_hit, n_hits);

            // Per-pathway seed derived from index to keep results deterministic.
            let pw_seed = base_seed
                .wrapping_add(pw_idx as u64)
                .wrapping_mul(0x9E3779B97F4A7C15);

            // Adaptive multilevel p-value.
            let p_multilevel = multilevel_pvalue(es, n_genes, n_hits, EPS, pw_seed);

            // Size-matched cache p-value.
            let p_cache = null_cache
                .pvalue_from_cache(n_hits, es)
                .unwrap_or(p_multilevel);

            let p_value = p_multilevel.min(p_cache);

            // NES: normalise against size-matched null cache mean.
            let nes = null_cache.normalize_es(n_hits, es);

            let leading_edge: Vec<String> = if es >= 0.0 {
                ranked_upper_refs[..=peak_idx]
                    .iter()
                    .zip(ranked_genes[..=peak_idx].iter())
                    .filter(|(&upper, _)| pathway_set.contains(upper))
                    .map(|(_, (g, _))| g.clone())
                    .collect()
            } else {
                ranked_upper_refs[peak_idx..]
                    .iter()
                    .zip(ranked_genes[peak_idx..].iter())
                    .filter(|(&upper, _)| pathway_set.contains(upper))
                    .map(|(_, (g, _))| g.clone())
                    .collect()
            };

            Some(GseaResult {
                pathway_id: pid.to_string(),
                pathway_name: pname.to_string(),
                es,
                nes,
                p_value,
                padj: f64::NAN,
                n_genes_pathway: n_hits,
                n_genes_ranked: n_genes,
                leading_edge,
                running_sum,
            })
        })
        .collect();

    if results.is_empty() {
        return Vec::new();
    }

    // BH FDR correction.
    let pvals: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    let padj_vals = benjamini_hochberg(&pvals);
    for (r, padj) in results.iter_mut().zip(padj_vals) {
        r.padj = padj;
    }

    // Sort by NES descending (most positively enriched first).
    results.sort_unstable_by(|a, b| {
        b.nes
            .partial_cmp(&a.nes)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_ranked(n: usize) -> Vec<(String, f64)> {
        (0..n)
            .map(|i| (format!("GENE{i:04}"), (n as f64) - (i as f64)))
            .collect()
    }

    #[test]
    fn test_perfect_pathway_at_top_is_enriched() {
        let n = 200;
        let ranked = make_ranked(n);
        let pathway_genes: Vec<&str> = ranked[..20].iter().map(|(g, _)| g.as_str()).collect();
        let pathways: Vec<(&str, &str, &[&str])> =
            vec![("PW_TOP", "Top-loaded pathway", pathway_genes.as_slice())];
        let results = gsea_preranked(&ranked, &pathways, 5, 500, 100);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert!(r.es > 0.0, "ES={}", r.es);
        assert!(r.nes > 1.0, "NES={}", r.nes);
        assert!(r.p_value <= 0.10, "p={}", r.p_value);
        assert!(!r.leading_edge.is_empty());
    }

    #[test]
    fn test_pathway_filtered_by_size_bounds() {
        let ranked = make_ranked(100);
        let tiny_genes: Vec<&str> = ranked[..2].iter().map(|(g, _)| g.as_str()).collect();
        let pathways = vec![("TINY", "Tiny pathway", tiny_genes.as_slice())];
        let results = gsea_preranked(&ranked, &pathways, 5, 500, 10);
        assert!(results.is_empty());
    }

    #[test]
    fn test_running_sum_length_matches_ranked() {
        let n = 50;
        let ranked = make_ranked(n);
        let pathway_genes: Vec<&str> = ranked[..10].iter().map(|(g, _)| g.as_str()).collect();
        let pathways = vec![("PW", "Test", pathway_genes.as_slice())];
        let results = gsea_preranked(&ranked, &pathways, 5, 500, 50);
        assert_eq!(results[0].running_sum.len(), n);
    }
}
