//! GSEA pre-ranked enrichment (Subramanian et al. 2005 classic GSEA).
//!
//! Implements the enrichment score (ES) walk, NES normalization, empirical
//! p-value estimation via permutation testing (using Rayon for parallelism),
//! and Benjamini-Hochberg FDR correction.
//!
//! # References
//! Subramanian A, et al. (2005). Gene set enrichment analysis: A knowledge-based
//! approach for interpreting genome-wide expression profiles.
//! PNAS 102(43):15545–15550. <https://doi.org/10.1073/pnas.0506580102>
//!
//! # Statistical note on p-values
//! This implementation uses empirical permutation-based p-values (default 1000
//! permutations). The permutation scheme permutes gene labels, recomputes ES,
//! and estimates p = fraction(|ES_perm| >= |ES_obs|). This is the same strategy
//! used in the original GSEA software. For very small pathways the resolution is
//! limited by `n_perm`. Publication-grade analyses should use n_perm ≥ 10 000.

use ahash::AHashSet;
use biomics_core::statistics::benjamini_hochberg;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

// ── Public types ──────────────────────────────────────────────────────────────

/// Result of GSEA pre-ranked enrichment for one gene set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GseaResult {
    /// Pathway / gene-set identifier (e.g. "hsa04110").
    pub pathway_id: String,
    /// Human-readable name.
    pub pathway_name: String,
    /// Enrichment score: maximum-deviation point of the running sum statistic.
    pub es: f64,
    /// Normalized enrichment score: ES / mean(|ES_null|).
    pub nes: f64,
    /// Empirical permutation p-value.
    pub p_value: f64,
    /// Benjamini-Hochberg adjusted p-value.
    pub padj: f64,
    /// Number of pathway genes that appear in the ranked list.
    pub n_genes_pathway: usize,
    /// Total number of genes in the ranked list.
    pub n_genes_ranked: usize,
    /// Leading-edge genes: genes contributing to the ES peak (for positive ES,
    /// genes appearing before the peak position in the ranked list).
    pub leading_edge: Vec<String>,
    /// Running sum values at every position in the ranked list (for plotting).
    pub running_sum: Vec<f64>,
}

// ── Core enrichment-score calculation ────────────────────────────────────────

/// Compute the enrichment score and the full running-sum trace for one gene set.
///
/// Returns `(es, running_sum, peak_index)` where `peak_index` is the index of
/// the maximum absolute deviation in the running sum.
///
/// `is_hit[i]` is `true` iff ranked gene `i` belongs to the pathway.
/// `n_hits` is the number of hits (pre-computed for efficiency).
fn compute_es(is_hit: &[bool], n_hits: usize) -> (f64, Vec<f64>, usize) {
    let n = is_hit.len();
    let n_miss = n - n_hits;

    if n_hits == 0 || n_miss == 0 {
        // Degenerate cases
        let rs = vec![0.0; n];
        return (0.0, rs, 0);
    }

    // Step contributions per position
    let hit_inc = ((n_miss as f64) / (n_hits as f64)).sqrt();
    let miss_dec = ((n_hits as f64) / (n_miss as f64)).sqrt();

    let mut running = 0.0_f64;
    let mut peak_val = 0.0_f64; // tracks maximum absolute value (with sign)
    let mut peak_idx = 0_usize;
    let mut rs = Vec::with_capacity(n);

    for (i, &hit) in is_hit.iter().enumerate() {
        if hit {
            running += hit_inc;
        } else {
            running -= miss_dec;
        }
        rs.push(running);

        if running.abs() > peak_val.abs() {
            peak_val = running;
            peak_idx = i;
        }
    }

    (peak_val, rs, peak_idx)
}

// ── Permutation null distribution ────────────────────────────────────────────

/// Compute `n_perm` permuted ES values by shuffling the hit indicator vector.
///
/// Uses a fast Fisher-Yates shuffle seeded per-permutation from the permutation
/// index to avoid requiring a shared RNG (no rayon-unsafe state).
fn permutation_null(is_hit: &[bool], n_hits: usize, n_perm: usize) -> Vec<f64> {
    (0..n_perm)
        .into_par_iter()
        .map(|perm_idx| {
            // Shuffle hit positions using a simple deterministic PRNG.
            // We scatter the n_hits hit positions uniformly across the n positions.
            let n = is_hit.len();
            let seed = (perm_idx as u64).wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let shuffled_positions = xorshift_sample(n, n_hits, seed);

            let mut perm_hit = vec![false; n];
            for pos in &shuffled_positions {
                perm_hit[*pos] = true;
            }

            let (es_perm, _, _) = compute_es(&perm_hit, n_hits);
            es_perm
        })
        .collect()
}

/// Sample `k` unique positions from [0, n) using Knuth reservoir sampling
/// driven by xorshift64 with the given seed.
fn xorshift_sample(n: usize, k: usize, seed: u64) -> Vec<usize> {
    let mut state = if seed == 0 { 1 } else { seed };

    let xorshift = |s: &mut u64| -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    };

    // Build a shuffled index array via partial Fisher-Yates
    let mut indices: Vec<usize> = (0..n).collect();
    for i in 0..k {
        let rand_val = xorshift(&mut state);
        let j = i + (rand_val as usize % (n - i));
        indices.swap(i, j);
    }
    indices[..k].to_vec()
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Run GSEA pre-ranked enrichment analysis.
///
/// # Arguments
/// - `ranked_genes`: genes sorted **descending** by enrichment metric
///   (e.g. sign(log2FC) × −log10(p)).
/// - `pathways`: list of `(id, name, genes)` tuples.
/// - `min_size`: minimum number of pathway genes present in the ranked list.
/// - `max_size`: maximum number of pathway genes present in the ranked list.
/// - `n_perm`: number of permutations for null ES estimation. Must be ≥ 1.
///   Use ≥ 1000 for valid p-values; ≥ 10 000 for publication quality.
///
/// # Returns
/// A `Vec<GseaResult>` sorted by NES descending (most positively enriched first).
/// BH FDR correction is applied across all returned results.
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

    let n_perm = n_perm.max(1);
    let n_genes = ranked_genes.len();

    // Pre-build the ranked gene lookup for quick membership testing
    // Compute results for each eligible pathway
    let mut results: Vec<GseaResult> = pathways
        .iter()
        .filter_map(|&(pid, pname, pgenes)| {
            // Build pathway gene set (uppercase for case-insensitive matching)
            let pathway_set: AHashSet<String> =
                pgenes.iter().map(|g| g.to_uppercase()).collect();

            // Count how many pathway genes appear in the ranked list
            let n_hits = ranked_genes
                .iter()
                .filter(|(g, _)| pathway_set.contains(g.as_str()))
                .count();

            if n_hits < min_size || n_hits > max_size {
                return None;
            }

            // Build hit indicator vector (true = gene at rank i is in pathway)
            let is_hit: Vec<bool> = ranked_genes
                .iter()
                .map(|(g, _)| pathway_set.contains(g.as_str()))
                .collect();

            // Observed enrichment score
            let (es, running_sum, peak_idx) = compute_es(&is_hit, n_hits);

            // Permutation null distribution
            let null_es = permutation_null(&is_hit, n_hits, n_perm);

            // NES: normalize by mean of absolute null ES values
            let null_abs_mean = {
                let s: f64 = null_es.iter().map(|x| x.abs()).sum();
                let m = s / null_es.len() as f64;
                if m < 1e-9 { 1.0 } else { m }
            };
            let nes = es / null_abs_mean;

            // Empirical p-value: fraction of permuted |ES| >= observed |ES|
            let obs_abs = es.abs();
            let p_value = {
                let extreme = null_es.iter().filter(|&&x| x.abs() >= obs_abs).count();
                // Use Laplace correction to avoid p=0 (prevents padj = 0)
                (extreme as f64 + 1.0) / (null_es.len() as f64 + 1.0)
            };

            // Leading edge: genes before (and including) the peak for positive ES,
            // or genes after (and including) the peak for negative ES.
            let leading_edge: Vec<String> = if es >= 0.0 {
                ranked_genes[..=peak_idx]
                    .iter()
                    .filter(|(g, _)| pathway_set.contains(g.as_str()))
                    .map(|(g, _)| g.clone())
                    .collect()
            } else {
                ranked_genes[peak_idx..]
                    .iter()
                    .filter(|(g, _)| pathway_set.contains(g.as_str()))
                    .map(|(g, _)| g.clone())
                    .collect()
            };

            Some(GseaResult {
                pathway_id: pid.to_string(),
                pathway_name: pname.to_string(),
                es,
                nes,
                p_value,
                padj: f64::NAN, // filled below
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

    // BH FDR correction
    let pvals: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    let padj_vals = benjamini_hochberg(&pvals);
    for (r, padj) in results.iter_mut().zip(padj_vals) {
        r.padj = padj;
    }

    // Sort by NES descending (most positively enriched first)
    results.sort_unstable_by(|a, b| {
        b.nes.partial_cmp(&a.nes).unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a ranked gene list of length `n`. The first `k` genes form a
    /// "perfect pathway" and receive the highest metric scores.
    fn make_ranked(n: usize) -> Vec<(String, f64)> {
        (0..n)
            .map(|i| {
                let id = format!("GENE{i:04}");
                // Linearly decreasing metric: highest rank at i=0
                let metric = (n as f64) - (i as f64);
                (id, metric)
            })
            .collect()
    }

    #[test]
    fn test_perfect_pathway_at_top_is_enriched() {
        // 200-gene ranked list; pathway = first 20 genes (all at top)
        let n = 200;
        let ranked = make_ranked(n);
        let pathway_genes: Vec<&str> = ranked[..20].iter().map(|(g, _)| g.as_str()).collect();

        let pathways: Vec<(&str, &str, &[&str])> =
            vec![("PW_TOP", "Top-loaded pathway", pathway_genes.as_slice())];

        let results = gsea_preranked(&ranked, &pathways, 5, 500, 100);

        assert_eq!(results.len(), 1, "should return exactly one result");
        let r = &results[0];

        assert!(r.es > 0.0, "ES should be positive for pathway at top: ES={}", r.es);
        assert!(r.nes > 1.0, "NES should exceed 1.0 for a strongly enriched pathway: NES={}", r.nes);
        assert!(r.p_value <= 0.10, "p-value should be low for perfect pathway: p={}", r.p_value);
        assert!(!r.leading_edge.is_empty(), "leading edge should be non-empty");
    }

    #[test]
    fn test_uniform_pathway_has_lower_es_than_enriched() {
        // A pathway uniformly spaced through the ranked list should produce a
        // substantially lower ES than one perfectly concentrated at the top.
        let n = 200_usize;
        let ranked = make_ranked(n);

        // Uniform pathway: every 10th gene (20 genes spread throughout)
        let uniform_genes: Vec<&str> =
            ranked.iter().step_by(10).map(|(g, _)| g.as_str()).collect();
        let uniform_pathways: Vec<(&str, &str, &[&str])> =
            vec![("PW_UNIFORM", "Uniform pathway", uniform_genes.as_slice())];

        // Enriched pathway: all 20 genes at the top
        let enriched_genes: Vec<&str> =
            ranked[..20].iter().map(|(g, _)| g.as_str()).collect();
        let enriched_pathways: Vec<(&str, &str, &[&str])> =
            vec![("PW_TOP", "Top pathway", enriched_genes.as_slice())];

        let uniform_res = gsea_preranked(&ranked, &uniform_pathways, 5, 500, 200);
        let enriched_res = gsea_preranked(&ranked, &enriched_pathways, 5, 500, 200);

        assert_eq!(uniform_res.len(), 1);
        assert_eq!(enriched_res.len(), 1);

        // The enriched pathway must score higher than the uniform one
        assert!(
            enriched_res[0].es > uniform_res[0].es,
            "top-loaded pathway ES ({}) should exceed uniform pathway ES ({})",
            enriched_res[0].es, uniform_res[0].es,
        );
    }

    #[test]
    fn test_pathway_filtered_by_size_bounds() {
        let ranked = make_ranked(100);
        // A pathway with only 2 genes should be filtered out with min_size=5
        let tiny_genes: Vec<&str> = ranked[..2].iter().map(|(g, _)| g.as_str()).collect();
        let pathways = vec![("TINY", "Tiny pathway", tiny_genes.as_slice())];
        let results = gsea_preranked(&ranked, &pathways, 5, 500, 10);
        assert!(results.is_empty(), "pathway below min_size should be filtered");
    }

    #[test]
    fn test_running_sum_length_matches_ranked() {
        let n = 50;
        let ranked = make_ranked(n);
        let pathway_genes: Vec<&str> = ranked[..10].iter().map(|(g, _)| g.as_str()).collect();
        let pathways = vec![("PW", "Test", pathway_genes.as_slice())];
        let results = gsea_preranked(&ranked, &pathways, 5, 500, 50);
        assert_eq!(results[0].running_sum.len(), n, "running sum length should equal ranked list length");
    }
}
