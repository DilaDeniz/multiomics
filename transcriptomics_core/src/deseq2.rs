//! DESeq2-style size-factor normalization and differential expression for raw
//! integer count matrices.
//!
//! Size factors are estimated by the median-of-ratios method (Anders & Huber 2010).
//! Differential expression uses Welch t-test on log₂(normalized + 0.5) followed
//! by Benjamini-Hochberg FDR correction.

use anyhow::{bail, Context, Result};
use biomics_core::statistics::{benjamini_hochberg, welch_t_test};

use crate::types::DiffExprResult;

// ── Public data structures ────────────────────────────────────────────────────

/// Per-sample size factors estimated by the median-of-ratios method.
#[derive(Debug, Clone)]
pub struct SizeFactors {
    pub sample_names: Vec<String>,
    /// One factor per sample, same order as `sample_names`.
    pub factors: Vec<f64>,
}

/// A full normalized count matrix together with the size factors and the
/// pooled dispersion estimate.
#[derive(Debug, Clone)]
pub struct NormalizedMatrix {
    pub gene_ids: Vec<String>,
    pub sample_names: Vec<String>,
    pub size_factors: SizeFactors,
    /// Raw counts — `counts[gene_index][sample_index]`.
    pub counts: Vec<Vec<f64>>,
    /// Normalized counts — `normalized[gene_index][sample_index]`.
    pub normalized: Vec<Vec<f64>>,
    /// Pooled (median) negative-binomial dispersion across all genes.
    pub global_dispersion: f64,
}

// ── Core estimation functions ─────────────────────────────────────────────────

/// Estimate per-sample size factors using the DESeq2 median-of-ratios method.
///
/// `counts[gene][sample]` — raw integer counts as f64.
pub fn estimate_size_factors(counts: &[Vec<f64>], sample_names: &[String]) -> Result<SizeFactors> {
    let n_genes = counts.len();
    if n_genes == 0 {
        bail!("estimate_size_factors: count matrix is empty");
    }
    let n_samples = sample_names.len();
    if n_samples == 0 {
        bail!("estimate_size_factors: no sample names provided");
    }
    // Validate row lengths
    for (g, row) in counts.iter().enumerate() {
        if row.len() != n_samples {
            bail!(
                "estimate_size_factors: gene {} has {} values but {} samples declared",
                g,
                row.len(),
                n_samples
            );
        }
    }

    // Step 1: geometric mean per gene in log-space.
    // Only include non-zero counts to avoid log(0).
    // Genes where every sample is zero are excluded from size-factor estimation
    // (their geometric mean would be 0, making ratios undefined).
    let geom_means: Vec<f64> = counts
        .iter()
        .map(|row| {
            let positive_logs: Vec<f64> =
                row.iter().filter(|&&c| c > 0.0).map(|&c| c.ln()).collect();
            if positive_logs.is_empty() {
                0.0 // sentinel: excluded below
            } else {
                let mean_log = positive_logs.iter().sum::<f64>() / positive_logs.len() as f64;
                mean_log.exp()
            }
        })
        .collect();

    // Step 2: for each sample, collect ratios count_gj / mu_g for genes where
    // mu_g > 0, then take the median.
    let mut factors = Vec::with_capacity(n_samples);
    for j in 0..n_samples {
        let mut ratios: Vec<f64> = counts
            .iter()
            .zip(geom_means.iter())
            .filter(|(_, &mu)| mu > 0.0)
            .map(|(row, &mu)| row[j] / mu)
            .collect();

        if ratios.is_empty() {
            bail!(
                "estimate_size_factors: no valid genes for size-factor estimation of sample {}",
                sample_names[j]
            );
        }

        ratios.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median = median_sorted(&ratios);

        // Guard against degenerate size factors
        if median <= 0.0 || !median.is_finite() {
            bail!(
                "estimate_size_factors: computed non-positive size factor ({}) for sample {}",
                median,
                sample_names[j]
            );
        }
        factors.push(median);
    }

    Ok(SizeFactors {
        sample_names: sample_names.to_vec(),
        factors,
    })
}

/// Normalize a count matrix and estimate global dispersion.
///
/// `counts[gene][sample]` — raw integer counts as f64.
pub fn normalize_counts(
    gene_ids: &[String],
    counts: &[Vec<f64>],
    sample_names: &[String],
) -> Result<NormalizedMatrix> {
    if gene_ids.len() != counts.len() {
        bail!(
            "normalize_counts: {} gene IDs but {} count rows",
            gene_ids.len(),
            counts.len()
        );
    }

    let size_factors = estimate_size_factors(counts, sample_names)
        .context("normalize_counts: size-factor estimation failed")?;

    let n_samples = sample_names.len();

    // Normalize: normalized_gj = count_gj / s_j
    let normalized: Vec<Vec<f64>> = counts
        .iter()
        .map(|row| {
            row.iter()
                .zip(size_factors.factors.iter())
                .map(|(&c, &s)| c / s)
                .collect()
        })
        .collect();

    // Dispersion estimation via method of moments.
    // For each gene: mean_g and var_g of normalized counts across samples.
    // alpha_g = max(0, (var_g - mean_g) / mean_g^2)   [NB dispersion]
    let mut dispersions: Vec<f64> = Vec::with_capacity(normalized.len());

    for norm_row in &normalized {
        if norm_row.is_empty() {
            continue;
        }
        let n = norm_row.len() as f64;
        let mean_g = norm_row.iter().sum::<f64>() / n;
        if mean_g <= 0.0 {
            dispersions.push(0.0);
            continue;
        }
        let var_g = if n_samples > 1 {
            norm_row.iter().map(|&x| (x - mean_g).powi(2)).sum::<f64>() / (n - 1.0)
        } else {
            0.0
        };
        let alpha = ((var_g - mean_g) / mean_g.powi(2)).max(0.0);
        dispersions.push(alpha);
    }

    // Global dispersion: median across all genes
    dispersions.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let global_dispersion = if dispersions.is_empty() {
        0.0
    } else {
        median_sorted(&dispersions)
    };

    Ok(NormalizedMatrix {
        gene_ids: gene_ids.to_vec(),
        sample_names: sample_names.to_vec(),
        size_factors,
        counts: counts.to_vec(),
        normalized,
        global_dispersion,
    })
}

/// Differential expression on a `NormalizedMatrix`.
///
/// Groups are split by position — first ⌈n/2⌉ samples vs. the remainder —
/// matching the convention in `diffexpr.rs`. Welch's t-test is applied on
/// log₂(normalized + 0.5). BH FDR is applied across genes with a valid p-value.
/// Results are sorted by ascending `padj` (NaN last), then descending |log2FC|.
pub fn deseq2_differential_expression(matrix: &NormalizedMatrix) -> Vec<DiffExprResult> {
    let n_samples = matrix.sample_names.len();
    if n_samples == 0 || matrix.normalized.is_empty() {
        return Vec::new();
    }

    let split = n_samples.div_ceil(2);
    let can_test = split >= 2 && (n_samples - split) >= 2;

    let mut results: Vec<DiffExprResult> = matrix
        .normalized
        .iter()
        .zip(matrix.gene_ids.iter())
        .map(|(norm_row, gene_id)| {
            // log₂(normalized + 0.5) to handle zeros and approach the log-normal
            let g1: Vec<f64> = norm_row[..split]
                .iter()
                .map(|&v| (v + 0.5_f64).log2())
                .collect();
            let g2: Vec<f64> = norm_row[split..]
                .iter()
                .map(|&v| (v + 0.5_f64).log2())
                .collect();

            let mean1 = g1.iter().sum::<f64>() / g1.len().max(1) as f64;
            let mean2 = g2.iter().sum::<f64>() / g2.len().max(1) as f64;
            let lfc = mean2 - mean1; // log2 scale → log2FC

            // Raw (non-log) means for reporting
            let mean_s1 = norm_row[..split].iter().sum::<f64>() / split.max(1) as f64;
            let mean_s2 = norm_row[split..].iter().sum::<f64>() / (n_samples - split).max(1) as f64;

            let (p_value, padj) = if can_test {
                let pval = welch_t_test(&g1, &g2).map(|(_, p)| p).unwrap_or(f64::NAN);
                (pval, f64::NAN) // padj filled below
            } else {
                (f64::NAN, f64::NAN)
            };

            DiffExprResult {
                gene_id: gene_id.clone(),
                log2_fold_change: lfc,
                mean_s1,
                mean_s2,
                p_value,
                padj,
            }
        })
        .collect();

    // BH FDR correction over the subset with valid p-values
    if can_test {
        let valid_indices: Vec<usize> = results
            .iter()
            .enumerate()
            .filter(|(_, r)| !r.p_value.is_nan())
            .map(|(i, _)| i)
            .collect();

        if !valid_indices.is_empty() {
            let pvals: Vec<f64> = valid_indices.iter().map(|&i| results[i].p_value).collect();
            let padj_vals = benjamini_hochberg(&pvals);
            for (vi, &orig_i) in valid_indices.iter().enumerate() {
                results[orig_i].padj = padj_vals[vi];
            }
        }
    }

    // Sort: smallest padj first (NaN last), break ties by descending |log2FC|
    results.sort_unstable_by(|a, b| match (a.padj.is_nan(), b.padj.is_nan()) {
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => a
            .padj
            .partial_cmp(&b.padj)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(
                b.log2_fold_change
                    .abs()
                    .partial_cmp(&a.log2_fold_change.abs())
                    .unwrap_or(std::cmp::Ordering::Equal),
            ),
    });

    results
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Median of a **sorted** slice (no allocation).
/// Panics if the slice is empty — callers must guard.
fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    debug_assert!(n > 0, "median_sorted called on empty slice");
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// A small 4-gene × 4-sample count matrix with known geometric structure.
    ///
    /// Sample library sizes: [100, 200, 150, 300]. After scaling each column by
    /// its total relative to a reference (sample 0), expected size factors are
    /// approximately [1.0, 2.0, 1.5, 3.0]. The median-of-ratios method produces
    /// values close to this but not identical; we allow ±20 %.
    fn make_test_counts() -> (Vec<String>, Vec<Vec<f64>>, Vec<String>) {
        // 6 genes, 4 samples, with a consistent 2× library-size difference
        // between (sample0, sample1) and (sample2, sample3) pairs
        let gene_ids: Vec<String> = (0..6).map(|i| format!("gene{i}")).collect();
        let sample_names: Vec<String> = vec![
            "s1".to_string(),
            "s2".to_string(),
            "s3".to_string(),
            "s4".to_string(),
        ];
        // counts[gene][sample]
        // s2 has 2× the counts of s1; s3 has 1.5×; s4 has 3×
        let counts: Vec<Vec<f64>> = vec![
            vec![10.0, 20.0, 15.0, 30.0],
            vec![20.0, 40.0, 30.0, 60.0],
            vec![30.0, 60.0, 45.0, 90.0],
            vec![40.0, 80.0, 60.0, 120.0],
            vec![50.0, 100.0, 75.0, 150.0],
            vec![60.0, 120.0, 90.0, 180.0],
        ];
        (gene_ids, counts, sample_names)
    }

    #[test]
    fn test_size_factors() {
        let (_, counts, sample_names) = make_test_counts();
        let sf = estimate_size_factors(&counts, &sample_names).expect("size factors");

        // Expected relative sizes: s1=1.0, s2=2.0, s3=1.5, s4=3.0
        // The median-of-ratios method anchors to geometric-mean rows,
        // so absolute values may shift, but the *ratios* between samples
        // must reflect the actual library-size differences.
        let f = &sf.factors;
        assert_eq!(f.len(), 4);

        // All factors must be positive and finite
        for (name, &factor) in sample_names.iter().zip(f.iter()) {
            assert!(
                factor > 0.0 && factor.is_finite(),
                "size factor for {name} is {factor}"
            );
        }

        // Ratios: f[1]/f[0] ≈ 2.0, f[2]/f[0] ≈ 1.5, f[3]/f[0] ≈ 3.0
        let tol = 0.20; // 20 % relative tolerance
        let ratio_12 = f[1] / f[0];
        assert!(
            (ratio_12 - 2.0).abs() < 2.0 * tol,
            "ratio s2/s1 = {ratio_12} not within 20% of 2.0"
        );
        let ratio_13 = f[2] / f[0];
        assert!(
            (ratio_13 - 1.5).abs() < 1.5 * tol,
            "ratio s3/s1 = {ratio_13} not within 20% of 1.5"
        );
        let ratio_14 = f[3] / f[0];
        assert!(
            (ratio_14 - 3.0).abs() < 3.0 * tol,
            "ratio s4/s1 = {ratio_14} not within 20% of 3.0"
        );
    }

    #[test]
    fn test_normalization() {
        let (gene_ids, counts, sample_names) = make_test_counts();
        let matrix = normalize_counts(&gene_ids, &counts, &sample_names).expect("normalize_counts");

        // After normalization, the per-sample sums of normalized counts should
        // be roughly equal (within 20 %) across samples, because the size
        // factors absorb the library-size differences.
        let col_sums: Vec<f64> = (0..sample_names.len())
            .map(|j| matrix.normalized.iter().map(|row| row[j]).sum::<f64>())
            .collect();

        let mean_sum = col_sums.iter().sum::<f64>() / col_sums.len() as f64;
        for (name, &s) in sample_names.iter().zip(col_sums.iter()) {
            let rel_dev = (s - mean_sum).abs() / mean_sum;
            assert!(
                rel_dev < 0.20,
                "sample {name}: normalized library size {s:.1} deviates {:.1}% from mean {mean_sum:.1}",
                rel_dev * 100.0
            );
        }

        // Global dispersion must be non-negative and finite
        assert!(
            matrix.global_dispersion >= 0.0 && matrix.global_dispersion.is_finite(),
            "global_dispersion = {}",
            matrix.global_dispersion
        );
    }

    #[test]
    fn test_differential_expression_no_panic_small() {
        // 2 groups of 2 → can_test = true (split=2, remainder=2)
        // Design: genes 0-3 all have ~10× higher raw counts in group 2 (b1,b2)
        // vs group 1 (a1,a2). After DESeq2 normalization the size factors will
        // absorb that global library-size difference, so per-gene LFC will be
        // near zero (all genes change proportionally). What we DO assert:
        //   - all p-values and padj are finite (no NaN),
        //   - results are sorted correctly (ascending padj),
        //   - there is at least one result per gene.
        let gene_ids: Vec<String> = (0..4).map(|i| format!("g{i}")).collect();
        let sample_names: Vec<String> = vec![
            "a1".to_string(),
            "a2".to_string(),
            "b1".to_string(),
            "b2".to_string(),
        ];
        // group 1 (a1,a2) has low raw counts; group 2 (b1,b2) has high raw counts.
        // After median-of-ratios normalization the library-size difference is
        // corrected, leaving only within-gene variance — so LFC ≈ 0 for all genes.
        let counts: Vec<Vec<f64>> = vec![
            vec![10.0, 12.0, 100.0, 110.0],
            vec![8.0, 9.0, 80.0, 90.0],
            vec![5.0, 6.0, 50.0, 55.0],
            vec![20.0, 22.0, 200.0, 220.0],
        ];
        let matrix = normalize_counts(&gene_ids, &counts, &sample_names).expect("normalize_counts");
        let de = deseq2_differential_expression(&matrix);
        assert_eq!(de.len(), 4);

        // All results must have finite p-values (4 samples → 2+2 split → can_test)
        for r in &de {
            assert!(r.p_value.is_finite(), "p_value is NaN for {}", r.gene_id);
            assert!(r.padj.is_finite(), "padj is NaN for {}", r.gene_id);
        }

        // After normalization, all genes are proportionally identical across groups,
        // so |log2FC| must be small (< 0.5 after normalization absorbs library size).
        for r in &de {
            assert!(
                r.log2_fold_change.abs() < 0.5,
                "expected |log2FC| < 0.5 after normalization for {}, got {}",
                r.gene_id,
                r.log2_fold_change
            );
        }

        // Test with a count matrix that has genuine differential expression:
        // gene A is truly up in group 2 (beyond library-size correction).
        let gene_ids2: Vec<String> = vec!["de_gene".to_string(), "null_gene".to_string()];
        let sample_names2: Vec<String> = vec![
            "x1".to_string(),
            "x2".to_string(),
            "y1".to_string(),
            "y2".to_string(),
        ];
        // de_gene: group 2 is 16× higher AFTER adjusting for 2× overall library size
        // null_gene: group 2 is 2× higher (pure library-size effect, corrected to 1×)
        let counts2: Vec<Vec<f64>> = vec![
            vec![10.0, 11.0, 160.0, 175.0], // de_gene: up ~8× after normalization
            vec![20.0, 22.0, 40.0, 44.0],   // null_gene: pure library-size, LFC ≈ 0
        ];
        let matrix2 = normalize_counts(&gene_ids2, &counts2, &sample_names2).unwrap();
        let de2 = deseq2_differential_expression(&matrix2);

        // de_gene should have a clearly positive LFC after normalization.
        // With pseudocount 0.5, log2(norm+0.5) attenuates the true 8× fold change,
        // so we use a conservative threshold of > 1.0.
        let de_gene_result = de2.iter().find(|r| r.gene_id == "de_gene").unwrap();
        assert!(
            de_gene_result.log2_fold_change > 1.0,
            "de_gene should have log2FC > 1 after normalization, got {}",
            de_gene_result.log2_fold_change
        );
    }
}
