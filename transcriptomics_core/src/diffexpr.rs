use biomics_core::statistics::{benjamini_hochberg, welch_t_test};

use crate::types::{DiffExprResult, GeneRecord};

/// Compute differential expression for all genes using Welch's t-test + BH FDR.
///
/// ## Grouping strategy
/// Samples are split by position: **group 1** = first ⌈n/2⌉ columns,
/// **group 2** = remaining columns. When n < 4, a t-test cannot be computed
/// (at least one group has fewer than 2 observations); in that case `p_value`
/// and `padj` are set to `f64::NAN` and only `log2_fold_change` is reported.
///
/// ## Log₂ transform
/// All TPM values are transformed as log₂(TPM + 0.5) before testing. The
/// pseudocount of 0.5 prevents –∞ for zero-expression genes while being
/// smaller than typical expressed-gene TPM levels (≥ 1).
///
/// ## Output
/// Results are BH-corrected across all genes with a valid p-value, then
/// sorted by ascending `padj` (NaN entries last).
pub fn differential_expression(records: &[GeneRecord]) -> Vec<DiffExprResult> {
    if records.is_empty() {
        return Vec::new();
    }

    let n_samples = records.first().map(|r| r.samples.len()).unwrap_or(0);
    let split = n_samples.div_ceil(2); // ceil(n/2) → group 1 size
    let can_test = split >= 2 && (n_samples - split) >= 2;

    // Compute log₂FC and, if possible, raw p-values
    let mut results: Vec<DiffExprResult> = records
        .iter()
        .filter(|r| r.samples.len() == n_samples)
        .map(|r| {
            let g1: Vec<f64> = r.samples[..split].iter().map(|&v| (v + 0.5_f64).log2()).collect();
            let g2: Vec<f64> = r.samples[split..].iter().map(|&v| (v + 0.5_f64).log2()).collect();

            let mean1 = g1.iter().sum::<f64>() / g1.len().max(1) as f64;
            let mean2 = g2.iter().sum::<f64>() / g2.len().max(1) as f64;
            let lfc = mean2 - mean1; // log2 scale difference = log2FC

            let mean_s1 = r.samples[..split].iter().sum::<f64>() / split.max(1) as f64;
            let mean_s2 =
                r.samples[split..].iter().sum::<f64>() / (n_samples - split).max(1) as f64;

            let (p_value, padj) = if can_test {
                let pval = welch_t_test(&g1, &g2).map(|(_, p)| p).unwrap_or(f64::NAN);
                (pval, f64::NAN) // padj filled in below
            } else {
                (f64::NAN, f64::NAN)
            };

            DiffExprResult { gene_id: r.gene_id.clone(), log2_fold_change: lfc, mean_s1, mean_s2, p_value, padj }
        })
        .collect();

    // Apply BH FDR correction to the subset with valid p-values
    if can_test {
        let valid_indices: Vec<usize> =
            results.iter().enumerate().filter(|(_, r)| !r.p_value.is_nan()).map(|(i, _)| i).collect();

        if !valid_indices.is_empty() {
            let pvals: Vec<f64> = valid_indices.iter().map(|&i| results[i].p_value).collect();
            let padj_vals = benjamini_hochberg(&pvals);
            for (vi, &orig_i) in valid_indices.iter().enumerate() {
                results[orig_i].padj = padj_vals[vi];
            }
        }
    }

    // Sort: significant first (by padj, NaN last), then by |log2FC|
    results.sort_unstable_by(|a, b| {
        match (a.padj.is_nan(), b.padj.is_nan()) {
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            _ => a.padj.partial_cmp(&b.padj)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(b.log2_fold_change.abs().partial_cmp(&a.log2_fold_change.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)),
        }
    });

    results
}

/// Return gene IDs that pass significance thresholds.
///
/// When `padj` is available, the filter is `padj < 0.05 AND |log2FC| ≥ 1`.
/// When `padj` is `NaN` (n < 4 samples), falls back to `|log2FC| ≥ threshold`
/// with the caller-supplied `min_tpm` guard.
pub fn significant_de_genes(results: &[DiffExprResult], threshold: f64, min_tpm: f64) -> Vec<String> {
    results
        .iter()
        .filter(|r| {
            let expressed = r.mean_s1 >= min_tpm || r.mean_s2 >= min_tpm;
            let sig = if r.padj.is_nan() {
                r.log2_fold_change.abs() >= threshold
            } else {
                r.padj < 0.05 && r.log2_fold_change.abs() >= 1.0
            };
            expressed && sig
        })
        .map(|r| r.gene_id.clone())
        .collect()
}
