use crate::types::{DiffExprResult, GeneRecord};
use biomics_core::stats::log2_fold_change;

/// Compute log₂ fold-change between sample index 0 and sample index 1.
///
/// Only meaningful when `records` come from a 2-sample expression matrix.
/// Genes where either sample has TPM = 0 are excluded (undefined fold change).
///
/// Results are sorted by absolute fold change, descending.
pub fn two_sample_diff_expr(records: &[GeneRecord]) -> Vec<DiffExprResult> {
    let mut results: Vec<DiffExprResult> = records
        .iter()
        .filter(|r| r.samples.len() >= 2)
        .filter_map(|r| {
            let s1 = r.samples[0];
            let s2 = r.samples[1];
            let lfc = log2_fold_change(s1, s2)?;
            Some(DiffExprResult {
                gene_id: r.gene_id.clone(),
                log2_fold_change: lfc,
                mean_s1: s1,
                mean_s2: s2,
            })
        })
        .collect();

    results.sort_unstable_by(|a, b| {
        b.log2_fold_change
            .abs()
            .partial_cmp(&a.log2_fold_change.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

/// Return gene IDs with |log₂FC| ≥ `threshold` and mean expression ≥ `min_tpm`
/// in at least one sample.
pub fn significant_de_genes(results: &[DiffExprResult], threshold: f64, min_tpm: f64) -> Vec<String> {
    results
        .iter()
        .filter(|r| {
            r.log2_fold_change.abs() >= threshold
                && (r.mean_s1 >= min_tpm || r.mean_s2 >= min_tpm)
        })
        .map(|r| r.gene_id.clone())
        .collect()
}
