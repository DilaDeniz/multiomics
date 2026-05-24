use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// One row from the expression matrix TSV.
///
/// The TSV header row provides sample names; subsequent rows are:
/// `gene_id\tsample1_tpm\tsample2_tpm\t...`
#[derive(Debug, Clone)]
pub struct GeneRecord {
    pub gene_id: String,
    /// TPM values, one per sample column, in header order.
    pub samples: Vec<f64>,
}

/// Per-gene descriptive statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneStats {
    pub mean: f64,
    pub std: f64,
    pub max: f64,
}

/// Differential expression result for one gene.
///
/// For n=2 samples, `p_value` and `padj` are `f64::NAN` (no within-group
/// variance to estimate). For n≥4 samples, Welch's t-test is used and BH
/// FDR correction is applied across all tested genes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiffExprResult {
    pub gene_id: String,
    /// log₂(mean_group2 / mean_group1) — positive = up in group 2.
    pub log2_fold_change: f64,
    pub mean_s1: f64,
    pub mean_s2: f64,
    /// Raw two-tailed Welch t-test p-value. NaN when n < 4.
    pub p_value: f64,
    /// Benjamini-Hochberg adjusted p-value. NaN when n < 4.
    pub padj: f64,
}

/// Final summary produced by `TranscriptomicsAccum::finalize`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TranscriptomicsSummary {
    pub total_genes: u64,
    /// Genes with mean TPM ≥ 1.0 across all samples.
    pub expressed_genes: u64,
    pub low_expression_genes: Vec<String>,
    pub gene_stats: HashMap<String, GeneStats>,
    /// Top 100 expressed genes by mean TPM, sorted descending.
    pub top_100_expressed: Vec<(String, f64)>,
    /// Differential expression results (populated when n_samples ≥ 2).
    pub diff_expr: Option<Vec<DiffExprResult>>,
    pub sample_count: usize,
    pub sample_names: Vec<String>,
}
