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

/// Two-sample differential expression result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiffExprResult {
    pub gene_id: String,
    /// log₂(sample2 / sample1) — positive = upregulated in sample2.
    pub log2_fold_change: f64,
    pub mean_s1: f64,
    pub mean_s2: f64,
}

/// Final summary produced by `TranscriptomicsAccum::finalize`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptomicsSummary {
    pub total_genes: u64,
    /// Genes with mean TPM ≥ 1.0 across all samples.
    pub expressed_genes: u64,
    pub low_expression_genes: Vec<String>,
    pub gene_stats: HashMap<String, GeneStats>,
    /// Top 100 expressed genes by mean TPM, sorted descending.
    pub top_100_expressed: Vec<(String, f64)>,
    /// Differential expression results (populated only when exactly 2 samples).
    pub diff_expr: Option<Vec<DiffExprResult>>,
    pub sample_count: usize,
    pub sample_names: Vec<String>,
}
