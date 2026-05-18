//! Public result type for GSEA analysis.

use serde::{Deserialize, Serialize};

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
    /// Empirical permutation p-value (multilevel adaptive estimate).
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
