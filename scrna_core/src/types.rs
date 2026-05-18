//! Shared summary types for the single-cell pipeline.

use crate::de::ClusterMarker;

/// High-level summary produced by [`crate::run_scrna_pipeline`].
#[derive(Debug, Clone, serde::Serialize)]
pub struct SingleCellSummary {
    /// Number of cells before QC filtering.
    pub n_cells_raw: usize,
    /// Number of cells retained after QC filtering.
    pub n_cells_after_qc: usize,
    /// Total number of features (genes) in the dataset.
    pub n_genes: usize,
    /// Number of highly variable genes selected.
    pub n_hvg: usize,
    /// Number of Leiden clusters identified.
    pub n_clusters: u32,
    /// Median number of detected genes per cell (post-QC).
    pub median_genes_per_cell: f64,
    /// Median total UMI counts per cell (post-QC).
    pub median_counts_per_cell: f64,
    /// Top 3 marker genes per cluster.
    pub top_markers: Vec<ClusterMarker>,
}
