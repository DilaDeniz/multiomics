use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// One record from a 6-column BED methylation file.
///
/// Standard ENCODE bisulfite BED format:
/// `chrom  start  end  name  score  strand`
/// where `score` is methylation percentage (0–1000, divide by 10 for %).
#[derive(Debug, Clone)]
pub struct MethylationRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// Methylation percentage in [0.0, 100.0].
    pub methylation: f64,
    /// Gene name from the name column (col4) if present and non-numeric.
    pub gene: Option<String>,
}

/// A detected CpG island.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CpGIsland {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub gc_percent: f64,
    pub mean_methylation: f64,
    /// Observed-to-expected CpG ratio (CpO/E). Values ≥ 0.6 define CpG islands
    /// under the Gardiner-Garden & Frommer (1987) criterion.
    pub cpoe: f64,
}

/// Per-chromosome methylation statistics.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ChromMethylation {
    pub total_sites: u64,
    pub sum_methylation: f64,
    /// Filled at finalize time.
    pub mean_methylation: f64,
}

/// Classification of a methylation region.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum RegionKind {
    /// Mean methylation > 80%.
    Hypermethylated,
    /// Mean methylation < 20%.
    Hypomethylated,
}

/// A genomic region with extreme methylation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MethylationRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub mean_methylation: f64,
    pub kind: RegionKind,
}

/// Final summary produced by `EpigenomicsAccum::finalize`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct EpigenomicsSummary {
    pub total_sites: u64,
    /// Mean methylation across all sites.
    pub global_methylation_pct: f64,
    pub per_chrom: HashMap<String, ChromMethylation>,
    pub cpg_islands: Vec<CpGIsland>,
    pub hypermethylated: Vec<MethylationRegion>,
    pub hypomethylated: Vec<MethylationRegion>,
    /// Mean methylation per gene, populated when the BED file contains gene names in col4.
    /// Empty when no gene annotations are present in the BED file.
    pub gene_methylation: HashMap<String, f64>,
}
