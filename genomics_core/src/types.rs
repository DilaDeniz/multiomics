use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::cancer::{HrdScore, KataegisLocus, LohChromosome, MsiResult, TmbResult, TumorPurityResult};

/// Whether a SNP is a transition or transversion, or an indel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TiTvClass {
    /// Purine↔purine or pyrimidine↔pyrimidine substitution.
    Transition,
    /// Purine↔pyrimidine substitution.
    Transversion,
    /// Insertion or deletion.
    Indel,
    /// Multi-allelic or ambiguous.
    Other,
}

/// A single parsed variant record, cheaply cloneable.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantRecord {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    /// Phred-scaled quality score from QUAL column.
    pub qual: f32,
    pub titv: TiTvClass,
    /// Allele frequency from INFO/AF field, if present.
    pub af: Option<f32>,
    /// Gene name from INFO/GENE or ANN field, if present.
    pub gene: Option<String>,
}

/// Per-chromosome variant counts.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ChromDensity {
    pub total: u64,
    pub snps: u64,
    pub indels: u64,
}

/// Final summary produced by `GenomicsAccum::finalize`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GenomicsSummary {
    pub total_variants: u64,
    pub snp_count: u64,
    pub indel_count: u64,
    /// Ratio of transitions to transversions (expected ~2.0–2.1 for WGS).
    pub titv_ratio: f64,
    pub per_chrom: HashMap<String, ChromDensity>,
    /// Variants with QUAL > 30.
    pub high_impact: Vec<VariantRecord>,
    /// Allele frequency histogram: 20 bins over [0.0, 1.0).
    pub af_histogram: Vec<u64>,
    /// Approximate count of unique (chrom, pos) positions.
    pub unique_positions: u64,
    /// All gene names appearing in high-impact variants.
    pub high_impact_genes: Vec<String>,
    /// Tumor purity estimate (None until set by integration layer).
    #[serde(default)]
    pub tumor_purity: Option<TumorPurityResult>,
    /// Kataegis (hypermutation) loci detected from all variants.
    #[serde(default)]
    pub kataegis_loci: Vec<KataegisLocus>,
    /// Homologous recombination deficiency score from indel spectrum.
    #[serde(default)]
    pub hrd: Option<HrdScore>,
    /// Per-chromosome loss of heterozygosity assessment.
    #[serde(default)]
    pub loh_chromosomes: Vec<LohChromosome>,
    /// Tumor mutational burden (filled in by runner after context detection).
    #[serde(default)]
    pub tmb: Option<TmbResult>,
    /// Microsatellite instability score derived from homopolymer indel fraction.
    #[serde(default)]
    pub msi: Option<MsiResult>,
}
