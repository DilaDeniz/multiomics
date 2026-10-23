use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::cancer::{
    HrdScore, KataegisLocus, LohChromosome, MsiResult, TmbResult, TumorPurityResult,
};
use crate::cosmic::MutationalSignatureResult;
use crate::prs::PrsResult;

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

/// Variant class of one ALT allele relative to REF.
///
/// Mirrors htslib's `bcf_set_variant_type`, i.e. the classification behind
/// `bcftools stats`, so per-type counts can be validated against bcftools.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum VariantKind {
    /// Single-nucleotide variant (after trimming shared flanking bases).
    #[default]
    Snv,
    /// Multi-nucleotide polymorphism: equal-length REF/ALT differing at >1 base.
    Mnp,
    /// Insertion or deletion.
    Indel,
    /// Symbolic allele (`<DEL>`, `<DUP>` …) or complex substitution.
    Other,
    /// Breakend (`G]chr2:5000]`, `[chr2:5000[G` …).
    Breakend,
}

/// Putative functional impact from SnpEff (`ANN`, field 3) or Ensembl VEP
/// (`CSQ`, `IMPACT` field). Ordered from least to most severe.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Impact {
    Modifier,
    Low,
    Moderate,
    High,
}

impl Impact {
    pub fn parse(s: &[u8]) -> Option<Self> {
        match s {
            b"HIGH" => Some(Impact::High),
            b"MODERATE" => Some(Impact::Moderate),
            b"LOW" => Some(Impact::Low),
            b"MODIFIER" => Some(Impact::Modifier),
            _ => None,
        }
    }
}

/// Bit flags describing a whole VCF record (site), stored on the first ALT
/// allele emitted for that record and zero on the others.
pub mod site_flags {
    /// The record is the first allele of its site.
    pub const FIRST: u8 = 1;
    pub const SNV: u8 = 1 << 1;
    pub const MNP: u8 = 1 << 2;
    pub const INDEL: u8 = 1 << 3;
    pub const OTHER: u8 = 1 << 4;
    /// The record had more than one ALT allele (including `<NON_REF>`/`*`).
    pub const MULTIALLELIC: u8 = 1 << 5;
}

/// A single biallelic variant, cheaply cloneable.
///
/// Multi-allelic VCF records are split into one `VariantRecord` per ALT allele
/// at parse time (like `bcftools norm -m-`), so downstream analyses can rely
/// on `alt_allele` holding exactly one allele.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantRecord {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    /// Phred-scaled quality score from QUAL column.
    pub qual: f32,
    /// Transition/transversion class for SNVs; `Indel`/`Other` otherwise.
    pub titv: TiTvClass,
    /// Allele frequency from INFO/AF field, if present.
    pub af: Option<f32>,
    /// Gene name from INFO/GENE, SnpEff ANN or VEP CSQ, if present.
    pub gene: Option<String>,
    /// Most severe SnpEff/VEP impact annotated for this ALT allele.
    #[serde(default)]
    pub impact: Option<Impact>,
    /// Copies of this ALT allele in the first sample's GT (0–2 for diploid
    /// calls); `None` for sites-only VCFs or missing genotypes.
    #[serde(default)]
    pub dosage: Option<u8>,
    /// Variant class of this ALT allele.
    #[serde(default)]
    pub kind: VariantKind,
    /// [`site_flags`] bitmask; non-zero only on the first allele of a record.
    #[serde(default)]
    pub site_flags: u8,
}

impl VariantRecord {
    /// Build a standalone biallelic record, classifying the allele pair.
    pub fn biallelic(
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
        qual: f32,
        af: Option<f32>,
        gene: Option<&str>,
    ) -> Self {
        let (kind, titv) =
            match crate::vcf::classify_allele(ref_allele.as_bytes(), alt_allele.as_bytes()) {
                crate::vcf::AlleleClass::Variant(kind, titv) => (kind, titv),
                crate::vcf::AlleleClass::NotVariant => (VariantKind::Other, TiTvClass::Other),
            };
        let flag = match kind {
            VariantKind::Snv => site_flags::SNV,
            VariantKind::Mnp => site_flags::MNP,
            VariantKind::Indel => site_flags::INDEL,
            VariantKind::Other => site_flags::OTHER,
            VariantKind::Breakend => 0,
        };
        Self {
            chrom: chrom.to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            qual,
            titv,
            af,
            gene: gene.map(str::to_string),
            impact: None,
            dosage: None,
            kind,
            site_flags: site_flags::FIRST | flag,
        }
    }

    #[inline]
    pub fn is_snv(&self) -> bool {
        self.kind == VariantKind::Snv
    }

    #[inline]
    pub fn is_indel(&self) -> bool {
        self.kind == VariantKind::Indel
    }

    /// SNVs, MNPs and indels — the variant classes counted towards TMB.
    #[inline]
    pub fn is_small_variant(&self) -> bool {
        matches!(
            self.kind,
            VariantKind::Snv | VariantKind::Mnp | VariantKind::Indel
        )
    }
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
    /// Records carrying at least one ALT allele (SNV, MNP, indel, symbolic or
    /// breakend). Records with no ALT (`.`, `<NON_REF>`, `<*>`) or only a
    /// spanning-deletion `*` allele are excluded.
    pub total_variants: u64,
    /// Records with at least one SNV allele (`bcftools stats` "number of SNPs").
    pub snp_count: u64,
    /// Records with at least one indel allele (`bcftools stats` "number of indels").
    pub indel_count: u64,
    /// Records with at least one MNP allele (`bcftools stats` "number of MNPs").
    #[serde(default)]
    pub mnp_count: u64,
    /// Records with a symbolic or complex allele (`bcftools stats` "number of others").
    #[serde(default)]
    pub other_count: u64,
    /// Records with more than one ALT allele (`bcftools stats` "number of multiallelic sites").
    #[serde(default)]
    pub multiallelic_sites: u64,
    /// SNV + MNP + indel ALT alleles; the numerator for TMB.
    #[serde(default)]
    pub small_variants: u64,
    /// Transition SNV alleles, counted per ALT allele as `bcftools stats` does.
    #[serde(default)]
    pub transitions: u64,
    /// Transversion SNV alleles, counted per ALT allele.
    #[serde(default)]
    pub transversions: u64,
    /// Transitions / transversions over all SNV ALT alleles (expected ~2.0–2.1 for WGS).
    pub titv_ratio: f64,
    pub per_chrom: HashMap<String, ChromDensity>,
    /// Quality-filtered variants (QUAL > 30), used for VAF-based purity and
    /// per-gene variant evidence.
    pub high_impact: Vec<VariantRecord>,
    /// Number of ALT alleles carrying a SnpEff/VEP impact annotation.
    #[serde(default)]
    pub annotated_variants: u64,
    /// Number of ALT alleles annotated IMPACT=HIGH.
    #[serde(default)]
    pub impact_high_variants: u64,
    /// How `high_impact_genes` was derived: "IMPACT=HIGH (SnpEff/VEP)" when the
    /// VCF is functionally annotated, otherwise "QUAL>30 (no functional annotation)".
    #[serde(default)]
    pub impact_source: String,
    /// Allele frequency histogram: 20 bins over [0.0, 1.0).
    pub af_histogram: Vec<u64>,
    /// Approximate count of unique (chrom, pos) positions.
    pub unique_positions: u64,
    /// Genes with an IMPACT=HIGH variant when the VCF carries SnpEff/VEP
    /// annotation; otherwise genes of QUAL > 30 variants (see `impact_source`).
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
    /// COSMIC mutational signature analysis from 6-channel SBS spectrum.
    #[serde(default)]
    pub cosmic_signatures: Option<MutationalSignatureResult>,
    /// Polygenic risk scores for cancer-relevant GWAS loci.
    #[serde(default)]
    pub prs_scores: Vec<PrsResult>,
    /// Complete variant list, retained in memory for reference-guided analyses
    /// (e.g. 96-channel COSMIC deconvolution) that run after the parse phase.
    /// Never serialized — it would bloat the JSON report with the full VCF.
    #[serde(skip)]
    pub all_variants: Vec<VariantRecord>,
}
