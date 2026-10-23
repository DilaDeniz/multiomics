use ahash::AHashMap;

use biomics_core::{BatchAccum, HyperLogLog};

use crate::cancer::{compute_hrd_score, compute_msi, detect_kataegis, detect_loh};
use crate::types::{site_flags, ChromDensity, GenomicsSummary, TiTvClass, VariantRecord};

/// Compact position key: FNV-1a hash of (chrom, pos).
///
/// Using `AHashSet<u64>` instead of `AHashSet<(String, u64)>` eliminates one
/// String heap-allocation per variant during the parallel fold.
#[inline(always)]
fn position_key(chrom: &str, pos: u64) -> u64 {
    // FNV-1a over chrom bytes, then mix pos
    let mut h: u64 = 14_695_981_039_346_656_037;
    for b in chrom.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(1_099_511_628_211);
    }
    h ^= pos;
    h = h.wrapping_mul(1_099_511_628_211);
    h ^ (h >> 17)
}

/// Lock-free accumulator for genomic variant statistics.
///
/// Each rayon worker holds an independent instance; `merge` is called once
/// per worker at finalization. No shared mutable state during `process`.
pub struct GenomicsAccum {
    /// VCF records carrying at least one variant allele.
    pub total: u64,
    /// Records with an SNV / indel / MNP / other allele (bcftools semantics).
    pub snps: u64,
    pub indels: u64,
    pub mnps: u64,
    pub others: u64,
    pub multiallelic: u64,
    /// SNV + MNP + indel alleles (TMB numerator).
    pub small_variants: u64,
    /// Alleles with a SnpEff/VEP impact annotation, and those with IMPACT=HIGH.
    pub annotated: u64,
    pub impact_high: u64,
    /// Genes carrying an IMPACT=HIGH allele.
    pub impact_high_genes: ahash::AHashSet<String>,
    /// Ti / Tv counted per SNV allele.
    pub transitions: u64,
    pub transversions: u64,
    pub per_chrom: AHashMap<String, ChromDensity>,
    pub high_impact: Vec<VariantRecord>,
    /// 20-bin allele frequency histogram over [0.0, 1.0).
    pub af_histogram: [u64; 20],
    /// HyperLogLog cardinality sketch — 16 KB fixed, ~0.81 % error.
    /// Replaces AHashSet<u64> which grew O(n) and consumed hundreds of MB on WGS.
    pub positions: HyperLogLog,
    /// Full variant list (all variants, not just high-impact) for cancer analyses.
    pub all_variants: Vec<VariantRecord>,
}

impl Default for GenomicsAccum {
    #[inline]
    fn default() -> Self {
        Self {
            total: 0,
            snps: 0,
            indels: 0,
            mnps: 0,
            others: 0,
            multiallelic: 0,
            small_variants: 0,
            annotated: 0,
            impact_high: 0,
            impact_high_genes: ahash::AHashSet::new(),
            transitions: 0,
            transversions: 0,
            per_chrom: AHashMap::new(),
            high_impact: Vec::new(),
            af_histogram: [0u64; 20],
            positions: HyperLogLog::new(),
            all_variants: Vec::new(),
        }
    }
}

impl BatchAccum for GenomicsAccum {
    type Record = VariantRecord;
    type Summary = GenomicsSummary;

    #[inline(always)]
    fn process(&mut self, r: &VariantRecord) -> anyhow::Result<()> {
        // Site-level counts: one per VCF record, carried by its first allele.
        let flags = r.site_flags;
        if flags & site_flags::FIRST != 0 {
            let has = |bit: u8| u64::from(flags & bit != 0);
            self.total += 1;
            self.snps += has(site_flags::SNV);
            self.indels += has(site_flags::INDEL);
            self.mnps += has(site_flags::MNP);
            self.others += has(site_flags::OTHER);
            self.multiallelic += has(site_flags::MULTIALLELIC);

            let entry = self.per_chrom.entry(r.chrom.clone()).or_default();
            entry.total += 1;
            entry.snps += has(site_flags::SNV);
            entry.indels += has(site_flags::INDEL);
        }

        // Allele-level counts.
        match r.titv {
            TiTvClass::Transition => self.transitions += 1,
            TiTvClass::Transversion => self.transversions += 1,
            TiTvClass::Indel | TiTvClass::Other => {}
        }
        if r.is_small_variant() {
            self.small_variants += 1;
        }
        if let Some(impact) = r.impact {
            self.annotated += 1;
            if impact == crate::types::Impact::High {
                self.impact_high += 1;
                if let Some(gene) = &r.gene {
                    self.impact_high_genes.insert(gene.clone());
                }
            }
        }

        if r.qual > 30.0 {
            self.high_impact.push(r.clone());
        }
        self.all_variants.push(r.clone());

        if let Some(af) = r.af {
            let bin = biomics_core::stats::histogram_bin(af as f64, 0.0, 1.0, 20);
            self.af_histogram[bin] += 1;
        }

        self.positions.insert_hashed(position_key(&r.chrom, r.pos));

        Ok(())
    }

    #[inline(always)]
    fn merge(&mut self, other: Self) {
        self.total += other.total;
        self.snps += other.snps;
        self.indels += other.indels;
        self.mnps += other.mnps;
        self.others += other.others;
        self.multiallelic += other.multiallelic;
        self.small_variants += other.small_variants;
        self.annotated += other.annotated;
        self.impact_high += other.impact_high;
        self.impact_high_genes.extend(other.impact_high_genes);
        self.transitions += other.transitions;
        self.transversions += other.transversions;

        for (chrom, density) in other.per_chrom {
            let entry = self.per_chrom.entry(chrom).or_default();
            entry.total += density.total;
            entry.snps += density.snps;
            entry.indels += density.indels;
        }

        self.high_impact.extend(other.high_impact);
        self.all_variants.extend(other.all_variants);

        for (i, count) in other.af_histogram.iter().enumerate() {
            self.af_histogram[i] += count;
        }

        self.positions.merge(&other.positions);
    }

    fn finalize(self) -> anyhow::Result<GenomicsSummary> {
        let titv_ratio = if self.transversions == 0 {
            0.0
        } else {
            self.transitions as f64 / self.transversions as f64
        };

        // With functional annotation, "high impact" means SnpEff/VEP IMPACT=HIGH;
        // without it, fall back to the genes of quality-filtered variants.
        let (mut high_impact_genes, impact_source): (Vec<String>, &str) = if self.annotated > 0 {
            (
                self.impact_high_genes.into_iter().collect(),
                "IMPACT=HIGH (SnpEff/VEP)",
            )
        } else {
            (
                self.high_impact
                    .iter()
                    .filter_map(|v| v.gene.clone())
                    .collect(),
                "QUAL>30 (no functional annotation)",
            )
        };
        high_impact_genes.sort_unstable();
        high_impact_genes.dedup();

        // Cancer analyses run on the full variant set
        let kataegis_loci = detect_kataegis(&self.all_variants);
        let hrd = Some(compute_hrd_score(&self.all_variants));
        let loh_chromosomes = detect_loh(&self.all_variants);
        let msi = Some(compute_msi(&self.all_variants));
        let cosmic_signatures = Some(crate::cosmic::compute_mutational_signatures(
            &self.all_variants,
        ));
        let prs_scores = crate::prs::compute_prs(&self.all_variants);

        Ok(GenomicsSummary {
            total_variants: self.total,
            snp_count: self.snps,
            indel_count: self.indels,
            mnp_count: self.mnps,
            other_count: self.others,
            multiallelic_sites: self.multiallelic,
            small_variants: self.small_variants,
            transitions: self.transitions,
            transversions: self.transversions,
            titv_ratio,
            per_chrom: self.per_chrom.into_iter().collect(),
            high_impact: self.high_impact,
            annotated_variants: self.annotated,
            impact_high_variants: self.impact_high,
            impact_source: impact_source.to_string(),
            af_histogram: self.af_histogram.to_vec(),
            unique_positions: self.positions.cardinality(),
            high_impact_genes,
            tumor_purity: None,
            kataegis_loci,
            hrd,
            loh_chromosomes,
            tmb: None, // filled in by runner.rs after context detection
            msi,
            cosmic_signatures,
            prs_scores,
            all_variants: self.all_variants,
        })
    }
}
