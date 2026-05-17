use ahash::AHashMap;

use biomics_core::{BatchAccum, HyperLogLog};

use crate::types::{ChromDensity, GenomicsSummary, TiTvClass, VariantRecord};

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
    pub total: u64,
    pub snps: u64,
    pub indels: u64,
    pub transitions: u64,
    pub transversions: u64,
    pub per_chrom: AHashMap<String, ChromDensity>,
    pub high_impact: Vec<VariantRecord>,
    /// 20-bin allele frequency histogram over [0.0, 1.0).
    pub af_histogram: [u64; 20],
    /// HyperLogLog cardinality sketch — 16 KB fixed, ~0.81 % error.
    /// Replaces AHashSet<u64> which grew O(n) and consumed hundreds of MB on WGS.
    pub positions: HyperLogLog,
}

impl Default for GenomicsAccum {
    #[inline]
    fn default() -> Self {
        Self {
            total: 0,
            snps: 0,
            indels: 0,
            transitions: 0,
            transversions: 0,
            per_chrom: AHashMap::new(),
            high_impact: Vec::new(),
            af_histogram: [0u64; 20],
            positions: HyperLogLog::new(),
        }
    }
}

impl BatchAccum for GenomicsAccum {
    type Record = VariantRecord;
    type Summary = GenomicsSummary;

    #[inline(always)]
    fn process(&mut self, r: &VariantRecord) -> anyhow::Result<()> {
        self.total += 1;

        let entry = self.per_chrom.entry(r.chrom.clone()).or_default();
        entry.total += 1;

        match r.titv {
            TiTvClass::Transition => {
                self.snps += 1;
                self.transitions += 1;
                entry.snps += 1;
            }
            TiTvClass::Transversion => {
                self.snps += 1;
                self.transversions += 1;
                entry.snps += 1;
            }
            TiTvClass::Indel => {
                self.indels += 1;
                entry.indels += 1;
            }
            TiTvClass::Other => {}
        }

        if r.qual > 30.0 {
            self.high_impact.push(r.clone());
        }

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
        self.transitions += other.transitions;
        self.transversions += other.transversions;

        for (chrom, density) in other.per_chrom {
            let entry = self.per_chrom.entry(chrom).or_default();
            entry.total += density.total;
            entry.snps += density.snps;
            entry.indels += density.indels;
        }

        self.high_impact.extend(other.high_impact);

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

        let mut high_impact_genes: Vec<String> = self
            .high_impact
            .iter()
            .filter_map(|v| v.gene.clone())
            .collect();
        high_impact_genes.sort_unstable();
        high_impact_genes.dedup();

        Ok(GenomicsSummary {
            total_variants: self.total,
            snp_count: self.snps,
            indel_count: self.indels,
            titv_ratio,
            per_chrom: self.per_chrom.into_iter().collect(),
            high_impact: self.high_impact,
            af_histogram: self.af_histogram.to_vec(),
            unique_positions: self.positions.cardinality(),
            high_impact_genes,
        })
    }
}
