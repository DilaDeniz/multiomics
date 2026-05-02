use ahash::AHashMap;

use biomics_core::BatchAccum;

use crate::cpg::detect_cpg_islands;
use crate::types::{ChromMethylation, EpigenomicsSummary, MethylationRecord, MethylationRegion, RegionKind};

/// Lock-free accumulator for methylation statistics.
///
/// `AHashMap` replaces `std::HashMap` for 3-5× faster hash table operations.
pub struct EpigenomicsAccum {
    pub total_sites: u64,
    pub sum_methylation: f64,
    pub per_chrom: AHashMap<String, ChromMethylation>,
    /// Per-chromosome site list for CpG island detection: (start, end, methylation_pct).
    pub chrom_sites: AHashMap<String, Vec<(u64, u64, f64)>>,
}

impl Default for EpigenomicsAccum {
    #[inline]
    fn default() -> Self {
        Self {
            total_sites: 0,
            sum_methylation: 0.0,
            per_chrom: AHashMap::new(),
            chrom_sites: AHashMap::new(),
        }
    }
}

impl BatchAccum for EpigenomicsAccum {
    type Record = MethylationRecord;
    type Summary = EpigenomicsSummary;

    #[inline(always)]
    fn process(&mut self, r: &MethylationRecord) -> anyhow::Result<()> {
        self.total_sites += 1;
        self.sum_methylation += r.methylation;

        let entry = self.per_chrom.entry(r.chrom.clone()).or_default();
        entry.total_sites += 1;
        entry.sum_methylation += r.methylation;

        self.chrom_sites
            .entry(r.chrom.clone())
            .or_default()
            .push((r.start, r.end, r.methylation));

        Ok(())
    }

    #[inline(always)]
    fn merge(&mut self, other: Self) {
        self.total_sites += other.total_sites;
        self.sum_methylation += other.sum_methylation;

        for (chrom, cm) in other.per_chrom {
            let entry = self.per_chrom.entry(chrom).or_default();
            entry.total_sites += cm.total_sites;
            entry.sum_methylation += cm.sum_methylation;
        }

        for (chrom, sites) in other.chrom_sites {
            self.chrom_sites.entry(chrom).or_default().extend(sites);
        }
    }

    fn finalize(mut self) -> anyhow::Result<EpigenomicsSummary> {
        let global_methylation_pct = if self.total_sites == 0 {
            0.0
        } else {
            self.sum_methylation / self.total_sites as f64
        };

        for cm in self.per_chrom.values_mut() {
            cm.mean_methylation = if cm.total_sites == 0 {
                0.0
            } else {
                cm.sum_methylation / cm.total_sites as f64
            };
        }

        // Sort each chromosome's sites by start position for CpG island detection
        for sites in self.chrom_sites.values_mut() {
            sites.sort_unstable_by_key(|s| s.0);
        }

        let mut cpg_islands = Vec::new();
        for (chrom, sites) in &self.chrom_sites {
            cpg_islands.extend(detect_cpg_islands(chrom, sites, 200, 50));
        }

        let mut hypermethylated = Vec::new();
        let mut hypomethylated = Vec::new();

        for (chrom, sites) in &self.chrom_sites {
            for region in sliding_window_regions(chrom, sites, 5) {
                if region.mean_methylation > 80.0 {
                    hypermethylated.push(region);
                } else if region.mean_methylation < 20.0 {
                    hypomethylated.push(region);
                }
            }
        }

        Ok(EpigenomicsSummary {
            total_sites: self.total_sites,
            global_methylation_pct,
            per_chrom: self.per_chrom.into_iter().collect(),
            cpg_islands,
            hypermethylated,
            hypomethylated,
        })
    }
}

/// Group consecutive sites into regions using a sliding window of `window_size` sites.
#[inline]
fn sliding_window_regions(
    chrom: &str,
    sites: &[(u64, u64, f64)],
    window_size: usize,
) -> Vec<MethylationRegion> {
    if sites.len() < window_size {
        return Vec::new();
    }
    let mut regions = Vec::with_capacity(sites.len() - window_size + 1);
    for window in sites.windows(window_size) {
        let start = window[0].0;
        let end = window[window.len() - 1].1;
        let mean = window.iter().map(|s| s.2).sum::<f64>() / window.len() as f64;
        let kind = if mean > 80.0 { RegionKind::Hypermethylated } else { RegionKind::Hypomethylated };
        regions.push(MethylationRegion { chrom: chrom.to_string(), start, end, mean_methylation: mean, kind });
    }
    regions
}
