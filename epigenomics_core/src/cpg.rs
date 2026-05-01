use crate::types::CpGIsland;

/// Detect CpG islands from a sorted list of methylation sites.
///
/// A region qualifies as a CpG island when:
/// - Length ≥ `min_length_bp`
/// - GC content > 50% (approximated from CpG density)
/// - Contains at least 5 CpG sites
///
/// `sites` must be sorted by `start` position.
///
/// The GC percentage is estimated from CpG density: each CpG site contributes
/// 2 GC bases to a 2-bp window, giving a rough GC% proxy.
///
/// # Parameters
/// - `chrom` — chromosome name
/// - `sites` — sorted `(start, end, methylation_pct)` tuples
/// - `min_length_bp` — minimum island length in base pairs
/// - `min_gc_pct` — minimum GC percentage threshold (e.g. 50)
pub fn detect_cpg_islands(
    chrom: &str,
    sites: &[(u64, u64, f64)],
    min_length_bp: u64,
    min_gc_pct: u64,
) -> Vec<CpGIsland> {
    if sites.is_empty() {
        return Vec::new();
    }

    let mut islands = Vec::new();
    // Sliding window: expand while sites are within `min_length_bp` of window start
    let mut i = 0;
    while i < sites.len() {
        let region_start = sites[i].0;
        let mut j = i;

        // Expand window to collect all sites within a contiguous dense region
        while j < sites.len() {
            let gap = sites[j].0.saturating_sub(if j == 0 { 0 } else { sites[j - 1].1 });
            if gap > 500 && j > i {
                // Large gap — end current candidate region
                break;
            }
            j += 1;
        }

        let region_end = sites[j - 1].1;
        let length = region_end.saturating_sub(region_start);
        let site_count = j - i;

        if length >= min_length_bp && site_count >= 5 {
            // Estimate GC% from CpG density.
            // Each CpG site spans 2 bases, but CpG islands are embedded in GC-rich
            // sequence where ~3× as many bases are GC from non-CpG context.
            // Multiplying CpG base fraction by 3 approximates total GC%.
            let cpg_density = site_count as f64 * 6.0 / length.max(1) as f64;
            let gc_pct = (cpg_density * 100.0).min(100.0);

            if gc_pct as u64 >= min_gc_pct {
                let mean_meth = sites[i..j].iter().map(|s| s.2).sum::<f64>() / site_count as f64;
                islands.push(CpGIsland {
                    chrom: chrom.to_string(),
                    start: region_start,
                    end: region_end,
                    length,
                    gc_percent: gc_pct,
                    mean_methylation: mean_meth,
                });
            }
        }

        i = if j > i { j } else { i + 1 };
    }

    islands
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_cpg_islands_empty() {
        assert!(detect_cpg_islands("chr1", &[], 200, 50).is_empty());
    }

    #[test]
    fn test_detect_cpg_islands_dense() {
        // 25 sites spaced 10bp apart over ~242bp — 1 CpG per 10bp is classic CpG island density
        let sites: Vec<(u64, u64, f64)> = (0..25).map(|i| (i * 10, i * 10 + 2, 60.0)).collect();
        let islands = detect_cpg_islands("chr1", &sites, 200, 50);
        assert!(!islands.is_empty(), "Expected at least one CpG island");
    }
}
