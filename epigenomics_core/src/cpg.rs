use crate::types::CpGIsland;

/// Detect CpG islands from a sorted list of methylation sites using an
/// overlapping sliding-window algorithm based on the Gardiner-Garden &
/// Frommer (1987) criteria adapted for WGBS site data.
///
/// ## Algorithm
/// 1. Slide a 200 bp window anchored at each site (two-pointer, O(n)).
/// 2. Compute CpO/E = `observed / (WINDOW × 0.0625)` where 0.0625 = p(C)×p(G)
///    assuming 50% GC content.
/// 3. A window qualifies when CpO/E ≥ 0.6 AND ≥ 5 sites are present.
/// 4. Adjacent / overlapping qualifying windows are merged.
/// 5. Merged regions shorter than `min_length_bp` are discarded.
///
/// ## Parameters
/// - `chrom`         — chromosome name (copied into each output island)
/// - `sites`         — sorted `(start, end, methylation_pct)` tuples
/// - `min_length_bp` — minimum island length in base pairs (typically 200)
/// - `_min_gc_pct`   — kept for API compatibility; threshold encoded in CpO/E
pub fn detect_cpg_islands(
    chrom: &str,
    sites: &[(u64, u64, f64)],
    min_length_bp: u64,
    _min_gc_pct: u64,
) -> Vec<CpGIsland> {
    const WINDOW: u64 = 200;
    // Expected CpG dinucleotide rate at 50% GC: p(C)×p(G) = 0.25×0.25 = 0.0625
    const EXPECTED_CPG_RATE: f64 = 0.0625;
    const CPOE_MIN: f64 = 0.6;
    const MIN_SITES: usize = 5;

    let n = sites.len();
    if n < MIN_SITES {
        return Vec::new();
    }

    // Two-pointer: both i and j only advance → O(n) total
    let mut qualifying: Vec<(u64, u64)> = Vec::new();
    let mut j = 0usize;

    for i in 0..n {
        let win_start = sites[i].0;
        let win_end = win_start + WINDOW;

        while j < n && sites[j].0 < win_end {
            j += 1;
        }

        let count = j - i;
        if count < MIN_SITES {
            continue;
        }

        let cpoe = count as f64 / (WINDOW as f64 * EXPECTED_CPG_RATE);
        if cpoe < CPOE_MIN {
            continue;
        }

        // Merge into the previous qualifying window if they overlap
        if let Some(last) = qualifying.last_mut() {
            if win_start <= last.1 {
                last.1 = last.1.max(win_end);
                continue;
            }
        }
        qualifying.push((win_start, win_end));
    }

    if qualifying.is_empty() {
        return Vec::new();
    }

    // Convert merged windows into CpGIsland structs
    let mut islands = Vec::with_capacity(qualifying.len());
    for (start, end) in qualifying {
        let length = end - start;
        if length < min_length_bp {
            continue;
        }

        // Binary search for the sites contained in this island
        let lo = sites.partition_point(|s| s.0 < start);
        let hi = sites.partition_point(|s| s.0 < end);
        let island_sites = &sites[lo..hi];

        if island_sites.len() < MIN_SITES {
            continue;
        }

        let cpoe = island_sites.len() as f64 / (length as f64 * EXPECTED_CPG_RATE);
        // GC% approximation: each CpG = 2 GC bases; multiply by 3 for flanking GC context
        let gc_pct = (island_sites.len() as f64 * 6.0 / length as f64 * 100.0).min(100.0);
        let mean_meth =
            island_sites.iter().map(|s| s.2).sum::<f64>() / island_sites.len() as f64;

        islands.push(CpGIsland {
            chrom: chrom.to_string(),
            start,
            end,
            length,
            gc_percent: gc_pct,
            mean_methylation: mean_meth,
            cpoe,
        });
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
        // 25 sites 10 bp apart; 200-bp window contains 20 sites
        // CpO/E = 20 / (200 × 0.0625) = 1.6 — well above threshold
        let sites: Vec<(u64, u64, f64)> = (0..25).map(|i| (i * 10, i * 10 + 2, 60.0)).collect();
        let islands = detect_cpg_islands("chr1", &sites, 200, 50);
        assert!(!islands.is_empty(), "Expected at least one CpG island");
        assert!(islands[0].cpoe >= 0.6, "CpO/E should meet the 0.6 threshold");
    }

    #[test]
    fn test_detect_cpg_islands_sparse() {
        // 5 sites 400 bp apart: CpO/E = 1/(200×0.0625) ≈ 0.08 → below threshold
        let sites: Vec<(u64, u64, f64)> = (0..5).map(|i| (i * 400, i * 400 + 2, 50.0)).collect();
        assert!(detect_cpg_islands("chr1", &sites, 200, 50).is_empty());
    }
}
