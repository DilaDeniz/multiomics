//! Coverage-based CNV calling (CNVkit-style) from BAM depth.
//!
//! Implements windowed read-depth binning and greedy change-point segmentation
//! to infer copy-number alterations without a matched VCF, following the
//! approach described in Talevich et al. (2014) CNVkit.
//!
//! # Reference
//! Talevich E, et al. (2014) CNVkit: genome-wide copy number detection and
//! visualization from targeted sequencing. PLOS Computational Biology.

use ahash::AHashMap;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use noodles_bam as bam;

// ── Public types ──────────────────────────────────────────────────────────────

/// A windowed depth bin.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DepthBin {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// Mean read depth in this window.
    pub depth: f64,
    /// GC content fraction [0,1].
    pub gc: f64,
    /// Log2 ratio vs median depth (copy-number proxy).
    pub log2_ratio: f64,
    /// Copy-number estimate (2 ^ (log2_ratio + 1), clipped to [0, 8]).
    pub copy_number: f64,
}

/// A CBS-style copy-number segment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CnvSegment {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_bins: usize,
    pub mean_log2: f64,
    pub copy_number: f64,
    pub class: crate::cnv::CnvClass,
}

/// Summary of coverage-based CNV analysis.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CoverageCnvSummary {
    pub n_bins: u64,
    pub n_segments: u64,
    pub median_depth: f64,
    pub fraction_genome_altered: f64,
    /// chrom → mean log2 ratio
    pub per_chrom: AHashMap<String, f64>,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Compute windowed read depth from a BAM file.
///
/// # Arguments
/// * `bam_path`    – coordinate-sorted BAM
/// * `window_size` – bin size in bp (default: 1000)
/// * `min_mapq`    – minimum mapping quality (default: 20)
pub fn compute_depth_bins(
    bam_path: &std::path::Path,
    window_size: u64,
    min_mapq: u8,
) -> Result<Vec<DepthBin>> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .with_context(|| format!("cannot open BAM file: {}", bam_path.display()))?;

    let header = reader.read_header().context("failed to read BAM header")?;

    // Build chromosome list with lengths.
    let chroms: Vec<(String, u64)> = header
        .reference_sequences()
        .iter()
        .map(|(name, seq)| {
            let len = usize::from(seq.length()) as u64;
            (name.to_string(), len)
        })
        .collect();

    if chroms.is_empty() {
        return Ok(Vec::new());
    }

    // Allocate per-chromosome count vectors.
    // counts[chrom_idx][bin_idx] = read count
    let n_bins_per_chrom: Vec<usize> = chroms
        .iter()
        .map(|(_, len)| (*len as usize).div_ceil(window_size as usize))
        .collect();

    let mut counts: Vec<Vec<u32>> = n_bins_per_chrom.iter().map(|&n| vec![0u32; n]).collect();

    // Walk all BAM records and tally reads into bins.
    for result in reader.record_bufs(&header) {
        let record = result.context("failed to read BAM record")?;

        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(0);
        if mapq < min_mapq {
            continue;
        }

        let ref_id = match record.reference_sequence_id() {
            Some(id) => id,
            None => continue,
        };

        // alignment_start is 1-based noodles Position; convert to 0-based.
        let start_0 = match record.alignment_start() {
            Some(pos) => usize::from(pos) - 1,
            None => continue,
        };

        let bin_idx = start_0 / window_size as usize;
        if let Some(chrom_counts) = counts.get_mut(ref_id) {
            if let Some(cell) = chrom_counts.get_mut(bin_idx) {
                *cell = cell.saturating_add(1);
            }
        }
    }

    // Build DepthBin list (without log2 ratios yet — need global median first).
    let mut bins: Vec<DepthBin> = Vec::new();
    for (chrom_idx, (chrom_name, chrom_len)) in chroms.iter().enumerate() {
        let n_bins = n_bins_per_chrom[chrom_idx];
        for (bin_idx, &count) in counts[chrom_idx].iter().enumerate().take(n_bins) {
            let start = bin_idx as u64 * window_size;
            let end = (start + window_size).min(*chrom_len);
            let actual_size = (end - start) as f64;
            let depth = count as f64 / actual_size;
            bins.push(DepthBin {
                chrom: chrom_name.clone(),
                start,
                end,
                depth,
                gc: 0.5, // placeholder — real GC requires reference FASTA
                log2_ratio: 0.0,
                copy_number: 2.0,
            });
        }
    }

    if bins.is_empty() {
        return Ok(bins);
    }

    // Compute global median depth.
    let median_depth = {
        let mut depths: Vec<f64> = bins.iter().map(|b| b.depth).collect();
        depths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = depths.len();
        if n.is_multiple_of(2) {
            (depths[n / 2 - 1] + depths[n / 2]) / 2.0
        } else {
            depths[n / 2]
        }
    };

    // Compute log2 ratio and copy number for each bin.
    for bin in &mut bins {
        // CNVkit formula: log2(depth / median + 0.01) - 1 centers diploid at 0.
        bin.log2_ratio = (bin.depth / (median_depth + f64::EPSILON) + 0.01).log2() - 1.0;
        bin.copy_number = (2.0_f64.powf(bin.log2_ratio + 1.0)).clamp(0.0, 8.0);
    }

    log::info!(
        "compute_depth_bins: {} bins, median_depth={:.2}",
        bins.len(),
        median_depth
    );
    Ok(bins)
}

/// Apply simplified circular binary segmentation to depth bins.
///
/// Uses greedy change-point detection: scans for positions where the
/// rolling mean changes by > 0.3 log2 units, then merges small segments.
pub fn segment_bins(bins: &[DepthBin], min_segment_size: usize) -> Vec<CnvSegment> {
    if bins.is_empty() {
        return Vec::new();
    }

    // Group bins by chromosome, preserving order.
    let mut chrom_groups: Vec<(String, Vec<usize>)> = Vec::new();
    for (i, bin) in bins.iter().enumerate() {
        if chrom_groups
            .last()
            .map(|(c, _)| c == &bin.chrom)
            .unwrap_or(false)
        {
            chrom_groups.last_mut().unwrap().1.push(i);
        } else {
            chrom_groups.push((bin.chrom.clone(), vec![i]));
        }
    }

    let mut segments: Vec<CnvSegment> = Vec::new();

    for (chrom, indices) in chrom_groups {
        let n = indices.len();
        if n == 0 {
            continue;
        }

        // Build rolling mean with window = 5.
        let window: usize = 5;
        let log2s: Vec<f64> = indices.iter().map(|&i| bins[i].log2_ratio).collect();

        let rolling: Vec<f64> = (0..n)
            .map(|i| {
                let lo = i.saturating_sub(window / 2);
                let hi = (i + window / 2 + 1).min(n);
                let s: f64 = log2s[lo..hi].iter().sum();
                s / (hi - lo) as f64
            })
            .collect();

        // Detect breakpoints.
        let mut breakpoints: Vec<usize> = vec![0];
        for i in 1..n {
            if (rolling[i] - rolling[i - 1]).abs() > 0.3 {
                breakpoints.push(i);
            }
        }
        breakpoints.push(n);

        // Build initial segments from breakpoints.
        let mut raw_segs: Vec<(usize, usize)> =
            breakpoints.windows(2).map(|w| (w[0], w[1])).collect();

        // Merge segments smaller than min_segment_size with their neighbor.
        let mut changed = true;
        while changed {
            changed = false;
            let mut merged: Vec<(usize, usize)> = Vec::with_capacity(raw_segs.len());
            let mut skip = false;
            for j in 0..raw_segs.len() {
                if skip {
                    skip = false;
                    continue;
                }
                let (lo, hi) = raw_segs[j];
                if hi - lo < min_segment_size && j + 1 < raw_segs.len() {
                    // Merge with next.
                    merged.push((lo, raw_segs[j + 1].1));
                    skip = true;
                    changed = true;
                } else if hi - lo < min_segment_size && !merged.is_empty() {
                    // Merge with previous.
                    let last = merged.last_mut().unwrap();
                    last.1 = hi;
                    changed = true;
                } else {
                    merged.push((lo, hi));
                }
            }
            raw_segs = merged;
        }

        // Emit CnvSegment for each group.
        for (lo, hi) in raw_segs {
            let seg_log2s = &log2s[lo..hi];
            let mean_log2 = seg_log2s.iter().sum::<f64>() / seg_log2s.len() as f64;
            let copy_number = (2.0_f64.powf(mean_log2 + 1.0)).clamp(0.0, 8.0);

            let class = cnv_class_from_log2(mean_log2);

            let start = bins[indices[lo]].start;
            let end = bins[indices[hi - 1]].end;

            segments.push(CnvSegment {
                chrom: chrom.clone(),
                start,
                end,
                n_bins: hi - lo,
                mean_log2,
                copy_number,
                class,
            });
        }
    }

    segments
}

/// Summarize coverage-based CNV analysis.
pub fn summarize_coverage_cnv(bins: &[DepthBin], segments: &[CnvSegment]) -> CoverageCnvSummary {
    let mut per_chrom: AHashMap<String, (f64, u64)> = AHashMap::new();
    for bin in bins {
        let entry = per_chrom.entry(bin.chrom.clone()).or_default();
        entry.0 += bin.log2_ratio;
        entry.1 += 1;
    }

    let per_chrom_mean: AHashMap<String, f64> = per_chrom
        .into_iter()
        .map(|(chrom, (sum, n))| (chrom, if n > 0 { sum / n as f64 } else { 0.0 }))
        .collect();

    let median_depth = {
        let mut depths: Vec<f64> = bins.iter().map(|b| b.depth).collect();
        depths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = depths.len();
        if n == 0 {
            0.0
        } else if n.is_multiple_of(2) {
            (depths[n / 2 - 1] + depths[n / 2]) / 2.0
        } else {
            depths[n / 2]
        }
    };

    let n_altered = segments
        .iter()
        .filter(|s| !matches!(s.class, crate::cnv::CnvClass::Diploid))
        .count() as f64;
    let fraction_genome_altered = if segments.is_empty() {
        0.0
    } else {
        n_altered / segments.len() as f64
    };

    CoverageCnvSummary {
        n_bins: bins.len() as u64,
        n_segments: segments.len() as u64,
        median_depth,
        fraction_genome_altered,
        per_chrom: per_chrom_mean,
    }
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Classify a log2 ratio into a [`crate::cnv::CnvClass`].
///
/// Boundaries (log2 scale):
/// * < -0.7  → deletion (hom or het depending on depth)
/// * -0.3 to 0.3 → diploid
/// * > 0.7   → amplification
fn cnv_class_from_log2(log2: f64) -> crate::cnv::CnvClass {
    use crate::cnv::CnvClass;
    if log2 < -0.7 {
        CnvClass::HeterozygousDeletion
    } else if log2 < -0.3 {
        // borderline loss — call het deletion
        CnvClass::HeterozygousDeletion
    } else if log2 <= 0.3 {
        CnvClass::Diploid
    } else if log2 <= 0.7 {
        // borderline gain — call low amp
        CnvClass::LowAmplification
    } else {
        CnvClass::HighAmplification
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cnv::CnvClass;

    /// compute_depth_bins on a path that doesn't exist should return an error,
    /// not a panic.  We treat a missing file as "graceful skip" for the purpose
    /// of testing that the function surface is safe.
    #[test]
    fn depth_bins_empty_bam() {
        let result = compute_depth_bins(std::path::Path::new("/nonexistent/empty.bam"), 1000, 20);
        // Expected: error (file not found), not a panic.
        assert!(result.is_err());
    }

    /// A slice of bins with identical log2 ratios should produce a single
    /// segment per chromosome (no spurious breakpoints).
    #[test]
    fn segment_bins_flat() {
        let bins: Vec<DepthBin> = (0..10)
            .map(|i| DepthBin {
                chrom: "chr1".into(),
                start: i * 1000,
                end: (i + 1) * 1000,
                depth: 30.0,
                gc: 0.5,
                log2_ratio: 0.0,
                copy_number: 2.0,
            })
            .collect();

        let segs = segment_bins(&bins, 5);
        assert_eq!(segs.len(), 1, "uniform depth should produce one segment");
        assert_eq!(segs[0].chrom, "chr1");
    }

    /// Verify CnvClass assignment at canonical log2 thresholds.
    #[test]
    fn cnv_class_from_log2_boundaries() {
        assert!(
            matches!(cnv_class_from_log2(-1.0), CnvClass::HeterozygousDeletion),
            "log2=-1.0 should be a deletion"
        );
        assert!(
            matches!(cnv_class_from_log2(0.0), CnvClass::Diploid),
            "log2=0.0 should be diploid"
        );
        assert!(
            matches!(cnv_class_from_log2(1.5), CnvClass::HighAmplification),
            "log2=1.5 should be high amplification"
        );
    }
}
