//! ATAC-seq de novo peak calling (MACS2-style) from sorted BAM.
//!
//! Implements Poisson-based peak calling from read pileups following the
//! approach described in Zhang et al. (2008) MACS.
//!
//! # Reference
//! Zhang Y, et al. (2008) Model-based analysis of ChIP-seq (MACS).
//! Genome Biology, 9(9):R137.

use ahash::AHashMap;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use noodles_bam as bam;

// ── Public types ──────────────────────────────────────────────────────────────

/// A called ATAC-seq peak.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalledPeak {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// Position of highest coverage within peak.
    pub summit: u64,
    /// Peak height (read depth at summit).
    pub summit_height: u32,
    /// -log10(Poisson p-value) at summit vs local background.
    pub neg_log10_pvalue: f64,
    /// Fold enrichment vs genome-wide lambda.
    pub fold_enrichment: f64,
}

/// Summary of called peaks.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PeakSummary {
    pub total_peaks: u64,
    pub median_peak_width: f64,
    pub median_fold_enrichment: f64,
    /// Fraction of reads in peaks (FRiP).
    pub fraction_reads_in_peaks: f64,
    /// chrom → number of peaks on that chromosome.
    pub per_chrom: AHashMap<String, u64>,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Call peaks from a sorted BAM file.
///
/// # Arguments
/// * `bam_path`       – coordinate-sorted BAM
/// * `shift_size`     – read shift in bp (default: 100 for ATAC, 0 to disable)
/// * `ext_size`       – read extension in bp (default: 200)
/// * `min_mapq`       – mapping quality filter (default: 20)
/// * `p_value_cutoff` – Poisson p-value threshold (default: 1e-5)
pub fn call_peaks(
    bam_path: &std::path::Path,
    shift_size: i64,
    ext_size: u64,
    min_mapq: u8,
    p_value_cutoff: f64,
) -> Result<Vec<CalledPeak>> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .with_context(|| format!("cannot open BAM file: {}", bam_path.display()))?;

    let header = reader.read_header().context("failed to read BAM header")?;

    // Build chromosome name/length list from header.
    let chroms: Vec<(String, u64)> = header
        .reference_sequences()
        .iter()
        .map(|(name, seq)| (name.to_string(), usize::from(seq.length()) as u64))
        .collect();

    let genome_size: u64 = chroms.iter().map(|(_, l)| l).sum();
    if genome_size == 0 {
        return Ok(Vec::new());
    }

    // Allocate per-chromosome pileup vectors.
    let mut pileups: Vec<Vec<u32>> = chroms
        .iter()
        .map(|(_, len)| vec![0u32; *len as usize])
        .collect();

    let mut total_tags: u64 = 0;

    // Walk all BAM records.
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
            Some(pos) => usize::from(pos).saturating_sub(1),
            None => continue,
        };

        let chrom_len = match chroms.get(ref_id) {
            Some((_, l)) => *l as usize,
            None => continue,
        };

        // Determine strand from flags.
        let is_reverse = flags.is_reverse_complemented();

        // Apply strand-aware shift.
        let shifted_start: i64 = if is_reverse {
            start_0 as i64 - shift_size
        } else {
            start_0 as i64 + shift_size
        };

        // Clamp to chromosome bounds after shift.
        let frag_start = shifted_start.max(0) as usize;
        let frag_end = (shifted_start + ext_size as i64)
            .max(frag_start as i64 + 1)
            .min(chrom_len as i64) as usize;

        if frag_start >= frag_end || frag_end > chrom_len {
            continue;
        }

        if let Some(pileup) = pileups.get_mut(ref_id) {
            for pos in frag_start..frag_end {
                if let Some(cell) = pileup.get_mut(pos) {
                    *cell = cell.saturating_add(1);
                }
            }
        }

        total_tags += 1;
    }

    if total_tags == 0 {
        return Ok(Vec::new());
    }

    // Genome-wide lambda = average tags per bp.
    let global_lambda = total_tags as f64 / genome_size as f64;
    let neg_log10_cutoff = -p_value_cutoff.log10();

    let window: u64 = 200;
    let mut peaks: Vec<CalledPeak> = Vec::new();

    for (chrom_idx, (chrom_name, chrom_len)) in chroms.iter().enumerate() {
        let pileup = match pileups.get(chrom_idx) {
            Some(p) => p,
            None => continue,
        };
        let clen = *chrom_len as usize;
        if clen == 0 {
            continue;
        }

        // Sliding-window scan.
        let mut in_peak = false;
        let mut peak_start = 0usize;
        let mut summit_pos = 0usize;
        let mut summit_depth = 0u32;

        let mut pos = 0usize;
        while pos + (window as usize) <= clen {
            let win_end = pos + window as usize;

            // Count tags in window.
            let window_depth: u32 = pileup[pos..win_end].iter().sum();
            let local_depth = window_depth;

            // Local background: use the window sum as the observed count.
            // lambda_local = global_lambda * window (expected tags in window).
            let lambda_local = global_lambda * window as f64;
            let nl10p = neg_log10_poisson(local_depth, lambda_local);

            if nl10p > neg_log10_cutoff {
                if !in_peak {
                    in_peak = true;
                    peak_start = pos;
                    summit_pos = pos;
                    summit_depth = local_depth;
                } else if local_depth > summit_depth {
                    summit_depth = local_depth;
                    summit_pos = pos;
                }
            } else if in_peak {
                // Close peak.
                let peak_end = pos + window as usize;
                let fold = if global_lambda > 0.0 {
                    summit_depth as f64 / lambda_local
                } else {
                    0.0
                };
                peaks.push(CalledPeak {
                    chrom: chrom_name.clone(),
                    start: peak_start as u64,
                    end: peak_end as u64,
                    summit: summit_pos as u64,
                    summit_height: summit_depth,
                    neg_log10_pvalue: neg_log10_poisson(summit_depth, lambda_local),
                    fold_enrichment: fold,
                });
                in_peak = false;
                summit_depth = 0;
            }

            pos += window as usize / 2; // 50 % overlap
        }

        // Close any open peak at chromosome end.
        if in_peak {
            let lambda_local = global_lambda * window as f64;
            let fold = if global_lambda > 0.0 {
                summit_depth as f64 / lambda_local
            } else {
                0.0
            };
            peaks.push(CalledPeak {
                chrom: chrom_name.clone(),
                start: peak_start as u64,
                end: clen as u64,
                summit: summit_pos as u64,
                summit_height: summit_depth,
                neg_log10_pvalue: neg_log10_poisson(summit_depth, lambda_local),
                fold_enrichment: fold,
            });
        }
    }

    // Merge overlapping peaks on the same chromosome.
    peaks = merge_peaks(peaks);

    log::info!(
        "call_peaks: {} peaks called, total_tags={}",
        peaks.len(),
        total_tags,
    );

    Ok(peaks)
}

/// Summarize peak calls.
pub fn summarize_peaks(peaks: &[CalledPeak], total_reads: u64) -> PeakSummary {
    if peaks.is_empty() {
        return PeakSummary::default();
    }

    let mut per_chrom: AHashMap<String, u64> = AHashMap::new();
    let mut widths: Vec<f64> = Vec::with_capacity(peaks.len());
    let mut fes: Vec<f64> = Vec::with_capacity(peaks.len());
    let mut reads_in_peaks: u64 = 0;

    for peak in peaks {
        *per_chrom.entry(peak.chrom.clone()).or_default() += 1;
        widths.push((peak.end.saturating_sub(peak.start)) as f64);
        fes.push(peak.fold_enrichment);
        reads_in_peaks += peak.summit_height as u64;
    }

    widths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    fes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let median_peak_width = median_f64(&widths);
    let median_fold_enrichment = median_f64(&fes);

    let fraction_reads_in_peaks = if total_reads > 0 {
        reads_in_peaks as f64 / total_reads as f64
    } else {
        0.0
    };

    PeakSummary {
        total_peaks: peaks.len() as u64,
        median_peak_width,
        median_fold_enrichment,
        fraction_reads_in_peaks,
        per_chrom,
    }
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Compute -log10(P(X >= k | lambda)) using log-space arithmetic.
///
/// For k == 0: p = 1.0, returns 0.0.
/// For small k (< 50): exact computation in log-space.
/// For large k: normal approximation z = (k - lambda) / sqrt(lambda).
fn neg_log10_poisson(k: u32, lambda: f64) -> f64 {
    if k == 0 {
        return 0.0;
    }
    if lambda <= 0.0 {
        return if k > 0 { f64::INFINITY } else { 0.0 };
    }

    if k < 50 {
        // Compute log P(X = k) = k*ln(lambda) - lambda - lgamma(k+1)
        // Then upper tail: P(X >= k) ≈ P(X = k) for large signal (conservative)
        // For better accuracy, use: P(X >= k) = 1 - P(X <= k-1)
        // Compute sum of P(X = i) for i = 0..k-1 in log-space.
        let log_lambda = lambda.ln();
        let mut log_cdf = f64::NEG_INFINITY; // log(0)

        let mut log_term = -lambda; // log P(X=0) = -lambda
        log_cdf = log_sum_exp(log_cdf, log_term);

        for i in 1..k {
            log_term += log_lambda - (i as f64).ln();
            log_cdf = log_sum_exp(log_cdf, log_term);
        }

        // P(X >= k) = 1 - CDF(k-1)
        let log_survival = log1m_exp(log_cdf);
        let neg_log10_p = -log_survival / std::f64::consts::LN_10;
        neg_log10_p.max(0.0)
    } else {
        // Normal approximation.
        let z = (k as f64 - lambda) / lambda.sqrt();
        // P(Z >= z) ≈ 0.5 * erfc(z / sqrt(2))
        let p = 0.5 * erfc(z / std::f64::consts::SQRT_2);
        if p <= 0.0 {
            return 300.0; // cap at 300
        }
        (-p.log10()).max(0.0)
    }
}

/// log(exp(a) + exp(b)) — numerically stable.
#[inline]
fn log_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (lo - hi).exp().ln_1p()
}

/// log(1 - exp(x)) for x < 0 — numerically stable.
#[inline]
fn log1m_exp(x: f64) -> f64 {
    if x >= 0.0 {
        return f64::NEG_INFINITY;
    }
    if x < -std::f64::consts::LN_2 {
        // |x| is large: use log(1 - exp(x)) directly.
        (-x.exp()).ln_1p()
    } else {
        // |x| is small: use log(-expm1(x)).
        (-x.exp_m1()).ln()
    }
}

/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    poly * (-x * x).exp()
}

/// Merge overlapping/adjacent peaks (sorted by chrom then start).
fn merge_peaks(mut peaks: Vec<CalledPeak>) -> Vec<CalledPeak> {
    if peaks.is_empty() {
        return peaks;
    }

    peaks.sort_unstable_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

    let mut merged: Vec<CalledPeak> = Vec::with_capacity(peaks.len());

    for peak in peaks {
        if let Some(last) = merged.last_mut() {
            if last.chrom == peak.chrom && peak.start <= last.end {
                // Merge: extend end, keep the higher summit.
                last.end = last.end.max(peak.end);
                if peak.summit_height > last.summit_height {
                    last.summit = peak.summit;
                    last.summit_height = peak.summit_height;
                    last.neg_log10_pvalue = peak.neg_log10_pvalue;
                }
                last.fold_enrichment = last.fold_enrichment.max(peak.fold_enrichment);
                continue;
            }
        }
        merged.push(peak);
    }

    merged
}

/// Compute median of a pre-sorted slice.
fn median_f64(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n.is_multiple_of(2) {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// High signal (k=100, lambda=5) should produce a very significant p-value.
    #[test]
    fn poisson_pvalue_high_signal() {
        let nl10p = neg_log10_poisson(100, 5.0);
        assert!(
            nl10p > 10.0,
            "expected neg_log10_pvalue > 10 for k=100, lambda=5, got {}",
            nl10p
        );
    }

    /// summarize_peaks on an empty slice should return all-zeros without panicking.
    #[test]
    fn peaks_summary_empty() {
        let summary = summarize_peaks(&[], 0);
        assert_eq!(summary.total_peaks, 0);
        assert_eq!(summary.median_peak_width, 0.0);
        assert_eq!(summary.median_fold_enrichment, 0.0);
        assert_eq!(summary.fraction_reads_in_peaks, 0.0);
        assert!(summary.per_chrom.is_empty());
    }

    /// CnvClass boundaries tested via the helper imported from genomics_core
    /// would be an integration test; here we just verify the Poisson helper
    /// at boundary values.
    #[test]
    fn segment_class_boundaries() {
        // k == 0 → p-value = 1.0 → -log10(1) = 0
        assert_eq!(neg_log10_poisson(0, 5.0), 0.0);
        // k == lambda (no enrichment) should be small
        let nl10p_baseline = neg_log10_poisson(5, 5.0);
        assert!(
            nl10p_baseline < 2.0,
            "expected near-baseline p-value for k≈lambda, got {}",
            nl10p_baseline
        );
        // Large k well above lambda → highly significant
        let nl10p_sig = neg_log10_poisson(50, 5.0);
        assert!(
            nl10p_sig > 10.0,
            "expected very significant p-value for k=50, lambda=5, got {}",
            nl10p_sig
        );
    }
}
