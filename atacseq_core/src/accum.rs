use ahash::AHashMap;

use biomics_core::BatchAccum;

use crate::types::{AtacPeak, AtacSummary, ChromPeakStats};

/// Maximum number of top peaks retained per thread during accumulation.
/// Trimmed to `TOP_PEAKS_FINAL` at finalize time.
const TOP_PEAKS_ACCUM: usize = 200;

/// Number of top peaks reported in the final summary.
const TOP_PEAKS_FINAL: usize = 100;

/// Lock-free accumulator for ATAC-seq narrowPeak statistics.
///
/// Uses `AHashMap` for 3–5× faster per-chromosome bookkeeping than `std::HashMap`.
/// Top-peak tracking uses a bounded vector: we keep up to `TOP_PEAKS_ACCUM` peaks
/// and trim only when the vector exceeds that bound, deferring the sort cost to
/// merge/finalize time.
pub struct AtacAccum {
    /// Per-chrom running totals: (peak_count, open_bp, signal_sum, width_sum).
    per_chrom: AHashMap<String, (u64, u64, f64, u64)>,
    /// Candidate top peaks — sorted and trimmed during merge and finalize.
    top_peaks: Vec<AtacPeak>,
    /// Minimum signal among current `top_peaks` candidates (used as fast guard).
    top_min_signal: f64,
    /// Total peaks seen by this accumulator.
    total_peaks: u64,
    /// Sum of signal values across all peaks.
    signal_sum: f64,
    /// Sum of all peak widths (open-chromatin base pairs).
    open_bp_total: u64,
    /// Number of peaks with width < 500 bp.
    narrow_peaks: u64,
    /// All signal values (for median computation at finalize time).
    signal_values: Vec<f64>,
}

impl Default for AtacAccum {
    #[inline]
    fn default() -> Self {
        Self {
            per_chrom: AHashMap::new(),
            top_peaks: Vec::with_capacity(TOP_PEAKS_ACCUM + 1),
            top_min_signal: f64::NEG_INFINITY,
            total_peaks: 0,
            signal_sum: 0.0,
            open_bp_total: 0,
            narrow_peaks: 0,
            signal_values: Vec::new(),
        }
    }
}

impl AtacAccum {
    /// Trim `top_peaks` to `TOP_PEAKS_ACCUM` by descending signal value.
    /// Only called when `top_peaks.len() > TOP_PEAKS_ACCUM`.
    #[inline]
    fn trim_top_peaks(&mut self) {
        // Partial sort: place the top-N by signal in the first N slots.
        self.top_peaks.sort_unstable_by(|a, b| {
            b.signal_value
                .partial_cmp(&a.signal_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        self.top_peaks.truncate(TOP_PEAKS_ACCUM);
        self.top_min_signal = self
            .top_peaks
            .last()
            .map(|p| p.signal_value)
            .unwrap_or(f64::NEG_INFINITY);
    }
}

impl BatchAccum for AtacAccum {
    type Record = AtacPeak;
    type Summary = AtacSummary;

    #[inline(always)]
    fn process(&mut self, r: &AtacPeak) -> anyhow::Result<()> {
        let width = r.width();

        // Global counters.
        self.total_peaks += 1;
        self.signal_sum += r.signal_value;
        self.open_bp_total += width;
        self.signal_values.push(r.signal_value);
        if r.is_narrow() {
            self.narrow_peaks += 1;
        }

        // Per-chromosome counters.
        let entry = self
            .per_chrom
            .entry(r.chrom.clone())
            .or_insert((0, 0, 0.0, 0));
        entry.0 += 1;
        entry.1 += width;
        entry.2 += r.signal_value;
        entry.3 += width;

        // Top-peak tracking: admit if there is room or if this peak beats the minimum.
        if self.top_peaks.len() < TOP_PEAKS_ACCUM || r.signal_value > self.top_min_signal {
            self.top_peaks.push(r.clone());
            // Trim only when the buffer overflows to avoid O(n) work every record.
            if self.top_peaks.len() > TOP_PEAKS_ACCUM {
                self.trim_top_peaks();
            } else if r.signal_value < self.top_min_signal
                || self.top_min_signal == f64::NEG_INFINITY
            {
                self.top_min_signal = r.signal_value;
            }
        }

        Ok(())
    }

    #[inline(always)]
    fn merge(&mut self, other: Self) {
        self.total_peaks += other.total_peaks;
        self.signal_sum += other.signal_sum;
        self.open_bp_total += other.open_bp_total;
        self.narrow_peaks += other.narrow_peaks;
        self.signal_values.extend(other.signal_values);

        // Merge per-chrom maps by accumulating into self.
        for (chrom, (peaks, open_bp, sig_sum, width_sum)) in other.per_chrom {
            let entry = self.per_chrom.entry(chrom).or_insert((0, 0, 0.0, 0));
            entry.0 += peaks;
            entry.1 += open_bp;
            entry.2 += sig_sum;
            entry.3 += width_sum;
        }

        // Merge top-peaks candidate lists, then re-trim.
        self.top_peaks.extend(other.top_peaks);
        if self.top_peaks.len() > TOP_PEAKS_ACCUM {
            self.trim_top_peaks();
        }
    }

    fn finalize(mut self) -> anyhow::Result<AtacSummary> {
        let total = self.total_peaks;

        let mean_peak_width = if total == 0 {
            0.0
        } else {
            self.open_bp_total as f64 / total as f64
        };

        let mean_signal_value = if total == 0 {
            0.0
        } else {
            self.signal_sum / total as f64
        };

        // Compute exact median by sorting all signal values.
        let median_signal_value = if self.signal_values.is_empty() {
            0.0
        } else {
            self.signal_values
                .sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let n = self.signal_values.len();
            if n % 2 == 1 {
                self.signal_values[n / 2]
            } else {
                (self.signal_values[n / 2 - 1] + self.signal_values[n / 2]) / 2.0
            }
        };

        let fraction_narrow_peaks = if total == 0 {
            0.0
        } else {
            self.narrow_peaks as f64 / total as f64
        };

        // Build per-chrom stats.
        let per_chrom: AHashMap<String, ChromPeakStats> = self
            .per_chrom
            .into_iter()
            .map(|(chrom, (n_peaks, open_bp, sig_sum, width_sum))| {
                let mean_signal = if n_peaks == 0 {
                    0.0
                } else {
                    sig_sum / n_peaks as f64
                };
                let mean_width = if n_peaks == 0 {
                    0.0
                } else {
                    width_sum as f64 / n_peaks as f64
                };
                (
                    chrom,
                    ChromPeakStats {
                        total_peaks: n_peaks,
                        total_open_bp: open_bp,
                        mean_signal,
                        mean_width,
                    },
                )
            })
            .collect();

        // Finalize top peaks: sort descending by signal and take top 100.
        self.top_peaks.sort_unstable_by(|a, b| {
            b.signal_value
                .partial_cmp(&a.signal_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        self.top_peaks.truncate(TOP_PEAKS_FINAL);

        Ok(AtacSummary {
            total_peaks: total,
            total_open_chromatin_bp: self.open_bp_total,
            mean_peak_width,
            mean_signal_value,
            median_signal_value,
            per_chrom,
            top_peaks: self.top_peaks,
            fraction_narrow_peaks,
        })
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_peak(chrom: &str, start: u64, end: u64, signal: f64) -> AtacPeak {
        AtacPeak {
            chrom: chrom.to_owned(),
            start,
            end,
            name: format!("peak_{}_{}", start, end),
            score: 500.0,
            strand: None,
            signal_value: signal,
            p_value_log10: 5.0,
            q_value_log10: 3.0,
            peak_offset: -1,
        }
    }

    #[test]
    fn test_empty_finalize() {
        let accum = AtacAccum::default();
        let summary = accum.finalize().expect("finalize");
        assert_eq!(summary.total_peaks, 0);
        assert_eq!(summary.total_open_chromatin_bp, 0);
        assert!((summary.mean_peak_width - 0.0).abs() < f64::EPSILON);
        assert!((summary.fraction_narrow_peaks - 0.0).abs() < f64::EPSILON);
        assert!(summary.top_peaks.is_empty());
    }

    #[test]
    fn test_single_peak() {
        let mut accum = AtacAccum::default();
        // width = 300, narrow (< 500)
        accum
            .process(&make_peak("chr1", 1000, 1300, 7.5))
            .expect("process");
        let summary = accum.finalize().expect("finalize");
        assert_eq!(summary.total_peaks, 1);
        assert_eq!(summary.total_open_chromatin_bp, 300);
        assert!((summary.mean_peak_width - 300.0).abs() < 1e-10);
        assert!((summary.mean_signal_value - 7.5).abs() < 1e-10);
        assert!((summary.median_signal_value - 7.5).abs() < 1e-10);
        assert!((summary.fraction_narrow_peaks - 1.0).abs() < f64::EPSILON);
        assert_eq!(summary.top_peaks.len(), 1);
    }

    #[test]
    fn test_per_chrom_aggregation() {
        let mut accum = AtacAccum::default();
        accum
            .process(&make_peak("chr1", 100, 600, 5.0))
            .expect("process"); // width 500, not narrow
        accum
            .process(&make_peak("chr1", 700, 900, 3.0))
            .expect("process"); // width 200, narrow
        accum
            .process(&make_peak("chr2", 0, 100, 9.0))
            .expect("process"); // width 100, narrow

        let summary = accum.finalize().expect("finalize");
        assert_eq!(summary.total_peaks, 3);

        let chr1 = summary.per_chrom.get("chr1").expect("chr1");
        assert_eq!(chr1.total_peaks, 2);
        assert_eq!(chr1.total_open_bp, 700); // 500 + 200

        let chr2 = summary.per_chrom.get("chr2").expect("chr2");
        assert_eq!(chr2.total_peaks, 1);
        assert_eq!(chr2.total_open_bp, 100);

        // 2 narrow out of 3 total
        assert!((summary.fraction_narrow_peaks - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_top_peaks_ordering() {
        let mut accum = AtacAccum::default();
        for i in 0..150u64 {
            accum
                .process(&make_peak("chr1", i * 1000, i * 1000 + 200, i as f64))
                .expect("process");
        }
        let summary = accum.finalize().expect("finalize");
        assert_eq!(summary.top_peaks.len(), TOP_PEAKS_FINAL);
        // Highest signal should be first
        assert!(summary.top_peaks[0].signal_value >= summary.top_peaks[1].signal_value);
        // All top peaks should have signal >= the 100th-highest
        let min_top = summary
            .top_peaks
            .last()
            .map(|p| p.signal_value)
            .unwrap_or(0.0);
        // signals were 0..149; top 100 should be 50..149
        assert!((min_top - 50.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_merge_equivalence() {
        // Processing all records in one accumulator must equal merging two halves.
        let peaks: Vec<AtacPeak> = (0..50u64)
            .map(|i| make_peak("chr1", i * 100, i * 100 + 50, i as f64))
            .collect();

        let mut single = AtacAccum::default();
        for p in &peaks {
            single.process(p).expect("process");
        }
        let s1 = single.finalize().expect("finalize");

        let mut left = AtacAccum::default();
        let mut right = AtacAccum::default();
        for p in &peaks[..25] {
            left.process(p).expect("process");
        }
        for p in &peaks[25..] {
            right.process(p).expect("process");
        }
        left.merge(right);
        let s2 = left.finalize().expect("finalize");

        assert_eq!(s1.total_peaks, s2.total_peaks);
        assert_eq!(s1.total_open_chromatin_bp, s2.total_open_chromatin_bp);
        assert!((s1.mean_signal_value - s2.mean_signal_value).abs() < 1e-10);
        assert!((s1.median_signal_value - s2.median_signal_value).abs() < 1e-10);
    }
}
