use ahash::AHashMap;
use serde::{Deserialize, Serialize};

/// One record from an ENCODE narrowPeak file (BED6+4).
///
/// Format:
/// ```text
/// chrom  start  end  name  score  strand  signalValue  pValue  qValue  peak
/// ```
/// - `score`: 0–1000 integer score (stored as f64 for uniformity).
/// - `strand`: `'+'`, `'-'`, or `'.'`; the dot case is represented as `None`.
/// - `peak_offset`: byte offset from `start` to the summit; `-1` in the file
///   encodes "not determined" and is stored as `-1i64` (sentinel).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtacPeak {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    /// Score in [0, 1000].
    pub score: f64,
    /// `None` when the strand field is `'.'`.
    pub strand: Option<char>,
    /// Fold enrichment or average signal over the peak.
    pub signal_value: f64,
    /// –log10(p-value); `-1.0` means not computed.
    pub p_value_log10: f64,
    /// –log10(q-value / FDR); `-1.0` means not computed.
    pub q_value_log10: f64,
    /// Offset (in bp) from `start` to the peak summit. `-1` = not determined.
    pub peak_offset: i64,
}

impl AtacPeak {
    /// Peak width in base pairs.
    #[inline]
    pub fn width(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    /// Whether this is a "narrow" peak (width < 500 bp).
    #[inline]
    pub fn is_narrow(&self) -> bool {
        self.width() < 500
    }
}

/// Per-chromosome ATAC-seq statistics.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ChromPeakStats {
    /// Number of peaks on this chromosome.
    pub total_peaks: u64,
    /// Sum of peak widths — the total base pairs of open chromatin.
    pub total_open_bp: u64,
    /// Mean signal value across all peaks on this chromosome.
    pub mean_signal: f64,
    /// Mean peak width (bp) on this chromosome.
    pub mean_width: f64,
}

/// Final report produced by [`crate::accum::AtacAccum::finalize`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtacSummary {
    /// Total number of peaks across all chromosomes.
    pub total_peaks: u64,
    /// Total base pairs of open chromatin (sum of all peak widths).
    pub total_open_chromatin_bp: u64,
    /// Mean peak width across all peaks.
    pub mean_peak_width: f64,
    /// Mean signal value across all peaks.
    pub mean_signal_value: f64,
    /// Median signal value across all peaks (exact, via sorted vector).
    pub median_signal_value: f64,
    /// Per-chromosome breakdown.
    pub per_chrom: AHashMap<String, ChromPeakStats>,
    /// Top 100 peaks by signal value (descending).
    pub top_peaks: Vec<AtacPeak>,
    /// Fraction of peaks with width < 500 bp.
    pub fraction_narrow_peaks: f64,
}
