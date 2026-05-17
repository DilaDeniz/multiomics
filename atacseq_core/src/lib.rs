pub mod accum;
pub mod narrowpeak;
pub mod types;

pub use accum::AtacAccum;
pub use narrowpeak::parse_narrowpeak;
pub use types::{AtacPeak, AtacSummary, ChromPeakStats};

use anyhow::Result;
use biomics_core::{parallel_fold, ProgressEvent};
use crossbeam_channel::Sender;
use std::path::Path;

/// Parse a narrowPeak file and compute open-chromatin statistics in parallel.
///
/// All peaks are loaded from `path`, then distributed across rayon worker threads
/// via [`parallel_fold`]. Progress events are sent on `progress_tx` when provided.
///
/// # Errors
///
/// Returns an error if the file cannot be opened, memory-mapped, or if the
/// parallel fold encounters an unrecoverable internal error (individual malformed
/// lines are logged and skipped, not propagated as errors).
pub fn analyze_narrowpeak(
    path: &Path,
    progress_tx: Option<&Sender<ProgressEvent>>,
) -> Result<AtacSummary> {
    let records = parse_narrowpeak(path)?;
    parallel_fold::<AtacAccum>(&records, "atac", progress_tx)
}
