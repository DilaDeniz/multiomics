pub mod accum;
pub mod bed;
pub mod cpg;
#[cfg(feature = "longread")]
pub mod longread;
pub mod types;

pub use accum::EpigenomicsAccum;
pub use bed::parse_bed;
pub use cpg::detect_cpg_islands;
#[cfg(feature = "longread")]
pub use longread::{longread_to_methylation_records, parse_longread_methylation, LongReadMethCall};
pub use types::{
    ChromMethylation, CpGIsland, EpigenomicsSummary, MethylationRecord, MethylationRegion,
    RegionKind,
};

use anyhow::Result;
use crossbeam_channel::Sender;
use std::path::Path;

use biomics_core::{parallel_fold, ProgressEvent};

/// Run the full epigenomics analysis pipeline on a BED methylation file.
///
/// Parses the BED file, distributes records across the rayon thread pool,
/// and returns a finalized `EpigenomicsSummary` including CpG island detection
/// and hypermethylated/hypomethylated region lists.
pub fn analyze_bed(
    path: &Path,
    progress_tx: Option<&Sender<ProgressEvent>>,
) -> Result<EpigenomicsSummary> {
    let records = parse_bed(path)?;
    parallel_fold::<EpigenomicsAccum>(&records, "epigenomics", progress_tx)
}
