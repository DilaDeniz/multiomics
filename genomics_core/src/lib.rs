pub mod accum;
pub mod cnv;
pub mod fastq;
pub mod types;
pub mod vcf;

pub use accum::GenomicsAccum;
pub use cnv::{parse_cnv_vcf, summarize_cnv, CnvRecord, CnvSummary};
pub use fastq::{parse_fastq, FastqSummary};
pub use types::{ChromDensity, GenomicsSummary, TiTvClass, VariantRecord};
pub use vcf::parse_vcf;

use anyhow::Result;
use crossbeam_channel::Sender;
use std::path::Path;

use biomics_core::{parallel_fold, ProgressEvent};

/// Run the full genomics analysis pipeline on a VCF file.
pub fn analyze_vcf(path: &Path, progress_tx: Option<&Sender<ProgressEvent>>) -> Result<GenomicsSummary> {
    let records = parse_vcf(path)?;
    parallel_fold::<GenomicsAccum>(&records, "genomics", progress_tx)
}
