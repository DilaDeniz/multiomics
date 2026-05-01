pub mod accum;
pub mod diffexpr;
pub mod tsv;
pub mod types;

pub use accum::TranscriptomicsAccum;
pub use diffexpr::{significant_de_genes, two_sample_diff_expr};
pub use tsv::parse_tsv;
pub use types::{DiffExprResult, GeneRecord, GeneStats, TranscriptomicsSummary};

use anyhow::Result;
use crossbeam_channel::Sender;
use std::path::Path;

use biomics_core::{parallel_fold, ProgressEvent};

/// Run the full transcriptomics analysis pipeline on an expression matrix TSV.
///
/// Parses the TSV, runs a parallel fold, injects sample names and (if 2 samples)
/// differential expression into the returned summary.
pub fn analyze_tsv(
    path: &Path,
    progress_tx: Option<&Sender<ProgressEvent>>,
) -> Result<TranscriptomicsSummary> {
    let (records, sample_names) = parse_tsv(path)?;

    let n_samples = sample_names.len();
    let mut summary =
        parallel_fold::<TranscriptomicsAccum>(&records, "transcriptomics", progress_tx)?;

    summary.sample_names = sample_names;
    summary.sample_count = n_samples;

    if n_samples == 2 {
        let de = two_sample_diff_expr(&records);
        summary.diff_expr = Some(de);
    }

    Ok(summary)
}
