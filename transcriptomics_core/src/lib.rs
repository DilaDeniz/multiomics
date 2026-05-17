pub mod accum;
pub mod deseq2;
pub mod diffexpr;
pub mod tsv;
pub mod types;

pub use accum::TranscriptomicsAccum;
pub use deseq2::{deseq2_differential_expression, normalize_counts, NormalizedMatrix, SizeFactors};
pub use diffexpr::{differential_expression, significant_de_genes};
pub use tsv::parse_tsv;
pub use types::{DiffExprResult, GeneRecord, GeneStats, TranscriptomicsSummary};

use anyhow::Result;
use crossbeam_channel::Sender;
use std::path::Path;

use biomics_core::{parallel_fold, ProgressEvent};

/// Run the full transcriptomics analysis pipeline on an expression matrix TSV.
///
/// Parses the TSV, runs a parallel fold to accumulate per-gene statistics,
/// injects sample names, then runs differential expression for any n ≥ 2:
/// - n = 2: log₂FC only (no valid within-group variance for t-test)
/// - n ≥ 4: Welch t-test + BH FDR correction
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

    if n_samples >= 2 {
        let de = differential_expression(&records);
        summary.diff_expr = Some(de);
    }

    Ok(summary)
}
