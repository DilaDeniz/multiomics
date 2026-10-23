pub mod accum;
pub mod deseq2;
pub mod diffexpr;
pub mod quantify;
pub mod tsv;
pub mod types;

pub use accum::TranscriptomicsAccum;
pub use deseq2::{deseq2_differential_expression, normalize_counts, NormalizedMatrix, SizeFactors};
pub use diffexpr::{differential_expression, significant_de_genes};
pub use quantify::{count_reads, parse_gtf, quantify, summarize_quant, GeneInterval, QuantSummary};
pub use tsv::parse_tsv;
pub use types::{DiffExprResult, GeneRecord, GeneStats, TranscriptomicsSummary};

use anyhow::{Context, Result};
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

/// Like [`analyze_tsv`] for a matrix of raw integer read counts: differential
/// expression is computed with the DESeq2 negative-binomial Wald test
/// ([`deseq2::run_deseq2`]) instead of the Welch screen.
///
/// Samples are split by position, as in [`differential_expression`]: the first
/// ⌈n/2⌉ columns form group 1 (reference) and the rest group 2. At least three
/// samples are required. `mean_s1`/`mean_s2` hold mean normalised counts.
pub fn analyze_count_tsv(
    path: &Path,
    progress_tx: Option<&Sender<ProgressEvent>>,
) -> Result<TranscriptomicsSummary> {
    let (records, sample_names) = parse_tsv(path)?;
    let mut summary =
        parallel_fold::<TranscriptomicsAccum>(&records, "transcriptomics", progress_tx)?;
    summary.sample_count = sample_names.len();
    let de = deseq2_de_results(&records, &sample_names).with_context(|| {
        format!(
            "DESeq2 on '{}' failed — raw integer read counts are required \
             (omit --raw-counts for TPM/FPKM input)",
            path.display()
        )
    })?;
    summary.diff_expr = Some(de);
    summary.sample_names = sample_names;
    Ok(summary)
}

/// Run DESeq2 on parsed count records with the positional two-group split and
/// convert the results to [`DiffExprResult`]s sorted by padj (NaN last).
pub fn deseq2_de_results(
    records: &[GeneRecord],
    sample_names: &[String],
) -> Result<Vec<DiffExprResult>> {
    let n = sample_names.len();
    if n < 3 {
        anyhow::bail!("DESeq2 needs at least 3 samples (found {n})");
    }
    let split = n.div_ceil(2);
    let group0: Vec<usize> = (0..split).collect();
    let group1: Vec<usize> = (split..n).collect();
    let ids: Vec<String> = records.iter().map(|r| r.gene_id.clone()).collect();
    let counts: Vec<Vec<f64>> = records.iter().map(|r| r.samples.clone()).collect();
    let analysis = deseq2::run_deseq2(&counts, &ids, &group0, &group1, sample_names)?;
    log::info!(
        "DESeq2: {} genes, {} with padj < 0.05, {} removed by independent filtering \
         (mean count < {:.3}), trend {:?}",
        analysis.results.len(),
        analysis.n_significant,
        analysis.n_filtered,
        analysis.filter_threshold,
        analysis.dispersion_trend.kind
    );
    let mut results: Vec<DiffExprResult> = analysis
        .results
        .into_iter()
        .map(|r| DiffExprResult {
            gene_id: r.gene_id,
            log2_fold_change: r.log2_fold_change,
            mean_s1: r.mean_group0,
            mean_s2: r.mean_group1,
            p_value: r.p_value,
            padj: r.padj,
        })
        .collect();
    results.sort_by(|a, b| match (a.padj.is_nan(), b.padj.is_nan()) {
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => a
            .padj
            .partial_cmp(&b.padj)
            .unwrap_or(std::cmp::Ordering::Equal),
    });
    Ok(results)
}
