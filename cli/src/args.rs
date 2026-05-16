use std::path::PathBuf;

use clap::Parser;

/// BioMultiOmics — parallel multi-omics analysis: VCF + TSV + BED → HTML/JSON report.
#[derive(Parser, Debug, Clone)]
#[command(
    name = "bioomics",
    version,
    about = "Parallel multi-omics analysis: VCF + TSV + BED → integrated HTML/JSON report",
    long_about = None,
)]
pub struct Cli {
    /// VCF file for genomics variant analysis
    #[arg(long, value_name = "FILE")]
    pub genomics: PathBuf,

    /// Expression matrix TSV or raw-count matrix (genes × samples, first row = sample names)
    #[arg(long, value_name = "FILE")]
    pub transcriptomics: PathBuf,

    /// BED methylation file for epigenomics analysis (ENCODE bisulfite or 4-column format)
    #[arg(long, value_name = "FILE")]
    pub epigenomics: PathBuf,

    /// ATAC-seq peaks in ENCODE narrowPeak (BED6+4) format
    #[arg(long, value_name = "FILE")]
    pub atac: Option<PathBuf>,

    /// VCF with CN INFO field for copy-number variation analysis
    #[arg(long, value_name = "FILE")]
    pub cnv: Option<PathBuf>,

    /// Optional FASTQ input for sequence-level QC (paired or single-end)
    #[arg(long, value_name = "FILE")]
    pub fastq: Option<PathBuf>,

    /// GMT pathway file for custom gene set enrichment (tab-delimited: name, desc, genes...)
    #[arg(long, value_name = "FILE")]
    pub gmt: Option<PathBuf>,

    /// Treat --transcriptomics input as raw counts and apply DESeq2 size-factor normalization
    #[arg(long, default_value_t = false)]
    pub raw_counts: bool,

    /// Output directory (created if it does not exist)
    #[arg(long, value_name = "DIR", default_value = "./bioomics_out")]
    pub output: PathBuf,

    /// Number of parallel worker threads (default: all logical cores)
    #[arg(long, value_name = "N")]
    pub threads: Option<usize>,

    /// Skip the ML integration layer (PCA and correlation matrix)
    #[arg(long, default_value_t = false)]
    pub no_ml: bool,

    /// Emit JSON output only — no TUI, no HTML report
    #[arg(long, default_value_t = false)]
    pub json: bool,
}
