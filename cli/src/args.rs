use std::path::PathBuf;

use clap::Parser;

/// BioMultiOmics — parallel multi-omics analysis: VCF + TSV + BED → HTML/JSON report.
///
/// Ingests genomics (VCF), transcriptomics (expression TSV), and epigenomics
/// (methylation BED) data simultaneously and produces integrated biological insights.
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

    /// Expression matrix TSV file (genes × samples, first row = header with sample names)
    #[arg(long, value_name = "FILE")]
    pub transcriptomics: PathBuf,

    /// BED methylation file for epigenomics analysis (ENCODE bisulfite or 4-column format)
    #[arg(long, value_name = "FILE")]
    pub epigenomics: PathBuf,

    /// Optional FASTQ input for sequence-level QC (paired or single-end)
    #[arg(long, value_name = "FILE")]
    pub fastq: Option<PathBuf>,

    /// Output directory (created if it does not exist)
    #[arg(long, value_name = "DIR", default_value = "./bioomics_out")]
    pub output: PathBuf,

    /// Number of parallel worker threads (default: all logical cores)
    #[arg(long, value_name = "N")]
    pub threads: Option<usize>,

    /// Skip the ML integration layer (PCA and Pearson correlation)
    #[arg(long, default_value_t = false)]
    pub no_ml: bool,

    /// Compare against a second set of input files (JSON with genomics/transcriptomics/epigenomics keys)
    #[arg(long, value_name = "FILE")]
    pub compare: Option<PathBuf>,

    /// Emit JSON output only — no TUI, no HTML report
    #[arg(long, default_value_t = false)]
    pub json: bool,
}
