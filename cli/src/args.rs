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
    // ── Primary inputs ────────────────────────────────────────────────────────

    /// VCF file for genomics variant analysis
    #[arg(long, value_name = "FILE")]
    pub genomics: Option<PathBuf>,

    /// Expression matrix TSV or raw-count matrix (genes × samples, first row = sample names)
    #[arg(long, value_name = "FILE")]
    pub transcriptomics: Option<PathBuf>,

    /// BED methylation file for epigenomics analysis (ENCODE bisulfite or 4-column format)
    #[arg(long, value_name = "FILE")]
    pub epigenomics: Option<PathBuf>,

    /// ATAC-seq peaks in ENCODE narrowPeak (BED6+4) format
    #[arg(long, value_name = "FILE")]
    pub atac: Option<PathBuf>,

    /// VCF with CN INFO field for copy-number variation analysis
    #[arg(long, value_name = "FILE")]
    pub cnv: Option<PathBuf>,

    /// Optional FASTQ input for sequence-level QC
    #[arg(long, value_name = "FILE")]
    pub fastq: Option<PathBuf>,

    // ── Comparison mode (tumor-vs-normal / treatment-vs-control) ─────────────

    /// Control/normal VCF for comparison mode (enables --compare-* flags)
    #[arg(long, value_name = "FILE", requires = "compare_transcriptomics")]
    pub compare_genomics: Option<PathBuf>,

    /// Control/normal expression TSV for comparison mode
    #[arg(long, value_name = "FILE")]
    pub compare_transcriptomics: Option<PathBuf>,

    /// Control/normal methylation BED for comparison mode
    #[arg(long, value_name = "FILE")]
    pub compare_epigenomics: Option<PathBuf>,

    /// Control ATAC-seq narrowPeak for comparison mode
    #[arg(long, value_name = "FILE")]
    pub compare_atac: Option<PathBuf>,

    // ── Enrichment ────────────────────────────────────────────────────────────

    /// GMT pathway file for custom gene set enrichment (name, desc, genes…)
    #[arg(long, value_name = "FILE")]
    pub gmt: Option<PathBuf>,

    // ── Normalization ─────────────────────────────────────────────────────────

    /// Treat --transcriptomics as raw counts; apply DESeq2 size-factor normalization
    #[arg(long, default_value_t = false)]
    pub raw_counts: bool,

    // ── Config ────────────────────────────────────────────────────────────────

    /// TOML configuration file (thresholds, output settings, performance tuning)
    #[arg(long, value_name = "FILE")]
    pub config: Option<PathBuf>,

    /// Print default configuration as TOML and exit
    #[arg(long, default_value_t = false)]
    pub dump_config: bool,

    // ── Output ────────────────────────────────────────────────────────────────

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

impl Cli {
    /// Returns true if all three compare inputs are provided.
    pub fn has_compare(&self) -> bool {
        self.compare_genomics.is_some()
            && self.compare_transcriptomics.is_some()
            && self.compare_epigenomics.is_some()
    }
}
