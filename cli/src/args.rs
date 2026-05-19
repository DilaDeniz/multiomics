use std::path::PathBuf;

use clap::Parser;

use crate::config;

/// Multiomics — parallel multi-omics analysis: VCF + TSV + BED → HTML/JSON report.
#[derive(Parser, Debug, Clone)]
#[command(
    name = "multiomics",
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

    // ── Proteomics ────────────────────────────────────────────────────────────
    /// One or more mzML files for proteomics database search (requires --fasta).
    /// Pass multiple files to run a multi-file experiment-level search.
    #[arg(long, value_name = "FILE", num_args = 1..)]
    pub proteomics: Vec<PathBuf>,

    /// Directory containing *.mzML files (alternative to listing files individually)
    #[arg(long, value_name = "DIR")]
    pub proteomics_dir: Option<PathBuf>,

    /// Protein database FASTA for proteomics search
    #[arg(long, value_name = "FILE")]
    pub fasta: Option<PathBuf>,

    /// FDR threshold for proteomics reporting [default: 0.01]
    #[arg(long, default_value_t = 0.01f64)]
    pub proteomics_fdr: f64,

    // ── Somatic variant calling ───────────────────────────────────────────────
    /// Tumor BAM for somatic variant calling (requires --normal-bam)
    #[arg(long, value_name = "FILE")]
    pub tumor_bam: Option<PathBuf>,

    /// Matched normal BAM for somatic variant calling
    #[arg(long, value_name = "FILE")]
    pub normal_bam: Option<PathBuf>,

    /// Minimum tumor log-odds score for somatic calls (default: 6.3)
    #[arg(long, default_value_t = 6.3)]
    pub somatic_min_lod: f64,

    // ── Reference-guided alignment ────────────────────────────────────────────
    /// Reference FASTA for aligning reads (enables --fastq processing)
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,

    // ── Single-cell ───────────────────────────────────────────────────────────
    /// 10x Genomics MEX directory for single-cell analysis
    #[arg(long, value_name = "DIR")]
    pub scrna: Option<PathBuf>,

    /// Number of UMAP neighbors (default: 15)
    #[arg(long, default_value_t = 15usize)]
    pub umap_neighbors: usize,

    // ── Gene quantification ───────────────────────────────────────────────────
    /// BAM file for gene quantification (requires --gtf)
    #[arg(long, value_name = "FILE")]
    pub bam: Option<PathBuf>,

    /// GTF/GFF3 annotation for gene quantification
    #[arg(long, value_name = "FILE")]
    pub gtf: Option<PathBuf>,

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

    /// Pre-configured threshold set for common use cases (cancer, plant, rna-seq, wgbs, atac, clinical)
    #[arg(long, value_enum)]
    pub preset: Option<config::Preset>,

    /// List all available presets and exit
    #[arg(long, default_value_t = false)]
    pub list_presets: bool,

    // ── Output ────────────────────────────────────────────────────────────────
    /// Output directory (created if it does not exist)
    #[arg(long, value_name = "DIR", default_value = "./multiomics_out")]
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
