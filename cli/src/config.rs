//! TOML configuration file support.
//!
//! Load with `--config bioomics.toml`. CLI flags override config values.
//! Run `bioomics --dump-config` to print defaults.

use std::path::Path;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

/// Master configuration. All fields have sensible defaults via `Default`.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct BioomicsConfig {
    pub genomics: GenomicsConfig,
    pub transcriptomics: TranscriptomicsConfig,
    pub epigenomics: EpigenomicsConfig,
    pub atac: AtacConfig,
    pub integration: IntegrationConfig,
    pub output: OutputConfig,
    pub performance: PerformanceConfig,
    pub compare: CompareConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct GenomicsConfig {
    /// Minimum QUAL score to consider a variant high-impact.
    pub high_impact_qual: f64,
    /// Minimum allele frequency to include in AF histogram.
    pub min_af: f64,
    /// Ti/Tv ratio below this triggers a WARNING insight.
    pub titv_warn_below: f64,
    /// Ti/Tv ratio above this triggers a WARNING insight.
    pub titv_warn_above: f64,
    /// Maximum high-impact variants to keep in memory.
    pub max_high_impact: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct TranscriptomicsConfig {
    /// TPM threshold for "expressed" classification.
    pub expressed_tpm: f64,
    /// Adjusted p-value cutoff for significant DE.
    pub padj_threshold: f64,
    /// Log₂ fold-change threshold for DE reporting.
    pub log2fc_threshold: f64,
    /// Maximum DE genes to store in the summary.
    pub max_de_genes: usize,
    /// Top N expressed genes to report.
    pub top_n_expressed: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct EpigenomicsConfig {
    /// Global methylation below this triggers a CRITICAL insight.
    pub global_meth_crit_below: f64,
    /// Global methylation below this triggers a WARNING insight.
    pub global_meth_warn_below: f64,
    /// Minimum CpG island length in bp (Gardiner-Garden criterion).
    pub cpg_island_min_len: u64,
    /// Minimum GC fraction for CpG island detection.
    pub cpg_island_min_gc: f64,
    /// Minimum CpO/E ratio for CpG island detection.
    pub cpg_island_min_cpoe: f64,
    /// Methylation above this is "hypermethylated" (per-region).
    pub hypermeth_threshold: f64,
    /// Methylation below this is "hypomethylated" (per-region).
    pub hypometh_threshold: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct AtacConfig {
    /// Minimum signal value to include a peak.
    pub min_signal: f64,
    /// Minimum q-value (−log₁₀) for a peak to be considered significant.
    pub min_qvalue: f64,
    /// Top N peaks to keep per chromosome.
    pub top_n_peaks: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct IntegrationConfig {
    /// Minimum gene overlap for pathway enrichment.
    pub min_pathway_overlap: usize,
    /// Maximum pathways to report.
    pub top_n_pathways: usize,
    /// Number of GSEA permutations.
    pub gsea_n_permutations: usize,
    /// Minimum gene set size for GSEA.
    pub gsea_min_size: usize,
    /// Maximum gene set size for GSEA.
    pub gsea_max_size: usize,
    /// Absolute correlation threshold to trigger a cross-modality insight.
    pub correlation_insight_threshold: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct OutputConfig {
    /// Title displayed in the HTML report header.
    pub report_title: String,
    /// Include raw per-record data arrays in JSON output.
    pub include_raw_data: bool,
    /// Colour scheme: "dark" (default) or "light".
    pub color_scheme: String,
    /// Maximum insights to show in TUI panel.
    pub max_tui_insights: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PerformanceConfig {
    /// Records per rayon chunk. Tune for your L2 cache size.
    pub chunk_size: usize,
    /// Stack size per rayon worker thread in bytes.
    pub thread_stack_bytes: usize,
    /// Pre-allocate this many records when parsing (0 = auto).
    pub preallocate_records: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CompareConfig {
    /// Label for the case/treatment group in comparison reports.
    pub case_label: String,
    /// Label for the control/normal group in comparison reports.
    pub control_label: String,
    /// Log₂FC threshold for variant-level comparison insights.
    pub variant_fc_threshold: f64,
    /// Methylation difference (absolute) to flag as differentially methylated.
    pub delta_meth_threshold: f64,
}


impl Default for GenomicsConfig {
    fn default() -> Self {
        Self {
            high_impact_qual: 30.0,
            min_af: 0.0,
            titv_warn_below: 1.8,
            titv_warn_above: 3.0,
            max_high_impact: 10_000,
        }
    }
}

impl Default for TranscriptomicsConfig {
    fn default() -> Self {
        Self {
            expressed_tpm: 1.0,
            padj_threshold: 0.05,
            log2fc_threshold: 1.0,
            max_de_genes: 50_000,
            top_n_expressed: 100,
        }
    }
}

impl Default for EpigenomicsConfig {
    fn default() -> Self {
        Self {
            global_meth_crit_below: 40.0,
            global_meth_warn_below: 55.0,
            cpg_island_min_len: 200,
            cpg_island_min_gc: 0.5,
            cpg_island_min_cpoe: 0.6,
            hypermeth_threshold: 80.0,
            hypometh_threshold: 20.0,
        }
    }
}

impl Default for AtacConfig {
    fn default() -> Self {
        Self {
            min_signal: 0.0,
            min_qvalue: 0.0,
            top_n_peaks: 200,
        }
    }
}

impl Default for IntegrationConfig {
    fn default() -> Self {
        Self {
            min_pathway_overlap: 1,
            top_n_pathways: 20,
            gsea_n_permutations: 1000,
            gsea_min_size: 10,
            gsea_max_size: 500,
            correlation_insight_threshold: 0.7,
        }
    }
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            report_title: "BioMultiOmics Report".into(),
            include_raw_data: false,
            color_scheme: "dark".into(),
            max_tui_insights: 12,
        }
    }
}

impl Default for PerformanceConfig {
    fn default() -> Self {
        Self {
            chunk_size: 65_536,
            thread_stack_bytes: 8 * 1024 * 1024,
            preallocate_records: 0,
        }
    }
}

impl Default for CompareConfig {
    fn default() -> Self {
        Self {
            case_label: "case".into(),
            control_label: "control".into(),
            variant_fc_threshold: 1.5,
            delta_meth_threshold: 15.0,
        }
    }
}

/// Load config from a TOML file, merging with defaults.
pub fn load_config(path: &Path) -> Result<BioomicsConfig> {
    let raw = std::fs::read_to_string(path)
        .with_context(|| format!("Cannot read config '{}'", path.display()))?;
    let cfg: BioomicsConfig = toml::from_str(&raw)
        .with_context(|| format!("Invalid TOML in '{}'", path.display()))?;
    Ok(cfg)
}

/// Serialize the default config to TOML for `--dump-config`.
pub fn dump_default_config() -> String {
    toml::to_string_pretty(&BioomicsConfig::default())
        .expect("default config is always serializable")
}
