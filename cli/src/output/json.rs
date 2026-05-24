use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use serde::Serialize;

use epigenomics_core::EpigenomicsSummary;
use genomics_core::GenomicsSummary;
use integration_layer::{EnrichmentResult, GeneParadox, GeneRegulatoryProfile, Insight, IntegrationSummary};
use transcriptomics_core::TranscriptomicsSummary;

/// MultiQC-compatible general stats for the top-level table.
#[derive(Debug, Serialize)]
pub struct GeneralStats {
    pub total_variants: u64,
    pub snp_count: u64,
    pub indel_count: u64,
    pub titv_ratio: f64,
    pub total_genes: u64,
    pub expressed_genes: u64,
    pub total_sites: u64,
    pub global_methylation_pct: f64,
    pub cpg_islands_detected: usize,
}

/// Per-column display metadata for MultiQC general stats table.
#[derive(Debug, Serialize)]
pub struct ColumnMeta {
    pub title: &'static str,
    pub format: &'static str,
    pub scale: &'static str,
}

#[derive(Debug, Serialize)]
pub struct JsonGenomicsSection {
    pub total_variants: u64,
    pub snp_count: u64,
    pub indel_count: u64,
    pub titv_ratio: f64,
    pub high_impact_count: usize,
    pub high_impact_genes: Vec<String>,
    pub unique_positions: u64,
    pub af_histogram: Vec<u64>,
    pub per_chrom: HashMap<String, ChromStats>,
}

#[derive(Debug, Serialize)]
pub struct ChromStats {
    pub total: u64,
    pub snps: u64,
    pub indels: u64,
}

#[derive(Debug, Serialize)]
pub struct JsonTranscriptomicsSection {
    pub total_genes: u64,
    pub expressed_genes: u64,
    pub sample_count: usize,
    pub sample_names: Vec<String>,
    pub top_expressed: Vec<(String, f64)>,
    pub diff_expr_count: Option<usize>,
}

#[derive(Debug, Serialize)]
pub struct JsonEpigenomicsSection {
    pub total_sites: u64,
    pub global_methylation_pct: f64,
    pub cpg_islands_detected: usize,
    pub hypermethylated_regions: usize,
    pub hypomethylated_regions: usize,
    pub per_chrom_methylation: HashMap<String, f64>,
}

#[derive(Debug, Serialize)]
pub struct JsonIntegrationSection {
    pub correlation_matrix: Vec<Vec<f64>>,
    pub pca_points: Vec<[f64; 2]>,
    pub pca_explained_variance: Vec<f64>,
    pub top_pathways: Vec<EnrichmentResult>,
    pub insights: Vec<Insight>,
    pub paradoxes: Vec<GeneParadox>,
    pub gene_states: Vec<GeneRegulatoryProfile>,
}

#[derive(Debug, Serialize)]
pub struct ReportMetadata {
    pub tool: &'static str,
    pub version: &'static str,
    pub generated_at: DateTime<Utc>,
    pub threads_used: usize,
    pub elapsed_seconds: u64,
}

/// Root MultiQC-compatible output object.
#[derive(Debug, Serialize)]
pub struct MultiQcOutput {
    pub report_general_stats_data: Vec<HashMap<String, GeneralStats>>,
    pub report_general_stats_headers: HashMap<String, ColumnMeta>,
    pub multiomics_genomics: JsonGenomicsSection,
    pub multiomics_transcriptomics: JsonTranscriptomicsSection,
    pub multiomics_epigenomics: JsonEpigenomicsSection,
    pub multiomics_integration: JsonIntegrationSection,
    pub metadata: ReportMetadata,
}

/// Build the full `MultiQcOutput` from all modality summaries.
pub fn build_multiqc_output(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    integration: &IntegrationSummary,
    threads_used: usize,
    elapsed_secs: u64,
) -> MultiQcOutput {
    let general = GeneralStats {
        total_variants: genomics.total_variants,
        snp_count: genomics.snp_count,
        indel_count: genomics.indel_count,
        titv_ratio: genomics.titv_ratio,
        total_genes: transcr.total_genes,
        expressed_genes: transcr.expressed_genes,
        total_sites: epigen.total_sites,
        global_methylation_pct: epigen.global_methylation_pct,
        cpg_islands_detected: epigen.cpg_islands.len(),
    };

    let mut general_map = HashMap::new();
    general_map.insert("multiomics".to_string(), general);

    let mut headers = HashMap::new();
    headers.insert(
        "total_variants".to_string(),
        ColumnMeta {
            title: "Variants",
            format: "{:.0}",
            scale: "Blues",
        },
    );
    headers.insert(
        "titv_ratio".to_string(),
        ColumnMeta {
            title: "Ti/Tv",
            format: "{:.2}",
            scale: "RdYlGn",
        },
    );
    headers.insert(
        "expressed_genes".to_string(),
        ColumnMeta {
            title: "Expressed Genes",
            format: "{:.0}",
            scale: "Greens",
        },
    );
    headers.insert(
        "global_methylation_pct".to_string(),
        ColumnMeta {
            title: "Global Meth %",
            format: "{:.1}",
            scale: "Oranges",
        },
    );

    let per_chrom: HashMap<String, ChromStats> = genomics
        .per_chrom
        .iter()
        .map(|(chrom, d)| {
            (
                chrom.clone(),
                ChromStats {
                    total: d.total,
                    snps: d.snps,
                    indels: d.indels,
                },
            )
        })
        .collect();

    let per_chrom_methylation: HashMap<String, f64> = epigen
        .per_chrom
        .iter()
        .map(|(chrom, cm)| (chrom.clone(), cm.mean_methylation))
        .collect();

    MultiQcOutput {
        report_general_stats_data: vec![general_map],
        report_general_stats_headers: headers,
        multiomics_genomics: JsonGenomicsSection {
            total_variants: genomics.total_variants,
            snp_count: genomics.snp_count,
            indel_count: genomics.indel_count,
            titv_ratio: genomics.titv_ratio,
            high_impact_count: genomics.high_impact.len(),
            high_impact_genes: genomics.high_impact_genes.clone(),
            unique_positions: genomics.unique_positions,
            af_histogram: genomics.af_histogram.clone(),
            per_chrom,
        },
        multiomics_transcriptomics: JsonTranscriptomicsSection {
            total_genes: transcr.total_genes,
            expressed_genes: transcr.expressed_genes,
            sample_count: transcr.sample_count,
            sample_names: transcr.sample_names.clone(),
            top_expressed: transcr.top_100_expressed.clone(),
            diff_expr_count: transcr.diff_expr.as_ref().map(|de| de.len()),
        },
        multiomics_epigenomics: JsonEpigenomicsSection {
            total_sites: epigen.total_sites,
            global_methylation_pct: epigen.global_methylation_pct,
            cpg_islands_detected: epigen.cpg_islands.len(),
            hypermethylated_regions: epigen.hypermethylated.len(),
            hypomethylated_regions: epigen.hypomethylated.len(),
            per_chrom_methylation,
        },
        multiomics_integration: JsonIntegrationSection {
            correlation_matrix: integration.correlation_matrix.clone(),
            pca_points: integration.pca.points.clone(),
            pca_explained_variance: integration.pca.explained_variance_ratio.clone(),
            top_pathways: integration.top_pathways.clone(),
            insights: integration.insights.clone(),
            paradoxes: integration.paradoxes.clone(),
            gene_states: integration.gene_states.clone(),
        },
        metadata: ReportMetadata {
            tool: "multiomics",
            version: env!("CARGO_PKG_VERSION"),
            generated_at: Utc::now(),
            threads_used,
            elapsed_seconds: elapsed_secs,
        },
    }
}

/// Write a `MultiQcOutput` to `{output_dir}/multiqc_multiomics.json`.
pub fn write_json(output: &MultiQcOutput, output_dir: &Path) -> Result<()> {
    let path = output_dir.join("multiqc_multiomics.json");
    let file = std::fs::File::create(&path)
        .with_context(|| format!("Cannot create JSON output '{}'", path.display()))?;
    serde_json::to_writer_pretty(file, output)
        .with_context(|| format!("Cannot serialize JSON output to '{}'", path.display()))?;
    log::info!("JSON report written to '{}'", path.display());
    Ok(())
}
