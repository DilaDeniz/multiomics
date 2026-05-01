use serde::{Deserialize, Serialize};

use genomics_core::GenomicsSummary;
use transcriptomics_core::TranscriptomicsSummary;
use epigenomics_core::EpigenomicsSummary;

use crate::pathway::EnrichmentResult;

/// Severity level of a discovered insight.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum InsightLevel {
    Info,
    Warning,
    Critical,
}

impl InsightLevel {
    pub fn as_str(self) -> &'static str {
        match self {
            InsightLevel::Info => "INFO",
            InsightLevel::Warning => "WARN",
            InsightLevel::Critical => "CRIT",
        }
    }

    pub fn color_hex(self) -> &'static str {
        match self {
            InsightLevel::Info => "#28a745",
            InsightLevel::Warning => "#ffc107",
            InsightLevel::Critical => "#dc3545",
        }
    }
}

/// Which modality generated this insight.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum InsightModality {
    Genomics,
    Transcriptomics,
    Epigenomics,
    Integration,
}

impl InsightModality {
    pub fn as_str(self) -> &'static str {
        match self {
            InsightModality::Genomics => "Genomics",
            InsightModality::Transcriptomics => "Transcriptomics",
            InsightModality::Epigenomics => "Epigenomics",
            InsightModality::Integration => "Integration",
        }
    }
}

/// A plain-English biological finding produced by the insight engine.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Insight {
    pub level: InsightLevel,
    pub modality: InsightModality,
    pub message: String,
}

impl Insight {
    fn info(modality: InsightModality, message: impl Into<String>) -> Self {
        Self { level: InsightLevel::Info, modality, message: message.into() }
    }

    fn warn(modality: InsightModality, message: impl Into<String>) -> Self {
        Self { level: InsightLevel::Warning, modality, message: message.into() }
    }

    fn crit(modality: InsightModality, message: impl Into<String>) -> Self {
        Self { level: InsightLevel::Critical, modality, message: message.into() }
    }
}

/// Derive plain-English biological insights from the combined analysis results.
///
/// Rules are evaluated in priority order; each rule produces at most one insight.
pub fn derive_insights(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    corr: &[[f64; 3]; 3],
    enrichment: &[EnrichmentResult],
) -> Vec<Insight> {
    let mut insights = Vec::new();

    // --- Genomics rules ---
    check_titv_ratio(genomics, &mut insights);
    check_variant_load(genomics, &mut insights);
    check_high_impact_genes(genomics, &mut insights);
    check_indel_burden(genomics, &mut insights);

    // --- Transcriptomics rules ---
    check_expression_breadth(transcr, &mut insights);
    check_differential_expression(transcr, &mut insights);
    check_low_expression(transcr, &mut insights);

    // --- Epigenomics rules ---
    check_global_methylation(epigen, &mut insights);
    check_cpg_islands(epigen, &mut insights);
    check_extreme_methylation(epigen, &mut insights);

    // --- Integration rules ---
    check_cross_modality_correlation(corr, &mut insights);
    check_pathway_enrichment(enrichment, &mut insights);
    check_silencing_signature(genomics, transcr, epigen, &mut insights);

    insights
}

fn check_titv_ratio(g: &GenomicsSummary, out: &mut Vec<Insight>) {
    if g.titv_ratio > 0.0 && g.titv_ratio < 1.5 {
        out.push(Insight::warn(
            InsightModality::Genomics,
            format!(
                "Ti/Tv ratio of {:.2} is below the expected 2.0–2.1 for whole-genome sequencing. \
                 Possible sequencing artifact, FFPE damage, or enrichment bias.",
                g.titv_ratio
            ),
        ));
    } else if g.titv_ratio >= 2.0 && g.titv_ratio <= 2.2 {
        out.push(Insight::info(
            InsightModality::Genomics,
            format!(
                "Ti/Tv ratio of {:.2} is within the expected range for whole-genome sequencing.",
                g.titv_ratio
            ),
        ));
    } else if g.titv_ratio > 2.5 {
        out.push(Insight::warn(
            InsightModality::Genomics,
            format!(
                "Ti/Tv ratio of {:.2} is elevated. This may indicate hypermutation or a \
                 CpG-driven mutational signature.",
                g.titv_ratio
            ),
        ));
    }
}

fn check_variant_load(g: &GenomicsSummary, out: &mut Vec<Insight>) {
    if g.total_variants == 0 {
        out.push(Insight::warn(
            InsightModality::Genomics,
            "No variants detected. Check that the VCF file is correctly formatted and non-empty."
                .to_string(),
        ));
    } else if g.total_variants > 5_000_000 {
        out.push(Insight::warn(
            InsightModality::Genomics,
            format!(
                "High variant load: {} variants detected. Typical WGS yields 3–4 M SNPs.",
                g.total_variants
            ),
        ));
    } else {
        out.push(Insight::info(
            InsightModality::Genomics,
            format!(
                "{} variants detected ({} SNPs, {} indels, Ti/Tv = {:.2}).",
                g.total_variants, g.snp_count, g.indel_count, g.titv_ratio
            ),
        ));
    }
}

fn check_high_impact_genes(g: &GenomicsSummary, out: &mut Vec<Insight>) {
    let n = g.high_impact.len();
    if n == 0 {
        return;
    }

    // Check for known tumor suppressor / oncogene hits
    let flagged: Vec<&str> = ["TP53", "KRAS", "BRCA1", "BRCA2", "APC", "PTEN", "RB1", "VHL", "MLH1", "MSH2"]
        .iter()
        .filter(|&&gene| g.high_impact_genes.iter().any(|g| g == gene))
        .copied()
        .collect();

    if !flagged.is_empty() {
        out.push(Insight::crit(
            InsightModality::Genomics,
            format!(
                "High-impact variants detected in cancer driver genes: {}. \
                 These may be pathogenic. Manual review recommended.",
                flagged.join(", ")
            ),
        ));
    } else {
        out.push(Insight::info(
            InsightModality::Genomics,
            format!("{} high-quality variants (QUAL > 30) identified.", n),
        ));
    }
}

fn check_indel_burden(g: &GenomicsSummary, out: &mut Vec<Insight>) {
    if g.snp_count > 0 {
        let indel_frac = g.indel_count as f64 / (g.snp_count + g.indel_count) as f64;
        if indel_frac > 0.3 {
            out.push(Insight::warn(
                InsightModality::Genomics,
                format!(
                    "High indel fraction ({:.1}%). This may indicate microsatellite instability \
                     (MSI) or BRCA-associated homologous recombination deficiency.",
                    indel_frac * 100.0
                ),
            ));
        }
    }
}

fn check_expression_breadth(t: &TranscriptomicsSummary, out: &mut Vec<Insight>) {
    if t.total_genes == 0 {
        out.push(Insight::warn(
            InsightModality::Transcriptomics,
            "No genes found in expression matrix. Check TSV format.".to_string(),
        ));
        return;
    }
    let pct = t.expressed_genes as f64 / t.total_genes as f64 * 100.0;
    if pct < 30.0 {
        out.push(Insight::warn(
            InsightModality::Transcriptomics,
            format!(
                "Only {:.1}% of genes ({}/{}) have TPM ≥ 1. Sample may have low RNA quality \
                 or represent a highly specialized cell type.",
                pct, t.expressed_genes, t.total_genes
            ),
        ));
    } else {
        out.push(Insight::info(
            InsightModality::Transcriptomics,
            format!(
                "{} of {} genes ({:.1}%) expressed at TPM ≥ 1.",
                t.expressed_genes, t.total_genes, pct
            ),
        ));
    }
}

fn check_differential_expression(t: &TranscriptomicsSummary, out: &mut Vec<Insight>) {
    if let Some(ref de) = t.diff_expr {
        let sig = de.iter().filter(|r| r.log2_fold_change.abs() >= 2.0).count();
        if sig > 500 {
            out.push(Insight::crit(
                InsightModality::Transcriptomics,
                format!(
                    "{} genes with |log₂FC| ≥ 2 between samples. Large-scale transcriptional \
                     reprogramming detected.",
                    sig
                ),
            ));
        } else if sig > 0 {
            out.push(Insight::info(
                InsightModality::Transcriptomics,
                format!(
                    "{} significantly differentially expressed genes (|log₂FC| ≥ 2) between the \
                     two samples.",
                    sig
                ),
            ));
        }
    }
}

fn check_low_expression(t: &TranscriptomicsSummary, out: &mut Vec<Insight>) {
    let n = t.low_expression_genes.len();
    if n > 10_000 {
        out.push(Insight::info(
            InsightModality::Transcriptomics,
            format!(
                "{} genes below TPM < 1 threshold (low/absent expression).",
                n
            ),
        ));
    }
}

fn check_global_methylation(e: &EpigenomicsSummary, out: &mut Vec<Insight>) {
    let pct = e.global_methylation_pct;
    if pct < 40.0 && e.total_sites > 0 {
        out.push(Insight::crit(
            InsightModality::Epigenomics,
            format!(
                "Global methylation of {:.1}% is severely hypomethylated. Values < 40% are \
                 associated with chromosomal instability and loss of imprinting.",
                pct
            ),
        ));
    } else if pct > 90.0 && e.total_sites > 0 {
        out.push(Insight::warn(
            InsightModality::Epigenomics,
            format!(
                "Global methylation of {:.1}% is unusually high. Possible technical artifact \
                 or hypermethylation phenotype (CIMP).",
                pct
            ),
        ));
    } else if e.total_sites > 0 {
        out.push(Insight::info(
            InsightModality::Epigenomics,
            format!(
                "Global methylation: {:.1}% across {} CpG sites.",
                pct, e.total_sites
            ),
        ));
    }
}

fn check_cpg_islands(e: &EpigenomicsSummary, out: &mut Vec<Insight>) {
    let n = e.cpg_islands.len();
    if n > 0 {
        let methylated = e
            .cpg_islands
            .iter()
            .filter(|i| i.mean_methylation > 70.0)
            .count();
        if methylated > n / 2 {
            out.push(Insight::warn(
                InsightModality::Epigenomics,
                format!(
                    "{} of {} CpG islands ({:.0}%) are hypermethylated (mean > 70%). \
                     CpG island methylator phenotype (CIMP) suspected.",
                    methylated,
                    n,
                    methylated as f64 / n as f64 * 100.0
                ),
            ));
        } else {
            out.push(Insight::info(
                InsightModality::Epigenomics,
                format!("{} CpG islands detected.", n),
            ));
        }
    }
}

fn check_extreme_methylation(e: &EpigenomicsSummary, out: &mut Vec<Insight>) {
    let hyper = e.hypermethylated.len();
    let hypo = e.hypomethylated.len();
    if hyper > 0 || hypo > 0 {
        out.push(Insight::info(
            InsightModality::Epigenomics,
            format!(
                "{} hypermethylated regions (>80%) and {} hypomethylated regions (<20%) detected.",
                hyper, hypo
            ),
        ));
    }
}

fn check_cross_modality_correlation(corr: &[[f64; 3]; 3], out: &mut Vec<Insight>) {
    let geno_epi = corr[0][2].abs();
    let geno_trans = corr[0][1].abs();
    let trans_epi = corr[1][2].abs();

    if geno_epi > 0.7 {
        out.push(Insight::info(
            InsightModality::Integration,
            format!(
                "Strong genomics–epigenomics correlation (r = {:.2}) suggests that variant \
                 burden is co-localised with methylation changes — possible variant-driven \
                 epigenetic remodelling.",
                corr[0][2]
            ),
        ));
    }

    if geno_trans > 0.7 {
        out.push(Insight::info(
            InsightModality::Integration,
            format!(
                "Strong genomics–transcriptomics correlation (r = {:.2}) indicates that \
                 chromosomal variant density is associated with altered gene expression.",
                corr[0][1]
            ),
        ));
    }

    if trans_epi > 0.7 {
        out.push(Insight::warn(
            InsightModality::Integration,
            format!(
                "Strong transcriptomics–epigenomics correlation (r = {:.2}) is consistent \
                 with epigenetic silencing of transcription.",
                corr[1][2]
            ),
        ));
    }
}

fn check_pathway_enrichment(enrichment: &[EnrichmentResult], out: &mut Vec<Insight>) {
    if let Some(top) = enrichment.first() {
        out.push(Insight::info(
            InsightModality::Integration,
            format!(
                "Top enriched pathway: '{}' (score = {:.3}, {} overlapping genes). \
                 This pathway may be functionally altered across modalities.",
                top.pathway_name, top.score, top.overlap
            ),
        ));
    }
}

fn check_silencing_signature(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    out: &mut Vec<Insight>,
) {
    // Look for genes that appear mutated (high-impact) AND are low-expressed
    let mutated_genes: std::collections::HashSet<&str> = genomics
        .high_impact_genes
        .iter()
        .map(|g| g.as_str())
        .collect();

    let low_expr_genes: std::collections::HashSet<&str> = transcr
        .low_expression_genes
        .iter()
        .map(|g| g.as_str())
        .collect();

    let silenced: Vec<&str> = mutated_genes
        .intersection(&low_expr_genes)
        .copied()
        .collect();

    if !silenced.is_empty() {
        let examples: Vec<&str> = silenced.iter().take(3).copied().collect();
        let has_hyper = !epigen.hypermethylated.is_empty();
        if has_hyper {
            out.push(Insight::crit(
                InsightModality::Integration,
                format!(
                    "Multi-omic silencing signature detected: gene(s) {} are mutated \
                     (genomics), underexpressed (transcriptomics), AND the sample shows \
                     widespread hypermethylation (epigenomics). This pattern is consistent \
                     with epigenetic silencing of tumour suppressors.",
                    examples.join(", ")
                ),
            ));
        } else {
            out.push(Insight::warn(
                InsightModality::Integration,
                format!(
                    "Gene(s) {} are mutated (genomics) AND underexpressed (transcriptomics). \
                     Consider epigenetic analysis to determine whether silencing is involved.",
                    examples.join(", ")
                ),
            ));
        }
    }
}
