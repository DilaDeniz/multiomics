use std::collections::{HashMap, HashSet};
use serde::{Deserialize, Serialize};

use epigenomics_core::EpigenomicsSummary;
use genomics_core::GenomicsSummary;
use transcriptomics_core::{DiffExprResult, TranscriptomicsSummary};

use crate::paradox::{GeneParadox, ParadoxKind};

/// Regulatory state of a gene derived from its multi-modal molecular profile.
///
/// Based on Roadmap Epigenomics chromatin state logic (Roadmap Epigenomics
/// Consortium, Nature 2015) simplified to three measurable modalities.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum GeneState {
    /// Expressed (TPM ≥ 10) with low or absent methylation (< 30%).
    /// Consistent with an active promoter and open chromatin.
    Active,
    /// Not expressed (TPM < 1) with high methylation (> 70%) or a silencing variant.
    /// Consistent with a repressed or imprinted locus.
    Silenced,
    /// Low expression (TPM 1–10) with low methylation (< 30%).
    /// Gene is accessible but not strongly transcribed; may be primed for activation.
    Poised,
    /// Low expression (TPM 1–10) with high methylation (> 70%).
    /// Conflicting marks reminiscent of bivalent stem-cell promoters.
    Bivalent,
    /// Carries a high-impact variant that co-occurs with significant differential
    /// expression (|log2FC| > 1, padj < 0.05). Variant likely drives expression change.
    VariantDriven,
    /// Multi-modal contradiction detected (e.g., methylated but expressed).
    /// See the paradox report for details.
    Paradoxical,
    /// Gene observed in at least one modality but state could not be determined.
    Unknown,
}

impl GeneState {
    pub fn as_str(&self) -> &'static str {
        match self {
            GeneState::Active => "ACTIVE",
            GeneState::Silenced => "SILENCED",
            GeneState::Poised => "POISED",
            GeneState::Bivalent => "BIVALENT",
            GeneState::VariantDriven => "VARIANT_DRIVEN",
            GeneState::Paradoxical => "PARADOXICAL",
            GeneState::Unknown => "UNKNOWN",
        }
    }

    /// CSS color for the HTML report.
    pub fn html_color(&self) -> &'static str {
        match self {
            GeneState::Active => "#27ae60",        // green
            GeneState::Silenced => "#7f8c8d",      // grey
            GeneState::Poised => "#f39c12",        // amber
            GeneState::Bivalent => "#8e44ad",      // purple
            GeneState::VariantDriven => "#e74c3c", // red
            GeneState::Paradoxical => "#e67e22",   // orange
            GeneState::Unknown => "#bdc3c7",       // light grey
        }
    }

    /// Sort priority: lower = higher priority (shown first in table).
    fn priority(&self) -> u8 {
        match self {
            GeneState::VariantDriven => 0,
            GeneState::Paradoxical => 1,
            GeneState::Silenced => 2,
            GeneState::Active => 3,
            GeneState::Bivalent => 4,
            GeneState::Poised => 5,
            GeneState::Unknown => 6,
        }
    }
}

/// The full regulatory profile of a single gene across all modalities.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneRegulatoryProfile {
    pub gene: String,
    pub state: GeneState,
    /// Mean TPM across samples (None if not observed in transcriptomics).
    pub mean_tpm: Option<f64>,
    /// Mean methylation at annotated CpG sites (None if not in BED annotation).
    pub mean_methylation: Option<f64>,
    /// Whether a high-impact variant (QUAL > 30) is present.
    pub has_variant: bool,
    /// Max QUAL among variants for this gene.
    pub max_variant_qual: Option<f32>,
    /// log2 fold-change if DE results are available.
    pub log2_fold_change: Option<f64>,
    /// BH-adjusted p-value for DE.
    pub padj: Option<f64>,
    /// Human-readable description of the state assignment.
    pub description: String,
}

/// Classify every gene that appears in at least one modality into a regulatory state.
///
/// Returns profiles sorted by state priority (VariantDriven, Paradoxical first),
/// then by mean_tpm descending within each state.
pub fn classify_gene_states(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    paradoxes: &[GeneParadox],
) -> Vec<GeneRegulatoryProfile> {
    // 1. Build set of paradoxical genes (MethylatedButExpressed only — clearest regulatory paradox)
    let paradox_genes: HashSet<&str> = paradoxes
        .iter()
        .filter(|p| p.kind == ParadoxKind::MethylatedButExpressed)
        .map(|p| p.gene.as_str())
        .collect();

    // 2. Build variant gene set and gene → max QUAL map from high_impact Vec
    let mut gene_max_qual: HashMap<&str, f32> = HashMap::new();
    for v in &genomics.high_impact {
        if let Some(ref g) = v.gene {
            let entry = gene_max_qual.entry(g.as_str()).or_insert(0.0_f32);
            if v.qual > *entry {
                *entry = v.qual;
            }
        }
    }
    // Also index high_impact_genes as a HashSet (may include genes without QUAL record)
    let variant_gene_set: HashSet<&str> = genomics
        .high_impact_genes
        .iter()
        .map(|s| s.as_str())
        .collect();

    // 3. Build HashMap<gene, &DiffExprResult> from transcr.diff_expr
    let de_map: HashMap<&str, &DiffExprResult> =
        if let Some(ref de) = transcr.diff_expr {
            de.iter().map(|r| (r.gene_id.as_str(), r)).collect()
        } else {
            HashMap::new()
        };

    // 4. Collect all genes seen in any modality
    let mut all_genes: HashSet<String> = HashSet::new();
    for g in transcr.gene_stats.keys() {
        all_genes.insert(g.clone());
    }
    for g in &genomics.high_impact_genes {
        all_genes.insert(g.clone());
    }
    for g in epigen.gene_methylation.keys() {
        all_genes.insert(g.clone());
    }

    // 5. Classify each gene
    let mut profiles: Vec<GeneRegulatoryProfile> = all_genes
        .into_iter()
        .map(|gene| {
            let mean_tpm = transcr.gene_stats.get(&gene).map(|s| s.mean);
            let mean_meth = epigen.gene_methylation.get(&gene).copied();
            let has_variant = variant_gene_set.contains(gene.as_str());
            let max_variant_qual = gene_max_qual.get(gene.as_str()).copied();
            let de = de_map.get(gene.as_str()).copied();
            let log2_fold_change = de.map(|d| d.log2_fold_change);
            let padj = de.map(|d| d.padj);

            let state = classify_single_gene(
                gene.as_str(),
                mean_tpm,
                mean_meth,
                has_variant,
                de,
                &paradox_genes,
            );

            let description = build_description(&state, mean_tpm, mean_meth, max_variant_qual, log2_fold_change, padj);

            GeneRegulatoryProfile {
                gene,
                state,
                mean_tpm,
                mean_methylation: mean_meth,
                has_variant,
                max_variant_qual,
                log2_fold_change,
                padj,
                description,
            }
        })
        .collect();

    // 7. Sort: by state priority ascending, then mean_tpm descending within each group
    profiles.sort_by(|a, b| {
        let pa = a.state.priority();
        let pb = b.state.priority();
        pa.cmp(&pb).then_with(|| {
            let ta = a.mean_tpm.unwrap_or(0.0);
            let tb = b.mean_tpm.unwrap_or(0.0);
            tb.partial_cmp(&ta)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });

    // 8. Cap at 500 genes
    profiles.truncate(500);

    profiles
}

/// Classify a single gene into a GeneState using the priority logic.
fn classify_single_gene(
    gene: &str,
    mean_tpm: Option<f64>,
    mean_meth: Option<f64>,
    has_variant: bool,
    de: Option<&DiffExprResult>,
    paradox_genes: &HashSet<&str>,
) -> GeneState {
    // VARIANT_DRIVEN: variant + significant DE (|log2FC| > 1, padj < 0.05)
    if has_variant
        && de
            .map(|d| {
                !d.padj.is_nan()
                    && d.padj < 0.05
                    && d.log2_fold_change.abs() > 1.0
            })
            .unwrap_or(false)
    {
        return GeneState::VariantDriven;
    }

    // PARADOXICAL: MethylatedButExpressed flagged
    if paradox_genes.contains(gene) {
        return GeneState::Paradoxical;
    }

    // SILENCED: TPM < 1 AND (high methylation > 70% OR has variant)
    if mean_tpm.map(|t| t < 1.0).unwrap_or(false)
        && (mean_meth.map(|m| m > 70.0).unwrap_or(false) || has_variant)
    {
        return GeneState::Silenced;
    }

    // ACTIVE: TPM >= 10 AND low or absent methylation (absent = assume active)
    if mean_tpm.map(|t| t >= 10.0).unwrap_or(false)
        && mean_meth.map(|m| m < 30.0).unwrap_or(true)
    {
        return GeneState::Active;
    }

    // BIVALENT: TPM 1–10 AND high methylation > 70%
    if mean_tpm.map(|t| (1.0..10.0).contains(&t)).unwrap_or(false)
        && mean_meth.map(|m| m > 70.0).unwrap_or(false)
    {
        return GeneState::Bivalent;
    }

    // POISED: TPM 1–10 AND low methylation < 30%
    if mean_tpm.map(|t| (1.0..10.0).contains(&t)).unwrap_or(false)
        && mean_meth.map(|m| m < 30.0).unwrap_or(false)
    {
        return GeneState::Poised;
    }

    // ACTIVE fallback: expressed regardless of methylation data
    if mean_tpm.map(|t| t >= 10.0).unwrap_or(false) {
        return GeneState::Active;
    }

    GeneState::Unknown
}

/// Build a human-readable description for a gene profile.
fn build_description(
    state: &GeneState,
    mean_tpm: Option<f64>,
    mean_meth: Option<f64>,
    max_variant_qual: Option<f32>,
    log2_fold_change: Option<f64>,
    padj: Option<f64>,
) -> String {
    match state {
        GeneState::Active => {
            let tpm_str = mean_tpm.map(|t| format!("TPM={t:.1}")).unwrap_or_else(|| "TPM=N/A".to_string());
            let meth_str = mean_meth
                .map(|m| format!(" with low methylation ({m:.1}%)"))
                .unwrap_or_else(|| " (no methylation data)".to_string());
            format!("Expressed ({tpm_str}){meth_str}")
        }
        GeneState::Silenced => {
            let tpm_str = mean_tpm.map(|t| format!("TPM={t:.2}")).unwrap_or_else(|| "TPM=N/A".to_string());
            let meth_str = mean_meth
                .map(|m| format!(" with hypermethylated promoter ({m:.1}%)"))
                .unwrap_or_default();
            format!("Silent ({tpm_str}){meth_str}")
        }
        GeneState::Poised => {
            let tpm_str = mean_tpm.map(|t| format!("{t:.1}")).unwrap_or_else(|| "N/A".to_string());
            let meth_str = mean_meth.map(|m| format!("{m:.1}%")).unwrap_or_else(|| "N/A".to_string());
            format!("Low expression (TPM={tpm_str}) with low methylation ({meth_str}) — gene accessible but not actively transcribed")
        }
        GeneState::Bivalent => {
            let tpm_str = mean_tpm.map(|t| format!("{t:.1}")).unwrap_or_else(|| "N/A".to_string());
            let meth_str = mean_meth.map(|m| format!("{m:.1}%")).unwrap_or_else(|| "N/A".to_string());
            format!("Low expression (TPM={tpm_str}) with high methylation ({meth_str}) — conflicting bivalent marks")
        }
        GeneState::VariantDriven => {
            let qual_str = max_variant_qual
                .map(|q| format!("{q:.0}"))
                .unwrap_or_else(|| ">30".to_string());
            let lfc_str = log2_fold_change
                .map(|l| format!("{l:.2}"))
                .unwrap_or_else(|| "N/A".to_string());
            let padj_str = padj
                .map(|p| format!("{p:.2e}"))
                .unwrap_or_else(|| "N/A".to_string());
            format!("Variant (QUAL={qual_str}) co-occurs with log2FC={lfc_str} (padj={padj_str})")
        }
        GeneState::Paradoxical => {
            let meth_str = mean_meth.map(|m| format!("{m:.1}%")).unwrap_or_else(|| "N/A".to_string());
            let tpm_str = mean_tpm.map(|t| format!("{t:.1}")).unwrap_or_else(|| "N/A".to_string());
            format!("Methylated ({meth_str}) but actively expressed (TPM={tpm_str}) — see paradox report")
        }
        GeneState::Unknown => {
            "State could not be determined from available multi-modal data".to_string()
        }
    }
}
