use serde::{Deserialize, Serialize};
use transcriptomics_core::TranscriptomicsSummary;

/// Checkpoint gene panel for immune evasion scoring.
/// Genes, canonical_aliases, weight
/// Reference: Chen & Mellman 2017 (Nature), Ribas & Wolchok 2018 (Science)
const CHECKPOINT_GENES: &[(&str, &[&str], f64)] = &[
    ("CD274",  &["PD-L1", "PDCD1LG1"],        1.0),  // PD-L1
    ("PDCD1",  &["PD-1"],                      0.9),  // PD-1 receptor
    ("CTLA4",  &["CTLA-4"],                    1.0),  // CTLA-4
    ("LAG3",   &["LAG-3"],                     0.8),  // LAG-3
    ("HAVCR2", &["TIM3", "TIM-3"],             0.8),  // TIM-3
    ("TIGIT",  &[],                             0.8),  // TIGIT
    ("FOXP3",  &[],                             0.6),  // Treg marker
];

const B2M: (&str, &[&str]) = ("B2M", &[]);

/// Immune evasion score derived from checkpoint gene expression.
///
/// Scoring formula (Chen & Mellman 2017):
/// - checkpoint_score = weighted mean of log2(TPM+1)/log2(101) for inhibitory checkpoint genes
/// - antigen_presentation_score = log2(B2M_TPM+1)/log2(101) (high = good, low = escape)
/// - immune_evasion_score = 0.6 × checkpoint_score + 0.4 × (1 - antigen_presentation_score)
///
/// Classification: HIGH > 0.40, MODERATE 0.20–0.40, LOW < 0.20
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImmuneEvasionScore {
    /// Weighted mean normalized checkpoint gene expression (0–1).
    pub checkpoint_score: f64,
    /// Antigen presentation score from B2M expression (0–1, high = intact).
    pub antigen_presentation_score: f64,
    /// Final composite immune evasion score (0–1).
    pub immune_evasion_score: f64,
    /// "HIGH" (>0.40), "MODERATE" (0.20–0.40), "LOW" (<0.20).
    pub evasion_class: String,
    /// Detected checkpoint gene expressions: (gene_symbol, mean_tpm).
    pub detected_genes: Vec<(String, f64)>,
    /// Checkpoint genes not found in expression data.
    pub missing_genes: Vec<String>,
    /// B2M mean TPM (None if not detected).
    pub b2m_tpm: Option<f64>,
    /// Note when detected_genes < 3 (insufficient data for reliable scoring).
    pub note: Option<String>,
}

/// Compute immune evasion score from transcriptomics gene stats.
pub fn compute_immune_evasion(transcr: &TranscriptomicsSummary) -> Option<ImmuneEvasionScore> {
    if transcr.gene_stats.is_empty() {
        return None;
    }

    let log2_101 = (101.0_f64).log2();

    let mut weighted_sum = 0.0_f64;
    let mut weight_total = 0.0_f64;
    let mut detected_genes: Vec<(String, f64)> = Vec::new();
    let mut missing_genes: Vec<String> = Vec::new();

    for (canonical, aliases, weight) in CHECKPOINT_GENES {
        // Try canonical name, then aliases
        let tpm_opt = transcr.gene_stats.get(*canonical).map(|s| s.mean).or_else(|| {
            aliases.iter().find_map(|alias| transcr.gene_stats.get(*alias).map(|s| s.mean))
        });

        match tpm_opt {
            Some(tpm) => {
                let norm = (tpm + 1.0).log2() / log2_101;
                weighted_sum += norm * weight;
                weight_total += weight;
                detected_genes.push((canonical.to_string(), tpm));
            }
            None => {
                missing_genes.push(canonical.to_string());
            }
        }
    }

    let checkpoint_score = if weight_total > 0.0 {
        (weighted_sum / weight_total).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // B2M (antigen presentation)
    let b2m_tpm = transcr.gene_stats.get(B2M.0).map(|s| s.mean).or_else(|| {
        B2M.1.iter().find_map(|alias| transcr.gene_stats.get(*alias).map(|s| s.mean))
    });
    let antigen_presentation_score = match b2m_tpm {
        Some(tpm) => ((tpm + 1.0).log2() / log2_101).clamp(0.0, 1.0),
        None => 0.5, // unknown: assume average
    };

    let immune_evasion_score = (0.6 * checkpoint_score + 0.4 * (1.0 - antigen_presentation_score)).clamp(0.0, 1.0);

    let evasion_class = if immune_evasion_score > 0.40 {
        "HIGH".to_string()
    } else if immune_evasion_score >= 0.20 {
        "MODERATE".to_string()
    } else {
        "LOW".to_string()
    };

    let note = if detected_genes.len() < 3 {
        Some(format!(
            "Only {}/{} checkpoint genes detected — immune evasion score may be unreliable",
            detected_genes.len(),
            CHECKPOINT_GENES.len()
        ))
    } else {
        None
    };

    Some(ImmuneEvasionScore {
        checkpoint_score,
        antigen_presentation_score,
        immune_evasion_score,
        evasion_class,
        detected_genes,
        missing_genes,
        b2m_tpm,
        note,
    })
}
