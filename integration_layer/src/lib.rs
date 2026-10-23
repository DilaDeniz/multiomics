pub mod correlation;
pub mod gene_state;
pub mod gmt;
pub mod gsea;
pub mod immune;
pub mod insights;
pub mod mofa;
pub mod paradox;
pub mod pathway;
pub mod pca;

pub use gene_state::{classify_gene_states, GeneRegulatoryProfile, GeneState};
pub use gmt::{gmt_enrichment_analysis, parse_gmt, GmtPathway};
pub use gsea::{gsea_preranked, GseaResult};
pub use insights::{derive_insights, Insight, InsightLevel, InsightModality};
pub use mofa::{run_mofa, MofaConfig, MofaResult};
pub use paradox::{detect_paradoxes, GeneParadox, ParadoxEvidence, ParadoxKind};
pub use pathway::{enrichment_analysis, EnrichmentResult, KeggPathway, KEGG_PATHWAYS};
pub use pca::{run_pca, PcaResult};

use anyhow::Result;
use serde::{Deserialize, Serialize};

use epigenomics_core::EpigenomicsSummary;
use genomics_core::GenomicsSummary;
use transcriptomics_core::{significant_de_genes, TranscriptomicsSummary};

/// All cross-modality analysis results.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntegrationSummary {
    /// 3×3 Pearson correlation matrix: [genomics, transcriptomics, epigenomics].
    pub correlation_matrix: Vec<Vec<f64>>,
    /// PCA projection of the three modalities to 2D.
    pub pca: PcaResult,
    /// MOFA+ joint factor analysis across all modalities (None when --no-ml).
    pub mofa: Option<MofaResult>,
    /// Top pathway enrichment results.
    pub top_pathways: Vec<EnrichmentResult>,
    /// Plain-English biological insights.
    pub insights: Vec<Insight>,
    /// Multi-modal biological paradoxes detected across all three modalities.
    pub paradoxes: Vec<GeneParadox>,
    /// Per-gene regulatory state classifications (Active/Silenced/Poised/Bivalent/VariantDriven/Paradoxical).
    pub gene_states: Vec<GeneRegulatoryProfile>,
    /// Cross-modal tumor purity estimate (VAF + methylation).
    #[serde(default)]
    pub tumor_purity: Option<genomics_core::cancer::TumorPurityResult>,
    /// Horvath epigenetic age clock result (cloned from epigenomics summary).
    #[serde(default)]
    pub methylation_age: Option<epigenomics_core::clock::MethylationAgeResult>,
    /// Tumor mutational burden result (cloned from genomics summary).
    #[serde(default)]
    pub tmb: Option<genomics_core::cancer::TmbResult>,
    /// Microsatellite instability result (cloned from genomics summary).
    #[serde(default)]
    pub msi: Option<genomics_core::cancer::MsiResult>,
    /// Immune evasion score derived from checkpoint gene expression.
    #[serde(default)]
    pub immune_evasion: Option<immune::ImmuneEvasionScore>,
}

impl IntegrationSummary {
    /// Produce an empty summary when `--no-ml` is requested.
    pub fn empty() -> Self {
        Self {
            correlation_matrix: vec![
                vec![1.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
                vec![0.0, 0.0, 1.0],
            ],
            pca: PcaResult {
                points: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                explained_variance_ratio: vec![1.0, 0.0],
            },
            mofa: None,
            top_pathways: Vec::new(),
            insights: Vec::new(),
            paradoxes: Vec::new(),
            gene_states: Vec::new(),
            tumor_purity: None,
            methylation_age: None,
            tmb: None,
            msi: None,
            immune_evasion: None,
        }
    }
}

/// Which modality summaries come from real input files. Summaries of absent
/// modalities are `Default` values and must not generate findings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Modalities {
    pub genomics: bool,
    pub transcriptomics: bool,
    pub epigenomics: bool,
}

impl Modalities {
    /// All three modalities were provided.
    pub const ALL: Self = Self {
        genomics: true,
        transcriptomics: true,
        epigenomics: true,
    };
}

/// Pearson correlation between per-chromosome variant counts and mean
/// methylation, over chromosomes present in both inputs (at least three).
/// Returns (r, number of chromosomes).
pub fn genomics_epigenomics_correlation(
    genomics: &GenomicsSummary,
    epigen: &EpigenomicsSummary,
) -> Option<(f64, usize)> {
    let mut shared: Vec<&String> = genomics
        .per_chrom
        .keys()
        .filter(|c| epigen.per_chrom.contains_key(*c))
        .collect();
    if shared.len() < 3 {
        return None;
    }
    shared.sort();
    let g = ndarray::Array1::from_iter(shared.iter().map(|c| genomics.per_chrom[*c].total as f64));
    let e =
        ndarray::Array1::from_iter(shared.iter().map(|c| epigen.per_chrom[*c].mean_methylation));
    Some((correlation::pearson_r(&g, &e), shared.len()))
}

/// Run the integration layer assuming all three modalities were provided.
pub fn run_integration(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    skip_ml: bool,
) -> Result<IntegrationSummary> {
    run_integration_with(genomics, transcr, epigen, skip_ml, Modalities::ALL)
}

/// Run the integration layer over the modality summaries that were provided.
///
/// The cross-modality correlation matrix is ordered [genomics,
/// transcriptomics, epigenomics]. Only the genomics–epigenomics entry is
/// defined: both are summarised per chromosome. The expression matrix carries
/// no genomic coordinates, so its entries are NaN (`null` in JSON). PCA and
/// MOFA+ need several samples measured in every modality and are not run on
/// these single-sample summaries (`pca` is empty, `mofa` is `None`).
pub fn run_integration_with(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    skip_ml: bool,
    provided: Modalities,
) -> Result<IntegrationSummary> {
    // Collect query genes from DE results and high-impact genomic variants
    let mut query_genes: Vec<String> = if provided.genomics {
        genomics.high_impact_genes.clone()
    } else {
        Vec::new()
    };
    if let Some(ref de) = transcr.diff_expr {
        let de_genes = significant_de_genes(de, 1.0, 1.0);
        query_genes.extend(de_genes);
    }
    query_genes.sort_unstable();
    query_genes.dedup();

    // Pathway enrichment
    let top_pathways = enrichment_analysis(&query_genes, 1);
    let top_pathways: Vec<EnrichmentResult> = top_pathways.into_iter().take(10).collect();

    let diag = |present: bool| if present { 1.0 } else { f64::NAN };
    let r_ge = if !skip_ml && provided.genomics && provided.epigenomics {
        genomics_epigenomics_correlation(genomics, epigen).map_or(f64::NAN, |(r, _)| r)
    } else {
        f64::NAN
    };
    let correlation_matrix = vec![
        vec![diag(provided.genomics), f64::NAN, r_ge],
        vec![f64::NAN, diag(provided.transcriptomics), f64::NAN],
        vec![r_ge, f64::NAN, diag(provided.epigenomics)],
    ];
    let pca = PcaResult {
        points: Vec::new(),
        explained_variance_ratio: Vec::new(),
    };
    let mofa_result: Option<MofaResult> = None;
    let corr_arr: [[f64; 3]; 3] = [
        [
            correlation_matrix[0][0],
            correlation_matrix[0][1],
            correlation_matrix[0][2],
        ],
        [
            correlation_matrix[1][0],
            correlation_matrix[1][1],
            correlation_matrix[1][2],
        ],
        [
            correlation_matrix[2][0],
            correlation_matrix[2][1],
            correlation_matrix[2][2],
        ],
    ];

    let mut insights = derive_insights(
        genomics,
        transcr,
        epigen,
        &corr_arr,
        &top_pathways,
        provided,
    );

    // Run paradox detection across all three modalities
    let paradoxes = detect_paradoxes(genomics, transcr, epigen);
    log::info!(
        "Multi-modal paradox detection: {} paradoxes found",
        paradoxes.len()
    );
    if !paradoxes.is_empty() {
        insights.push(insights::Insight {
            level: insights::InsightLevel::Info,
            modality: insights::InsightModality::Integration,
            message: format!(
                "[INFO] {} multi-modal paradoxes detected — see report for details",
                paradoxes.len()
            ),
        });
    }

    // Run per-gene regulatory state classification
    let gene_states = classify_gene_states(genomics, transcr, epigen, &paradoxes);
    let count_active = gene_states
        .iter()
        .filter(|g| g.state == GeneState::Active)
        .count();
    let count_silenced = gene_states
        .iter()
        .filter(|g| g.state == GeneState::Silenced)
        .count();
    let count_paradoxical = gene_states
        .iter()
        .filter(|g| g.state == GeneState::Paradoxical)
        .count();
    log::info!(
        "Gene state classification: {} genes classified ({} active, {} silenced, {} paradoxical)",
        gene_states.len(),
        count_active,
        count_silenced,
        count_paradoxical
    );

    // Warn if many silenced genes carry high-impact variants (tumor suppressor candidates)
    let silenced_with_variants = gene_states
        .iter()
        .filter(|g| g.state == GeneState::Silenced && g.has_variant)
        .count();
    if silenced_with_variants >= 2 {
        insights.push(insights::Insight {
            level: insights::InsightLevel::Warning,
            modality: insights::InsightModality::Integration,
            message: format!(
                "[WARN] {} genes silenced with high-impact variants — possible tumor suppressor candidates",
                silenced_with_variants
            ),
        });
    }

    // Tumor purity cross-modal estimation
    let tumor_purity_result = genomics_core::cancer::estimate_tumor_purity(
        &genomics.high_impact,
        epigen.global_methylation_pct,
    );

    // Cancer-specific insights
    if let Some(consensus) = tumor_purity_result.consensus_purity {
        if consensus > 0.7 {
            insights.push(insights::Insight {
                level: insights::InsightLevel::Info,
                modality: insights::InsightModality::Integration,
                message: format!(
                    "[INFO] High tumor purity estimated: {:.0}% (VAF-based)",
                    consensus * 100.0
                ),
            });
        }
    }
    if tumor_purity_result.discordant {
        let vaf_pct = tumor_purity_result
            .vaf_purity
            .map(|v| v * 100.0)
            .unwrap_or(0.0);
        let meth_pct = tumor_purity_result
            .methylation_purity
            .map(|v| v * 100.0)
            .unwrap_or(0.0);
        insights.push(insights::Insight {
            level: insights::InsightLevel::Warning,
            modality: insights::InsightModality::Integration,
            message: format!(
                "[WARN] Purity estimates discordant between VAF ({:.0}%) and methylation ({:.0}%) — possible tumor heterogeneity",
                vaf_pct, meth_pct
            ),
        });
    }
    if !genomics.kataegis_loci.is_empty() {
        insights.push(insights::Insight {
            level: insights::InsightLevel::Warning,
            modality: insights::InsightModality::Integration,
            message: format!(
                "[WARN] {} kataegis loci detected — APOBEC/AID mutagenesis signature",
                genomics.kataegis_loci.len()
            ),
        });
    }
    if let Some(ref hrd) = genomics.hrd {
        if hrd.hrd_class == "HRD-HIGH" {
            insights.push(insights::Insight {
                level: insights::InsightLevel::Warning,
                modality: insights::InsightModality::Integration,
                message: "[WARN] HRD-HIGH indel signature — may indicate BRCA1/2 deficiency; consider PARP inhibitor sensitivity".to_string(),
            });
        }
    }

    // Epigenetic age clock insights
    if let Some(ref ma) = epigen.methylation_age {
        insights.push(insights::Insight {
            level: insights::InsightLevel::Info,
            modality: insights::InsightModality::Epigenomics,
            message: format!(
                "[INFO] Epigenetic age: {:.1} years (coverage: {}/{} clock CpGs)",
                ma.biological_age, ma.cpgs_found, ma.cpgs_total
            ),
        });
        if ma.age_accelerated == Some(true) {
            if let Some(delta) = ma.age_delta {
                insights.push(insights::Insight {
                    level: insights::InsightLevel::Warning,
                    modality: insights::InsightModality::Epigenomics,
                    message: format!(
                        "[WARN] Epigenetic age acceleration detected (+{:.1} years) — associated with cancer and disease risk",
                        delta
                    ),
                });
            }
        }
    }

    // Immune evasion score
    let immune_evasion = immune::compute_immune_evasion(transcr);
    if let Some(ref ie) = immune_evasion {
        if ie.evasion_class == "HIGH" {
            insights.push(insights::Insight {
                level: insights::InsightLevel::Warning,
                modality: insights::InsightModality::Integration,
                message: "[WARN] HIGH immune evasion score — active checkpoint inhibition detected; consider immune checkpoint therapy".to_string(),
            });
        }
    }

    // TMB / MSI insights
    let tmb_high = genomics
        .tmb
        .as_ref()
        .map(|t| t.tmb_class == "TMB-H")
        .unwrap_or(false);
    let msi_high = genomics
        .msi
        .as_ref()
        .map(|m| m.msi_class == "MSI-H")
        .unwrap_or(false);

    if tmb_high && msi_high {
        insights.push(insights::Insight {
            level: insights::InsightLevel::Warning,
            modality: insights::InsightModality::Integration,
            message: "[WARN] TMB-H + MSI-H — strong immunotherapy candidate".to_string(),
        });
    } else {
        if tmb_high {
            if let Some(ref tmb) = genomics.tmb {
                insights.push(insights::Insight {
                    level: insights::InsightLevel::Warning,
                    modality: insights::InsightModality::Integration,
                    message: format!(
                        "[WARN] TMB-H: {:.1} mut/Mb — may benefit from immune checkpoint therapy (pembrolizumab)",
                        tmb.tmb
                    ),
                });
            }
        }
        if msi_high {
            if let Some(ref msi) = genomics.msi {
                insights.push(insights::Insight {
                    level: insights::InsightLevel::Warning,
                    modality: insights::InsightModality::Integration,
                    message: format!(
                        "[WARN] MSI-H detected (score={:.2}) — mismatch repair deficient; pembrolizumab indicated",
                        msi.msi_score
                    ),
                });
            }
        }
    }

    Ok(IntegrationSummary {
        correlation_matrix,
        pca,
        mofa: mofa_result,
        top_pathways,
        insights,
        paradoxes,
        gene_states,
        tumor_purity: Some(tumor_purity_result),
        methylation_age: epigen.methylation_age.clone(),
        tmb: genomics.tmb.clone(),
        msi: genomics.msi.clone(),
        immune_evasion,
    })
}
