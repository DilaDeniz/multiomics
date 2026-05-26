pub mod correlation;
pub mod gene_state;
pub mod gmt;
pub mod gsea;
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

use correlation::{build_cross_modality_matrix, pearson_correlation_matrix};

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
        }
    }
}

/// Run the full integration pipeline over all three modality summaries.
///
/// When `skip_ml` is true, the correlation matrix and PCA are skipped and
/// only insights from simple rule checks are computed.
pub fn run_integration(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
    skip_ml: bool,
) -> Result<IntegrationSummary> {
    // Collect query genes from DE results and high-impact genomic variants
    let mut query_genes: Vec<String> = genomics.high_impact_genes.clone();
    if let Some(ref de) = transcr.diff_expr {
        let de_genes = significant_de_genes(de, 1.0, 1.0);
        query_genes.extend(de_genes);
    }
    query_genes.sort_unstable();
    query_genes.dedup();

    // Pathway enrichment
    let top_pathways = enrichment_analysis(&query_genes, 1);
    let top_pathways: Vec<EnrichmentResult> = top_pathways.into_iter().take(10).collect();

    // Correlation + PCA + MOFA+
    let (correlation_matrix, pca, mofa_result) = if skip_ml {
        let identity = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let pca_fallback = PcaResult {
            points: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            explained_variance_ratio: vec![1.0, 0.0],
        };
        (identity, pca_fallback, None)
    } else {
        let matrix = build_cross_modality_matrix(genomics, transcr, epigen);
        let corr = pearson_correlation_matrix(&matrix)?;
        let corr_vec: Vec<Vec<f64>> = (0..corr.nrows())
            .map(|i| corr.row(i).iter().copied().collect())
            .collect();

        let pca_result = run_pca(&matrix, 2)?;

        // MOFA+: one row per modality, features = the same cross-modality columns.
        // With 3 samples (rows) and K=2 factors we extract the maximal joint structure.
        let mofa = if matrix.nrows() >= 2 && matrix.ncols() >= 2 {
            let view_names = ["genomics", "transcriptomics", "epigenomics"];
            let views: Vec<(&str, &ndarray::Array2<f64>)> = view_names
                .iter()
                .zip(std::iter::repeat(&matrix))
                .map(|(&n, m)| (n, m))
                .collect();
            let cfg = MofaConfig {
                n_factors: 2,
                max_iter: 500,
                tol: 1e-5,
                ..Default::default()
            };
            match run_mofa(&views, &cfg) {
                Ok(r) => {
                    log::info!(
                        "MOFA+ converged in {} iterations (ELBO={:.2})",
                        r.n_iter,
                        r.elbo
                    );
                    Some(r)
                }
                Err(e) => {
                    log::warn!("MOFA+ failed: {e:#}");
                    None
                }
            }
        } else {
            None
        };

        (corr_vec, pca_result, mofa)
    };

    // Convert correlation matrix to fixed-size array for insight engine
    let corr_arr: [[f64; 3]; 3] = [
        [
            correlation_matrix
                .first()
                .and_then(|r| r.first())
                .copied()
                .unwrap_or(1.0),
            correlation_matrix
                .first()
                .and_then(|r| r.get(1))
                .copied()
                .unwrap_or(0.0),
            correlation_matrix
                .first()
                .and_then(|r| r.get(2))
                .copied()
                .unwrap_or(0.0),
        ],
        [
            correlation_matrix
                .get(1)
                .and_then(|r| r.first())
                .copied()
                .unwrap_or(0.0),
            correlation_matrix
                .get(1)
                .and_then(|r| r.get(1))
                .copied()
                .unwrap_or(1.0),
            correlation_matrix
                .get(1)
                .and_then(|r| r.get(2))
                .copied()
                .unwrap_or(0.0),
        ],
        [
            correlation_matrix
                .get(2)
                .and_then(|r| r.first())
                .copied()
                .unwrap_or(0.0),
            correlation_matrix
                .get(2)
                .and_then(|r| r.get(1))
                .copied()
                .unwrap_or(0.0),
            correlation_matrix
                .get(2)
                .and_then(|r| r.get(2))
                .copied()
                .unwrap_or(1.0),
        ],
    ];

    let mut insights = derive_insights(genomics, transcr, epigen, &corr_arr, &top_pathways);

    // Run paradox detection across all three modalities
    let paradoxes = detect_paradoxes(genomics, transcr, epigen);
    log::info!("Multi-modal paradox detection: {} paradoxes found", paradoxes.len());
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
    let count_active = gene_states.iter().filter(|g| g.state == GeneState::Active).count();
    let count_silenced = gene_states.iter().filter(|g| g.state == GeneState::Silenced).count();
    let count_paradoxical = gene_states.iter().filter(|g| g.state == GeneState::Paradoxical).count();
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
    })
}
