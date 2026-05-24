pub mod correlation;
pub mod gmt;
pub mod gsea;
pub mod insights;
pub mod mofa;
pub mod pathway;
pub mod pca;

pub use gmt::{gmt_enrichment_analysis, parse_gmt, GmtPathway};
pub use gsea::{gsea_preranked, GseaResult};
pub use insights::{derive_insights, Insight, InsightLevel, InsightModality};
pub use mofa::{run_mofa, MofaConfig, MofaResult};
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

    let insights = derive_insights(genomics, transcr, epigen, &corr_arr, &top_pathways);

    Ok(IntegrationSummary {
        correlation_matrix,
        pca,
        mofa: mofa_result,
        top_pathways,
        insights,
    })
}
