use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

use epigenomics_core::EpigenomicsSummary;
use genomics_core::GenomicsSummary;
use transcriptomics_core::TranscriptomicsSummary;

/// Classification of why a gene is flagged as paradoxical or multi-hit.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ParadoxKind {
    /// Gene has high CpG methylation (>70% mean at annotated sites) but is
    /// actively expressed (mean TPM > 10). Promoter methylation should silence.
    /// Possible: alternate promoter, enhancer bypass, NMD escape.
    MethylatedButExpressed,

    /// Gene carries a high-quality variant (QUAL > 30) and is in the top 10%
    /// expressed genes. Variant may affect a functionally active protein.
    /// Possible: gain-of-function, dominant negative, haploinsufficiency.
    VariantInActiveGene,

    /// Gene carries a high-quality variant but is very lowly expressed
    /// (mean TPM < 1). Variant may have caused transcriptional silencing,
    /// or the variant affects a tissue-specific gene not expressed here.
    VariantInSilentGene,

    /// Gene is significantly differentially expressed (|log2FC| > 2, padj < 0.05)
    /// but has no high-impact variant detected → epigenetic or post-transcriptional
    /// mechanism is the likely driver (e.g., promoter methylation, miRNA, splicing).
    DifferentialWithoutVariant,

    /// Gene has a high-impact variant but shows no differential expression
    /// (|log2FC| < 0.5 when DE data is present). Variant may be a passenger
    /// or compensated by the other allele.
    VariantWithoutExpression,

    /// Gene is abnormal in 2+ modalities simultaneously (any combination of
    /// variant + expression change + methylation change). Based on the
    /// Knudson two-hit hypothesis; these are highest-priority candidates.
    MultiHit { n_modalities: u8 },
}

/// A gene flagged with a biological paradox or multi-modal convergence signal.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneParadox {
    pub gene: String,
    pub kind: ParadoxKind,
    /// Supporting evidence values (TPM, methylation%, QUAL, log2FC, etc.)
    pub evidence: ParadoxEvidence,
    /// Plain-English one-line summary for the report.
    pub summary: String,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ParadoxEvidence {
    pub mean_tpm: Option<f64>,
    pub mean_methylation: Option<f64>,
    pub max_variant_qual: Option<f32>,
    pub log2_fold_change: Option<f64>,
    pub padj: Option<f64>,
    pub n_modalities_hit: u8,
}

/// Detect multi-modal paradoxes across all three modality summaries.
///
/// Returns a list of paradoxical genes sorted by priority:
/// MultiHit first, then by number of modalities, then alphabetically.
pub fn detect_paradoxes(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
) -> Vec<GeneParadox> {
    // 1. Build set of high-impact variant genes
    let variant_genes: HashSet<String> = genomics.high_impact_genes.iter().cloned().collect();

    // Build a map from gene → max QUAL among high-impact variants
    let mut gene_max_qual: HashMap<&str, f32> = HashMap::new();
    for v in &genomics.high_impact {
        if let Some(ref g) = v.gene {
            let entry = gene_max_qual.entry(g.as_str()).or_insert(0.0_f32);
            if v.qual > *entry {
                *entry = v.qual;
            }
        }
    }

    // 2. Compute 90th-percentile TPM threshold from gene_stats
    let mut all_means: Vec<f64> = transcr
        .gene_stats
        .values()
        .map(|s| s.mean)
        .collect();
    all_means.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p90_tpm = if all_means.len() < 10 {
        // Use median if fewer than 10 genes
        if all_means.is_empty() {
            0.0
        } else {
            all_means[all_means.len() / 2]
        }
    } else {
        let idx = (0.9 * all_means.len() as f64) as usize;
        let idx = idx.min(all_means.len() - 1);
        all_means[idx]
    };

    // 3. Build map of significantly DE genes: gene → DiffExprResult reference
    // Criteria: padj < 0.05 && |log2FC| > 2.0 && !padj.is_nan()
    let de_sig: HashMap<String, (f64, f64)> = if let Some(ref de) = transcr.diff_expr {
        de.iter()
            .filter(|r| {
                !r.padj.is_nan()
                    && r.padj < 0.05
                    && r.log2_fold_change.abs() > 2.0
            })
            .map(|r| (r.gene_id.clone(), (r.log2_fold_change, r.padj)))
            .collect()
    } else {
        HashMap::new()
    };

    // Also build a full DE map for VariantWithoutExpression check (|log2FC| < 0.5)
    let de_all: HashMap<String, f64> = if let Some(ref de) = transcr.diff_expr {
        de.iter()
            .filter(|r| !r.log2_fold_change.is_nan())
            .map(|r| (r.gene_id.clone(), r.log2_fold_change))
            .collect()
    } else {
        HashMap::new()
    };

    // 4. Methylation lookup
    let meth_map = &epigen.gene_methylation;

    // 5. Collect all genes across all data sources
    let mut all_genes: HashSet<String> = HashSet::new();
    for g in &genomics.high_impact_genes {
        all_genes.insert(g.clone());
    }
    for g in transcr.gene_stats.keys() {
        all_genes.insert(g.clone());
    }
    for g in de_sig.keys() {
        all_genes.insert(g.clone());
    }
    for g in meth_map.keys() {
        all_genes.insert(g.clone());
    }

    // 6. Collect all paradoxes
    let mut paradoxes: Vec<GeneParadox> = Vec::new();
    // Track genes already flagged as MultiHit to avoid duplicates
    let mut multi_hit_genes: HashSet<String> = HashSet::new();

    for gene in &all_genes {
        let mean_tpm = transcr.gene_stats.get(gene).map(|s| s.mean);
        let mean_meth = meth_map.get(gene).copied();
        let max_qual = gene_max_qual.get(gene.as_str()).copied();
        let in_variants = variant_genes.contains(gene);
        let de_result = de_sig.get(gene);
        let de_log2fc = de_result.map(|(l, _)| *l);
        let de_padj = de_result.map(|(_, p)| *p);

        // Compute n_modalities for MultiHit
        let mut n_modalities: u8 = 0;
        if in_variants {
            n_modalities += 1;
        }
        let tpm_active = mean_tpm.map(|t| t > 10.0).unwrap_or(false);
        if tpm_active || de_result.is_some() {
            n_modalities += 1;
        }
        if mean_meth.map(|m| m > 70.0).unwrap_or(false) {
            n_modalities += 1;
        }

        let base_evidence = ParadoxEvidence {
            mean_tpm,
            mean_methylation: mean_meth,
            max_variant_qual: max_qual,
            log2_fold_change: de_log2fc,
            padj: de_padj,
            n_modalities_hit: n_modalities,
        };

        // MethylatedButExpressed
        if let (Some(meth), Some(tpm)) = (mean_meth, mean_tpm) {
            if meth > 70.0 && tpm > 10.0 {
                let summary = format!(
                    "{gene}: mean methylation={meth:.1}% but mean TPM={tpm:.1} \
                     — promoter CpG methylation should silence expression",
                );
                paradoxes.push(GeneParadox {
                    gene: gene.clone(),
                    kind: ParadoxKind::MethylatedButExpressed,
                    evidence: base_evidence.clone(),
                    summary,
                });
            }
        }

        // VariantInActiveGene
        if in_variants {
            if let Some(tpm) = mean_tpm {
                if tpm >= p90_tpm {
                    let pct_rank = 10u8; // top 10%
                    let qual_str = max_qual.map(|q| format!("{q:.0}")).unwrap_or_else(|| ">30".to_string());
                    let summary = format!(
                        "{gene}: high-impact variant (QUAL={qual_str}) in actively expressed gene \
                         (TPM={tpm:.1}, top {pct_rank}%)",
                    );
                    paradoxes.push(GeneParadox {
                        gene: gene.clone(),
                        kind: ParadoxKind::VariantInActiveGene,
                        evidence: base_evidence.clone(),
                        summary,
                    });
                }

                // VariantInSilentGene
                if tpm < 1.0 {
                    let qual_str = max_qual.map(|q| format!("{q:.0}")).unwrap_or_else(|| ">30".to_string());
                    let summary = format!(
                        "{gene}: high-impact variant (QUAL={qual_str}) but gene is silenced \
                         (mean TPM={tpm:.2}) — variant may have caused transcriptional silencing",
                    );
                    paradoxes.push(GeneParadox {
                        gene: gene.clone(),
                        kind: ParadoxKind::VariantInSilentGene,
                        evidence: base_evidence.clone(),
                        summary,
                    });
                }
            }
        }

        // DifferentialWithoutVariant
        if let Some((lfc, pj)) = de_result {
            if !in_variants {
                let summary = format!(
                    "{gene}: |log2FC|={lfc:.2} (padj={pj:.2e}) but no high-impact variant detected \
                     — epigenetic or post-transcriptional driver likely",
                    lfc = lfc.abs(),
                );
                paradoxes.push(GeneParadox {
                    gene: gene.clone(),
                    kind: ParadoxKind::DifferentialWithoutVariant,
                    evidence: base_evidence.clone(),
                    summary,
                });
            }
        }

        // VariantWithoutExpression: gene in high_impact AND DE data exists AND |log2FC| < 0.5
        if in_variants && !de_all.is_empty() {
            if let Some(&lfc) = de_all.get(gene) {
                if lfc.abs() < 0.5 {
                    let qual_str = max_qual.map(|q| format!("{q:.0}")).unwrap_or_else(|| ">30".to_string());
                    let summary = format!(
                        "{gene}: high-impact variant (QUAL={qual_str}) but no differential \
                         expression (|log2FC|={lfc:.2}) — variant may be a passenger or \
                         functionally compensated",
                        lfc = lfc.abs(),
                    );
                    paradoxes.push(GeneParadox {
                        gene: gene.clone(),
                        kind: ParadoxKind::VariantWithoutExpression,
                        evidence: base_evidence.clone(),
                        summary,
                    });
                }
            }
        }

        // MultiHit: n_modalities >= 2
        if n_modalities >= 2 && !multi_hit_genes.contains(gene) {
            multi_hit_genes.insert(gene.clone());
            let summary = format!(
                "{gene}: abnormal in {n_modalities} modalities \
                 — high-priority multi-hit candidate",
            );
            paradoxes.push(GeneParadox {
                gene: gene.clone(),
                kind: ParadoxKind::MultiHit { n_modalities },
                evidence: base_evidence,
                summary,
            });
        }
    }

    // 8. Sort: MultiHit first (descending n_modalities), then others, then alphabetical
    paradoxes.sort_by(|a, b| {
        let a_multi = matches!(a.kind, ParadoxKind::MultiHit { .. });
        let b_multi = matches!(b.kind, ParadoxKind::MultiHit { .. });

        match (a_multi, b_multi) {
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            (true, true) => {
                // Both MultiHit — sort by descending n_modalities, then gene name
                let a_n = if let ParadoxKind::MultiHit { n_modalities } = a.kind {
                    n_modalities
                } else {
                    0
                };
                let b_n = if let ParadoxKind::MultiHit { n_modalities } = b.kind {
                    n_modalities
                } else {
                    0
                };
                b_n.cmp(&a_n).then_with(|| a.gene.cmp(&b.gene))
            }
            (false, false) => a.gene.cmp(&b.gene),
        }
    });

    paradoxes
}
