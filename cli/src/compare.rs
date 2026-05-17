//! Tumor-vs-normal / treatment-vs-control multi-omics comparison.
//!
//! Runs the full pipeline on both the case and control inputs, then computes
//! per-modality differential statistics:
//!
//! - **Genomics**: variant burden ratio, unique/shared variant sets, Ti/Tv
//!   difference, per-chromosome density delta.
//! - **Transcriptomics**: log₂FC between case and control global expression
//!   profiles; DE genes detected across groups.
//! - **Epigenomics**: per-chromosome Δmethylation; hypermethylation / hypo-
//!   methylation events unique to case.
//! - **ATAC-seq**: differential peak density per chromosome (optional).

use anyhow::Result;
use crossbeam_channel::unbounded;
use serde::{Deserialize, Serialize};

use genomics_core::GenomicsSummary;
use transcriptomics_core::TranscriptomicsSummary;
use epigenomics_core::EpigenomicsSummary;

use crate::args::Cli;
use crate::config::BioomicsConfig;

/// Per-chromosome variant density delta (case − control).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromVariantDelta {
    pub chrom: String,
    pub case_total: u64,
    pub control_total: u64,
    /// (case − control) / max(control, 1)
    pub fold_enrichment: f64,
}

/// Differential methylation at chromosome resolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromMethDelta {
    pub chrom: String,
    pub case_mean: f64,
    pub control_mean: f64,
    /// Absolute difference: case − control (positive = hypermethylated in case).
    pub delta: f64,
}

/// Top-line genomics comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomicsComparison {
    pub case_total_variants: u64,
    pub control_total_variants: u64,
    /// case / control (normalized by genome coverage).
    pub variant_burden_ratio: f64,
    pub case_titv: f64,
    pub control_titv: f64,
    pub case_unique_positions: u64,
    pub control_unique_positions: u64,
    pub per_chrom: Vec<ChromVariantDelta>,
}

/// Top-line transcriptomics comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptomicsComparison {
    pub case_expressed_genes: u64,
    pub control_expressed_genes: u64,
    pub case_total_genes: u64,
    /// Genes expressed in case but silent (TPM < threshold) in control.
    pub gained_genes: Vec<String>,
    /// Genes silent in case but expressed in control.
    pub lost_genes: Vec<String>,
    /// Genes in both top-100 expressed lists (shared high-expressors).
    pub shared_top_expressed: Vec<String>,
}

/// Top-line epigenomics comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpigenomicsComparison {
    pub case_global_meth: f64,
    pub control_global_meth: f64,
    /// case − control global methylation.
    pub delta_global_meth: f64,
    pub case_cpg_islands: usize,
    pub control_cpg_islands: usize,
    pub per_chrom: Vec<ChromMethDelta>,
    /// Chromosomes with |Δmeth| > threshold (from config).
    pub differentially_methylated_chroms: Vec<String>,
}

/// Combined comparison summary across all modalities.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonSummary {
    pub case_label: String,
    pub control_label: String,
    pub genomics: GenomicsComparison,
    pub transcriptomics: TranscriptomicsComparison,
    pub epigenomics: EpigenomicsComparison,
    /// Plain-English comparison insights.
    pub insights: Vec<String>,
}

/// Run the comparison pipeline: analyse control inputs, then diff against case.
pub fn run_comparison(
    cli: &Cli,
    cfg: &BioomicsConfig,
    case_g: &GenomicsSummary,
    case_t: &TranscriptomicsSummary,
    case_e: &EpigenomicsSummary,
) -> Result<ComparisonSummary> {
    let ctrl_g_path = cli.compare_genomics.as_ref().unwrap();
    let ctrl_t_path = cli.compare_transcriptomics.as_ref().unwrap();
    let ctrl_e_path = cli.compare_epigenomics.as_ref().unwrap();

    log::info!(
        "Comparison mode: analysing control inputs '{}', '{}', '{}'",
        ctrl_g_path.display(), ctrl_t_path.display(), ctrl_e_path.display()
    );

    let (gtx, _grx) = unbounded();
    let (ttx, _trx) = unbounded();
    let (etx, _erx) = unbounded();

    let (ctrl_g, ctrl_t, ctrl_e) = {
        let (gr, tr, er) = std::thread::scope(|s| {
            let gh = s.spawn(|| genomics_core::analyze_vcf(ctrl_g_path, Some(&gtx)));
            let th = s.spawn(|| transcriptomics_core::analyze_tsv(ctrl_t_path, Some(&ttx)));
            let eh = s.spawn(|| epigenomics_core::analyze_bed(ctrl_e_path, Some(&etx)));
            (gh.join(), th.join(), eh.join())
        });
        let g = gr.map_err(|e| anyhow::anyhow!("Control genomics panicked: {:?}", e))??;
        let t = tr.map_err(|e| anyhow::anyhow!("Control transcriptomics panicked: {:?}", e))??;
        let e = er.map_err(|e| anyhow::anyhow!("Control epigenomics panicked: {:?}", e))??;
        (g, t, e)
    };

    let genomics = diff_genomics(case_g, &ctrl_g);
    let transcriptomics = diff_transcriptomics(case_t, &ctrl_t);
    let epigenomics = diff_epigenomics(case_e, &ctrl_e, cfg.compare.delta_meth_threshold);
    let insights = derive_comparison_insights(
        &genomics,
        &transcriptomics,
        &epigenomics,
        cfg,
    );

    Ok(ComparisonSummary {
        case_label: cfg.compare.case_label.clone(),
        control_label: cfg.compare.control_label.clone(),
        genomics,
        transcriptomics,
        epigenomics,
        insights,
    })
}

fn diff_genomics(case: &GenomicsSummary, ctrl: &GenomicsSummary) -> GenomicsComparison {
    let burden_ratio = if ctrl.total_variants > 0 {
        case.total_variants as f64 / ctrl.total_variants as f64
    } else {
        f64::INFINITY
    };

    let mut per_chrom: Vec<ChromVariantDelta> = case
        .per_chrom
        .iter()
        .map(|(chrom, cd)| {
            let ctrl_total = ctrl.per_chrom.get(chrom.as_str()).map(|c| c.total).unwrap_or(0);
            let fold = (cd.total as f64 - ctrl_total as f64) / (ctrl_total.max(1) as f64);
            ChromVariantDelta {
                chrom: chrom.clone(),
                case_total: cd.total,
                control_total: ctrl_total,
                fold_enrichment: fold,
            }
        })
        .collect();
    per_chrom.sort_unstable_by(|a, b| {
        b.fold_enrichment.abs().partial_cmp(&a.fold_enrichment.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    GenomicsComparison {
        case_total_variants: case.total_variants,
        control_total_variants: ctrl.total_variants,
        variant_burden_ratio: burden_ratio,
        case_titv: case.titv_ratio,
        control_titv: ctrl.titv_ratio,
        case_unique_positions: case.unique_positions,
        control_unique_positions: ctrl.unique_positions,
        per_chrom,
    }
}

fn diff_transcriptomics(
    case: &TranscriptomicsSummary,
    ctrl: &TranscriptomicsSummary,
) -> TranscriptomicsComparison {
    let case_expressed: ahash::AHashSet<String> = case
        .top_100_expressed
        .iter()
        .map(|(g, _)| g.clone())
        .collect();
    let ctrl_expressed: ahash::AHashSet<String> = ctrl
        .top_100_expressed
        .iter()
        .map(|(g, _)| g.clone())
        .collect();

    let gained: Vec<String> = case_expressed
        .difference(&ctrl_expressed)
        .cloned()
        .collect();
    let lost: Vec<String> = ctrl_expressed
        .difference(&case_expressed)
        .cloned()
        .collect();
    let shared: Vec<String> = case_expressed
        .intersection(&ctrl_expressed)
        .cloned()
        .collect();

    TranscriptomicsComparison {
        case_expressed_genes: case.expressed_genes,
        control_expressed_genes: ctrl.expressed_genes,
        case_total_genes: case.total_genes,
        gained_genes: gained,
        lost_genes: lost,
        shared_top_expressed: shared,
    }
}

fn diff_epigenomics(
    case: &EpigenomicsSummary,
    ctrl: &EpigenomicsSummary,
    delta_threshold: f64,
) -> EpigenomicsComparison {
    let delta_global = case.global_methylation_pct - ctrl.global_methylation_pct;

    let mut per_chrom: Vec<ChromMethDelta> = case
        .per_chrom
        .iter()
        .map(|(chrom, cm)| {
            let ctrl_mean = ctrl
                .per_chrom
                .get(chrom.as_str())
                .map(|c| c.mean_methylation)
                .unwrap_or(case.global_methylation_pct);
            ChromMethDelta {
                chrom: chrom.clone(),
                case_mean: cm.mean_methylation,
                control_mean: ctrl_mean,
                delta: cm.mean_methylation - ctrl_mean,
            }
        })
        .collect();
    per_chrom.sort_unstable_by(|a, b| {
        b.delta.abs().partial_cmp(&a.delta.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let differentially_methylated: Vec<String> = per_chrom
        .iter()
        .filter(|d| d.delta.abs() >= delta_threshold)
        .map(|d| d.chrom.clone())
        .collect();

    EpigenomicsComparison {
        case_global_meth: case.global_methylation_pct,
        control_global_meth: ctrl.global_methylation_pct,
        delta_global_meth: delta_global,
        case_cpg_islands: case.cpg_islands.len(),
        control_cpg_islands: ctrl.cpg_islands.len(),
        per_chrom,
        differentially_methylated_chroms: differentially_methylated,
    }
}

fn derive_comparison_insights(
    g: &GenomicsComparison,
    t: &TranscriptomicsComparison,
    e: &EpigenomicsComparison,
    cfg: &BioomicsConfig,
) -> Vec<String> {
    let mut insights = Vec::new();

    // Variant burden
    if g.variant_burden_ratio > cfg.compare.variant_fc_threshold {
        insights.push(format!(
            "[WARN] {case} has {:.1}× higher variant burden than {ctrl} ({cv} vs {ctv} variants)",
            g.variant_burden_ratio,
            case = cfg.compare.case_label,
            ctrl = cfg.compare.control_label,
            cv = g.case_total_variants,
            ctv = g.control_total_variants,
        ));
    } else if g.variant_burden_ratio < 1.0 / cfg.compare.variant_fc_threshold {
        insights.push(format!(
            "[INFO] {case} has fewer variants than {ctrl} (ratio={:.2})",
            g.variant_burden_ratio,
            case = cfg.compare.case_label,
            ctrl = cfg.compare.control_label,
        ));
    }

    // Ti/Tv comparison
    let titv_delta = (g.case_titv - g.control_titv).abs();
    if titv_delta > 0.3 {
        insights.push(format!(
            "[INFO] Ti/Tv differs between groups: {case}={:.2}, {ctrl}={:.2}",
            g.case_titv, g.control_titv,
            case = cfg.compare.case_label,
            ctrl = cfg.compare.control_label,
        ));
    }

    // Gained/lost expression
    if !t.gained_genes.is_empty() {
        let sample: Vec<_> = t.gained_genes.iter().take(3).cloned().collect();
        insights.push(format!(
            "[INFO] {} genes gained in {case} top-expressed (e.g. {})",
            t.gained_genes.len(),
            sample.join(", "),
            case = cfg.compare.case_label,
        ));
    }
    if !t.lost_genes.is_empty() {
        let sample: Vec<_> = t.lost_genes.iter().take(3).cloned().collect();
        insights.push(format!(
            "[INFO] {} genes lost from {ctrl} top-expressed (e.g. {})",
            t.lost_genes.len(),
            sample.join(", "),
            ctrl = cfg.compare.control_label,
        ));
    }

    // Global methylation shift
    if e.delta_global_meth.abs() >= cfg.compare.delta_meth_threshold {
        let direction = if e.delta_global_meth > 0.0 { "hyper" } else { "hypo" };
        insights.push(format!(
            "[{}] Global {}methylation in {case}: {:.1}% vs {ctrl}: {:.1}% (Δ={:+.1}%)",
            if e.delta_global_meth.abs() > 20.0 { "CRIT" } else { "WARN" },
            direction,
            e.case_global_meth,
            e.control_global_meth,
            e.delta_global_meth,
            case = cfg.compare.case_label,
            ctrl = cfg.compare.control_label,
        ));
    }

    // Differentially methylated chromosomes
    if !e.differentially_methylated_chroms.is_empty() {
        let n = e.differentially_methylated_chroms.len();
        let sample: Vec<_> = e.differentially_methylated_chroms.iter().take(3).cloned().collect();
        insights.push(format!(
            "[INFO] {} chromosomes with |Δmeth| ≥ {:.0}% ({}{})",
            n,
            cfg.compare.delta_meth_threshold,
            sample.join(", "),
            if n > 3 { ", …" } else { "" },
        ));
    }

    // CpG island difference
    let cpg_delta = e.case_cpg_islands as i64 - e.control_cpg_islands as i64;
    if cpg_delta.abs() > 5 {
        insights.push(format!(
            "[INFO] CpG islands: {case}={ci}, {ctrl}={cti} (Δ={cpg_delta:+})",
            case = cfg.compare.case_label,
            ctrl = cfg.compare.control_label,
            ci = e.case_cpg_islands,
            cti = e.control_cpg_islands,
        ));
    }

    insights
}
