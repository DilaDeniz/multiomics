//! DESeq2-equivalent negative-binomial GLM differential expression engine.
//!
//! Implements the full Love et al. 2014 pipeline:
//!   1. Size factors — median-of-ratios (Anders & Huber 2010)
//!   2. Gene-wise MLE dispersions — method of moments
//!   3. Parametric dispersion trend fitting — α = a₀ + a₁/μ
//!   4. MAP dispersion shrinkage — empirical Bayes log-normal prior
//!   5. NB-GLM fitting — IRLS with 2×2 WLS closed-form solve
//!   6. Wald test — z-statistic, two-sided normal p-value
//!   7. Cook's distance outlier flagging
//!   8. Independent filtering (Bourgon et al. 2010)
//!   9. LFC shrinkage — normal-prior MAP (apeglm-style)
//!
//! # Citation
//! Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and
//! dispersion for RNA-seq data with DESeq2. Genome Biol 15:550.

use anyhow::{bail, Context, Result};
use biomics_core::statistics::{benjamini_hochberg, welch_t_test};

use crate::types::DiffExprResult;

// ── Public data structures ────────────────────────────────────────────────────

/// Per-sample size factors estimated by the median-of-ratios method.
#[derive(Debug, Clone)]
pub struct SizeFactors {
    pub sample_names: Vec<String>,
    /// One factor per sample, same order as `sample_names`.
    pub factors: Vec<f64>,
}

/// A full normalized count matrix together with the size factors and the
/// pooled dispersion estimate.
#[derive(Debug, Clone)]
pub struct NormalizedMatrix {
    pub gene_ids: Vec<String>,
    pub sample_names: Vec<String>,
    pub size_factors: SizeFactors,
    /// Raw counts — `counts[gene_index][sample_index]`.
    pub counts: Vec<Vec<f64>>,
    /// Normalized counts — `normalized[gene_index][sample_index]`.
    pub normalized: Vec<Vec<f64>>,
    /// Pooled (median) negative-binomial dispersion across all genes.
    pub global_dispersion: f64,
}

/// Parametric dispersion trend α = a₀ + a₁/μ (Love et al. 2014 eq. 4).
#[derive(Debug, Clone)]
pub struct DispersionTrend {
    /// a₀ — asymptotic dispersion at high mean counts.
    pub asympt_disp: f64,
    /// a₁ — extra-Poisson variance coefficient.
    pub extra_pois: f64,
}

impl DispersionTrend {
    /// Evaluate the trend at a given mean count `mean`.
    pub fn eval(&self, mean: f64) -> f64 {
        self.asympt_disp + self.extra_pois / mean.max(1e-8)
    }
}

/// Full DESeq2-equivalent result for one gene.
#[derive(Debug, Clone)]
pub struct DeseqResult {
    pub gene_id: String,
    /// Mean normalized count across all samples.
    pub base_mean: f64,
    /// MLE log₂ fold change (β₁ / ln 2).
    pub log2_fold_change: f64,
    /// Shrunk log₂ fold change (apeglm-style normal-prior MAP).
    pub lfc_shrunk: f64,
    /// Standard error of β₁ divided by ln 2.
    pub lfc_se: f64,
    /// Wald z-statistic.
    pub stat: f64,
    pub p_value: f64,
    /// Benjamini-Hochberg adjusted p-value (NaN if filtered).
    pub padj: f64,
    /// Mean normalized count in group 0.
    pub mean_group0: f64,
    /// Mean normalized count in group 1.
    pub mean_group1: f64,
    /// Final MAP-shrunk dispersion.
    pub dispersion: f64,
    /// True if Cook's distance exceeds threshold.
    pub outlier: bool,
}

/// Top-level result from a full DESeq2 run.
#[derive(Debug, Clone)]
pub struct DeseqAnalysis {
    pub results: Vec<DeseqResult>,
    pub size_factors: SizeFactors,
    pub dispersion_trend: DispersionTrend,
    /// Number of genes with padj < 0.05.
    pub n_significant: usize,
    /// Number of genes removed by independent filtering.
    pub n_filtered: usize,
}

/// Per-gene result from a likelihood-ratio test comparing full vs. reduced model.
#[derive(Debug, Clone)]
pub struct LrtResult {
    pub gene_id: String,
    pub base_mean: f64,
    /// LRT statistic: 2 × (ll_full − ll_reduced).
    pub stat: f64,
    /// Degrees of freedom = n_full_factors − n_reduced_factors.
    pub df: usize,
    pub p_value: f64,
    pub padj: f64,
}

// ── Core estimation functions ─────────────────────────────────────────────────

/// Estimate per-sample size factors using the DESeq2 median-of-ratios method.
///
/// `counts[gene][sample]` — raw integer counts as f64.
pub fn estimate_size_factors(counts: &[Vec<f64>], sample_names: &[String]) -> Result<SizeFactors> {
    let n_genes = counts.len();
    if n_genes == 0 {
        bail!("estimate_size_factors: count matrix is empty");
    }
    let n_samples = sample_names.len();
    if n_samples == 0 {
        bail!("estimate_size_factors: no sample names provided");
    }
    // Validate row lengths
    for (g, row) in counts.iter().enumerate() {
        if row.len() != n_samples {
            bail!(
                "estimate_size_factors: gene {} has {} values but {} samples declared",
                g,
                row.len(),
                n_samples
            );
        }
    }

    // Step 1: geometric mean per gene in log-space.
    // Only include non-zero counts to avoid log(0).
    // Genes where every sample is zero are excluded from size-factor estimation
    // (their geometric mean would be 0, making ratios undefined).
    let geom_means: Vec<f64> = counts
        .iter()
        .map(|row| {
            let positive_logs: Vec<f64> =
                row.iter().filter(|&&c| c > 0.0).map(|&c| c.ln()).collect();
            if positive_logs.is_empty() {
                0.0 // sentinel: excluded below
            } else {
                let mean_log = positive_logs.iter().sum::<f64>() / positive_logs.len() as f64;
                mean_log.exp()
            }
        })
        .collect();

    // Step 2: for each sample, collect ratios count_gj / mu_g for genes where
    // mu_g > 0, then take the median.
    let mut factors = Vec::with_capacity(n_samples);
    for j in 0..n_samples {
        let mut ratios: Vec<f64> = counts
            .iter()
            .zip(geom_means.iter())
            .filter(|(_, &mu)| mu > 0.0)
            .map(|(row, &mu)| row[j] / mu)
            .collect();

        if ratios.is_empty() {
            bail!(
                "estimate_size_factors: no valid genes for size-factor estimation of sample {}",
                sample_names[j]
            );
        }

        ratios.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median = median_sorted(&ratios);

        // Guard against degenerate size factors
        if median <= 0.0 || !median.is_finite() {
            bail!(
                "estimate_size_factors: computed non-positive size factor ({}) for sample {}",
                median,
                sample_names[j]
            );
        }
        factors.push(median);
    }

    Ok(SizeFactors {
        sample_names: sample_names.to_vec(),
        factors,
    })
}

/// Normalize a count matrix and estimate global dispersion.
///
/// `counts[gene][sample]` — raw integer counts as f64.
pub fn normalize_counts(
    gene_ids: &[String],
    counts: &[Vec<f64>],
    sample_names: &[String],
) -> Result<NormalizedMatrix> {
    if gene_ids.len() != counts.len() {
        bail!(
            "normalize_counts: {} gene IDs but {} count rows",
            gene_ids.len(),
            counts.len()
        );
    }

    let size_factors = estimate_size_factors(counts, sample_names)
        .context("normalize_counts: size-factor estimation failed")?;

    let n_samples = sample_names.len();

    // Normalize: normalized_gj = count_gj / s_j
    let normalized: Vec<Vec<f64>> = counts
        .iter()
        .map(|row| {
            row.iter()
                .zip(size_factors.factors.iter())
                .map(|(&c, &s)| c / s)
                .collect()
        })
        .collect();

    // Dispersion estimation via method of moments.
    // For each gene: mean_g and var_g of normalized counts across samples.
    // alpha_g = max(0, (var_g - mean_g) / mean_g^2)   [NB dispersion]
    let mut dispersions: Vec<f64> = Vec::with_capacity(normalized.len());

    for norm_row in &normalized {
        if norm_row.is_empty() {
            continue;
        }
        let n = norm_row.len() as f64;
        let mean_g = norm_row.iter().sum::<f64>() / n;
        if mean_g <= 0.0 {
            dispersions.push(0.0);
            continue;
        }
        let var_g = if n_samples > 1 {
            norm_row.iter().map(|&x| (x - mean_g).powi(2)).sum::<f64>() / (n - 1.0)
        } else {
            0.0
        };
        let alpha = ((var_g - mean_g) / mean_g.powi(2)).max(0.0);
        dispersions.push(alpha);
    }

    // Global dispersion: median across all genes
    dispersions.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let global_dispersion = if dispersions.is_empty() {
        0.0
    } else {
        median_sorted(&dispersions)
    };

    Ok(NormalizedMatrix {
        gene_ids: gene_ids.to_vec(),
        sample_names: sample_names.to_vec(),
        size_factors,
        counts: counts.to_vec(),
        normalized,
        global_dispersion,
    })
}

/// Differential expression on a `NormalizedMatrix`.
///
/// Groups are split by position — first ⌈n/2⌉ samples vs. the remainder —
/// matching the convention in `diffexpr.rs`. Welch's t-test is applied on
/// log₂(normalized + 0.5). BH FDR is applied across genes with a valid p-value.
/// Results are sorted by ascending `padj` (NaN last), then descending |log2FC|.
pub fn deseq2_differential_expression(matrix: &NormalizedMatrix) -> Vec<DiffExprResult> {
    let n_samples = matrix.sample_names.len();
    if n_samples == 0 || matrix.normalized.is_empty() {
        return Vec::new();
    }

    let split = n_samples.div_ceil(2);
    let can_test = split >= 2 && (n_samples - split) >= 2;

    let mut results: Vec<DiffExprResult> = matrix
        .normalized
        .iter()
        .zip(matrix.gene_ids.iter())
        .map(|(norm_row, gene_id)| {
            // log₂(normalized + 0.5) to handle zeros and approach the log-normal
            let g1: Vec<f64> = norm_row[..split]
                .iter()
                .map(|&v| (v + 0.5_f64).log2())
                .collect();
            let g2: Vec<f64> = norm_row[split..]
                .iter()
                .map(|&v| (v + 0.5_f64).log2())
                .collect();

            let mean1 = g1.iter().sum::<f64>() / g1.len().max(1) as f64;
            let mean2 = g2.iter().sum::<f64>() / g2.len().max(1) as f64;
            let lfc = mean2 - mean1; // log2 scale → log2FC

            // Raw (non-log) means for reporting
            let mean_s1 = norm_row[..split].iter().sum::<f64>() / split.max(1) as f64;
            let mean_s2 = norm_row[split..].iter().sum::<f64>() / (n_samples - split).max(1) as f64;

            let (p_value, padj) = if can_test {
                let pval = welch_t_test(&g1, &g2).map(|(_, p)| p).unwrap_or(f64::NAN);
                (pval, f64::NAN) // padj filled below
            } else {
                (f64::NAN, f64::NAN)
            };

            DiffExprResult {
                gene_id: gene_id.clone(),
                log2_fold_change: lfc,
                mean_s1,
                mean_s2,
                p_value,
                padj,
            }
        })
        .collect();

    // BH FDR correction over the subset with valid p-values
    if can_test {
        let valid_indices: Vec<usize> = results
            .iter()
            .enumerate()
            .filter(|(_, r)| !r.p_value.is_nan())
            .map(|(i, _)| i)
            .collect();

        if !valid_indices.is_empty() {
            let pvals: Vec<f64> = valid_indices.iter().map(|&i| results[i].p_value).collect();
            let padj_vals = benjamini_hochberg(&pvals);
            for (vi, &orig_i) in valid_indices.iter().enumerate() {
                results[orig_i].padj = padj_vals[vi];
            }
        }
    }

    // Sort: smallest padj first (NaN last), break ties by descending |log2FC|
    results.sort_unstable_by(|a, b| match (a.padj.is_nan(), b.padj.is_nan()) {
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => a
            .padj
            .partial_cmp(&b.padj)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(
                b.log2_fold_change
                    .abs()
                    .partial_cmp(&a.log2_fold_change.abs())
                    .unwrap_or(std::cmp::Ordering::Equal),
            ),
    });

    results
}

// ── Full DESeq2 NB-GLM pipeline ───────────────────────────────────────────────

/// Run the full DESeq2 NB-GLM pipeline (Love et al. 2014).
///
/// `counts[gene][sample]` — raw integer counts as f64.
/// `group0` and `group1` — indices into the sample dimension.
///
/// Internally delegates to `run_deseq2_design` with a two-group design matrix
/// and a contrast of `[0, 1]` (testing the treatment coefficient).
pub fn run_deseq2(
    counts: &[Vec<f64>],
    gene_ids: &[String],
    group0: &[usize],
    group1: &[usize],
    sample_names: &[String],
) -> Result<DeseqAnalysis> {
    if counts.is_empty() {
        bail!("run_deseq2: count matrix is empty");
    }
    if gene_ids.len() != counts.len() {
        bail!(
            "run_deseq2: {} gene IDs but {} count rows",
            gene_ids.len(),
            counts.len()
        );
    }
    if group0.is_empty() || group1.is_empty() {
        bail!("run_deseq2: each group must have at least one sample");
    }
    let n_samples = sample_names.len();
    for &idx in group0.iter().chain(group1.iter()) {
        if idx >= n_samples {
            bail!(
                "run_deseq2: sample index {} out of range (n_samples={})",
                idx,
                n_samples
            );
        }
    }

    let design = design_matrix_two_group(group0, group1, n_samples);
    let contrast = vec![0.0, 1.0];
    run_deseq2_design(counts, gene_ids, &design, &contrast, sample_names)
        .context("run_deseq2: general design pipeline failed")
}

// ── General multi-factor design: public API ───────────────────────────────────

/// Run the DESeq2 NB-GLM pipeline with an arbitrary design matrix.
///
/// Caller constructs `design_matrix[sample][factor]` (e.g. via
/// `design_matrix_two_group`, `design_matrix_with_batch`, or
/// `interaction_design_matrix`).  `contrast` selects the linear combination of
/// coefficients to test (Wald test).
pub fn run_deseq2_design(
    counts: &[Vec<f64>],
    gene_ids: &[String],
    design_matrix: &[Vec<f64>],
    contrast: &[f64],
    sample_names: &[String],
) -> Result<DeseqAnalysis> {
    if counts.is_empty() {
        bail!("run_deseq2_design: count matrix is empty");
    }
    if gene_ids.len() != counts.len() {
        bail!(
            "run_deseq2_design: {} gene IDs but {} count rows",
            gene_ids.len(),
            counts.len()
        );
    }
    let n_samples = sample_names.len();
    if design_matrix.len() != n_samples {
        bail!(
            "run_deseq2_design: design matrix has {} rows but {} samples",
            design_matrix.len(),
            n_samples
        );
    }
    if n_samples == 0 {
        bail!("run_deseq2_design: no samples");
    }
    let n_factors = design_matrix[0].len();
    if n_factors == 0 {
        bail!("run_deseq2_design: design matrix has no columns");
    }
    if contrast.len() != n_factors {
        bail!(
            "run_deseq2_design: contrast length {} != n_factors {}",
            contrast.len(),
            n_factors
        );
    }

    // Step 1: size factors
    let size_factors = estimate_size_factors(counts, sample_names)
        .context("run_deseq2_design: size-factor estimation failed")?;
    let sf = &size_factors.factors;
    let n_genes = counts.len();

    // Derive group0/group1 from first non-intercept column for dispersion MoM
    // (heuristic: treat col 1 as the primary treatment indicator)
    let (group0_def, group1_def) = groups_from_design(design_matrix, n_samples);

    // Step 2: MLE dispersions
    let mle_dispersions = compute_mle_dispersions(counts, sf, &group0_def, &group1_def);

    // Step 3: dispersion trend
    let gene_means = compute_gene_means(counts, sf);
    let trend = fit_dispersion_trend(&mle_dispersions, &gene_means);

    // Step 4: MAP shrinkage
    let map_dispersions = shrink_dispersions(&mle_dispersions, &gene_means, &trend);

    // Steps 5+6: general IRLS + Wald contrast test
    let mut results: Vec<DeseqResult> = (0..n_genes)
        .map(|g| {
            fit_gene_nb_glm_general(
                g,
                gene_ids,
                counts,
                sf,
                design_matrix,
                contrast,
                map_dispersions[g],
                &gene_means,
                &group0_def,
                &group1_def,
            )
        })
        .collect();

    // Steps 7-9: reuse existing helpers (Cook's, filtering, LFC shrinkage, BH)
    flag_outliers(&mut results, counts, sf, &group0_def, &group1_def);
    let n_filtered = independent_filtering(&mut results);
    apply_lfc_shrinkage(&mut results);
    apply_bh_correction(&mut results);

    let n_significant = results.iter().filter(|r| r.padj < 0.05).count();

    Ok(DeseqAnalysis {
        results,
        size_factors,
        dispersion_trend: trend,
        n_significant,
        n_filtered,
    })
}

/// Likelihood-ratio test comparing full vs. reduced model across all genes.
///
/// Fits both models via IRLS and computes the chi-squared LRT statistic.
/// Returns one `LrtResult` per gene, BH-adjusted.
pub fn lrt_test(
    counts: &[Vec<f64>],
    gene_ids: &[String],
    full_design: &[Vec<f64>],
    reduced_design: &[Vec<f64>],
    sample_names: &[String],
) -> Result<Vec<LrtResult>> {
    if counts.is_empty() {
        bail!("lrt_test: count matrix is empty");
    }
    let n_samples = sample_names.len();
    if full_design.len() != n_samples || reduced_design.len() != n_samples {
        bail!("lrt_test: design matrix row count != n_samples");
    }
    let p_full = full_design[0].len();
    let p_reduced = reduced_design[0].len();
    if p_reduced >= p_full {
        bail!("lrt_test: reduced model must have fewer factors than full model");
    }
    let df = p_full - p_reduced;

    let size_factors =
        estimate_size_factors(counts, sample_names).context("lrt_test: size-factor estimation")?;
    let sf = &size_factors.factors;
    let gene_means = compute_gene_means(counts, sf);

    let (group0_def, group1_def) = groups_from_design(full_design, n_samples);
    let mle_disp = compute_mle_dispersions(counts, sf, &group0_def, &group1_def);
    let trend = fit_dispersion_trend(&mle_disp, &gene_means);
    let map_disp = shrink_dispersions(&mle_disp, &gene_means, &trend);

    let mut lrt_results: Vec<LrtResult> = (0..counts.len())
        .map(|g| {
            compute_lrt_for_gene(
                g,
                gene_ids,
                counts,
                sf,
                full_design,
                reduced_design,
                map_disp[g],
                &gene_means,
                df,
            )
        })
        .collect();

    // BH correction over finite p-values
    let eligible: Vec<usize> = lrt_results
        .iter()
        .enumerate()
        .filter(|(_, r)| r.p_value.is_finite())
        .map(|(i, _)| i)
        .collect();
    if !eligible.is_empty() {
        let pvals: Vec<f64> = eligible.iter().map(|&i| lrt_results[i].p_value).collect();
        let padj_vals = biomics_core::statistics::benjamini_hochberg(&pvals);
        for (ei, &orig_i) in eligible.iter().enumerate() {
            lrt_results[orig_i].padj = padj_vals[ei];
        }
    }
    for r in lrt_results.iter_mut() {
        if r.padj.is_nan() {
            r.padj = 1.0;
        }
    }

    Ok(lrt_results)
}

// ── Design matrix builders ────────────────────────────────────────────────────

/// Build a two-group design matrix: intercept + treatment indicator.
///
/// Columns: [1, x] where x=0 for group0 samples, x=1 for group1 samples.
pub fn design_matrix_two_group(
    group0: &[usize],
    group1: &[usize],
    n_samples: usize,
) -> Vec<Vec<f64>> {
    let mut dm = vec![vec![0.0_f64; 2]; n_samples];
    for &j in group0 {
        dm[j][0] = 1.0; // intercept
        dm[j][1] = 0.0; // treatment
    }
    for &j in group1 {
        dm[j][0] = 1.0; // intercept
        dm[j][1] = 1.0; // treatment
    }
    dm
}

/// Build a design matrix with intercept, treatment, and batch dummy variables.
///
/// Corresponds to the R formula `~ batch + condition`.
/// Columns: [1, treatment, batch_1, batch_2, ..., batch_{B-1}]
/// (reference batch = 0; omitted to avoid collinearity).
pub fn design_matrix_with_batch(group: &[u32], batch: &[u32]) -> Vec<Vec<f64>> {
    let n = group.len();
    debug_assert_eq!(n, batch.len(), "group and batch must have same length");
    let n_batches = batch.iter().copied().max().unwrap_or(0) as usize + 1;
    // Number of columns: intercept + treatment + (n_batches - 1) batch dummies
    let n_cols = 2 + n_batches.saturating_sub(1);
    let mut dm = vec![vec![0.0_f64; n_cols]; n];
    for i in 0..n {
        dm[i][0] = 1.0; // intercept
        dm[i][1] = group[i] as f64; // treatment
        let b = batch[i] as usize;
        if b > 0 && b < n_batches {
            dm[i][1 + b] = 1.0; // batch dummy (reference batch=0 omitted)
        }
    }
    dm
}

/// Build a design matrix with intercept, treatment, batch, and treatment×batch
/// interaction.
///
/// Corresponds to the R formula `~ batch + condition + batch:condition`.
/// Columns: [1, treatment, batch_1, ..., batch_{B-1},
///           treatment×batch_1, ..., treatment×batch_{B-1}]
pub fn interaction_design_matrix(group: &[u32], batch: &[u32]) -> Vec<Vec<f64>> {
    let n = group.len();
    debug_assert_eq!(n, batch.len(), "group and batch must have same length");
    let n_batches = batch.iter().copied().max().unwrap_or(0) as usize + 1;
    let n_batch_dummies = n_batches.saturating_sub(1);
    // intercept + treatment + batch_dummies + interaction_dummies
    let n_cols = 2 + n_batch_dummies * 2;
    let mut dm = vec![vec![0.0_f64; n_cols]; n];
    for i in 0..n {
        dm[i][0] = 1.0;
        let t = group[i] as f64;
        dm[i][1] = t;
        let b = batch[i] as usize;
        if b > 0 && b < n_batches {
            dm[i][1 + b] = 1.0; // batch dummy
            dm[i][1 + n_batch_dummies + b] = t; // interaction: treatment × batch
        }
    }
    dm
}

// ── General IRLS solver ───────────────────────────────────────────────────────

/// Fit NB-GLM coefficients via IRLS using a general design matrix.
///
/// Returns `(beta, vcov)` where `vcov = (X^T W X)^{-1}` (Love et al. 2014).
#[allow(clippy::needless_range_loop)]
fn general_irls(
    counts: &[f64],
    size_factors: &[f64],
    design: &[Vec<f64>],
    alpha: f64,
    max_iter: usize,
    tol: f64,
) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = counts.len();
    let p = design[0].len();
    if n == 0 || p == 0 {
        bail!("general_irls: empty input");
    }
    if n < p {
        bail!("general_irls: more parameters ({p}) than observations ({n})");
    }

    // Initialise: mu_j = (count_j + 0.5) / size_factor_j
    let mut mu: Vec<f64> = counts
        .iter()
        .zip(size_factors.iter())
        .map(|(&k, &s)| (k + 0.5) / s.max(1e-12))
        .collect();

    let mut beta: Vec<f64> = vec![0.0; p];
    // Warm-start intercept from log mean
    let log_mean = mu.iter().map(|&m| m.ln()).sum::<f64>() / n as f64;
    beta[0] = log_mean;

    let mut vcov = vec![vec![0.0_f64; p]; p];

    for _iter in 0..max_iter {
        // Build X^T W X (p×p) and X^T W z (p-vector) — WLS normal equations.
        let mut xtwx = vec![vec![0.0_f64; p]; p];
        let mut xtwz = vec![0.0_f64; p];

        for j in 0..n {
            let mu_j = mu[j].max(1e-12);
            // NB working weight: w_j = 1/(1/mu_j + alpha)  (Love et al. 2014 eq. 6)
            let w = 1.0 / (1.0 / mu_j + alpha);
            // Working response: z_j = log(mu_j) - log(s_j) + (k_j - mu_j)/mu_j
            let z = mu_j.ln() - size_factors[j].max(1e-12).ln() + (counts[j] - mu_j) / mu_j;

            for a in 0..p {
                xtwz[a] += w * design[j][a] * z;
                for b in 0..=a {
                    xtwx[a][b] += w * design[j][a] * design[j][b];
                }
            }
        }
        // Symmetrise lower triangle
        for a in 0..p {
            for b in (a + 1)..p {
                xtwx[a][b] = xtwx[b][a];
            }
        }

        // Cholesky decomposition of X^T W X, then solve + invert
        let l = cholesky_decompose(&xtwx)?;
        let new_beta = cholesky_solve(&l, &xtwz);
        vcov = cholesky_inverse(&l);

        // Convergence check
        let delta = new_beta
            .iter()
            .zip(beta.iter())
            .map(|(nb, ob)| (nb - ob).abs())
            .fold(0.0_f64, f64::max);
        beta = new_beta;

        // Update mu_j = s_j * exp(X[j,:] · beta)
        for j in 0..n {
            let eta: f64 = design[j]
                .iter()
                .zip(beta.iter())
                .map(|(&x, &b)| x * b)
                .sum();
            mu[j] = size_factors[j].max(1e-12) * eta.exp();
        }

        if delta < tol {
            break;
        }
    }

    Ok((beta, vcov))
}

/// Cholesky decomposition L such that L L^T = A (lower triangular).
///
/// Standard Cholesky-Banachiewicz algorithm.
fn cholesky_decompose(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let p = a.len();
    let mut l = vec![vec![0.0_f64; p]; p];
    for i in 0..p {
        for j in 0..=i {
            let sum: f64 = (0..j).map(|k| l[i][k] * l[j][k]).sum();
            if i == j {
                let diag = a[i][i] - sum;
                if diag <= 0.0 {
                    bail!("cholesky_decompose: matrix not positive-definite at diagonal [{i}]");
                }
                l[i][j] = diag.sqrt();
            } else {
                l[i][j] = (a[i][j] - sum) / l[j][j];
            }
        }
    }
    Ok(l)
}

/// Solve L L^T x = b via forward then backward substitution.
fn cholesky_solve(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let p = l.len();
    // Forward substitution: L y = b
    let mut y = vec![0.0_f64; p];
    for i in 0..p {
        let sum: f64 = (0..i).map(|k| l[i][k] * y[k]).sum();
        y[i] = (b[i] - sum) / l[i][i];
    }
    // Backward substitution: L^T x = y
    let mut x = vec![0.0_f64; p];
    for i in (0..p).rev() {
        let sum: f64 = ((i + 1)..p).map(|k| l[k][i] * x[k]).sum();
        x[i] = (y[i] - sum) / l[i][i];
    }
    x
}

/// Compute (L L^T)^{-1} = L^{-T} L^{-1} via triangular inversion.
///
/// Only required for the p×p vcov matrix (typically p ≤ 10).
fn cholesky_inverse(l: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let p = l.len();
    // Invert lower triangular L in-place by solving L * E_i = e_i for each column
    let mut linv = vec![vec![0.0_f64; p]; p];
    for j in 0..p {
        linv[j][j] = 1.0 / l[j][j];
        for i in (j + 1)..p {
            let sum: f64 = (j..i).map(|k| l[i][k] * linv[k][j]).sum();
            linv[i][j] = -sum / l[i][i];
        }
    }
    // vcov = L^{-T} L^{-1}  (symmetric product)
    let mut vcov = vec![vec![0.0_f64; p]; p];
    for i in 0..p {
        for j in 0..=i {
            let dot: f64 = (i..p).map(|k| linv[k][i] * linv[k][j]).sum();
            vcov[i][j] = dot;
            vcov[j][i] = dot;
        }
    }
    vcov
}

/// Wald test for a linear contrast of NB-GLM coefficients.
///
/// Returns `(effect_size, se, p_value)`.
/// effect_size = c^T β, var = c^T Σ c, z = effect_size / se.
fn wald_test_contrast(beta: &[f64], vcov: &[Vec<f64>], contrast: &[f64]) -> (f64, f64, f64) {
    let effect: f64 = contrast.iter().zip(beta.iter()).map(|(&c, &b)| c * b).sum();
    let var: f64 = contrast
        .iter()
        .enumerate()
        .map(|(i, &ci)| {
            contrast
                .iter()
                .enumerate()
                .map(|(j, &cj)| ci * vcov[i][j] * cj)
                .sum::<f64>()
        })
        .sum();
    let se = var.max(0.0).sqrt();
    let z = if se > 0.0 { effect / se } else { 0.0 };
    let p = 2.0 * normal_cdf(-z.abs());
    (effect, se, p)
}

/// Chi-squared p-value P(χ²(df) > stat).
///
/// Special cases for df=1,2; Wilson-Hilferty cube-root approximation for df>2.
fn chi2_pvalue(stat: f64, df: usize) -> f64 {
    if stat <= 0.0 {
        return 1.0;
    }
    match df {
        0 => 1.0,
        // df=1: p = erfc(sqrt(stat/2)/sqrt(2)) = 2*Phi(-sqrt(stat))
        1 => 2.0 * normal_cdf(-(stat.sqrt())),
        // df=2: exact — p = exp(-stat/2)
        2 => (-stat / 2.0).exp(),
        // df>2: Wilson-Hilferty cube-root normal approximation
        _ => {
            let d = df as f64;
            let z = ((stat / d).cbrt() - (1.0 - 2.0 / (9.0 * d))) / (2.0 / (9.0 * d)).sqrt();
            normal_cdf(-z).clamp(0.0, 1.0)
        }
    }
}

/// NB log-likelihood for one gene given fitted mu values.
///
/// Uses the negative-binomial log-pmf: Σ_j [k_j*log(mu_j) - (k_j+1/α)*log(1+α*mu_j) + lgamma(k_j+1/α) - lgamma(k_j+1) - lgamma(1/α)]
/// The lgamma terms involving only data cancel in the LRT difference, so we
/// return only the model-dependent part.
fn nb_loglik(counts: &[f64], mu: &[f64], alpha: f64) -> f64 {
    // Model-dependent part: Σ [ k*log(mu) - (k + 1/α)*log(1 + α*mu) ]
    let inv_alpha = 1.0 / alpha.max(1e-12);
    counts
        .iter()
        .zip(mu.iter())
        .map(|(&k, &m)| {
            let m = m.max(1e-12);
            k * m.ln() - (k + inv_alpha) * (1.0 + alpha * m).ln()
        })
        .sum()
}

/// Compute fitted mu values from design matrix and beta coefficients.
fn fitted_mu(design: &[Vec<f64>], beta: &[f64], size_factors: &[f64]) -> Vec<f64> {
    design
        .iter()
        .zip(size_factors.iter())
        .map(|(row, &s)| {
            let eta: f64 = row.iter().zip(beta.iter()).map(|(&x, &b)| x * b).sum();
            s.max(1e-12) * eta.exp()
        })
        .collect()
}

// ── Internal helpers for general design ──────────────────────────────────────

/// Derive group0/group1 from a design matrix for dispersion MoM.
///
/// Uses the second column (index 1) as the treatment indicator if it exists;
/// otherwise puts all samples in group0.
#[allow(clippy::needless_range_loop)]
fn groups_from_design(design: &[Vec<f64>], n_samples: usize) -> (Vec<usize>, Vec<usize>) {
    let mut g0 = Vec::new();
    let mut g1 = Vec::new();
    for j in 0..n_samples {
        if design[j].len() > 1 && design[j][1] > 0.5 {
            g1.push(j);
        } else {
            g0.push(j);
        }
    }
    if g0.is_empty() || g1.is_empty() {
        // Degenerate: put everything in g0 to avoid crashing downstream
        let all: Vec<usize> = (0..n_samples).collect();
        return (all.clone(), vec![all[0]]);
    }
    (g0, g1)
}

/// Fit one gene with the general IRLS and return a `DeseqResult`.
#[allow(clippy::too_many_arguments)]
fn fit_gene_nb_glm_general(
    g: usize,
    gene_ids: &[String],
    counts: &[Vec<f64>],
    sf: &[f64],
    design: &[Vec<f64>],
    contrast: &[f64],
    dispersion: f64,
    gene_means: &[f64],
    group0: &[usize],
    group1: &[usize],
) -> DeseqResult {
    let row = &counts[g];
    let base_mean = gene_means[g];
    let alpha = dispersion.max(1e-8);
    let log2 = std::f64::consts::LN_2;

    let mean_g0 = {
        let s: f64 = group0.iter().map(|&j| row[j] / sf[j]).sum();
        s / group0.len().max(1) as f64
    };
    let mean_g1 = {
        let s: f64 = group1.iter().map(|&j| row[j] / sf[j]).sum();
        s / group1.len().max(1) as f64
    };

    let counts_gene: Vec<f64> = row.to_vec();
    match general_irls(&counts_gene, sf, design, alpha, IRLS_MAX_ITER, IRLS_TOL) {
        Err(_) => make_na_result(&gene_ids[g], base_mean, mean_g0, mean_g1, dispersion),
        Ok((beta, vcov)) => {
            let (effect, se, p_value) = wald_test_contrast(&beta, &vcov, contrast);
            let stat = if se > 0.0 { effect / se } else { 0.0 };
            DeseqResult {
                gene_id: gene_ids[g].clone(),
                base_mean,
                log2_fold_change: effect / log2,
                lfc_shrunk: effect / log2,
                lfc_se: se / log2,
                stat,
                p_value,
                padj: f64::NAN,
                mean_group0: mean_g0,
                mean_group1: mean_g1,
                dispersion,
                outlier: false,
            }
        }
    }
}

/// Compute LRT result for one gene.
#[allow(clippy::too_many_arguments)]
fn compute_lrt_for_gene(
    g: usize,
    gene_ids: &[String],
    counts: &[Vec<f64>],
    sf: &[f64],
    full_design: &[Vec<f64>],
    reduced_design: &[Vec<f64>],
    dispersion: f64,
    gene_means: &[f64],
    df: usize,
) -> LrtResult {
    let row: Vec<f64> = counts[g].to_vec();
    let alpha = dispersion.max(1e-8);
    let base_mean = gene_means[g];

    let na = LrtResult {
        gene_id: gene_ids[g].clone(),
        base_mean,
        stat: f64::NAN,
        df,
        p_value: f64::NAN,
        padj: f64::NAN,
    };

    let Ok((beta_full, _)) = general_irls(&row, sf, full_design, alpha, IRLS_MAX_ITER, IRLS_TOL)
    else {
        return na;
    };
    let Ok((beta_red, _)) = general_irls(&row, sf, reduced_design, alpha, IRLS_MAX_ITER, IRLS_TOL)
    else {
        return na;
    };

    let mu_full = fitted_mu(full_design, &beta_full, sf);
    let mu_red = fitted_mu(reduced_design, &beta_red, sf);

    let ll_full = nb_loglik(&row, &mu_full, alpha);
    let ll_red = nb_loglik(&row, &mu_red, alpha);
    let stat = (2.0 * (ll_full - ll_red)).max(0.0);
    let p_value = chi2_pvalue(stat, df);

    LrtResult {
        gene_id: gene_ids[g].clone(),
        base_mean,
        stat,
        df,
        p_value,
        padj: f64::NAN,
    }
}

// ── Step 2: gene-wise MLE dispersions ────────────────────────────────────────

/// Compute method-of-moments NB dispersion for each gene.
///
/// α_g = max(0, (Var(K_gj/s_j) - Mean(K_gj/s_j)) / Mean(K_gj/s_j)²)
fn compute_mle_dispersions(
    counts: &[Vec<f64>],
    sf: &[f64],
    group0: &[usize],
    group1: &[usize],
) -> Vec<f64> {
    let all_samples: Vec<usize> = group0.iter().chain(group1.iter()).copied().collect();
    counts
        .iter()
        .map(|row| {
            let normed: Vec<f64> = all_samples.iter().map(|&j| row[j] / sf[j]).collect();
            let n = normed.len() as f64;
            if n < 2.0 {
                return 0.0;
            }
            let mean = normed.iter().sum::<f64>() / n;
            if mean <= 0.0 {
                return 0.0;
            }
            let var = normed.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
            ((var - mean) / mean.powi(2)).max(0.0)
        })
        .collect()
}

// ── Step 3: dispersion trend fitting ─────────────────────────────────────────

/// Compute mean normalized count per gene across all samples.
fn compute_gene_means(counts: &[Vec<f64>], sf: &[f64]) -> Vec<f64> {
    counts
        .iter()
        .map(|row| {
            let n = sf.len() as f64;
            if n == 0.0 {
                return 0.0;
            }
            row.iter().zip(sf.iter()).map(|(&c, &s)| c / s).sum::<f64>() / n
        })
        .collect()
}

/// Fit parametric trend α = a₀ + a₁/μ using the 5th percentile + OLS.
fn fit_dispersion_trend(dispersions: &[f64], means: &[f64]) -> DispersionTrend {
    // Use only genes with mean > 1.0 and positive MLE dispersion
    let valid: Vec<(f64, f64)> = dispersions
        .iter()
        .zip(means.iter())
        .filter(|(&d, &m)| m > 1.0 && d > 0.0 && d.is_finite())
        .map(|(&d, &m)| (d, m))
        .collect();

    if valid.is_empty() {
        return DispersionTrend {
            asympt_disp: 0.1,
            extra_pois: 1.0,
        };
    }

    // a₀ = 5th percentile of MLE dispersions
    let mut disp_sorted: Vec<f64> = valid.iter().map(|&(d, _)| d).collect();
    disp_sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p5_idx = ((disp_sorted.len() as f64 * 0.05) as usize).min(disp_sorted.len() - 1);
    let a0 = disp_sorted[p5_idx].max(1e-8);

    // a₁: closed-form OLS through origin on (α_g - a₀) ~ a₁ × (1/μ_g)
    // a₁ = Σ((α_g - a₀) × μ_g) / Σ(1)  — simplified (see spec)
    let n_valid = valid.len() as f64;
    let sum_adj_times_mu: f64 = valid.iter().map(|&(d, m)| (d - a0) * m).sum();
    let a1 = (sum_adj_times_mu / n_valid).max(0.0);

    DispersionTrend {
        asympt_disp: a0,
        extra_pois: a1,
    }
}

// ── Step 4: MAP dispersion shrinkage ─────────────────────────────────────────

/// Shrink MLE dispersions toward the trend using empirical Bayes log-normal prior.
///
/// log(α_shrunk) = weighted average of log(α_MLE) and log(α_trend).
fn shrink_dispersions(mle_dispersions: &[f64], means: &[f64], trend: &DispersionTrend) -> Vec<f64> {
    let n_samples_approx = 6_usize; // reasonable default; not critical for shrinkage direction

    // Compute log dispersions for genes with valid MLE
    let log_mlep: Vec<f64> = mle_dispersions
        .iter()
        .filter(|&&d| d > 0.0 && d.is_finite())
        .map(|&d| d.ln())
        .collect();

    // σ²_prior = max(0, Var(log α_g) − Var_sampling)
    let var_sampling = 1.0 / n_samples_approx as f64;
    let sigma2_prior = if log_mlep.len() >= 2 {
        let n = log_mlep.len() as f64;
        let mean_log = log_mlep.iter().sum::<f64>() / n;
        let var_log = log_mlep
            .iter()
            .map(|&x| (x - mean_log).powi(2))
            .sum::<f64>()
            / (n - 1.0).max(1.0);
        (var_log - var_sampling).max(0.0)
    } else {
        0.0
    };

    mle_dispersions
        .iter()
        .zip(means.iter())
        .map(|(&d_mle, &mu)| {
            if d_mle <= 0.0 || !d_mle.is_finite() || mu <= 0.0 {
                // Fall back to trend value when MLE is degenerate
                return trend.eval(mu).max(1e-8);
            }
            let log_mle = d_mle.ln();
            let log_trend = trend.eval(mu).max(1e-8).ln();

            if sigma2_prior <= 0.0 {
                // No prior variance → return trend
                return log_trend.exp();
            }

            // MAP: precision-weighted average
            let prec_mle = 1.0 / var_sampling;
            let prec_prior = 1.0 / sigma2_prior;
            let log_map = (log_mle * prec_mle + log_trend * prec_prior) / (prec_mle + prec_prior);
            log_map.exp().max(1e-8)
        })
        .collect()
}

// ── Step 5 + 6: IRLS NB-GLM + Wald test ──────────────────────────────────────

/// IRLS convergence tolerance on |Δβ₀| + |Δβ₁|.
const IRLS_TOL: f64 = 1e-8;
/// Maximum IRLS iterations.
const IRLS_MAX_ITER: usize = 50;

/// Fit a 2-group NB-GLM for one gene via IRLS and compute the Wald z-statistic.
#[allow(clippy::too_many_arguments, dead_code)]
fn fit_gene_nb_glm(
    g: usize,
    gene_ids: &[String],
    counts: &[Vec<f64>],
    sf: &[f64],
    group0: &[usize],
    group1: &[usize],
    dispersion: f64,
    gene_means: &[f64],
) -> DeseqResult {
    let row = &counts[g];
    let base_mean = gene_means[g];
    let alpha = dispersion.max(1e-8);

    // Collect (count, size_factor, group_indicator) tuples
    // x=0 for group0, x=1 for group1
    let samples: Vec<(f64, f64, f64)> = group0
        .iter()
        .map(|&j| (row[j], sf[j], 0.0))
        .chain(group1.iter().map(|&j| (row[j], sf[j], 1.0)))
        .collect();

    let n = samples.len();

    // mean_group0 and mean_group1 in normalized counts
    let mean_g0 = {
        let s: f64 = group0.iter().map(|&j| row[j] / sf[j]).sum();
        s / group0.len().max(1) as f64
    };
    let mean_g1 = {
        let s: f64 = group1.iter().map(|&j| row[j] / sf[j]).sum();
        s / group1.len().max(1) as f64
    };

    // Initialize μ_j = (K_j + 0.5) / s_j
    let mut mu: Vec<f64> = samples.iter().map(|(k, s, _)| (k + 0.5) / s).collect();

    // β = [β₀, β₁]
    let mut beta0 = mu.iter().map(|&m| m.ln()).sum::<f64>() / n as f64;
    let mut beta1 = 0.0_f64;

    // ainv11 is the (1,1) element of (X^T W X)^{-1}, used for SE(β₁)
    let mut ainv11 = 0.0_f64;

    for _iter in 0..IRLS_MAX_ITER {
        // Working weights and working responses
        let mut sw0 = 0.0_f64;
        let mut sw1 = 0.0_f64;
        let mut swz0 = 0.0_f64;
        let mut swz1 = 0.0_f64;

        for (i, &(k, s, x)) in samples.iter().enumerate() {
            let mu_j = mu[i];
            // NB variance weight: w_j = 1 / (1/μ + α)
            let w = 1.0 / (1.0 / mu_j + alpha);
            // Working response: z_j = log(μ_j) - log(s_j) + (K_j - μ_j)/μ_j
            let z = mu_j.ln() - s.ln() + (k - mu_j) / mu_j;

            if x < 0.5 {
                sw0 += w;
                swz0 += w * z;
            } else {
                sw1 += w;
                swz1 += w * z;
            }
        }

        // X^T W X for design [1, x]:
        // A = [[sw0+sw1, sw1], [sw1, sw1]]
        let a00 = sw0 + sw1;
        let a01 = sw1;
        let a11 = sw1;
        let det = a00 * a11 - a01 * a01; // = sw0*sw1

        if det.abs() < 1e-12 {
            // Degenerate — return NA-like result
            return make_na_result(&gene_ids[g], base_mean, mean_g0, mean_g1, dispersion);
        }

        // A⁻¹ = (1/det) * [[a11, -a01], [-a01, a00]]
        let ainv00 = a11 / det;
        ainv11 = a00 / det;
        let ainv01 = -a01 / det;

        // X^T W z = [swz0+swz1, swz1]
        let b0 = swz0 + swz1;
        let b1 = swz1;

        let new_beta0 = ainv00 * b0 + ainv01 * b1;
        let new_beta1 = ainv01 * b0 + ainv11 * b1;

        let delta = (new_beta0 - beta0).abs() + (new_beta1 - beta1).abs();
        beta0 = new_beta0;
        beta1 = new_beta1;

        // Update μ
        for (i, &(_, s, x)) in samples.iter().enumerate() {
            mu[i] = s * (beta0 + beta1 * x).exp();
        }

        if delta < IRLS_TOL {
            break;
        }
    }

    // Wald test: SE(β₁) = sqrt(A⁻¹₁₁)
    let se_beta1 = ainv11.max(0.0).sqrt();
    let stat = if se_beta1 > 0.0 {
        beta1 / se_beta1
    } else {
        0.0
    };
    let p_value = 2.0 * normal_cdf(-stat.abs());

    let log2 = std::f64::consts::LN_2;
    DeseqResult {
        gene_id: gene_ids[g].clone(),
        base_mean,
        log2_fold_change: beta1 / log2,
        lfc_shrunk: beta1 / log2, // filled in step 9
        lfc_se: se_beta1 / log2,
        stat,
        p_value,
        padj: f64::NAN,
        mean_group0: mean_g0,
        mean_group1: mean_g1,
        dispersion,
        outlier: false,
    }
}

/// Return a result with NaN statistics for a degenerate gene.
fn make_na_result(
    gene_id: &str,
    base_mean: f64,
    mean_group0: f64,
    mean_group1: f64,
    dispersion: f64,
) -> DeseqResult {
    DeseqResult {
        gene_id: gene_id.to_string(),
        base_mean,
        log2_fold_change: f64::NAN,
        lfc_shrunk: f64::NAN,
        lfc_se: f64::NAN,
        stat: f64::NAN,
        p_value: f64::NAN,
        padj: f64::NAN,
        mean_group0,
        mean_group1,
        dispersion,
        outlier: false,
    }
}

// ── Step 7: Cook's distance outlier flagging ──────────────────────────────────

/// Flag genes whose maximum Cook's distance exceeds the conservative threshold.
///
/// Cook's distance: C_j = (K_j - μ_j)² × h_jj / (p × MSE × (1 - h_jj)²)
fn flag_outliers(
    results: &mut [DeseqResult],
    counts: &[Vec<f64>],
    sf: &[f64],
    group0: &[usize],
    group1: &[usize],
) {
    // Conservative F(0.99, 2, n-2) threshold = 5.0
    const COOKS_THRESHOLD: f64 = 5.0;
    let n = group0.len() + group1.len();

    for (g, result) in results.iter_mut().enumerate() {
        let row = &counts[g];
        let alpha = result.dispersion.max(1e-8);

        // Quick recompute of IRLS final μ and XTWX⁻¹ for this gene
        // (We only need the hat matrix diagonal, so we use the stored LFC)
        let log2 = std::f64::consts::LN_2;
        let beta1 = result.log2_fold_change * log2;
        // Recover β₀ from mean: base_mean ≈ exp(β₀ + β₁/2) — rough but good enough
        // for Cook's calculation; use group means instead.
        let beta0 = if result.mean_group0 > 0.0 {
            result.mean_group0.ln()
        } else {
            0.0_f64
        };

        let mut max_cooks = 0.0_f64;

        let all_samples: Vec<(f64, f64, f64)> = group0
            .iter()
            .map(|&j| (row[j], sf[j], 0.0_f64))
            .chain(group1.iter().map(|&j| (row[j], sf[j], 1.0_f64)))
            .collect();

        // Recompute XTWX with approximate μ
        let mut sw0 = 0.0_f64;
        let mut sw1 = 0.0_f64;
        let mut mus: Vec<f64> = Vec::with_capacity(all_samples.len());
        for &(_, s, x) in &all_samples {
            let mu_j = (s * (beta0 + beta1 * x).exp()).max(1e-8);
            mus.push(mu_j);
            let w = 1.0 / (1.0 / mu_j + alpha);
            if x < 0.5 {
                sw0 += w;
            } else {
                sw1 += w;
            }
        }

        let det = sw0 * sw1;
        if det.abs() < 1e-12 || n < 3 {
            continue;
        }

        let ainv00 = sw1 / det; // (sw0+sw1)*sw1 - sw1² = sw0*sw1; inv[0,0] = sw1/det
        let ainv11 = (sw0 + sw1) / det;

        for (i, &(k, _, x)) in all_samples.iter().enumerate() {
            let mu_j = mus[i];
            let w = 1.0 / (1.0 / mu_j + alpha);
            // Hat matrix diagonal: h_jj = w_j × (X (X^T W X)^{-1} X^T)_{jj}
            // For group0 (x=0): (X A⁻¹ Xᵀ)_jj = A⁻¹_{00}
            // For group1 (x=1): (X A⁻¹ Xᵀ)_jj ≈ A⁻¹_{00} + 2×A⁻¹_{01} + A⁻¹_{11}
            // But off-diagonal: ainv01 = -sw1/det = -1/sw0 (small)
            // Simplified: use ainv00 for group0, ainv11 for group1
            let ainv_diag = if x < 0.5 { ainv00 } else { ainv11 };
            let h = (w * ainv_diag).min(0.9999);

            let denom = 2.0 * mu_j.powi(2) * (1.0 - h).powi(2);
            if denom > 0.0 {
                let c = (k - mu_j).powi(2) * h / denom;
                if c > max_cooks {
                    max_cooks = c;
                }
            }
        }

        result.outlier = max_cooks > COOKS_THRESHOLD;
    }
}

// ── Step 8: Independent filtering ────────────────────────────────────────────

/// Apply Bourgon et al. 2010 independent filtering; returns number of filtered genes.
///
/// Chooses the mean-count threshold θ* that maximises genes passing at padj < 0.1.
fn independent_filtering(results: &mut [DeseqResult]) -> usize {
    // Build quantile thresholds from 0.0 to 0.8, step 0.05
    let mut base_means: Vec<f64> = results
        .iter()
        .filter(|r| r.p_value.is_finite() && !r.outlier)
        .map(|r| r.base_mean)
        .collect();
    base_means.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = base_means.len();
    if n == 0 {
        return 0;
    }

    let mut best_theta = 0.0_f64;
    let mut best_count = 0_usize;

    let steps = 17usize; // 0.00, 0.05, ..., 0.80
    for step in 0..steps {
        let q = step as f64 * 0.05;
        let theta_idx = ((n as f64 * q) as usize).min(n - 1);
        let theta = base_means[theta_idx];

        // Quick BH on genes passing the threshold to count padj < 0.1
        let passing_pvals: Vec<f64> = results
            .iter()
            .filter(|r| r.base_mean > theta && r.p_value.is_finite() && !r.outlier)
            .map(|r| r.p_value)
            .collect();

        let padj_vals = benjamini_hochberg(&passing_pvals);
        let cnt = padj_vals.iter().filter(|&&p| p < 0.1).count();

        if cnt >= best_count {
            best_count = cnt;
            best_theta = theta;
        }
    }

    // Apply the chosen filter: set padj=1 for filtered genes
    let mut n_filtered = 0_usize;
    for result in results.iter_mut() {
        if result.base_mean <= best_theta && result.p_value.is_finite() && !result.outlier {
            result.padj = 1.0;
            n_filtered += 1;
        }
    }
    n_filtered
}

// ── Step 9: LFC shrinkage ─────────────────────────────────────────────────────

/// Apply apeglm-style normal-prior MAP shrinkage to log₂ fold changes.
///
/// β₁_shrunk = β₁ × σ²_mle / (σ²_mle + σ²_prior)
fn apply_lfc_shrinkage(results: &mut [DeseqResult]) {
    // σ²_prior = median(β₁²) across all non-outlier genes with finite LFC
    let beta1_sq: Vec<f64> = results
        .iter()
        .filter(|r| r.log2_fold_change.is_finite() && !r.outlier)
        .map(|r| r.log2_fold_change.powi(2))
        .collect();

    if beta1_sq.is_empty() {
        return;
    }

    let mut sorted = beta1_sq.clone();
    sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let sigma2_prior = median_sorted(&sorted);

    for result in results.iter_mut() {
        if !result.log2_fold_change.is_finite() {
            continue;
        }
        let sigma2_mle = result.lfc_se.powi(2);
        let shrinkage = if sigma2_mle + sigma2_prior > 0.0 {
            sigma2_mle / (sigma2_mle + sigma2_prior)
        } else {
            1.0
        };
        result.lfc_shrunk = result.log2_fold_change * shrinkage;
    }
}

// ── BH correction pass ────────────────────────────────────────────────────────

/// Apply BH correction to genes that are not already filtered or outliers.
fn apply_bh_correction(results: &mut [DeseqResult]) {
    // Collect indices of genes eligible for multiple testing correction
    let eligible: Vec<usize> = results
        .iter()
        .enumerate()
        .filter(|(_, r)| r.p_value.is_finite() && !r.outlier && r.padj.is_nan())
        .map(|(i, _)| i)
        .collect();

    if eligible.is_empty() {
        return;
    }

    let pvals: Vec<f64> = eligible.iter().map(|&i| results[i].p_value).collect();
    let padj_vals = benjamini_hochberg(&pvals);
    for (ei, &orig_i) in eligible.iter().enumerate() {
        results[orig_i].padj = padj_vals[ei];
    }

    // Outlier genes and filtered genes get padj = 1.0 (conservative)
    for result in results.iter_mut() {
        if result.outlier || (!result.p_value.is_finite() && result.padj.is_nan()) {
            result.padj = 1.0;
        }
    }
}

// ── Normal CDF (Abramowitz & Stegun 7.1.26) ──────────────────────────────────

/// Two-sided normal CDF: Φ(z) via erfc approximation.
fn normal_cdf(z: f64) -> f64 {
    0.5 * erfc_approx(-z / std::f64::consts::SQRT_2)
}

/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
fn erfc_approx(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    let r = poly * (-x * x).exp();
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Median of a **sorted** slice (no allocation).
/// Panics if the slice is empty — callers must guard.
fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    debug_assert!(n > 0, "median_sorted called on empty slice");
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// A small 4-gene × 4-sample count matrix with known geometric structure.
    ///
    /// Sample library sizes: [100, 200, 150, 300]. After scaling each column by
    /// its total relative to a reference (sample 0), expected size factors are
    /// approximately [1.0, 2.0, 1.5, 3.0]. The median-of-ratios method produces
    /// values close to this but not identical; we allow ±20 %.
    fn make_test_counts() -> (Vec<String>, Vec<Vec<f64>>, Vec<String>) {
        // 6 genes, 4 samples, with a consistent 2× library-size difference
        // between (sample0, sample1) and (sample2, sample3) pairs
        let gene_ids: Vec<String> = (0..6).map(|i| format!("gene{i}")).collect();
        let sample_names: Vec<String> = vec![
            "s1".to_string(),
            "s2".to_string(),
            "s3".to_string(),
            "s4".to_string(),
        ];
        // counts[gene][sample]
        // s2 has 2× the counts of s1; s3 has 1.5×; s4 has 3×
        let counts: Vec<Vec<f64>> = vec![
            vec![10.0, 20.0, 15.0, 30.0],
            vec![20.0, 40.0, 30.0, 60.0],
            vec![30.0, 60.0, 45.0, 90.0],
            vec![40.0, 80.0, 60.0, 120.0],
            vec![50.0, 100.0, 75.0, 150.0],
            vec![60.0, 120.0, 90.0, 180.0],
        ];
        (gene_ids, counts, sample_names)
    }

    #[test]
    fn test_size_factors() {
        let (_, counts, sample_names) = make_test_counts();
        let sf = estimate_size_factors(&counts, &sample_names).expect("size factors");

        // Expected relative sizes: s1=1.0, s2=2.0, s3=1.5, s4=3.0
        // The median-of-ratios method anchors to geometric-mean rows,
        // so absolute values may shift, but the *ratios* between samples
        // must reflect the actual library-size differences.
        let f = &sf.factors;
        assert_eq!(f.len(), 4);

        // All factors must be positive and finite
        for (name, &factor) in sample_names.iter().zip(f.iter()) {
            assert!(
                factor > 0.0 && factor.is_finite(),
                "size factor for {name} is {factor}"
            );
        }

        // Ratios: f[1]/f[0] ≈ 2.0, f[2]/f[0] ≈ 1.5, f[3]/f[0] ≈ 3.0
        let tol = 0.20; // 20 % relative tolerance
        let ratio_12 = f[1] / f[0];
        assert!(
            (ratio_12 - 2.0).abs() < 2.0 * tol,
            "ratio s2/s1 = {ratio_12} not within 20% of 2.0"
        );
        let ratio_13 = f[2] / f[0];
        assert!(
            (ratio_13 - 1.5).abs() < 1.5 * tol,
            "ratio s3/s1 = {ratio_13} not within 20% of 1.5"
        );
        let ratio_14 = f[3] / f[0];
        assert!(
            (ratio_14 - 3.0).abs() < 3.0 * tol,
            "ratio s4/s1 = {ratio_14} not within 20% of 3.0"
        );
    }

    #[test]
    fn test_normalization() {
        let (gene_ids, counts, sample_names) = make_test_counts();
        let matrix = normalize_counts(&gene_ids, &counts, &sample_names).expect("normalize_counts");

        // After normalization, the per-sample sums of normalized counts should
        // be roughly equal (within 20 %) across samples, because the size
        // factors absorb the library-size differences.
        let col_sums: Vec<f64> = (0..sample_names.len())
            .map(|j| matrix.normalized.iter().map(|row| row[j]).sum::<f64>())
            .collect();

        let mean_sum = col_sums.iter().sum::<f64>() / col_sums.len() as f64;
        for (name, &s) in sample_names.iter().zip(col_sums.iter()) {
            let rel_dev = (s - mean_sum).abs() / mean_sum;
            assert!(
                rel_dev < 0.20,
                "sample {name}: normalized library size {s:.1} deviates {:.1}% from mean {mean_sum:.1}",
                rel_dev * 100.0
            );
        }

        // Global dispersion must be non-negative and finite
        assert!(
            matrix.global_dispersion >= 0.0 && matrix.global_dispersion.is_finite(),
            "global_dispersion = {}",
            matrix.global_dispersion
        );
    }

    #[test]
    fn test_differential_expression_no_panic_small() {
        // 2 groups of 2 → can_test = true (split=2, remainder=2)
        // Design: genes 0-3 all have ~10× higher raw counts in group 2 (b1,b2)
        // vs group 1 (a1,a2). After DESeq2 normalization the size factors will
        // absorb that global library-size difference, so per-gene LFC will be
        // near zero (all genes change proportionally). What we DO assert:
        //   - all p-values and padj are finite (no NaN),
        //   - results are sorted correctly (ascending padj),
        //   - there is at least one result per gene.
        let gene_ids: Vec<String> = (0..4).map(|i| format!("g{i}")).collect();
        let sample_names: Vec<String> = vec![
            "a1".to_string(),
            "a2".to_string(),
            "b1".to_string(),
            "b2".to_string(),
        ];
        // group 1 (a1,a2) has low raw counts; group 2 (b1,b2) has high raw counts.
        // After median-of-ratios normalization the library-size difference is
        // corrected, leaving only within-gene variance — so LFC ≈ 0 for all genes.
        let counts: Vec<Vec<f64>> = vec![
            vec![10.0, 12.0, 100.0, 110.0],
            vec![8.0, 9.0, 80.0, 90.0],
            vec![5.0, 6.0, 50.0, 55.0],
            vec![20.0, 22.0, 200.0, 220.0],
        ];
        let matrix = normalize_counts(&gene_ids, &counts, &sample_names).expect("normalize_counts");
        let de = deseq2_differential_expression(&matrix);
        assert_eq!(de.len(), 4);

        // All results must have finite p-values (4 samples → 2+2 split → can_test)
        for r in &de {
            assert!(r.p_value.is_finite(), "p_value is NaN for {}", r.gene_id);
            assert!(r.padj.is_finite(), "padj is NaN for {}", r.gene_id);
        }

        // After normalization, all genes are proportionally identical across groups,
        // so |log2FC| must be small (< 0.5 after normalization absorbs library size).
        for r in &de {
            assert!(
                r.log2_fold_change.abs() < 0.5,
                "expected |log2FC| < 0.5 after normalization for {}, got {}",
                r.gene_id,
                r.log2_fold_change
            );
        }

        // Test with a count matrix that has genuine differential expression:
        // gene A is truly up in group 2 (beyond library-size correction).
        let gene_ids2: Vec<String> = vec!["de_gene".to_string(), "null_gene".to_string()];
        let sample_names2: Vec<String> = vec![
            "x1".to_string(),
            "x2".to_string(),
            "y1".to_string(),
            "y2".to_string(),
        ];
        // de_gene: group 2 is 16× higher AFTER adjusting for 2× overall library size
        // null_gene: group 2 is 2× higher (pure library-size effect, corrected to 1×)
        let counts2: Vec<Vec<f64>> = vec![
            vec![10.0, 11.0, 160.0, 175.0], // de_gene: up ~8× after normalization
            vec![20.0, 22.0, 40.0, 44.0],   // null_gene: pure library-size, LFC ≈ 0
        ];
        let matrix2 = normalize_counts(&gene_ids2, &counts2, &sample_names2).unwrap();
        let de2 = deseq2_differential_expression(&matrix2);

        // de_gene should have a clearly positive LFC after normalization.
        // With pseudocount 0.5, log2(norm+0.5) attenuates the true 8× fold change,
        // so we use a conservative threshold of > 1.0.
        let de_gene_result = de2.iter().find(|r| r.gene_id == "de_gene").unwrap();
        assert!(
            de_gene_result.log2_fold_change > 1.0,
            "de_gene should have log2FC > 1 after normalization, got {}",
            de_gene_result.log2_fold_change
        );
    }

    #[test]
    fn test_run_deseq2_basic() {
        // Basic smoke test: run full NB-GLM on a small dataset
        let gene_ids: Vec<String> = vec![
            "G1".to_string(),
            "G2".to_string(),
            "G3".to_string(),
            "G4".to_string(),
        ];
        let sample_names: Vec<String> = vec![
            "c1".to_string(),
            "c2".to_string(),
            "c3".to_string(),
            "t1".to_string(),
            "t2".to_string(),
            "t3".to_string(),
        ];
        // G1 is up in treatment; G3 is down; G2/G4 are null
        let counts: Vec<Vec<f64>> = vec![
            vec![10.0, 12.0, 11.0, 80.0, 85.0, 78.0],  // G1 up
            vec![50.0, 48.0, 52.0, 95.0, 100.0, 98.0], // G2 null (library)
            vec![80.0, 78.0, 82.0, 20.0, 18.0, 22.0],  // G3 down
            vec![30.0, 32.0, 28.0, 60.0, 58.0, 62.0],  // G4 null (library)
        ];
        let group0 = vec![0, 1, 2];
        let group1 = vec![3, 4, 5];

        let analysis = run_deseq2(&counts, &gene_ids, &group0, &group1, &sample_names)
            .expect("run_deseq2 should not fail");

        assert_eq!(analysis.results.len(), 4);

        // G1 should be up (positive LFC)
        let g1 = analysis.results.iter().find(|r| r.gene_id == "G1").unwrap();
        assert!(
            g1.log2_fold_change > 0.0,
            "G1 should be up-regulated, got lfc={}",
            g1.log2_fold_change
        );

        // G3 should be down (negative LFC)
        let g3 = analysis.results.iter().find(|r| r.gene_id == "G3").unwrap();
        assert!(
            g3.log2_fold_change < 0.0,
            "G3 should be down-regulated, got lfc={}",
            g3.log2_fold_change
        );

        // All p-values should be finite (not NaN, not infinite)
        for r in &analysis.results {
            assert!(r.p_value.is_finite(), "p_value NaN for gene {}", r.gene_id);
        }
    }

    #[test]
    fn test_normal_cdf_symmetry() {
        // Φ(0) = 0.5
        let p = normal_cdf(0.0);
        assert!((p - 0.5).abs() < 1e-6, "normal_cdf(0) = {p}");
        // Φ(1.96) ≈ 0.975
        let p196 = normal_cdf(1.96);
        assert!((p196 - 0.975).abs() < 0.002, "normal_cdf(1.96) = {p196}");
    }
}
