use ndarray::{Array1, Array2};

use biomics_core::statistics::spearman_r;
use epigenomics_core::EpigenomicsSummary;
use genomics_core::GenomicsSummary;
use transcriptomics_core::TranscriptomicsSummary;

/// Build a feature vector from genomics data.
///
/// Returns per-chromosome variant counts, normalized to variants per thousand sites,
/// padded or trimmed to a fixed length for correlation computation.
pub fn genomics_feature_vec(summary: &GenomicsSummary, chrom_order: &[String]) -> Array1<f64> {
    let mut v: Vec<f64> = chrom_order
        .iter()
        .map(|c| {
            summary
                .per_chrom
                .get(c)
                .map(|d| d.total as f64)
                .unwrap_or(0.0)
        })
        .collect();

    // Normalize by total
    let total: f64 = v.iter().sum();
    if total > 0.0 {
        v.iter_mut().for_each(|x| *x /= total);
    }

    Array1::from_vec(v)
}

/// Build a feature vector from transcriptomics top-N gene expression.
pub fn transcriptomics_feature_vec(summary: &TranscriptomicsSummary, top_n: usize) -> Array1<f64> {
    let mut v: Vec<f64> = summary
        .top_100_expressed
        .iter()
        .take(top_n)
        .map(|(_, tpm)| *tpm)
        .collect();

    // Pad to top_n if fewer genes available
    while v.len() < top_n {
        v.push(0.0);
    }

    // Normalize
    let total: f64 = v.iter().sum();
    if total > 0.0 {
        v.iter_mut().for_each(|x| *x /= total);
    }

    Array1::from_vec(v)
}

/// Build a feature vector from epigenomics per-chromosome mean methylation.
pub fn epigenomics_feature_vec(
    summary: &EpigenomicsSummary,
    chrom_order: &[String],
) -> Array1<f64> {
    let v: Vec<f64> = chrom_order
        .iter()
        .map(|c| {
            summary
                .per_chrom
                .get(c)
                .map(|cm| cm.mean_methylation)
                .unwrap_or(0.0)
        })
        .collect();

    // Methylation is already in [0, 100]; normalize to [0, 1]
    let v: Vec<f64> = v.iter().map(|&x| x / 100.0).collect();
    Array1::from_vec(v)
}

/// Compute Pearson correlation between two 1D arrays of equal length.
///
/// Uses a single-pass numerically stable algorithm (Chan et al. 1979):
/// accumulates Σ(x-x̄)(y-ȳ), Σ(x-x̄)², Σ(y-ȳ)² in one loop without
/// pre-computing the means first.
pub fn pearson_r(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    // Welford-style single-pass: track means and cross-products simultaneously.
    let mut mean_a = 0.0f64;
    let mut mean_b = 0.0f64;
    let mut cov = 0.0f64;
    let mut var_a = 0.0f64;
    let mut var_b = 0.0f64;

    for (k, (&x, &y)) in a.iter().zip(b.iter()).enumerate() {
        let k1 = (k + 1) as f64;
        let dx = x - mean_a;
        let dy = y - mean_b;
        mean_a += dx / k1;
        mean_b += dy / k1;
        cov += dx * (y - mean_b);
        var_a += dx * (x - mean_a);
        var_b += dy * (y - mean_b);
    }

    let denom = var_a.sqrt() * var_b.sqrt();
    if denom < 1e-12 {
        0.0
    } else {
        (cov / denom).clamp(-1.0, 1.0)
    }
}

/// Compute an N×N Pearson correlation matrix for rows of `data`.
///
/// Each row represents one observation vector. Returns an `Array2<f64>` of
/// shape `(n_rows, n_rows)` where `result[[i, j]]` is the Pearson r between
/// row i and row j.
pub fn pearson_correlation_matrix(data: &Array2<f64>) -> anyhow::Result<Array2<f64>> {
    let n = data.nrows();
    let mut corr = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            let row_i = data.row(i).to_owned();
            let row_j = data.row(j).to_owned();
            corr[[i, j]] = pearson_r(&row_i, &row_j);
        }
    }
    Ok(corr)
}

/// Compute an N×N Spearman rank correlation matrix for rows of `data`.
///
/// Returns `Array2<f64>` of shape `(n_rows, n_rows)` where `result[[i, j]]`
/// is the Spearman r between row i and row j. More robust than Pearson for
/// non-Gaussian or monotone-but-nonlinear relationships.
pub fn spearman_correlation_matrix(data: &Array2<f64>) -> anyhow::Result<Array2<f64>> {
    let n = data.nrows();
    let mut corr = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            let row_i: Vec<f64> = data.row(i).iter().copied().collect();
            let row_j: Vec<f64> = data.row(j).iter().copied().collect();
            corr[[i, j]] = spearman_r(&row_i, &row_j);
        }
    }
    Ok(corr)
}

/// Build the 3×3 cross-modality feature matrix with rows:
/// [genomics_features, transcriptomics_features, epigenomics_features]
/// all aligned to the same feature dimension.
pub fn build_cross_modality_matrix(
    genomics: &GenomicsSummary,
    transcr: &TranscriptomicsSummary,
    epigen: &EpigenomicsSummary,
) -> Array2<f64> {
    // Use chromosomes present in any modality as the shared feature space
    let mut chroms: Vec<String> = genomics
        .per_chrom
        .keys()
        .chain(epigen.per_chrom.keys())
        .cloned()
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();
    chroms.sort();

    let dim = chroms.len().max(1);

    let gvec = genomics_feature_vec(genomics, &chroms);
    let evec = epigenomics_feature_vec(epigen, &chroms);

    // Transcriptomics uses top-N genes, resized to match `dim`
    let mut tvec_raw = transcriptomics_feature_vec(transcr, dim);
    tvec_raw.slice_mut(ndarray::s![..dim.min(tvec_raw.len())]);

    // Ensure all vectors are the same length
    let len = gvec.len().min(evec.len()).min(tvec_raw.len());
    let gvec = gvec.slice(ndarray::s![..len]).to_owned();
    let evec = evec.slice(ndarray::s![..len]).to_owned();
    let tvec = tvec_raw.slice(ndarray::s![..len]).to_owned();

    let mut matrix = Array2::<f64>::zeros((3, len));
    matrix.row_mut(0).assign(&gvec);
    matrix.row_mut(1).assign(&tvec);
    matrix.row_mut(2).assign(&evec);
    matrix
}
