//! Differential expression via Wilcoxon rank-sum test.
//!
//! For each cluster, tests every gene against the remaining cells using the
//! Wilcoxon-Mann-Whitney U statistic with normal approximation for large
//! samples. Adjusted p-values use Benjamini-Hochberg FDR correction.

use ndarray::Array2;

/// Marker gene result for one cluster × gene pair.
#[derive(Debug, Clone, serde::Serialize)]
pub struct ClusterMarker {
    /// Gene name.
    pub gene_id: String,
    /// Cluster index.
    pub cluster: u32,
    /// log2 fold-change (cluster mean / rest mean, pseudocount 1).
    pub log2_fold_change: f64,
    /// Mean normalised expression in the cluster.
    pub mean_expr_cluster: f32,
    /// Mean normalised expression in all other cells.
    pub mean_expr_rest: f32,
    /// Two-sided Wilcoxon p-value.
    pub p_value: f64,
    /// Benjamini-Hochberg adjusted p-value (within cluster).
    pub padj: f64,
    /// Area under ROC curve: U / (n1 × n2).
    pub auc: f64,
}

/// Find marker genes for each cluster using the Wilcoxon rank-sum test.
///
/// Only genes expressed in at least `min_pct` fraction of cells in the
/// cluster are tested, reducing the multiple-testing burden.
///
/// `norm_matrix`: `[n_cells × n_genes]`.
pub fn find_cluster_markers(
    norm_matrix: &Array2<f32>,
    clusters: &[u32],
    feature_names: &[String],
    min_pct: f32,
) -> Vec<ClusterMarker> {
    let (n_cells, n_genes) = norm_matrix.dim();
    if n_cells == 0 || n_genes == 0 {
        return Vec::new();
    }

    let cluster_ids: Vec<u32> = {
        let mut ids: Vec<u32> = clusters.to_vec();
        ids.sort_unstable();
        ids.dedup();
        ids
    };

    let mut all_markers: Vec<ClusterMarker> = Vec::new();

    for &cid in &cluster_ids {
        let in_cluster: Vec<usize> = (0..n_cells).filter(|&i| clusters[i] == cid).collect();
        let in_rest: Vec<usize> = (0..n_cells).filter(|&i| clusters[i] != cid).collect();

        if in_cluster.is_empty() || in_rest.is_empty() {
            continue;
        }

        let n1 = in_cluster.len();
        let n2 = in_rest.len();

        let mut markers: Vec<ClusterMarker> = Vec::new();

        for g in 0..n_genes {
            // Expression vectors
            let group1: Vec<f32> = in_cluster.iter().map(|&i| norm_matrix[[i, g]]).collect();
            let group2: Vec<f32> = in_rest.iter().map(|&i| norm_matrix[[i, g]]).collect();

            // Minimum percent expressed filter
            let pct_expr = group1.iter().filter(|&&v| v > 0.0).count() as f32 / n1 as f32;
            if pct_expr < min_pct {
                continue;
            }

            let mean_cluster = group1.iter().sum::<f32>() / n1 as f32;
            let mean_rest = group2.iter().sum::<f32>() / n2 as f32;

            let lfc = ((mean_cluster as f64 + 1.0) / (mean_rest as f64 + 1.0)).log2();

            let (u_stat, p_value) = wilcoxon_ranksum(&group1, &group2);
            let auc = u_stat / (n1 as f64 * n2 as f64);

            let gene_name = feature_names
                .get(g)
                .cloned()
                .unwrap_or_else(|| format!("gene_{g}"));

            markers.push(ClusterMarker {
                gene_id: gene_name,
                cluster: cid,
                log2_fold_change: lfc,
                mean_expr_cluster: mean_cluster,
                mean_expr_rest: mean_rest,
                p_value,
                padj: p_value, // will be overwritten by BH below
                auc,
            });
        }

        // Benjamini-Hochberg FDR correction within cluster
        bh_correct(&mut markers);
        all_markers.extend(markers);
    }

    all_markers
}

/// Wilcoxon rank-sum test (two-sided, normal approximation for n1+n2 > 20).
///
/// Returns `(U_statistic, p_value)`.
fn wilcoxon_ranksum(group1: &[f32], group2: &[f32]) -> (f64, f64) {
    let n1 = group1.len();
    let n2 = group2.len();
    if n1 == 0 || n2 == 0 {
        return (0.0, 1.0);
    }

    // Combine and rank
    let mut combined: Vec<(f32, usize)> = group1
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i))
        .chain(group2.iter().enumerate().map(|(i, &v)| (v, n1 + i)))
        .collect();
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Average ranks for ties
    let n_total = n1 + n2;
    let mut ranks = vec![0.0f64; n_total];
    let mut i = 0;
    while i < n_total {
        let val = combined[i].0;
        let mut j = i + 1;
        while j < n_total && combined[j].0 == val {
            j += 1;
        }
        let avg_rank = (i + j + 1) as f64 / 2.0; // 1-indexed average
        for k in i..j {
            ranks[combined[k].1] = avg_rank;
        }
        i = j;
    }

    let r1: f64 = ranks[..n1].iter().sum();
    let u1 = r1 - n1 as f64 * (n1 as f64 + 1.0) / 2.0;

    let p = if n1 + n2 > 20 {
        let mean_u = n1 as f64 * n2 as f64 / 2.0;
        let var_u = n1 as f64 * n2 as f64 * (n1 + n2 + 1) as f64 / 12.0;
        let z = (u1 - mean_u) / var_u.sqrt();
        2.0 * normal_sf(z.abs())
    } else {
        // Small sample: approximate with p=1 (conservative)
        1.0
    };

    (u1, p.clamp(0.0, 1.0))
}

/// Survival function of the standard normal: P(Z > z).
///
/// Uses the relation to erfc: 0.5 × erfc(z / √2).
/// erfc approximated via Abramowitz & Stegun 7.1.26.
fn normal_sf(z: f64) -> f64 {
    0.5 * erfc_approx(z / std::f64::consts::SQRT_2)
}

/// Approximate complementary error function (Abramowitz & Stegun 7.1.26).
fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    poly * (-x * x).exp()
}

/// Benjamini-Hochberg FDR correction in-place (sorted by p-value).
fn bh_correct(markers: &mut [ClusterMarker]) {
    if markers.is_empty() {
        return;
    }
    let m = markers.len();
    // Sort by p_value ascending; keep original index to restore order
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| {
        markers[a]
            .p_value
            .partial_cmp(&markers[b].p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut padj = vec![1.0f64; m];
    let mut cummin = f64::INFINITY;
    // BH: iterate from largest p to smallest
    for rank in (0..m).rev() {
        let i = order[rank];
        let bh = markers[i].p_value * m as f64 / (rank + 1) as f64;
        cummin = cummin.min(bh);
        padj[i] = cummin.min(1.0);
    }

    for (i, p) in padj.into_iter().enumerate() {
        markers[i].padj = p;
    }
}
