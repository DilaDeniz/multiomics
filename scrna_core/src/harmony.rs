//! Harmony batch correction for single-cell data.
//!
//! Reference: Korsunsky et al. (2019) "Fast, sensitive and accurate integration of
//! single-cell data with Harmony" Nature Methods 16:1289.
//!
//! Iteratively clusters cells into soft clusters and applies per-batch, per-cluster
//! corrections to the PCA embedding until convergence.

use anyhow::Result;
use ndarray::Array2;

/// Result of Harmony batch integration.
pub struct HarmonyResult {
    /// Batch-corrected PCA embedding `[n_cells, n_pcs]`.
    pub embedding: Array2<f32>,
    /// Soft cluster assignments `[n_cells, n_clusters]`.
    pub soft_clusters: Array2<f32>,
    /// Number of iterations until convergence.
    pub n_iter: usize,
}

/// Run Harmony batch correction on a PCA embedding.
///
/// `batch_labels`: batch ID per cell (0-indexed, contiguous).
/// `n_clusters`: number of soft clusters (default: `min(100, n_cells / 10)`).
/// `theta`: diversity penalty (default 1.0).
/// `sigma`: bandwidth for soft cluster assignments (default 0.1).
pub fn harmony_integrate(
    pca_embedding: &Array2<f32>,
    batch_labels: &[u32],
    n_clusters: usize,
    max_iter: usize,
    theta: f32,
    sigma: f32,
) -> Result<HarmonyResult> {
    let (n_cells, n_pcs) = pca_embedding.dim();
    if n_cells == 0 {
        anyhow::bail!("empty embedding");
    }
    if batch_labels.len() != n_cells {
        anyhow::bail!(
            "batch_labels length ({}) != n_cells ({})",
            batch_labels.len(),
            n_cells
        );
    }

    let n_batches = batch_labels.iter().max().copied().unwrap_or(0) as usize + 1;
    let n_clust = if n_clusters == 0 {
        n_clusters.max(1).max((n_cells / 10).min(100))
    } else {
        n_clusters
    };
    let n_clust = n_clust.min(n_cells);

    // Count cells per batch
    let mut batch_counts = vec![0usize; n_batches];
    for &b in batch_labels {
        batch_counts[b as usize] += 1;
    }

    // Z: working corrected embedding [n_cells, n_pcs]
    let mut z = pca_embedding.to_owned();

    // Initialize centroids C [n_pcs, n_clusters] via k-means++
    let mut c = kmeans_pp_init(&z, n_clust);

    let mut r = Array2::<f32>::zeros((n_cells, n_clust)); // soft assignments
    let mut n_iter = 0usize;

    for iter in 0..max_iter {
        let z_prev = z.clone();

        // ---- Step 1: soft cluster assignment --------------------------------
        // dist[i,k] = ||Z[i,:] - C[:,k]||^2
        for i in 0..n_cells {
            for k in 0..n_clust {
                let mut d2 = 0.0f32;
                for p in 0..n_pcs {
                    let diff = z[[i, p]] - c[[p, k]];
                    d2 += diff * diff;
                }
                r[[i, k]] = (-d2 / sigma).exp();
            }
            // Row-normalize
            let row_sum: f32 = (0..n_clust).map(|k| r[[i, k]]).sum::<f32>().max(1e-12);
            for k in 0..n_clust {
                r[[i, k]] /= row_sum;
            }
        }

        // Diversity penalty
        // O[b,k] = sum_{i in batch b} R[i,k]
        // E[b,k] = n_cells_in_batch_b * sum_i(R[i,k]) / n_cells
        let mut o = vec![vec![0.0f32; n_clust]; n_batches];
        let mut cluster_sum = vec![0.0f32; n_clust];
        for i in 0..n_cells {
            let b = batch_labels[i] as usize;
            for k in 0..n_clust {
                o[b][k] += r[[i, k]];
                cluster_sum[k] += r[[i, k]];
            }
        }
        let mut e = vec![vec![0.0f32; n_clust]; n_batches];
        for b in 0..n_batches {
            for k in 0..n_clust {
                e[b][k] = batch_counts[b] as f32 * cluster_sum[k] / n_cells as f32;
            }
        }

        // Reweight R: R[i,k] *= (E[b_i,k] / O[b_i,k])^theta
        for i in 0..n_cells {
            let b = batch_labels[i] as usize;
            for k in 0..n_clust {
                let ratio = e[b][k] / o[b][k].max(1e-12);
                r[[i, k]] *= ratio.powf(theta);
            }
            // Re-normalize
            let row_sum: f32 = (0..n_clust).map(|k| r[[i, k]]).sum::<f32>().max(1e-12);
            for k in 0..n_clust {
                r[[i, k]] /= row_sum;
            }
        }

        // ---- Step 2: compute correction ------------------------------------
        // phi[b,k] = weighted centroid of batch b in cluster k
        // phi shape: [n_batches, n_clust, n_pcs]
        let mut phi = vec![vec![vec![0.0f32; n_pcs]; n_clust]; n_batches];
        let mut phi_w = vec![vec![0.0f32; n_clust]; n_batches]; // sum of weights
        for i in 0..n_cells {
            let b = batch_labels[i] as usize;
            for k in 0..n_clust {
                let w = r[[i, k]];
                phi_w[b][k] += w;
                for p in 0..n_pcs {
                    phi[b][k][p] += w * z[[i, p]];
                }
            }
        }
        for b in 0..n_batches {
            for k in 0..n_clust {
                let denom = phi_w[b][k].max(1e-12);
                for val in phi[b][k].iter_mut() {
                    *val /= denom;
                }
            }
        }

        // Update centroids: C[:,k] = sum_b(phi[b,k] * weight_b_k) / sum(R[:,k])
        for k in 0..n_clust {
            let total_w = cluster_sum[k].max(1e-12);
            for p in 0..n_pcs {
                let mut val = 0.0f32;
                for b in 0..n_batches {
                    val += phi[b][k][p] * phi_w[b][k];
                }
                c[[p, k]] = val / total_w;
            }
        }

        // Correct each cell: Z_new[i,:] = Z[i,:] - sum_k R[i,k]*(phi[b_i,k] - C[:,k])
        let mut z_new = z.clone();
        for i in 0..n_cells {
            let b = batch_labels[i] as usize;
            for k in 0..n_clust {
                let w = r[[i, k]];
                for p in 0..n_pcs {
                    z_new[[i, p]] -= w * (phi[b][k][p] - c[[p, k]]);
                }
            }
        }

        // Convergence check
        let delta = frobenius_relative(&z_new, &z_prev);
        z = z_new;
        n_iter = iter + 1;
        if delta < 1e-4 {
            break;
        }
    }

    Ok(HarmonyResult {
        embedding: z,
        soft_clusters: r,
        n_iter,
    })
}

/// k-means++ initialization of `k` centroids from `data` `[n_cells, n_pcs]`.
/// Returns `[n_pcs, k]`.
pub fn kmeans_pp_init(data: &Array2<f32>, k: usize) -> Array2<f32> {
    let (n_cells, n_pcs) = data.dim();
    let k = k.min(n_cells);
    let mut centroids = Array2::<f32>::zeros((n_pcs, k));
    if k == 0 || n_cells == 0 {
        return centroids;
    }

    // xorshift RNG
    let mut seed: u64 = 0xabcd_ef01_2345_6789;
    let mut rng = move || -> f64 {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };

    // First centroid: random cell
    let first = (rng() * n_cells as f64) as usize % n_cells;
    for p in 0..n_pcs {
        centroids[[p, 0]] = data[[first, p]];
    }

    let mut min_dists = vec![f32::INFINITY; n_cells];

    for ck in 1..k {
        // Update min distances to nearest already-chosen centroid
        for i in 0..n_cells {
            let mut d2 = 0.0f32;
            for p in 0..n_pcs {
                let diff = data[[i, p]] - centroids[[p, ck - 1]];
                d2 += diff * diff;
            }
            if d2 < min_dists[i] {
                min_dists[i] = d2;
            }
        }
        // Sample next centroid proportional to min_dists^2
        let total: f64 = min_dists.iter().map(|&d| d as f64).sum::<f64>();
        let target = rng() * total;
        let mut cumsum = 0.0f64;
        let mut chosen = 0usize;
        for (i, &d) in min_dists.iter().enumerate() {
            cumsum += d as f64;
            if cumsum >= target {
                chosen = i;
                break;
            }
        }
        for p in 0..n_pcs {
            centroids[[p, ck]] = data[[chosen, p]];
        }
    }

    centroids
}

/// Relative Frobenius norm: ||A - B||_F / ||B||_F.
fn frobenius_relative(a: &Array2<f32>, b: &Array2<f32>) -> f32 {
    let diff_sq: f32 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum();
    let b_sq: f32 = b.iter().map(|x| x * x).sum::<f32>().max(1e-12);
    (diff_sq / b_sq).sqrt()
}
