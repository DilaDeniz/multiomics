//! Doublet detection via the Scrublet algorithm.
//!
//! Reference: Wolock, Lopez & Klein (2019) "Scrublet: Computational Identification of
//! Cell Doublets in Single-Cell Transcriptomic Data" Cell Systems 8(4):281.
//!
//! Simulates artificial doublets by combining pairs of real cells, embeds real and
//! simulated cells together with PCA, and scores each real cell by the fraction of
//! its kNN that are simulated doublets.

use ndarray::Array2;

/// Per-cell doublet scores and classification.
pub struct DoubletScores {
    /// Doublet score per cell (0.0 = singlet, 1.0 = doublet).
    pub scores: Vec<f32>,
    /// Score threshold used for classification.
    pub threshold: f32,
    /// Which cells are predicted doublets.
    pub is_doublet: Vec<bool>,
    /// Expected doublet rate used for reporting.
    pub expected_rate: f32,
}

/// Detect doublets in a normalized expression matrix.
///
/// `norm_matrix`: `[n_cells, n_genes]` log-normalized counts.
/// `n_simulated`: number of synthetic doublets to generate (default: `min(10_000, n_cells * 2)`).
/// `n_neighbors`: kNN size for scoring (default: 30).
/// `n_pcs`: PCA components for embedding (default: 30).
pub fn detect_doublets(
    norm_matrix: &Array2<f32>,
    n_simulated: usize,
    n_neighbors: usize,
    n_pcs: usize,
) -> DoubletScores {
    let (n_cells, n_genes) = norm_matrix.dim();

    // --- 1. Simulate doublets via xorshift RNG seeded at 42 -----------------
    let mut rng_state: u64 = 42 ^ 0x9e37_79b9_7f4a_7c15;
    let mut xorshift = move || -> u64 {
        rng_state ^= rng_state << 13;
        rng_state ^= rng_state >> 7;
        rng_state ^= rng_state << 17;
        rng_state
    };

    let n_sim = if n_simulated == 0 {
        (n_cells * 2).min(10_000)
    } else {
        n_simulated
    };

    let mut sim = Array2::<f32>::zeros((n_sim, n_genes));
    for s in 0..n_sim {
        let a = (xorshift() as usize) % n_cells;
        let b = (xorshift() as usize) % n_cells;
        for g in 0..n_genes {
            sim[[s, g]] = norm_matrix[[a, g]] + norm_matrix[[b, g]];
        }
    }

    // --- 2. Combine real + simulated ----------------------------------------
    let n_total = n_cells + n_sim;
    let mut combined = Array2::<f32>::zeros((n_total, n_genes));
    for i in 0..n_cells {
        for g in 0..n_genes {
            combined[[i, g]] = norm_matrix[[i, g]];
        }
    }
    for s in 0..n_sim {
        for g in 0..n_genes {
            combined[[n_cells + s, g]] = sim[[s, g]];
        }
    }

    // --- 3. Center on real cells, project combined into PCA space -----------
    let k = n_pcs.min(n_cells).min(n_genes);
    let embedding = pca_project_combined(&combined, n_cells, k);

    // --- 4. kNN in PCA space (real cells vs all cells) ----------------------
    let k_nn = n_neighbors.min(n_total.saturating_sub(1));
    let scores_real = knn_doublet_scores(&embedding, n_cells, n_total, k_nn);

    // scores for simulated cells (same computation)
    let scores_sim = knn_doublet_scores_sim(&embedding, n_cells, n_total, k_nn);

    // --- 6. Threshold from simulated scores ---------------------------------
    let sim_mean = if scores_sim.is_empty() {
        0.5f32
    } else {
        scores_sim.iter().sum::<f32>() / scores_sim.len() as f32
    };
    let sim_var = if scores_sim.is_empty() {
        0.0f32
    } else {
        scores_sim
            .iter()
            .map(|&x| (x - sim_mean) * (x - sim_mean))
            .sum::<f32>()
            / scores_sim.len() as f32
    };
    let sim_std = sim_var.sqrt();
    let threshold = (sim_mean + 2.0 * sim_std).clamp(0.2, 0.8);

    let is_doublet: Vec<bool> = scores_real.iter().map(|&s| s > threshold).collect();

    let expected_rate = n_sim as f32 / n_total as f32;

    DoubletScores {
        scores: scores_real,
        threshold,
        is_doublet,
        expected_rate,
    }
}

/// Fit PCA on real cells (first `n_real` rows of `combined`), center both on real
/// means, then return the projected embedding `[n_total, k]`.
fn pca_project_combined(combined: &Array2<f32>, n_real: usize, k: usize) -> Array2<f32> {
    let (n_total, n_genes) = combined.dim();
    if k == 0 || n_total == 0 || n_genes == 0 {
        return Array2::zeros((n_total, 0));
    }

    // Column means from real cells only
    let mut means = vec![0.0f32; n_genes];
    for g in 0..n_genes {
        let mut s = 0.0f32;
        for i in 0..n_real {
            s += combined[[i, g]];
        }
        means[g] = s / n_real as f32;
    }

    // Center all rows
    let mut centered = combined.to_owned();
    for i in 0..n_total {
        for g in 0..n_genes {
            centered[[i, g]] -= means[g];
        }
    }

    // Power iteration PCA fitted on real cells only
    let real_view = centered.slice(ndarray::s![..n_real, ..]).to_owned();
    let components = power_iter_pca(&real_view, k); // [n_genes, k]

    // Project all cells: embedding = centered @ components  [n_total, k]
    mat_mul_f32(&centered, &components)
}

/// Power-iteration PCA. Returns the top-k right singular vectors `[n_genes, k]`.
fn power_iter_pca(matrix: &Array2<f32>, k: usize) -> Array2<f32> {
    let (n_cells, n_genes) = matrix.dim();
    let k = k.min(n_cells).min(n_genes);
    if k == 0 {
        return Array2::zeros((n_genes, 0));
    }

    // Initialize V [n_genes, k] with xorshift
    let mut seed: u64 = 0x1234_5678_9abc_def0;
    let mut v = Array2::<f32>::zeros((n_genes, k));
    for g in 0..n_genes {
        for j in 0..k {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            v[[g, j]] = ((seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5) as f32;
        }
    }
    qr_cols_inplace(&mut v);

    // Power iteration: V = (M^T M) V, then QR
    let mt = matrix.t().to_owned(); // [n_genes, n_cells]
    for _ in 0..10 {
        // tmp = M @ V  [n_cells, k]
        let tmp = mat_mul_f32(matrix, &v);
        // V = M^T @ tmp  [n_genes, k]
        v = mat_mul_f32(&mt, &tmp);
        qr_cols_inplace(&mut v);
    }
    v
}

/// Compute doublet scores for real cells (indices 0..n_real) by counting
/// how many of their k nearest neighbors among all cells are simulated.
fn knn_doublet_scores(
    embedding: &Array2<f32>,
    n_real: usize,
    n_total: usize,
    k: usize,
) -> Vec<f32> {
    (0..n_real)
        .map(|i| score_cell(embedding, i, n_real, n_total, k))
        .collect()
}

/// Compute doublet scores for simulated cells (indices n_real..n_total).
fn knn_doublet_scores_sim(
    embedding: &Array2<f32>,
    n_real: usize,
    n_total: usize,
    k: usize,
) -> Vec<f32> {
    (n_real..n_total)
        .map(|i| score_cell(embedding, i, n_real, n_total, k))
        .collect()
}

/// Score cell `i`: fraction of its k nearest neighbors (among all cells except itself)
/// that are simulated (index >= n_real).
fn score_cell(embedding: &Array2<f32>, i: usize, n_real: usize, n_total: usize, k: usize) -> f32 {
    let row_i = embedding.row(i);
    let slice_i = row_i.as_slice().unwrap_or(&[]);

    // Compute distances to all other cells
    let mut dists: Vec<(usize, f32)> = (0..n_total)
        .filter(|&j| j != i)
        .map(|j| {
            let row_j = embedding.row(j);
            let d = euclidean_dist(slice_i, row_j.as_slice().unwrap_or(&[]));
            (j, d)
        })
        .collect();

    dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    dists.truncate(k);

    let n_sim_neighbors = dists.iter().filter(|&&(idx, _)| idx >= n_real).count();
    n_sim_neighbors as f32 / k as f32
}

fn euclidean_dist(a: &[f32], b: &[f32]) -> f32 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f32>()
        .sqrt()
}

fn mat_mul_f32(a: &Array2<f32>, b: &Array2<f32>) -> Array2<f32> {
    let m = a.nrows();
    let k = a.ncols();
    let n = b.ncols();
    let mut c = Array2::<f32>::zeros((m, n));
    for i in 0..m {
        for p in 0..k {
            let aip = a[[i, p]];
            for j in 0..n {
                c[[i, j]] += aip * b[[p, j]];
            }
        }
    }
    c
}

/// Modified Gram-Schmidt orthonormalization of columns of `v` in-place.
fn qr_cols_inplace(v: &mut Array2<f32>) {
    let n = v.nrows();
    let m = v.ncols();
    for j in 0..m {
        let norm: f32 = (0..n).map(|i| v[[i, j]] * v[[i, j]]).sum::<f32>().sqrt();
        if norm > 1e-12 {
            for i in 0..n {
                v[[i, j]] /= norm;
            }
        }
        for l in (j + 1)..m {
            let dot: f32 = (0..n).map(|i| v[[i, j]] * v[[i, l]]).sum();
            for i in 0..n {
                let sub = dot * v[[i, j]];
                v[[i, l]] -= sub;
            }
        }
    }
}
