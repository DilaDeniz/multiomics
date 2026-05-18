//! Diffusion Pseudotime (DPT) trajectory inference.
//!
//! Reference: Haghverdi et al. (2016) "Diffusion pseudotime robustly reconstructs
//! lineage branching" Nature Methods 13:845.
//!
//! Orders cells along a developmental trajectory by computing their distance
//! to a root cell in diffusion map space.

use anyhow::Result;
use ndarray::Array2;

use crate::graph::KnnGraph;

/// Result of diffusion pseudotime computation.
pub struct PseudotimeResult {
    /// Pseudotime value per cell (0.0 = root, 1.0 = most distant).
    pub pseudotime: Vec<f64>,
    /// Diffusion map coordinates `[n_cells, n_components]`.
    pub diffusion_map: Array2<f64>,
    /// Index of the root cell.
    pub root_cell: usize,
}

/// Compute diffusion pseudotime from a PCA embedding and kNN graph.
///
/// `root_cell`: if `None`, the cell most distant from the centroid is used as root.
/// `n_diffusion_components`: number of diffusion map dimensions (default 10).
pub fn compute_pseudotime(
    pca_embedding: &Array2<f32>,
    knn_graph: &KnnGraph,
    root_cell: Option<usize>,
    n_diffusion_components: usize,
) -> Result<PseudotimeResult> {
    let n_cells = pca_embedding.nrows();
    if n_cells == 0 {
        anyhow::bail!("empty embedding: no cells");
    }
    if knn_graph.n_cells != n_cells {
        anyhow::bail!(
            "knn_graph.n_cells ({}) != pca_embedding rows ({})",
            knn_graph.n_cells,
            n_cells
        );
    }

    let n_comp = n_diffusion_components.min(n_cells.saturating_sub(1)).max(1);

    // --- 1. Build affinity matrix W (sparse adjacency list) -----------------
    // sigma = median of kNN distances
    let sigma = compute_sigma(pca_embedding, knn_graph);
    let sigma2 = (sigma * sigma).max(1e-12);

    // W[i] = list of (j, weight) where weight = exp(-dist^2/sigma^2), symmetric
    let mut w: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_cells];
    for i in 0..n_cells {
        let row_i = pca_embedding.row(i);
        let slice_i = row_i.as_slice().unwrap_or(&[]);
        for &nb in &knn_graph.neighbors[i] {
            let j = nb as usize;
            let row_j = pca_embedding.row(j);
            let d2 = euclidean_dist_sq(slice_i, row_j.as_slice().unwrap_or(&[]));
            let weight = (-(d2 as f64) / sigma2).exp();
            w[i].push((j, weight));
            w[j].push((i, weight)); // ensure symmetry
        }
    }

    // Symmetrize: average duplicate entries
    let w = symmetrize(w, n_cells);

    // --- 2. Diffusion operator T = D^{-1} W (row-normalized) ----------------
    let mut row_sums = vec![0.0f64; n_cells];
    for (i, row) in w.iter().enumerate() {
        for &(_, wt) in row {
            row_sums[i] += wt;
        }
    }

    // t[i] = [(j, w[i,j] / row_sum[i])]
    let t: Vec<Vec<(usize, f64)>> = w
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let d = row_sums[i].max(1e-300);
            row.iter().map(|&(j, wt)| (j, wt / d)).collect()
        })
        .collect();

    // --- 3. Power iteration for top eigenvectors of T -----------------------
    // We need n_comp + 1 eigenvectors (first is trivial), use n_comp+1 columns
    let n_vecs = (n_comp + 1).min(n_cells);
    let mut v = init_random_matrix(n_cells, n_vecs);
    qr_cols_f64(&mut v);

    for _ in 0..20 {
        v = sparse_matvec(&t, &v);
        qr_cols_f64(&mut v);
    }

    // Column norms of T@V give approximate eigenvalues
    // Skip the first eigenvector (constant), use columns 1..n_comp+1
    let n_use = (n_vecs - 1).max(1);
    let mut diffusion_map = Array2::<f64>::zeros((n_cells, n_use));
    let mut eigenvalues = vec![1.0f64; n_use];

    for k in 0..n_use {
        let col_idx = k + 1; // skip first
        if col_idx >= n_vecs {
            break;
        }
        // Eigenvalue estimate: ||T v_k|| / ||v_k|| ≈ column norm of T@V / 1 (already normalized)
        // After convergence, T@V ≈ lambda * V, so lambda ≈ (T@V)[:,k] dot V[:,k]
        let tv_col = sparse_matvec_single(&t, &v, col_idx);
        let lam: f64 = tv_col
            .iter()
            .zip(v.column(col_idx).iter())
            .map(|(a, b)| a * b)
            .sum();
        eigenvalues[k] = lam.abs().max(1e-12);
        for i in 0..n_cells {
            diffusion_map[[i, k]] = v[[i, col_idx]];
        }
    }

    // --- 4. Root cell selection ---------------------------------------------
    let root = match root_cell {
        Some(r) if r < n_cells => r,
        _ => most_peripheral_cell(pca_embedding),
    };

    // --- 5. DPT distance from root ------------------------------------------
    let mut dpt = vec![0.0f64; n_cells];
    for i in 0..n_cells {
        let mut d2 = 0.0f64;
        for k in 0..n_use {
            let diff = diffusion_map[[i, k]] - diffusion_map[[root, k]];
            d2 += diff * diff / eigenvalues[k].powi(2);
        }
        dpt[i] = d2.sqrt();
    }

    // Normalize to [0, 1]
    let min_d = dpt.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_d = dpt.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = (max_d - min_d).max(1e-12);
    let pseudotime: Vec<f64> = dpt.iter().map(|&d| (d - min_d) / range).collect();

    Ok(PseudotimeResult {
        pseudotime,
        diffusion_map,
        root_cell: root,
    })
}

/// Compute sigma as the median of all kNN distances.
fn compute_sigma(pca: &Array2<f32>, knn: &KnnGraph) -> f64 {
    let mut dists: Vec<f64> = Vec::new();
    for i in 0..knn.n_cells {
        let row_i = pca.row(i);
        let slice_i = row_i.as_slice().unwrap_or(&[]);
        for &nb in &knn.neighbors[i] {
            let j = nb as usize;
            let row_j = pca.row(j);
            let d = euclidean_dist_sq(slice_i, row_j.as_slice().unwrap_or(&[])).sqrt() as f64;
            dists.push(d);
        }
    }
    if dists.is_empty() {
        return 1.0;
    }
    dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = dists.len() / 2;
    if dists.len() % 2 == 1 {
        dists[mid]
    } else {
        (dists[mid - 1] + dists[mid]) / 2.0
    }
}

/// Symmetrize adjacency list by averaging duplicate (i,j) entries.
fn symmetrize(mut w: Vec<Vec<(usize, f64)>>, _n_cells: usize) -> Vec<Vec<(usize, f64)>> {
    // Sort and merge duplicates, then average
    for row in w.iter_mut() {
        row.sort_by_key(|&(j, _)| j);
        let mut merged: Vec<(usize, f64)> = Vec::new();
        for (j, wt) in row.drain(..) {
            if let Some(last) = merged.last_mut() {
                if last.0 == j {
                    last.1 = (last.1 + wt) / 2.0;
                    continue;
                }
            }
            merged.push((j, wt));
        }
        *row = merged;
    }
    w
}

/// Sparse matrix-vector product: result[i,:] = sum_j T[i,j] * v[j,:].
pub fn sparse_matvec(adjacency: &[Vec<(usize, f64)>], v: &Array2<f64>) -> Array2<f64> {
    let n = adjacency.len();
    let cols = v.ncols();
    let mut out = Array2::<f64>::zeros((n, cols));
    for i in 0..n {
        for &(j, w) in &adjacency[i] {
            for c in 0..cols {
                out[[i, c]] += w * v[[j, c]];
            }
        }
    }
    out
}

/// Sparse matrix-vector product for a single column of v.
fn sparse_matvec_single(adjacency: &[Vec<(usize, f64)>], v: &Array2<f64>, col: usize) -> Vec<f64> {
    let n = adjacency.len();
    let mut out = vec![0.0f64; n];
    for i in 0..n {
        for &(j, w) in &adjacency[i] {
            out[i] += w * v[[j, col]];
        }
    }
    out
}

/// Return the cell index most distant from the centroid of `pca`.
fn most_peripheral_cell(pca: &Array2<f32>) -> usize {
    let n = pca.nrows();
    let d = pca.ncols();
    if n == 0 {
        return 0;
    }
    let mut centroid = vec![0.0f32; d];
    for i in 0..n {
        for g in 0..d {
            centroid[g] += pca[[i, g]];
        }
    }
    for val in centroid.iter_mut() {
        *val /= n as f32;
    }
    let mut best = 0usize;
    let mut best_dist = f32::NEG_INFINITY;
    for i in 0..n {
        let row = pca.row(i);
        let d2: f32 = row
            .iter()
            .zip(centroid.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum();
        if d2 > best_dist {
            best_dist = d2;
            best = i;
        }
    }
    best
}

fn euclidean_dist_sq(a: &[f32], b: &[f32]) -> f32 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// Initialize a random `[n, m]` matrix with xorshift RNG.
fn init_random_matrix(n: usize, m: usize) -> Array2<f64> {
    let mut seed: u64 = 0xfeed_face_dead_beef;
    let mut out = Array2::<f64>::zeros((n, m));
    for i in 0..n {
        for j in 0..m {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            out[[i, j]] = (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
        }
    }
    out
}

/// Modified Gram-Schmidt orthonormalization of columns (in-place, f64 version).
fn qr_cols_f64(v: &mut Array2<f64>) {
    let n = v.nrows();
    let m = v.ncols();
    for j in 0..m {
        let norm: f64 = (0..n).map(|i| v[[i, j]] * v[[i, j]]).sum::<f64>().sqrt();
        if norm > 1e-14 {
            for i in 0..n {
                v[[i, j]] /= norm;
            }
        }
        for l in (j + 1)..m {
            let dot: f64 = (0..n).map(|i| v[[i, j]] * v[[i, l]]).sum();
            for i in 0..n {
                let sub = dot * v[[i, j]];
                v[[i, l]] -= sub;
            }
        }
    }
}
