//! UMAP (Uniform Manifold Approximation and Projection) dimensionality reduction.
//!
//! Reference: McInnes, Healy & Melville (2018) "UMAP: Uniform Manifold Approximation
//! and Projection for Dimension Reduction" arXiv:1802.03426.
//!
//! Implements Phase 1 (fuzzy simplicial set construction via KNN + sigma search)
//! and Phase 2 (SGD-based embedding optimization with attraction + repulsion).

use anyhow::Result;
use ndarray::Array2;
use rayon::prelude::*;

// ---------------------------------------------------------------------------
// Public result type
// ---------------------------------------------------------------------------

/// Output of a UMAP run.
pub struct UmapResult {
    /// Low-dimensional embedding `[n_cells, 2]`.
    pub embedding: Array2<f64>,
    /// Number of SGD epochs executed.
    pub n_epochs: usize,
}

// ---------------------------------------------------------------------------
// Phase 1: Fuzzy simplicial set
// ---------------------------------------------------------------------------

/// Build the fuzzy high-dimensional graph (Phase 1 of UMAP).
///
/// Returns `(adjacency, rho)` where:
/// - `adjacency[i]` is a list of `(neighbor_idx, weight)` pairs for cell `i`
///   (symmetrized: `w_sym = w_ij + w_ji - w_ij * w_ji`).
/// - `rho[i]` is the distance to the nearest neighbor of cell `i`.
pub fn compute_fuzzy_graph(
    data: &Array2<f64>,
    n_neighbors: usize,
) -> (Vec<Vec<(usize, f64)>>, Vec<f64>) {
    let n_cells = data.nrows();
    let k = n_neighbors.min(n_cells.saturating_sub(1));

    // Step 1: KNN
    let knn = knn_all(data, k);

    // Step 2: per-cell rho and sigma
    let mut rho = vec![0.0f64; n_cells];
    let mut sigma = vec![1.0f64; n_cells];

    for i in 0..n_cells {
        if knn[i].is_empty() {
            continue;
        }
        rho[i] = knn[i][0].1; // distance to nearest neighbour
        sigma[i] = find_sigma(&knn[i], rho[i], k);
    }

    // Step 3: directed weights
    // w_ij = exp(-(d_ij - rho_i) / sigma_i)
    let mut directed: Vec<Vec<(usize, f64)>> = vec![vec![]; n_cells];
    for i in 0..n_cells {
        for &(j, d) in &knn[i] {
            let w = if sigma[i] > 0.0 {
                (-(d - rho[i]).max(0.0) / sigma[i]).exp()
            } else {
                if (d - rho[i]).abs() < 1e-12 {
                    1.0
                } else {
                    0.0
                }
            };
            directed[i].push((j, w));
        }
    }

    // Step 4: symmetrize w_sym = w_ij + w_ji - w_ij * w_ji
    // Build a lookup: for each (i,j) pair we need both w_ij and w_ji.
    // Use a flat upper-triangle map via sorted pair keys.
    use std::collections::HashMap;
    // Store directed weights in a map keyed by (min, max) -> (w_from_min, w_from_max)
    let mut pair_map: HashMap<(usize, usize), (f64, f64)> = HashMap::new();

    for (i, dir_i) in directed.iter().enumerate() {
        for &(j, w) in dir_i {
            let key = if i < j { (i, j) } else { (j, i) };
            let entry = pair_map.entry(key).or_insert((0.0, 0.0));
            if i <= j {
                entry.0 = w; // w from smaller-index cell
            } else {
                entry.1 = w; // w from larger-index cell
            }
        }
    }

    let mut adjacency: Vec<Vec<(usize, f64)>> = vec![vec![]; n_cells];
    for (&(a, b), &(w_ab, w_ba)) in &pair_map {
        let w_sym = w_ab + w_ba - w_ab * w_ba;
        adjacency[a].push((b, w_sym));
        adjacency[b].push((a, w_sym));
    }

    // Sort each adjacency list by neighbour index for deterministic ordering.
    for nbrs in &mut adjacency {
        nbrs.sort_unstable_by_key(|&(j, _)| j);
    }

    (adjacency, rho)
}

/// Binary search for sigma such that sum_j exp(-(d_ij - rho_i)/sigma) ≈ log2(k).
fn find_sigma(neighbors: &[(usize, f64)], rho: f64, k: usize) -> f64 {
    let target = (k as f64).log2();
    let mut lo = 1e-10f64;
    let mut hi = 1e10f64;

    for _ in 0..64 {
        let mid = (lo + hi) / 2.0;
        let s: f64 = neighbors
            .iter()
            .map(|&(_, d)| (-(d - rho).max(0.0) / mid).exp())
            .sum();
        if s < target {
            hi = mid;
        } else {
            lo = mid;
        }
        if (hi - lo) / lo.max(1e-20) < 1e-8 {
            break;
        }
    }
    (lo + hi) / 2.0
}

// ---------------------------------------------------------------------------
// KNN (exact flat scan for small datasets, random-projection buckets otherwise)
// ---------------------------------------------------------------------------

/// Compute k nearest neighbours for every cell.
/// Returns `Vec<Vec<(idx, dist)>>` sorted by ascending distance (excluding self).
fn knn_all(data: &Array2<f64>, k: usize) -> Vec<Vec<(usize, f64)>> {
    let n_cells = data.nrows();
    if n_cells <= 10_000 {
        knn_flat(data, k)
    } else {
        knn_rp(data, k)
    }
}

/// Exact KNN via flat pairwise scan — O(n²) but straightforward.
fn knn_flat(data: &Array2<f64>, k: usize) -> Vec<Vec<(usize, f64)>> {
    let n_cells = data.nrows();
    (0..n_cells)
        .into_par_iter()
        .map(|i| {
            let mut dists: Vec<(usize, f64)> = (0..n_cells)
                .filter(|&j| j != i)
                .map(|j| (j, euclidean_sq(data, i, j).sqrt()))
                .collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            dists.truncate(k);
            dists
        })
        .collect()
}

/// Approximate KNN via random projection partitioning for large datasets.
fn knn_rp(data: &Array2<f64>, k: usize) -> Vec<Vec<(usize, f64)>> {
    let n_cells = data.nrows();
    let n_dims = data.ncols();
    let n_trees = 4usize;
    let max_leaf = (k * 10).max(50).min(n_cells);

    // For each cell, accumulate candidate neighbour sets
    let mut candidates: Vec<std::collections::HashSet<usize>> = (0..n_cells)
        .map(|_| std::collections::HashSet::new())
        .collect();

    let mut seed: u64 = 0x1234_5678_abcd_ef01;
    let mut xorshift = move || -> u64 {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };

    for _tree in 0..n_trees {
        // Random projection vector
        let proj: Vec<f64> = (0..n_dims)
            .map(|_| {
                let u = xorshift();
                // Box-Muller pair, take just one
                let u1 = (u >> 32) as f64 / u32::MAX as f64 + 1e-300;
                let u2 = (u & 0xFFFF_FFFF) as f64 / u32::MAX as f64;
                (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
            })
            .collect();

        // Project all cells
        let mut projections: Vec<(usize, f64)> = (0..n_cells)
            .map(|i| {
                let p: f64 = (0..n_dims).map(|d| data[[i, d]] * proj[d]).sum();
                (i, p)
            })
            .collect();
        projections.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        // Slide a window of max_leaf and add mutual candidates
        let window = max_leaf.min(n_cells);
        for start in (0..n_cells).step_by(window / 2 + 1) {
            let end = (start + window).min(n_cells);
            let bucket: Vec<usize> = projections[start..end]
                .iter()
                .map(|&(idx, _)| idx)
                .collect();
            for &ci in &bucket {
                for &cj in &bucket {
                    if ci != cj {
                        candidates[ci].insert(cj);
                    }
                }
            }
        }
    }

    // Exact KNN within candidate sets
    (0..n_cells)
        .into_par_iter()
        .map(|i| {
            let cands = &candidates[i];
            let mut dists: Vec<(usize, f64)> = cands
                .iter()
                .map(|&j| (j, euclidean_sq(data, i, j).sqrt()))
                .collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            dists.truncate(k);
            dists
        })
        .collect()
}

#[inline]
fn euclidean_sq(data: &Array2<f64>, i: usize, j: usize) -> f64 {
    let row_i = data.row(i);
    let row_j = data.row(j);
    row_i
        .iter()
        .zip(row_j.iter())
        .map(|(a, b)| (a - b) * (a - b))
        .sum()
}

// ---------------------------------------------------------------------------
// Phase 2: SGD embedding
// ---------------------------------------------------------------------------

/// Run UMAP dimensionality reduction.
///
/// # Arguments
/// - `data` — `[n_cells, n_dims]` input (e.g. PCA embedding).
/// - `n_neighbors` — number of nearest neighbours (typically 15).
/// - `n_epochs` — SGD iterations (typically 200).
/// - `min_dist` — minimum distance in low-dim space (typically 0.1).
/// - `learning_rate` — initial SGD step size (typically 1.0).
/// - `seed` — RNG seed for reproducible results.
pub fn run_umap(
    data: &Array2<f64>,
    n_neighbors: usize,
    n_epochs: usize,
    min_dist: f64,
    learning_rate: f64,
    seed: u64,
) -> Result<UmapResult> {
    let n_cells = data.nrows();
    if n_cells < 2 {
        anyhow::bail!("need at least 2 cells for UMAP");
    }
    if data.ncols() == 0 {
        anyhow::bail!("input has zero dimensions");
    }

    // Phase 1
    let (adjacency, _rho) = compute_fuzzy_graph(data, n_neighbors);

    // Precompute a, b from min_dist
    let (a, b) = ab_params(min_dist);

    // Initialise embedding with scaled random normal
    let mut emb = init_embedding(n_cells, seed);

    // Collect all edges (i, j, w) — only i < j to avoid double-counting.
    // Sort for deterministic SGD update order regardless of HashMap iteration order.
    let mut edges: Vec<(usize, usize, f64)> = Vec::new();
    for (i, nbrs) in adjacency.iter().enumerate() {
        for &(j, w) in nbrs {
            if i < j {
                edges.push((i, j, w));
            }
        }
    }
    edges.sort_unstable_by_key(|&(i, j, _)| (i, j));

    let n_neg = 5usize; // negative samples per positive edge
    let epsilon = 1e-4f64;

    // SGD loop — sequential over epochs, parallel within epoch where possible
    let mut rng_state = seed ^ 0xdeadbeef;

    for epoch in 0..n_epochs {
        let lr = learning_rate * (1.0 - epoch as f64 / n_epochs as f64);

        // Positive (attraction) edges
        // Sequential to allow in-place mutation without locks
        for &(i, j, _w) in &edges {
            let dy0 = emb[[i, 0]] - emb[[j, 0]];
            let dy1 = emb[[i, 1]] - emb[[j, 1]];
            let d2 = (dy0 * dy0 + dy1 * dy1).max(1e-12);
            let d = d2.sqrt();

            // Attraction gradient: -2ab d^(2b-2) / (1 + a d^(2b))
            let d_2b = a * d.powf(2.0 * b);
            let grad_coeff = -2.0 * a * b * d.powf(2.0 * b - 2.0) / (1.0 + d_2b);

            let g0 = (grad_coeff * dy0).clamp(-4.0, 4.0);
            let g1 = (grad_coeff * dy1).clamp(-4.0, 4.0);

            emb[[i, 0]] += lr * g0;
            emb[[i, 1]] += lr * g1;
            emb[[j, 0]] -= lr * g0;
            emb[[j, 1]] -= lr * g1;
        }

        // Negative (repulsion) samples
        for &(i, _, _) in &edges {
            for _ in 0..n_neg {
                let j = xorshift_next(&mut rng_state) as usize % n_cells;
                if j == i {
                    continue;
                }
                let dy0 = emb[[i, 0]] - emb[[j, 0]];
                let dy1 = emb[[i, 1]] - emb[[j, 1]];
                let d2 = dy0 * dy0 + dy1 * dy1;
                let d_2b = a * d2.powf(b);

                // Repulsion gradient: 2b / (epsilon + d^2) / (1 + a d^(2b))
                let grad_coeff = 2.0 * b / (epsilon + d2) / (1.0 + d_2b);

                let g0 = (grad_coeff * dy0).clamp(-4.0, 4.0);
                let g1 = (grad_coeff * dy1).clamp(-4.0, 4.0);

                emb[[i, 0]] += lr * g0;
                emb[[i, 1]] += lr * g1;
            }
        }
    }

    Ok(UmapResult {
        embedding: emb,
        n_epochs,
    })
}

/// Convenience wrapper: run UMAP with default hyperparameters on a PCA embedding.
pub fn umap_from_pca(
    pca_embedding: &Array2<f64>,
    n_neighbors: usize,
    n_epochs: usize,
) -> Result<UmapResult> {
    run_umap(pca_embedding, n_neighbors, n_epochs, 0.1, 1.0, 42)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Compute UMAP curve parameters a and b from `min_dist`.
///
/// Fits ψ(d) = 1 / (1 + a·d^(2b)) to the step function:
///   ψ(d) = 1 if d ≤ min_dist else exp(-(d - min_dist))
///
/// Uses hardcoded defaults for min_dist=0.1 (a≈1.929, b≈0.7915) and a
/// binary-search refinement for other values.
fn ab_params(min_dist: f64) -> (f64, f64) {
    // For the canonical min_dist=0.1 use the reference values from the paper.
    if (min_dist - 0.1).abs() < 1e-9 {
        return (1.929_001_525_678_12, 0.791_494_504_589_076_3);
    }

    // Otherwise do a simple Gauss-Newton-style fit over a grid of d values.
    // We fit b first (1-D search), then compute a from b analytically.
    // Target: step(d) = 1 if d <= min_dist else exp(-(d - min_dist))
    let n_pts = 300usize;
    let d_max = 3.0f64;
    let ds: Vec<f64> = (0..n_pts)
        .map(|k| d_max * k as f64 / (n_pts - 1) as f64)
        .collect();
    let ys: Vec<f64> = ds
        .iter()
        .map(|&d| {
            if d <= min_dist {
                1.0
            } else {
                (-(d - min_dist)).exp()
            }
        })
        .collect();

    // Search b in [0.1, 2.0]
    let best_b = {
        let mut best_b = 0.7915f64;
        let mut best_err = f64::MAX;
        let steps = 200usize;
        for step in 0..=steps {
            let b_try = 0.1 + 1.9 * step as f64 / steps as f64;
            // For fixed b, a = sum(y * d^(2b)) / sum(d^(4b)) over points where d > min_dist
            // derived from minimising sum (y - 1/(1+a*d^(2b)))^2 analytically → non-trivial,
            // so instead use a simpler moment match: a s.t. average psi matches average y.
            let a_try = compute_a(&ds, &ys, b_try);
            let err: f64 = ds
                .iter()
                .zip(ys.iter())
                .map(|(&d, &y)| {
                    let pred = 1.0 / (1.0 + a_try * d.powf(2.0 * b_try));
                    (y - pred) * (y - pred)
                })
                .sum();
            if err < best_err {
                best_err = err;
                best_b = b_try;
            }
        }
        best_b
    };

    let best_a = compute_a(&ds, &ys, best_b);
    (best_a.max(0.001), best_b)
}

fn compute_a(ds: &[f64], ys: &[f64], b: f64) -> f64 {
    // Minimise sum (y - 1/(1+a*x))^2 w.r.t. a where x = d^(2b).
    // Setting derivative to zero and approximating: a ≈ sum((1-y)/y * 1/x) / n
    // Better: use the least-squares normal equation for the linearised form.
    // Linear approx: if y ~ 1/(1+a*x) then a*x ~ (1-y)/y => a = mean((1-y)/y / x)
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for (&d, &y) in ds.iter().zip(ys.iter()) {
        if d < 1e-12 || y < 1e-9 {
            continue;
        }
        let x = d.powf(2.0 * b);
        // Weighted residual approach
        num += (1.0 - y) * x;
        den += x * x;
    }
    if den < 1e-30 {
        1.929
    } else {
        (num / den).max(0.001)
    }
}

/// Initialize a `[n_cells, 2]` embedding with small random normal values.
fn init_embedding(n_cells: usize, seed: u64) -> Array2<f64> {
    let mut s = seed ^ 0x9e37_79b9_7f4a_7c15u64;
    let mut rng = move || -> f64 {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };

    let mut emb = Array2::<f64>::zeros((n_cells, 2));
    let mut buf = 0.0f64;
    let mut buf_ready = false;
    for i in 0..n_cells {
        for c in 0..2 {
            let v = if buf_ready {
                buf_ready = false;
                buf
            } else {
                let u1 = rng().max(1e-300);
                let u2 = rng();
                let r = (-2.0 * u1.ln()).sqrt() * 0.0001;
                let theta = std::f64::consts::TAU * u2;
                buf = r * theta.sin();
                buf_ready = true;
                r * theta.cos()
            };
            emb[[i, c]] = v;
        }
    }
    emb
}

#[inline]
fn xorshift_next(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    fn random_data(n_cells: usize, n_dims: usize, seed: u64) -> Array2<f64> {
        let mut s = seed;
        let mut rng = move || -> f64 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        Array2::from_shape_fn((n_cells, n_dims), |_| rng())
    }

    #[test]
    fn umap_output_shape() {
        let data = random_data(50, 10, 1);
        let result = run_umap(&data, 10, 50, 0.1, 1.0, 42).expect("umap failed");
        assert_eq!(result.embedding.nrows(), 50);
        assert_eq!(result.embedding.ncols(), 2);
        assert_eq!(result.n_epochs, 50);
    }

    #[test]
    fn umap_is_reproducible() {
        let data = random_data(50, 10, 7);
        let r1 = run_umap(&data, 10, 30, 0.1, 1.0, 99).expect("umap failed");
        let r2 = run_umap(&data, 10, 30, 0.1, 1.0, 99).expect("umap failed");
        for i in 0..50 {
            for c in 0..2 {
                assert!(
                    (r1.embedding[[i, c]] - r2.embedding[[i, c]]).abs() < 1e-15,
                    "embedding not reproducible at [{i},{c}]"
                );
            }
        }
    }

    #[test]
    fn fuzzy_graph_weights_bounded() {
        let data = random_data(40, 8, 3);
        let (adjacency, _rho) = compute_fuzzy_graph(&data, 10);
        for neighbors in &adjacency {
            for &(_j, w) in neighbors {
                assert!((0.0..=1.0).contains(&w), "weight {w} out of [0,1]");
            }
        }
    }
}
