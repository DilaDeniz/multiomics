//! KNN graph construction via random-projection approximate nearest neighbours.
//!
//! Uses 20 random unit projections to find candidate neighbours, then
//! resolves exact Euclidean distances among candidates.

use ndarray::Array2;

/// Sparse KNN graph over cells.
pub struct KnnGraph {
    /// Number of cells.
    pub n_cells: usize,
    /// Number of neighbours per cell.
    pub k: usize,
    /// `n_cells × k` neighbour lists (cell indices).
    pub neighbors: Vec<Vec<u32>>,
}

/// Build an approximate KNN graph from a cell embedding.
///
/// Uses random-projection LSH (20 projections, ±15 position window) to
/// identify candidates, then ranks by exact Euclidean distance.
///
/// `embedding`: `[n_cells, n_components]`.
pub fn build_knn_graph(embedding: &Array2<f32>, k: usize) -> KnnGraph {
    let (n_cells, n_components) = embedding.dim();
    if n_cells == 0 || n_components == 0 {
        return KnnGraph {
            n_cells,
            k,
            neighbors: vec![Vec::new(); n_cells],
        };
    }

    const N_PROJ: usize = 20;
    const WINDOW: usize = 15;

    // Generate random unit projection vectors via xorshift + Box-Muller
    let projections = random_projections(n_components, N_PROJ);

    // Project each cell onto each projection vector: proj_values[p][i]
    let proj_values: Vec<Vec<f32>> = projections
        .iter()
        .map(|proj| {
            (0..n_cells)
                .map(|i| {
                    let row = embedding.row(i);
                    row.iter().zip(proj.iter()).map(|(a, b)| a * b).sum()
                })
                .collect()
        })
        .collect();

    // For each projection, sort cells by projection value and record sorted order
    let sorted_orders: Vec<Vec<usize>> = proj_values
        .iter()
        .map(|pv| {
            let mut order: Vec<usize> = (0..n_cells).collect();
            order.sort_by(|&a, &b| {
                pv[a]
                    .partial_cmp(&pv[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            order
        })
        .collect();

    // For each cell, invert to find its position in each sorted order
    let mut pos_in_order: Vec<Vec<usize>> = vec![vec![0usize; N_PROJ]; n_cells];
    for (p, order) in sorted_orders.iter().enumerate() {
        for (rank, &cell) in order.iter().enumerate() {
            pos_in_order[cell][p] = rank;
        }
    }

    // Build candidate sets per cell
    let k_actual = k.min(n_cells.saturating_sub(1));
    let mut neighbors: Vec<Vec<u32>> = Vec::with_capacity(n_cells);

    for (i, cell_positions) in pos_in_order.iter().enumerate() {
        let mut candidates = ahash::AHashSet::new();
        for (p, &rank) in cell_positions.iter().enumerate() {
            let lo = rank.saturating_sub(WINDOW);
            let hi = (rank + WINDOW + 1).min(n_cells);
            for &c in &sorted_orders[p][lo..hi] {
                if c != i {
                    candidates.insert(c);
                }
            }
        }

        // Compute exact distances to all candidates
        let row_i = embedding.row(i);
        let mut dists: Vec<(u32, f32)> = candidates
            .into_iter()
            .map(|c| {
                let row_c = embedding.row(c);
                let d = euclidean_dist(
                    row_i.as_slice().unwrap_or(&[]),
                    row_c.as_slice().unwrap_or(&[]),
                );
                (c as u32, d)
            })
            .collect();

        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        dists.truncate(k_actual);
        neighbors.push(dists.into_iter().map(|(c, _)| c).collect());
    }

    KnnGraph {
        n_cells,
        k: k_actual,
        neighbors,
    }
}

/// Euclidean distance between two slices of equal length.
fn euclidean_dist(a: &[f32], b: &[f32]) -> f32 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f32>()
        .sqrt()
}

/// Generate `n_proj` random unit vectors of dimension `dim` using
/// xorshift64 + Box-Muller transform.
fn random_projections(dim: usize, n_proj: usize) -> Vec<Vec<f32>> {
    let mut seed: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut rng = move || -> f64 {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };

    let mut out = Vec::with_capacity(n_proj);
    for _ in 0..n_proj {
        let mut v: Vec<f32> = Vec::with_capacity(dim);
        let mut i = 0;
        while i < dim {
            // Box-Muller transform
            let u1 = rng().max(1e-300);
            let u2 = rng();
            let r = (-2.0 * u1.ln()).sqrt();
            let theta = std::f64::consts::TAU * u2;
            v.push((r * theta.cos()) as f32);
            if i + 1 < dim {
                v.push((r * theta.sin()) as f32);
            }
            i += 2;
        }
        // Normalise
        let norm: f32 = v.iter().map(|x| x * x).sum::<f32>().sqrt();
        if norm > 1e-12 {
            for x in &mut v {
                *x /= norm;
            }
        }
        out.push(v);
    }
    out
}
