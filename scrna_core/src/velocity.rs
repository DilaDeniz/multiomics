//! RNA velocity estimation (steady-state model).
//!
//! Reference: La Manno et al. (2018) "RNA velocity of single cells"
//! Nature 560:494–498.
//!
//! Estimates per-gene degradation rates (γ) from the steady-state
//! relationship `spliced = γ × unspliced`, then computes residual
//! velocity per cell: `v_i = s_i − γ·u_i`.

use anyhow::Result;
use ndarray::Array2;
use serde::{Deserialize, Serialize};

/// Per-gene velocity estimates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneVelocity {
    pub gene_id: String,
    /// Transcription rate (gamma * steady_state_ratio).
    pub gamma: f64,
    /// R² of the steady-state fit.
    pub r_squared: f64,
    /// Velocity residual per cell: v_i = ds/dt - γ·u (positive = increasing).
    pub velocity: Vec<f64>,
}

/// Full RNA velocity result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VelocityResult {
    pub gene_velocities: Vec<GeneVelocity>,
    /// Cell-level velocity vectors in PCA space (n_cells × 2).
    pub velocity_embedding: Array2<f64>,
    pub n_cells: usize,
    pub n_velocity_genes: usize,
}

/// Compute RNA velocity using the steady-state model (La Manno 2018).
///
/// For each gene:
/// 1. Estimate γ as the median ratio `spliced / unspliced` among cells in
///    the top 5% of unspliced expression (upper-quantile regression).
/// 2. Compute residual velocity: `v_i = s_i − γ·u_i`.
/// 3. Keep gene if R² > 0.01 and expressed in > 5% of cells.
///
/// If `pca_loadings` (`n_genes × 2`) is provided, project cell velocities
/// to 2D PCA space; otherwise the embedding is zeros.
pub fn compute_rna_velocity(
    spliced: &Array2<f32>,   // n_cells × n_genes
    unspliced: &Array2<f32>, // n_cells × n_genes
    gene_ids: &[String],
    pca_loadings: Option<&Array2<f64>>, // n_genes × 2
) -> Result<VelocityResult> {
    let (n_cells, n_genes) = spliced.dim();
    if unspliced.dim() != (n_cells, n_genes) {
        anyhow::bail!(
            "spliced ({n_cells}×{n_genes}) and unspliced ({:?}) dimensions differ",
            unspliced.dim()
        );
    }
    if gene_ids.len() != n_genes {
        anyhow::bail!("gene_ids length {} != n_genes {n_genes}", gene_ids.len());
    }
    if n_cells == 0 || n_genes == 0 {
        return Ok(VelocityResult {
            gene_velocities: Vec::new(),
            velocity_embedding: Array2::zeros((n_cells, 2)),
            n_cells,
            n_velocity_genes: 0,
        });
    }

    let min_pct: f64 = 0.05;
    let min_r2: f64 = 0.01;

    // Per-gene full velocity matrix (n_cells × n_genes), filled for kept genes.
    // We accumulate the loadings dot product directly into the 2-D embedding.
    let mut gene_velocities: Vec<GeneVelocity> = Vec::new();

    // velocity_embedding[cell, dim] accumulated
    let mut vel_embed = Array2::<f64>::zeros((n_cells, 2));

    for g in 0..n_genes {
        let s: Vec<f64> = (0..n_cells).map(|i| spliced[[i, g]] as f64).collect();
        let u: Vec<f64> = (0..n_cells).map(|i| unspliced[[i, g]] as f64).collect();

        // Expression filter: > min_pct of cells must express the gene (s > 0)
        let n_expressed = s.iter().filter(|&&v| v > 0.0).count();
        if (n_expressed as f64 / n_cells as f64) < min_pct {
            continue;
        }

        // Estimate γ via upper-quantile regression:
        // take cells in top 5% of unspliced, use median(s/u) as γ.
        let gamma = estimate_gamma(&s, &u);

        // Compute velocity residuals
        let velocity: Vec<f64> = s
            .iter()
            .zip(u.iter())
            .map(|(&si, &ui)| si - gamma * ui)
            .collect();

        // R² of the steady-state fit (spliced ~ gamma * unspliced)
        let r2 = compute_r2(&s, &u, gamma);
        if r2 < min_r2 {
            continue;
        }

        // Project onto PCA loadings if available
        if let Some(loadings) = pca_loadings {
            if g < loadings.nrows() {
                let l0 = loadings[[g, 0]];
                let l1 = if loadings.ncols() > 1 {
                    loadings[[g, 1]]
                } else {
                    0.0
                };
                for i in 0..n_cells {
                    vel_embed[[i, 0]] += velocity[i] * l0;
                    vel_embed[[i, 1]] += velocity[i] * l1;
                }
            }
        }

        gene_velocities.push(GeneVelocity {
            gene_id: gene_ids[g].clone(),
            gamma,
            r_squared: r2,
            velocity,
        });
    }

    let n_velocity_genes = gene_velocities.len();

    Ok(VelocityResult {
        gene_velocities,
        velocity_embedding: vel_embed,
        n_cells,
        n_velocity_genes,
    })
}

/// Estimate γ as the median of `s_i / u_i` for cells in the top 5 % of
/// unspliced expression.  γ is clamped to [0.01, 100].
fn estimate_gamma(s: &[f64], u: &[f64]) -> f64 {
    let n = s.len();
    if n == 0 {
        return 1.0;
    }

    // Sort by unspliced expression to find the 95th-percentile threshold.
    let mut u_sorted = u.to_vec();
    u_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let q95_idx = ((n as f64 * 0.95) as usize).min(n.saturating_sub(1));
    let threshold = u_sorted[q95_idx];

    // Collect ratios for cells above the threshold.
    let mut ratios: Vec<f64> = s
        .iter()
        .zip(u.iter())
        .filter(|(_, &ui)| ui >= threshold && ui > 1e-10)
        .map(|(&si, &ui)| si / ui)
        .collect();

    if ratios.is_empty() {
        // Fall back: use all cells with u > 0.
        ratios = s
            .iter()
            .zip(u.iter())
            .filter(|(_, &ui)| ui > 1e-10)
            .map(|(&si, &ui)| si / ui)
            .collect();
    }

    if ratios.is_empty() {
        return 1.0;
    }

    // Median of ratios.
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = ratios.len() / 2;
    let gamma = if ratios.len() % 2 == 1 {
        ratios[mid]
    } else {
        (ratios[mid - 1] + ratios[mid]) / 2.0
    };

    gamma.clamp(0.01, 100.0)
}

/// R² of the linear fit `s = γ · u` (no intercept, as in the steady-state model).
fn compute_r2(s: &[f64], u: &[f64], gamma: f64) -> f64 {
    let n = s.len();
    if n == 0 {
        return 0.0;
    }
    let s_mean = s.iter().sum::<f64>() / n as f64;
    let ss_tot: f64 = s.iter().map(|&si| (si - s_mean).powi(2)).sum();
    if ss_tot < 1e-30 {
        return 0.0;
    }
    let ss_res: f64 = s
        .iter()
        .zip(u.iter())
        .map(|(&si, &ui)| (si - gamma * ui).powi(2))
        .sum();
    (1.0 - ss_res / ss_tot).max(0.0)
}

/// Compute cosine similarity between velocity vectors and displacement vectors
/// for each cell pair.
///
/// For cell *i* and cell *j*:
/// - displacement: `umap[j] - umap[i]`
/// - velocity: `velocity_embedding[i]`
/// - similarity: cosine(velocity_i, displacement_ij)
///
/// Returns an `n_cells × n_cells` cosine similarity matrix.
/// Only the `n_neighbors` nearest UMAP neighbours of each cell are evaluated;
/// all other entries remain 0.
pub fn velocity_graph(
    velocity_embedding: &Array2<f64>,
    umap_embedding: &Array2<f64>,
    n_neighbors: usize,
) -> Array2<f64> {
    let n_cells = velocity_embedding.nrows();
    let umap_dim = umap_embedding.ncols();
    let vel_dim = velocity_embedding.ncols();
    let mut out = Array2::<f64>::zeros((n_cells, n_cells));

    if n_cells == 0 || n_neighbors == 0 {
        return out;
    }

    let k = n_neighbors.min(n_cells.saturating_sub(1));

    for i in 0..n_cells {
        // Find k nearest neighbours in UMAP space (brute-force).
        let mut dists: Vec<(usize, f64)> = (0..n_cells)
            .filter(|&j| j != i)
            .map(|j| {
                let d2: f64 = (0..umap_dim)
                    .map(|d| {
                        let diff = umap_embedding[[i, d]] - umap_embedding[[j, d]];
                        diff * diff
                    })
                    .sum();
                (j, d2)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let vel_i: Vec<f64> = (0..vel_dim).map(|d| velocity_embedding[[i, d]]).collect();
        let vel_norm: f64 = vel_i.iter().map(|v| v * v).sum::<f64>().sqrt();

        for &(j, _) in dists.iter().take(k) {
            // Displacement vector: umap[j] - umap[i]
            let disp: Vec<f64> = (0..umap_dim.min(vel_dim))
                .map(|d| umap_embedding[[j, d]] - umap_embedding[[i, d]])
                .collect();
            let disp_norm: f64 = disp.iter().map(|v| v * v).sum::<f64>().sqrt();

            if vel_norm < 1e-14 || disp_norm < 1e-14 {
                continue;
            }
            let dot: f64 = vel_i.iter().zip(disp.iter()).map(|(a, b)| a * b).sum();
            out[[i, j]] = dot / (vel_norm * disp_norm);
        }
    }

    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    /// Build simple spliced/unspliced matrices where spliced > unspliced
    /// so γ should come out > 1.
    #[test]
    fn velocity_gamma_positive() {
        let n_cells = 20;
        let n_genes = 3;
        // spliced = 2 * unspliced  →  γ ≈ 2
        let spliced =
            Array2::<f32>::from_shape_fn((n_cells, n_genes), |(i, _)| (i as f32 + 1.0) * 2.0);
        let unspliced = Array2::<f32>::from_shape_fn((n_cells, n_genes), |(i, _)| i as f32 + 1.0);
        let gene_ids: Vec<String> = (0..n_genes).map(|i| format!("g{i}")).collect();

        let result = compute_rna_velocity(&spliced, &unspliced, &gene_ids, None).unwrap();
        assert!(!result.gene_velocities.is_empty(), "no velocity genes kept");
        for gv in &result.gene_velocities {
            assert!(
                gv.gamma > 1.0,
                "expected gamma > 1 for gene {}, got {}",
                gv.gene_id,
                gv.gamma
            );
        }
    }

    /// Output velocity_embedding should have shape n_cells × 2.
    #[test]
    fn velocity_result_shape() {
        let n_cells = 15;
        let n_genes = 5;
        let spliced =
            Array2::<f32>::from_shape_fn((n_cells, n_genes), |(i, g)| ((i + 1) * (g + 2)) as f32);
        let unspliced =
            Array2::<f32>::from_shape_fn((n_cells, n_genes), |(i, g)| ((i + 1) * (g + 1)) as f32);
        let gene_ids: Vec<String> = (0..n_genes).map(|i| format!("g{i}")).collect();

        // Provide identity-like PCA loadings (n_genes × 2)
        let mut loadings = Array2::<f64>::zeros((n_genes, 2));
        for g in 0..n_genes.min(2) {
            loadings[[g, g]] = 1.0;
        }

        let result =
            compute_rna_velocity(&spliced, &unspliced, &gene_ids, Some(&loadings)).unwrap();
        assert_eq!(
            result.velocity_embedding.dim(),
            (n_cells, 2),
            "embedding shape mismatch"
        );
        assert_eq!(result.n_cells, n_cells);
    }

    /// velocity_graph should return an n_cells × n_cells matrix.
    #[test]
    fn velocity_graph_shape() {
        let n_cells = 10;
        let vel = Array2::<f64>::from_shape_fn((n_cells, 2), |(i, j)| (i + j) as f64 * 0.1);
        let umap = Array2::<f64>::from_shape_fn((n_cells, 2), |(i, j)| (i * 2 + j) as f64);
        let graph = velocity_graph(&vel, &umap, 3);
        assert_eq!(graph.dim(), (n_cells, n_cells), "graph shape mismatch");
    }
}
