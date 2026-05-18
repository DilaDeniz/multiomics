//! Normalization routines for scRNA-seq count data.
//!
//! Implements a simplified pooling normalization inspired by scran
//! (Lun, Bach & Marioni 2016, Genome Biology 17:75) and standard
//! log-normalization.

use anyhow::Result;
use ndarray::Array2;

use crate::io::mex::CsrMatrix;

/// Compute per-cell size factors using scran-inspired pooling (Lun et al. 2016).
///
/// Cells are sorted by library size; sliding windows of each requested pool
/// size are used to construct sum-equations, which are solved for individual
/// size factors via 5 rounds of iterative refinement.
pub fn scran_size_factors(matrix: &CsrMatrix, pool_sizes: &[usize]) -> Result<Vec<f64>> {
    let n = matrix.n_cols;
    anyhow::ensure!(n >= 2, "need at least 2 cells to compute size factors");

    // --- Step 1: library sizes and sort order ---
    let lib_sizes: Vec<u64> = (0..n)
        .map(|j| {
            let mut s = 0u64;
            for row in 0..matrix.n_rows {
                for (col, val) in matrix.row(row) {
                    if col as usize == j {
                        s += val as u64;
                    }
                }
            }
            s
        })
        .collect();

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| lib_sizes[i]);

    // Compute pseudo-reference: geometric mean of library sizes (log-space mean)
    let log_ref: f64 = lib_sizes
        .iter()
        .map(|&l| (l as f64 + 1.0).ln())
        .sum::<f64>()
        / n as f64;
    let geo_ref = log_ref.exp();

    // Per-cell reference ratio (library size / geo_mean)
    let ref_ratio: Vec<f64> = lib_sizes.iter().map(|&l| l as f64 / geo_ref).collect();

    // --- Step 2: build pool equations ---
    // For each pool_size p and each starting cell offset o (circular):
    //   pool_sf = median-ratio of pooled library / pooled reference
    // Equation: sum of s_i for cells in window = pool_sf

    let mut pool_membership: Vec<Vec<usize>> = Vec::new(); // which cells (sorted order) in each pool
    let mut pool_sf: Vec<f64> = Vec::new(); // target sum for that pool

    for &p in pool_sizes {
        let p = p.min(n);
        for o in 0..n {
            let window: Vec<usize> = (0..p).map(|k| order[(o + k) % n]).collect();
            let pooled_lib: f64 = window.iter().map(|&i| lib_sizes[i] as f64).sum();
            let pooled_ref: f64 = window.iter().map(|&i| ref_ratio[i]).sum();
            let sf = if pooled_ref > 0.0 {
                pooled_lib / (geo_ref * pooled_ref)
            } else {
                1.0
            };
            pool_membership.push(window);
            pool_sf.push(sf);
        }
    }

    // --- Step 3: iterative refinement ---
    // Initialize size factors as average pool_sf for pools containing each cell
    let mut sf = vec![0.0f64; n];
    let mut cell_pool_count = vec![0usize; n];

    for (pm, &psf) in pool_membership.iter().zip(pool_sf.iter()) {
        for &i in pm {
            sf[i] += psf;
            cell_pool_count[i] += 1;
        }
    }
    for i in 0..n {
        sf[i] = if cell_pool_count[i] > 0 {
            sf[i] / cell_pool_count[i] as f64
        } else {
            1.0
        };
    }

    // 5 iterations of refinement: s_i = mean over pools containing i of
    // (pool_sf / sum_{k in pool, k != i} s_k)
    for _iter in 0..5 {
        let mut new_sf = vec![0.0f64; n];
        let mut counts = vec![0usize; n];
        for (pm, &psf) in pool_membership.iter().zip(pool_sf.iter()) {
            let pool_sum: f64 = pm.iter().map(|&i| sf[i]).sum();
            for &i in pm {
                let rest = pool_sum - sf[i];
                if rest > 1e-12 {
                    new_sf[i] += psf / rest;
                    counts[i] += 1;
                }
            }
        }
        for i in 0..n {
            sf[i] = if counts[i] > 0 {
                new_sf[i] / counts[i] as f64
            } else {
                sf[i]
            };
        }
    }

    // --- Step 4: floor non-positive size factors ---
    for s in &mut sf {
        if *s < 0.01 {
            *s = 0.01;
        }
    }

    anyhow::ensure!(
        sf.len() == n,
        "size factor count mismatch: expected {} got {}",
        n,
        sf.len()
    );
    Ok(sf)
}

/// Log-normalize counts using pre-computed size factors.
///
/// Returns a `[n_cells × n_genes]` `f32` matrix where each entry is
/// `log1p(count / size_factor × 10_000)` — the standard Seurat logNormalize.
pub fn log_normalize(matrix: &CsrMatrix, size_factors: &[f64]) -> Result<Array2<f32>> {
    anyhow::ensure!(
        size_factors.len() == matrix.n_cols,
        "size_factors length {} != n_cols {}",
        size_factors.len(),
        matrix.n_cols
    );

    let mut out = Array2::<f32>::zeros((matrix.n_cols, matrix.n_rows));

    for row in 0..matrix.n_rows {
        for (col, val) in matrix.row(row) {
            let j = col as usize;
            let sf = size_factors[j];
            let norm = (val as f64 / sf * 10_000.0).ln_1p() as f32;
            out[[j, row]] = norm;
        }
    }

    if !out.iter().all(|v| v.is_finite()) {
        anyhow::bail!("log_normalize produced non-finite values");
    }

    Ok(out)
}
