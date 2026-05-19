//! CITE-seq multimodal analysis with Weighted Nearest Neighbor (WNN) integration.
//!
//! Implements the WNN algorithm from Hao et al. 2021 (Seurat v4):
//! per-cell modality weights derived from within-modality prediction scores,
//! weighted KNN graph construction, Leiden clustering, and UMAP embedding.
//!
//! # References
//! - Hao Y, et al. (2021) Cell 184:3573–3587 (Seurat v4 WNN)

use anyhow::Result;
use ndarray::Array2;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Result of a WNN integration of two modalities.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WnnResult {
    /// Combined embedding (n_cells × 2) from weighted KNN graph.
    pub embedding: Array2<f64>,
    /// Per-cell RNA modality weight (complement is protein weight).
    pub rna_weights: Vec<f64>,
    /// Per-cell cluster labels from Leiden on WNN graph.
    pub clusters: Vec<u32>,
    pub n_cells: usize,
}

/// A CITE-seq dataset with RNA and ADT (protein) modalities.
pub struct CiteseqData {
    /// RNA count matrix (n_cells × n_rna_genes), normalized log counts.
    pub rna: Array2<f32>,
    pub rna_genes: Vec<String>,
    /// ADT count matrix (n_cells × n_proteins), CLR-normalized.
    pub adt: Array2<f32>,
    pub adt_names: Vec<String>,
    pub barcodes: Vec<String>,
}

// ---------------------------------------------------------------------------
// CLR normalization
// ---------------------------------------------------------------------------

/// Centered log-ratio normalization for ADT counts.
///
/// CLR(x_i) = log(x_i + 1) − log(geometric_mean_i) for each cell independently.
/// This is the standard normalization for CITE-seq protein data.
///
/// For each row (cell): geometric_mean = exp(mean(log(x + 1)))
/// CLR[i,j] = log(adt[i,j] + 1) - log(geometric_mean_i)
pub fn clr_normalize_adt(adt: &Array2<f32>) -> Array2<f32> {
    let (n_cells, n_proteins) = adt.dim();
    let mut out = Array2::<f32>::zeros((n_cells, n_proteins));

    for i in 0..n_cells {
        if n_proteins == 0 {
            continue;
        }
        // Geometric mean in log space: mean of log(x+1)
        let log_mean: f32 = (0..n_proteins)
            .map(|j| (adt[[i, j]] + 1.0_f32).ln())
            .sum::<f32>()
            / n_proteins as f32;

        for j in 0..n_proteins {
            out[[i, j]] = (adt[[i, j]] + 1.0_f32).ln() - log_mean;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Parse ADT matrix
// ---------------------------------------------------------------------------

/// Parse ADT count matrix from a CSV/TSV file.
///
/// Format: rows = cells (barcodes), cols = ADT names, with a header row.
/// Auto-detects delimiter (tab vs comma) from the first line.
///
/// Returns `(barcodes, adt_names, matrix)`.
pub fn parse_adt_matrix(data: &[u8]) -> Result<(Vec<String>, Vec<String>, Array2<f32>)> {
    let text = std::str::from_utf8(data)
        .map_err(|e| anyhow::anyhow!("ADT file is not valid UTF-8: {e}"))?;

    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("ADT file is empty"))?;

    // Auto-detect delimiter
    let delim = if header.contains('\t') { '\t' } else { ',' };

    let adt_names: Vec<String> = header
        .splitn(usize::MAX, delim)
        .skip(1) // first col is barcode
        .map(|s| s.trim().to_owned())
        .collect();

    let n_proteins = adt_names.len();
    if n_proteins == 0 {
        anyhow::bail!("ADT file has no protein columns");
    }

    let mut barcodes: Vec<String> = Vec::new();
    let mut values: Vec<f32> = Vec::new();

    for (line_no, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let mut cols = line.splitn(usize::MAX, delim);
        let barcode = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("line {} has no barcode column", line_no + 2))?
            .trim()
            .to_owned();
        barcodes.push(barcode);

        for (col_idx, token) in cols.enumerate() {
            if col_idx >= n_proteins {
                break;
            }
            let v: f32 = token.trim().parse().map_err(|_| {
                anyhow::anyhow!("line {}: invalid float '{}'", line_no + 2, token.trim())
            })?;
            values.push(v);
        }
        // Pad with zeros if row is shorter than expected
        let filled = values.len() - (barcodes.len() - 1) * n_proteins;
        values.extend(std::iter::repeat_n(0.0_f32, n_proteins - filled));
    }

    let n_cells = barcodes.len();
    if n_cells == 0 {
        anyhow::bail!("ADT file has no data rows");
    }

    let matrix = Array2::from_shape_vec((n_cells, n_proteins), values)
        .map_err(|e| anyhow::anyhow!("failed to build ADT matrix: {e}"))?;

    Ok((barcodes, adt_names, matrix))
}

// ---------------------------------------------------------------------------
// WNN integration
// ---------------------------------------------------------------------------

/// Compute WNN integration of RNA and protein modalities.
///
/// Algorithm (Hao et al. 2021):
/// 1. Compute KNN graphs separately for RNA (k=n_neighbors) and ADT in PCA space.
/// 2. For each cell i, compute within-modality prediction score:
///    - s_RNA(i) = mean cosine similarity of i's RNA neighbors to their own RNA neighbors
///    - s_ADT(i) = mean cosine similarity of i's ADT neighbors to their own ADT neighbors
/// 3. Per-cell modality weight: w_RNA(i) = s_RNA(i) / (s_RNA(i) + s_ADT(i))
/// 4. Build weighted KNN graph: for each cell, pool RNA and ADT neighbors,
///    weight each edge by modality weight × cosine similarity.
/// 5. Run Leiden clustering on weighted graph.
/// 6. UMAP embedding on weighted graph.
pub fn run_wnn(
    data: &CiteseqData,
    n_neighbors: usize,
    n_pca_rna: usize,
    n_pca_adt: usize,
    seed: u64,
) -> Result<WnnResult> {
    let n_cells = data.rna.nrows();
    if n_cells == 0 {
        anyhow::bail!("CiteseqData has zero cells");
    }
    if data.adt.nrows() != n_cells {
        anyhow::bail!(
            "RNA and ADT matrices have different cell counts: {} vs {}",
            n_cells,
            data.adt.nrows()
        );
    }

    let k = n_neighbors.min(n_cells.saturating_sub(1)).max(1);

    // --- Step 1: PCA for each modality ---
    let rna_f64 = f32_to_f64(&data.rna);
    let adt_clr = clr_normalize_adt(&data.adt);
    let adt_f64 = f32_to_f64(&adt_clr);

    let k_rna = n_pca_rna.min(data.rna.ncols()).min(n_cells);
    let k_adt = n_pca_adt.min(data.adt.ncols()).min(n_cells);

    let (rna_pca, _) = truncated_svd(&rna_f64, k_rna);
    let (adt_pca, _) = truncated_svd(&adt_f64, k_adt);

    // --- Step 2: KNN graphs in PCA space ---
    let rna_knn = knn_brute(&rna_pca, k); // Vec<Vec<usize>>
    let adt_knn = knn_brute(&adt_pca, k);

    // --- Step 3: Per-cell modality weights ---
    // s_RNA(i) = mean cosine-sim of cell i's RNA neighbors to each other's RNA neighbors
    let rna_normed = row_normalize(&rna_pca);
    let adt_normed = row_normalize(&adt_pca);

    let s_rna = prediction_scores(&rna_normed, &rna_knn);
    let s_adt = prediction_scores(&adt_normed, &adt_knn);

    let rna_weights: Vec<f64> = (0..n_cells)
        .map(|i| {
            let sr = s_rna[i].max(0.0);
            let sa = s_adt[i].max(0.0);
            let denom = sr + sa;
            if denom < 1e-15 {
                0.5
            } else {
                sr / denom
            }
        })
        .collect();

    // --- Step 4: Build weighted KNN graph ---
    // For each cell, pool RNA and ADT neighbors; weight by modality weight × cosine-sim.
    // Build as adjacency list with weights for Leiden and UMAP.
    let wnn_adj = build_wnn_graph(&rna_normed, &adt_normed, &rna_knn, &adt_knn, &rna_weights);

    // --- Step 5: Leiden clustering ---
    let knn_graph = wnn_adj_to_knn_graph(&wnn_adj, n_cells);
    let clusters = crate::clustering::leiden_cluster(&knn_graph, 1.0);

    // --- Step 6: UMAP on WNN graph ---
    // Build a synthetic PCA-like coordinate by concatenating scaled RNA+ADT PCA scores,
    // then run UMAP on the weighted adjacency via run_umap.
    // Per spec: call run_umap from umap module.
    let combined_pca = concat_pca(&rna_pca, &adt_pca, &rna_weights);
    let n_umap_neighbors = k.min(n_cells.saturating_sub(1)).max(2);
    let umap_result = crate::umap::run_umap(&combined_pca, n_umap_neighbors, 100, 0.1, 1.0, seed)?;

    Ok(WnnResult {
        embedding: umap_result.embedding,
        rna_weights,
        clusters,
        n_cells,
    })
}

// ---------------------------------------------------------------------------
// Truncated SVD (power iteration)
// ---------------------------------------------------------------------------

/// Truncated SVD via power iteration.
///
/// Returns `(scores: n_cells × k, singular_values: Vec<f64>)`.
/// Uses power iteration: start with random init, iterate X X^T v, 3 iterations.
fn truncated_svd(x: &Array2<f64>, k: usize) -> (Array2<f64>, Vec<f64>) {
    let (n_cells, n_features) = x.dim();
    if k == 0 || n_cells == 0 || n_features == 0 {
        return (Array2::zeros((n_cells, 0)), Vec::new());
    }
    let k = k.min(n_cells).min(n_features);

    // Center each column
    let mut centered = x.to_owned();
    for j in 0..n_features {
        let mean = centered.column(j).sum() / n_cells as f64;
        for i in 0..n_cells {
            centered[[i, j]] -= mean;
        }
    }

    // Power iteration for top-k singular vectors via deflation
    let mut scores = Array2::<f64>::zeros((n_cells, k));
    let mut svals = vec![0.0f64; k];
    let mut residual = centered.clone();

    for comp in 0..k {
        // Random init for right singular vector (n_features)
        let mut v = xorshift_vec_f64(
            n_features,
            (comp as u64 + 1).wrapping_mul(0x9e3779b97f4a7c15),
        );
        normalize_f64(&mut v);

        // 3 power iterations: u = X v / |X v|, v = X^T u / |X^T u|
        for _ in 0..3 {
            // u = residual @ v  (n_cells,)
            let mut u = vec![0.0f64; n_cells];
            #[allow(clippy::needless_range_loop)]
            for i in 0..n_cells {
                for j in 0..n_features {
                    u[i] += residual[[i, j]] * v[j];
                }
            }
            normalize_f64(&mut u);

            // v = residual^T @ u  (n_features,)
            let mut v_new = vec![0.0f64; n_features];
            #[allow(clippy::needless_range_loop)]
            for j in 0..n_features {
                for i in 0..n_cells {
                    v_new[j] += residual[[i, j]] * u[i];
                }
            }
            normalize_f64(&mut v_new);
            v = v_new;
        }

        // Final u and singular value
        let mut u = vec![0.0f64; n_cells];
        #[allow(clippy::needless_range_loop)]
        for i in 0..n_cells {
            for j in 0..n_features {
                u[i] += residual[[i, j]] * v[j];
            }
        }
        let sigma = l2_norm(&u);
        if sigma > 1e-15 {
            for val in u.iter_mut() {
                *val /= sigma;
            }
        }

        svals[comp] = sigma;
        for i in 0..n_cells {
            scores[[i, comp]] = u[i] * sigma;
        }

        // Deflate: residual -= sigma * u * v^T
        #[allow(clippy::needless_range_loop)]
        for i in 0..n_cells {
            for j in 0..n_features {
                residual[[i, j]] -= sigma * u[i] * v[j];
            }
        }
    }

    (scores, svals)
}

// ---------------------------------------------------------------------------
// KNN (brute force, for N < 5000)
// ---------------------------------------------------------------------------

/// Brute-force KNN using L2 distance. Returns top-k neighbor indices per cell.
fn knn_brute(pca: &Array2<f64>, k: usize) -> Vec<Vec<usize>> {
    let n_cells = pca.nrows();
    let k_actual = k.min(n_cells.saturating_sub(1));

    (0..n_cells)
        .map(|i| {
            let mut dists: Vec<(usize, f64)> = (0..n_cells)
                .filter(|&j| j != i)
                .map(|j| {
                    let d = (0..pca.ncols())
                        .map(|d| {
                            let diff = pca[[i, d]] - pca[[j, d]];
                            diff * diff
                        })
                        .sum::<f64>();
                    (j, d)
                })
                .collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            dists.truncate(k_actual);
            dists.into_iter().map(|(idx, _)| idx).collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Prediction scores (within-modality)
// ---------------------------------------------------------------------------

/// Compute per-cell within-modality prediction scores.
///
/// s(i) = mean cosine similarity between cell i and the neighbors of its neighbors
/// (the "second-order" neighborhood consistency).
fn prediction_scores(normed: &Array2<f64>, knn: &[Vec<usize>]) -> Vec<f64> {
    let n_cells = normed.nrows();
    let mut scores = vec![0.0f64; n_cells];

    for i in 0..n_cells {
        if knn[i].is_empty() {
            scores[i] = 0.0;
            continue;
        }
        let mut total = 0.0f64;
        let mut count = 0usize;
        for &nbr in &knn[i] {
            // cosine similarity between i and each of nbr's neighbors
            for &nbr2 in &knn[nbr] {
                if nbr2 != i {
                    total += cosine_sim_rows(normed, i, nbr2);
                    count += 1;
                }
            }
        }
        scores[i] = if count > 0 { total / count as f64 } else { 0.0 };
    }
    scores
}

/// Cosine similarity between row `i` and row `j` in a row-normalized matrix.
/// Since rows are already unit vectors, this is just the dot product.
#[inline]
fn cosine_sim_rows(normed: &Array2<f64>, i: usize, j: usize) -> f64 {
    let n = normed.ncols();
    let mut dot = 0.0f64;
    for d in 0..n {
        dot += normed[[i, d]] * normed[[j, d]];
    }
    dot
}

// ---------------------------------------------------------------------------
// WNN graph construction
// ---------------------------------------------------------------------------

/// Build the weighted nearest-neighbor adjacency list.
///
/// For each cell i, pools RNA and ADT neighbors, weights each edge by
/// modality_weight × cosine_similarity.
fn build_wnn_graph(
    rna_normed: &Array2<f64>,
    adt_normed: &Array2<f64>,
    rna_knn: &[Vec<usize>],
    adt_knn: &[Vec<usize>],
    rna_weights: &[f64],
) -> Vec<Vec<(usize, f64)>> {
    let n_cells = rna_normed.nrows();
    let mut adj: Vec<ahash::AHashMap<usize, f64>> =
        (0..n_cells).map(|_| ahash::AHashMap::new()).collect();

    for i in 0..n_cells {
        let w_rna = rna_weights[i];
        let w_adt = 1.0 - w_rna;

        // RNA neighbors
        for &j in &rna_knn[i] {
            let sim = cosine_sim_rows(rna_normed, i, j).max(0.0);
            let edge_w = w_rna * sim;
            *adj[i].entry(j).or_insert(0.0) += edge_w;
            *adj[j].entry(i).or_insert(0.0) += edge_w;
        }

        // ADT neighbors
        for &j in &adt_knn[i] {
            let sim = cosine_sim_rows(adt_normed, i, j).max(0.0);
            let edge_w = w_adt * sim;
            *adj[i].entry(j).or_insert(0.0) += edge_w;
            *adj[j].entry(i).or_insert(0.0) += edge_w;
        }
    }

    adj.into_iter()
        .map(|m| {
            let mut v: Vec<(usize, f64)> = m.into_iter().collect();
            v.sort_unstable_by_key(|&(j, _)| j);
            v
        })
        .collect()
}

/// Convert WNN adjacency list to a `KnnGraph` for Leiden clustering.
fn wnn_adj_to_knn_graph(adj: &[Vec<(usize, f64)>], n_cells: usize) -> crate::graph::KnnGraph {
    let k = adj.iter().map(|v| v.len()).max().unwrap_or(0);
    let neighbors: Vec<Vec<u32>> = adj
        .iter()
        .map(|nbrs| nbrs.iter().map(|&(j, _)| j as u32).collect())
        .collect();
    crate::graph::KnnGraph {
        n_cells,
        k,
        neighbors,
    }
}

// ---------------------------------------------------------------------------
// Combined PCA for UMAP input
// ---------------------------------------------------------------------------

/// Concatenate RNA and ADT PCA scores, weighting each by per-cell modality weights.
///
/// Returns an `[n_cells, n_rna_pca + n_adt_pca]` matrix.
fn concat_pca(rna_pca: &Array2<f64>, adt_pca: &Array2<f64>, rna_weights: &[f64]) -> Array2<f64> {
    let n_cells = rna_pca.nrows();
    let n_rna = rna_pca.ncols();
    let n_adt = adt_pca.ncols();
    let n_total = n_rna + n_adt;

    let mut out = Array2::<f64>::zeros((n_cells, n_total));
    for i in 0..n_cells {
        let w_rna = rna_weights[i];
        let w_adt = 1.0 - w_rna;
        for j in 0..n_rna {
            out[[i, j]] = rna_pca[[i, j]] * w_rna;
        }
        for j in 0..n_adt {
            out[[i, n_rna + j]] = adt_pca[[i, j]] * w_adt;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Convert an `Array2<f32>` to `Array2<f64>`.
fn f32_to_f64(x: &Array2<f32>) -> Array2<f64> {
    x.mapv(|v| v as f64)
}

/// Row-normalize a matrix so each row has unit L2 norm.
fn row_normalize(x: &Array2<f64>) -> Array2<f64> {
    let (n_rows, n_cols) = x.dim();
    let mut out = x.to_owned();
    for i in 0..n_rows {
        let norm = (0..n_cols)
            .map(|j| out[[i, j]] * out[[i, j]])
            .sum::<f64>()
            .sqrt();
        if norm > 1e-15 {
            for j in 0..n_cols {
                out[[i, j]] /= norm;
            }
        }
    }
    out
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn normalize_f64(v: &mut [f64]) {
    let n = l2_norm(v);
    if n > 1e-15 {
        for x in v.iter_mut() {
            *x /= n;
        }
    }
}

/// Generate a pseudorandom f64 vector using xorshift.
fn xorshift_vec_f64(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed ^ 0x6c62272e07bb0142;
    (0..n)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            // Map to [-0.5, 0.5]
            (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    /// Build a synthetic CiteseqData with random values.
    fn synthetic_data(n_cells: usize, n_rna: usize, n_adt: usize, seed: u64) -> CiteseqData {
        let mut s = seed ^ 0xdeadbeef;
        let mut rng = move || -> f32 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            ((s >> 11) as f64 / (1u64 << 53) as f64) as f32
        };

        let rna = Array2::from_shape_fn((n_cells, n_rna), |_| rng() * 5.0);
        let adt = Array2::from_shape_fn((n_cells, n_adt), |_| (rng() * 100.0).floor());

        CiteseqData {
            rna,
            rna_genes: (0..n_rna).map(|i| format!("gene_{i}")).collect(),
            adt,
            adt_names: (0..n_adt).map(|i| format!("protein_{i}")).collect(),
            barcodes: (0..n_cells).map(|i| format!("cell_{i}")).collect(),
        }
    }

    #[test]
    fn clr_normalize_adt_basic() {
        // 2 cells × 2 proteins
        let adt = Array2::from_shape_vec((2, 2), vec![1.0f32, 3.0, 5.0, 7.0]).unwrap();
        let clr = clr_normalize_adt(&adt);
        assert_eq!(clr.nrows(), 2);
        assert_eq!(clr.ncols(), 2);

        // Each row should sum to approximately 0
        for i in 0..2 {
            let row_sum: f32 = clr.row(i).iter().sum();
            assert!(
                row_sum.abs() < 1e-5,
                "CLR row {i} does not sum to 0: {row_sum}"
            );
        }
    }

    #[test]
    fn wnn_output_shape() {
        let data = synthetic_data(20, 10, 5, 42);
        let result = run_wnn(&data, 5, 5, 3, 42).expect("WNN failed");
        assert_eq!(result.embedding.nrows(), 20, "embedding rows");
        assert_eq!(result.embedding.ncols(), 2, "embedding cols");
        assert_eq!(result.n_cells, 20, "n_cells");
        assert_eq!(result.rna_weights.len(), 20, "rna_weights length");
        assert_eq!(result.clusters.len(), 20, "clusters length");
    }

    #[test]
    fn rna_weights_bounded() {
        let data = synthetic_data(20, 10, 5, 7);
        let result = run_wnn(&data, 5, 5, 3, 7).expect("WNN failed");
        for (i, &w) in result.rna_weights.iter().enumerate() {
            assert!(w > 0.0 && w < 1.0, "rna_weight[{i}] = {w} not in (0, 1)");
        }
    }
}
