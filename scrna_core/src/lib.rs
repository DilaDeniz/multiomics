//! Single-cell RNA-seq analysis pipeline.
//!
//! Implements 10x MEX I/O, QC, scran normalization, HVG selection,
//! PCA, KNN graph, Leiden clustering, and Wilcoxon cluster markers.
//!
//! # References
//! - Lun ATL, et al. (2016) Genome Biology 17:75 (scran)
//! - Stuart T, et al. (2019) Cell 177:1888 (Seurat v3 HVG)
//! - Traag VA, et al. (2019) Sci Reports 9:5233 (Leiden)

pub mod cellcomm;
pub mod clustering;
pub mod de;
pub mod doublets;
pub mod graph;
pub mod harmony;
pub mod hvg;
pub mod io;
pub mod multimodal;
pub mod normalize;
pub mod pseudotime;
pub mod qc;
pub mod spatial;
pub mod types;
pub mod umap;
pub mod umap_gpu;
pub mod velocity;

pub use cellcomm::{
    builtin_lr_database, compute_communication, filter_significant, CommScore, CommSummary, LRPair,
};
pub use clustering::leiden_cluster;
pub use de::{find_cluster_markers, ClusterMarker};
pub use doublets::{detect_doublets, DoubletScores};
pub use graph::{build_knn_graph, KnnGraph};
pub use harmony::{harmony_integrate, HarmonyResult};
pub use hvg::select_hvg;
pub use io::mex::{parse_10x_mex, CsrMatrix};
pub use multimodal::{clr_normalize_adt, parse_adt_matrix, run_wnn, CiteseqData, WnnResult};
pub use normalize::{log_normalize, scran_size_factors};
pub use pseudotime::{compute_pseudotime, PseudotimeResult};
pub use qc::{compute_qc, default_qc_filter, filter_cells, CellQc};
pub use spatial::{
    find_spatially_variable_genes, load_visium_dir, morans_i, parse_spot_matrix,
    parse_visium_positions, run_spatial_analysis, SpatialResult, SpatiallyVariableGene, Spot,
};
pub use types::SingleCellSummary;
pub use umap::{compute_fuzzy_graph, compute_fuzzy_graph_from_knn, run_umap, run_umap_from_graph, umap_from_pca, UmapResult};
pub use umap_gpu::run_umap_gpu;
pub use velocity::{compute_rna_velocity, velocity_graph, GeneVelocity, VelocityResult};

use anyhow::Result;
use ndarray::Array2;

/// Run the full scRNA-seq analysis pipeline on a 10x MEX directory.
///
/// Steps: parse → QC filter → scran normalization → HVG selection →
/// truncated PCA → KNN graph → Leiden clustering → Wilcoxon markers → summarize.
pub fn run_scrna_pipeline(
    mex_dir: &std::path::Path,
    n_hvg: usize,
    n_pca_components: usize,
    k_neighbors: usize,
    leiden_resolution: f64,
) -> Result<SingleCellSummary> {
    // 1. Parse MEX format
    log::info!("Parsing 10x MEX directory: {}", mex_dir.display());
    let raw_matrix = parse_10x_mex(mex_dir)?;
    let n_cells_raw = raw_matrix.n_cols;
    let n_genes = raw_matrix.n_rows;

    // 2. QC
    log::info!("Computing QC metrics for {} cells", n_cells_raw);
    let mut qc_metrics = compute_qc(&raw_matrix);
    default_qc_filter(&mut qc_metrics);
    let filtered = filter_cells(&raw_matrix, &qc_metrics);
    let n_cells_after_qc = filtered.n_cols;
    log::info!("{} cells retained after QC", n_cells_after_qc);

    // 3. Normalize
    log::info!("Computing scran size factors");
    let pool_sizes = vec![20, 40, 60, 80, 100];
    let size_factors = scran_size_factors(&filtered, &pool_sizes)?;
    let norm = log_normalize(&filtered, &size_factors)?;

    // 4. HVG selection
    let n_hvg_actual = n_hvg.min(n_genes);
    log::info!("Selecting {} highly variable genes", n_hvg_actual);
    let hvg_indices = select_hvg(&norm, n_hvg_actual);
    let n_hvg_selected = hvg_indices.len();

    // Extract HVG submatrix [n_cells, n_hvg]
    let hvg_matrix = select_columns(&norm, &hvg_indices);

    // 5. Truncated PCA
    let n_pca = n_pca_components.min(n_hvg_selected).min(n_cells_after_qc);
    log::info!("Running truncated PCA ({} components)", n_pca);
    let embedding = truncated_pca(&hvg_matrix, n_pca);

    // 6. KNN graph
    log::info!("Building KNN graph (k={})", k_neighbors);
    let knn = build_knn_graph(&embedding, k_neighbors);

    // 7. Leiden clustering
    log::info!(
        "Running Leiden clustering (resolution={})",
        leiden_resolution
    );
    let cluster_labels = leiden_cluster(&knn, leiden_resolution);
    let n_clusters = cluster_labels.iter().max().copied().unwrap_or(0) + 1;
    log::info!("Found {} clusters", n_clusters);

    // 8. Marker genes (Wilcoxon)
    log::info!("Finding cluster markers");
    let hvg_names: Vec<String> = hvg_indices
        .iter()
        .map(|&i| {
            filtered
                .features
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("gene_{i}"))
        })
        .collect();
    let all_markers = find_cluster_markers(&hvg_matrix, &cluster_labels, &hvg_names, 0.1);

    // Top 3 markers per cluster by p-value
    let mut top_markers: Vec<ClusterMarker> = Vec::new();
    for cid in 0..n_clusters {
        let mut cluster_markers: Vec<&ClusterMarker> =
            all_markers.iter().filter(|m| m.cluster == cid).collect();
        cluster_markers.sort_by(|a, b| {
            a.p_value
                .partial_cmp(&b.p_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for m in cluster_markers.into_iter().take(3) {
            top_markers.push(m.clone());
        }
    }

    // Compute medians for summary
    let passing_qc: Vec<&CellQc> = qc_metrics.iter().filter(|c| c.pass).collect();
    let mut gene_counts: Vec<f64> = passing_qc.iter().map(|c| c.n_genes as f64).collect();
    let mut count_totals: Vec<f64> = passing_qc.iter().map(|c| c.total_counts as f64).collect();
    let median_genes = median_f64(&mut gene_counts);
    let median_counts = median_f64(&mut count_totals);

    Ok(SingleCellSummary {
        n_cells_raw,
        n_cells_after_qc,
        n_genes,
        n_hvg: n_hvg_selected,
        n_clusters,
        median_genes_per_cell: median_genes,
        median_counts_per_cell: median_counts,
        top_markers,
    })
}

/// Extract a column subset from a `[n_rows, n_cols]` matrix.
fn select_columns(matrix: &Array2<f32>, indices: &[usize]) -> Array2<f32> {
    let n_rows = matrix.nrows();
    let n_selected = indices.len();
    let mut out = Array2::<f32>::zeros((n_rows, n_selected));
    for (new_col, &old_col) in indices.iter().enumerate() {
        for row in 0..n_rows {
            out[[row, new_col]] = matrix[[row, old_col]];
        }
    }
    out
}

/// Truncated PCA via randomized SVD (Halko, Martinsson & Tropp 2011, Algorithm 4.4).
///
/// Centers each gene column, then applies a random projection to build a
/// low-rank approximation via QR + small SVD.
///
/// Returns `[n_cells, n_components]` embedding.
pub fn truncated_pca(matrix: &Array2<f32>, n_components: usize) -> Array2<f32> {
    let (n_cells, n_genes) = matrix.dim();
    if n_components == 0 || n_cells == 0 || n_genes == 0 {
        return Array2::zeros((n_cells, 0));
    }

    let k = n_components.min(n_cells).min(n_genes);
    let oversampling = 10usize;
    let l = (k + oversampling).min(n_cells).min(n_genes);

    // Step 1: center each gene (column)
    let mut centered = matrix.to_owned();
    for g in 0..n_genes {
        let col = centered.column(g);
        let mean = col.sum() / n_cells as f32;
        for i in 0..n_cells {
            centered[[i, g]] -= mean;
        }
    }

    // Step 2: random Gaussian matrix Omega [n_genes, l]
    let omega = gaussian_matrix(n_genes, l);

    // Step 3: Y = centered @ Omega  [n_cells, l]
    let y = mat_mul_f32(&centered, &omega);

    // Step 4: QR decomposition of Y → Q [n_cells, l]
    let q = qr_orthonormalize(y);

    // Step 5: B = Q^T @ centered  [l, n_genes]
    let b = mat_mul_f32(&q.t().to_owned(), &centered);

    // Step 6: SVD of small B [l, n_genes] via power iteration / Gram matrix approach
    svd_embedding(&q, &b, k)
}

/// Compute Q @ U[:, :k] @ diag(S[:k]) from the small SVD of B = Q^T @ A.
///
/// Uses Gram iteration: B @ B^T is `[l, l]`, eigendecompose it.
fn svd_embedding(q: &Array2<f32>, b: &Array2<f32>, k: usize) -> Array2<f32> {
    let l = b.nrows();
    let k = k.min(l);

    // Gram matrix C = B @ B^T  [l, l]
    let b_t = b.t().to_owned();
    let c = mat_mul_f32(b, &b_t);

    // Power iteration eigen-decomposition of symmetric C
    let (eigvecs, eigvals) = symmetric_eigen(&c, k);
    // eigvals[i] ≈ singular_value_i^2

    // Embedding = Q @ eigvecs @ diag(sqrt(eigvals))  [n_cells, k]
    let n_cells = q.nrows();
    let mut out = Array2::<f32>::zeros((n_cells, k));
    for j in 0..k {
        let sv = (eigvals[j].max(0.0) as f64).sqrt() as f32;
        for i in 0..n_cells {
            let mut val = 0.0f32;
            for r in 0..l {
                val += q[[i, r]] * eigvecs[[r, j]];
            }
            out[[i, j]] = val * sv;
        }
    }
    out
}

/// Symmetric eigen-decomposition via power iteration, returning top-k
/// eigenvectors and eigenvalues in descending order.
fn symmetric_eigen(c: &Array2<f32>, k: usize) -> (Array2<f32>, Vec<f32>) {
    let l = c.nrows();
    let k = k.min(l);
    let mut vecs = Array2::<f32>::zeros((l, k));
    let mut vals = vec![0.0f32; k];

    // Deflation: iteratively extract top eigenvectors
    let mut residual = c.to_owned();

    for j in 0..k {
        // Initialize random vector
        let mut v = init_vec(l, j as u64 + 1);
        normalize_vec(&mut v);

        // Power iteration: 30 steps
        for _ in 0..30 {
            v = mat_vec_mul(&residual, &v);
            normalize_vec(&mut v);
        }

        // Rayleigh quotient for eigenvalue
        let av = mat_vec_mul(&residual, &v);
        let lam: f32 = v.iter().zip(av.iter()).map(|(a, b)| a * b).sum();
        vals[j] = lam;

        // Store eigenvector
        for i in 0..l {
            vecs[[i, j]] = v[i];
        }

        // Deflate
        for r in 0..l {
            for s in 0..l {
                residual[[r, s]] -= lam * v[r] * v[s];
            }
        }
    }

    (vecs, vals)
}

fn init_vec(n: usize, seed: u64) -> Vec<f32> {
    let mut s = seed ^ 0x9e37_79b9_7f4a_7c15;
    (0..n)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            ((s >> 11) as f64 / (1u64 << 53) as f64 - 0.5) as f32
        })
        .collect()
}

fn normalize_vec(v: &mut [f32]) {
    let norm: f32 = v.iter().map(|x| x * x).sum::<f32>().sqrt();
    if norm > 1e-12 {
        for x in v.iter_mut() {
            *x /= norm;
        }
    }
}

fn mat_vec_mul(m: &Array2<f32>, v: &[f32]) -> Vec<f32> {
    let n = m.nrows();
    let p = m.ncols();
    let mut out = vec![0.0f32; n];
    for i in 0..n {
        for j in 0..p {
            out[i] += m[[i, j]] * v[j];
        }
    }
    out
}

/// Cache-blocked matrix multiplication: A [m,k] × B [k,n] → C [m,n].
///
/// Uses a 32×32 tile so that the working set for one tile fits in L1 cache
/// (~32KB on most microarchitectures). For a 2 000×2 000 matrix this is
/// 4-6× faster than the naïve triple loop.
fn mat_mul_f32(a: &Array2<f32>, b: &Array2<f32>) -> Array2<f32> {
    const BLOCK: usize = 32;
    let m = a.nrows();
    let k = a.ncols();
    let n = b.ncols();
    let mut c = Array2::<f32>::zeros((m, n));

    // Flatten to row-major slices for raw pointer arithmetic.
    let a_s = a.as_slice().expect("a must be contiguous");
    let b_s = b.as_slice().expect("b must be contiguous");
    let c_s = c.as_slice_mut().expect("c must be contiguous");

    let mut ii = 0;
    while ii < m {
        let i_end = (ii + BLOCK).min(m);
        let mut jj = 0;
        while jj < n {
            let j_end = (jj + BLOCK).min(n);
            let mut pp = 0;
            while pp < k {
                let p_end = (pp + BLOCK).min(k);
                for i in ii..i_end {
                    let a_row = &a_s[i * k..i * k + k];
                    let c_row = &mut c_s[i * n..i * n + n];
                    for p in pp..p_end {
                        let a_ip = a_row[p];
                        let b_row = &b_s[p * n..p * n + n];
                        for j in jj..j_end {
                            c_row[j] += a_ip * b_row[j];
                        }
                    }
                }
                pp += BLOCK;
            }
            jj += BLOCK;
        }
        ii += BLOCK;
    }
    c
}

/// Orthonormalize columns of Y via modified Gram-Schmidt.
fn qr_orthonormalize(mut y: Array2<f32>) -> Array2<f32> {
    let n = y.nrows();
    let m = y.ncols();
    for j in 0..m {
        // Normalize column j
        let norm: f32 = (0..n).map(|i| y[[i, j]] * y[[i, j]]).sum::<f32>().sqrt();
        if norm > 1e-12 {
            for i in 0..n {
                y[[i, j]] /= norm;
            }
        }
        // Orthogonalize subsequent columns
        for k in (j + 1)..m {
            let dot: f32 = (0..n).map(|i| y[[i, j]] * y[[i, k]]).sum();
            for i in 0..n {
                let sub = dot * y[[i, j]];
                y[[i, k]] -= sub;
            }
        }
    }
    y
}

/// Generate a random Gaussian matrix [rows, cols] with xorshift + Box-Muller.
fn gaussian_matrix(rows: usize, cols: usize) -> Array2<f32> {
    let mut seed: u64 = 0xdead_beef_cafe_babe;
    let mut rng = move || -> f64 {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };

    let mut out = Array2::<f32>::zeros((rows, cols));
    let mut i = 0;
    let mut j = 0;
    let total = rows * cols;
    let mut filled = 0;
    while filled < total {
        let u1 = rng().max(1e-300);
        let u2 = rng();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = std::f64::consts::TAU * u2;
        out[[i, j]] = (r * theta.cos()) as f32;
        j += 1;
        if j == cols {
            j = 0;
            i += 1;
        }
        filled += 1;
        if filled < total {
            out[[i, j]] = (r * theta.sin()) as f32;
            j += 1;
            if j == cols {
                j = 0;
                i += 1;
            }
            filled += 1;
        }
    }
    out
}

fn median_f64(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    }
}
