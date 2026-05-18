//! Per-cell QC metrics and adaptive MAD-based filtering.

use crate::io::mex::CsrMatrix;

/// QC metrics for a single cell.
#[derive(Debug, Clone)]
pub struct CellQc {
    /// Cell barcode string.
    pub barcode: String,
    /// Number of genes with at least one count.
    pub n_genes: u32,
    /// Total UMI counts.
    pub total_counts: u32,
    /// Fraction of counts mapping to mitochondrial genes (0–100).
    pub pct_mito: f32,
    /// Whether the cell passes QC filters.
    pub pass: bool,
}

/// Compute per-cell QC metrics from a count matrix.
///
/// Mitochondrial genes are identified by the `MT-` prefix (case-insensitive).
pub fn compute_qc(matrix: &CsrMatrix) -> Vec<CellQc> {
    let mito_flags: Vec<bool> = matrix
        .features
        .iter()
        .map(|f| f.to_ascii_uppercase().starts_with("MT-"))
        .collect();

    let mut n_genes = vec![0u32; matrix.n_cols];
    let mut total_counts = vec![0u32; matrix.n_cols];
    let mut mito_counts = vec![0u32; matrix.n_cols];

    for (row, &is_mito) in mito_flags.iter().enumerate() {
        for (col, val) in matrix.row(row) {
            let c = col as usize;
            n_genes[c] += 1;
            total_counts[c] += val;
            if is_mito {
                mito_counts[c] += val;
            }
        }
    }

    (0..matrix.n_cols)
        .map(|j| {
            let tot = total_counts[j];
            let pct_mito = if tot > 0 {
                mito_counts[j] as f32 / tot as f32 * 100.0
            } else {
                0.0
            };
            CellQc {
                barcode: matrix.barcodes[j].clone(),
                n_genes: n_genes[j],
                total_counts: tot,
                pct_mito,
                pass: true,
            }
        })
        .collect()
}

/// Apply adaptive MAD-based filtering to `qc`, setting `pass = false` for outliers.
///
/// Cells are flagged when:
/// - `n_genes < median − 3 × MAD` (low-complexity / empty droplets)
/// - `n_genes > median + 5 × MAD` (potential doublets)
/// - `pct_mito > 20.0` (high mitochondrial content)
pub fn default_qc_filter(qc: &mut [CellQc]) {
    let genes: Vec<f64> = qc.iter().map(|c| c.n_genes as f64).collect();
    let low = mad_threshold(&genes, 3.0);
    let high = mad_upper(&genes, 5.0);

    for cell in qc.iter_mut() {
        let g = cell.n_genes as f64;
        if g < low || g > high || cell.pct_mito > 20.0 {
            cell.pass = false;
        }
    }
}

/// Rebuild the matrix keeping only cells where `qc[j].pass == true`.
pub fn filter_cells(matrix: &CsrMatrix, qc: &[CellQc]) -> CsrMatrix {
    let keep: Vec<usize> = (0..matrix.n_cols).filter(|&j| qc[j].pass).collect();

    // Build a mapping old_col → new_col (u32::MAX means dropped)
    let mut col_map = vec![u32::MAX; matrix.n_cols];
    for (new, &old) in keep.iter().enumerate() {
        col_map[old] = new as u32;
    }

    let new_n_cols = keep.len();
    let mut indptr = vec![0u32; matrix.n_rows + 1];
    let mut indices: Vec<u32> = Vec::new();
    let mut data: Vec<u32> = Vec::new();

    for row in 0..matrix.n_rows {
        for (old_col, val) in matrix.row(row) {
            let new_col = col_map[old_col as usize];
            if new_col != u32::MAX {
                indices.push(new_col);
                data.push(val);
            }
        }
        indptr[row + 1] = indices.len() as u32;
    }

    let barcodes = keep.iter().map(|&j| matrix.barcodes[j].clone()).collect();

    CsrMatrix {
        n_rows: matrix.n_rows,
        n_cols: new_n_cols,
        indptr,
        indices,
        data,
        barcodes,
        features: matrix.features.clone(),
    }
}

/// Lower MAD threshold: `median − n_mad × MAD`.
fn mad_threshold(values: &[f64], n_mad: f64) -> f64 {
    let med = median_f64(values);
    let abs_devs: Vec<f64> = values.iter().map(|&x| (x - med).abs()).collect();
    let mad = median_f64(&abs_devs);
    med - n_mad * mad
}

/// Upper MAD threshold: `median + n_mad × MAD`.
fn mad_upper(values: &[f64], n_mad: f64) -> f64 {
    let med = median_f64(values);
    let abs_devs: Vec<f64> = values.iter().map(|&x| (x - med).abs()).collect();
    let mad = median_f64(&abs_devs);
    med + n_mad * mad
}

fn median_f64(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}
