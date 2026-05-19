//! Spatial transcriptomics analysis: Visium I/O, Moran's I, spatially variable genes.
//!
//! Supports 10x Visium output format (tissue_positions.csv, MEX matrix).
//! Implements Moran's I for spatial autocorrelation with a permutation test.

use ahash::AHashMap;
use anyhow::{Context, Result};
use ndarray::Array2;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// A spatial spot (Visium bead or Slide-seq bead).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Spot {
    pub barcode: String,
    /// Physical x coordinate (microns or array position).
    pub x: f64,
    /// Physical y coordinate.
    pub y: f64,
    /// Total UMI count.
    pub total_counts: u64,
    /// Leiden cluster label (set after clustering).
    pub cluster: Option<u32>,
}

/// A spatially variable gene result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatiallyVariableGene {
    pub gene_id: String,
    /// Moran's I statistic \[-1, 1\]. Values near +1 = strong spatial clustering.
    pub morans_i: f64,
    /// Pseudo p-value from permutation test (100 permutations).
    pub p_value: f64,
}

/// Full result of a spatial transcriptomics analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialResult {
    pub spots: Vec<Spot>,
    pub n_spots: usize,
    pub n_genes: usize,
    pub spatially_variable_genes: Vec<SpatiallyVariableGene>,
    pub mean_counts_per_spot: f64,
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

/// Parse a Visium-style spot positions file.
///
/// Format: CSV or TSV with columns:
/// `barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres`
/// (standard 10x Visium tissue_positions.csv format).
///
/// Only spots where `in_tissue == 1` are returned.
pub fn parse_visium_positions(data: &[u8]) -> Result<Vec<Spot>> {
    let text = std::str::from_utf8(data).context("tissue_positions is not valid UTF-8")?;
    let mut spots = Vec::new();

    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        // Detect separator (comma or tab)
        let sep = if line.contains('\t') { '\t' } else { ',' };
        let parts: Vec<&str> = line.splitn(7, sep).collect();

        // Skip header row (first field is "barcode" literally)
        if line_no == 0 && parts.first().map(|s| s.trim()) == Some("barcode") {
            continue;
        }

        if parts.len() < 6 {
            anyhow::bail!(
                "tissue_positions line {}: expected ≥6 fields, got {}",
                line_no + 1,
                parts.len()
            );
        }

        let barcode = parts[0].trim().to_owned();
        let in_tissue: u8 = parts[1]
            .trim()
            .parse()
            .with_context(|| format!("line {}: parsing in_tissue", line_no + 1))?;
        if in_tissue == 0 {
            continue;
        }

        // Use pixel coordinates (columns 4 and 5) as spatial coords.
        let pxl_row: f64 = parts[4]
            .trim()
            .parse()
            .with_context(|| format!("line {}: parsing pxl_row_in_fullres", line_no + 1))?;
        let pxl_col: f64 = parts[5]
            .trim()
            .parse()
            .with_context(|| format!("line {}: parsing pxl_col_in_fullres", line_no + 1))?;

        spots.push(Spot {
            barcode,
            x: pxl_col,
            y: pxl_row,
            total_counts: 0,
            cluster: None,
        });
    }

    Ok(spots)
}

/// Parse a simple spot × gene count matrix.
///
/// Format: rows = spots (barcodes), columns = genes.
/// First row = header (gene names), first column = barcode.
///
/// Returns `(barcodes, gene_ids, count_matrix [n_spots × n_genes])`.
pub fn parse_spot_matrix(data: &[u8]) -> Result<(Vec<String>, Vec<String>, Array2<f32>)> {
    let text = std::str::from_utf8(data).context("spot matrix is not valid UTF-8")?;
    let mut lines = text.lines().peekable();

    // Parse header row
    let header = lines.next().context("spot matrix is empty")?;
    let sep = if header.contains('\t') { '\t' } else { ',' };
    let gene_ids: Vec<String> = header
        .splitn(usize::MAX, sep)
        .skip(1) // skip barcode column
        .map(|s| s.trim().to_owned())
        .collect();
    let n_genes = gene_ids.len();

    let mut barcodes = Vec::new();
    let mut rows: Vec<Vec<f32>> = Vec::new();

    for (line_no, line) in lines.enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.splitn(n_genes + 2, sep).collect();
        if parts.is_empty() {
            continue;
        }
        let barcode = parts[0].trim().to_owned();
        barcodes.push(barcode);

        let mut row = Vec::with_capacity(n_genes);
        for (col_idx, &field) in parts.iter().skip(1).enumerate() {
            if col_idx >= n_genes {
                break;
            }
            let val: f32 = field
                .trim()
                .parse()
                .with_context(|| format!("line {}: column {}", line_no + 2, col_idx + 2))?;
            row.push(val);
        }
        // Pad with zeros if row is short
        while row.len() < n_genes {
            row.push(0.0);
        }
        rows.push(row);
    }

    let n_spots = barcodes.len();
    let mut matrix = Array2::<f32>::zeros((n_spots, n_genes));
    for (i, row) in rows.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            matrix[[i, j]] = val;
        }
    }

    Ok((barcodes, gene_ids, matrix))
}

/// Load spatial data from a 10x Visium output directory.
///
/// Expects:
/// - `tissue_positions.csv` (or `tissue_positions_list.csv`)
/// - `matrix.mtx[.gz]`, `barcodes.tsv[.gz]`, `features.tsv[.gz]`
///
/// Returns `(spots_in_tissue, csr_matrix)`.
pub fn load_visium_dir(dir: &std::path::Path) -> Result<(Vec<Spot>, crate::io::mex::CsrMatrix)> {
    // Try tissue_positions.csv first, then tissue_positions_list.csv
    let positions_path = {
        let p1 = dir.join("tissue_positions.csv");
        let p2 = dir.join("tissue_positions_list.csv");
        if p1.exists() {
            p1
        } else if p2.exists() {
            p2
        } else {
            anyhow::bail!(
                "No tissue_positions.csv or tissue_positions_list.csv found in {}",
                dir.display()
            );
        }
    };

    let positions_bytes = std::fs::read(&positions_path)
        .with_context(|| format!("reading {}", positions_path.display()))?;
    let mut spots = parse_visium_positions(&positions_bytes)?;

    // Build barcode → spot index map (owned keys to avoid lifetime issues)
    let barcode_to_spot: AHashMap<String, usize> = spots
        .iter()
        .enumerate()
        .map(|(i, spot)| (spot.barcode.clone(), i))
        .collect();

    // Parse MEX matrix
    let matrix = crate::io::mex::parse_10x_mex(dir)?;

    // Fill total_counts from the matrix for barcodes that are in-tissue spots
    for (col_idx, barcode) in matrix.barcodes.iter().enumerate() {
        if let Some(&spot_idx) = barcode_to_spot.get(barcode.as_str()) {
            let total: u64 = (0..matrix.n_rows)
                .flat_map(|row| matrix.row(row))
                .filter(|&(c, _)| c as usize == col_idx)
                .map(|(_, v)| v as u64)
                .sum();
            spots[spot_idx].total_counts = total;
        }
    }

    Ok((spots, matrix))
}

// ---------------------------------------------------------------------------
// Spatial weight matrix helpers
// ---------------------------------------------------------------------------

/// Build a k-nearest-neighbour weight matrix (row-normalized).
///
/// Returns `weights[i][j]` = normalised weight for pair (i, j) where j is a
/// neighbour of i, plus the total sum W = Σ wᵢⱼ.
fn build_knn_weights(spots: &[Spot], n_neighbors: usize) -> (Vec<Vec<(usize, f64)>>, f64) {
    let n = spots.len();
    if n == 0 {
        return (vec![], 0.0);
    }
    let k = n_neighbors.min(n.saturating_sub(1));

    // For each spot, find the k nearest neighbours by Euclidean distance.
    let mut weights: Vec<Vec<(usize, f64)>> = Vec::with_capacity(n);
    for i in 0..n {
        let xi = spots[i].x;
        let yi = spots[i].y;

        let mut dists: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| {
                let dx = xi - spots[j].x;
                let dy = yi - spots[j].y;
                (j, dx * dx + dy * dy)
            })
            .collect();

        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        dists.truncate(k);

        let w = if k > 0 { 1.0 / k as f64 } else { 0.0 };
        let row: Vec<(usize, f64)> = dists.into_iter().map(|(j, _)| (j, w)).collect();
        weights.push(row);
    }

    let total_w: f64 = weights
        .iter()
        .map(|row| row.iter().map(|(_, w)| w).sum::<f64>())
        .sum();
    (weights, total_w)
}

// ---------------------------------------------------------------------------
// Moran's I
// ---------------------------------------------------------------------------

/// Compute Moran's I for a gene's expression across spots.
///
/// Moran's I = (N / W) × (Σᵢ Σⱼ wᵢⱼ(xᵢ − x̄)(xⱼ − x̄)) / (Σᵢ(xᵢ − x̄)²)
///
/// Spatial weight matrix W: wᵢⱼ = 1 if spot j is among the `n_neighbors`
/// nearest neighbours of i, else 0; W is row-normalised.
///
/// A permutation test (`n_permutations` shuffles) yields the p-value.
pub fn morans_i(
    expression: &[f64],
    spots: &[Spot],
    n_neighbors: usize,
    n_permutations: usize,
    seed: u64,
) -> (f64, f64) {
    let n = expression.len();
    if n < 3 {
        return (0.0, 1.0);
    }

    let (weights, total_w) = build_knn_weights(spots, n_neighbors);
    if total_w < 1e-12 {
        return (0.0, 1.0);
    }

    let observed = compute_morans_i_stat(expression, &weights, total_w);

    // Permutation test
    let mut perm_expr = expression.to_vec();
    let mut count_extreme = 0usize;
    let mut rng_state = seed ^ 0x9e37_79b9_7f4a_7c15u64;

    for _ in 0..n_permutations {
        // Fisher-Yates shuffle with xorshift64 RNG
        for i in (1..n).rev() {
            rng_state = xorshift64(rng_state);
            let j = (rng_state as usize) % (i + 1);
            perm_expr.swap(i, j);
        }
        let perm_i = compute_morans_i_stat(&perm_expr, &weights, total_w);
        if perm_i >= observed {
            count_extreme += 1;
        }
    }

    let p_value = (count_extreme + 1) as f64 / (n_permutations + 1) as f64;
    (observed, p_value)
}

fn xorshift64(mut state: u64) -> u64 {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    state
}

fn compute_morans_i_stat(expression: &[f64], weights: &[Vec<(usize, f64)>], total_w: f64) -> f64 {
    let n = expression.len() as f64;
    let mean = expression.iter().sum::<f64>() / n;

    let deviations: Vec<f64> = expression.iter().map(|&x| x - mean).collect();

    let numerator: f64 = weights
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .map(|&(j, w)| w * deviations[i] * deviations[j])
                .sum::<f64>()
        })
        .sum();

    let denominator: f64 = deviations.iter().map(|&d| d * d).sum::<f64>();

    if denominator < 1e-12 {
        return 0.0;
    }

    (n / total_w) * (numerator / denominator)
}

// ---------------------------------------------------------------------------
// Spatially variable gene detection
// ---------------------------------------------------------------------------

/// Find spatially variable genes using Moran's I.
///
/// Returns the top `n_top` genes sorted by Moran's I descending.
///
/// Steps:
/// 1. Compute log1p-normalised expression for each gene.
/// 2. Filter to genes expressed in > 10% of spots.
/// 3. Compute Moran's I for each gene in parallel (rayon).
pub fn find_spatially_variable_genes(
    count_matrix: &Array2<f32>,
    gene_ids: &[String],
    spots: &[Spot],
    n_top: usize,
    n_neighbors: usize,
) -> Vec<SpatiallyVariableGene> {
    let (n_spots, n_genes) = count_matrix.dim();
    if n_spots == 0 || n_genes == 0 {
        return Vec::new();
    }

    let min_spots = (n_spots as f64 * 0.10).ceil() as usize;

    // Pre-build weight matrix (shared across genes)
    let (weights, total_w) = build_knn_weights(spots, n_neighbors);
    let weights = std::sync::Arc::new(weights);

    let mut results: Vec<SpatiallyVariableGene> = (0..n_genes)
        .into_par_iter()
        .filter_map(|g| {
            let gene_id = gene_ids.get(g)?.clone();

            // Build log1p-normalised expression vector
            let expr: Vec<f64> = (0..n_spots)
                .map(|s| (count_matrix[[s, g]] as f64 + 1.0).ln())
                .collect();

            // Filter: must be expressed in > min_spots spots
            let n_expressed = expr.iter().filter(|&&v| v > 0.0).count();
            if n_expressed <= min_spots {
                return None;
            }

            let (mi, pv) =
                compute_morans_i_with_permutations(&expr, &weights, total_w, 100, g as u64 + 42);

            Some(SpatiallyVariableGene {
                gene_id,
                morans_i: mi,
                p_value: pv,
            })
        })
        .collect();

    // Sort by Moran's I descending
    results.sort_by(|a, b| {
        b.morans_i
            .partial_cmp(&a.morans_i)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    results.truncate(n_top);
    results
}

/// Internal: compute Moran's I with permutation test using pre-built weights.
fn compute_morans_i_with_permutations(
    expression: &[f64],
    weights: &[Vec<(usize, f64)>],
    total_w: f64,
    n_permutations: usize,
    seed: u64,
) -> (f64, f64) {
    let n = expression.len();
    if n < 3 || total_w < 1e-12 {
        return (0.0, 1.0);
    }

    let observed = compute_morans_i_stat(expression, weights, total_w);

    let mut perm_expr = expression.to_vec();
    let mut count_extreme = 0usize;
    let mut rng_state = seed ^ 0x9e37_79b9_7f4a_7c15u64;

    for _ in 0..n_permutations {
        for i in (1..n).rev() {
            rng_state = xorshift64(rng_state);
            let j = (rng_state as usize) % (i + 1);
            perm_expr.swap(i, j);
        }
        let perm_i = compute_morans_i_stat(&perm_expr, weights, total_w);
        if perm_i >= observed {
            count_extreme += 1;
        }
    }

    let p_value = (count_extreme + 1) as f64 / (n_permutations + 1) as f64;
    (observed, p_value)
}

// ---------------------------------------------------------------------------
// Full pipeline
// ---------------------------------------------------------------------------

/// Run the full spatial transcriptomics pipeline.
pub fn run_spatial_analysis(
    spots: &mut [Spot],
    count_matrix: &Array2<f32>,
    gene_ids: &[String],
    n_top_svg: usize,
) -> Result<SpatialResult> {
    let n_spots = spots.len();
    let n_genes = gene_ids.len();

    log::info!(
        "Running spatial analysis: {} spots, {} genes",
        n_spots,
        n_genes
    );

    if n_spots == 0 {
        anyhow::bail!("No spots provided for spatial analysis");
    }

    // Fill total_counts from count_matrix if not already set
    for (s, spot) in spots.iter_mut().enumerate() {
        if spot.total_counts == 0 && s < count_matrix.nrows() {
            spot.total_counts = count_matrix.row(s).iter().map(|&v| v as u64).sum();
        }
    }

    let mean_counts_per_spot =
        spots.iter().map(|s| s.total_counts as f64).sum::<f64>() / n_spots as f64;

    log::info!("Finding spatially variable genes (top {})", n_top_svg);
    let spatially_variable_genes =
        find_spatially_variable_genes(count_matrix, gene_ids, spots, n_top_svg, 6);

    Ok(SpatialResult {
        spots: spots.to_owned(),
        n_spots,
        n_genes,
        spatially_variable_genes,
        mean_counts_per_spot,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// 10 spots in a line with expression increasing left→right → Moran's I > 0.3
    #[test]
    fn morans_i_positive_clustering() {
        let spots: Vec<Spot> = (0..10)
            .map(|i| Spot {
                barcode: format!("bc{i}"),
                x: i as f64 * 100.0,
                y: 0.0,
                total_counts: i as u64 * 10,
                cluster: None,
            })
            .collect();

        // Expression increases monotonically left→right
        let expression: Vec<f64> = (0..10).map(|i| i as f64).collect();

        let (mi, _pv) = morans_i(&expression, &spots, 6, 100, 42);
        assert!(
            mi > 0.3,
            "Expected Moran's I > 0.3 for linearly increasing expression, got {mi}"
        );
    }

    /// Parse a minimal CSV string and verify spot coordinates.
    #[test]
    fn parse_visium_positions_basic() {
        let csv = b"barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n\
                    AAACAAGTATCTCCCA-1,1,0,16,1000,2000\n\
                    AAACACCAATAACTGC-1,0,1,17,1050,2050\n\
                    AAACAGAGCGACTCCT-1,1,2,18,1100,2100\n";

        let spots = parse_visium_positions(csv).expect("parse should succeed");
        // Only in_tissue == 1 spots are returned
        assert_eq!(spots.len(), 2);

        assert_eq!(spots[0].barcode, "AAACAAGTATCTCCCA-1");
        assert!((spots[0].x - 2000.0).abs() < 1e-6, "x should be pxl_col");
        assert!((spots[0].y - 1000.0).abs() < 1e-6, "y should be pxl_row");

        assert_eq!(spots[1].barcode, "AAACAGAGCGACTCCT-1");
        assert!((spots[1].x - 2100.0).abs() < 1e-6);
        assert!((spots[1].y - 1100.0).abs() < 1e-6);
    }

    /// Verify `find_spatially_variable_genes` returns exactly `n_top` genes (or fewer).
    #[test]
    fn spatially_variable_genes_returns_n() {
        // Build a small synthetic dataset: 20 spots in a 4×5 grid, 15 genes
        let n_spots = 20usize;
        let n_genes = 15usize;
        let n_top = 5usize;

        let spots: Vec<Spot> = (0..n_spots)
            .map(|i| Spot {
                barcode: format!("bc{i}"),
                x: (i % 4) as f64 * 100.0,
                y: (i / 4) as f64 * 100.0,
                total_counts: 100,
                cluster: None,
            })
            .collect();

        // Build count matrix with enough non-zero values to pass the 10% filter
        let mut matrix = Array2::<f32>::zeros((n_spots, n_genes));
        // Make every gene expressed in all spots (guaranteed to pass 10% filter)
        for s in 0..n_spots {
            for g in 0..n_genes {
                // Simple increasing pattern for half the genes to create variation
                matrix[[s, g]] = ((s + g) % 5 + 1) as f32;
            }
        }

        let gene_ids: Vec<String> = (0..n_genes).map(|i| format!("gene_{i}")).collect();

        let svgs = find_spatially_variable_genes(&matrix, &gene_ids, &spots, n_top, 6);

        assert!(
            svgs.len() <= n_top,
            "Should return at most n_top={n_top} genes, got {}",
            svgs.len()
        );
        assert!(
            !svgs.is_empty(),
            "Should return at least some SVGs from this dataset"
        );

        // Verify sorted by Moran's I descending
        for w in svgs.windows(2) {
            assert!(
                w[0].morans_i >= w[1].morans_i,
                "SVGs should be sorted by Moran's I descending"
            );
        }
    }
}
