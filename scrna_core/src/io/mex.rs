//! 10x Genomics MEX format reader.
//!
//! Parses `matrix.mtx`, `barcodes.tsv`, and `features.tsv` (plain or `.gz`).

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Compressed-sparse-row count matrix from a 10x MEX directory.
pub struct CsrMatrix {
    /// Number of features (genes) — rows in the original MTX.
    pub n_rows: usize,
    /// Number of cells — columns in the original MTX.
    pub n_cols: usize,
    /// Row pointers (length `n_rows + 1`).
    pub indptr: Vec<u32>,
    /// Column indices (0-indexed), sorted within each row.
    pub indices: Vec<u32>,
    /// Non-zero count values.
    pub data: Vec<u32>,
    /// Cell barcodes (one per column).
    pub barcodes: Vec<String>,
    /// Feature / gene names (one per row).
    pub features: Vec<String>,
}

impl CsrMatrix {
    /// Iterate over `(column_index, value)` pairs for row `i`.
    pub fn row(&self, i: usize) -> impl Iterator<Item = (u32, u32)> + '_ {
        let start = self.indptr[i] as usize;
        let end = self.indptr[i + 1] as usize;
        self.indices[start..end]
            .iter()
            .zip(self.data[start..end].iter())
            .map(|(&col, &val)| (col, val))
    }

    /// Total UMI counts for cell `j` (sums over all features). O(nnz).
    pub fn col_sum(&self, j: usize) -> u32 {
        let mut sum = 0u32;
        for row in 0..self.n_rows {
            for (col, val) in self.row(row) {
                if col as usize == j {
                    sum += val;
                }
            }
        }
        sum
    }

    /// Number of non-zero features per cell (length `n_cols`). O(nnz).
    pub fn nnz_per_col(&self) -> Vec<u32> {
        let mut counts = vec![0u32; self.n_cols];
        for &col in &self.indices {
            counts[col as usize] += 1;
        }
        counts
    }
}

/// Open a file, preferring `<name>.gz` when both exist, falling back to plain.
fn open_reader(dir: &Path, stem: &str) -> Result<Box<dyn BufRead>> {
    let gz_path = dir.join(format!("{stem}.gz"));
    let plain_path = dir.join(stem);
    if gz_path.exists() {
        let f = File::open(&gz_path).with_context(|| format!("opening {}", gz_path.display()))?;
        Ok(Box::new(BufReader::new(GzDecoder::new(f))))
    } else {
        let f =
            File::open(&plain_path).with_context(|| format!("opening {}", plain_path.display()))?;
        Ok(Box::new(BufReader::new(f)))
    }
}

fn read_lines_trimmed(reader: Box<dyn BufRead>) -> Result<Vec<String>> {
    let mut out = Vec::new();
    for line in reader.lines() {
        let l = line.context("reading line")?;
        let t = l.trim().to_owned();
        if !t.is_empty() {
            out.push(t);
        }
    }
    Ok(out)
}

/// Parse a 10x MEX directory, supporting `.gz` compressed files.
///
/// Expects `matrix.mtx[.gz]`, `barcodes.tsv[.gz]`, and `features.tsv[.gz]`.
pub fn parse_10x_mex(dir: &Path) -> Result<CsrMatrix> {
    let barcodes = {
        let rdr = open_reader(dir, "barcodes.tsv")?;
        read_lines_trimmed(rdr).context("reading barcodes.tsv")?
    };

    let features: Vec<String> = {
        let rdr = open_reader(dir, "features.tsv")?;
        let lines = read_lines_trimmed(rdr).context("reading features.tsv")?;
        // features.tsv has tab-separated columns; take the second column (gene symbol)
        // if present, otherwise the first.
        lines
            .into_iter()
            .map(|line| {
                let mut parts = line.splitn(3, '\t');
                let first = parts.next().unwrap_or("").to_owned();
                let second = parts.next().map(|s| s.to_owned());
                second.filter(|s| !s.is_empty()).unwrap_or(first)
            })
            .collect()
    };

    let (n_rows, n_cols, coo) = parse_mtx(dir)?;

    if n_rows != features.len() {
        anyhow::bail!(
            "matrix.mtx has {} features but features.tsv has {}",
            n_rows,
            features.len()
        );
    }
    if n_cols != barcodes.len() {
        anyhow::bail!(
            "matrix.mtx has {} cells but barcodes.tsv has {}",
            n_cols,
            barcodes.len()
        );
    }

    let (indptr, indices, data) = coo_to_csr(n_rows, coo);

    Ok(CsrMatrix {
        n_rows,
        n_cols,
        indptr,
        indices,
        data,
        barcodes,
        features,
    })
}

/// COO entry: (row, col, value), all 0-indexed.
type CooEntry = (u32, u32, u32);

/// Parse matrix.mtx, returning `(n_features, n_barcodes, COO entries)`.
fn parse_mtx(dir: &Path) -> Result<(usize, usize, Vec<CooEntry>)> {
    let rdr = open_reader(dir, "matrix.mtx")?;
    let mut lines = rdr.lines();

    // Skip comment/header lines beginning with '%'
    let mut dim_line = String::new();
    for line in lines.by_ref() {
        let l = line.context("reading matrix.mtx")?;
        if l.starts_with('%') {
            continue;
        }
        dim_line = l;
        break;
    }

    let mut parts = dim_line.split_whitespace();
    let n_rows: usize = parts
        .next()
        .context("missing n_rows in matrix.mtx header")?
        .parse()
        .context("parsing n_rows")?;
    let n_cols: usize = parts
        .next()
        .context("missing n_cols in matrix.mtx header")?
        .parse()
        .context("parsing n_cols")?;
    let nnz: usize = parts
        .next()
        .context("missing nnz in matrix.mtx header")?
        .parse()
        .context("parsing nnz")?;

    let mut coo: Vec<(u32, u32, u32)> = Vec::with_capacity(nnz);
    for line in lines {
        let l = line.context("reading matrix.mtx entry")?;
        let l = l.trim();
        if l.is_empty() || l.starts_with('%') {
            continue;
        }
        let mut p = l.split_whitespace();
        let row: u32 = p
            .next()
            .context("missing row in entry")?
            .parse::<u32>()
            .context("parsing row")?
            .checked_sub(1)
            .context("row index 0 is invalid (MTX is 1-indexed)")?;
        let col: u32 = p
            .next()
            .context("missing col in entry")?
            .parse::<u32>()
            .context("parsing col")?
            .checked_sub(1)
            .context("col index 0 is invalid (MTX is 1-indexed)")?;
        let val: u32 = p
            .next()
            .context("missing value in entry")?
            .parse()
            .context("parsing value")?;
        coo.push((row, col, val));
    }

    Ok((n_rows, n_cols, coo))
}

/// Convert COO triples (sorted stable by row) into CSR arrays.
fn coo_to_csr(n_rows: usize, mut coo: Vec<(u32, u32, u32)>) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    coo.sort_by_key(|&(r, c, _)| (r, c));

    let nnz = coo.len();
    let mut indptr = vec![0u32; n_rows + 1];
    let mut indices = Vec::with_capacity(nnz);
    let mut data = Vec::with_capacity(nnz);

    for &(row, col, val) in &coo {
        indptr[row as usize + 1] += 1;
        indices.push(col);
        data.push(val);
    }
    // Prefix-sum to build indptr
    for i in 1..=n_rows {
        indptr[i] += indptr[i - 1];
    }

    (indptr, indices, data)
}
