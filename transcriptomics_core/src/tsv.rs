use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::types::GeneRecord;

/// Parse a genes-by-samples expression matrix TSV file.
///
/// Expected format:
/// ```text
/// gene_id\tsample1\tsample2\t...
/// GAPDH\t1234.5\t987.6\t...
/// ```
/// The first column is gene identifier; subsequent columns are TPM values.
///
/// Returns both the parsed records and the ordered sample names extracted
/// from the header row.
///
/// # Errors
/// Returns an error when the file cannot be opened, the header is missing,
/// or a data row has fewer columns than the header.
pub fn parse_tsv(path: &Path) -> Result<(Vec<GeneRecord>, Vec<String>)> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open expression matrix TSV '{}'", path.display()))?;

    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Read and parse header
    let header_line = lines
        .next()
        .with_context(|| {
            format!(
                "Expression matrix TSV '{}' is empty",
                path.display()
            )
        })?
        .with_context(|| format!("Read error in '{}'", path.display()))?;

    let header_cols: Vec<&str> = header_line.split('\t').collect();
    if header_cols.len() < 2 {
        bail!(
            "Expression matrix TSV '{}' must have at least one sample column (found {})",
            path.display(),
            header_cols.len()
        );
    }

    // First column is gene_id label; the rest are sample names
    let sample_names: Vec<String> = header_cols[1..].iter().map(|s| s.trim().to_string()).collect();
    let n_samples = sample_names.len();

    let mut records = Vec::new();
    let mut line_no = 1usize;

    for line in lines {
        line_no += 1;
        let line =
            line.with_context(|| format!("Read error at line {} of '{}'", line_no, path.display()))?;

        if line.trim().is_empty() {
            continue;
        }

        let cols: Vec<&str> = line.splitn(n_samples + 2, '\t').collect();
        if cols.len() < 2 {
            log::warn!(
                "Skipping line {} in '{}': too few columns",
                line_no,
                path.display()
            );
            continue;
        }

        let gene_id = cols[0].trim().to_string();
        let mut samples = Vec::with_capacity(n_samples);
        for i in 1..=n_samples {
            let val_str = cols.get(i).unwrap_or(&"0").trim();
            let val: f64 = val_str.parse().unwrap_or_else(|_| {
                log::warn!(
                    "Cannot parse TPM '{}' at line {} col {} of '{}', using 0",
                    val_str, line_no, i, path.display()
                );
                0.0
            });
            samples.push(val);
        }

        records.push(GeneRecord { gene_id, samples });
    }

    log::info!(
        "Parsed {} gene records × {} samples from '{}'",
        records.len(),
        n_samples,
        path.display()
    );

    Ok((records, sample_names))
}
