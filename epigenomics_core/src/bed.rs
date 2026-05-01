use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use memmap2::Mmap;

use crate::types::MethylationRecord;

/// Parse a BED methylation file into a `Vec<MethylationRecord>`.
///
/// Supports both 4-column (chrom, start, end, methylation_pct) and
/// 6-column ENCODE bisulfite format (chrom, start, end, name, score, strand)
/// where `score` is in [0, 1000] (divide by 10 for percentage).
///
/// Lines starting with `#` or `track` are treated as headers and skipped.
/// Malformed lines are logged and skipped.
///
/// # Errors
/// Returns an error when the file cannot be opened or memory-mapped.
pub fn parse_bed(path: &Path) -> Result<Vec<MethylationRecord>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open BED file '{}'", path.display()))?;

    // SAFETY: the file is not modified while this process holds the Mmap.
    let mmap = unsafe { Mmap::map(&file) }
        .with_context(|| format!("Cannot mmap BED file '{}'", path.display()))?;

    let reader = BufReader::new(mmap.as_ref());
    let mut records = Vec::new();
    let mut line_no = 0usize;

    for line in reader.lines() {
        line_no += 1;
        let line = line
            .with_context(|| format!("Read error at line {} of '{}'", line_no, path.display()))?;

        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }

        match parse_bed_line(trimmed) {
            Some(record) => records.push(record),
            None => {
                log::warn!(
                    "Skipping unparseable BED line {} in '{}': {}",
                    line_no,
                    path.display(),
                    &trimmed[..trimmed.len().min(80)]
                );
            }
        }
    }

    log::info!(
        "Parsed {} methylation records from '{}'",
        records.len(),
        path.display()
    );

    Ok(records)
}

/// Parse a single BED data line.
fn parse_bed_line(line: &str) -> Option<MethylationRecord> {
    let mut cols = line.split('\t');
    let chrom = cols.next()?.to_string();
    let start: u64 = cols.next()?.parse().ok()?;
    let end: u64 = cols.next()?.parse().ok()?;

    // Try to determine methylation value from column 4 or 5
    let col4 = cols.next();
    let col5 = cols.next();

    let methylation = if let Some(score_str) = col5 {
        // 6-column ENCODE format: col5 is score in [0, 1000]
        let score: f64 = score_str.parse().ok()?;
        score / 10.0 // convert to percentage
    } else if let Some(val_str) = col4 {
        // 4-column format: col4 is direct methylation percentage
        let val: f64 = val_str.parse().ok()?;
        // Heuristic: if value > 1.0, treat as percentage; if ≤ 1.0, scale to %
        if val > 1.0 {
            val.min(100.0)
        } else {
            val * 100.0
        }
    } else {
        return None;
    };

    Some(MethylationRecord {
        chrom,
        start,
        end,
        methylation: methylation.clamp(0.0, 100.0),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bed_line_6col() {
        let line = "chr1\t100\t101\tCpG\t850\t+";
        let rec = parse_bed_line(line).unwrap();
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.start, 100);
        assert!((rec.methylation - 85.0).abs() < 1e-6);
    }

    #[test]
    fn test_parse_bed_line_4col() {
        let line = "chr1\t200\t201\t72.5";
        let rec = parse_bed_line(line).unwrap();
        assert!((rec.methylation - 72.5).abs() < 1e-6);
    }
}
