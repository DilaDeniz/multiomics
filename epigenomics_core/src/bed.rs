use std::path::Path;

use anyhow::{Context, Result};
use memmap2::Mmap;

use crate::types::MethylationRecord;
use biomics_core::parse::{parse_f64, parse_u64, trim_bytes, ByteLines, TabFields};

/// Parse a BED methylation file into a `Vec<MethylationRecord>` using
/// memory-mapped zero-allocation byte-level parsing.
///
/// Supports:
/// - 4-column: `chrom  start  end  methylation_pct`
/// - 6-column ENCODE bisulfite: `chrom  start  end  name  score  strand`
///   where `score` ∈ [0, 1000] → divide by 10 for percentage.
///
/// `madvise(SEQUENTIAL)` is applied for kernel read-ahead.
pub fn parse_bed(path: &Path) -> Result<Vec<MethylationRecord>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open BED file '{}'", path.display()))?;

    // SAFETY: the file is not modified while this process holds the Mmap.
    let mmap = unsafe { Mmap::map(&file) }
        .with_context(|| format!("Cannot mmap BED file '{}'", path.display()))?;

    let _ = mmap.advise(memmap2::Advice::Sequential);

    let data = mmap.as_ref();
    // Estimate: a typical 6-col BED line is ~40 bytes
    let mut records = Vec::with_capacity(data.len() / 40);

    for line in ByteLines::new(data) {
        let line = trim_bytes(line);
        if line.is_empty()
            || line[0] == b'#'
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }
        match parse_bed_line(line) {
            Some(record) => records.push(record),
            None => log::warn!(
                "Skipping unparseable BED line: {}",
                String::from_utf8_lossy(&line[..line.len().min(80)])
            ),
        }
    }

    log::info!(
        "Parsed {} methylation records from '{}'",
        records.len(),
        path.display()
    );
    Ok(records)
}

/// Parse a single BED data line from a byte slice.
#[inline]
fn parse_bed_line(line: &[u8]) -> Option<MethylationRecord> {
    let mut cols = TabFields::new(line);

    let chrom = std::str::from_utf8(cols.next()?).ok()?.to_string();
    let start = parse_u64(cols.next()?)?;
    let end = parse_u64(cols.next()?)?;

    let col4 = cols.next();
    let col5 = cols.next();

    // Determine gene name: col4 is a gene name if it exists and cannot be parsed as f64.
    let gene: Option<String> = col4.and_then(|b| {
        let s = std::str::from_utf8(trim_bytes(b)).ok()?;
        // If col5 exists, col4 is the name field in 6-col ENCODE format.
        // If col5 doesn't exist, col4 is the methylation value — treat as numeric.
        if col5.is_some() {
            // col4 is name; only keep if non-numeric
            if s.parse::<f64>().is_err() && !s.is_empty() {
                Some(s.to_string())
            } else {
                None
            }
        } else {
            // col4 is the methylation value (4-col format) — no gene name
            None
        }
    });

    let methylation = if let Some(score_b) = col5 {
        // 6-column ENCODE: score in [0, 1000]
        let score: f64 = parse_f64(trim_bytes(score_b))?;
        score / 10.0
    } else if let Some(val_b) = col4 {
        let val: f64 = parse_f64(trim_bytes(val_b))?;
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
        gene,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bed_line_6col() {
        let line = b"chr1\t100\t101\tCpG\t850\t+";
        let rec = parse_bed_line(line).unwrap();
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.start, 100);
        assert!((rec.methylation - 85.0).abs() < 1e-6);
        assert_eq!(rec.gene, Some("CpG".to_string()));
    }

    #[test]
    fn test_parse_bed_line_6col_with_gene() {
        let line = b"chr1\t100\t101\tBRCA1\t750\t+";
        let rec = parse_bed_line(line).unwrap();
        assert_eq!(rec.gene, Some("BRCA1".to_string()));
        assert!((rec.methylation - 75.0).abs() < 1e-6);
    }

    #[test]
    fn test_parse_bed_line_4col() {
        let line = b"chr1\t200\t201\t72.5";
        let rec = parse_bed_line(line).unwrap();
        assert!((rec.methylation - 72.5).abs() < 1e-6);
        assert_eq!(rec.gene, None);
    }
}
