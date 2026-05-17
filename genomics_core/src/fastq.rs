use std::path::Path;

use anyhow::{Context, Result};

/// Summary statistics from an optional FASTQ input.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FastqSummary {
    pub total_reads: u64,
    pub total_bases: u64,
    pub mean_read_length: f64,
    pub gc_content_pct: f64,
    /// Percentage of bases with Phred quality ≥ 30.
    pub q30_pct: f64,
}

/// Parse a FASTQ file and return basic QC statistics.
///
/// Uses `needletail` for zero-copy FASTQ/FASTA parsing.
///
/// # Errors
/// Returns an error when the file cannot be opened or parsed.
pub fn parse_fastq(path: &Path) -> Result<FastqSummary> {
    use needletail::parse_fastx_file;

    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("Cannot open FASTQ file '{}'", path.display()))?;

    let mut total_reads = 0u64;
    let mut total_bases = 0u64;
    let mut gc_bases = 0u64;
    let mut q30_bases = 0u64;

    while let Some(record) = reader.next() {
        let rec = record.with_context(|| format!("Parse error in FASTQ '{}'", path.display()))?;

        total_reads += 1;
        let seq = rec.seq();
        total_bases += seq.len() as u64;

        for &b in seq.iter() {
            if matches!(b.to_ascii_uppercase(), b'G' | b'C') {
                gc_bases += 1;
            }
        }

        if let Some(quals) = rec.qual() {
            for &q in quals.iter() {
                // Illumina offset 33: q - 33 gives Phred score
                if q.saturating_sub(33) >= 30 {
                    q30_bases += 1;
                }
            }
        }
    }

    let mean_read_length = if total_reads == 0 {
        0.0
    } else {
        total_bases as f64 / total_reads as f64
    };

    let gc_content_pct = if total_bases == 0 {
        0.0
    } else {
        gc_bases as f64 / total_bases as f64 * 100.0
    };

    let q30_pct = if total_bases == 0 {
        0.0
    } else {
        q30_bases as f64 / total_bases as f64 * 100.0
    };

    log::info!(
        "FASTQ '{}': {} reads, {} bases, GC={:.1}%, Q30={:.1}%",
        path.display(),
        total_reads,
        total_bases,
        gc_content_pct,
        q30_pct
    );

    Ok(FastqSummary {
        total_reads,
        total_bases,
        mean_read_length,
        gc_content_pct,
        q30_pct,
    })
}
