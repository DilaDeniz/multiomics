use std::path::Path;

use anyhow::{bail, Context, Result};
use memmap2::Mmap;

use crate::types::GeneRecord;
use biomics_core::parse::{parse_f64, trim_bytes, ByteLines, TabFields};

/// Parse a genes-by-samples expression matrix TSV file using memory-mapped
/// zero-allocation byte-level parsing.
///
/// Expected format:
/// ```text
/// gene_id\tsample1\tsample2\t...
/// GAPDH\t1234.5\t987.6\t...
/// ```
///
/// Returns both the parsed records and the ordered sample names from the header.
pub fn parse_tsv(path: &Path) -> Result<(Vec<GeneRecord>, Vec<String>)> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open expression matrix TSV '{}'", path.display()))?;

    // SAFETY: the file is not modified while this process holds the Mmap.
    let mmap = unsafe { Mmap::map(&file) }
        .with_context(|| format!("Cannot mmap TSV file '{}'", path.display()))?;

    let _ = mmap.advise(memmap2::Advice::Sequential);

    let data = mmap.as_ref();
    let mut lines = ByteLines::new(data);

    // ── Header ────────────────────────────────────────────────────────────────
    let header = lines
        .next()
        .with_context(|| format!("Expression matrix TSV '{}' is empty", path.display()))?;

    let header_fields: Vec<&[u8]> = TabFields::new(header).collect();
    if header_fields.len() < 2 {
        bail!(
            "Expression matrix TSV '{}' must have at least one sample column (found {})",
            path.display(),
            header_fields.len()
        );
    }

    let sample_names: Vec<String> = header_fields[1..]
        .iter()
        .map(|f| String::from_utf8_lossy(trim_bytes(f)).into_owned())
        .collect();
    let n_samples = sample_names.len();

    // ── Data rows ─────────────────────────────────────────────────────────────
    // Estimate: ~50 bytes per row average (short gene IDs + float values)
    let mut records = Vec::with_capacity(data.len() / 50);

    for line in lines {
        let line = trim_bytes(line);
        if line.is_empty() {
            continue;
        }

        let mut fields = TabFields::new(line);
        let gene_id_bytes = match fields.next() {
            Some(b) if !b.is_empty() => b,
            _ => continue,
        };
        let gene_id = String::from_utf8_lossy(trim_bytes(gene_id_bytes)).into_owned();

        let mut samples = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            let val = fields
                .next()
                .and_then(|b| parse_f64(trim_bytes(b)))
                .unwrap_or(0.0);
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
