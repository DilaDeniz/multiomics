//! Long-read (Nanopore / PacBio) methylation calling via MM/ML BAM auxiliary tags.
//!
//! # MM/ML tag format
//!
//! The SAM/BAM specification (v1.7+) defines two auxiliary tags for base
//! modification information:
//!
//! - **MM** (`MM:Z:…`) — Modified bases descriptor. Each semicolon-separated
//!   clause has the form `<base><strand><type>[?][,<skip>…]`, e.g.
//!   `C+m?,3,1,0;` means 5-methylcytosine (5mC) on the forward strand, with
//!   skip counts 3, 1, 0 between consecutive modified positions. The `?` flag
//!   signals that undetected positions should be treated as having unknown
//!   methylation status rather than being inferred as unmethylated.
//!
//! - **ML** (`ML:B:C,…`) — Probability array. Each `uint8` entry corresponds
//!   to one modification entry in MM, in the same order. Divide by 255 to
//!   obtain a probability in \[0.0, 1.0\].
//!
//! # Scope and limitations
//!
//! This module currently handles **5-methylcytosine (C+m) only**. Other
//! modification codes (e.g. `C+h` for 5-hydroxymethylcytosine, `A+a` for
//! 6-methyladenine) present in the MM tag are silently skipped.
//!
//! CIGAR operations supported for reference-coordinate mapping:
//! `M`, `=`, `X` (consuming both query and reference), `D`/`N` (reference
//! only), `I`/`S` (query only). Hard clips (`H`) and padding (`P`) are
//! handled correctly (neither query nor reference consumed).
//!
//! # Feature flag
//!
//! This module is gated behind the **`longread`** feature to avoid mandatory
//! dependencies on noodles-bam / noodles-sam. Enable with:
//! ```toml
//! epigenomics_core = { path = "…", features = ["longread"] }
//! ```

use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::path::Path;

use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
// The record::Cigar trait must be in scope to call `.iter()` on record_buf::Cigar.
use noodles_sam::alignment::record::Cigar as CigarTrait;
use noodles_sam::alignment::record_buf::data::field::value::Array as BufArray;
use noodles_sam::alignment::record_buf::data::field::Value as BufValue;

use crate::types::MethylationRecord;

// ── Public types ──────────────────────────────────────────────────────────────

/// A single methylation call derived from MM/ML tags in a long-read BAM record.
#[derive(Debug, Clone)]
pub struct LongReadMethCall {
    /// Reference chromosome / contig name.
    pub chrom: String,
    /// 0-based reference position of the modified cytosine.
    pub position: u64,
    /// Methylation probability in \[0.0, 1.0\] decoded from the ML tag
    /// (`raw_byte / 255.0`).
    pub methylation_prob: f64,
    /// `true` when `methylation_prob >= min_prob` supplied to
    /// [`parse_longread_methylation`].
    pub is_methylated: bool,
    /// `'+'` for forward-strand reads, `'-'` for reverse-strand reads.
    pub strand: char,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Parse MM/ML methylation tags from a BAM file.
///
/// Iterates all records in the file sequentially and, for each record:
/// 1. Skips unmapped records and records with MAPQ below `min_mapq`.
/// 2. Extracts the `MM:Z` auxiliary tag and locates the `C+m` (5mC) clause.
/// 3. Parses comma-separated skip counts to find query positions of modified
///    cytosines.
/// 4. Uses the CIGAR string to project each query position to a reference
///    coordinate.
/// 5. Reads the corresponding probability byte from the `ML:B:C` array and
///    divides by 255 to yield a probability in \[0.0, 1.0\].
///
/// Records that lack MM or ML tags are silently skipped.
///
/// # Arguments
/// * `path`     – Filesystem path to the BAM file.
/// * `min_mapq` – Minimum mapping quality; reads below this threshold are
///                discarded.
/// * `min_prob` – Probability threshold for calling a site methylated.
///                Typically 0.5.
///
/// # Returns
/// A `Vec<LongReadMethCall>` sorted by `(chrom, position)` ascending.
///
/// # Errors
/// Returns an [`anyhow::Error`] if the BAM file cannot be opened, the header
/// cannot be parsed, or a record is structurally invalid.
pub fn parse_longread_methylation(
    path: &Path,
    min_mapq: u8,
    min_prob: f64,
) -> Result<Vec<LongReadMethCall>> {
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("cannot open BAM file: {}", path.display()))?;

    let header = reader.read_header().context("failed to read BAM header")?;

    let mut calls: Vec<LongReadMethCall> = Vec::new();

    for result in reader.record_bufs(&header) {
        let record = result.context("failed to read BAM record")?;

        // Skip unmapped records (SAM flag 0x4).
        if record.flags().is_unmapped() {
            continue;
        }

        // Apply mapping-quality filter.
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(0);
        if mapq < min_mapq {
            continue;
        }

        // Resolve chromosome name from the reference sequence dictionary.
        let ref_id = match record.reference_sequence_id() {
            Some(id) => id,
            None => continue,
        };
        let chrom = match header.reference_sequences().get_index(ref_id) {
            Some((name, _)) => name.to_string(),
            None => continue,
        };

        // Reference start position — noodles `Position` is 1-based; convert to 0-based u64.
        let ref_start = match record.alignment_start() {
            Some(pos) => usize::from(pos) as u64 - 1,
            None => continue,
        };

        let strand = if record.flags().is_reverse_complemented() {
            '-'
        } else {
            '+'
        };

        // Retrieve MM and ML auxiliary tags.
        let data = record.data();

        // Data::get() accepts any K where Tag: Borrow<K> and K: Eq.
        // [u8; 2] satisfies this via Tag: Borrow<[u8; 2]>.
        let mm_str = match data.get(b"MM") {
            Some(BufValue::String(s)) => String::from_utf8_lossy(s).into_owned(),
            _ => continue, // No MM tag — skip this record.
        };

        let ml_raw: Vec<u8> = match data.get(b"ML") {
            Some(BufValue::Array(BufArray::UInt8(bytes))) => bytes.clone(),
            _ => continue, // No ML tag — skip this record.
        };

        // Parse query-level 5mC positions and map to reference coordinates.
        let cigar = record.cigar();
        let record_calls =
            extract_5mc_calls(&chrom, ref_start, strand, &mm_str, &ml_raw, cigar, min_prob)
                .with_context(|| {
                    format!(
                        "failed to extract 5mC calls for record at {}:{}",
                        chrom, ref_start
                    )
                })?;

        calls.extend(record_calls);
    }

    // Sort by (chrom, position) for deterministic, merge-friendly output.
    calls.sort_unstable_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| a.position.cmp(&b.position))
    });

    Ok(calls)
}

/// Convert a slice of [`LongReadMethCall`] records to [`MethylationRecord`]
/// format for use with [`crate::accum::EpigenomicsAccum`].
///
/// Methylation probability is scaled to the \[0.0, 100.0\] percentage range
/// used by the ENCODE BED format: `prob × 100.0`. All sites (methylated or
/// not) are included so that the accumulator can compute per-site statistics.
pub fn longread_to_methylation_records(calls: &[LongReadMethCall]) -> Vec<MethylationRecord> {
    calls
        .iter()
        .map(|c| MethylationRecord {
            chrom: c.chrom.clone(),
            start: c.position,
            end: c.position + 1,
            methylation: c.methylation_prob * 100.0,
        })
        .collect()
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Parse the MM tag string, locate the first `C+m` clause, extract
/// query-relative modification positions, then map each to a reference
/// coordinate via the CIGAR string.
///
/// Returns `Vec<LongReadMethCall>` — one per modified cytosine found.
fn extract_5mc_calls(
    chrom: &str,
    ref_start: u64,
    strand: char,
    mm_str: &str,
    ml_raw: &[u8],
    cigar: &noodles_sam::alignment::record_buf::Cigar,
    min_prob: f64,
) -> Result<Vec<LongReadMethCall>> {
    // The MM tag may contain multiple semicolon-separated modification clauses.
    // Find the first `C+m` clause (5mC, forward strand).
    let cm_clause = mm_str.split(';').find(|clause| {
        let c = clause.trim_start();
        c.starts_with("C+m") || c.starts_with("c+m")
    });

    let Some(clause) = cm_clause else {
        return Ok(Vec::new()); // No 5mC modification in this record.
    };

    // Strip the type prefix (`C+m` or `C+m?`) and parse the skip counts that
    // follow the comma separator.
    //
    // Format: `C+m?,<skip0>,<skip1>,...`
    // The skip count at position i is the number of same-base (C) positions in
    // the read that are skipped *before* the i-th modified base.
    let after_prefix = strip_cm_prefix(clause)?;

    if after_prefix.is_empty() {
        return Ok(Vec::new());
    }

    let skip_counts: Vec<u64> = after_prefix
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.trim()
                .parse::<u64>()
                .with_context(|| format!("invalid skip count in MM tag: {s:?}"))
        })
        .collect::<Result<Vec<_>>>()?;

    if skip_counts.is_empty() {
        return Ok(Vec::new());
    }

    // Determine the index into ml_raw for the first C+m entry, accounting for
    // any modification clauses that appear before C+m in the MM string.
    let ml_offset = ml_offset_for_cm(mm_str)?;

    // Build a query-position → reference-position lookup table from the CIGAR.
    let q_to_ref = query_to_ref_map(ref_start, cigar)?;
    let query_len = cigar_query_len(cigar)?;

    // Walk through skip counts to recover each modified cytosine's query
    // position.  Skip counts are cumulative between consecutive modified Cs:
    // after landing on modified C[i], advance by skip[i+1]+1 to reach C[i+1].
    let mut calls: Vec<LongReadMethCall> = Vec::with_capacity(skip_counts.len());
    let mut query_pos: u64 = 0;

    for (i, &skip) in skip_counts.iter().enumerate() {
        query_pos += skip + 1;

        if query_pos as usize > query_len {
            break; // Ran off the end of the read — malformed record, stop.
        }

        // Map 0-based query position to reference coordinate.
        let Some(&ref_pos) = q_to_ref.get(&(query_pos - 1)) else {
            // Position falls inside an insertion or soft-clip — no ref coord.
            continue;
        };

        // Fetch the corresponding ML probability byte.
        let ml_idx = ml_offset + i;
        let prob_byte = ml_raw.get(ml_idx).copied().unwrap_or(0);
        let methylation_prob = prob_byte as f64 / 255.0;

        calls.push(LongReadMethCall {
            chrom: chrom.to_string(),
            position: ref_pos,
            methylation_prob,
            is_methylated: methylation_prob >= min_prob,
            strand,
        });
    }

    Ok(calls)
}

/// Strip the `C+m` (or `C+m?`) prefix from an MM tag clause and return the
/// portion after the first comma (the skip count list).
///
/// Returns an empty string when the clause has no skip counts (zero modified
/// bases). Returns an error when the prefix is unrecognised.
fn strip_cm_prefix(clause: &str) -> Result<&str> {
    // Possible prefixes, longest-first to avoid prefix ambiguity.
    for prefix in &["C+m?,", "c+m?,", "C+m,", "c+m,"] {
        if let Some(rest) = clause.strip_prefix(prefix) {
            return Ok(rest);
        }
    }
    // Clause with no skip counts at all (e.g. exactly "C+m" or "C+m?").
    if matches!(clause, "C+m" | "C+m?" | "c+m" | "c+m?") {
        return Ok("");
    }
    bail!("unexpected C+m clause format: {:?}", clause);
}

/// Count how many ML entries are consumed by modification clauses that appear
/// *before* the first `C+m` clause in the MM tag string.
///
/// Each non-C+m clause contributes as many ML entries as it has skip-count
/// values (i.e., the count of commas after the type prefix).
fn ml_offset_for_cm(mm_str: &str) -> Result<usize> {
    let mut offset = 0_usize;
    for clause in mm_str.split(';') {
        let c = clause.trim_start();
        if c.starts_with("C+m") || c.starts_with("c+m") {
            break;
        }
        // Count the skip-count values: everything after the first comma.
        if let Some(rest) = c.splitn(2, ',').nth(1) {
            offset += rest.split(',').filter(|s| !s.is_empty()).count();
        }
    }
    Ok(offset)
}

/// Build a `HashMap<query_pos, ref_pos>` from a CIGAR string and a reference
/// start position.
///
/// Both positions are 0-based. Only query positions that align to the reference
/// (CIGAR ops `M`, `=`, `X`) are inserted. Insertions and soft-clips advance
/// the query counter without a reference entry. Deletions and skips advance
/// only the reference counter.
fn query_to_ref_map(
    ref_start: u64,
    cigar: &noodles_sam::alignment::record_buf::Cigar,
) -> Result<HashMap<u64, u64>> {
    let mut map = HashMap::new();
    let mut q: u64 = 0;
    let mut r: u64 = ref_start;

    for op_result in cigar.iter() {
        let op = op_result.context("invalid CIGAR operation")?;
        let len = op.len() as u64;
        match op.kind() {
            CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                for k in 0..len {
                    map.insert(q + k, r + k);
                }
                q += len;
                r += len;
            }
            CigarKind::Insertion | CigarKind::SoftClip => {
                q += len;
            }
            CigarKind::Deletion | CigarKind::Skip => {
                r += len;
            }
            CigarKind::HardClip | CigarKind::Pad => {
                // Hard clips and padding consume neither query nor reference.
            }
        }
    }

    Ok(map)
}

/// Return the total query-consuming length from a CIGAR (M, I, S, =, X ops).
fn cigar_query_len(cigar: &noodles_sam::alignment::record_buf::Cigar) -> Result<usize> {
    let mut total = 0_usize;
    for op_result in cigar.iter() {
        let op = op_result.context("invalid CIGAR operation")?;
        match op.kind() {
            CigarKind::Match
            | CigarKind::Insertion
            | CigarKind::SoftClip
            | CigarKind::SequenceMatch
            | CigarKind::SequenceMismatch => {
                total += op.len();
            }
            _ => {}
        }
    }
    Ok(total)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_longread_to_methylation_records_conversion() {
        let calls = vec![
            LongReadMethCall {
                chrom: "chr1".to_string(),
                position: 1000,
                methylation_prob: 0.9,
                is_methylated: true,
                strand: '+',
            },
            LongReadMethCall {
                chrom: "chr1".to_string(),
                position: 2000,
                methylation_prob: 0.1,
                is_methylated: false,
                strand: '-',
            },
        ];

        let records = longread_to_methylation_records(&calls);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].start, 1000);
        assert_eq!(records[0].end, 1001);
        assert!((records[0].methylation - 90.0).abs() < 1e-9);
        assert!((records[1].methylation - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_ml_offset_no_prior_clauses() {
        // C+m is the only clause — offset should be 0.
        let offset = ml_offset_for_cm("C+m?,3,1,0").unwrap();
        assert_eq!(offset, 0);
    }

    #[test]
    fn test_ml_offset_with_prior_clause() {
        // A prior clause with 3 skip-count values appears before C+m.
        let offset = ml_offset_for_cm("A+a,1,2,3;C+m?,0,1").unwrap();
        assert_eq!(offset, 3);
    }

    #[test]
    fn test_strip_cm_prefix_variants() {
        assert_eq!(strip_cm_prefix("C+m?,3,1,0").unwrap(), "3,1,0");
        assert_eq!(strip_cm_prefix("C+m,3,1").unwrap(), "3,1");
        assert_eq!(strip_cm_prefix("c+m?,0").unwrap(), "0");
        assert_eq!(strip_cm_prefix("C+m").unwrap(), "");
        assert_eq!(strip_cm_prefix("C+m?").unwrap(), "");
    }

    #[test]
    fn test_strip_cm_prefix_invalid() {
        assert!(strip_cm_prefix("G+m?,1,2").is_err());
    }
}
