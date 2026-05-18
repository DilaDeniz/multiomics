//! BAM pileup: converts sorted BAM to per-position read stacks.
//!
//! Reads a sorted BAM file using noodles-bam and produces [`PileupColumn`]
//! records. Each column contains all reads overlapping a genomic position
//! after applying base quality and mapping quality filters.
//!
//! # Reference
//! Li H, et al. (2009) The Sequence Alignment/Map format and SAMtools.
//! Bioinformatics 25(16):2078–2079.

use anyhow::{Context, Result};
use std::path::Path;

use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles_sam::alignment::record::Cigar as CigarTrait;

// ── Public types ──────────────────────────────────────────────────────────────

/// A single base observation at one genomic position from one read.
#[derive(Debug, Clone)]
pub struct PileupBase {
    /// Observed base: one of b'A', b'C', b'G', b'T', b'N'.
    pub base: u8,
    /// Phred-scaled base quality score.
    pub base_qual: u8,
    /// Mapping quality of the parent read.
    pub mapq: u8,
    /// `true` when the read is on the reverse strand.
    pub is_rev: bool,
}

/// All bases from reads that overlap a single reference position.
#[derive(Debug, Clone)]
pub struct PileupColumn {
    /// Reference sequence (chromosome) name.
    pub chrom: String,
    /// 0-based reference position.
    pub pos: u64,
    /// Reference base at this position (`b'N'` if unknown).
    pub ref_base: u8,
    /// Per-read base observations at this position.
    pub bases: Vec<PileupBase>,
}

impl PileupColumn {
    /// Number of reads covering this position.
    pub fn depth(&self) -> usize {
        self.bases.len()
    }
}

// ── Internal tuple used during collection ─────────────────────────────────────

struct BaseTuple {
    chrom_idx: u32,
    pos: u64,
    base: u8,
    base_qual: u8,
    mapq: u8,
    is_rev: bool,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Build a pileup from a sorted BAM file.
///
/// Each returned [`PileupColumn`] contains all filtered bases at one reference
/// position. Positions with zero depth after filtering are omitted.
///
/// # Arguments
/// * `bam_path`      – Path to a coordinate-sorted BAM file.
/// * `min_base_qual` – Minimum Phred base quality; bases below this are dropped.
/// * `min_mapq`      – Minimum mapping quality; reads below this are skipped.
/// * `max_depth`     – Maximum bases per column (SAMtools default: 8000).
///
/// # Errors
/// Returns an error if the BAM file cannot be opened or a record is malformed.
pub fn build_pileup(
    bam_path: &Path,
    min_base_qual: u8,
    min_mapq: u8,
    max_depth: usize,
) -> Result<Vec<PileupColumn>> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .with_context(|| format!("cannot open BAM file: {}", bam_path.display()))?;

    let header = reader.read_header().context("failed to read BAM header")?;

    // Collect all (chrom_idx, pos, base, qual, mapq, strand) tuples.
    let mut tuples: Vec<BaseTuple> = Vec::new();

    for result in reader.record_bufs(&header) {
        let record = result.context("failed to read BAM record")?;

        // Skip unmapped, secondary, supplementary, and duplicate reads.
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
        {
            continue;
        }

        // Apply mapping quality filter.
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(0);
        if mapq < min_mapq {
            continue;
        }

        // Resolve chromosome index.
        let ref_id = match record.reference_sequence_id() {
            Some(id) => id as u32,
            None => continue,
        };

        // Reference start (noodles Position is 1-based; convert to 0-based).
        let ref_start = match record.alignment_start() {
            Some(pos) => usize::from(pos) as u64 - 1,
            None => continue,
        };

        let is_rev = flags.is_reverse_complemented();

        let cigar = record.cigar();
        let seq = record.sequence();
        let qual_bytes: &[u8] = record.quality_scores().as_ref();

        // Walk CIGAR to emit (ref_pos, query_idx) pairs.
        let mut q_idx: usize = 0;
        let mut r_pos: u64 = ref_start;

        for op_result in cigar.iter() {
            let op = op_result.context("invalid CIGAR operation")?;
            let len = op.len();

            match op.kind() {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    for k in 0..len {
                        let qi = q_idx + k;
                        // Sequence returns raw IUPAC bytes; normalize to ACGTN.
                        let base_byte = normalize_base(seq.get(qi).unwrap_or(b'N'));
                        let bq = qual_bytes.get(qi).copied().unwrap_or(0);

                        if bq >= min_base_qual && base_byte != b'N' {
                            tuples.push(BaseTuple {
                                chrom_idx: ref_id,
                                pos: r_pos + k as u64,
                                base: base_byte,
                                base_qual: bq,
                                mapq,
                                is_rev,
                            });
                        }
                    }
                    q_idx += len;
                    r_pos += len as u64;
                }
                CigarKind::Insertion | CigarKind::SoftClip => {
                    q_idx += len;
                }
                CigarKind::Deletion | CigarKind::Skip => {
                    r_pos += len as u64;
                }
                CigarKind::HardClip | CigarKind::Pad => {}
            }
        }
    }

    // Sort by (chrom_idx, pos) to group into columns.
    tuples.sort_unstable_by_key(|t| (t.chrom_idx, t.pos));

    // Build a chrom index → name lookup from the header.
    let chrom_names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    // Group tuples into PileupColumn records.
    let mut columns: Vec<PileupColumn> = Vec::new();
    let mut i = 0;

    while i < tuples.len() {
        let chrom_idx = tuples[i].chrom_idx;
        let pos = tuples[i].pos;

        // Find the end of this (chrom, pos) group.
        let mut j = i;
        while j < tuples.len() && tuples[j].chrom_idx == chrom_idx && tuples[j].pos == pos {
            j += 1;
        }

        let chrom = chrom_names
            .get(chrom_idx as usize)
            .cloned()
            .unwrap_or_else(|| format!("ref{}", chrom_idx));

        // Apply max_depth cap (take the first max_depth entries).
        let end = j.min(i + max_depth);
        let bases: Vec<PileupBase> = tuples[i..end]
            .iter()
            .map(|t| PileupBase {
                base: t.base,
                base_qual: t.base_qual,
                mapq: t.mapq,
                is_rev: t.is_rev,
            })
            .collect();

        if !bases.is_empty() {
            columns.push(PileupColumn {
                chrom,
                pos,
                ref_base: b'N',
                bases,
            });
        }

        i = j;
    }

    Ok(columns)
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Normalize an IUPAC base byte to uppercase A/C/G/T, or N for anything else.
fn normalize_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'A',
        b'C' => b'C',
        b'G' => b'G',
        b'T' => b'T',
        _ => b'N',
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pileup_column_depth() {
        let col = PileupColumn {
            chrom: "chr1".into(),
            pos: 100,
            ref_base: b'A',
            bases: vec![
                PileupBase {
                    base: b'A',
                    base_qual: 30,
                    mapq: 60,
                    is_rev: false,
                },
                PileupBase {
                    base: b'G',
                    base_qual: 25,
                    mapq: 60,
                    is_rev: true,
                },
            ],
        };
        assert_eq!(col.depth(), 2);
    }
}
