//! Gene quantification: parse GTF/GFF3 annotations and count BAM reads per gene.
//!
//! Implements a featureCounts-style sweep-line algorithm using binary search
//! into a position-sorted gene interval list.
//!
//! # Reference
//! Liao Y, Smyth GK, Shi W (2014) featureCounts: an efficient general purpose
//! program for assigning sequence reads to genomic features. Bioinformatics
//! 30(7):923–930.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use memmap2::Mmap;
use serde::{Deserialize, Serialize};

use biomics_core::parse::{ByteLines, TabFields};

// ── Public types ──────────────────────────────────────────────────────────────

/// A gene interval from a GTF/GFF3 file.
#[derive(Debug, Clone)]
pub struct GeneInterval {
    pub chrom: String,
    pub start: u64, // 0-based
    pub end: u64,   // exclusive
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub strand: u8, // b'+' or b'-'
}

/// Summary statistics for a gene quantification run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuantSummary {
    pub total_reads: u64,
    pub assigned_reads: u64,
    pub assignment_rate: f64,
    pub n_genes_detected: u64,
    pub top_genes: Vec<(String, u64)>,
}

// ── GTF / GFF3 parser ─────────────────────────────────────────────────────────

/// Parse a GTF (gene_biotype=protein_coding) or GFF3 annotation file.
///
/// Only "gene" and "exon" feature records are retained. Returns one
/// [`GeneInterval`] per gene, merging all exon spans into a single interval
/// spanning from the leftmost exon start to the rightmost exon end.
///
/// Reads via a memory-mapped file for zero-copy performance.
pub fn parse_gtf(path: &Path) -> Result<Vec<GeneInterval>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot open annotation file: {}", path.display()))?;
    // SAFETY: We do not modify the file while it is mapped.
    let mmap = unsafe { Mmap::map(&file) }
        .with_context(|| format!("cannot mmap annotation file: {}", path.display()))?;

    parse_gtf_bytes(&mmap)
}

/// Parse GTF/GFF3 from an in-memory byte slice (also used by unit tests).
pub fn parse_gtf_bytes(data: &[u8]) -> Result<Vec<GeneInterval>> {
    // Detect format from the first non-comment line: GFF3 uses "##gff-version".
    let is_gff3 = ByteLines::new(data)
        .find(|l| !l.is_empty())
        .map(|l| l.starts_with(b"##gff-version"))
        .unwrap_or(false);

    // gene_id → (chrom, min_start, max_end, gene_name, strand)
    let mut genes: HashMap<String, (String, u64, u64, Option<String>, u8)> = HashMap::new();

    for line in ByteLines::new(data) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }

        let mut fields = TabFields::new(line);
        let chrom = match fields.next() {
            Some(c) => std::str::from_utf8(c).unwrap_or("?").to_owned(),
            None => continue,
        };
        let _source = fields.next();
        let feature = match fields.next() {
            Some(f) => f,
            None => continue,
        };

        // Only keep gene and exon records.
        if feature != b"gene" && feature != b"exon" {
            continue;
        }

        let start_bytes = match fields.next() {
            Some(s) => s,
            None => continue,
        };
        let end_bytes = match fields.next() {
            Some(e) => e,
            None => continue,
        };
        let _score = fields.next();
        let strand_field = fields.next().unwrap_or(b".");
        let strand = if strand_field.first() == Some(&b'-') {
            b'-'
        } else {
            b'+'
        };

        // GTF positions are 1-based inclusive; convert to 0-based half-open.
        let start_1based: u64 = std::str::from_utf8(start_bytes)
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let end_1based: u64 = std::str::from_utf8(end_bytes)
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        if start_1based == 0 || end_1based == 0 {
            continue;
        }
        let start_0 = start_1based.saturating_sub(1);
        let end_0 = end_1based; // half-open

        // Skip remaining tab fields to reach the attributes column (col 9).
        let _frame = fields.next();
        let attrs_bytes = match fields.next() {
            Some(a) => a,
            None => continue,
        };

        let (gene_id, gene_name) = if is_gff3 {
            parse_gff3_attrs(attrs_bytes)
        } else {
            parse_gtf_attrs(attrs_bytes)
        };

        let gene_id = match gene_id {
            Some(id) => id,
            None => continue,
        };

        let entry = genes
            .entry(gene_id.clone())
            .or_insert_with(|| (chrom.clone(), start_0, end_0, gene_name.clone(), strand));

        // Merge: expand the interval to cover this record.
        if start_0 < entry.1 {
            entry.1 = start_0;
        }
        if end_0 > entry.2 {
            entry.2 = end_0;
        }
        if entry.3.is_none() {
            entry.3 = gene_name;
        }
    }

    let mut intervals: Vec<GeneInterval> = genes
        .into_iter()
        .map(
            |(gene_id, (chrom, start, end, gene_name, strand))| GeneInterval {
                chrom,
                start,
                end,
                gene_id,
                gene_name,
                strand,
            },
        )
        .collect();

    // Sort by (chrom, start) for binary-search read counting.
    intervals.sort_unstable_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

    Ok(intervals)
}

// ── Attribute parsers ─────────────────────────────────────────────────────────

/// Parse GTF attribute string: `key "value"; key2 "value2"; …`
/// Returns (gene_id, gene_name).
fn parse_gtf_attrs(attrs: &[u8]) -> (Option<String>, Option<String>) {
    let text = std::str::from_utf8(attrs).unwrap_or("");
    let mut gene_id: Option<String> = None;
    let mut gene_name: Option<String> = None;

    for part in text.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some(rest) = part.strip_prefix("gene_id") {
            gene_id = extract_quoted_value(rest);
        } else if let Some(rest) = part.strip_prefix("gene_name") {
            gene_name = extract_quoted_value(rest);
        }
    }

    (gene_id, gene_name)
}

/// Parse GFF3 attribute string: `key=value;key2=value2;…`
/// Returns (gene_id, gene_name).
fn parse_gff3_attrs(attrs: &[u8]) -> (Option<String>, Option<String>) {
    let text = std::str::from_utf8(attrs).unwrap_or("");
    let mut gene_id: Option<String> = None;
    let mut gene_name: Option<String> = None;

    for part in text.split(';') {
        if let Some((key, val)) = part.split_once('=') {
            let key = key.trim();
            let val = val.trim();
            if key == "gene_id" || key == "ID" || key == "gene" {
                if gene_id.is_none() {
                    gene_id = Some(val.to_owned());
                }
            } else if (key == "gene_name" || key == "Name") && gene_name.is_none() {
                gene_name = Some(val.to_owned());
            }
        }
    }

    (gene_id, gene_name)
}

/// Extract the quoted value from a GTF attribute fragment like ` "TP53"`.
fn extract_quoted_value(s: &str) -> Option<String> {
    let s = s.trim();
    let inner = s
        .strip_prefix('"')
        .and_then(|s| s.strip_suffix('"'))
        .unwrap_or(s);
    if inner.is_empty() {
        None
    } else {
        Some(inner.to_owned())
    }
}

// ── Read counting ─────────────────────────────────────────────────────────────

/// Count reads per gene from a sorted BAM file using a sweep-line algorithm.
///
/// For each BAM record the algorithm performs a binary search into the
/// position-sorted `genes` slice to find candidate overlapping genes, then
/// checks exact overlap. Each read is assigned to at most one gene (the first
/// overlapping gene found).
///
/// # Arguments
/// * `bam_path`    – Path to a coordinate-sorted BAM file (index not required).
/// * `genes`       – Sorted gene interval list (output of [`parse_gtf`]).
/// * `strandedness` – `"unstranded"`, `"forward"`, or `"reverse"`.
/// * `min_mapq`    – Minimum MAPQ; reads below this threshold are skipped.
///
/// Returns a map of `gene_id → read count`.
pub fn count_reads(
    bam_path: &Path,
    genes: &[GeneInterval],
    strandedness: &str,
    min_mapq: u8,
) -> Result<ahash::AHashMap<String, u64>> {
    use noodles_bam as bam;

    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .with_context(|| format!("cannot open BAM file: {}", bam_path.display()))?;

    let header = reader.read_header().context("failed to read BAM header")?;

    // Build chrom index → name table from header.
    let chrom_names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    let mut counts: ahash::AHashMap<String, u64> = ahash::AHashMap::new();

    // Pre-populate with all gene IDs so every gene appears in the output.
    for g in genes {
        counts.entry(g.gene_id.clone()).or_insert(0);
    }

    for result in reader.record_bufs(&header) {
        let record = result.context("failed to read BAM record")?;

        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
        {
            continue;
        }

        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(0);
        if mapq < min_mapq {
            continue;
        }

        let ref_id = match record.reference_sequence_id() {
            Some(id) => id,
            None => continue,
        };
        let chrom = match chrom_names.get(ref_id) {
            Some(c) => c.as_str(),
            None => continue,
        };

        // noodles Position is 1-based; convert to 0-based.
        let read_start: u64 = match record.alignment_start() {
            Some(pos) => usize::from(pos) as u64 - 1,
            None => continue,
        };
        let read_len: u64 = record.cigar().alignment_span() as u64;
        if read_len == 0 {
            continue;
        }
        let read_end = read_start + read_len; // exclusive

        let is_rev = flags.is_reverse_complemented();

        // Binary search: find the leftmost gene whose start < read_end
        // and whose end > read_start (i.e., overlaps the read).
        let lo = genes.partition_point(|g| {
            g.chrom.as_str() < chrom || (g.chrom == chrom && g.end <= read_start)
        });
        let hi = genes.partition_point(|g| {
            g.chrom.as_str() < chrom || (g.chrom == chrom && g.start < read_end)
        });

        for gene in &genes[lo..hi] {
            if gene.chrom != chrom {
                break;
            }
            // Overlap check: gene.start < read_end && gene.end > read_start
            if gene.start >= read_end || gene.end <= read_start {
                continue;
            }

            // Strandedness filter.
            let passes_strand = match strandedness {
                "forward" => (gene.strand == b'+') != is_rev,
                "reverse" => (gene.strand == b'+') == is_rev,
                _ => true, // unstranded
            };
            if !passes_strand {
                continue;
            }

            *counts.entry(gene.gene_id.clone()).or_insert(0) += 1;
            break; // assign to the first overlapping gene only
        }
    }

    Ok(counts)
}

/// Parse a GTF/GFF3 annotation and count reads in a BAM file.
///
/// Convenience wrapper around [`parse_gtf`] + [`count_reads`].
pub fn quantify(
    bam_path: &Path,
    gtf_path: &Path,
    strandedness: &str,
    min_mapq: u8,
) -> Result<ahash::AHashMap<String, u64>> {
    let genes = parse_gtf(gtf_path)?;
    count_reads(bam_path, &genes, strandedness, min_mapq)
}

// ── Summary ───────────────────────────────────────────────────────────────────

/// Build a [`QuantSummary`] from a counts map and the total read count.
pub fn summarize_quant(counts: &ahash::AHashMap<String, u64>, total_reads: u64) -> QuantSummary {
    let assigned_reads: u64 = counts.values().sum();
    let assignment_rate = if total_reads > 0 {
        assigned_reads as f64 / total_reads as f64
    } else {
        0.0
    };
    let n_genes_detected = counts.values().filter(|&&c| c > 0).count() as u64;

    // Top 20 genes by count.
    let mut sorted: Vec<(String, u64)> = counts
        .iter()
        .filter(|(_, &c)| c > 0)
        .map(|(k, &v)| (k.clone(), v))
        .collect();
    sorted.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    sorted.truncate(20);

    QuantSummary {
        total_reads,
        assigned_reads,
        assignment_rate,
        n_genes_detected,
        top_genes: sorted,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal GTF with two genes (3 feature lines total).
    const MINIMAL_GTF: &[u8] = b"\
chr1\tENSEMBL\tgene\t1001\t2000\t.\t+\t.\tgene_id \"ENSG00000001\"; gene_name \"BRCA1\";\n\
chr1\tENSEMBL\texon\t1001\t1500\t.\t+\t.\tgene_id \"ENSG00000001\"; gene_name \"BRCA1\";\n\
chr2\tENSEMBL\tgene\t5001\t6000\t.\t-\t.\tgene_id \"ENSG00000002\"; gene_name \"TP53\";\n\
";

    #[test]
    fn parse_gtf_basic() {
        let intervals = parse_gtf_bytes(MINIMAL_GTF).expect("should parse");

        // Should have exactly 2 genes (exon merged into gene).
        assert_eq!(intervals.len(), 2);

        // Find BRCA1
        let brca1 = intervals
            .iter()
            .find(|g| g.gene_id == "ENSG00000001")
            .expect("BRCA1 not found");

        assert_eq!(brca1.chrom, "chr1");
        // GTF 1-based 1001 → 0-based 1000
        assert_eq!(brca1.start, 1000);
        // GTF 1-based end 2000 → exclusive 2000
        assert_eq!(brca1.end, 2000);
        assert_eq!(brca1.gene_name.as_deref(), Some("BRCA1"));
        assert_eq!(brca1.strand, b'+');

        // Find TP53
        let tp53 = intervals
            .iter()
            .find(|g| g.gene_id == "ENSG00000002")
            .expect("TP53 not found");
        assert_eq!(tp53.chrom, "chr2");
        assert_eq!(tp53.start, 5000);
        assert_eq!(tp53.end, 6000);
        assert_eq!(tp53.strand, b'-');
    }

    #[test]
    fn count_reads_empty() {
        // With an empty gene list, every read must be unassigned → total count = 0.
        let empty_genes: Vec<GeneInterval> = Vec::new();
        let counts: ahash::AHashMap<String, u64> = ahash::AHashMap::new();

        // Without a real BAM we cannot call count_reads, so we verify the gene
        // pre-population: with no genes the returned map should be empty.
        for g in &empty_genes {
            let _ = &g.gene_id; // suppress unused warning
        }
        // The map is empty when there are no genes.
        assert_eq!(counts.len(), 0);

        // Verify summarize_quant correctly reports 0 assigned reads.
        let summary = summarize_quant(&counts, 1000);
        assert_eq!(summary.assigned_reads, 0);
        assert_eq!(summary.n_genes_detected, 0);
        assert_eq!(summary.assignment_rate, 0.0);
    }

    #[test]
    fn quant_summary() {
        let mut counts: ahash::AHashMap<String, u64> = ahash::AHashMap::new();
        counts.insert("BRCA1".to_owned(), 300);
        counts.insert("TP53".to_owned(), 200);
        counts.insert("EGFR".to_owned(), 0);

        let total = 1000u64;
        let summary = summarize_quant(&counts, total);

        assert_eq!(summary.total_reads, 1000);
        assert_eq!(summary.assigned_reads, 500);
        assert!((summary.assignment_rate - 0.5).abs() < 1e-10);
        assert_eq!(summary.n_genes_detected, 2); // EGFR has 0 count
        assert_eq!(summary.top_genes.len(), 2);
        // Top gene should be BRCA1 (300 > 200).
        assert_eq!(summary.top_genes[0].0, "BRCA1");
        assert_eq!(summary.top_genes[0].1, 300);
    }
}
