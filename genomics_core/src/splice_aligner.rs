//! Splice-aware RNA-seq aligner inspired by STAR (Dobin et al. 2013).
//!
//! Algorithm overview (two-pass):
//! Pass 1 – Seed finding and unspliced chaining; gaps > 20 bp in reference space
//!           but 0 bp in query space are treated as candidate introns.
//!           Candidate introns validated by GT..AG (or GC..AG) donor/acceptor motifs.
//! Pass 2 – Re-align using expanded junction database (known + novel from pass 1).
//!           Smith-Waterman per exon segment; introns get a fixed GT-AG bonus.
//!           CIGAR assembled with M (match) and N (intron skip) operators.

use ahash::AHashMap;
use anyhow::Result;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

// ── Scoring constants ──────────────────────────────────────────────────────────

const MATCH_SCORE: i32 = 2;
const MISMATCH_SCORE: i32 = -1;
const GAP_OPEN: i32 = -3;
const GAP_EXTEND: i32 = -1;
const SW_BAND: usize = 20;
/// Bonus added to the alignment score for each validated GT-AG intron.
const GTAG_BONUS: i32 = 5;
/// Minimum intron length to consider a reference-space gap as a candidate intron.
const MIN_INTRON: usize = 20;
/// k-mer size for junction seeds.
const SPLICE_K: usize = 16;

// ── Public types ───────────────────────────────────────────────────────────────

/// A splice junction (intron boundary) discovered or loaded from annotation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct SpliceJunction {
    pub chrom: String,
    pub donor: u64,    // 5′ splice site (last exon base, 0-based)
    pub acceptor: u64, // 3′ splice site (first exon base, 0-based)
    pub strand: u8,    // b'+', b'-', or b'.'
    pub novel: bool,   // true = discovered in pass 1; false = from GTF
}

/// Result of aligning one RNA-seq read.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpliceAlignment {
    pub chrom: String,
    pub ref_start: u64,
    pub ref_end: u64,
    /// CIGAR with N operators for introns, e.g. "50M100000N50M"
    pub cigar: String,
    pub score: i32,
    pub mapq: u8,
    pub is_reverse: bool,
    /// Junctions spanned by this read.
    pub junctions: Vec<SpliceJunction>,
    pub n_mismatches: u32,
}

// ── Junction index ─────────────────────────────────────────────────────────────

pub struct SpliceIndex {
    /// Reference genome sequence (flattened chromosomes).
    genome: Vec<u8>,
    chrom_offsets: Vec<(String, usize)>,
    /// k-mer index of the genome (k=SPLICE_K for junction seeds).
    kmer_index: AHashMap<u64, Vec<u32>>,
    k: usize,
    /// Known and discovered junctions keyed by (chrom_idx, donor).
    junctions: AHashMap<(u32, u64), Vec<SpliceJunction>>,
}

impl SpliceIndex {
    /// Build index from a FASTA reference and optional GTF junction list.
    /// Pass `gtf_junctions = &[]` for de novo discovery.
    pub fn build(fasta_bytes: &[u8], gtf_junctions: &[SpliceJunction]) -> Result<Self> {
        let chroms = parse_fasta(fasta_bytes);
        anyhow::ensure!(!chroms.is_empty(), "No chromosomes found in FASTA input");

        let mut genome: Vec<u8> = Vec::new();
        let mut chrom_offsets: Vec<(String, usize)> = Vec::new();

        for (name, seq) in &chroms {
            chrom_offsets.push((name.clone(), genome.len()));
            genome.extend_from_slice(seq);
        }

        let k = SPLICE_K;
        let mask: u64 = if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };

        let mut kmer_index: AHashMap<u64, Vec<u32>> = AHashMap::new();

        for (_name, start_offset) in &chrom_offsets {
            let offset = *start_offset;
            let chrom_seq = &genome[offset..];
            let chrom_len = chrom_seq.len();
            if chrom_len < k {
                continue;
            }

            let first_hash = kmer_hash_dna(&chrom_seq[..k], k);
            let mut valid_window = first_hash.is_some();
            let mut hash = first_hash.unwrap_or(0);
            let mut bad_run = if valid_window { 0usize } else { k };

            if valid_window {
                kmer_index.entry(hash).or_default().push(offset as u32);
            }

            for i in 1..=(chrom_len - k) {
                let new_base = chrom_seq[i + k - 1];
                if base_bits(new_base).is_none() {
                    bad_run = k;
                    valid_window = false;
                    continue;
                }
                if bad_run > 0 {
                    bad_run -= 1;
                    if bad_run == 0 {
                        hash = match kmer_hash_dna(&chrom_seq[i..i + k], k) {
                            Some(h) => h,
                            None => {
                                bad_run = k;
                                valid_window = false;
                                continue;
                            }
                        };
                        valid_window = true;
                    } else {
                        continue;
                    }
                } else {
                    match kmer_roll(hash, new_base, mask) {
                        Some(h) => {
                            hash = h;
                            valid_window = true;
                        }
                        None => {
                            bad_run = k;
                            valid_window = false;
                            continue;
                        }
                    }
                }

                if valid_window {
                    kmer_index
                        .entry(hash)
                        .or_default()
                        .push((offset + i) as u32);
                }
            }
        }

        // Populate junction table from GTF-provided junctions.
        let mut junctions: AHashMap<(u32, u64), Vec<SpliceJunction>> = AHashMap::new();
        for junc in gtf_junctions {
            let chrom_idx = chrom_offsets
                .iter()
                .position(|(n, _)| n == &junc.chrom)
                .unwrap_or(0) as u32;
            junctions
                .entry((chrom_idx, junc.donor))
                .or_default()
                .push(junc.clone());
        }

        Ok(Self {
            genome,
            chrom_offsets,
            kmer_index,
            k,
            junctions,
        })
    }

    /// Parse splice junctions from a GTF file (exon features only).
    /// Adjacent exons on the same transcript define a junction.
    pub fn junctions_from_gtf(gtf_bytes: &[u8]) -> Vec<SpliceJunction> {
        // Parse exon records keyed by transcript_id → sorted list of (start, end, chrom, strand).
        let mut transcripts: AHashMap<String, Vec<(u64, u64, String, u8)>> = AHashMap::new();

        for raw_line in gtf_bytes.split(|&b| b == b'\n') {
            let line = std::str::from_utf8(raw_line).unwrap_or("").trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.splitn(9, '\t').collect();
            if fields.len() < 9 {
                continue;
            }
            let feature = fields[2];
            if feature != "exon" {
                continue;
            }
            let chrom = fields[0].to_owned();
            let start: u64 = fields[3].parse::<u64>().unwrap_or(0).saturating_sub(1); // GTF is 1-based
            let end: u64 = fields[4].parse::<u64>().unwrap_or(0).saturating_sub(1);
            let strand_ch = fields[6].as_bytes().first().copied().unwrap_or(b'.');
            let strand = match strand_ch {
                b'+' => b'+',
                b'-' => b'-',
                _ => b'.',
            };
            let attrs = fields[8];
            // Extract transcript_id "..." from attributes.
            let tx_id = if let Some(pos) = attrs.find("transcript_id") {
                let rest = &attrs[pos + 13..].trim_start_matches(' ');
                let rest = rest.trim_start_matches('"');
                let end_q = rest.find('"').unwrap_or(rest.len());
                rest[..end_q].to_owned()
            } else {
                continue;
            };
            transcripts
                .entry(tx_id)
                .or_default()
                .push((start, end, chrom, strand));
        }

        let mut junctions: Vec<SpliceJunction> = Vec::new();
        for (_tx, mut exons) in transcripts {
            exons.sort_unstable_by_key(|e| e.0);
            for i in 0..exons.len().saturating_sub(1) {
                let (_, end_a, ref chrom_a, strand_a) = exons[i];
                let (start_b, _, ref chrom_b, _) = exons[i + 1];
                if chrom_a != chrom_b {
                    continue;
                }
                if start_b <= end_a {
                    continue; // overlapping exons — skip
                }
                junctions.push(SpliceJunction {
                    chrom: chrom_a.clone(),
                    donor: end_a,
                    acceptor: start_b,
                    strand: strand_a,
                    novel: false,
                });
            }
        }
        junctions
    }

    /// Two-pass alignment:
    /// Pass 1: unspliced alignment → discover novel junctions from gap evidence
    /// Pass 2: re-align using expanded junction database
    pub fn align_rna(&self, read: &[u8]) -> Option<SpliceAlignment> {
        let (fwd, fwd_rev) = self.align_rna_strand(read, false);
        let rc = rev_comp(read);
        let (rev, rev_rev) = self.align_rna_strand(&rc, true);

        // Choose best alignment by score.
        let candidates = [fwd, rev];
        let second_candidates = [fwd_rev, rev_rev];

        let best = candidates.into_iter().flatten().max_by_key(|a| a.score)?;

        let second_best_score = second_candidates
            .into_iter()
            .flatten()
            .map(|a| a.score)
            .max();

        let mapq = splice_mapq(best.score, second_best_score);
        Some(SpliceAlignment { mapq, ..best })
    }

    /// Parallel batch alignment.
    pub fn align_batch_rna(&self, reads: &[Vec<u8>]) -> Vec<Option<SpliceAlignment>> {
        reads.par_iter().map(|r| self.align_rna(r)).collect()
    }

    // ── Internal alignment logic ───────────────────────────────────────────────

    /// Returns (primary_alignment, second_best_score_holder) for one strand.
    fn align_rna_strand(
        &self,
        read: &[u8],
        is_reverse: bool,
    ) -> (Option<SpliceAlignment>, Option<SpliceAlignment>) {
        let k = self.k;
        if read.len() < k {
            return (None, None);
        }

        // ── Pass 1: collect seeds ─────────────────────────────────────────────
        let seeds = self.collect_seeds(read);
        if seeds.is_empty() {
            return (None, None);
        }

        // ── Pass 1: chain seeds; detect candidate introns ─────────────────────
        let chain = chain_seeds_splice(&seeds, k);
        if chain.is_empty() {
            return (None, None);
        }

        // Detect candidate introns from chain gaps.
        let mut novel_junctions: Vec<SpliceJunction> = Vec::new();
        let chrom_idx = self.global_chrom_idx(chain[0].rpos as usize) as u32;
        let chrom_name = self.chrom_offsets[chrom_idx as usize].0.clone();

        for w in chain.windows(2) {
            let (a, b) = (w[0], w[1]);
            let ref_gap = (b.rpos as i64) - (a.rpos as i64 + k as i64);
            let qry_gap = (b.qpos as i64) - (a.qpos as i64 + k as i64);
            // Candidate intron: large ref gap, tiny query gap.
            if ref_gap > MIN_INTRON as i64 && qry_gap.abs() <= 2 {
                let donor_pos = a.rpos as usize + k; // first base of intron
                let acceptor_pos = b.rpos as usize; // first base of exon after intron
                if is_gt_ag(&self.genome, donor_pos, acceptor_pos) {
                    novel_junctions.push(SpliceJunction {
                        chrom: chrom_name.clone(),
                        donor: (donor_pos.saturating_sub(1)) as u64,
                        acceptor: acceptor_pos as u64,
                        strand: b'+', // GT-AG strand assignment (simplified)
                        novel: true,
                    });
                }
            }
        }

        // ── Pass 2: splice-aware alignment ────────────────────────────────────
        // Build combined junction set (known + novel).
        let mut local_junctions: Vec<SpliceJunction> = novel_junctions;
        if let Some(known) = self.junctions.get(&(chrom_idx, 0)) {
            // This key doesn't exist; iterate whole map for this chrom.
            local_junctions.extend_from_slice(known);
        }
        // Collect all known junctions for this chromosome.
        for ((cidx, _), jvec) in &self.junctions {
            if *cidx == chrom_idx {
                for j in jvec {
                    if !local_junctions.contains(j) {
                        local_junctions.push(j.clone());
                    }
                }
            }
        }

        // Determine exon segments from the chain.
        let (segments, intron_lengths) = chain_to_exon_segments(&chain, k, read.len());

        if segments.is_empty() {
            return (None, None);
        }

        // Align each exon segment with banded SW.
        let chrom_offset = self.chrom_offsets[chrom_idx as usize].1;
        let mut total_score = 0i32;
        let mut cigar_parts: Vec<String> = Vec::new();
        let mut spanned_junctions: Vec<SpliceJunction> = Vec::new();
        let mut n_mismatches = 0u32;
        let mut aln_ref_start: Option<u64> = None;
        let mut aln_ref_end: u64 = 0;

        for (seg_idx, seg) in segments.iter().enumerate() {
            let q_slice = &read[seg.q_start..seg.q_end];
            if q_slice.is_empty() {
                continue;
            }
            let r_global_start = chrom_offset + seg.r_win_start;
            let r_global_end = (chrom_offset + seg.r_win_end).min(self.genome.len());
            if r_global_start >= self.genome.len() || r_global_start >= r_global_end {
                continue;
            }
            let r_slice = &self.genome[r_global_start..r_global_end];

            let (score, cigar, sw_r_start, sw_r_end) =
                banded_smith_waterman(q_slice, r_slice, SW_BAND);
            if score <= 0 {
                continue;
            }
            total_score += score;
            cigar_parts.push(cigar.clone());

            let actual_r_start = seg.r_win_start + sw_r_start;
            let actual_r_end = seg.r_win_start + sw_r_end;

            if aln_ref_start.is_none() {
                aln_ref_start = Some(actual_r_start as u64);
            }
            aln_ref_end = actual_r_end as u64;

            // Count mismatches from the aligned segment.
            let mm_ref_end = (r_global_start + sw_r_end).min(self.genome.len());
            if r_global_start + sw_r_start < mm_ref_end {
                n_mismatches += count_mismatches(
                    q_slice,
                    &self.genome[r_global_start + sw_r_start..mm_ref_end],
                    &cigar,
                );
            }

            // If this isn't the last segment, add an N operator for the intron.
            if seg_idx < intron_lengths.len() {
                let intron_len = intron_lengths[seg_idx];
                if intron_len > 0 {
                    cigar_parts.push(format!("{}N", intron_len));

                    // Intron donor = seed end of current exon; acceptor = seed start of next exon.
                    let donor = seg.r_seed_end.saturating_sub(1) as u64;
                    let acceptor = (seg.r_seed_end + intron_len) as u64;

                    // GT-AG bonus.
                    let donor_global = chrom_offset + seg.r_seed_end;
                    let acceptor_global = chrom_offset + seg.r_seed_end + intron_len;
                    if is_gt_ag(&self.genome, donor_global, acceptor_global) {
                        total_score += GTAG_BONUS;
                    }

                    // Find matching junction.
                    let found_junc = local_junctions
                        .iter()
                        .find(|j| j.donor == donor && j.acceptor == acceptor)
                        .cloned()
                        .unwrap_or_else(|| SpliceJunction {
                            chrom: chrom_name.clone(),
                            donor,
                            acceptor,
                            strand: b'.',
                            novel: true,
                        });
                    spanned_junctions.push(found_junc);
                }
            }
        }

        if aln_ref_start.is_none() || total_score <= 0 {
            return (None, None);
        }

        let combined_cigar = cigar_parts.concat();

        let aln = SpliceAlignment {
            chrom: chrom_name,
            ref_start: aln_ref_start.unwrap(),
            ref_end: aln_ref_end,
            cigar: combined_cigar,
            score: total_score,
            mapq: 0, // filled by caller
            is_reverse,
            junctions: spanned_junctions,
            n_mismatches,
        };

        (Some(aln), None)
    }

    /// Collect k-mer seeds for a read.
    fn collect_seeds(&self, read: &[u8]) -> Vec<Seed> {
        let k = self.k;
        let mask: u64 = if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };

        let mut seeds: Vec<Seed> = Vec::new();
        let mut hash = kmer_hash_dna(&read[..k], k);

        if let Some(h) = hash {
            if let Some(positions) = self.kmer_index.get(&h) {
                for &rpos in positions.iter().take(50) {
                    seeds.push(Seed { qpos: 0, rpos });
                }
            }
        }

        for i in 1..=(read.len() - k) {
            hash = match hash {
                Some(h) => kmer_roll(h, read[i + k - 1], mask),
                None => kmer_hash_dna(&read[i..i + k], k),
            };
            if let Some(h) = hash {
                if let Some(positions) = self.kmer_index.get(&h) {
                    for &rpos in positions.iter().take(50) {
                        seeds.push(Seed {
                            qpos: i as u32,
                            rpos,
                        });
                    }
                }
            }
        }

        seeds
    }

    /// Return the chromosome index for a global position.
    fn global_chrom_idx(&self, global: usize) -> usize {
        let idx = self.chrom_offsets.partition_point(|o| o.1 <= global);
        if idx == 0 {
            0
        } else {
            idx - 1
        }
    }
}

// ── Helpers ────────────────────────────────────────────────────────────────────

/// Check for GT..AG (or GC..AG) donor-acceptor motifs.
///
/// `donor_pos` = first base of the intron (genome index).
/// `acceptor_pos` = first base of the downstream exon (genome index).
fn is_gt_ag(genome: &[u8], donor_pos: usize, acceptor_pos: usize) -> bool {
    if donor_pos + 2 > genome.len() || acceptor_pos < 2 || acceptor_pos > genome.len() {
        return false;
    }
    let d1 = genome[donor_pos].to_ascii_uppercase();
    let d2 = genome[donor_pos + 1].to_ascii_uppercase();
    let a1 = genome[acceptor_pos - 2].to_ascii_uppercase();
    let a2 = genome[acceptor_pos - 1].to_ascii_uppercase();

    let donor_ok = d1 == b'G' && (d2 == b'T' || d2 == b'C');
    let acceptor_ok = a1 == b'A' && a2 == b'G';
    donor_ok && acceptor_ok
}

/// Encode a k-mer into a u64 using 2-bit encoding.
pub(crate) fn kmer_hash_dna(seq: &[u8], k: usize) -> Option<u64> {
    if seq.len() < k {
        return None;
    }
    let mut h: u64 = 0;
    for &b in seq[..k].iter() {
        let bits = base_bits(b)?;
        h = (h << 2) | bits;
    }
    Some(h)
}

/// 2-bit encoding of a nucleotide.
#[inline]
fn base_bits(b: u8) -> Option<u64> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Rolling k-mer hash update.
#[inline]
fn kmer_roll(prev: u64, new_base: u8, mask: u64) -> Option<u64> {
    let bits = base_bits(new_base)?;
    Some(((prev << 2) | bits) & mask)
}

/// Parse a (possibly multi-line) FASTA byte stream.
fn parse_fasta(bytes: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut result: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for raw_line in bytes.split(|&b| b == b'\n') {
        let line = if raw_line.last() == Some(&b'\r') {
            &raw_line[..raw_line.len() - 1]
        } else {
            raw_line
        };
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(name) = current_name.take() {
                result.push((name, std::mem::take(&mut current_seq)));
            }
            let header = &line[1..];
            let name_end = header
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .unwrap_or(header.len());
            let name = String::from_utf8_lossy(&header[..name_end]).into_owned();
            current_name = Some(name);
            current_seq = Vec::new();
        } else if current_name.is_some() {
            for &b in line {
                let ub = b.to_ascii_uppercase();
                if matches!(ub, b'A' | b'C' | b'G' | b'T') {
                    current_seq.push(ub);
                } else {
                    current_seq.push(b'N');
                }
            }
        }
    }
    if let Some(name) = current_name {
        result.push((name, current_seq));
    }
    result
}

/// Reverse complement of a DNA sequence.
fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[inline]
fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

/// Compute MAPQ from best and (optionally) second-best score.
fn splice_mapq(best_score: i32, second_score: Option<i32>) -> u8 {
    match second_score {
        None => 60,
        Some(s) => {
            let diff = (best_score - s).max(0) as f64;
            ((diff * 10.0 / 2.0_f64.ln()) as u8).min(60)
        }
    }
}

// ── Seed / chain types ─────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct Seed {
    qpos: u32,
    rpos: u32,
}

/// Chain seeds, allowing large ref gaps (introns) while keeping query-collinear.
fn chain_seeds_splice(seeds: &[Seed], k: usize) -> Vec<Seed> {
    if seeds.is_empty() {
        return vec![];
    }

    let mut sorted = seeds.to_vec();
    sorted.sort_unstable_by_key(|s| (s.rpos, s.qpos));

    let n = sorted.len();
    let mut dp = vec![k as i64; n];
    let mut prev = vec![usize::MAX; n];

    for j in 1..n {
        for i in 0..j {
            if sorted[j].rpos <= sorted[i].rpos {
                continue;
            }
            if sorted[j].qpos <= sorted[i].qpos {
                continue;
            }
            // Gap sizes: how many bases are skipped between the END of seed i and START of seed j.
            let dr_gap = (sorted[j].rpos as i64) - (sorted[i].rpos as i64 + k as i64);
            let dq_gap = (sorted[j].qpos as i64) - (sorted[i].qpos as i64 + k as i64);

            // Allow large ref gaps (introns) but keep query gap small.
            // A valid intron: large ref gap, effectively zero query gap.
            let is_intron_gap = dr_gap > MIN_INTRON as i64 && dq_gap.abs() <= 2;

            if !is_intron_gap {
                // For normal (non-intron) extension, require diagonal consistency.
                if dq_gap < 0 || dq_gap > k as i64 * 4 {
                    continue;
                }
                if (dr_gap - dq_gap).abs() > k as i64 * 2 {
                    continue;
                }
            }

            let score = dp[i] + k as i64;
            if score > dp[j] {
                dp[j] = score;
                prev[j] = i;
            }
        }
    }

    let best_idx = match dp
        .iter()
        .enumerate()
        .max_by_key(|&(_, &s)| s)
        .map(|(i, _)| i)
    {
        Some(i) => i,
        None => return vec![],
    };

    let mut chain = Vec::new();
    let mut cur = best_idx;
    loop {
        chain.push(sorted[cur]);
        if prev[cur] == usize::MAX {
            break;
        }
        cur = prev[cur];
    }
    chain.reverse();
    chain
}

/// An exon segment: query range and ref SW-window in chrom-relative coordinates.
/// `r_seed_start` / `r_seed_end` are the tight seed-derived boundaries (before padding).
struct ExonSegment {
    q_start: usize,
    q_end: usize,
    /// Padded reference window for Smith-Waterman.
    r_win_start: usize,
    r_win_end: usize,
    /// Seed-derived reference boundary (tight, no padding) — used for intron length calc.
    _r_seed_start: usize,
    r_seed_end: usize,
}

/// Convert a seed chain to a list of exon segments with intron gap sizes.
///
/// Returns `(segments, intron_lengths)` where `intron_lengths[i]` is the number
/// of reference bases skipped between `segments[i]` and `segments[i+1]`.
fn chain_to_exon_segments(
    chain: &[Seed],
    k: usize,
    read_len: usize,
) -> (Vec<ExonSegment>, Vec<usize>) {
    if chain.is_empty() {
        return (vec![], vec![]);
    }

    let mut segments: Vec<ExonSegment> = Vec::new();
    let mut intron_lengths: Vec<usize> = Vec::new();
    let mut block_start = 0usize;

    for i in 1..=chain.len() {
        let is_last = i == chain.len();
        let is_split = if is_last {
            true
        } else {
            let dr = (chain[i].rpos as i64) - (chain[i - 1].rpos as i64 + k as i64);
            let dq = (chain[i].qpos as i64) - (chain[i - 1].qpos as i64 + k as i64);
            dr > MIN_INTRON as i64 && dq.abs() <= 2
        };

        if is_split {
            let first_seed = chain[block_start];
            let last_seed = chain[i - 1];

            let q_start = first_seed.qpos as usize;
            let q_end = (last_seed.qpos as usize + k).min(read_len);

            // Tight seed boundaries.
            let r_seed_start = first_seed.rpos as usize;
            let r_seed_end = last_seed.rpos as usize + k;

            // Padded SW window: expand by unaligned read flanks.
            let r_win_start = r_seed_start.saturating_sub(q_start);
            let r_win_end = r_seed_end + (read_len - q_end);

            // Intron length = distance from this segment's seed end to next segment's seed start.
            if !is_last {
                let next_r_seed_start = chain[i].rpos as usize;
                let intron = next_r_seed_start.saturating_sub(r_seed_end);
                intron_lengths.push(intron);
            }

            segments.push(ExonSegment {
                q_start,
                q_end,
                r_win_start,
                r_win_end,
                _r_seed_start: r_seed_start,
                r_seed_end,
            });

            block_start = i;
        }
    }

    // No intron detected: whole read is one segment.
    if segments.is_empty() {
        let first = chain[0];
        let last = chain[chain.len() - 1];
        let q_start = first.qpos as usize;
        let q_end = (last.qpos as usize + k).min(read_len);
        let r_seed_start = first.rpos as usize;
        let r_seed_end = last.rpos as usize + k;
        let r_win_start = r_seed_start.saturating_sub(q_start);
        let r_win_end = r_seed_end + (read_len - q_end);
        segments.push(ExonSegment {
            q_start: 0,
            q_end: read_len,
            r_win_start,
            r_win_end,
            _r_seed_start: r_seed_start,
            r_seed_end,
        });
    }

    (segments, intron_lengths)
}

// ── Banded Smith-Waterman ──────────────────────────────────────────────────────

/// Banded Smith-Waterman; returns (score, cigar, ref_start_in_window, ref_end_in_window).
fn banded_smith_waterman(
    query: &[u8],
    reference: &[u8],
    band: usize,
) -> (i32, String, usize, usize) {
    let qlen = query.len();
    let rlen = reference.len();

    if qlen == 0 || rlen == 0 {
        return (0, String::new(), 0, 0);
    }

    let cols = rlen + 1;
    let rows = qlen + 1;

    let mut h = vec![0i32; rows * cols];
    let mut e = vec![i32::MIN / 2; rows * cols];
    let mut f = vec![i32::MIN / 2; rows * cols];
    let mut tb = vec![0u8; rows * cols];

    let idx = |i: usize, j: usize| i * cols + j;

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    #[allow(clippy::needless_range_loop)]
    for i in 1..rows {
        let j_lo = if i > band { i - band } else { 1 };
        let j_hi = (i + band).min(rlen);

        for j in j_lo..=j_hi {
            let diag = h[idx(i - 1, j - 1)];
            let score_m = if query[i - 1].eq_ignore_ascii_case(&reference[j - 1]) {
                diag + MATCH_SCORE
            } else {
                diag + MISMATCH_SCORE
            };

            let e_open = h[idx(i - 1, j)].saturating_add(GAP_OPEN + GAP_EXTEND);
            let e_ext = e[idx(i - 1, j)].saturating_add(GAP_EXTEND);
            e[idx(i, j)] = e_open.max(e_ext);

            let f_open = h[idx(i, j - 1)].saturating_add(GAP_OPEN + GAP_EXTEND);
            let f_ext = f[idx(i, j - 1)].saturating_add(GAP_EXTEND);
            f[idx(i, j)] = f_open.max(f_ext);

            let cell = score_m.max(e[idx(i, j)]).max(f[idx(i, j)]).max(0);
            h[idx(i, j)] = cell;

            tb[idx(i, j)] = if cell == 0 {
                0
            } else if cell == score_m {
                1
            } else if cell == e[idx(i, j)] {
                2
            } else {
                3
            };

            if cell > best_score {
                best_score = cell;
                best_i = i;
                best_j = j;
            }
        }
    }

    if best_score <= 0 {
        return (0, String::new(), 0, 0);
    }

    let ref_end_in_window = best_j;
    let (cigar, ref_start_1based) = traceback(&tb, &h, rows, cols, best_i, best_j);
    let ref_start_in_window = ref_start_1based.saturating_sub(1);

    (best_score, cigar, ref_start_in_window, ref_end_in_window)
}

/// SW traceback → CIGAR string + 1-based ref start.
fn traceback(
    tb: &[u8],
    h: &[i32],
    _rows: usize,
    cols: usize,
    start_i: usize,
    start_j: usize,
) -> (String, usize) {
    let idx = |i: usize, j: usize| i * cols + j;
    let mut ops: Vec<(char, u32)> = Vec::new();
    let mut i = start_i;
    let mut j = start_j;

    while i > 0 && j > 0 && h[idx(i, j)] > 0 {
        match tb[idx(i, j)] {
            1 => {
                push_op(&mut ops, 'M');
                i -= 1;
                j -= 1;
            }
            2 => {
                push_op(&mut ops, 'D');
                i -= 1;
            }
            3 => {
                push_op(&mut ops, 'I');
                j -= 1;
            }
            _ => break,
        }
    }

    ops.reverse();
    let cigar = ops
        .iter()
        .map(|(op, count)| format!("{}{}", count, op))
        .collect();
    (cigar, j + 1)
}

#[inline]
fn push_op(ops: &mut Vec<(char, u32)>, op: char) {
    if let Some(last) = ops.last_mut() {
        if last.0 == op {
            last.1 += 1;
            return;
        }
    }
    ops.push((op, 1));
}

/// Count mismatches between a query slice and aligned reference bases (M ops only).
fn count_mismatches(query: &[u8], reference: &[u8], cigar: &str) -> u32 {
    let mut mismatches = 0u32;
    let mut qi = 0usize;
    let mut ri = 0usize;
    let mut num_buf = String::new();

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_buf.push(ch);
        } else {
            let n: usize = num_buf.parse().unwrap_or(0);
            num_buf.clear();
            match ch {
                'M' => {
                    for _ in 0..n {
                        if qi < query.len() && ri < reference.len() {
                            if !query[qi].eq_ignore_ascii_case(&reference[ri]) {
                                mismatches += 1;
                            }
                            qi += 1;
                            ri += 1;
                        }
                    }
                }
                'I' => qi += n,
                'D' => ri += n,
                _ => {}
            }
        }
    }
    mismatches
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a FASTA with a known intron structure.
    ///
    /// Layout: exon1 (200 bp) | intron (1000 bp, starts GT, ends AG) | exon2 (200 bp)
    fn make_spliced_fasta() -> (Vec<u8>, Vec<u8>, usize, usize) {
        // Generate distinct pseudo-random (but deterministic) exon sequences.
        // Using an LCG avoids any external RNG dependency; different seeds produce
        // sequences that share no k=16 substrings with each other or with the intron.
        fn lcg_seq(seed: u64, len: usize) -> Vec<u8> {
            let bases: &[u8] = b"ACGT";
            let mut s = seed;
            (0..len)
                .map(|_| {
                    s = s
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    bases[((s >> 33) as usize) % 4]
                })
                .collect()
        }
        let exon1 = lcg_seq(0xDEAD_BEEF_1234_5678, 200);
        let exon2 = lcg_seq(0xCAFE_BABE_8765_4321, 200);

        // Intron: starts with GT, ends with AG, rest is filler.
        let mut intron = Vec::with_capacity(1000);
        intron.extend_from_slice(b"GT");
        for i in 2..998 {
            intron.push(b"TTTTCCCC"[i % 8]);
        }
        intron.extend_from_slice(b"AG");
        assert_eq!(intron.len(), 1000);

        let exon1_len = exon1.len();
        let exon2_start = exon1_len + intron.len();

        let mut genome = Vec::new();
        genome.extend_from_slice(&exon1);
        genome.extend_from_slice(&intron);
        genome.extend_from_slice(&exon2);

        let mut fasta = b">chr1\n".to_vec();
        fasta.extend_from_slice(&genome);
        fasta.push(b'\n');

        (fasta, genome, exon1_len, exon2_start)
    }

    // ── splice_junction_gtag_check ─────────────────────────────────────────────

    #[test]
    fn splice_junction_gtag_check() {
        // Construct a tiny genome with a GT..AG intron at positions 10-19.
        // Genome: ACGTACGTAC | GT XXXXXX AG | ACGTACGTAC
        //         0          10             20
        let mut genome = b"ACGTACGTAC".to_vec(); // 10 bp exon
        genome.extend_from_slice(b"GTTTTTTTAG"); // 10 bp intron (GT..AG)
        genome.extend_from_slice(b"ACGTACGTAC"); // 10 bp exon

        // donor_pos = 10 (first base of intron), acceptor_pos = 20 (first base of exon 2)
        assert!(is_gt_ag(&genome, 10, 20), "should detect GT-AG motif");

        // Negative case: GA..AG (not a valid donor).
        let mut bad = b"ACGTACGTAC".to_vec();
        bad.extend_from_slice(b"GATTTTTTAG");
        bad.extend_from_slice(b"ACGTACGTAC");
        assert!(!is_gt_ag(&bad, 10, 20), "GA donor should not be detected");
    }

    // ── align_unspliced_read ───────────────────────────────────────────────────

    #[test]
    fn align_unspliced_read() {
        let (fasta, genome, _exon1_len, _exon2_start) = make_spliced_fasta();
        let idx = SpliceIndex::build(&fasta, &[]).expect("build index");

        // Read fully within exon1 (position 50..100).
        let read = genome[50..100].to_vec();
        let aln = idx.align_rna(&read).expect("should align");

        // No intron → no N in CIGAR.
        assert!(
            !aln.cigar.contains('N'),
            "unspliced read should have no N in CIGAR, got: {}",
            aln.cigar
        );
        assert!(aln.score > 0, "score should be positive");
    }

    // ── align_spliced_read ─────────────────────────────────────────────────────

    #[test]
    fn align_spliced_read() {
        let (fasta, genome, exon1_len, exon2_start) = make_spliced_fasta();
        let idx = SpliceIndex::build(&fasta, &[]).expect("build index");

        // 100 bp read: 50 bp from end of exon1 + 50 bp from start of exon2.
        let mut read = Vec::new();
        read.extend_from_slice(&genome[exon1_len - 50..exon1_len]);
        read.extend_from_slice(&genome[exon2_start..exon2_start + 50]);
        assert_eq!(read.len(), 100);

        let aln = idx.align_rna(&read).expect("should align spliced read");

        // Must contain an N operator.
        assert!(
            aln.cigar.contains('N'),
            "spliced read should have N in CIGAR, got: {}",
            aln.cigar
        );

        // The intron span encoded in the N op should be ~1000.
        let n_len = extract_n_length(&aln.cigar);
        assert!(
            n_len >= 900 && n_len <= 1100,
            "N length should be ~1000, got {n_len} (CIGAR: {})",
            aln.cigar
        );
    }

    // ── junction_batch_parallel ────────────────────────────────────────────────

    #[test]
    fn junction_batch_parallel() {
        let (fasta, genome, _exon1_len, _exon2_start) = make_spliced_fasta();
        let idx = SpliceIndex::build(&fasta, &[]).expect("build index");

        // 10 reads, each 30 bp, fully within exon1 at various offsets.
        let reads: Vec<Vec<u8>> = (0..10)
            .map(|i| genome[i * 5..i * 5 + 30].to_vec())
            .collect();

        let results = idx.align_batch_rna(&reads);
        assert_eq!(results.len(), 10, "should return 10 results");

        let aligned = results.iter().filter(|r| r.is_some()).count();
        assert!(
            aligned >= 8,
            "at least 8/10 unspliced reads should align, got {aligned}/10"
        );
    }

    // ── Helper: extract the length from the first N operator in a CIGAR ────────

    fn extract_n_length(cigar: &str) -> u64 {
        let mut num_buf = String::new();
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                num_buf.push(ch);
            } else {
                let n: u64 = num_buf.parse().unwrap_or(0);
                num_buf.clear();
                if ch == 'N' {
                    return n;
                }
            }
        }
        0
    }
}
