//! Seed-chain-extend short-read aligner inspired by BWA-MEM (Li 2013).
//!
//! Algorithm overview:
//! 1. Build a k-mer hash index over the reference genome.
//! 2. For each read, collect k-mer seed hits from the index.
//! 3. Chain seeds with a collinear DP.
//! 4. Extend the best chain with banded Smith-Waterman.
//! 5. Compute MAPQ from the ratio of best vs second-best alignment score.

use ahash::AHashMap;
use anyhow::{Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

// ── Scoring constants ─────────────────────────────────────────────────────────

const MATCH_SCORE: i32 = 2;
const MISMATCH_SCORE: i32 = -1;
const GAP_OPEN: i32 = -3;
const GAP_EXTEND: i32 = -1;
const SW_BAND: usize = 20;

// ── Default k-mer size ────────────────────────────────────────────────────────

const DEFAULT_K: usize = 19;

// ── Types ─────────────────────────────────────────────────────────────────────

/// A reference index built from a FASTA-format reference.
pub struct ReferenceIndex {
    /// Flattened sequence of all chromosomes (uppercase ASCII).
    sequence: Vec<u8>,
    /// Start offset of each chromosome in `sequence`.
    chrom_offsets: Vec<(String, usize)>, // (name, offset)
    /// k-mer → list of global positions in the reference.
    kmer_index: AHashMap<u64, Vec<u32>>,
    k: usize,
}

/// A single alignment result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Alignment {
    pub chrom: String,
    pub ref_start: u64,
    pub ref_end: u64,
    pub query_start: usize,
    pub query_end: usize,
    /// Number of matching bases.
    pub n_matches: u32,
    /// CIGAR string, e.g. "50M" or "10M2D38M".
    pub cigar: String,
    /// Alignment score (SW score).
    pub score: i32,
    /// Mapping quality (PHRED-scaled, capped at 60).
    pub mapq: u8,
    pub is_reverse: bool,
}

// ── k-mer hash ────────────────────────────────────────────────────────────────

/// Encode a k-mer into a u64 using 2-bit encoding.
///
/// Returns `None` if any character is not in `{A, C, G, T}` (upper or lower).
fn kmer_hash(seq: &[u8], k: usize) -> Option<u64> {
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

/// Map a nucleotide byte to its 2-bit representation.
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

/// Rolling update: slide the k-mer window one position to the right.
///
/// `mask` must equal `(1u64 << (2*k)) - 1`.
#[inline]
fn kmer_hash_roll(prev: u64, new_base: u8, mask: u64) -> Option<u64> {
    let bits = base_bits(new_base)?;
    Some(((prev << 2) | bits) & mask)
}

/// Return `true` if the k-mer is low-complexity (>50 % dominated by one base).
fn is_low_complexity(seq: &[u8], k: usize) -> bool {
    let mut counts = [0u32; 4];
    for &b in seq[..k].iter() {
        match b.to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    let threshold = (k / 2) as u32; // > 50 %
    counts.iter().any(|&c| c > threshold)
}

// ── FASTA parser ──────────────────────────────────────────────────────────────

/// Parse a (possibly multi-line) FASTA byte stream.
///
/// Returns a vector of `(chromosome_name, sequence_bytes)`.
fn parse_fasta(bytes: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut result: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in bytes.split(|&b| b == b'\n') {
        // Trim carriage returns.
        let line = if line.last() == Some(&b'\r') {
            &line[..line.len() - 1]
        } else {
            line
        };
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            // Flush previous record.
            if let Some(name) = current_name.take() {
                result.push((name, std::mem::take(&mut current_seq)));
            }
            // Parse chromosome name: everything up to the first whitespace.
            let header = &line[1..];
            let name_end = header
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .unwrap_or(header.len());
            let name = String::from_utf8_lossy(&header[..name_end]).into_owned();
            current_name = Some(name);
            current_seq = Vec::new();
        } else if current_name.is_some() {
            // Sequence line — keep only ACGT (uppercase).
            for &b in line {
                let ub = b.to_ascii_uppercase();
                if matches!(ub, b'A' | b'C' | b'G' | b'T') {
                    current_seq.push(ub);
                } else {
                    // Ambiguous base — push 'N' as sentinel but we'll skip it during hashing.
                    current_seq.push(b'N');
                }
            }
        }
    }
    // Flush final record.
    if let Some(name) = current_name {
        result.push((name, current_seq));
    }
    result
}

// ── Index construction ────────────────────────────────────────────────────────

impl ReferenceIndex {
    /// Build a k-mer index from a reference genome already loaded in memory.
    ///
    /// `k = 19` is a good default for 150 bp reads.
    pub fn build(fasta_bytes: &[u8], k: usize) -> Result<Self> {
        let chroms = parse_fasta(fasta_bytes);
        if chroms.is_empty() {
            anyhow::bail!("No chromosomes found in FASTA input");
        }

        let mut sequence: Vec<u8> = Vec::new();
        let mut chrom_offsets: Vec<(String, usize)> = Vec::new();

        for (name, seq) in &chroms {
            chrom_offsets.push((name.clone(), sequence.len()));
            sequence.extend_from_slice(seq);
        }

        let mask: u64 = if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };

        let mut kmer_index: AHashMap<u64, Vec<u32>> = AHashMap::new();

        // Index all chromosomes.
        for (_name, start_offset) in &chrom_offsets {
            let offset = *start_offset;
            let chrom_seq = &sequence[offset..];
            let chrom_len = chrom_seq.len();
            if chrom_len < k {
                continue;
            }

            // Compute hash of the very first k-mer.
            let first_hash = kmer_hash(&chrom_seq[..k], k);
            // Track whether the current window is valid (all ACGT).
            let mut valid_window = first_hash.is_some();
            // Use 0 as a sentinel for "invalid" — guarded by `valid_window`.
            let mut hash = first_hash.unwrap_or(0);

            if valid_window && !is_low_complexity(&chrom_seq[..k], k) {
                let global_pos = offset as u32; // position 0 of this chrom
                kmer_index.entry(hash).or_default().push(global_pos);
            }

            let mut bad_run = 0usize; // how many consecutive bad bases we've seen
            if !valid_window {
                bad_run = k; // pretend we just saw k bad bases
            }

            for i in 1..=(chrom_len - k) {
                let new_base = chrom_seq[i + k - 1];
                if base_bits(new_base).is_none() {
                    // Non-ACGT base — invalidate window.
                    bad_run = k;
                    valid_window = false;
                    // We'll need to rebuild the hash once we have k good bases again.
                    continue;
                }
                if bad_run > 0 {
                    bad_run -= 1;
                    if bad_run == 0 {
                        // Window is good again; recompute hash from scratch.
                        hash = match kmer_hash(&chrom_seq[i..i + k], k) {
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
                    // Normal rolling update.
                    match kmer_hash_roll(hash, new_base, mask) {
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

                if valid_window && !is_low_complexity(&chrom_seq[i..i + k], k) {
                    let global_pos = (offset + i) as u32;
                    kmer_index.entry(hash).or_default().push(global_pos);
                }
            }
        }

        Ok(Self {
            sequence,
            chrom_offsets,
            kmer_index,
            k,
        })
    }

    /// Build a k-mer index from a FASTA file at `path` (memory-maps it).
    pub fn from_path(path: &std::path::Path) -> Result<Self> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("Cannot open reference FASTA: {}", path.display()))?;
        // Safety: we do not mutate the file while it is mapped.
        let mmap = unsafe { memmap2::Mmap::map(&file) }
            .with_context(|| format!("Cannot mmap: {}", path.display()))?;
        Self::build(&mmap[..], DEFAULT_K)
    }
}

// ── Alignment ─────────────────────────────────────────────────────────────────

/// An internal seed hit: (query_pos, ref_global_pos).
#[derive(Clone, Copy)]
struct Seed {
    qpos: u32,
    rpos: u32,
}

impl ReferenceIndex {
    /// Align a single read to the reference.  Returns the best alignment or `None`.
    pub fn align(&self, read: &[u8]) -> Option<Alignment> {
        let (best, second_best) = self.align_both_strands(read);
        let mapq = compute_mapq(best.as_ref(), second_best.as_ref());
        best.map(|mut a| {
            a.mapq = mapq;
            a
        })
    }

    /// Align all reads in parallel using Rayon.
    pub fn align_batch(&self, reads: &[Vec<u8>]) -> Vec<Option<Alignment>> {
        reads.par_iter().map(|r| self.align(r)).collect()
    }

    // ── Internal helpers ──────────────────────────────────────────────────────

    fn align_both_strands(&self, read: &[u8]) -> (Option<Alignment>, Option<Alignment>) {
        let fwd = self.align_one_strand(read, false);
        let rc = rev_comp(read);
        let rev = self.align_one_strand(&rc, true);

        // Pick best by score.
        let (best, second) = match (&fwd, &rev) {
            (Some(f), Some(r)) => {
                if f.score >= r.score {
                    (fwd.clone(), rev)
                } else {
                    (rev.clone(), fwd)
                }
            }
            (Some(_), None) => (fwd, rev),
            (None, Some(_)) => (rev, fwd),
            (None, None) => (None, None),
        };
        (best, second)
    }

    fn align_one_strand(&self, read: &[u8], is_reverse: bool) -> Option<Alignment> {
        let k = self.k;
        if read.len() < k {
            return None;
        }

        // ── 1. Collect seeds ──────────────────────────────────────────────────
        let mut seeds: Vec<Seed> = Vec::new();
        let mask: u64 = if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };

        // Compute initial hash.
        let mut hash = kmer_hash(&read[..k], k);
        if let Some(h) = hash {
            if let Some(positions) = self.kmer_index.get(&h) {
                for &rpos in positions.iter().take(50) {
                    // cap hits per seed
                    seeds.push(Seed { qpos: 0, rpos });
                }
            }
        }

        for i in 1..=(read.len() - k) {
            hash = match hash {
                Some(h) => kmer_hash_roll(h, read[i + k - 1], mask),
                None => kmer_hash(&read[i..i + k], k),
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

        if seeds.is_empty() {
            return None;
        }

        // ── 2. Chain seeds ────────────────────────────────────────────────────
        let chain = chain_seeds(&seeds, k);
        if chain.is_empty() {
            return None;
        }

        // Determine reference window from the chain.
        let first = chain[0];
        let last = chain[chain.len() - 1];

        // Query span covered by the chain.
        let q_start = first.qpos as usize;
        let q_end = (last.qpos as usize) + k;

        // Reference window: pad by query-start on the left and leftover on the right.
        let ref_chain_start = first.rpos as usize;
        let ref_chain_end = (last.rpos as usize) + k;

        // Extend the reference window to cover the full read.
        let ref_start = ref_chain_start.saturating_sub(q_start);
        let overhang_right = read.len().saturating_sub(q_end);
        let ref_end_raw = ref_chain_end + overhang_right;
        let ref_end = ref_end_raw.min(self.sequence.len());

        if ref_start >= self.sequence.len() || ref_start >= ref_end {
            return None;
        }

        let ref_slice = &self.sequence[ref_start..ref_end];

        // ── 3. Banded Smith-Waterman ──────────────────────────────────────────
        let (score, cigar, sw_ref_start, sw_ref_end) =
            banded_smith_waterman(read, ref_slice, SW_BAND);
        if score <= 0 {
            return None;
        }

        // Translate SW window-relative positions to global positions.
        let aln_global_start = ref_start + sw_ref_start;
        let aln_global_end = ref_start + sw_ref_end;

        // Count matches from CIGAR.
        let n_matches = count_matches_from_cigar(&cigar);

        // Chromosome lookup.
        let (chrom, chrom_start) = global_to_chrom(&self.chrom_offsets, aln_global_start);

        let chrom_ref_start = (aln_global_start - chrom_start) as u64;
        let chrom_ref_end = (aln_global_end - chrom_start) as u64;

        Some(Alignment {
            chrom,
            ref_start: chrom_ref_start,
            ref_end: chrom_ref_end,
            query_start: 0,
            query_end: read.len(),
            n_matches,
            cigar,
            score,
            mapq: 0, // filled in by caller
            is_reverse,
        })
    }
}

/// Collinear chaining of seeds.
///
/// DP: `dp[i]` = best chain score ending at seed `i`.
/// Transition: `dp[j] = max(dp[i] + k)` for all `i < j` where
/// `seeds[j].rpos > seeds[i].rpos` and the "diagonal" is roughly preserved.
fn chain_seeds(seeds: &[Seed], k: usize) -> Vec<Seed> {
    if seeds.is_empty() {
        return vec![];
    }

    // Sort by reference position then query position.
    let mut sorted = seeds.to_vec();
    sorted.sort_unstable_by_key(|s| (s.rpos, s.qpos));

    let n = sorted.len();
    let mut dp = vec![k as i64; n]; // score of chain ending here
    let mut prev = vec![usize::MAX; n]; // back-pointer

    for j in 1..n {
        for i in 0..j {
            // Must advance on both ref and query.
            if sorted[j].rpos <= sorted[i].rpos {
                continue;
            }
            if sorted[j].qpos <= sorted[i].qpos {
                continue;
            }
            // Diagonal consistency: |Δref - Δquery| <= k (allow small gaps/indels)
            let dr = (sorted[j].rpos - sorted[i].rpos) as i64;
            let dq = (sorted[j].qpos - sorted[i].qpos) as i64;
            if (dr - dq).abs() > k as i64 * 2 {
                continue;
            }
            let score = dp[i] + k as i64;
            if score > dp[j] {
                dp[j] = score;
                prev[j] = i;
            }
        }
    }

    // Traceback best chain.
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

// ── Banded Smith-Waterman ─────────────────────────────────────────────────────

/// Banded Smith-Waterman alignment.
///
/// Uses affine gap penalties: match = +2, mismatch = −1, gap_open = −3, gap_extend = −1.
/// The band is centred on the main diagonal.
///
/// Returns `(score, cigar_string, ref_start_in_window, ref_end_in_window)`.
/// `ref_start_in_window` and `ref_end_in_window` are 0-based half-open offsets
/// into the `reference` slice.
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

    // Three DP matrices: H (best), E (gap in ref / deletion), F (gap in query / insertion).
    let mut h = vec![0i32; rows * cols];
    let mut e = vec![i32::MIN / 2; rows * cols];
    let mut f = vec![i32::MIN / 2; rows * cols];

    // Track-back matrix: 0 = none, 1 = match/mismatch (diag), 2 = deletion (up), 3 = insertion (left).
    let mut tb = vec![0u8; rows * cols];

    let idx = |i: usize, j: usize| i * cols + j;

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    #[allow(clippy::needless_range_loop)]
    for i in 1..rows {
        // Band limits on j (reference axis).
        let j_lo = if i > band { i - band } else { 1 };
        let j_hi = (i + band).min(rlen);

        for j in j_lo..=j_hi {
            // Match / mismatch.
            let diag = h[idx(i - 1, j - 1)];
            let score_m = if query[i - 1].eq_ignore_ascii_case(&reference[j - 1]) {
                diag + MATCH_SCORE
            } else {
                diag + MISMATCH_SCORE
            };

            // Gap in reference (deletion in query).
            let e_open = h[idx(i - 1, j)].saturating_add(GAP_OPEN + GAP_EXTEND);
            let e_ext = e[idx(i - 1, j)].saturating_add(GAP_EXTEND);
            e[idx(i, j)] = e_open.max(e_ext);

            // Gap in query (insertion in query).
            let f_open = h[idx(i, j - 1)].saturating_add(GAP_OPEN + GAP_EXTEND);
            let f_ext = f[idx(i, j - 1)].saturating_add(GAP_EXTEND);
            f[idx(i, j)] = f_open.max(f_ext);

            let cell = score_m.max(e[idx(i, j)]).max(f[idx(i, j)]).max(0);
            h[idx(i, j)] = cell;

            // Track back source.
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

    // ref_end (exclusive) = best_j (1-based → convert to 0-based end = best_j).
    let ref_end_in_window = best_j;

    // Traceback: returns the cigar and the ref start position (1-based).
    let (cigar, ref_start_1based) = traceback(&tb, &h, rows, cols, best_i, best_j);
    // Convert to 0-based.
    let ref_start_in_window = ref_start_1based.saturating_sub(1);

    (best_score, cigar, ref_start_in_window, ref_end_in_window)
}

/// Traceback the SW DP to produce a CIGAR string and return the 1-based reference start.
///
/// Returns `(cigar_string, ref_start_1based)`.
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
    // j is now the 0-based index just before the alignment start in the reference;
    // the 1-based start = j + 1.
    (cigar, j + 1)
}

/// Push an operation onto a run-length encoded CIGAR list.
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

// ── Utilities ─────────────────────────────────────────────────────────────────

/// Reverse complement of a DNA sequence (ASCII).
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

/// Convert a global reference position to `(chrom_name, chrom_offset)`.
fn global_to_chrom(offsets: &[(String, usize)], global: usize) -> (String, usize) {
    // Binary search for the last chromosome whose offset <= global.
    let idx = offsets.partition_point(|o| o.1 <= global);
    if idx == 0 {
        return (offsets[0].0.clone(), 0);
    }
    let idx = idx - 1;
    (offsets[idx].0.clone(), offsets[idx].1)
}

/// Count 'M' operations in a CIGAR string.
fn count_matches_from_cigar(cigar: &str) -> u32 {
    let mut total = 0u32;
    let mut num_buf = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_buf.push(ch);
        } else {
            let n: u32 = num_buf.parse().unwrap_or(0);
            num_buf.clear();
            if ch == 'M' {
                total += n;
            }
        }
    }
    total
}

/// Compute MAPQ from the best and second-best alignment scores.
///
/// MAPQ = −10 · log10(p_wrong), capped at 60.
fn compute_mapq(best: Option<&Alignment>, second: Option<&Alignment>) -> u8 {
    match (best, second) {
        (None, _) => 0,
        (Some(_), None) => 60,
        (Some(b), Some(s)) => {
            if b.score <= 0 {
                return 0;
            }
            // Probability of wrong mapping ~ exp(s2 - s1).
            let diff = (b.score - s.score).max(0) as f64;
            let mapq = (diff * 10.0 / 2.0_f64.ln()).min(60.0) as u8; // rough heuristic
            mapq.min(60)
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a simple 100-bp reference FASTA.
    fn make_ref_fasta(seq: &[u8]) -> Vec<u8> {
        let mut out = b">chr1\n".to_vec();
        out.extend_from_slice(seq);
        out.push(b'\n');
        out
    }

    /// Build a simple repetitive-free reference of length `len`.
    fn acgt_cycle(len: usize) -> Vec<u8> {
        let bases = b"ACGTACGTACGTACGTACGT";
        (0..len).map(|i| bases[i % bases.len()]).collect()
    }

    // ── kmer_hash_consistency ─────────────────────────────────────────────────

    #[test]
    fn kmer_hash_consistency() {
        let k = 10;
        let seq = b"ACGTACGTAC";
        let h1 = kmer_hash(seq, k);
        let h2 = kmer_hash(seq, k);
        assert!(h1.is_some());
        assert_eq!(h1, h2, "same k-mer must always hash to the same value");

        // Reverse complement should generally differ.
        let rc = rev_comp(seq);
        let h_rc = kmer_hash(&rc, k);
        assert!(h_rc.is_some());
        // RC of ACGTACGTAC is GTACGTACGT — different bit pattern.
        assert_ne!(h1, h_rc, "reverse complement should hash differently");
    }

    // ── align_exact_match ─────────────────────────────────────────────────────

    #[test]
    fn align_exact_match() {
        // Use a 200bp reference with k=10 so the 20bp read has seeds that map.
        // Also make the reference longer than one period to reduce ambiguity.
        let ref_seq = acgt_cycle(200);
        let fasta = make_ref_fasta(&ref_seq);
        // Use a smaller k so that a 20bp read yields more than 1 seed.
        let index = ReferenceIndex::build(&fasta, 10).expect("build index");

        // 20-bp read from position 5.
        let read = ref_seq[5..25].to_vec();
        let aln = index.align(&read).expect("should align");

        // Read is 20 bp, max SW score with match=+2 is 40.
        assert_eq!(aln.score, 40, "exact match score should be 40");
        // The alignment should cover exactly 20 bases in the reference.
        assert_eq!(
            aln.ref_end - aln.ref_start,
            20,
            "alignment should span 20 reference bases"
        );
        assert!(!aln.is_reverse, "should align on forward strand");
    }

    // ── align_with_mismatch ───────────────────────────────────────────────────

    #[test]
    fn align_with_mismatch() {
        // Use a 50bp read so k=19 seeds from outside the mismatch position can be found.
        let ref_seq = acgt_cycle(200);
        let fasta = make_ref_fasta(&ref_seq);
        let index = ReferenceIndex::build(&fasta, DEFAULT_K).expect("build index");

        // Take a 50bp read from position 0.
        let mut read = ref_seq[0..50].to_vec();
        // Introduce one mismatch near the start (position 2), leaving 19+ bp seeds at
        // positions 3..22, 4..23, … that don't contain the mismatch.
        read[2] = if read[2] == b'A' { b'T' } else { b'A' };

        let aln = index
            .align(&read)
            .expect("should still align with 1 mismatch");
        // Score should be less than 100 (50 × +2) but positive.
        assert!(aln.score > 0, "score should be positive");
        assert!(aln.score < 100, "score should be less than perfect match");
    }

    // ── align_batch_parallel ──────────────────────────────────────────────────

    #[test]
    fn align_batch_parallel() {
        let ref_seq = acgt_cycle(200);
        let fasta = make_ref_fasta(&ref_seq);
        let index = ReferenceIndex::build(&fasta, DEFAULT_K).expect("build index");

        // 10 reads, each 20 bp, from different positions.
        let reads: Vec<Vec<u8>> = (0..10)
            .map(|i| ref_seq[i * 5..i * 5 + 20].to_vec())
            .collect();

        let results = index.align_batch(&reads);
        assert_eq!(results.len(), 10);

        let aligned_count = results.iter().filter(|r| r.is_some()).count();
        assert_eq!(
            aligned_count, 10,
            "all 10 reads should align (got {aligned_count}/10)"
        );
    }
}
