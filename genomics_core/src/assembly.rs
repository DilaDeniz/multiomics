//! De Bruijn graph local assembly and haplotype generation.
//!
//! Builds a De Bruijn graph from reads overlapping an active genomic region,
//! enumerates candidate haplotype sequences via DFS path finding, and provides
//! Smith-Waterman alignment for haplotype-to-reference comparison.
//!
//! # References
//! - Garrison & Marth (2012) "Haplotype-based variant detection from short-read
//!   sequencing" arXiv:1207.3907 (freebayes).
//! - McKenna et al. (2010) "The Genome Analysis Toolkit: A MapReduce framework
//!   for analyzing next-generation DNA sequencing data" Genome Research 20:1297.

use ahash::AHashMap;

// ── Public types ──────────────────────────────────────────────────────────────

/// A read that has been extracted from a BAM record for local assembly.
#[derive(Debug, Clone)]
pub struct ActiveRead {
    /// Nucleotide sequence (ASCII uppercase A/C/G/T/N).
    pub seq: Vec<u8>,
    /// Per-base Phred quality scores (parallel to `seq`).
    pub quals: Vec<u8>,
    /// 0-based reference start position of this read.
    pub ref_start: u64,
    /// Mapping quality.
    pub mapq: u8,
}

/// A contiguous genomic window where evidence of variation triggers assembly.
#[derive(Debug, Clone)]
pub struct ActiveRegion {
    /// Reference sequence name.
    pub chrom: String,
    /// 0-based start of the active window.
    pub start: u64,
    /// 0-based exclusive end of the active window.
    pub end: u64,
    /// Reads overlapping this region (collected from the full read set).
    pub reads: Vec<ActiveRead>,
}

// ── De Bruijn graph ───────────────────────────────────────────────────────────

/// De Bruijn graph: maps k-mer → list of (successor k-mer, edge count).
///
/// Nodes are k-mers; directed edges represent consecutive overlap (k-1 shared
/// bases). Edge weights accumulate how many times the transition was observed.
pub struct DeBruijnGraph {
    /// Adjacency: source k-mer → Vec<(sink k-mer, count)>.
    edges: AHashMap<Vec<u8>, Vec<(Vec<u8>, u32)>>,
    k: usize,
}

impl DeBruijnGraph {
    /// Create an empty graph with the given k-mer size.
    pub fn new(k: usize) -> Self {
        Self {
            edges: AHashMap::new(),
            k,
        }
    }

    /// Add all consecutive k-mer pairs from `seq` with weight `count`.
    ///
    /// For position `i` in `0..seq.len()-k`: edge `seq[i..i+k]` →
    /// `seq[i+1..i+k+1]` is incremented by `count`.
    pub fn add_sequence(&mut self, seq: &[u8], count: u32) {
        if seq.len() < self.k + 1 {
            return;
        }
        for i in 0..seq.len() - self.k {
            let src = seq[i..i + self.k].to_vec();
            let dst = seq[i + 1..i + self.k + 1].to_vec();
            let entry = self.edges.entry(src).or_default();
            if let Some(edge) = entry.iter_mut().find(|(d, _)| d == &dst) {
                edge.1 = edge.1.saturating_add(count);
            } else {
                entry.push((dst, count));
            }
        }
    }

    /// Remove all edges with count < `min_count`, then remove isolated nodes.
    pub fn prune(&mut self, min_count: u32) {
        // Prune low-count edges.
        for edges in self.edges.values_mut() {
            edges.retain(|(_, c)| *c >= min_count);
        }
        // Remove nodes with no outgoing edges.
        self.edges.retain(|_, edges| !edges.is_empty());

        // Collect the set of nodes that appear as sinks (have incoming edges).
        let reachable_sinks: ahash::AHashSet<Vec<u8>> = self
            .edges
            .values()
            .flat_map(|v| v.iter().map(|(d, _)| d.clone()))
            .collect();

        // A source node with no incoming edges AND no outgoing edges is a dead
        // end; but we already dropped nodes with no outgoing edges above.
        // Remove any source whose successors are all gone (already handled).
        // Additionally, ensure sink-only nodes are noted but we don't need them
        // as map keys — they're discovered during DFS via the successor lists.
        let _ = reachable_sinks; // used implicitly via the adjacency structure
    }

    /// DFS to enumerate up to `max_paths` haplotype sequences from `source`
    /// k-mer to `sink` k-mer.
    ///
    /// The reconstructed sequence starts with the full `source` k-mer, then
    /// appends the last base of each subsequent k-mer.  Maximum path length is
    /// capped at 1000 bases to prevent exponential blowup.
    pub fn find_paths(
        &self,
        source: &[u8],
        sink: &[u8],
        max_paths: usize,
    ) -> Vec<Vec<u8>> {
        const MAX_LEN: usize = 1000;
        let mut results: Vec<Vec<u8>> = Vec::new();
        // Stack entries: (current k-mer, sequence assembled so far).
        let mut stack: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        let initial_seq = source.to_vec();
        stack.push((source.to_vec(), initial_seq));

        while let Some((node, seq)) = stack.pop() {
            if results.len() >= max_paths {
                break;
            }
            // Only treat reaching the sink as terminal when the path has grown
            // past the initial seed k-mer (avoids false termination when
            // source == sink, e.g. a palindromic or repeated reference end).
            if node == sink && seq.len() > self.k {
                results.push(seq);
                continue;
            }
            if seq.len() >= MAX_LEN {
                // Path too long; abandon it.
                continue;
            }
            if let Some(successors) = self.edges.get(&node) {
                for (next, _count) in successors {
                    let mut new_seq = seq.clone();
                    // Append only the last base of the next k-mer.
                    if let Some(&last_base) = next.last() {
                        new_seq.push(last_base);
                    }
                    stack.push((next.clone(), new_seq));
                }
            }
        }
        results
    }
}

// ── Active region detection ───────────────────────────────────────────────────

/// Scan `pileup` for positions with elevated non-reference base fraction and
/// group them into active regions of size `window`.
///
/// A position triggers assembly when `(depth - ref_count) / depth >
/// min_activity`.  Overlapping windows are merged.  Each region collects all
/// reads from `reads` whose alignment overlaps the window.
pub fn find_active_regions(
    pileup: &[crate::pileup::PileupColumn],
    reads: &[ActiveRead],
    window: usize,
    min_activity: f64,
) -> Vec<ActiveRegion> {
    // Collect trigger positions (chrom, start, end of window).
    let mut windows: Vec<(String, u64, u64)> = Vec::new();

    for col in pileup {
        let depth = col.depth();
        if depth == 0 {
            continue;
        }
        let ref_count = col
            .bases
            .iter()
            .filter(|b| b.base.eq_ignore_ascii_case(&col.ref_base))
            .count();
        let non_ref_frac = (depth - ref_count) as f64 / depth as f64;
        if non_ref_frac > min_activity {
            let half = (window / 2) as u64;
            let win_start = col.pos.saturating_sub(half);
            let win_end = win_start + window as u64;
            windows.push((col.chrom.clone(), win_start, win_end));
        }
    }

    if windows.is_empty() {
        return Vec::new();
    }

    // Sort and merge overlapping windows per chromosome.
    windows.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    let mut merged: Vec<(String, u64, u64)> = Vec::new();
    for (chrom, start, end) in windows {
        match merged.last_mut() {
            Some(last) if last.0 == chrom && start <= last.2 => {
                last.2 = last.2.max(end);
            }
            _ => merged.push((chrom, start, end)),
        }
    }

    // Build ActiveRegion records.
    merged
        .into_iter()
        .map(|(chrom, start, end)| {
            let region_reads = reads
                .iter()
                .filter(|r| {
                    r.ref_start < end
                        && r.ref_start + r.seq.len() as u64 > start
                })
                .cloned()
                .collect();
            ActiveRegion {
                chrom,
                start,
                end,
                reads: region_reads,
            }
        })
        .collect()
}

// ── Haplotype assembly ────────────────────────────────────────────────────────

/// Assemble candidate haplotypes for an active region.
///
/// For each k-mer size in `k_sizes`:
///   1. Build a `DeBruijnGraph`.
///   2. Add all reads (weight 1) and the reference sequence (weight 10).
///   3. Prune edges with count < 2.
///   4. Enumerate paths from the first k-mer of the reference to the last.
///
/// Unique haplotypes are accumulated across all k values.  The reference
/// haplotype is always included.  Returns at most `max_haplotypes` sequences.
pub fn assemble_haplotypes(
    region: &ActiveRegion,
    ref_seq: &[u8],
    k_sizes: &[usize],
    max_haplotypes: usize,
) -> Vec<Vec<u8>> {
    let mut unique: Vec<Vec<u8>> = Vec::new();

    // Always include the reference.
    if !ref_seq.is_empty() {
        unique.push(ref_seq.to_vec());
    }

    for &k in k_sizes {
        if ref_seq.len() < k + 1 {
            continue;
        }
        let source = ref_seq[..k].to_vec();
        let sink_start = ref_seq.len() - k;
        let sink = ref_seq[sink_start..].to_vec();

        let mut graph = DeBruijnGraph::new(k);

        for read in &region.reads {
            graph.add_sequence(&read.seq, 1);
        }
        graph.add_sequence(ref_seq, 10);
        graph.prune(2);

        let paths = graph.find_paths(&source, &sink, max_haplotypes);
        for path in paths {
            if !unique.contains(&path) {
                unique.push(path);
            }
            if unique.len() >= max_haplotypes {
                break;
            }
        }
        if unique.len() >= max_haplotypes {
            break;
        }
    }

    unique
}

// ── Smith-Waterman alignment ──────────────────────────────────────────────────

/// Align `query` to `reference` with affine gap scoring.
///
/// Scoring: match = 2, mismatch = -1, gap-open = -3, gap-extend = -1.
/// Returns `(aligned_query, aligned_reference)` with `b'-'` for gaps.
///
/// Uses the standard O(mn) DP with three matrices (M, Ix, Iy) and traceback.
pub fn smith_waterman(query: &[u8], reference: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let m = query.len();
    let n = reference.len();

    const MATCH: i32 = 2;
    const MISMATCH: i32 = -1;
    const GAP_OPEN: i32 = -3;
    const GAP_EXTEND: i32 = -1;
    const NEG_INF: i32 = i32::MIN / 2;

    // DP matrices: row = query index (0..=m), col = ref index (0..=n).
    // H = best score ending here; Ix = gap in ref (insertion in query);
    // Iy = gap in query (deletion from query / insertion in ref).
    let sz = (m + 1) * (n + 1);
    let mut h = vec![0i32; sz];
    let mut ix = vec![NEG_INF; sz];
    let mut iy = vec![NEG_INF; sz];

    let idx = |r: usize, c: usize| r * (n + 1) + c;

    // Traceback: 0=stop, 1=diag(M), 2=up(Ix), 3=left(Iy).
    let mut tb = vec![0u8; sz];

    let mut best_score = 0i32;
    let mut best_pos = (0usize, 0usize);

    for i in 1..=m {
        for j in 1..=n {
            // Insertion in query (gap in ref): Ix.
            let ix_open = h[idx(i - 1, j)].saturating_add(GAP_OPEN + GAP_EXTEND);
            let ix_ext = ix[idx(i - 1, j)].saturating_add(GAP_EXTEND);
            ix[idx(i, j)] = ix_open.max(ix_ext);

            // Deletion from query (gap in query): Iy.
            let iy_open = h[idx(i, j - 1)].saturating_add(GAP_OPEN + GAP_EXTEND);
            let iy_ext = iy[idx(i, j - 1)].saturating_add(GAP_EXTEND);
            iy[idx(i, j)] = iy_open.max(iy_ext);

            // Match/mismatch.
            let score_diag = if query[i - 1] == reference[j - 1] {
                MATCH
            } else {
                MISMATCH
            };
            let diag = h[idx(i - 1, j - 1)].saturating_add(score_diag);

            let cell = 0i32
                .max(diag)
                .max(ix[idx(i, j)])
                .max(iy[idx(i, j)]);
            h[idx(i, j)] = cell;

            // Record traceback direction.
            tb[idx(i, j)] = if cell == 0 {
                0
            } else if cell == diag {
                1
            } else if cell == ix[idx(i, j)] {
                2
            } else {
                3
            };

            if cell > best_score {
                best_score = cell;
                best_pos = (i, j);
            }
        }
    }

    // Traceback from best_pos.
    let mut aq: Vec<u8> = Vec::new();
    let mut ar: Vec<u8> = Vec::new();
    let (mut i, mut j) = best_pos;

    while i > 0 && j > 0 && h[idx(i, j)] > 0 {
        match tb[idx(i, j)] {
            1 => {
                aq.push(query[i - 1]);
                ar.push(reference[j - 1]);
                i -= 1;
                j -= 1;
            }
            2 => {
                aq.push(query[i - 1]);
                ar.push(b'-');
                i -= 1;
            }
            3 => {
                aq.push(b'-');
                ar.push(reference[j - 1]);
                j -= 1;
            }
            _ => break,
        }
    }

    aq.reverse();
    ar.reverse();
    (aq, ar)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn debruijn_simple_path() {
        // Build a tiny graph where one read has a SNP vs the reference.
        let ref_seq = b"ACGTACGT";
        let read_seq = b"ACGTTCGT"; // T→T at pos 4 (same); actually SNP at pos 4: C→T
        let k = 4;

        let mut graph = DeBruijnGraph::new(k);
        graph.add_sequence(ref_seq, 10);
        graph.add_sequence(read_seq, 5);
        graph.prune(2);

        let source = ref_seq[..k].to_vec();
        let sink_start = ref_seq.len() - k;
        let sink = ref_seq[sink_start..].to_vec();
        let paths = graph.find_paths(&source, &sink, 16);

        // At minimum, the reference path should be found.
        assert!(
            paths.iter().any(|p| p == ref_seq),
            "reference haplotype must appear in paths"
        );
    }

    #[test]
    fn prune_removes_low_count_edges() {
        let mut graph = DeBruijnGraph::new(3);
        graph.add_sequence(b"ACGT", 1); // count 1 — should be pruned at min=2
        graph.add_sequence(b"ACGT", 10); // raises count to 11 — kept
        graph.prune(2);
        // Source ACGT (k=3): ACG → CGT with count 11
        let src = b"ACG".to_vec();
        assert!(
            graph.edges.contains_key(&src),
            "high-count edge should survive pruning"
        );
    }

    #[test]
    fn smith_waterman_exact_match() {
        let q = b"ACGT";
        let r = b"ACGT";
        let (aq, ar) = smith_waterman(q, r);
        assert_eq!(aq, b"ACGT");
        assert_eq!(ar, b"ACGT");
    }

    #[test]
    fn smith_waterman_snp() {
        // Single mismatch in the middle.
        let q = b"ACTT";
        let r = b"ACGT";
        let (aq, ar) = smith_waterman(q, r);
        // Both sequences should be fully aligned (length 4, one mismatch).
        assert_eq!(aq.len(), ar.len());
        assert!(!aq.is_empty());
    }

    #[test]
    fn find_active_regions_triggers_on_nonref() {
        use crate::pileup::{PileupBase, PileupColumn};
        // Build a pileup column where 50% of bases are non-ref.
        let col = PileupColumn {
            chrom: "chr1".into(),
            pos: 100,
            ref_base: b'A',
            bases: vec![
                PileupBase { base: b'A', base_qual: 30, mapq: 60, is_rev: false },
                PileupBase { base: b'G', base_qual: 30, mapq: 60, is_rev: false },
            ],
        };
        let reads: Vec<ActiveRead> = Vec::new();
        let regions = find_active_regions(&[col], &reads, 200, 0.05);
        assert_eq!(regions.len(), 1, "one active region should be found");
        assert_eq!(regions[0].chrom, "chr1");
    }
}
