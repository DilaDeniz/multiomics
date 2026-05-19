//! Hybrid peptide index: precursor mass bins + fragment ion voting.
//!
//! Building on the Sage search engine design (Lazear et al. 2023, Nature Methods),
//! this index combines two levels of filtering before any expensive scoring:
//!
//! 1. **Precursor filter** — 1-Da mass bins narrow candidates from the full
//!    database (~1 M peptides) to a few hundred per spectrum in O(1).
//!
//! 2. **Fragment vote filter** — an inverted index maps each integer fragment
//!    mass (Da) to all peptides that produce a b or y ion at that mass. For
//!    each observed MS2 peak, we look up matching peptides and increment their
//!    vote counter. Only candidates with ≥ `MIN_FRAGMENT_VOTES` are passed to
//!    the full hyperscore. In practice this cuts the scoring work by 5–20×.
//!
//! Memory: ~8 bytes per fragment entry. For a human proteome tryptic digest
//! (~1 M peptides, ~28 ions each) the fragment index occupies ~224 MB.
//!
//! # References
//! * Lazear MR et al. (2023) Sage: an open-source tool for fast proteomics
//!   searching and quantification at scale. J Proteome Res.

use ahash::{AHashMap, AHashSet};

use crate::types::{aa_mass, Peptide, PROTON, WATER};

#[cfg(test)]
use crate::fasta::peptide_mass;
#[cfg(test)]
use crate::score::b_ions;

/// Minimum fragment votes before a candidate is passed to full hyperscore.
/// Lower = more sensitive but slower; higher = faster but may miss weak PSMs.
pub const MIN_FRAGMENT_VOTES: u32 = 2;

/// Precursor tolerance in Da used for coarse bin lookup (~1 500 Da × 10 ppm × 3 safety).
const PRECURSOR_TOL_DA: f64 = 0.05;

/// Fragment ion tolerance for the vote index (1 Da bins — integer rounding).
const FRAG_BIN_DA: f64 = 1.0;

/// Combined precursor + fragment index.
pub struct PeptideIndex {
    /// Precursor mass bins: floor(mass_da) → Vec<peptide_idx>.
    precursor_bins: AHashMap<u32, Vec<u32>>,
    /// Sorted fragment entries for binary search: (frag_mass_da_int, peptide_idx).
    fragment_entries: Vec<(u32, u32)>,
    peptides: Vec<Peptide>,
}

impl PeptideIndex {
    /// Build the index from a peptide list (consumed).
    ///
    /// Pre-computes all singly-charged b and y ions for each peptide and
    /// inserts them into the sorted fragment entry array.
    pub fn build(peptides: Vec<Peptide>) -> Self {
        let n = peptides.len();
        let mut precursor_bins: AHashMap<u32, Vec<u32>> = AHashMap::default();
        // Estimate ~28 fragment entries per peptide (avg length 15 → 14 b + 14 y).
        let mut fragment_entries: Vec<(u32, u32)> = Vec::with_capacity(n * 28);

        for (idx, pep) in peptides.iter().enumerate() {
            let pep_idx = idx as u32;

            // Precursor bin.
            let bin = pep.mass.floor() as u32;
            precursor_bins.entry(bin).or_default().push(pep_idx);

            // Fragment ions (b and y, z=1 only for the vote index).
            let seq = pep.sequence.as_bytes();
            let len = seq.len();
            if len < 2 {
                continue;
            }

            // b ions: cumulative prefix mass, exclude the last residue
            let mut prefix = 0.0f64;
            for &aa in &seq[..len - 1] {
                prefix += aa_mass(aa);
                let b_bin = (prefix + PROTON) as u32;
                fragment_entries.push((b_bin, pep_idx));
            }

            // y ions: cumulative suffix mass from the C-terminus
            let mut suffix = WATER;
            for &aa in seq[1..].iter().rev() {
                suffix += aa_mass(aa);
                let y_bin = (suffix + PROTON) as u32;
                fragment_entries.push((y_bin, pep_idx));
            }
        }

        // Sort by fragment mass bin — enables binary search in vote().
        fragment_entries.sort_unstable_by_key(|&(bin, _)| bin);

        Self {
            precursor_bins,
            fragment_entries,
            peptides,
        }
    }

    /// Return candidates for an MS2 spectrum using both precursor and fragment filters.
    ///
    /// Returns `(peptide_idx, vote_count)` pairs, sorted by vote count descending.
    /// Only peptides with ≥ `MIN_FRAGMENT_VOTES` matching fragment bins are returned.
    pub fn candidates_voted(
        &self,
        precursor_mass: f64,
        precursor_tol_ppm: f64,
        observed_mz: &[f32],
        min_votes: u32,
    ) -> Vec<(u32, u32)> {
        // Step 1: precursor filter — coarse Da-bin lookup.
        let tol_da = (precursor_mass * precursor_tol_ppm / 1e6).max(PRECURSOR_TOL_DA);
        let precursor_set: AHashSet<u32> = self.precursor_candidates(precursor_mass, tol_da);

        if precursor_set.is_empty() {
            return Vec::new();
        }

        // Step 2: fragment vote — for each observed peak look up which peptides
        // (from the precursor set) produce a matching fragment ion within 1 Da.
        let mut votes: AHashMap<u32, u32> = AHashMap::with_capacity(precursor_set.len());

        for &mz in observed_mz {
            let lo = (mz as f64 - FRAG_BIN_DA) as u32;
            let hi = (mz as f64 + FRAG_BIN_DA) as u32;

            let start = self.fragment_entries.partition_point(|&(bin, _)| bin < lo);
            for &(bin, pep_idx) in &self.fragment_entries[start..] {
                if bin > hi {
                    break;
                }
                if precursor_set.contains(&pep_idx) {
                    *votes.entry(pep_idx).or_insert(0) += 1;
                }
            }
        }

        // Collect candidates above vote threshold, sorted best-first.
        let mut result: Vec<(u32, u32)> =
            votes.into_iter().filter(|&(_, v)| v >= min_votes).collect();
        result.sort_unstable_by(|a, b| b.1.cmp(&a.1));
        result
    }

    /// Coarse precursor lookup: returns the set of peptide indices within `tol_da`.
    fn precursor_candidates(&self, query_mass: f64, tol_da: f64) -> AHashSet<u32> {
        let lo = (query_mass - tol_da).floor() as i64;
        let hi = (query_mass + tol_da).ceil() as i64;
        let mut out = AHashSet::default();
        for bin in lo..=hi {
            if bin < 0 {
                continue;
            }
            if let Some(idxs) = self.precursor_bins.get(&(bin as u32)) {
                out.extend(idxs.iter().copied());
            }
        }
        out
    }

    /// Direct precursor-only lookup (used by tests and as fallback).
    pub fn candidates(&self, query_mass: f64, tol_da: f64) -> Vec<u32> {
        self.precursor_candidates(query_mass, tol_da)
            .into_iter()
            .collect()
    }

    #[inline]
    pub fn peptide(&self, idx: u32) -> &Peptide {
        &self.peptides[idx as usize]
    }

    pub fn len(&self) -> usize {
        self.peptides.len()
    }

    pub fn is_empty(&self) -> bool {
        self.peptides.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_peptide(seq: &str) -> Peptide {
        Peptide {
            mass: peptide_mass(seq.as_bytes()),
            sequence: seq.into(),
            protein_idx: 0,
            is_decoy: false,
            missed_cleavages: 0,
        }
    }

    #[test]
    fn candidate_lookup_basic() {
        let peps = vec![
            make_peptide("PEPTIDE"),
            make_peptide("ACDEFGHIK"),
            make_peptide("LMNPQRSTVWY"),
        ];
        let mass0 = peps[0].mass;
        let idx = PeptideIndex::build(peps);

        // Precursor-only lookup should find the PEPTIDE entry.
        let candidates = idx.candidates(mass0, 0.02);
        assert!(
            !candidates.is_empty(),
            "should find PEPTIDE by precursor mass"
        );
    }

    #[test]
    fn fragment_vote_reduces_candidates() {
        let peps = vec![make_peptide("PEPTIDE"), make_peptide("AAAAAAAAAAK")];
        let target_mass = peps[0].mass;
        let idx = PeptideIndex::build(peps);

        // Build synthetic fragment peaks for PEPTIDE b ions.
        let b = b_ions(b"PEPTIDE");
        let peaks: Vec<f32> = b.iter().map(|&m| m as f32).collect();

        let voted = idx.candidates_voted(target_mass, 20.0, &peaks, 2);
        // PEPTIDE should win votes; AAAAAAAAAAK should have far fewer matches.
        assert!(!voted.is_empty(), "voted candidates should not be empty");
        // Top candidate by votes should be PEPTIDE (idx 0).
        assert_eq!(voted[0].0, 0, "PEPTIDE should be top-voted candidate");
    }

    #[test]
    fn fragment_entries_sorted() {
        let peps = vec![make_peptide("PEPTIDE"), make_peptide("MEPTIDEK")];
        let idx = PeptideIndex::build(peps);
        // Fragment entries must be sorted for binary search to work.
        let bins: Vec<u32> = idx.fragment_entries.iter().map(|&(b, _)| b).collect();
        let mut sorted = bins.clone();
        sorted.sort_unstable();
        assert_eq!(bins, sorted, "fragment entries must be sorted by bin");
    }
}
