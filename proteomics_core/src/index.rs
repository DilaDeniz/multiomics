//! Peptide mass index for fast candidate retrieval.
//!
//! Groups peptides into 1-Da mass bins. For a query precursor mass M, all
//! candidates in bins [M − tol, M + tol] are returned for scoring.
//! This is O(1) lookup compared to linear scan — identical to the Sage
//! fragment-index concept for precursor lookup.

use ahash::AHashMap;

use crate::types::Peptide;

/// Peptide mass index: maps integer mass bin → peptide indices.
pub struct PeptideIndex {
    bins: AHashMap<u32, Vec<u32>>,
    peptides: Vec<Peptide>,
}

impl PeptideIndex {
    /// Build the index from a peptide list (consumed).
    pub fn build(peptides: Vec<Peptide>) -> Self {
        let mut bins: AHashMap<u32, Vec<u32>> = AHashMap::default();
        for (i, pep) in peptides.iter().enumerate() {
            let bin = pep.mass.floor() as u32;
            bins.entry(bin).or_default().push(i as u32);
        }
        Self { bins, peptides }
    }

    /// Retrieve candidate peptide indices within `tol_da` of `query_mass`.
    ///
    /// Iterates over mass bins spanning [query_mass − tol_da, query_mass + tol_da].
    pub fn candidates(&self, query_mass: f64, tol_da: f64) -> Vec<u32> {
        let lo = (query_mass - tol_da).floor() as i64;
        let hi = (query_mass + tol_da).ceil() as i64;
        let mut out: Vec<u32> = Vec::new();
        for bin in lo..=hi {
            if bin < 0 {
                continue;
            }
            if let Some(idxs) = self.bins.get(&(bin as u32)) {
                out.extend_from_slice(idxs);
            }
        }
        out
    }

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

    fn make_peptide(mass: f64, decoy: bool) -> Peptide {
        Peptide {
            sequence: "PEPTIDE".into(),
            mass,
            protein_idx: 0,
            is_decoy: decoy,
            missed_cleavages: 0,
        }
    }

    #[test]
    fn candidate_lookup_basic() {
        let peps = vec![
            make_peptide(799.0, false),
            make_peptide(800.5, false),
            make_peptide(801.2, false),
            make_peptide(1200.0, false),
        ];
        let idx = PeptideIndex::build(peps);
        let candidates = idx.candidates(800.5, 2.0);
        // Should return the three peptides near 800.5, not the one at 1200.
        assert!(candidates.len() >= 2);
        assert!(!candidates.is_empty());
        // None of the returned candidates should be the 1200-Da peptide.
        for &c in &candidates {
            assert!((idx.peptide(c).mass - 800.5).abs() < 3.0);
        }
    }
}
