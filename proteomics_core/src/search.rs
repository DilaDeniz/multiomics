//! Core database search: score all MS2 spectra against the peptide index.
//!
//! Each spectrum is processed independently. For each MS2:
//!   1. Compute precursor neutral mass from precursor m/z and charge.
//!   2. Two-stage candidate filter: precursor mass bin + fragment vote index.
//!   3. Score each passing candidate with the hyperscore.
//!   4. Keep the best and second-best scores → winner PSM + delta score.
//!
//! The fragment vote pre-filter (Sage-inspired) eliminates 5–20× more
//! candidates before the expensive hyperscore step, dramatically cutting
//! total search time on large databases.
//!
//! All spectra are processed in parallel via rayon.

use ahash::{AHashMap, AHashSet};
use rayon::prelude::*;

use crate::fdr::assign_qvalues;
use crate::index::{PeptideIndex, MIN_FRAGMENT_VOTES};
use crate::score::hyperscore;
use crate::types::{Psm, Spectrum};

/// Precursor mass tolerance in ppm (10 ppm matches most high-res instruments).
pub const PRECURSOR_TOL_PPM: f64 = 10.0;

/// Run the full database search on a set of MS2 spectra.
///
/// Returns all PSMs after target-decoy competition and q-value assignment
/// (both targets and decoys are returned; callers filter on `q_value ≤ 0.01`
/// and `!is_decoy`).
pub fn search_spectra(
    spectra: &[Spectrum],
    index: &PeptideIndex,
    protein_names: &[String],
    frag_tol_ppm: f64,
    precursor_tol_ppm: f64,
) -> Vec<Psm> {
    let ms2: Vec<&Spectrum> = spectra.iter().filter(|s| s.ms_level == 2).collect();

    let mut psms: Vec<Psm> = ms2
        .par_iter()
        .filter_map(|spec| {
            if spec.precursor_z == 0 || spec.precursor_mz == 0.0 {
                return None;
            }
            let precursor_mass = spec.precursor_mass();

            // Two-stage filter: precursor mass bin + fragment ion votes.
            // candidates_voted returns (pep_idx, vote_count) sorted by votes desc.
            let voted = index.candidates_voted(
                precursor_mass,
                precursor_tol_ppm,
                &spec.mz,
                MIN_FRAGMENT_VOTES,
            );
            if voted.is_empty() {
                return None;
            }

            let mut best_score = 0.0f64;
            let mut second_score = 0.0f64;
            let mut best_pep_idx: Option<u32> = None;
            let mut best_nb = 0u32;
            let mut best_ny = 0u32;

            for &(cand_idx, _votes) in &voted {
                let pep = index.peptide(cand_idx);
                // Precise PPM filter on the candidate mass.
                let mass_err_ppm = ((pep.mass - precursor_mass) / precursor_mass * 1e6).abs();
                if mass_err_ppm > precursor_tol_ppm {
                    continue;
                }

                let (score, nb, ny) = hyperscore(spec, pep, frag_tol_ppm);
                if score > best_score {
                    second_score = best_score;
                    best_score = score;
                    best_pep_idx = Some(cand_idx);
                    best_nb = nb;
                    best_ny = ny;
                } else if score > second_score {
                    second_score = score;
                }
            }

            let pep_idx = best_pep_idx?;
            if best_score <= 0.0 {
                return None;
            }

            let pep = index.peptide(pep_idx);
            let prot_name = protein_names
                .get(pep.protein_idx as usize)
                .cloned()
                .unwrap_or_else(|| "unknown".into());

            let mass_error_ppm = ((pep.mass - precursor_mass) / precursor_mass * 1e6).abs();

            Some(Psm {
                scan: spec.scan,
                rt: spec.rt,
                peptide: pep.sequence.clone(),
                protein: prot_name,
                charge: spec.precursor_z,
                hyperscore: best_score,
                delta_score: best_score - second_score,
                mass_error_ppm,
                n_matched_b: best_nb,
                n_matched_y: best_ny,
                q_value: 1.0,
                is_decoy: pep.is_decoy,
            })
        })
        .collect();

    assign_qvalues(&mut psms);
    psms
}

/// Protein inference: group PSMs by protein and compute protein-level q-values.
pub fn infer_proteins(psms: &[Psm]) -> Vec<crate::types::ProteinGroup> {
    let mut groups: AHashMap<&str, (u32, AHashSet<&str>, f64, bool)> = AHashMap::default();

    for psm in psms {
        let e = groups
            .entry(psm.protein.as_str())
            .or_insert_with(|| (0, AHashSet::default(), 0.0, psm.is_decoy));
        e.0 += 1;
        e.1.insert(psm.peptide.as_str());
        if psm.hyperscore > e.2 {
            e.2 = psm.hyperscore;
        }
    }

    let mut proteins: Vec<crate::types::ProteinGroup> = groups
        .iter()
        .map(
            |(&name, &(n_psms, ref peps, top_score, decoy))| crate::types::ProteinGroup {
                protein: name.to_string(),
                n_psms,
                n_unique_peptides: peps.len() as u32,
                top_score,
                q_value: 1.0,
                is_decoy: decoy,
            },
        )
        .collect();

    // Sort by top score descending, assign protein q-values.
    proteins.sort_unstable_by(|a, b| {
        b.top_score
            .partial_cmp(&a.top_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut n_target = 0u64;
    let mut n_decoy = 0u64;
    let n = proteins.len();
    let mut fdr = Vec::with_capacity(n);
    for p in &proteins {
        if p.is_decoy {
            n_decoy += 1;
        } else {
            n_target += 1;
        }
        fdr.push(((n_decoy as f64 + 1.0) / n_target.max(1) as f64).min(1.0));
    }
    // Backward rolling minimum.
    let mut min_fdr = 1.0f64;
    let mut q_vals: Vec<f64> = fdr
        .iter()
        .rev()
        .map(|&f| {
            min_fdr = min_fdr.min(f);
            min_fdr
        })
        .collect();
    q_vals.reverse();
    for (p, q) in proteins.iter_mut().zip(q_vals.iter()) {
        p.q_value = *q;
    }

    proteins
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::{digest, parse_fasta};
    use crate::score::FRAG_TOL_PPM;

    #[test]
    fn search_empty_spectra() {
        let fasta = b">PROT\nPEPTIDEKAAAAR\n";
        let prots = parse_fasta(fasta);
        let peps = digest(&prots, 1, 5, 50);
        let names: Vec<String> = prots.iter().map(|p| p.header.clone()).collect();
        let index = PeptideIndex::build(peps);
        let results = search_spectra(&[], &index, &names, FRAG_TOL_PPM, PRECURSOR_TOL_PPM);
        assert!(results.is_empty());
    }
}
