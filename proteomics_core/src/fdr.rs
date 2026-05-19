//! Target-decoy FDR estimation and q-value computation.
//!
//! Uses the standard target-decoy competition model:
//! for each spectrum the best-scoring target PSM competes with the best-scoring
//! decoy PSM; only the winner is kept. PSMs are then sorted by score descending
//! and q-values computed as a rolling minimum of the estimated FDR.
//!
//! Reference: Elias & Gygi (2007) Nature Methods 4:207-214.

use crate::types::Psm;

/// Assign q-values to a list of PSMs in-place.
///
/// The input slice must contain the competition-winner PSM per spectrum
/// (both targets and decoys). PSMs are sorted by `hyperscore` descending.
/// After this call, `psm.q_value` is the minimum-FDR estimate at each
/// score threshold.
pub fn assign_qvalues(psms: &mut [Psm]) {
    if psms.is_empty() {
        return;
    }

    // Sort by score descending.
    psms.sort_unstable_by(|a, b| {
        b.hyperscore
            .partial_cmp(&a.hyperscore)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let n = psms.len();
    let mut n_target = 0u64;
    let mut n_decoy = 0u64;

    let mut fdr: Vec<f64> = Vec::with_capacity(n);
    for psm in psms.iter() {
        if psm.is_decoy {
            n_decoy += 1;
        } else {
            n_target += 1;
        }
        // Laplace-corrected: (decoy + 1) / target
        let est = (n_decoy as f64 + 1.0) / n_target.max(1) as f64;
        fdr.push(est.min(1.0));
    }

    // Backward pass: enforce monotone non-increasing q-value (rolling minimum).
    let mut q = Vec::with_capacity(n);
    let mut running_min = 1.0f64;
    for &f in fdr.iter().rev() {
        running_min = running_min.min(f);
        q.push(running_min);
    }
    q.reverse();

    for (psm, qv) in psms.iter_mut().zip(q.iter()) {
        psm.q_value = *qv;
    }
}

/// Filter PSMs to those passing `fdr_threshold` (e.g. 0.01 for 1 %).
pub fn filter_psms(psms: &[Psm], fdr_threshold: f64) -> Vec<Psm> {
    psms.iter()
        .filter(|p| !p.is_decoy && p.q_value <= fdr_threshold)
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_psm(score: f64, decoy: bool) -> Psm {
        Psm {
            scan: 0,
            rt: 0.0,
            peptide: String::new(),
            protein: String::new(),
            charge: 2,
            hyperscore: score,
            delta_score: 0.0,
            mass_error_ppm: 0.0,
            n_matched_b: 0,
            n_matched_y: 0,
            q_value: 1.0,
            is_decoy: decoy,
        }
    }

    #[test]
    fn qvalues_monotone() {
        let mut psms = vec![
            make_psm(10.0, false),
            make_psm(9.0, false),
            make_psm(8.5, false),
            make_psm(8.0, true),
            make_psm(7.0, false),
            make_psm(6.0, true),
        ];
        assign_qvalues(&mut psms);
        // q-values must be non-decreasing along the sorted PSM list.
        let qvs: Vec<f64> = psms.iter().map(|p| p.q_value).collect();
        for i in 1..qvs.len() {
            assert!(qvs[i] >= qvs[i - 1], "q-values not monotone: {qvs:?}");
        }
    }

    #[test]
    fn filter_removes_decoys_and_highfdr() {
        let mut psms = vec![
            make_psm(10.0, false),
            make_psm(5.0, true),
            make_psm(3.0, false),
        ];
        assign_qvalues(&mut psms);
        let filtered = filter_psms(&psms, 0.01);
        // Only high-confidence targets should survive.
        assert!(filtered.iter().all(|p| !p.is_decoy));
    }
}
