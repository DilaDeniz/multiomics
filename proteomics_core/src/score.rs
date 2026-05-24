//! Fragment ion generation and PSM scoring.
//!
//! Implements the hyperscore from X!Tandem (Craig & Beavis 2004):
//!
//!   hyperscore = ln(dot_b + 1) + ln(dot_y + 1)
//!                + ln(n_b!) + ln(n_y!)
//!
//! where dot_b / dot_y are summed intensities of matched b / y ions and
//! n_b / n_y are counts of matched b / y ions.
//!
//! Fragment masses are computed for singly and doubly charged b/y ions.
//! Matching uses binary search on the sorted observed m/z array with a
//! configurable PPM tolerance (default: 20 ppm).

use crate::types::{aa_mass, Peptide, Spectrum, PROTON, WATER};

/// Default fragment matching tolerance in ppm.
pub const FRAG_TOL_PPM: f64 = 20.0;

/// Score a single peptide against an MS2 spectrum.
///
/// Returns `(hyperscore, n_matched_b, n_matched_y)`.
pub fn hyperscore(spectrum: &Spectrum, peptide: &Peptide, tol_ppm: f64) -> (f64, u32, u32) {
    let seq = peptide.sequence.as_bytes();
    let n = seq.len();
    if n < 2 {
        return (0.0, 0, 0);
    }

    // Build per-position modification deltas.
    // mod_delta[i] = total mass shift of all modifications AT position i.
    let mut mod_delta = vec![0.0f64; n];
    for m in &peptide.modifications {
        if m.position < n {
            mod_delta[m.position] += m.mass_delta;
        }
    }

    // Prefix cumulative modification mass (for b ions: mods at 0..=i contribute to b_i).
    let mut mod_prefix_cum = vec![0.0f64; n];
    mod_prefix_cum[0] = mod_delta[0];
    for i in 1..n {
        mod_prefix_cum[i] = mod_prefix_cum[i - 1] + mod_delta[i];
    }

    // Suffix cumulative modification mass (for y ions: mods at i..n contribute to y starting at i).
    let mut mod_suffix_cum = vec![0.0f64; n];
    mod_suffix_cum[n - 1] = mod_delta[n - 1];
    for i in (0..n - 1).rev() {
        mod_suffix_cum[i] = mod_suffix_cum[i + 1] + mod_delta[i];
    }

    // Pre-compute prefix residue masses for b ions and suffix masses for y ions.
    let mut prefix = vec![0.0f64; n];
    let mut suffix = vec![0.0f64; n];

    prefix[0] = aa_mass(seq[0]);
    for i in 1..n {
        prefix[i] = prefix[i - 1] + aa_mass(seq[i]);
    }
    suffix[n - 1] = aa_mass(seq[n - 1]);
    for i in (0..n - 1).rev() {
        suffix[i] = suffix[i + 1] + aa_mass(seq[i]);
    }

    let mut dot_b = 0.0f64;
    let mut dot_y = 0.0f64;
    let mut n_b = 0u32;
    let mut n_y = 0u32;

    // b ions: b_i covers residues 0..=i, including mods at those positions.
    for (i, &pm) in prefix[..n - 1].iter().enumerate() {
        let b_mass = pm + mod_prefix_cum[i];
        let b1 = b_mass + PROTON;
        let b2 = (b_mass + 2.0 * PROTON) / 2.0;
        if let Some(int) = find_peak(spectrum, b1 as f32, tol_ppm) {
            dot_b += int as f64;
            n_b += 1;
        } else if let Some(int) = find_peak(spectrum, b2 as f32, tol_ppm) {
            dot_b += int as f64;
            n_b += 1;
        }
    }

    // y ions: y_i covers residues i..n, including mods at those positions.
    for (i, &sm) in suffix[1..].iter().enumerate() {
        let y_mass = sm + WATER + mod_suffix_cum[i + 1];
        let y1 = y_mass + PROTON;
        let y2 = (y_mass + 2.0 * PROTON) / 2.0;
        if let Some(int) = find_peak(spectrum, y1 as f32, tol_ppm) {
            dot_y += int as f64;
            n_y += 1;
        } else if let Some(int) = find_peak(spectrum, y2 as f32, tol_ppm) {
            dot_y += int as f64;
            n_y += 1;
        }
    }

    let score = (dot_b + 1.0).ln() + (dot_y + 1.0).ln() + log_factorial(n_b) + log_factorial(n_y);

    (score, n_b, n_y)
}

/// Binary-search the sorted observed m/z array for a theoretical fragment.
/// Returns the normalized intensity if found within `tol_ppm`.
#[inline]
fn find_peak(spectrum: &Spectrum, mz_target: f32, tol_ppm: f64) -> Option<f32> {
    let mz = &spectrum.mz;
    if mz.is_empty() {
        return None;
    }
    let tol = (mz_target as f64 * tol_ppm / 1e6) as f32;
    let lo = mz_target - tol;
    let hi = mz_target + tol;

    // Binary search for the leftmost peak ≥ lo.
    let start = mz.partition_point(|&x| x < lo);
    // Scan right to find the closest peak within [lo, hi].
    let mut best: Option<(f32, f32)> = None; // (delta, intensity)
    for (mz_val, int_val) in mz[start..].iter().zip(spectrum.intensity[start..].iter()) {
        if *mz_val > hi {
            break;
        }
        let delta = (*mz_val - mz_target).abs();
        if best.is_none_or(|(d, _)| delta < d) {
            best = Some((delta, *int_val));
        }
    }
    best.map(|(_, int)| int)
}

/// ln(n!) via Stirling for n > 20, exact table for n ≤ 20.
fn log_factorial(n: u32) -> f64 {
    const TABLE: [f64; 21] = [
        0.0,
        0.0,
        std::f64::consts::LN_2,
        1.791_759_47,
        3.178_053_83,
        4.787_491_74,
        6.579_251_21,
        8.525_161_36,
        10.604_602_9,
        12.801_827_5,
        15.104_412_6,
        17.502_307_8,
        19.987_214_5,
        22.552_163_9,
        25.191_221_2,
        27.899_271_4,
        30.671_860_1,
        33.505_073_4,
        36.395_445_2,
        39.339_884_2,
        42.335_616_5,
    ];
    if n <= 20 {
        TABLE[n as usize]
    } else {
        let x = n as f64;
        x * x.ln() - x + 0.5 * (2.0 * std::f64::consts::PI * x).ln()
    }
}

/// Compute theoretical b-ion masses for a peptide (z=1, not including water/proton overhead).
pub fn b_ions(seq: &[u8]) -> Vec<f64> {
    let mut masses = Vec::with_capacity(seq.len().saturating_sub(1));
    let mut sum = 0.0;
    for &aa in &seq[..seq.len().saturating_sub(1)] {
        sum += aa_mass(aa);
        masses.push(sum + PROTON);
    }
    masses
}

/// Compute theoretical y-ion masses for a peptide (z=1).
pub fn y_ions(seq: &[u8]) -> Vec<f64> {
    let mut masses = Vec::with_capacity(seq.len().saturating_sub(1));
    let mut sum = WATER;
    for &aa in seq[1..].iter().rev() {
        sum += aa_mass(aa);
        masses.push(sum + PROTON);
    }
    masses
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_spectrum(peaks: &[(f32, f32)]) -> Spectrum {
        let mut mz: Vec<f32> = peaks.iter().map(|&(m, _)| m).collect();
        let mut int: Vec<f32> = peaks.iter().map(|&(_, i)| i).collect();
        // Sort by m/z (required for binary search).
        let mut pairs: Vec<(f32, f32)> = mz.iter().copied().zip(int.iter().copied()).collect();
        pairs.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        mz = pairs.iter().map(|p| p.0).collect();
        int = pairs.iter().map(|p| p.1).collect();
        Spectrum {
            scan: 1,
            ms_level: 2,
            rt: 30.0,
            precursor_mz: 500.0,
            precursor_z: 2,
            mz,
            intensity: int,
        }
    }

    fn make_peptide(seq: &str) -> Peptide {
        use crate::fasta::peptide_mass;
        Peptide {
            mass: peptide_mass(seq.as_bytes()),
            sequence: seq.into(),
            protein_idx: 0,
            is_decoy: false,
            missed_cleavages: 0,
            modifications: vec![],
        }
    }

    #[test]
    fn b_y_ion_count() {
        // PEPTIDE has 7 residues → 6 b ions and 6 y ions.
        let b = b_ions(b"PEPTIDE");
        let y = y_ions(b"PEPTIDE");
        assert_eq!(b.len(), 6);
        assert_eq!(y.len(), 6);
    }

    #[test]
    fn hyperscore_perfect_match() {
        // Build a spectrum whose peaks ARE the b/y ions of PEPTIDE.
        let seq = b"PEPTIDE";
        let mut peaks: Vec<(f32, f32)> = Vec::new();
        for m in b_ions(seq) {
            peaks.push((m as f32, 1.0));
        }
        for m in y_ions(seq) {
            peaks.push((m as f32, 1.0));
        }
        let spectrum = make_spectrum(&peaks);
        let peptide = make_peptide("PEPTIDE");
        let (score, nb, ny) = hyperscore(&spectrum, &peptide, 10.0);
        assert!(score > 0.0, "score should be positive: {score}");
        assert!(nb > 0, "should match b ions");
        assert!(ny > 0, "should match y ions");
    }

    #[test]
    fn hyperscore_no_match() {
        // Spectrum peaks far from PEPTIDE fragment ions.
        let peaks: Vec<(f32, f32)> = (0..5).map(|i| (50.0 + i as f32, 1.0)).collect();
        let spectrum = make_spectrum(&peaks);
        let peptide = make_peptide("PEPTIDE");
        let (score, nb, ny) = hyperscore(&spectrum, &peptide, 10.0);
        // With no matches, dot_b=0 and dot_y=0 → ln(1)+ln(1)+0+0 = 0.
        assert_eq!(score, 0.0);
        assert_eq!(nb, 0);
        assert_eq!(ny, 0);
    }
}
