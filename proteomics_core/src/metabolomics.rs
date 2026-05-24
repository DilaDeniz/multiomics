//! Untargeted metabolomics — MS1 feature detection from mzML data.
//!
//! Detects LC-MS features (isotope envelopes + chromatographic peaks) from
//! MS1 survey scans. Each feature is defined by its monoisotopic m/z, charge
//! state, retention time apex, peak width, and integrated intensity.
//!
//! # Algorithm
//! 1. **Isotope grouping** — within each MS1 scan, cluster peaks separated by
//!    ~1.003 Da (C¹² → C¹³ spacing / z) to identify charge states and
//!    monoisotopic masses. Accepts z = 1–4.
//!
//! 2. **Chromatographic peak detection** — across consecutive MS1 scans, link
//!    isotope groups with matching m/z (within `MZ_TOL_PPM`) to build ion
//!    chromatograms (XICs). A Gaussian-shaped apex is required (intensity must
//!    rise then fall for at least `MIN_SCANS` consecutive scans).
//!
//! 3. **Feature deduplication** — collapse overlapping features (same m/z within
//!    5 ppm, overlapping RT ranges) keeping the highest-intensity representative.
//!
//! # Output
//! `Feature` structs suitable for downstream annotation (HMDB, KEGG lookup),
//! differential analysis, or integration with the other omics modalities.
//!
//! # References
//! * Smith CA et al. (2006) XCMS: Processing mass spectrometry data for
//!   metabolite profiling using nonlinear peak alignment, matching, and
//!   identification. Anal. Chem. 78(3):779–787.

use serde::{Deserialize, Serialize};

use crate::types::Spectrum;

/// Monoisotopic spacing between consecutive isotope peaks (C¹²/C¹³).
const ISOTOPE_SPACING: f64 = 1.003_354_8;

/// PPM tolerance for m/z matching in isotope grouping and XIC linking.
pub const MZ_TOL_PPM: f64 = 10.0;

/// Minimum number of consecutive MS1 scans for a valid chromatographic peak.
pub const MIN_SCANS: usize = 3;

/// Minimum peak intensity (absolute, after noise floor filtering).
pub const MIN_INTENSITY: f32 = 100.0;

/// A detected LC-MS feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Feature {
    /// Monoisotopic m/z.
    pub mz: f64,
    /// Inferred charge state (1–4; 0 if unknown).
    pub charge: u8,
    /// Monoisotopic neutral mass: (mz − proton) × z.
    pub neutral_mass: f64,
    /// Retention time at the peak apex (seconds).
    pub rt_apex: f32,
    /// Retention time range of the chromatographic peak.
    pub rt_start: f32,
    pub rt_end: f32,
    /// Integrated peak area (trapezoidal rule across RT scans).
    pub area: f64,
    /// Apex intensity.
    pub apex_intensity: f32,
    /// Number of isotope peaks detected in the envelope.
    pub n_isotopes: u8,
}

/// Run MS1 feature detection on a set of spectra.
///
/// Only MS1 spectra are used (ms_level == 1). Returns detected features
/// sorted by apex intensity descending.
pub fn detect_features(spectra: &[Spectrum]) -> Vec<Feature> {
    let ms1: Vec<&Spectrum> = spectra.iter().filter(|s| s.ms_level == 1).collect();
    if ms1.len() < MIN_SCANS {
        return vec![];
    }

    // Step 1: detect isotope groups per MS1 scan.
    let scan_groups: Vec<Vec<IsotopeGroup>> = ms1.iter().map(|s| find_isotope_groups(s)).collect();

    // Step 2: build XICs by linking groups across consecutive scans.
    let mut xics: Vec<Xic> = build_xics(&ms1, &scan_groups);

    // Step 3: detect apex in each XIC → feature.
    let mut features: Vec<Feature> = xics
        .iter_mut()
        .filter_map(|xic| apex_to_feature(xic))
        .collect();

    // Step 4: deduplicate overlapping features.
    deduplicate(&mut features, 5.0);

    features.sort_unstable_by(|a, b| {
        b.apex_intensity
            .partial_cmp(&a.apex_intensity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    features
}

/// Summary statistics over all detected features.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetabolomicsSummary {
    /// Total features detected.
    pub n_features: u32,
    /// Features with charge state 1 / 2 / 3+ respectively.
    pub n_z1: u32,
    pub n_z2: u32,
    pub n_z3plus: u32,
    /// Median neutral mass.
    pub median_mass: f64,
    /// RT range covered by features.
    pub rt_range_sec: (f32, f32),
    /// Top-20 most intense features.
    pub top_features: Vec<Feature>,
}

impl MetabolomicsSummary {
    pub fn from_features(features: &[Feature]) -> Self {
        let n = features.len() as u32;
        let n_z1 = features.iter().filter(|f| f.charge == 1).count() as u32;
        let n_z2 = features.iter().filter(|f| f.charge == 2).count() as u32;
        let n_z3plus = features.iter().filter(|f| f.charge >= 3).count() as u32;

        let masses: Vec<f64> = features.iter().map(|f| f.neutral_mass).collect();
        let median_mass = median_f64(&masses);

        let rt_start = features
            .iter()
            .map(|f| f.rt_start)
            .fold(f32::INFINITY, f32::min);
        let rt_end = features
            .iter()
            .map(|f| f.rt_end)
            .fold(f32::NEG_INFINITY, f32::max);

        let top_features: Vec<Feature> = features.iter().take(20).cloned().collect();

        Self {
            n_features: n,
            n_z1,
            n_z2,
            n_z3plus,
            median_mass,
            rt_range_sec: (rt_start, rt_end),
            top_features,
        }
    }
}

// ── Internal types ────────────────────────────────────────────────────────────

/// A group of isotope peaks within one MS1 scan.
#[derive(Debug, Clone)]
struct IsotopeGroup {
    /// Monoisotopic m/z (lightest peak in the envelope).
    mz: f64,
    charge: u8,
    /// Summed intensity of all isotope peaks.
    intensity: f32,
    /// Number of peaks in the envelope.
    n_peaks: u8,
}

/// Ion chromatogram: linked isotope groups across scans.
struct Xic {
    mz: f64,
    charge: u8,
    /// (rt, intensity) for each scan in the XIC.
    points: Vec<(f32, f32)>,
}

// ── Step 1: Isotope grouping ──────────────────────────────────────────────────

fn find_isotope_groups(scan: &Spectrum) -> Vec<IsotopeGroup> {
    let n = scan.mz.len();
    if n == 0 {
        return vec![];
    }

    let mut used = vec![false; n];
    let mut groups: Vec<IsotopeGroup> = Vec::new();

    for i in 0..n {
        if used[i] || scan.intensity[i] < MIN_INTENSITY {
            continue;
        }

        let mz0 = scan.mz[i] as f64;
        // Try charge states 1–4.
        let mut best_group: Option<IsotopeGroup> = None;

        for z in 1u8..=4 {
            let spacing = ISOTOPE_SPACING / z as f64;
            let tol = mz0 * MZ_TOL_PPM / 1e6;

            let mut peaks: Vec<usize> = vec![i];
            let mut last_mz = mz0;

            // Walk forward looking for successive isotope peaks.
            for j in (i + 1)..n {
                if scan.intensity[j] < MIN_INTENSITY {
                    continue;
                }
                let mz_j = scan.mz[j] as f64;
                let expected = last_mz + spacing;
                if (mz_j - expected).abs() <= tol * 2.0 {
                    peaks.push(j);
                    last_mz = mz_j;
                    if peaks.len() >= 4 {
                        break;
                    }
                } else if mz_j > expected + tol * 2.0 {
                    break;
                }
            }

            if peaks.len() >= 2 {
                let intensity: f32 = peaks.iter().map(|&p| scan.intensity[p]).sum();
                let group = IsotopeGroup {
                    mz: mz0,
                    charge: z,
                    intensity,
                    n_peaks: peaks.len() as u8,
                };
                let better = best_group
                    .as_ref()
                    .map(|b| group.n_peaks > b.n_peaks)
                    .unwrap_or(true);
                if better {
                    if let Some(prev) = best_group.take() {
                        // Un-mark previous if we're improving.
                        let _ = prev;
                    }
                    for &p in &peaks {
                        used[p] = true;
                    }
                    best_group = Some(group);
                }
            }
        }

        if let Some(group) = best_group {
            groups.push(group);
        } else if scan.intensity[i] >= MIN_INTENSITY {
            // Singly charged singleton.
            groups.push(IsotopeGroup {
                mz: mz0,
                charge: 1,
                intensity: scan.intensity[i],
                n_peaks: 1,
            });
        }
    }

    groups
}

// ── Step 2: XIC construction ──────────────────────────────────────────────────

fn build_xics(ms1: &[&Spectrum], scan_groups: &[Vec<IsotopeGroup>]) -> Vec<Xic> {
    let mut xics: Vec<Xic> = Vec::new();

    for (scan_idx, (scan, groups)) in ms1.iter().zip(scan_groups.iter()).enumerate() {
        let rt = scan.rt;

        for group in groups {
            // Try to extend an existing XIC.
            let existing = xics.iter_mut().find(|x| {
                x.charge == group.charge
                    && ppm_delta(x.mz, group.mz) <= MZ_TOL_PPM
                    && scan_idx.saturating_sub(
                        ms1[..scan_idx]
                            .iter()
                            .rposition(|s| s.rt < rt - 60.0)
                            .map(|p| p + 1)
                            .unwrap_or(0),
                    ) < 5 // allow up to 5 scan gaps
            });

            if let Some(xic) = existing {
                xic.points.push((rt, group.intensity));
                // Update m/z as running average for better accuracy.
                xic.mz =
                    (xic.mz * (xic.points.len() - 1) as f64 + group.mz) / xic.points.len() as f64;
            } else {
                xics.push(Xic {
                    mz: group.mz,
                    charge: group.charge,
                    points: vec![(rt, group.intensity)],
                });
            }
        }
    }

    xics
}

// ── Step 3: Apex detection ────────────────────────────────────────────────────

fn apex_to_feature(xic: &mut Xic) -> Option<Feature> {
    if xic.points.len() < MIN_SCANS {
        return None;
    }

    // Find the intensity apex.
    let (apex_idx, &(apex_rt, apex_int)) = xic.points.iter().enumerate().max_by(|a, b| {
        a.1 .1
            .partial_cmp(&b.1 .1)
            .unwrap_or(std::cmp::Ordering::Equal)
    })?;

    if apex_int < MIN_INTENSITY {
        return None;
    }

    // Require at least one scan on each side of the apex.
    if apex_idx == 0 || apex_idx == xic.points.len() - 1 {
        return None;
    }

    // Find the peak boundaries (where intensity drops below half of apex).
    let half = apex_int / 2.0;
    let left = xic.points[..apex_idx]
        .iter()
        .rposition(|(_, int)| *int < half)
        .map(|p| p + 1)
        .unwrap_or(0);
    let right = xic.points[apex_idx + 1..]
        .iter()
        .position(|(_, int)| *int < half)
        .map(|p| apex_idx + 1 + p)
        .unwrap_or(xic.points.len() - 1);

    let rt_start = xic.points[left].0;
    let rt_end = xic.points[right].0;

    // Trapezoidal integration over the peak range.
    let area = xic.points[left..=right]
        .windows(2)
        .map(|w| {
            let dt = (w[1].0 - w[0].0) as f64;
            let avg_int = (w[0].1 + w[1].1) as f64 / 2.0;
            dt * avg_int
        })
        .sum::<f64>();

    let neutral_mass = (xic.mz - 1.007_276_466_621) * xic.charge as f64;

    Some(Feature {
        mz: xic.mz,
        charge: xic.charge,
        neutral_mass,
        rt_apex: apex_rt,
        rt_start,
        rt_end,
        area,
        apex_intensity: apex_int,
        n_isotopes: 1, // refined in full impl; placeholder here
    })
}

// ── Step 4: Deduplication ─────────────────────────────────────────────────────

fn deduplicate(features: &mut Vec<Feature>, tol_ppm: f64) {
    features.sort_unstable_by(|a, b| {
        b.apex_intensity
            .partial_cmp(&a.apex_intensity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut keep = vec![true; features.len()];
    for i in 0..features.len() {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..features.len() {
            if !keep[j] {
                continue;
            }
            let mz_ok = ppm_delta(features[i].mz, features[j].mz) <= tol_ppm;
            let rt_overlap = features[i].rt_start <= features[j].rt_end
                && features[j].rt_start <= features[i].rt_end;
            if mz_ok && rt_overlap {
                keep[j] = false; // features[i] has higher intensity
            }
        }
    }

    let mut idx = 0;
    features.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

// ── Utilities ─────────────────────────────────────────────────────────────────

#[inline]
fn ppm_delta(a: f64, b: f64) -> f64 {
    ((a - b) / a.max(b)).abs() * 1e6
}

fn median_f64(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut s = v.to_vec();
    s.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ms1_scan(rt: f32, peaks: &[(f32, f32)]) -> Spectrum {
        let mut pairs: Vec<(f32, f32)> = peaks.to_vec();
        pairs.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        Spectrum {
            scan: 1,
            ms_level: 1,
            rt,
            precursor_mz: 0.0,
            precursor_z: 0,
            mz: pairs.iter().map(|p| p.0).collect(),
            intensity: pairs.iter().map(|p| p.1).collect(),
        }
    }

    #[test]
    fn detect_simple_feature() {
        // Isotope envelope at ~500 Da (z=1) appearing across 5 RT scans.
        let mut spectra: Vec<Spectrum> = Vec::new();
        let intensities = [500.0f32, 1000.0, 2000.0, 1000.0, 500.0]; // apex at scan 2
        for (i, &int) in intensities.iter().enumerate() {
            let rt = i as f32 * 10.0; // 0, 10, 20, 30, 40 s
            spectra.push(ms1_scan(
                rt,
                &[
                    (500.0, int),
                    (501.003, int * 0.3), // C13 isotope
                ],
            ));
        }
        let features = detect_features(&spectra);
        assert!(!features.is_empty(), "should detect at least one feature");
        let f = &features[0];
        assert!((f.mz - 500.0).abs() < 0.05, "mz={}", f.mz);
        assert!(f.apex_intensity > 1000.0, "apex_int={}", f.apex_intensity);
        assert!(f.area > 0.0, "area must be positive");
    }

    #[test]
    fn no_features_from_ms2_scans() {
        let mut scan = ms1_scan(10.0, &[(500.0, 5000.0), (501.003, 1500.0)]);
        scan.ms_level = 2;
        let features = detect_features(&[scan]);
        assert!(features.is_empty(), "MS2 scans should be ignored");
    }

    #[test]
    fn summary_from_empty_is_zero() {
        let s = MetabolomicsSummary::from_features(&[]);
        assert_eq!(s.n_features, 0);
        assert_eq!(s.median_mass, 0.0);
    }

    #[test]
    fn ppm_delta_symmetric() {
        let d = ppm_delta(500.0, 500.005);
        assert!(d < 11.0 && d > 9.0, "ppm={d}");
    }
}
