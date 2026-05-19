//! Label-free quantification via MS1 extracted-ion chromatograms (XIC).
//!
//! For each identified peptide, integrates MS1 peak area over a ±5 ppm m/z
//! window and a ±30 s retention-time window around the identified PSM's RT.
//! Peak area is computed by trapezoidal integration.

use crate::types::{Psm, Spectrum, PROTON};

/// Peptide quantification result.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PeptideQuant {
    pub peptide: String,
    pub protein: String,
    pub charge: u8,
    pub rt_apex: f32,
    pub xic_area: f64,
    pub n_points: usize,
}

/// Extract XIC peak areas for a list of PSMs from MS1 spectra.
///
/// `ms1_spectra` must already be filtered to MS level 1 and sorted by RT.
/// `rt_window_sec`: half-width of RT window (default 30 s).
/// `mz_tol_ppm`: m/z tolerance for XIC extraction (default 5 ppm).
pub fn quantify_psms(
    psms: &[Psm],
    ms1_spectra: &[Spectrum],
    rt_window_sec: f32,
    mz_tol_ppm: f64,
) -> Vec<PeptideQuant> {
    psms.iter()
        .map(|psm| {
            let mz_center = (psm.charge as f64 + PROTON * psm.charge as f64) / psm.charge as f64;
            let area = xic_area(
                ms1_spectra,
                mz_center as f32,
                psm.rt,
                rt_window_sec,
                mz_tol_ppm,
            );
            PeptideQuant {
                peptide: psm.peptide.clone(),
                protein: psm.protein.clone(),
                charge: psm.charge,
                rt_apex: psm.rt,
                xic_area: area,
                n_points: count_points(ms1_spectra, psm.rt, rt_window_sec),
            }
        })
        .collect()
}

fn xic_area(ms1: &[Spectrum], mz_center: f32, rt_center: f32, rt_window: f32, tol_ppm: f64) -> f64 {
    let tol = (mz_center as f64 * tol_ppm / 1e6) as f32;
    let rt_lo = rt_center - rt_window;
    let rt_hi = rt_center + rt_window;

    let mut points: Vec<(f32, f64)> = Vec::new(); // (rt, intensity)

    for spec in ms1 {
        if spec.rt < rt_lo || spec.rt > rt_hi {
            continue;
        }
        let int = extract_intensity(&spec.mz, &spec.intensity, mz_center, tol);
        points.push((spec.rt, int));
    }

    if points.len() < 2 {
        return points.first().map(|p| p.1).unwrap_or(0.0);
    }
    points.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Trapezoidal integration.
    let mut area = 0.0f64;
    for w in points.windows(2) {
        let dt = (w[1].0 - w[0].0) as f64;
        area += 0.5 * (w[0].1 + w[1].1) * dt;
    }
    area
}

fn extract_intensity(mz: &[f32], intensity: &[f32], center: f32, tol: f32) -> f64 {
    let lo = center - tol;
    let hi = center + tol;
    let start = mz.partition_point(|&x| x < lo);
    let mut sum = 0.0f64;
    for i in start..mz.len() {
        if mz[i] > hi {
            break;
        }
        sum += intensity[i] as f64;
    }
    sum
}

fn count_points(ms1: &[Spectrum], rt_center: f32, rt_window: f32) -> usize {
    ms1.iter()
        .filter(|s| (s.rt - rt_center).abs() <= rt_window)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ms1_spec(rt: f32, mz: f32, int: f32) -> Spectrum {
        Spectrum {
            scan: 0,
            ms_level: 1,
            rt,
            precursor_mz: 0.0,
            precursor_z: 0,
            mz: vec![mz],
            intensity: vec![int],
        }
    }

    #[test]
    fn xic_area_triangle() {
        let ms1 = vec![
            ms1_spec(0.0, 500.0, 0.0),
            ms1_spec(10.0, 500.0, 1000.0),
            ms1_spec(20.0, 500.0, 0.0),
        ];
        // Triangle: base 20 s, height 1000 → area = 10_000.
        let area = xic_area(&ms1, 500.0, 10.0, 30.0, 10.0);
        assert!((area - 10_000.0).abs() < 1.0, "area: {area}");
    }
}
