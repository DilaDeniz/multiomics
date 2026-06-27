use serde::{Deserialize, Serialize};
use crate::types::VariantRecord;

/// 6-channel SBS mutation spectrum (no reference required).
///
/// Channels in pyrimidine context: C>A, C>G, C>T, T>A, T>C, T>G.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SbsSpectrum6 {
    pub c_to_a: u64,
    pub c_to_g: u64,
    pub c_to_t: u64,
    pub t_to_a: u64,
    pub t_to_c: u64,
    pub t_to_g: u64,
    pub total_snvs: u64,
    /// Fraction of each channel (array of 6, same order as fields above).
    pub fractions: [f64; 6],
}

/// Detected mutational signature and its estimated contribution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SignatureMatch {
    /// COSMIC SBS signature name (e.g., "SBS1", "SBS2/13").
    pub signature: String,
    /// Estimated fractional contribution [0, 1].
    pub weight: f64,
    /// Associated etiology.
    pub etiology: String,
}

/// Result of COSMIC mutational signature analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MutationalSignatureResult {
    /// 6-channel SBS spectrum.
    pub spectrum_6ch: SbsSpectrum6,
    /// Putative dominant signatures detected from 6-channel pattern.
    pub dominant_signatures: Vec<SignatureMatch>,
    /// Overall mutagenic process summary.
    pub summary: String,
    /// True if APOBEC (SBS2/SBS13) pattern is enriched.
    pub apobec_enriched: bool,
    /// True if tobacco smoking signature (SBS4) detected.
    pub tobacco_signature: bool,
    /// True if MMR deficiency (SBS6) or MSI-associated pattern.
    pub mismatch_repair_deficiency: bool,
    /// True if UV exposure (SBS7) pattern detected (C>T > 65%).
    pub uv_signature: bool,
    /// Note when total_snvs < 50 (low-confidence spectrum).
    pub note: Option<String>,
}

/// Normalize a single-base allele to pyrimidine context.
/// A, G (purine) → their pyrimidine complements T, C.
fn to_pyrimidine(base: u8) -> Option<u8> {
    match base {
        b'C' | b'T' => Some(base),
        b'A' => Some(b'T'),
        b'G' => Some(b'C'),
        _ => None,
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b => b,
    }
}

/// Compute the 6-channel SBS spectrum from SNV variant records.
pub fn compute_sbs6_spectrum(variants: &[VariantRecord]) -> SbsSpectrum6 {
    let mut c_to_a = 0u64;
    let mut c_to_g = 0u64;
    let mut c_to_t = 0u64;
    let mut t_to_a = 0u64;
    let mut t_to_c = 0u64;
    let mut t_to_g = 0u64;

    for v in variants {
        if v.ref_allele.len() != 1 || v.alt_allele.len() != 1 {
            continue; // skip indels and MNPs
        }
        let r = v.ref_allele.as_bytes()[0].to_ascii_uppercase();
        let a = v.alt_allele.as_bytes()[0].to_ascii_uppercase();
        if r == a {
            continue;
        }
        // Normalize: convert to pyrimidine context
        let (ref_py, alt_py) = match to_pyrimidine(r) {
            Some(rp) if rp == r => {
                // Already pyrimidine
                (r, a)
            }
            Some(_) => {
                // Was purine: take complement of both
                (complement(r), complement(a))
            }
            None => continue,
        };

        match (ref_py, alt_py) {
            (b'C', b'A') => c_to_a += 1,
            (b'C', b'G') => c_to_g += 1,
            (b'C', b'T') => c_to_t += 1,
            (b'T', b'A') => t_to_a += 1,
            (b'T', b'C') => t_to_c += 1,
            (b'T', b'G') => t_to_g += 1,
            _ => {}
        }
    }

    let total_snvs = c_to_a + c_to_g + c_to_t + t_to_a + t_to_c + t_to_g;
    let fractions = if total_snvs > 0 {
        let n = total_snvs as f64;
        [
            c_to_a as f64 / n,
            c_to_g as f64 / n,
            c_to_t as f64 / n,
            t_to_a as f64 / n,
            t_to_c as f64 / n,
            t_to_g as f64 / n,
        ]
    } else {
        [0.0; 6]
    };

    SbsSpectrum6 { c_to_a, c_to_g, c_to_t, t_to_a, t_to_c, t_to_g, total_snvs, fractions }
}

/// Detect dominant COSMIC mutational signatures from 6-channel SBS spectrum.
///
/// Pattern-based detection from the 6-channel distribution.
/// Reference: Alexandrov et al. 2020 (Nature), COSMIC v3.3.
pub fn detect_signatures_from_6ch(spec: &SbsSpectrum6) -> MutationalSignatureResult {
    let [fc_a, fc_g, fc_t, ft_a, ft_c, ft_g] = spec.fractions;
    let mut signatures: Vec<SignatureMatch> = Vec::new();
    let mut summary_parts: Vec<&str> = Vec::new();

    let apobec_enriched = fc_t > 0.35 && fc_g > 0.08;
    let tobacco_signature = fc_a > 0.28;
    let uv_signature = fc_t > 0.65;
    let mismatch_repair_deficiency = {
        // MMR-deficient tumors show relatively flat spectrum with elevated C>T and T>C
        let max_frac = [fc_a, fc_g, fc_t, ft_a, ft_c, ft_g]
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        max_frac < 0.45 && fc_t > 0.20 && ft_c > 0.15
    };

    // UV (SBS7a/b): extreme C>T enrichment (skin cancer)
    if uv_signature {
        signatures.push(SignatureMatch {
            signature: "SBS7a/7b".to_string(),
            weight: fc_t,
            etiology: "UV light exposure (skin cancer)".to_string(),
        });
        summary_parts.push("UV mutagenesis (SBS7)");
    }

    // APOBEC (SBS2 + SBS13): C>T + C>G elevated
    if apobec_enriched && !uv_signature {
        let apobec_weight = (fc_t + fc_g) / 2.0;
        signatures.push(SignatureMatch {
            signature: "SBS2/SBS13".to_string(),
            weight: apobec_weight,
            etiology: "APOBEC3A/B cytidine deaminase activity".to_string(),
        });
        summary_parts.push("APOBEC mutagenesis (SBS2/SBS13)");
    }

    // Tobacco (SBS4): high C>A
    if tobacco_signature {
        signatures.push(SignatureMatch {
            signature: "SBS4".to_string(),
            weight: fc_a,
            etiology: "Tobacco smoking / polycyclic aromatic hydrocarbons".to_string(),
        });
        summary_parts.push("Tobacco smoking (SBS4)");
    }

    // MMR deficiency (SBS6/14/15/20): flat spectrum
    if mismatch_repair_deficiency {
        signatures.push(SignatureMatch {
            signature: "SBS6".to_string(),
            weight: 0.5,
            etiology: "Mismatch repair deficiency (MMR-D / MSI)".to_string(),
        });
        summary_parts.push("MMR deficiency (SBS6)");
    }

    // Aging (SBS1 + SBS5): C>T dominant, not APOBEC pattern
    if fc_t > 0.30 && !apobec_enriched && !uv_signature {
        let aging_weight = fc_t * 0.7;
        signatures.push(SignatureMatch {
            signature: "SBS1/SBS5".to_string(),
            weight: aging_weight,
            etiology: "Aging / spontaneous CpG deamination".to_string(),
        });
        summary_parts.push("Aging / CpG deamination (SBS1/SBS5)");
    }

    // ROS / oxidative stress (SBS18): high C>A without tobacco
    if fc_a > 0.20 && ft_c > 0.20 && !tobacco_signature {
        signatures.push(SignatureMatch {
            signature: "SBS18".to_string(),
            weight: (fc_a + ft_c) / 2.0,
            etiology: "Reactive oxygen species (ROS) / oxidative stress".to_string(),
        });
        summary_parts.push("Oxidative stress (SBS18)");
    }

    // 5-FU treatment (SBS17): elevated T>G
    if ft_g > 0.18 {
        signatures.push(SignatureMatch {
            signature: "SBS17a/17b".to_string(),
            weight: ft_g,
            etiology: "Prior 5-fluorouracil (5-FU) chemotherapy".to_string(),
        });
        summary_parts.push("5-FU treatment (SBS17)");
    }

    // Sort by weight descending
    signatures.sort_by(|a, b| b.weight.partial_cmp(&a.weight).unwrap_or(std::cmp::Ordering::Equal));

    let summary = if summary_parts.is_empty() {
        "No dominant mutational signature detected (low SNV count or mixed pattern)".to_string()
    } else {
        summary_parts.join("; ")
    };

    let note = if spec.total_snvs < 50 {
        Some(format!(
            "Only {} SNVs in spectrum — COSMIC signature attribution may be unreliable (≥50 recommended)",
            spec.total_snvs
        ))
    } else {
        None
    };

    MutationalSignatureResult {
        spectrum_6ch: spec.clone(),
        dominant_signatures: signatures,
        summary,
        apobec_enriched,
        tobacco_signature,
        mismatch_repair_deficiency,
        uv_signature,
        note,
    }
}

/// Compute mutational signature result from variant list.
pub fn compute_mutational_signatures(variants: &[VariantRecord]) -> MutationalSignatureResult {
    let spec = compute_sbs6_spectrum(variants);
    detect_signatures_from_6ch(&spec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::TiTvClass;

    fn snv(chrom: &str, pos: u64, r: &str, a: &str) -> VariantRecord {
        VariantRecord {
            chrom: chrom.to_string(),
            pos,
            ref_allele: r.to_string(),
            alt_allele: a.to_string(),
            qual: 50.0,
            titv: TiTvClass::Transition,
            af: None,
            gene: None,
        }
    }

    #[test]
    fn empty_variants_yield_zero_spectrum() {
        let spec = compute_sbs6_spectrum(&[]);
        assert_eq!(spec.total_snvs, 0);
        assert_eq!(spec.fractions, [0.0; 6]);
    }

    #[test]
    fn indels_and_mnps_are_excluded_from_spectrum() {
        let variants = vec![
            snv("chr1", 100, "C", "CA"), // insertion
            snv("chr1", 200, "AT", "A"), // deletion
            snv("chr1", 300, "AC", "GT"), // MNP
        ];
        let spec = compute_sbs6_spectrum(&variants);
        assert_eq!(spec.total_snvs, 0);
    }

    #[test]
    fn purine_substitutions_are_normalized_to_pyrimidine_context() {
        // G>T on the reference strand is the complement of C>A.
        let variants = vec![snv("chr1", 100, "G", "T")];
        let spec = compute_sbs6_spectrum(&variants);
        assert_eq!(spec.total_snvs, 1);
        assert_eq!(spec.c_to_a, 1);
        assert_eq!(spec.c_to_g, 0);
    }

    #[test]
    fn identical_ref_and_alt_are_skipped() {
        let variants = vec![snv("chr1", 100, "C", "C")];
        let spec = compute_sbs6_spectrum(&variants);
        assert_eq!(spec.total_snvs, 0);
    }

    #[test]
    fn fractions_sum_to_one_when_spectrum_nonempty() {
        let variants = vec![
            snv("chr1", 100, "C", "A"),
            snv("chr1", 200, "C", "T"),
            snv("chr1", 300, "T", "C"),
            snv("chr1", 400, "T", "G"),
        ];
        let spec = compute_sbs6_spectrum(&variants);
        assert_eq!(spec.total_snvs, 4);
        let sum: f64 = spec.fractions.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9);
    }

    #[test]
    fn uv_signature_detected_from_extreme_c_to_t_enrichment() {
        let mut variants = Vec::new();
        for i in 0..80 {
            variants.push(snv("chr1", i, "C", "T"));
        }
        for i in 80..100 {
            variants.push(snv("chr1", i, "C", "A"));
        }
        let result = compute_mutational_signatures(&variants);
        assert!(result.uv_signature);
        assert!(result.dominant_signatures.iter().any(|s| s.signature == "SBS7a/7b"));
    }

    #[test]
    fn tobacco_signature_detected_from_high_c_to_a() {
        let mut variants = Vec::new();
        for i in 0..40 {
            variants.push(snv("chr1", i, "C", "A"));
        }
        for i in 40..100 {
            variants.push(snv("chr1", i, "T", "C"));
        }
        let result = compute_mutational_signatures(&variants);
        assert!(result.tobacco_signature);
        assert!(result.dominant_signatures.iter().any(|s| s.signature == "SBS4"));
    }

    #[test]
    fn low_snv_count_produces_a_confidence_note() {
        let variants = vec![snv("chr1", 100, "C", "A"); 10];
        let result = compute_mutational_signatures(&variants);
        assert!(result.note.is_some());
        assert!(result.note.unwrap().contains("10 SNVs"));
    }

    #[test]
    fn no_dominant_pattern_gives_default_summary() {
        let spec = SbsSpectrum6 {
            c_to_a: 0,
            c_to_g: 0,
            c_to_t: 0,
            t_to_a: 0,
            t_to_c: 0,
            t_to_g: 0,
            total_snvs: 0,
            fractions: [0.0; 6],
        };
        let result = detect_signatures_from_6ch(&spec);
        assert!(result.dominant_signatures.is_empty());
        assert!(result.summary.contains("No dominant"));
    }

    #[test]
    fn signatures_are_sorted_by_weight_descending() {
        let mut variants = Vec::new();
        for i in 0..50 {
            variants.push(snv("chr1", i, "T", "G")); // SBS17
        }
        for i in 50..70 {
            variants.push(snv("chr1", i, "C", "A")); // contributes to SBS18/SBS4
        }
        for i in 70..90 {
            variants.push(snv("chr1", i, "T", "C"));
        }
        let result = compute_mutational_signatures(&variants);
        for pair in result.dominant_signatures.windows(2) {
            assert!(pair[0].weight >= pair[1].weight);
        }
    }
}
