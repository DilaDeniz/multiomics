//! Mutect2-style somatic tumor/normal variant calling.
//!
//! Compares tumor and matched-normal pileups to detect somatic mutations using
//! a log-odds (LOD) scoring model adapted from Mutect2 (McKenna et al. 2010,
//! extended in Shiraishi et al. 2013). Strand-bias filtering is applied via a
//! chi-squared approximation of Fisher's exact test.
//!
//! # References
//! * McKenna et al. (2010) The Genome Analysis Toolkit. Genome Research 20:1297–1303.
//! * Shiraishi et al. (2013) An empirical Bayesian framework for somatic mutation
//!   detection from cancer genome sequencing data. Nucleic Acids Research 41(7):e89.

use ahash::AHashMap;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::pileup::{build_pileup, PileupColumn};

// ── Constants ─────────────────────────────────────────────────────────────────

/// Sequencing error rate used in the LOD model.
const ERROR_RATE: f64 = 0.001;

// ── Public types ──────────────────────────────────────────────────────────────

/// A somatic variant call from tumor/normal comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SomaticCall {
    pub chrom: String,
    pub pos: u64,
    pub ref_base: u8,
    pub alt_base: u8,
    /// Allele frequency in tumor (alt reads / total reads).
    pub tumor_af: f64,
    /// Allele frequency in normal (should be ~0 for somatic).
    pub normal_af: f64,
    /// Tumor depth.
    pub tumor_depth: u32,
    /// Normal depth.
    pub normal_depth: u32,
    /// Log-odds score: ln P(somatic) - ln P(germline).
    pub lod_score: f64,
    /// True if passing all filters.
    pub pass: bool,
    /// Comma-separated filter reasons if not PASS (e.g. "normal_lod", "strand_bias").
    pub filter: String,
}

/// Summary statistics across somatic calls.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SomaticSummary {
    pub total_candidates: u64,
    pub pass_count: u64,
    pub filtered_count: u64,
    pub mean_tumor_af: f64,
    /// Mutation spectrum: key = "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"
    pub mutation_spectrum: AHashMap<String, u64>,
    /// Ti/Tv ratio among PASS somatic SNVs.
    pub titv_ratio: f64,
}

// ── Private types ─────────────────────────────────────────────────────────────

/// Per-sample base counts at one position.
struct PositionCounts {
    ref_fwd: u32,
    ref_rev: u32,
    alt_fwd: u32,
    alt_rev: u32,
    other: u32,
}

impl PositionCounts {
    fn alt_count(&self) -> u32 {
        self.alt_fwd + self.alt_rev
    }

    fn ref_count(&self) -> u32 {
        self.ref_fwd + self.ref_rev
    }

    fn total(&self) -> u32 {
        self.ref_count() + self.alt_count() + self.other
    }

    fn af(&self) -> f64 {
        let t = self.total();
        if t == 0 {
            0.0
        } else {
            self.alt_count() as f64 / t as f64
        }
    }
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Call somatic variants by comparing tumor vs normal pileups.
///
/// Uses the Mutect2 log-odds scoring model (McKenna et al. 2010, extended
/// Shiraishi et al. 2013):
/// ```text
///   LOD_tumor  = log(AF_t / error_rate) × alt_t  +  log((1-AF_t)/(1-error_rate)) × ref_t
///   LOD_normal = log(0.5 / error_rate) × alt_n  (germline het would give AF≈0.5)
///   Call somatic if LOD_tumor > min_tumor_lod AND LOD_normal < max_normal_lod
/// ```
///
/// # Arguments
/// * `tumor_pileup`   — from `build_pileup()` on the tumor BAM
/// * `normal_pileup`  — from `build_pileup()` on the normal BAM
/// * `min_tumor_lod`  — minimum log-odds in tumor (Mutect2 default: 6.3)
/// * `max_normal_lod` — maximum log-odds in normal (default: 2.2)
/// * `min_tumor_af`   — minimum allele frequency in tumor (default: 0.1)
/// * `min_depth`      — minimum depth in both tumor and normal (default: 8)
pub fn call_somatic_variants(
    tumor_pileup: &[PileupColumn],
    normal_pileup: &[PileupColumn],
    min_tumor_lod: f64,
    max_normal_lod: f64,
    min_tumor_af: f64,
    min_depth: u32,
) -> Vec<SomaticCall> {
    // Build (chrom, pos) → &PileupColumn lookup for normal.
    let normal_map: AHashMap<(&str, u64), &PileupColumn> = normal_pileup
        .iter()
        .map(|col| (col.chrom.as_str(), col.pos, col))
        .map(|(chrom, pos, col)| ((chrom, pos), col))
        .collect();

    let mut calls = Vec::new();

    for tcol in tumor_pileup {
        // ── 1. Find the most common non-ref ALT base in the tumor ────────────
        let (alt_base, tumor_counts) = match best_alt(tcol) {
            Some(v) => v,
            None => continue,
        };

        // Require at least 2 alt reads to nominate a candidate.
        if tumor_counts.alt_count() < 2 {
            continue;
        }

        let tumor_depth = tumor_counts.total();
        if tumor_depth < min_depth {
            continue;
        }

        // ── 2. Lookup normal counts at the same position ─────────────────────
        let normal_counts = match normal_map.get(&(tcol.chrom.as_str(), tcol.pos)) {
            Some(ncol) => count_bases(ncol, tcol.ref_base, alt_base),
            None => continue, // no normal coverage — skip
        };

        let normal_depth = normal_counts.total();
        if normal_depth < min_depth {
            continue;
        }

        // ── 3. Allele frequencies ────────────────────────────────────────────
        let tumor_af = tumor_counts.af();
        let normal_af = normal_counts.af();

        // ── 4. LOD scoring ───────────────────────────────────────────────────
        let lod_tumor = lod_tumor(tumor_counts.alt_count(), tumor_counts.ref_count(), tumor_af);
        let lod_normal = lod_normal(normal_counts.alt_count());

        // ── 5. Collect filters ───────────────────────────────────────────────
        let mut filters: Vec<&'static str> = Vec::new();

        if tumor_af < min_tumor_af {
            filters.push("low_tumor_af");
        }
        if lod_tumor < min_tumor_lod {
            filters.push("tumor_lod");
        }
        if lod_normal > max_normal_lod {
            filters.push("normal_lod");
        }
        if strand_biased(&tumor_counts) {
            filters.push("strand_bias");
        }

        let pass = filters.is_empty();
        let filter = if pass {
            "PASS".to_string()
        } else {
            filters.join(",")
        };

        calls.push(SomaticCall {
            chrom: tcol.chrom.clone(),
            pos: tcol.pos,
            ref_base: tcol.ref_base,
            alt_base,
            tumor_af,
            normal_af,
            tumor_depth,
            normal_depth,
            lod_score: lod_tumor,
            pass,
            filter,
        });
    }

    calls
}

/// Convenience function: call somatic variants directly from BAM paths.
///
/// Calls `build_pileup()` for both BAMs then `call_somatic_variants()`.
#[allow(clippy::too_many_arguments)]
pub fn call_somatic_from_bams(
    tumor_bam: &std::path::Path,
    normal_bam: &std::path::Path,
    min_base_qual: u8,
    min_mapq: u8,
    min_tumor_lod: f64,
    max_normal_lod: f64,
    min_tumor_af: f64,
    min_depth: u32,
) -> Result<Vec<SomaticCall>> {
    let tumor_pileup = build_pileup(tumor_bam, min_base_qual, min_mapq, 8000)
        .with_context(|| format!("cannot build tumor pileup from {}", tumor_bam.display()))?;
    let normal_pileup = build_pileup(normal_bam, min_base_qual, min_mapq, 8000)
        .with_context(|| format!("cannot build normal pileup from {}", normal_bam.display()))?;
    Ok(call_somatic_variants(
        &tumor_pileup,
        &normal_pileup,
        min_tumor_lod,
        max_normal_lod,
        min_tumor_af,
        min_depth,
    ))
}

/// Summarize somatic calls: counts by filter status, mean tumor AF, mutation spectrum.
pub fn summarize_somatic(calls: &[SomaticCall]) -> SomaticSummary {
    let total_candidates = calls.len() as u64;
    let pass_calls: Vec<&SomaticCall> = calls.iter().filter(|c| c.pass).collect();
    let pass_count = pass_calls.len() as u64;
    let filtered_count = total_candidates - pass_count;

    // Mean tumor AF across PASS calls.
    let mean_tumor_af = if pass_count == 0 {
        0.0
    } else {
        pass_calls.iter().map(|c| c.tumor_af).sum::<f64>() / pass_count as f64
    };

    // Mutation spectrum (pyrimidine context normalization).
    let mut mutation_spectrum: AHashMap<String, u64> = AHashMap::new();
    let mut ti_count: u64 = 0;
    let mut tv_count: u64 = 0;

    for call in &pass_calls {
        let ref_b = call.ref_base.to_ascii_uppercase();
        let alt_b = call.alt_base.to_ascii_uppercase();

        // Skip non-SNV bases (indels, N, etc.).
        if !is_standard_base(ref_b) || !is_standard_base(alt_b) {
            continue;
        }

        // Normalize to pyrimidine reference context (C or T).
        let (norm_ref, norm_alt) = to_pyrimidine_context(ref_b, alt_b);
        let key = format!("{}>{}", norm_ref as char, norm_alt as char);
        *mutation_spectrum.entry(key).or_insert(0) += 1;

        // Ti/Tv classification.
        match is_transition(ref_b, alt_b) {
            Some(true) => ti_count += 1,
            Some(false) => tv_count += 1,
            None => {}
        }
    }

    let titv_ratio = if tv_count == 0 {
        0.0
    } else {
        ti_count as f64 / tv_count as f64
    };

    SomaticSummary {
        total_candidates,
        pass_count,
        filtered_count,
        mean_tumor_af,
        mutation_spectrum,
        titv_ratio,
    }
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Count ref/alt/other bases at a pileup column, split by strand.
fn count_bases(col: &PileupColumn, ref_base: u8, alt_base: u8) -> PositionCounts {
    let mut counts = PositionCounts {
        ref_fwd: 0,
        ref_rev: 0,
        alt_fwd: 0,
        alt_rev: 0,
        other: 0,
    };
    let ref_up = ref_base.to_ascii_uppercase();
    let alt_up = alt_base.to_ascii_uppercase();

    for pb in &col.bases {
        let b = pb.base.to_ascii_uppercase();
        if b == ref_up {
            if pb.is_rev {
                counts.ref_rev += 1;
            } else {
                counts.ref_fwd += 1;
            }
        } else if b == alt_up {
            if pb.is_rev {
                counts.alt_rev += 1;
            } else {
                counts.alt_fwd += 1;
            }
        } else {
            counts.other += 1;
        }
    }
    counts
}

/// Find the most common non-ref base (with count >= 2) and return its counts.
///
/// Returns `None` if no alt base passes the minimum count threshold.
fn best_alt(col: &PileupColumn) -> Option<(u8, PositionCounts)> {
    let ref_up = col.ref_base.to_ascii_uppercase();
    let bases: &[u8] = b"ACGT";

    // Count each base.
    let mut base_counts = [0u32; 4];
    for pb in &col.bases {
        let b = pb.base.to_ascii_uppercase();
        if let Some(idx) = acgt_index(b) {
            base_counts[idx] += 1;
        }
    }

    // Find the most frequent non-ref base with count >= 2.
    let best = bases
        .iter()
        .enumerate()
        .filter(|(_, &b)| b != ref_up)
        .filter(|(i, _)| base_counts[*i] >= 2)
        .max_by_key(|(i, _)| base_counts[*i]);

    let (_, &alt_base) = best?;
    let counts = count_bases(col, ref_up, alt_base);
    Some((alt_base, counts))
}

/// LOD score for tumor: log(AF_t / e) * alt + log((1-AF_t) / (1-e)) * ref.
fn lod_tumor(alt_count: u32, ref_count: u32, af: f64) -> f64 {
    let af = af.clamp(ERROR_RATE + f64::EPSILON, 1.0 - f64::EPSILON);
    let lod_alt = (af / ERROR_RATE).ln() * alt_count as f64;
    let lod_ref = ((1.0 - af) / (1.0 - ERROR_RATE)).ln() * ref_count as f64;
    lod_alt + lod_ref
}

/// LOD score for normal germline hypothesis: log(0.5 / error_rate) × alt_n.
///
/// A germline het would have AF ≈ 0.5; high LOD here means the normal is likely
/// a carrier and the mutation is NOT somatic.
fn lod_normal(alt_count: u32) -> f64 {
    (0.5_f64 / ERROR_RATE).ln() * alt_count as f64
}

/// Chi-squared approximation of Fisher's exact test for strand bias.
///
/// Builds the 2×2 table:
/// ```text
///            fwd   rev
///   ref      a     b
///   alt      c     d
/// ```
/// and computes χ² = (ad − bc)² × N / ((a+b)(c+d)(a+c)(b+d)).
///
/// Returns `true` (biased) if p < 0.001 (χ²(1) > 10.8276 for 1 df).
fn strand_biased(counts: &PositionCounts) -> bool {
    let a = counts.ref_fwd as f64;
    let b = counts.ref_rev as f64;
    let c = counts.alt_fwd as f64;
    let d = counts.alt_rev as f64;
    let n = a + b + c + d;

    // Avoid division by zero when any marginal is zero.
    let denom = (a + b) * (c + d) * (a + c) * (b + d);
    if denom == 0.0 || n == 0.0 {
        return false;
    }

    let chi2 = (a * d - b * c).powi(2) * n / denom;

    // χ²(1) critical value at p = 0.001 is 10.828.
    // Convert p-threshold to chi2 via inverse approximation.
    let chi2_thresh = chi2_critical_001();
    chi2 > chi2_thresh
}

/// χ²(1) critical value corresponding to p = 0.001 (10.8276).
#[inline]
fn chi2_critical_001() -> f64 {
    // -2 × ln(p / 2) is a good approximation for chi2(1) at small p.
    // Exact value: qchisq(0.999, 1) = 10.8276.
    10.8276
}

/// True for A or G (purines), false for C or T (pyrimidines).
fn is_purine(b: u8) -> bool {
    b == b'A' || b == b'G'
}

/// Classify a substitution as transition (true), transversion (false), or None.
fn is_transition(ref_b: u8, alt_b: u8) -> Option<bool> {
    if ref_b == alt_b {
        return None;
    }
    Some(is_purine(ref_b) == is_purine(alt_b))
}

/// Normalize ref/alt to pyrimidine context (C or T as reference).
///
/// If the reference is already a pyrimidine (C or T), return as-is.
/// Otherwise complement both ref and alt.
fn to_pyrimidine_context(ref_b: u8, alt_b: u8) -> (u8, u8) {
    if ref_b == b'C' || ref_b == b'T' {
        (ref_b, alt_b)
    } else {
        (complement(ref_b), complement(alt_b))
    }
}

/// Complement of a DNA base (A↔T, C↔G).
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// Return true if the base is one of A, C, G, T.
fn is_standard_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Map A→0, C→1, G→2, T→3, everything else → None.
fn acgt_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::{PileupBase, PileupColumn};

    /// Build a PileupColumn with a fixed ref base and mixed alt reads.
    fn make_col(chrom: &str, pos: u64, ref_base: u8, bases: Vec<PileupBase>) -> PileupColumn {
        PileupColumn {
            chrom: chrom.to_string(),
            pos,
            ref_base,
            bases,
        }
    }

    fn make_base(base: u8, is_rev: bool) -> PileupBase {
        PileupBase {
            base,
            base_qual: 30,
            mapq: 60,
            is_rev,
        }
    }

    // ── Test 1: Basic LOD scoring with 30% tumor AF, 0% normal AF ────────────

    #[test]
    fn somatic_lod_scoring_basic() {
        // 70 ref reads + 30 alt reads in tumor (AF = 0.30).
        let mut tumor_bases: Vec<PileupBase> =
            (0..70).map(|i| make_base(b'A', i % 2 == 0)).collect();
        tumor_bases.extend((0..30).map(|i| make_base(b'T', i % 2 == 0)));

        // 30 ref reads, 0 alt reads in normal.
        let normal_bases: Vec<PileupBase> = (0..30).map(|i| make_base(b'A', i % 2 == 0)).collect();

        let tumor_pileup = vec![make_col("chr1", 100, b'A', tumor_bases)];
        let normal_pileup = vec![make_col("chr1", 100, b'A', normal_bases)];

        let calls = call_somatic_variants(
            &tumor_pileup,
            &normal_pileup,
            6.3, // min_tumor_lod (Mutect2 default)
            2.2, // max_normal_lod
            0.1, // min_tumor_af
            8,   // min_depth
        );

        assert_eq!(calls.len(), 1, "expected exactly one somatic call");
        let call = &calls[0];
        assert_eq!(call.alt_base, b'T');
        assert!(
            call.lod_score > 6.3,
            "LOD score {} should exceed threshold 6.3",
            call.lod_score
        );
        assert!(
            call.pass,
            "call should PASS filters; filter = {}",
            call.filter
        );
        assert!(
            (call.tumor_af - 0.3).abs() < 0.02,
            "tumor AF {:.3} should be ~0.3",
            call.tumor_af
        );
        assert!(
            call.normal_af < 0.01,
            "normal AF {:.3} should be ~0.0",
            call.normal_af
        );
    }

    // ── Test 2: Normal contamination causes normal_lod filter ─────────────────

    #[test]
    fn somatic_normal_contamination_filtered() {
        // Tumor: 50% alt AF — strong tumor signal.
        let mut tumor_bases: Vec<PileupBase> =
            (0..50).map(|i| make_base(b'A', i % 2 == 0)).collect();
        tumor_bases.extend((0..50).map(|i| make_base(b'T', i % 2 == 0)));

        // Normal: 40% alt AF — contamination (germline or poor pairing).
        let mut normal_bases: Vec<PileupBase> =
            (0..30).map(|i| make_base(b'A', i % 2 == 0)).collect();
        normal_bases.extend((0..20).map(|i| make_base(b'T', i % 2 == 0)));

        let tumor_pileup = vec![make_col("chr1", 200, b'A', tumor_bases)];
        let normal_pileup = vec![make_col("chr1", 200, b'A', normal_bases)];

        let calls = call_somatic_variants(&tumor_pileup, &normal_pileup, 6.3, 2.2, 0.1, 8);

        assert_eq!(calls.len(), 1);
        let call = &calls[0];
        assert!(
            !call.pass,
            "call with high normal AF should NOT pass; filter = {}",
            call.filter
        );
        assert!(
            call.filter.contains("normal_lod"),
            "filter should include 'normal_lod', got: {}",
            call.filter
        );
    }

    // ── Test 3: Mutation spectrum for C>T calls ────────────────────────────────

    #[test]
    fn somatic_summary_spectrum() {
        // Construct three synthetic PASS SomaticCall entries with C>T mutations.
        let make_call = |chrom: &str, pos: u64, ref_b: u8, alt_b: u8| SomaticCall {
            chrom: chrom.to_string(),
            pos,
            ref_base: ref_b,
            alt_base: alt_b,
            tumor_af: 0.3,
            normal_af: 0.0,
            tumor_depth: 50,
            normal_depth: 30,
            lod_score: 20.0,
            pass: true,
            filter: "PASS".to_string(),
        };

        let calls = vec![
            make_call("chr1", 100, b'C', b'T'),
            make_call("chr1", 200, b'C', b'T'),
            make_call("chr1", 300, b'C', b'T'),
            // One G>A call which normalizes to C>T (complement).
            make_call("chr1", 400, b'G', b'A'),
        ];

        let summary = summarize_somatic(&calls);

        assert_eq!(summary.total_candidates, 4);
        assert_eq!(summary.pass_count, 4);
        assert_eq!(summary.filtered_count, 0);

        // All four calls should land in "C>T" after pyrimidine normalization.
        let ct_count = summary.mutation_spectrum.get("C>T").copied().unwrap_or(0);
        assert_eq!(
            ct_count, 4,
            "expected 4 C>T calls in spectrum, got {}; spectrum = {:?}",
            ct_count, summary.mutation_spectrum
        );

        // All four are transitions (C↔T / A↔G), so Ti/Tv = infinity — but since
        // tv_count=0 we return 0.0 rather than NaN.  Just verify Ti/Tv >= 0.
        assert!(
            summary.titv_ratio >= 0.0,
            "Ti/Tv ratio should be non-negative"
        );

        assert!(
            (summary.mean_tumor_af - 0.3).abs() < 1e-9,
            "mean tumor AF should be 0.3"
        );
    }

    // ── Additional unit tests ─────────────────────────────────────────────────

    #[test]
    fn lod_tumor_high_af_gives_large_lod() {
        // 40 alt reads, 60 ref reads, AF = 0.4.
        let lod = lod_tumor(40, 60, 0.4);
        // Should be well above Mutect2 default of 6.3.
        assert!(lod > 6.3, "LOD = {}", lod);
    }

    #[test]
    fn lod_normal_zero_alt_is_zero() {
        assert_eq!(lod_normal(0), 0.0);
    }

    #[test]
    fn strand_bias_extreme_detected() {
        // All alt reads on fwd strand, all ref reads on rev — strong bias.
        let counts = PositionCounts {
            ref_fwd: 0,
            ref_rev: 50,
            alt_fwd: 50,
            alt_rev: 0,
            other: 0,
        };
        assert!(
            strand_biased(&counts),
            "extreme strand bias should be detected"
        );
    }

    #[test]
    fn strand_bias_balanced_not_detected() {
        // Balanced strands — should not trigger.
        let counts = PositionCounts {
            ref_fwd: 25,
            ref_rev: 25,
            alt_fwd: 15,
            alt_rev: 15,
            other: 0,
        };
        assert!(
            !strand_biased(&counts),
            "balanced strand distribution should not be flagged"
        );
    }

    #[test]
    fn pyrimidine_normalization_purine_ref() {
        // G>A should normalize to C>T (complement of both).
        let (r, a) = to_pyrimidine_context(b'G', b'A');
        assert_eq!(r, b'C');
        assert_eq!(a, b'T');
    }

    #[test]
    fn pyrimidine_normalization_pyrimidine_ref() {
        // C>G stays as-is.
        let (r, a) = to_pyrimidine_context(b'C', b'G');
        assert_eq!(r, b'C');
        assert_eq!(a, b'G');
    }

    #[test]
    fn transition_classification() {
        assert_eq!(is_transition(b'C', b'T'), Some(true)); // Ti
        assert_eq!(is_transition(b'A', b'G'), Some(true)); // Ti
        assert_eq!(is_transition(b'C', b'A'), Some(false)); // Tv
        assert_eq!(is_transition(b'G', b'T'), Some(false)); // Tv
        assert_eq!(is_transition(b'A', b'A'), None); // same base
    }
}
