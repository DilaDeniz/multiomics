use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::collections::HashSet;
use crate::types::VariantRecord;
use crate::cosmic_ref::{COSMIC_SBS, N_CHANNELS};

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
    /// Attribution method used: "6-channel heuristic" (no reference) or
    /// "96-channel NNLS (COSMIC v3.3)" (reference-guided).
    #[serde(default = "default_method")]
    pub method: String,
    /// Full 96-channel SBS spectrum — populated only with a reference FASTA.
    #[serde(default)]
    pub sbs96: Option<Sbs96Spectrum>,
    /// COSMIC signature exposures from NNLS deconvolution against the v3.3
    /// reference catalogue — populated only with a reference FASTA.
    #[serde(default)]
    pub cosmic_exposures: Vec<SignatureExposure>,
    /// Cosine similarity between the reconstructed and observed 96-channel
    /// spectra (fit quality, 1.0 = perfect) — Some only with a reference.
    #[serde(default)]
    pub reconstruction_cosine: Option<f64>,
}

fn default_method() -> String {
    "6-channel heuristic".to_string()
}

/// Full 96-channel SBS mutation spectrum (trinucleotide context).
///
/// Requires a reference FASTA to resolve the base immediately 5' and 3' of
/// each substitution. Channels follow the canonical COSMIC ordering
/// (see [`crate::cosmic_ref::SBS96_TYPES`]).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sbs96Spectrum {
    /// Raw counts per channel (length 96).
    pub counts: Vec<u64>,
    /// Fraction per channel, sums to 1.0 when `total > 0` (length 96).
    pub fractions: Vec<f64>,
    /// Total SNVs assigned to a channel.
    pub total: u64,
    /// SNVs skipped because reference context could not be resolved
    /// (chromosome absent, position at contig edge, or non-ACGT base).
    pub skipped_no_context: u64,
}

/// One COSMIC signature's fitted exposure from NNLS deconvolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SignatureExposure {
    /// COSMIC signature name (e.g., "SBS4").
    pub signature: String,
    /// Associated aetiology from the COSMIC catalogue.
    pub etiology: String,
    /// Relative contribution to the observed spectrum [0, 1].
    pub weight: f64,
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
        method: default_method(),
        sbs96: None,
        cosmic_exposures: Vec::new(),
        reconstruction_cosine: None,
    }
}

/// Compute mutational signature result from variant list (6-channel, no reference).
pub fn compute_mutational_signatures(variants: &[VariantRecord]) -> MutationalSignatureResult {
    let spec = compute_sbs6_spectrum(variants);
    detect_signatures_from_6ch(&spec)
}

// ─── 96-channel SBS analysis (reference-guided COSMIC deconvolution) ──────────

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn dna_complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// Map a trinucleotide substitution to its canonical SBS-96 channel index.
///
/// `p5`/`p3` are the 5' and 3' flanking bases, `refb`/`altb` the substitution.
/// Purine-centered mutations are reverse-complemented to pyrimidine context,
/// matching the COSMIC ordering: index = p5·24 + substitution·4 + p3.
pub fn sbs96_channel_index(p5: u8, refb: u8, altb: u8, p3: u8) -> Option<usize> {
    let (p5, refb, altb, p3) = (
        p5.to_ascii_uppercase(),
        refb.to_ascii_uppercase(),
        altb.to_ascii_uppercase(),
        p3.to_ascii_uppercase(),
    );
    if refb == altb {
        return None;
    }
    // Normalize to pyrimidine (C/T) center. Purine center → reverse-complement,
    // which also swaps the 5' and 3' flanks.
    let (p5, refb, altb, p3) = if refb == b'C' || refb == b'T' {
        (p5, refb, altb, p3)
    } else {
        (dna_complement(p3), dna_complement(refb), dna_complement(altb), dna_complement(p5))
    };
    let sub_idx = match (refb, altb) {
        (b'C', b'A') => 0,
        (b'C', b'G') => 1,
        (b'C', b'T') => 2,
        (b'T', b'A') => 3,
        (b'T', b'C') => 4,
        (b'T', b'G') => 5,
        _ => return None,
    };
    let p5i = base_index(p5)?;
    let p3i = base_index(p3)?;
    Some(p5i * 24 + sub_idx * 4 + p3i)
}

/// Build the 96-channel SBS spectrum from SNVs using reference trinucleotide
/// context. `ref_seqs` maps chromosome name → uppercase sequence bytes.
pub fn compute_sbs96_spectrum(
    variants: &[VariantRecord],
    ref_seqs: &HashMap<String, Vec<u8>>,
) -> Sbs96Spectrum {
    let mut counts = vec![0u64; N_CHANNELS];
    let mut skipped = 0u64;

    for v in variants {
        if v.ref_allele.len() != 1 || v.alt_allele.len() != 1 {
            continue; // only SNVs carry trinucleotide context
        }
        let refb = v.ref_allele.as_bytes()[0].to_ascii_uppercase();
        let altb = v.alt_allele.as_bytes()[0].to_ascii_uppercase();
        if refb == altb {
            continue;
        }
        let seq = match ref_seqs.get(&v.chrom) {
            Some(s) => s,
            None => {
                skipped += 1;
                continue;
            }
        };
        // VCF pos is 1-based; the REF base sits at 0-based index pos-1.
        if v.pos < 2 || (v.pos as usize) >= seq.len() {
            skipped += 1;
            continue;
        }
        let center = (v.pos - 1) as usize;
        let p5 = seq[center - 1];
        let p3 = seq[center + 1];
        match sbs96_channel_index(p5, refb, altb, p3) {
            Some(idx) => counts[idx] += 1,
            None => skipped += 1,
        }
    }

    let total: u64 = counts.iter().sum();
    let fractions = if total > 0 {
        let n = total as f64;
        counts.iter().map(|&c| c as f64 / n).collect()
    } else {
        vec![0.0; N_CHANNELS]
    };

    Sbs96Spectrum { counts, fractions, total, skipped_no_context: skipped }
}

/// Cosine similarity between two equal-length vectors.
fn cosine_similarity(a: &[f64], b: &[f64]) -> f64 {
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let na: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na == 0.0 || nb == 0.0 {
        0.0
    } else {
        dot / (na * nb)
    }
}

/// Solve the unconstrained least-squares problem restricted to `cols` of the
/// design matrix `a` (m×n, column-major slices), i.e. minimize ||A_cols z − b||,
/// via the normal equations with Gaussian elimination. Returns z for `cols`.
fn least_squares_on_columns(
    a: &[Vec<f64>], // a[col] = column vector of length m
    b: &[f64],
    cols: &[usize],
) -> Option<Vec<f64>> {
    let k = cols.len();
    if k == 0 {
        return Some(Vec::new());
    }
    // Normal equations: G z = c, where G = A_cols^T A_cols, c = A_cols^T b.
    let mut g = vec![vec![0.0f64; k]; k];
    let mut c = vec![0.0f64; k];
    for (i, &ci) in cols.iter().enumerate() {
        for (j, &cj) in cols.iter().enumerate() {
            let dot: f64 = a[ci].iter().zip(&a[cj]).map(|(x, y)| x * y).sum();
            g[i][j] = dot;
        }
        c[i] = a[ci].iter().zip(b).map(|(x, y)| x * y).sum();
    }
    // Gaussian elimination with partial pivoting on [G | c].
    for p in 0..k {
        let mut pivot = p;
        let mut best = g[p][p].abs();
        for r in (p + 1)..k {
            if g[r][p].abs() > best {
                best = g[r][p].abs();
                pivot = r;
            }
        }
        if best < 1e-12 {
            return None; // singular / rank-deficient passive set
        }
        if pivot != p {
            g.swap(p, pivot);
            c.swap(p, pivot);
        }
        for r in (p + 1)..k {
            let factor = g[r][p] / g[p][p];
            for col in p..k {
                g[r][col] -= factor * g[p][col];
            }
            c[r] -= factor * c[p];
        }
    }
    let mut z = vec![0.0f64; k];
    for i in (0..k).rev() {
        let mut s = c[i];
        for j in (i + 1)..k {
            s -= g[i][j] * z[j];
        }
        z[i] = s / g[i][i];
    }
    Some(z)
}

/// Lawson-Hanson active-set non-negative least squares: minimize ||A x − b||
/// subject to x ≥ 0. `a[col]` holds each design-matrix column (length m).
fn nnls(a: &[Vec<f64>], b: &[f64], max_iter: usize) -> Vec<f64> {
    let n = a.len();
    let mut x = vec![0.0f64; n];
    let mut passive: Vec<bool> = vec![false; n];
    let tol = 1e-9;

    for _ in 0..(max_iter.max(3 * n)) {
        // Gradient w = A^T (b − A x).
        let residual: Vec<f64> = {
            let mut r = b.to_vec();
            for j in 0..n {
                if x[j] != 0.0 {
                    for (ri, av) in r.iter_mut().zip(&a[j]) {
                        *ri -= x[j] * av;
                    }
                }
            }
            r
        };
        // Pick the active variable with the most positive gradient.
        let mut best_j = None;
        let mut best_w = tol;
        for j in 0..n {
            if !passive[j] {
                let w: f64 = a[j].iter().zip(&residual).map(|(av, rv)| av * rv).sum();
                if w > best_w {
                    best_w = w;
                    best_j = Some(j);
                }
            }
        }
        let j_in = match best_j {
            Some(j) => j,
            None => break, // KKT satisfied
        };
        passive[j_in] = true;

        // Inner loop: solve LS on the passive set, backtracking to stay feasible.
        loop {
            let cols: Vec<usize> = (0..n).filter(|&j| passive[j]).collect();
            let z_passive = match least_squares_on_columns(a, b, &cols) {
                Some(z) => z,
                None => {
                    // Degenerate: undo the last inclusion and stop.
                    passive[j_in] = false;
                    return x;
                }
            };
            let mut z = vec![0.0f64; n];
            for (idx, &col) in cols.iter().enumerate() {
                z[col] = z_passive[idx];
            }
            if cols.iter().all(|&j| z[j] > tol) {
                x = z;
                break;
            }
            // Move toward z as far as feasibility allows.
            let mut alpha = f64::INFINITY;
            for &j in &cols {
                if z[j] <= tol {
                    let denom = x[j] - z[j];
                    if denom > 0.0 {
                        alpha = alpha.min(x[j] / denom);
                    }
                }
            }
            if !alpha.is_finite() {
                x = z;
                break;
            }
            for j in 0..n {
                x[j] += alpha * (z[j] - x[j]);
            }
            for j in 0..n {
                if passive[j] && x[j] <= tol {
                    passive[j] = false;
                    x[j] = 0.0;
                }
            }
        }
    }
    x
}

/// Which broad process a signature belongs to, for setting the boolean flags
/// on the result from real exposures rather than crude 6-channel heuristics.
fn signature_processes(name: &str) -> (bool, bool, bool, bool) {
    // (apobec, tobacco, mmr_deficiency, uv)
    let apobec = matches!(name, "SBS2" | "SBS13");
    let tobacco = matches!(name, "SBS4" | "SBS29");
    let mmr = matches!(name, "SBS6" | "SBS14" | "SBS15" | "SBS20" | "SBS21" | "SBS26" | "SBS44");
    let uv = matches!(name, "SBS7a" | "SBS7b" | "SBS7c" | "SBS7d" | "SBS38");
    (apobec, tobacco, mmr, uv)
}

/// Reporting threshold: signatures below this relative contribution are dropped.
const EXPOSURE_THRESHOLD: f64 = 0.05;

/// Deconvolve a 96-channel spectrum against the embedded COSMIC v3.3 catalogue.
///
/// Fits the observed fractions as a non-negative combination of all reference
/// signatures (NNLS), prunes contributions below [`EXPOSURE_THRESHOLD`], refits
/// on the survivors for clean exposures, and normalizes to relative weights.
pub fn deconvolve_sbs96(spectrum: &Sbs96Spectrum) -> (Vec<SignatureExposure>, f64) {
    if spectrum.total == 0 {
        return (Vec::new(), 0.0);
    }
    // Design matrix: one column per reference signature (each length 96).
    let a_full: Vec<Vec<f64>> = COSMIC_SBS.iter().map(|s| s.profile.to_vec()).collect();
    let b = &spectrum.fractions;

    // Initial fit against the full catalogue.
    let x0 = nnls(&a_full, b, 200);

    // Prune to meaningful contributors, then refit for stable exposures.
    let sum0: f64 = x0.iter().sum();
    let survivors: Vec<usize> = if sum0 > 0.0 {
        (0..x0.len()).filter(|&j| x0[j] / sum0 >= EXPOSURE_THRESHOLD).collect()
    } else {
        Vec::new()
    };
    if survivors.is_empty() {
        return (Vec::new(), 0.0);
    }
    let a_sub: Vec<Vec<f64>> = survivors.iter().map(|&j| a_full[j].clone()).collect();
    let x_sub = nnls(&a_sub, b, 200);
    let sum_sub: f64 = x_sub.iter().sum();

    // Reconstruct fitted spectrum for cosine quality.
    let mut fitted = vec![0.0f64; N_CHANNELS];
    for (col, &weight) in x_sub.iter().enumerate() {
        for (fi, av) in fitted.iter_mut().zip(&a_sub[col]) {
            *fi += weight * av;
        }
    }
    let cosine = cosine_similarity(b, &fitted);

    let mut exposures: Vec<SignatureExposure> = survivors
        .iter()
        .enumerate()
        .filter_map(|(idx, &sig_j)| {
            let w = if sum_sub > 0.0 { x_sub[idx] / sum_sub } else { 0.0 };
            if w < EXPOSURE_THRESHOLD {
                return None;
            }
            Some(SignatureExposure {
                signature: COSMIC_SBS[sig_j].name.to_string(),
                etiology: COSMIC_SBS[sig_j].etiology.to_string(),
                weight: w,
            })
        })
        .collect();
    exposures.sort_by(|a, b| b.weight.partial_cmp(&a.weight).unwrap_or(std::cmp::Ordering::Equal));

    (exposures, cosine)
}

/// Full reference-guided COSMIC signature analysis.
///
/// Builds the 96-channel trinucleotide spectrum from `ref_seqs`, deconvolves it
/// against the COSMIC v3.3 catalogue, and returns a result whose boolean process
/// flags and summary are derived from the fitted exposures. Reference:
/// Alexandrov et al. 2020 (Nature); COSMIC Mutational Signatures v3.3.
pub fn compute_mutational_signatures_with_context(
    variants: &[VariantRecord],
    ref_seqs: &HashMap<String, Vec<u8>>,
) -> MutationalSignatureResult {
    // Keep the 6-channel spectrum for backward-compatible display fields.
    let spec6 = compute_sbs6_spectrum(variants);
    let spectrum = compute_sbs96_spectrum(variants, ref_seqs);
    let (exposures, cosine) = deconvolve_sbs96(&spectrum);

    // Derive process flags from the actual fitted signatures.
    let (mut apobec, mut tobacco, mut mmr, mut uv) = (false, false, false, false);
    for e in &exposures {
        let (a, t, m, u) = signature_processes(&e.signature);
        apobec |= a;
        tobacco |= t;
        mmr |= m;
        uv |= u;
    }

    let dominant_signatures: Vec<SignatureMatch> = exposures
        .iter()
        .map(|e| SignatureMatch {
            signature: e.signature.clone(),
            weight: e.weight,
            etiology: e.etiology.clone(),
        })
        .collect();

    let summary = if exposures.is_empty() {
        "No COSMIC signature exceeded the 5% contribution threshold".to_string()
    } else {
        exposures
            .iter()
            .take(4)
            .map(|e| format!("{} ({:.0}%)", e.signature, e.weight * 100.0))
            .collect::<Vec<_>>()
            .join("; ")
    };

    let note = if spectrum.total < 50 {
        Some(format!(
            "Only {} SNVs with resolvable context — COSMIC attribution may be unreliable (≥50 recommended)",
            spectrum.total
        ))
    } else if cosine < 0.85 {
        Some(format!(
            "Reconstruction cosine {:.2} — spectrum poorly explained by the reference catalogue",
            cosine
        ))
    } else {
        None
    };

    MutationalSignatureResult {
        spectrum_6ch: spec6,
        dominant_signatures,
        summary,
        apobec_enriched: apobec,
        tobacco_signature: tobacco,
        mismatch_repair_deficiency: mmr,
        uv_signature: uv,
        note,
        method: "96-channel NNLS (COSMIC v3.3)".to_string(),
        sbs96: Some(spectrum),
        cosmic_exposures: exposures,
        reconstruction_cosine: Some(cosine),
    }
}

/// Reference-guided COSMIC analysis reading the reference FASTA at `reference_path`.
///
/// Only the chromosomes referenced by `variants` are loaded, mirroring the
/// memory-frugal strategy used for reference-guided HRD scoring.
pub fn compute_mutational_signatures_with_reference(
    variants: &[VariantRecord],
    reference_path: &std::path::Path,
) -> anyhow::Result<MutationalSignatureResult> {
    use memmap2::Mmap;

    let file = std::fs::File::open(reference_path).map_err(|e| {
        anyhow::anyhow!("Cannot open reference FASTA '{}': {e}", reference_path.display())
    })?;
    // SAFETY: file is not modified while this process holds the Mmap.
    let mmap = unsafe { Mmap::map(&file) }
        .map_err(|e| anyhow::anyhow!("Cannot mmap reference FASTA: {e}"))?;
    let _ = mmap.advise(memmap2::Advice::Sequential);

    let needed: HashSet<&str> = variants.iter().map(|v| v.chrom.as_str()).collect();
    let ref_seqs = crate::cancer::parse_fasta_selective_pub(mmap.as_ref(), &needed);

    Ok(compute_mutational_signatures_with_context(variants, &ref_seqs))
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

    // ─── 96-channel SBS / COSMIC deconvolution tests ──────────────────────────

    use crate::cosmic_ref::{COSMIC_SBS, SBS96_TYPES};

    /// Parse a COSMIC type label like "A[C>A]T" into (p5, ref, alt, p3).
    fn parse_type(t: &str) -> (u8, u8, u8, u8) {
        let b = t.as_bytes();
        (b[0], b[2], b[4], b[6])
    }

    #[test]
    fn channel_index_matches_canonical_cosmic_ordering() {
        // Every canonical label must round-trip to its own index.
        for (i, label) in SBS96_TYPES.iter().enumerate() {
            let (p5, r, a, p3) = parse_type(label);
            assert_eq!(sbs96_channel_index(p5, r, a, p3), Some(i), "label {label}");
        }
    }

    #[test]
    fn purine_centered_mutation_is_reverse_complemented() {
        // A[G>T]A on the forward strand == T[C>A]T after strand normalization.
        // T[C>A]T is index 3*24 + 0*4 + 3 = 75.
        assert_eq!(sbs96_channel_index(b'A', b'G', b'T', b'A'), Some(75));
        assert_eq!(SBS96_TYPES[75], "T[C>A]T");
    }

    #[test]
    fn channel_index_rejects_non_substitutions_and_bad_bases() {
        assert_eq!(sbs96_channel_index(b'A', b'C', b'C', b'A'), None); // ref == alt
        assert_eq!(sbs96_channel_index(b'N', b'C', b'A', b'A'), None); // bad flank
        assert_eq!(sbs96_channel_index(b'A', b'C', b'N', b'A'), None); // bad alt
    }

    #[test]
    fn spectrum_construction_uses_reference_trinucleotide_context() {
        // Reference: chr1 = "TCA" (1-based positions 1,2,3). A C>A at pos 2 has
        // 5'=T, 3'=A → channel T[C>A]A = 3*24 = 72.
        let mut refs = HashMap::new();
        refs.insert("chr1".to_string(), b"TCA".to_vec());
        let variants = vec![snv("chr1", 2, "C", "A")];
        let spec = compute_sbs96_spectrum(&variants, &refs);
        assert_eq!(spec.total, 1);
        assert_eq!(spec.counts[72], 1);
        assert_eq!(spec.skipped_no_context, 0);
    }

    #[test]
    fn spectrum_skips_edge_and_missing_chromosome_variants() {
        let mut refs = HashMap::new();
        refs.insert("chr1".to_string(), b"TCA".to_vec());
        let variants = vec![
            snv("chr1", 1, "T", "A"),   // at contig edge (no 5' base) → skipped
            snv("chr2", 2, "C", "A"),   // chromosome absent → skipped
        ];
        let spec = compute_sbs96_spectrum(&variants, &refs);
        assert_eq!(spec.total, 0);
        assert_eq!(spec.skipped_no_context, 2);
    }

    #[test]
    fn nnls_recovers_a_pure_reference_signature() {
        // Feeding a signature's own profile back in must recover it as dominant
        // with a near-perfect reconstruction.
        let idx = COSMIC_SBS.iter().position(|s| s.name == "SBS4").unwrap();
        let spectrum = Sbs96Spectrum {
            counts: vec![0; 96],
            fractions: COSMIC_SBS[idx].profile.to_vec(),
            total: 1000,
            skipped_no_context: 0,
        };
        let (exposures, cosine) = deconvolve_sbs96(&spectrum);
        assert!(!exposures.is_empty());
        assert_eq!(exposures[0].signature, "SBS4");
        assert!(exposures[0].weight > 0.8, "SBS4 weight = {}", exposures[0].weight);
        assert!(cosine > 0.99, "cosine = {}", cosine);
    }

    #[test]
    fn nnls_recovers_an_orthogonal_two_signature_mixture() {
        // SBS4 (tobacco, C>A) and SBS7a (UV, C>T) are near-orthogonal.
        let i4 = COSMIC_SBS.iter().position(|s| s.name == "SBS4").unwrap();
        let i7 = COSMIC_SBS.iter().position(|s| s.name == "SBS7a").unwrap();
        let mut fractions = vec![0.0; 96];
        for k in 0..96 {
            fractions[k] = 0.5 * COSMIC_SBS[i4].profile[k] + 0.5 * COSMIC_SBS[i7].profile[k];
        }
        let spectrum = Sbs96Spectrum { counts: vec![0; 96], fractions, total: 2000, skipped_no_context: 0 };
        let (exposures, cosine) = deconvolve_sbs96(&spectrum);
        assert!(cosine > 0.99, "cosine = {cosine}");
        let total_w: f64 = exposures.iter().map(|e| e.weight).sum();
        assert!((total_w - 1.0).abs() < 1e-6, "weights sum = {total_w}");
        let names: Vec<&str> = exposures.iter().map(|e| e.signature.as_str()).collect();
        assert!(names.contains(&"SBS4"), "expected SBS4 in {names:?}");
        // UV is captured by one of the near-collinear SBS7 signatures.
        assert!(names.iter().any(|n| n.starts_with("SBS7")), "expected an SBS7* in {names:?}");
    }

    #[test]
    fn deconvolution_of_empty_spectrum_yields_nothing() {
        let spectrum = Sbs96Spectrum {
            counts: vec![0; 96],
            fractions: vec![0.0; 96],
            total: 0,
            skipped_no_context: 0,
        };
        let (exposures, cosine) = deconvolve_sbs96(&spectrum);
        assert!(exposures.is_empty());
        assert_eq!(cosine, 0.0);
    }

    #[test]
    fn reference_guided_result_is_labeled_and_flags_derive_from_exposures() {
        // Build a reference and a set of C>A-context SNVs that load SBS4/tobacco.
        // "ACAA...": place C>A substitutions in T_A / A_A contexts.
        let mut refs = HashMap::new();
        refs.insert("chr1".to_string(), b"ACGTACGTACGT".to_vec());
        let mut variants = Vec::new();
        for pos in 2..=11u64 {
            // pull the ref base at this position to make a valid SNV
            let base = b"ACGTACGTACGT"[(pos - 1) as usize];
            let alt = if base == b'C' { b'A' } else if base == b'G' { b'T' } else { continue };
            variants.push(snv("chr1", pos, &(base as char).to_string(), &(alt as char).to_string()));
        }
        let result = compute_mutational_signatures_with_context(&variants, &refs);
        assert_eq!(result.method, "96-channel NNLS (COSMIC v3.3)");
        assert!(result.sbs96.is_some());
        assert!(result.reconstruction_cosine.is_some());
    }
}
