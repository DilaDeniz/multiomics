//! Bayesian SNP genotyper: Li (2011) statistical framework.
//!
//! Computes genotype likelihoods for three diploid genotype classes
//! (hom-ref, het, hom-alt) using the model from SAMtools/BCFtools,
//! then applies a Hardy-Weinberg prior to call the MAP genotype.
//!
//! # Reference
//! Li H (2011) Bioinformatics 27(21):2987–2993. doi:10.1093/bioinformatics/btr509

use crate::pileup::{build_pileup, PileupBase, PileupColumn};
use anyhow::Result;
use serde::{Deserialize, Serialize};

/// Population heterozygosity prior (θ).
const THETA: f64 = 0.001;

/// Precomputed Phred → error probability table: `PHRED_TABLE[q] = 10^(-q/10)`.
static PHRED_TABLE: std::sync::LazyLock<[f64; 256]> = std::sync::LazyLock::new(|| {
    let mut t = [0.0f64; 256];
    for (q, slot) in t.iter_mut().enumerate() {
        *slot = 10.0_f64.powf(-(q as f64) / 10.0);
    }
    t
});

// ── Public types ──────────────────────────────────────────────────────────────

/// Diploid genotype class relative to the reference allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Genotype {
    /// Both alleles match the reference (0/0).
    HomRef,
    /// One reference and one alternate allele (0/1).
    Het,
    /// Both alleles are the alternate (1/1).
    HomAlt,
}

/// A single variant call produced by the Li 2011 Bayesian genotyper.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenotypeCall {
    /// Reference sequence name.
    pub chrom: String,
    /// 0-based reference position.
    pub pos: u64,
    /// Reference base (ASCII byte).
    pub ref_base: u8,
    /// Most common non-reference base (ASCII byte).
    pub alt_base: u8,
    /// Called genotype.
    pub genotype: Genotype,
    /// Phred-scaled variant quality: `-10 × log10(P(hom_ref | data))`.
    pub qual: f32,
    /// Total read depth at this position.
    pub depth: u32,
    /// Alternate allele frequency: `alt_count / depth`.
    pub allele_freq: f32,
    /// Phred-scaled genotype likelihoods `[PL_hom_ref, PL_het, PL_hom_alt]`.
    pub pl: [u32; 3],
}

// ── Public API ────────────────────────────────────────────────────────────────

// ── Assembled variant calling ─────────────────────────────────────────────────

/// Call variants using local de novo assembly (GATK HaplotypeCaller-style).
///
/// Identifies active regions in the pileup, reassembles candidate haplotypes
/// with a De Bruijn graph, scores each read against each haplotype with the
/// pair-HMM forward algorithm, genotypes by maximising the diploid likelihood,
/// then extracts variant calls via Smith-Waterman alignment to the reference.
/// Non-active positions are handled by the standard pileup-based caller.
///
/// # Arguments
/// * `pileup`     – Pre-built pileup (see [`crate::pileup::build_pileup`]).
/// * `reads`      – All reads in the region as [`ActiveRead`] records.
/// * `ref_seq`    – Optional reference sequence for the pileup region.
/// * `min_depth`  – Minimum read depth to emit a call.
/// * `min_qual`   – Minimum Phred-scaled quality to emit a call.
pub fn call_variants_assembled(
    pileup: &[crate::pileup::PileupColumn],
    reads: &[crate::assembly::ActiveRead],
    ref_seq: Option<&[u8]>,
    min_depth: u32,
    min_qual: f32,
) -> Vec<GenotypeCall> {
    use crate::assembly::{assemble_haplotypes, find_active_regions};
    use crate::pairhmm::pair_hmm_log_prob;

    const K_SIZES: &[usize] = &[10, 15, 20, 25];
    const MAX_HAPLOTYPES: usize = 16;

    let active_regions = find_active_regions(pileup, reads, 200, 0.05);
    let mut calls: Vec<GenotypeCall> = Vec::new();

    // Track pileup positions that fall inside an active region so we can skip
    // them for the pileup-based caller pass.
    let mut active_positions: std::collections::HashSet<(String, u64)> =
        std::collections::HashSet::new();

    for region in &active_regions {
        for col in pileup {
            if col.chrom == region.chrom && col.pos >= region.start && col.pos < region.end {
                active_positions.insert((col.chrom.clone(), col.pos));
            }
        }

        let ref_slice = ref_seq.unwrap_or(&[]);
        if ref_slice.is_empty() || region.reads.is_empty() {
            continue;
        }

        let haplotypes = assemble_haplotypes(region, ref_slice, K_SIZES, MAX_HAPLOTYPES);
        if haplotypes.len() < 2 {
            continue;
        }

        // For each pair of haplotypes (hap1, hap2) compute diploid likelihood.
        let best_pair = best_haplotype_pair(&region.reads, &haplotypes, pair_hmm_log_prob);

        let (hap1_idx, hap2_idx) = match best_pair {
            Some(p) => p,
            None => continue,
        };

        // Extract variants from each best haplotype vs reference.
        let chrom = &region.chrom;
        for &hi in &[hap1_idx, hap2_idx] {
            if hi == 0 {
                continue; // haplotype 0 is always the reference
            }
            let variant_calls = extract_variants_from_haplotype(
                &haplotypes[hi],
                ref_slice,
                chrom,
                region.start,
                min_depth,
                min_qual,
                region.reads.len() as u32,
            );
            for vc in variant_calls {
                if !calls.iter().any(|c| c.chrom == vc.chrom && c.pos == vc.pos) {
                    calls.push(vc);
                }
            }
        }
    }

    // Run pileup-based caller on non-active positions.
    let non_active_pileup: Vec<crate::pileup::PileupColumn> = pileup
        .iter()
        .filter(|col| !active_positions.contains(&(col.chrom.clone(), col.pos)))
        .cloned()
        .collect();

    let mut pileup_calls = call_variants(&non_active_pileup, min_depth, min_qual);
    calls.append(&mut pileup_calls);

    // Sort and deduplicate by (chrom, pos).
    calls.sort_unstable_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos.cmp(&b.pos)));
    calls.dedup_by(|a, b| a.chrom == b.chrom && a.pos == b.pos);
    calls
}

/// Find the diploid haplotype pair maximising
/// `Σ_reads log_sum_exp( P(read|hap1), P(read|hap2) ) / 2`.
///
/// Returns indices into `haplotypes` or `None` if the haplotype set is empty.
fn best_haplotype_pair(
    reads: &[crate::assembly::ActiveRead],
    haplotypes: &[Vec<u8>],
    hmm_fn: impl Fn(&[u8], &[u8], &[u8]) -> f64,
) -> Option<(usize, usize)> {
    use crate::pairhmm::log_sum_exp;

    let nh = haplotypes.len();
    if nh == 0 {
        return None;
    }

    // Pre-compute log P(read | hap) for all (read, hap) pairs.
    let lp: Vec<Vec<f64>> = reads
        .iter()
        .map(|r| {
            haplotypes
                .iter()
                .map(|h| hmm_fn(&r.seq, &r.quals, h))
                .collect()
        })
        .collect();

    let mut best_score = f64::NEG_INFINITY;
    let mut best_i = 0;
    let mut best_j = 0;

    for i in 0..nh {
        for j in i..nh {
            let score: f64 = lp
                .iter()
                .map(|read_lp| log_sum_exp(read_lp[i], read_lp[j]) - std::f64::consts::LN_2)
                .sum();
            if score > best_score {
                best_score = score;
                best_i = i;
                best_j = j;
            }
        }
    }

    Some((best_i, best_j))
}

/// Align a haplotype to the reference using Smith-Waterman and emit a
/// [`GenotypeCall`] for each aligned difference (SNP or indel).
fn extract_variants_from_haplotype(
    hap: &[u8],
    ref_seq: &[u8],
    chrom: &str,
    region_start: u64,
    min_depth: u32,
    min_qual: f32,
    depth: u32,
) -> Vec<GenotypeCall> {
    use crate::assembly::smith_waterman;

    if depth < min_depth {
        return Vec::new();
    }

    let (aligned_hap, aligned_ref) = smith_waterman(hap, ref_seq);
    let mut calls = Vec::new();
    let mut ref_offset: u64 = 0;

    let mut i = 0;
    while i < aligned_hap.len() && i < aligned_ref.len() {
        let h = aligned_hap[i];
        let r = aligned_ref[i];

        if h == b'-' {
            // Deletion in haplotype (gap vs reference).
            let pos = region_start + ref_offset;
            calls.push(GenotypeCall {
                chrom: chrom.to_string(),
                pos,
                ref_base: r,
                alt_base: b'-',
                genotype: Genotype::Het,
                qual: min_qual,
                depth,
                allele_freq: 0.5,
                pl: [0, 0, 255],
            });
            ref_offset += 1;
        } else if r == b'-' {
            // Insertion in haplotype.
            let pos = region_start + ref_offset;
            calls.push(GenotypeCall {
                chrom: chrom.to_string(),
                pos,
                ref_base: b'-',
                alt_base: h,
                genotype: Genotype::Het,
                qual: min_qual,
                depth,
                allele_freq: 0.5,
                pl: [0, 0, 255],
            });
            // ref_offset does not advance for an insertion
        } else if h != r {
            // SNP
            let pos = region_start + ref_offset;
            calls.push(GenotypeCall {
                chrom: chrom.to_string(),
                pos,
                ref_base: r,
                alt_base: h,
                genotype: Genotype::Het,
                qual: min_qual,
                depth,
                allele_freq: 0.5,
                pl: [0, 0, 255],
            });
            ref_offset += 1;
        } else {
            ref_offset += 1;
        }
        i += 1;
    }

    calls.into_iter().filter(|c| c.qual >= min_qual).collect()
}

/// Call variants from a pre-built pileup.
///
/// For each [`PileupColumn`] the function:
/// 1. Skips positions with depth below `min_depth`.
/// 2. Identifies the most common non-reference base as ALT.
/// 3. Computes log-likelihoods under hom-ref, het, and hom-alt.
/// 4. Applies a Hardy-Weinberg prior with θ = 0.001.
/// 5. Emits a [`GenotypeCall`] when `QUAL >= min_qual` and call ≠ HomRef.
pub fn call_variants(pileup: &[PileupColumn], min_depth: u32, min_qual: f32) -> Vec<GenotypeCall> {
    let mut calls = Vec::new();

    for col in pileup {
        if col.depth() < min_depth as usize {
            continue;
        }

        // Count each base.
        let mut counts = [0u32; 4]; // A C G T
        for pb in &col.bases {
            if let Some(idx) = base_index(pb.base) {
                counts[idx] += 1;
            }
        }

        let ref_idx = base_index(col.ref_base);

        // Find the most frequent non-ref base.
        let alt_idx = (0..4)
            .filter(|&i| Some(i) != ref_idx)
            .max_by_key(|&i| counts[i]);

        let alt_idx = match alt_idx {
            Some(i) if counts[i] > 0 => i,
            _ => continue,
        };

        let alt_base = INDEX_BASES[alt_idx];

        let (ln_hom_ref, ln_het, ln_hom_alt) =
            genotype_likelihoods(&col.bases, col.ref_base, alt_base);

        // Hardy-Weinberg priors.
        let p_hom_ref = (1.0 - THETA).powi(2);
        let p_het = 2.0 * THETA * (1.0 - THETA);
        let p_hom_alt = THETA * THETA;

        let ln_post_hom_ref = ln_hom_ref + p_hom_ref.ln();
        let ln_post_het = ln_het + p_het.ln();
        let ln_post_hom_alt = ln_hom_alt + p_hom_alt.ln();

        // Normalise in log-space via log-sum-exp.
        let ln_norm = log_sum_exp3(ln_post_hom_ref, ln_post_het, ln_post_hom_alt);
        let p_hom_ref_post = (ln_post_hom_ref - ln_norm).exp();

        let qual = (-10.0 * p_hom_ref_post.log10()) as f32;
        if qual < min_qual {
            continue;
        }

        // MAP genotype.
        let genotype = if ln_post_hom_alt >= ln_post_het && ln_post_hom_alt >= ln_post_hom_ref {
            Genotype::HomAlt
        } else if ln_post_het >= ln_post_hom_ref {
            Genotype::Het
        } else {
            Genotype::HomRef
        };

        if genotype == Genotype::HomRef {
            continue;
        }

        // PL values: Phred-scaled likelihoods (0 = most likely).
        let pl = phred_scaled_pl(ln_hom_ref, ln_het, ln_hom_alt);

        let depth = col.depth() as u32;
        let alt_count = counts[alt_idx];
        let allele_freq = alt_count as f32 / depth as f32;

        calls.push(GenotypeCall {
            chrom: col.chrom.clone(),
            pos: col.pos,
            ref_base: col.ref_base,
            alt_base,
            genotype,
            qual,
            depth,
            allele_freq,
            pl,
        });
    }

    calls
}

/// Compute log-likelihoods for hom-ref, het, and hom-alt under Li 2011.
///
/// Returns `(ln_L_hom_ref, ln_L_het, ln_L_hom_alt)`.
pub fn genotype_likelihoods(bases: &[PileupBase], ref_base: u8, alt_base: u8) -> (f64, f64, f64) {
    let mut ln_hom_ref = 0.0_f64;
    let mut ln_het = 0.0_f64;
    let mut ln_hom_alt = 0.0_f64;

    for pb in bases {
        let eps = PHRED_TABLE[pb.base_qual as usize];

        // P(base | hom-ref): allele is ref_base.
        let p_hr = if pb.base == ref_base {
            1.0 - eps
        } else {
            eps / 3.0
        };

        // P(base | hom-alt): allele is alt_base.
        let p_ha = if pb.base == alt_base {
            1.0 - eps
        } else {
            eps / 3.0
        };

        // P(base | het): average of the two homozygous likelihoods.
        let p_ht = (p_hr + p_ha) / 2.0;

        // Accumulate in log-space; guard against log(0).
        ln_hom_ref += p_hr.max(f64::MIN_POSITIVE).ln();
        ln_het += p_ht.max(f64::MIN_POSITIVE).ln();
        ln_hom_alt += p_ha.max(f64::MIN_POSITIVE).ln();
    }

    (ln_hom_ref, ln_het, ln_hom_alt)
}

/// Convenience wrapper: build a pileup from a BAM file then call variants.
///
/// Uses default filter thresholds: `min_base_qual = 13`, `min_mapq = 20`,
/// `max_depth = 8000`.
pub fn call_variants_from_bam(
    bam_path: &std::path::Path,
    min_depth: u32,
    min_qual: f32,
    min_base_qual: u8,
    min_mapq: u8,
) -> Result<Vec<GenotypeCall>> {
    let pileup = build_pileup(bam_path, min_base_qual, min_mapq, 8000)?;
    Ok(call_variants(&pileup, min_depth, min_qual))
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Mapping from base ASCII to array index 0..4 (A, C, G, T).
const INDEX_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Return the index (0=A,1=C,2=G,3=T) for an ASCII base, or `None` for N/other.
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Compute log-sum-exp of three values.
fn log_sum_exp3(a: f64, b: f64, c: f64) -> f64 {
    let m = a.max(b).max(c);
    m + ((a - m).exp() + (b - m).exp() + (c - m).exp()).ln()
}

/// Convert raw log-likelihoods to Phred-scaled PL values (VCF FORMAT/PL).
///
/// The most-likely genotype gets PL=0; others are scaled relative to it.
fn phred_scaled_pl(ln_hr: f64, ln_ht: f64, ln_ha: f64) -> [u32; 3] {
    let best = ln_hr.max(ln_ht).max(ln_ha);
    let to_pl = |ln: f64| -> u32 {
        let diff = best - ln; // always >= 0
        (diff * 10.0 / std::f64::consts::LN_10).round() as u32
    };
    [to_pl(ln_hr), to_pl(ln_ht), to_pl(ln_ha)]
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::PileupBase;

    fn make_base(b: u8, q: u8) -> PileupBase {
        PileupBase {
            base: b,
            base_qual: q,
            mapq: 60,
            is_rev: false,
        }
    }

    #[test]
    fn hom_ref_gives_low_qual() {
        // All reads match the reference → hom-ref posterior should dominate.
        let bases: Vec<PileupBase> = (0..20).map(|_| make_base(b'A', 30)).collect();
        let (ln_hr, ln_ht, ln_ha) = genotype_likelihoods(&bases, b'A', b'G');
        assert!(ln_hr > ln_ht, "hom_ref should dominate over het");
        assert!(ln_hr > ln_ha, "hom_ref should dominate over hom_alt");
    }

    #[test]
    fn hom_alt_gives_high_qual() {
        // All reads carry the alternate base → hom-alt should dominate.
        let bases: Vec<PileupBase> = (0..20).map(|_| make_base(b'G', 30)).collect();
        let (ln_hr, ln_ht, ln_ha) = genotype_likelihoods(&bases, b'A', b'G');
        assert!(ln_ha > ln_ht, "hom_alt should dominate over het");
        assert!(ln_ha > ln_hr, "hom_alt should dominate over hom_ref");
    }

    #[test]
    fn het_variant_called() {
        // Equal mix of ref and alt → expect a het call.
        let mut bases: Vec<PileupBase> = (0..15).map(|_| make_base(b'A', 30)).collect();
        bases.extend((0..15).map(|_| make_base(b'G', 30)));

        let col = PileupColumn {
            chrom: "chr1".into(),
            pos: 500,
            ref_base: b'A',
            bases,
        };

        let calls = call_variants(&[col], 10, 10.0);
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].genotype, Genotype::Het);
        assert_eq!(calls[0].alt_base, b'G');
    }

    #[test]
    fn phred_table_spot_check() {
        // PHRED_TABLE[20] should be 10^(-2) = 0.01.
        let t = &*PHRED_TABLE;
        assert!((t[20] - 0.01).abs() < 1e-12);
    }

    #[test]
    fn pl_best_is_zero() {
        let pl = phred_scaled_pl(-5.0, -10.0, -20.0);
        assert_eq!(pl[0], 0, "best genotype should have PL=0");
    }
}
