//! Cancer-specific genomic analyses: tumor purity, kataegis, HRD, LOH.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::types::VariantRecord;

// ── Tumor Purity ──────────────────────────────────────────────────────────────

/// Tumor purity estimates derived from VAF distribution and methylation data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TumorPurityResult {
    /// Estimated tumor cell fraction [0.0, 1.0] from VAF distribution.
    /// Formula: purity = 2 × mode_AF for diploid het somatic variants.
    /// None when fewer than 20 heterozygous variants are available.
    pub vaf_purity: Option<f64>,
    /// Estimated tumor fraction from global methylation depletion.
    /// Based on: normal tissue ≈ 70% methylated; tumor shows LINE-1 hypomethylation.
    /// Formula: clamp(1.0 - global_meth_pct / 70.0, 0.0, 1.0).
    /// None when methylation data is absent (global_meth_pct == 0.0).
    pub methylation_purity: Option<f64>,
    /// Consensus estimate: weighted mean (VAF weight=0.7, meth weight=0.3).
    /// If only one source available, uses that source alone.
    pub consensus_purity: Option<f64>,
    /// "LOW" (<40%), "MODERATE" (40–70%), "HIGH" (>70%).
    pub purity_class: String,
    /// Cross-modal discordance flag: true when |vaf_purity - methylation_purity| > 0.25.
    /// Indicates sample heterogeneity or technical artifact.
    pub discordant: bool,
}

/// Estimate tumor purity from allele frequencies and optional methylation data.
///
/// Scientific basis: for diploid heterozygous somatic mutations, the expected VAF
/// is purity/2. The modal VAF peak in a pure tumor sample is therefore ~0.5 for
/// fully clonal variants, and the purity = 2 × mode_AF.
///
/// Reference: Carter et al. 2012 (Nature Biotechnology), ABSOLUTE algorithm.
pub fn estimate_tumor_purity(
    variants: &[VariantRecord],
    global_meth_pct: f64,
) -> TumorPurityResult {
    // Collect AFs in the heterozygous somatic window [0.1, 0.85)
    let afs: Vec<f64> = variants
        .iter()
        .filter_map(|v| v.af.map(|af| af as f64))
        .filter(|&af| af > 0.1 && af < 0.85)
        .collect();

    let vaf_purity = if afs.len() < 20 {
        None
    } else {
        // Bin into 20 bins over [0.1, 0.85)
        const N_BINS: usize = 20;
        const LO: f64 = 0.1;
        const HI: f64 = 0.85;
        const RANGE: f64 = HI - LO;

        let mut bins = [0u64; N_BINS];
        for &af in &afs {
            let bin = ((af - LO) / RANGE * N_BINS as f64) as usize;
            let bin = bin.min(N_BINS - 1);
            bins[bin] += 1;
        }

        // Find the bin with the maximum count
        let (mode_bin, _) = bins
            .iter()
            .enumerate()
            .max_by_key(|(_, &count)| count)
            .unwrap_or((0, &0));

        // Use bin center as mode_af
        let mode_af = LO + (mode_bin as f64 + 0.5) / N_BINS as f64 * RANGE;
        Some((mode_af * 2.0).clamp(0.0, 1.0))
    };

    let methylation_purity = if global_meth_pct > 0.0 {
        Some((1.0 - global_meth_pct / 70.0).clamp(0.0, 1.0))
    } else {
        None
    };

    let consensus_purity = match (vaf_purity, methylation_purity) {
        (Some(v), Some(m)) => Some(v * 0.7 + m * 0.3),
        (Some(v), None) => Some(v),
        (None, Some(m)) => Some(m),
        (None, None) => None,
    };

    let purity_for_class = consensus_purity.or(vaf_purity).or(methylation_purity);
    let purity_class = match purity_for_class {
        Some(p) if p > 0.7 => "HIGH".to_string(),
        Some(p) if p >= 0.4 => "MODERATE".to_string(),
        Some(_) => "LOW".to_string(),
        None => "UNKNOWN".to_string(),
    };

    let discordant = match (vaf_purity, methylation_purity) {
        (Some(v), Some(m)) => (v - m).abs() > 0.25,
        _ => false,
    };

    TumorPurityResult {
        vaf_purity,
        methylation_purity,
        consensus_purity,
        purity_class,
        discordant,
    }
}

// ── Kataegis ──────────────────────────────────────────────────────────────────

/// A hypermutation focus (kataegis locus) identified by clustered somatic mutations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KataegisLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_mutations: usize,
    /// Geometric mean inter-mutation distance (bp). Kataegis threshold: < 1000 bp.
    pub geometric_mean_imd: f64,
    /// Most common substitution type in this focus (e.g. "C>T").
    pub dominant_change: String,
}

/// Detect kataegis loci using the Alexandrov et al. 2013 IMD-based algorithm.
///
/// Kataegis is characterised by clusters of mutations with geometric mean
/// inter-mutation distance < 1000 bp. The APOBEC family of cytidine deaminases
/// is the predominant cause, producing characteristic C>T and C>G changes at
/// TpC context (COSMIC SBS2/SBS13).
///
/// Reference: Alexandrov et al. 2013 (Nature), Roberts et al. 2013 (Nature Genetics).
pub fn detect_kataegis(variants: &[VariantRecord]) -> Vec<KataegisLocus> {
    // Group variants by chromosome
    let mut by_chrom: HashMap<&str, Vec<(u64, &VariantRecord)>> = HashMap::new();
    for v in variants {
        by_chrom.entry(&v.chrom).or_default().push((v.pos, v));
    }

    let mut loci: Vec<KataegisLocus> = Vec::new();

    for (chrom, mut positions) in by_chrom {
        if positions.len() < 6 {
            continue;
        }

        // Sort by position
        positions.sort_by_key(|(pos, _)| *pos);

        // Compute IMD array (inter-mutation distances)
        let n = positions.len();
        let imd: Vec<u64> = (0..n - 1)
            .map(|i| positions[i + 1].0.saturating_sub(positions[i].0))
            .collect();

        // Scan for windows of >= 6 consecutive mutations where geometric mean IMD < 1000
        // Use a sliding window approach
        let mut start_idx = 0;
        while start_idx + 5 < n {
            // Window covers variants [start_idx .. start_idx+win_size]
            // with (win_size - 1) IMDs. We need at least 6 mutations => 5 IMDs.
            // Try to extend the window as far as possible while geomean < 1000.

            // Check if a window of exactly 6 starting at start_idx qualifies
            let window_imd: Vec<f64> = (start_idx..start_idx + 5)
                .map(|i| (imd[i] as f64 + 1.0).ln())
                .collect();
            let geomean_6 = (window_imd.iter().sum::<f64>() / window_imd.len() as f64).exp();

            if geomean_6 >= 1000.0 {
                start_idx += 1;
                continue;
            }

            // Extend the window greedily
            let mut end_idx = start_idx + 5; // inclusive: variants[start_idx..=end_idx]
            // imd indices: start_idx .. end_idx (end_idx - start_idx IMDs)
            let mut sum_ln_imd: f64 = (start_idx..end_idx)
                .map(|i| (imd[i] as f64 + 1.0).ln())
                .sum();
            let mut n_imd = (end_idx - start_idx) as f64;

            while end_idx + 1 < n {
                let new_imd_ln = (imd[end_idx] as f64 + 1.0).ln();
                let new_sum = sum_ln_imd + new_imd_ln;
                let new_n = n_imd + 1.0;
                let new_geomean = (new_sum / new_n).exp();
                if new_geomean < 1000.0 {
                    sum_ln_imd = new_sum;
                    n_imd = new_n;
                    end_idx += 1;
                } else {
                    break;
                }
            }

            let n_mutations = end_idx - start_idx + 1;
            let geometric_mean_imd = (sum_ln_imd / n_imd).exp();

            // Compute dominant substitution change among SNPs in this window
            let mut change_counts: HashMap<String, usize> = HashMap::new();
            for (_, v) in positions.iter().take(end_idx + 1).skip(start_idx) {
                // Only SNPs (single-base substitutions)
                if v.ref_allele.len() == 1 && v.alt_allele.len() == 1 {
                    let change = format!(
                        ">{}",
                        v.alt_allele.chars().next().unwrap_or('N')
                    );
                    let ref_char = v.ref_allele.chars().next().unwrap_or('N');
                    let key = format!("{ref_char}{change}");
                    *change_counts.entry(key).or_insert(0) += 1;
                }
            }
            let dominant_change = change_counts
                .into_iter()
                .max_by_key(|(_, count)| *count)
                .map(|(k, _)| k)
                .unwrap_or_else(|| "N/A".to_string());

            loci.push(KataegisLocus {
                chrom: chrom.to_string(),
                start: positions[start_idx].0,
                end: positions[end_idx].0,
                n_mutations,
                geometric_mean_imd,
                dominant_change,
            });

            // Advance past this window
            start_idx = end_idx + 1;
        }
    }

    // Sort by chrom then start
    loci.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    loci
}

// ── HRD Score ─────────────────────────────────────────────────────────────────

/// Homologous Recombination Deficiency score derived from indel size distribution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HrdScore {
    pub total_indels: u64,
    /// Fraction of deletions with size 1 bp (often NHEJ/replication slip).
    pub del_1bp_frac: f64,
    /// Fraction of deletions with size 2–5 bp.
    pub del_2_5bp_frac: f64,
    /// Fraction of deletions with size 6–50 bp (enriched in HRD tumors).
    pub del_6_50bp_frac: f64,
    /// Fraction of insertions with size > 3 bp.
    pub ins_gt3bp_frac: f64,
    /// Composite HRD-indel score [0.0, 1.0].
    /// Elevated del_6_50bp_frac and ins_gt3bp_frac are HRD markers.
    /// Score = (del_6_50bp_frac * 0.6 + ins_gt3bp_frac * 0.4).
    pub hrd_indel_score: f64,
    /// "HRD-HIGH" (score > 0.25), "HRD-INTERMEDIATE" (0.1–0.25), "HRD-LOW" (< 0.1).
    pub hrd_class: String,
    /// Note when total_indels < 50: "Low indel count — result may be unreliable".
    pub note: Option<String>,
}

/// Compute the HRD-indel score from variant indel size distribution.
///
/// Scientific basis: HRD tumors (BRCA1/2 deficient) show enrichment of deletions
/// at microhomology sequences (6–50 bp deletions, COSMIC signature ID8) and large
/// insertions. Without flanking reference sequence, deletion LENGTH is the best
/// available proxy.
///
/// References: Watkins et al. 2020 (Nature Communications),
/// Chan et al. 2015 (Nature Genetics).
pub fn compute_hrd_score(variants: &[VariantRecord]) -> HrdScore {
    let mut del_1bp: u64 = 0;
    let mut del_2_5bp: u64 = 0;
    let mut del_6_50bp: u64 = 0;
    let mut ins_gt3bp: u64 = 0;
    let mut total_indels: u64 = 0;

    for v in variants {
        let ref_len = v.ref_allele.len();
        let alt_len = v.alt_allele.len();
        if ref_len == alt_len {
            continue; // SNP or MNP — skip
        }
        total_indels += 1;
        if ref_len > alt_len {
            // Deletion
            let size = ref_len - alt_len;
            match size {
                1 => del_1bp += 1,
                2..=5 => del_2_5bp += 1,
                6..=50 => del_6_50bp += 1,
                _ => {} // >50 bp deletions not categorized in primary bins
            }
        } else {
            // Insertion
            let size = alt_len - ref_len;
            if size > 3 {
                ins_gt3bp += 1;
            }
        }
    }

    let (del_1bp_frac, del_2_5bp_frac, del_6_50bp_frac, ins_gt3bp_frac) =
        if total_indels == 0 {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let n = total_indels as f64;
            (
                del_1bp as f64 / n,
                del_2_5bp as f64 / n,
                del_6_50bp as f64 / n,
                ins_gt3bp as f64 / n,
            )
        };

    let hrd_indel_score = del_6_50bp_frac * 0.6 + ins_gt3bp_frac * 0.4;

    let hrd_class = if hrd_indel_score > 0.25 {
        "HRD-HIGH".to_string()
    } else if hrd_indel_score >= 0.1 {
        "HRD-INTERMEDIATE".to_string()
    } else {
        "HRD-LOW".to_string()
    };

    let note = if total_indels < 50 {
        Some("Low indel count — result may be unreliable".to_string())
    } else {
        None
    };

    HrdScore {
        total_indels,
        del_1bp_frac,
        del_2_5bp_frac,
        del_6_50bp_frac,
        ins_gt3bp_frac,
        hrd_indel_score,
        hrd_class,
        note,
    }
}

// ── LOH ───────────────────────────────────────────────────────────────────────

/// Per-chromosome loss of heterozygosity assessment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LohChromosome {
    pub chrom: String,
    /// Number of heterozygous variants used (AF in 0.15–0.85).
    pub n_het_variants: usize,
    /// Median |AF - 0.5| for heterozygous variants on this chromosome.
    /// Expected ~0.05 for balanced diploid; elevated (>0.15) suggests allelic imbalance.
    pub median_af_deviation: f64,
    /// Fraction of het variants with AF < 0.35 or AF > 0.65.
    pub skewed_fraction: f64,
    /// True when median_af_deviation > 0.15 AND n_het_variants >= 10.
    pub loh_flagged: bool,
}

/// Detect loss of heterozygosity per chromosome from allele frequency skewing.
///
/// In normal diploid tissue, heterozygous SNP allele frequencies cluster tightly
/// around 0.5. LOH causes one allele to be lost or amplified, shifting AF toward
/// 0.0 or 1.0. This method uses the median absolute deviation of AF from 0.5 as
/// a chromosome-level LOH indicator.
///
/// Reference: Van Loo et al. 2010 (PNAS), ASCAT algorithm.
pub fn detect_loh(variants: &[VariantRecord]) -> Vec<LohChromosome> {
    // Group variants by chromosome, keep only het variants (AF in [0.15, 0.85])
    let mut by_chrom: HashMap<String, Vec<f64>> = HashMap::new();
    for v in variants {
        if let Some(af) = v.af {
            let af = af as f64;
            if (0.15..=0.85).contains(&af) {
                by_chrom.entry(v.chrom.clone()).or_default().push(af);
            }
        }
    }

    let mut results: Vec<LohChromosome> = by_chrom
        .into_iter()
        .filter(|(_, afs)| afs.len() >= 10)
        .map(|(chrom, afs)| {
            let n = afs.len();

            // Compute |AF - 0.5| deviations
            let mut deviations: Vec<f64> = afs.iter().map(|&af| (af - 0.5).abs()).collect();
            deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

            // Median of deviations
            let median_af_deviation = if n % 2 == 0 {
                (deviations[n / 2 - 1] + deviations[n / 2]) / 2.0
            } else {
                deviations[n / 2]
            };

            // Fraction of variants with AF < 0.35 or AF > 0.65
            let skewed_count = afs.iter().filter(|&&af| !(0.35..=0.65).contains(&af)).count();
            let skewed_fraction = skewed_count as f64 / n as f64;

            let loh_flagged = median_af_deviation > 0.15;

            LohChromosome {
                chrom,
                n_het_variants: n,
                median_af_deviation,
                skewed_fraction,
                loh_flagged,
            }
        })
        .collect();

    // Sort by skewed_fraction descending (most skewed first)
    results.sort_by(|a, b| {
        b.skewed_fraction
            .partial_cmp(&a.skewed_fraction)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::TiTvClass;

    fn make_variant(chrom: &str, pos: u64, ref_a: &str, alt_a: &str, af: Option<f32>) -> VariantRecord {
        VariantRecord {
            chrom: chrom.to_string(),
            pos,
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
            qual: 50.0,
            titv: TiTvClass::Transition,
            af,
            gene: None,
        }
    }

    #[test]
    fn test_purity_insufficient_variants() {
        let variants: Vec<VariantRecord> = (0..10)
            .map(|i| make_variant("chr1", i * 1000, "A", "G", Some(0.45)))
            .collect();
        let result = estimate_tumor_purity(&variants, 0.0);
        assert!(result.vaf_purity.is_none());
        assert!(result.methylation_purity.is_none());
    }

    #[test]
    fn test_purity_with_methylation() {
        let result = estimate_tumor_purity(&[], 35.0);
        // meth_purity = clamp(1.0 - 35/70, 0, 1) = 0.5
        assert!((result.methylation_purity.unwrap() - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_hrd_score_empty() {
        let result = compute_hrd_score(&[]);
        assert_eq!(result.total_indels, 0);
        assert_eq!(result.hrd_class, "HRD-LOW");
    }

    #[test]
    fn test_hrd_score_del_classification() {
        let mut variants = Vec::new();
        // 1 bp deletion
        variants.push(make_variant("chr1", 100, "AG", "A", None));
        // 3 bp deletion
        variants.push(make_variant("chr1", 200, "ACGT", "A", None));
        // 8 bp deletion (HRD marker)
        variants.push(make_variant("chr1", 300, "ACGTACGTA", "A", None));
        let result = compute_hrd_score(&variants);
        assert_eq!(result.total_indels, 3);
        assert!((result.del_1bp_frac - 1.0 / 3.0).abs() < 1e-9);
        assert!((result.del_6_50bp_frac - 1.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_detect_kataegis_empty() {
        let loci = detect_kataegis(&[]);
        assert!(loci.is_empty());
    }

    #[test]
    fn test_detect_loh_empty() {
        let loci = detect_loh(&[]);
        assert!(loci.is_empty());
    }
}
