//! Cancer-specific genomic analyses: tumor purity, kataegis, HRD, LOH, TMB, MSI.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::types::VariantRecord;

// ── Tumor Mutational Burden (TMB) ─────────────────────────────────────────────

/// TMB result: total somatic mutations normalised to megabases of sequenced genome.
///
/// FDA approved pembrolizumab for TMB-H solid tumors (≥10 mut/Mb) in 2020.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TmbResult {
    /// Total variant count used for TMB calculation.
    pub total_variants: u64,
    /// Effective genome size in megabases used for normalization.
    pub genome_mb: f64,
    /// TMB = total_variants / genome_mb.
    pub tmb: f64,
    /// "TMB-H" (≥10), "TMB-L" (1–<10), "TMB-ZERO" (<1).
    pub tmb_class: String,
    /// Source of genome_mb: "WGS (auto)", "WES (auto)", "user-specified".
    pub genome_mb_source: String,
}

/// Compute TMB from a variant count and an effective genome size.
///
/// Reference: Chalmers et al. 2017 (Genome Medicine); FDA approval 2020.
pub fn compute_tmb(total_variants: u64, genome_mb: f64, genome_mb_source: &str) -> TmbResult {
    let tmb = total_variants as f64 / genome_mb;
    let tmb_class = if tmb >= 10.0 {
        "TMB-H".to_string()
    } else if tmb >= 1.0 {
        "TMB-L".to_string()
    } else {
        "TMB-ZERO".to_string()
    };
    TmbResult {
        total_variants,
        genome_mb,
        tmb,
        tmb_class,
        genome_mb_source: genome_mb_source.to_string(),
    }
}

// ── Microsatellite Instability (MSI) ─────────────────────────────────────────

/// MSI result derived from the homopolymer indel fraction and short-indel burden.
///
/// FDA approved pembrolizumab for MSI-H tumors regardless of histology (2017).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MsiResult {
    pub total_variants: u64,
    pub total_indels: u64,
    /// Number of indels classified as homopolymer-context (≥3 same bases deleted/inserted).
    pub homopolymer_indels: u64,
    /// Fraction of indels in homopolymer context.
    pub homopolymer_indel_frac: f64,
    /// Fraction of ALL variants that are short (1–4 bp) indels.
    pub short_indel_frac: f64,
    /// Composite MSI score: homopolymer_frac × 0.6 + short_indel_frac × 0.4.
    pub msi_score: f64,
    /// "MSI-H" (>0.30), "MSI-L" (0.10–0.30), "MSS" (<0.10).
    pub msi_class: String,
    /// Note when total_indels < 30: result may be unreliable.
    pub note: Option<String>,
}

/// Returns true if the indel (ref, alt) involves a homopolymer run of ≥3 identical bases.
///
/// For deletions: checks whether the deleted sequence is a single-base repeat.
/// For insertions: checks whether the inserted sequence is a single-base repeat.
fn is_homopolymer_indel(ref_allele: &str, alt_allele: &str) -> bool {
    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();
    if ref_len == alt_len {
        return false;
    }
    let seq = if ref_len > alt_len {
        // Deletion: deleted sequence = ref_allele[alt_len..]
        &ref_allele[alt_len..]
    } else {
        // Insertion: inserted sequence = alt_allele[ref_len..]
        &alt_allele[ref_len..]
    };
    if seq.len() < 3 {
        return false;
    }
    let first = seq.as_bytes()[0];
    seq.bytes().all(|b| b == first)
}

/// Compute MSI score from the full variant list.
///
/// References: Bonneville et al. 2017 (JCO Precision Oncology),
/// Cortes-Ciriano et al. 2017 (Nature Communications).
pub fn compute_msi(variants: &[VariantRecord]) -> MsiResult {
    let total_variants = variants.len() as u64;
    let mut total_indels: u64 = 0;
    let mut homopolymer_indels: u64 = 0;
    let mut short_indels: u64 = 0; // 1–4 bp indels

    for v in variants {
        let ref_len = v.ref_allele.len();
        let alt_len = v.alt_allele.len();
        if ref_len == alt_len {
            continue; // SNP or MNP — not an indel
        }
        total_indels += 1;
        let indel_size = ref_len.abs_diff(alt_len);
        if indel_size >= 1 && indel_size <= 4 {
            short_indels += 1;
        }
        if is_homopolymer_indel(&v.ref_allele, &v.alt_allele) {
            homopolymer_indels += 1;
        }
    }

    let homopolymer_indel_frac = if total_indels == 0 {
        0.0
    } else {
        homopolymer_indels as f64 / total_indels as f64
    };

    let short_indel_frac = if total_variants == 0 {
        0.0
    } else {
        short_indels as f64 / total_variants as f64
    };

    let msi_score = homopolymer_indel_frac * 0.6 + short_indel_frac * 0.4;

    let msi_class = if msi_score > 0.30 {
        "MSI-H".to_string()
    } else if msi_score >= 0.10 {
        "MSI-L".to_string()
    } else {
        "MSS".to_string()
    };

    let note = if total_indels < 30 {
        Some("Low indel count — MSI result may be unreliable".to_string())
    } else {
        None
    };

    MsiResult {
        total_variants,
        total_indels,
        homopolymer_indels,
        homopolymer_indel_frac,
        short_indel_frac,
        msi_score,
        msi_class,
        note,
    }
}

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
    /// Fraction of deletions with microhomology >= 1 bp at the deletion site.
    /// Only populated when a reference FASTA is provided.
    #[serde(default)]
    pub del_with_mh_frac: f64,
    /// Composite HRD-indel score [0.0, 1.0].
    /// Without reference: Score = (del_6_50bp_frac * 0.6 + ins_gt3bp_frac * 0.4).
    /// With reference:    Score = (del_with_mh_frac * 0.7 + del_6_50bp_frac * 0.3).
    pub hrd_indel_score: f64,
    /// "HRD-HIGH" (score > 0.25), "HRD-INTERMEDIATE" (0.1–0.25), "HRD-LOW" (< 0.1).
    pub hrd_class: String,
    /// Note when total_indels < 50: "Low indel count — result may be unreliable".
    pub note: Option<String>,
    /// True when microhomology-based scoring was used (reference FASTA was provided).
    #[serde(default)]
    pub reference_used: bool,
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
        del_with_mh_frac: 0.0,
        hrd_indel_score,
        hrd_class,
        note,
        reference_used: false,
    }
}

// ── Microhomology-based HRD ───────────────────────────────────────────────────

/// Compute microhomology length at a deletion site.
///
/// `deleted_seq`: the deleted bases (ref bases beyond the first shared base).
/// `left_flank`: reference bases immediately left of the deletion.
/// `right_flank`: reference bases immediately right of the deletion end.
///
/// Returns the maximum of left-flank and right-flank microhomology lengths.
fn microhomology_length(
    deleted_seq: &[u8],
    left_flank: &[u8],
    right_flank: &[u8],
) -> usize {
    let del_len = deleted_seq.len();
    // Left MH: compare end of left_flank with deleted_seq
    let left_mh = (1..=del_len.min(left_flank.len()))
        .rev()
        .find(|&k| left_flank[left_flank.len() - k..] == deleted_seq[..k])
        .unwrap_or(0);
    // Right MH: compare start of right_flank with end of deleted_seq
    let right_mh = (1..=del_len.min(right_flank.len()))
        .rev()
        .find(|&k| right_flank[..k] == deleted_seq[del_len - k..])
        .unwrap_or(0);
    left_mh.max(right_mh)
}

/// Public wrapper around [`parse_fasta_selective`] for reuse by other modules
/// (e.g. reference-guided COSMIC signature analysis).
pub fn parse_fasta_selective_pub(
    data: &[u8],
    needed_chroms: &std::collections::HashSet<&str>,
) -> HashMap<String, Vec<u8>> {
    parse_fasta_selective(data, needed_chroms)
}

/// Parse a FASTA file into a map of chromosome name → uppercase sequence bytes.
/// Only chromosomes present in `needed_chroms` are retained to save memory.
fn parse_fasta_selective(
    data: &[u8],
    needed_chroms: &std::collections::HashSet<&str>,
) -> HashMap<String, Vec<u8>> {
    let mut map: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in data.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            // Flush previous
            if let Some(name) = current_name.take() {
                if needed_chroms.contains(name.as_str()) {
                    map.insert(name, std::mem::take(&mut current_seq));
                } else {
                    current_seq.clear();
                }
            }
            // Parse header: take first whitespace-delimited token after '>'
            let header = &line[1..];
            let name_bytes = header
                .iter()
                .position(|&b| b == b' ' || b == b'\t' || b == b'\r')
                .map_or(header, |n| &header[..n]);
            let name = String::from_utf8_lossy(name_bytes).into_owned();
            current_name = Some(name);
        } else {
            // Strip carriage returns and accumulate as uppercase
            let trimmed = if line.last() == Some(&b'\r') {
                &line[..line.len() - 1]
            } else {
                line
            };
            current_seq.extend(trimmed.iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    // Flush last record
    if let Some(name) = current_name {
        if needed_chroms.contains(name.as_str()) {
            map.insert(name, current_seq);
        }
    }
    map
}

/// Compute the HRD-indel score using reference-guided microhomology detection.
///
/// Requires a reference FASTA file. For each deletion variant, looks up flanking
/// sequence and computes microhomology length. Deletions with MH >= 1 bp are
/// characteristic of HR deficiency (COSMIC signature ID8).
///
/// The enhanced score formula:
///   `hrd_indel_score = del_with_mh_frac * 0.7 + del_6_50bp_frac * 0.3`
///
/// References: Watkins et al. 2020 (Nature Communications),
/// Chan et al. 2015 (Nature Genetics).
pub fn compute_hrd_score_with_reference(
    variants: &[VariantRecord],
    reference_path: &std::path::Path,
) -> anyhow::Result<HrdScore> {
    use memmap2::Mmap;

    let file = std::fs::File::open(reference_path)
        .map_err(|e| anyhow::anyhow!("Cannot open reference FASTA '{}': {e}", reference_path.display()))?;

    // SAFETY: file is not modified while this process holds the Mmap.
    let mmap = unsafe { Mmap::map(&file) }
        .map_err(|e| anyhow::anyhow!("Cannot mmap reference FASTA: {e}"))?;

    let _ = mmap.advise(memmap2::Advice::Sequential);

    // Collect chromosome names needed
    let needed_chroms: std::collections::HashSet<&str> = variants
        .iter()
        .map(|v| v.chrom.as_str())
        .collect();

    let ref_seqs = parse_fasta_selective(mmap.as_ref(), &needed_chroms);

    let mut del_1bp: u64 = 0;
    let mut del_2_5bp: u64 = 0;
    let mut del_6_50bp: u64 = 0;
    let mut ins_gt3bp: u64 = 0;
    let mut del_with_mh: u64 = 0;
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
                _ => {}
            }
            // Microhomology detection
            // pos is 1-based in VCF; the ref allele starts at pos.
            // The deleted bases are ref_allele[alt_len..] (after the shared anchor base(s)).
            if let Some(seq) = ref_seqs.get(&v.chrom) {
                // VCF convention: REF includes anchor base at pos (1-based).
                // deleted bases = ref_allele[alt_len..] (the actual removed sequence)
                let deleted_bytes = &v.ref_allele.as_bytes()[alt_len..];
                let del_len = deleted_bytes.len();

                // del_start_0 in 0-based coords: first deleted base is at (v.pos - 1 + alt_len)
                let anchor_offset = alt_len; // number of shared bases (usually 1)
                let del_start_0 = v.pos as usize - 1 + anchor_offset;

                if del_start_0 < seq.len() {
                    // Left flank: up to del_len bases immediately left of deletion
                    let left_start = del_start_0.saturating_sub(del_len);
                    let left_flank = &seq[left_start..del_start_0];

                    // Right flank: up to del_len bases immediately right of deletion end
                    let right_end_0 = del_start_0 + del_len;
                    let right_end = (right_end_0 + del_len).min(seq.len());
                    let right_flank = if right_end_0 < seq.len() {
                        &seq[right_end_0..right_end]
                    } else {
                        &seq[0..0]
                    };

                    let mh = microhomology_length(deleted_bytes, left_flank, right_flank);
                    if mh >= 1 {
                        del_with_mh += 1;
                    }
                }
            }
        } else {
            // Insertion
            let size = alt_len - ref_len;
            if size > 3 {
                ins_gt3bp += 1;
            }
        }
    }

    let (del_1bp_frac, del_2_5bp_frac, del_6_50bp_frac, ins_gt3bp_frac, del_with_mh_frac) =
        if total_indels == 0 {
            (0.0, 0.0, 0.0, 0.0, 0.0)
        } else {
            let n = total_indels as f64;
            (
                del_1bp as f64 / n,
                del_2_5bp as f64 / n,
                del_6_50bp as f64 / n,
                ins_gt3bp as f64 / n,
                del_with_mh as f64 / n,
            )
        };

    // Enhanced score: microhomology fraction weighted more heavily
    let hrd_indel_score = del_with_mh_frac * 0.7 + del_6_50bp_frac * 0.3;

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

    Ok(HrdScore {
        total_indels,
        del_1bp_frac,
        del_2_5bp_frac,
        del_6_50bp_frac,
        ins_gt3bp_frac,
        del_with_mh_frac,
        hrd_indel_score,
        hrd_class,
        note,
        reference_used: true,
    })
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
