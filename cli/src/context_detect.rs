//! Automatic sample context detection — Feature #3.
//!
//! Scans the first N records of each input file to infer biological context
//! (species, assay type, mutation signatures, suggested preset) without
//! loading the entire file into memory.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use biomics_core::parse::{ByteLines, TabFields};
use serde::{Deserialize, Serialize};

/// Maximum non-header VCF / BED lines to scan.
const MAX_GENOMIC_LINES: usize = 10_000;

/// Maximum TSV lines to scan (header + data).
const MAX_TSV_LINES: usize = 100;

// ── Public types ──────────────────────────────────────────────────────────────

/// Detected biological species from file characteristics.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum DetectedSpecies {
    Human,
    Mouse,
    Plant, // e.g. Arabidopsis, rice
    Unknown,
}

impl DetectedSpecies {
    pub fn as_str(&self) -> &'static str {
        match self {
            DetectedSpecies::Human => "human",
            DetectedSpecies::Mouse => "mouse",
            DetectedSpecies::Plant => "plant",
            DetectedSpecies::Unknown => "unknown",
        }
    }
}

/// Detected sequencing assay type from file characteristics.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum DetectedAssay {
    WholeGenomeSeq,
    WholeExomeSeq,
    WholeGenomeBisulfite,
    ReducedRepresentationBisulfite,
    RnaSeq,
    Unknown,
}

impl DetectedAssay {
    pub fn as_str(&self) -> &'static str {
        match self {
            DetectedAssay::WholeGenomeSeq => "WGS",
            DetectedAssay::WholeExomeSeq => "WES",
            DetectedAssay::WholeGenomeBisulfite => "WGBS",
            DetectedAssay::ReducedRepresentationBisulfite => "RRBS",
            DetectedAssay::RnaSeq => "RNA-seq",
            DetectedAssay::Unknown => "unknown",
        }
    }
}

/// Suggested preset with a confidence score.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PresetSuggestion {
    /// e.g. "cancer", "plant", "wgbs"
    pub preset: String,
    /// 0.0–1.0
    pub confidence: f64,
    /// Human-readable evidence.
    pub reasons: Vec<String>,
}

/// Cross-modality concordance check result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConcordanceCheck {
    /// Whether the VCF and BED appear to be from the same species.
    pub species_concordant: bool,
    /// Whether Ti/Tv in VCF is compatible with expression characteristics.
    pub titv_expression_concordant: bool,
    /// Any warnings generated.
    pub warnings: Vec<String>,
}

/// Full sample context inferred from file scanning.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleContext {
    pub species: DetectedSpecies,
    pub genomics_assay: Option<DetectedAssay>,
    pub epigenomics_assay: Option<DetectedAssay>,
    /// Estimated somatic mutation burden (variants/Mb); None if insufficient data.
    pub somatic_burden_per_mb: Option<f64>,
    /// Dominant mutation signature hint (e.g. "C>T predominant (UV/aging)").
    pub mutation_signature_hint: Option<String>,
    /// Ti/Tv ratio from the first N variants.
    pub titv_ratio: Option<f64>,
    /// Suggested preset and confidence.
    pub suggested_preset: Option<PresetSuggestion>,
    /// Cross-modality concordance.
    pub concordance: ConcordanceCheck,
    /// All warnings generated during detection.
    pub warnings: Vec<String>,
}

// ── Internal scan results ─────────────────────────────────────────────────────

struct VcfScan {
    species: DetectedSpecies,
    assay: DetectedAssay,
    titv_ratio: f64,
    somatic_burden_per_mb: f64,
    mutation_signature_hint: Option<String>,
    warnings: Vec<String>,
}

struct BedScan {
    species: DetectedSpecies,
    assay: DetectedAssay,
    warnings: Vec<String>,
}

struct TsvScan {
    n_samples: usize,
}

// ── Entry point ───────────────────────────────────────────────────────────────

/// Scan input files and return a best-effort `SampleContext`.
///
/// Never fails — all IO errors are absorbed into `SampleContext::warnings`.
pub fn detect_sample_context(
    vcf_path: Option<&Path>,
    tsv_path: Option<&Path>,
    bed_path: Option<&Path>,
    atac_path: Option<&Path>,
) -> SampleContext {
    let mut global_warnings: Vec<String> = Vec::new();

    let vcf_scan = vcf_path.and_then(|p| {
        match scan_vcf(p) {
            Ok(s) => Some(s),
            Err(e) => {
                global_warnings.push(format!("VCF scan failed ({}): {}", p.display(), e));
                None
            }
        }
    });

    let bed_scan = bed_path.and_then(|p| {
        match scan_bed(p) {
            Ok(s) => Some(s),
            Err(e) => {
                global_warnings.push(format!("BED scan failed ({}): {}", p.display(), e));
                None
            }
        }
    });

    let tsv_scan = tsv_path.and_then(|p| {
        match scan_tsv(p) {
            Ok(s) => Some(s),
            Err(e) => {
                global_warnings.push(format!("TSV scan failed ({}): {}", p.display(), e));
                None
            }
        }
    });

    // Collect all modality warnings
    if let Some(ref v) = vcf_scan {
        global_warnings.extend(v.warnings.iter().cloned());
    }
    if let Some(ref b) = bed_scan {
        global_warnings.extend(b.warnings.iter().cloned());
    }

    // Determine species: VCF wins, then BED, then Unknown
    let species = vcf_scan
        .as_ref()
        .map(|v| v.species.clone())
        .or_else(|| bed_scan.as_ref().map(|b| b.species.clone()))
        .unwrap_or(DetectedSpecies::Unknown);

    let genomics_assay = vcf_scan.as_ref().map(|v| v.assay.clone());
    let epigenomics_assay = bed_scan.as_ref().map(|b| b.assay.clone());
    let titv_ratio = vcf_scan.as_ref().map(|v| v.titv_ratio);
    let somatic_burden_per_mb = vcf_scan.as_ref().map(|v| v.somatic_burden_per_mb);
    let mutation_signature_hint =
        vcf_scan.as_ref().and_then(|v| v.mutation_signature_hint.clone());

    // Concordance check
    let concordance = build_concordance(
        vcf_scan.as_ref(),
        bed_scan.as_ref(),
        titv_ratio,
        &mut global_warnings,
    );

    // Preset suggestion
    let suggested_preset = suggest_preset(
        &species,
        genomics_assay.as_ref(),
        epigenomics_assay.as_ref(),
        somatic_burden_per_mb,
        titv_ratio,
        tsv_scan.as_ref(),
        vcf_path,
        bed_path,
        atac_path,
    );

    // Warn about single-sample TSV
    if let Some(ref t) = tsv_scan {
        if t.n_samples == 1 {
            global_warnings
                .push("Single-sample expression detected — DE analysis requires ≥2 groups".into());
        }
    }

    SampleContext {
        species,
        genomics_assay,
        epigenomics_assay,
        somatic_burden_per_mb,
        mutation_signature_hint,
        titv_ratio,
        suggested_preset,
        concordance,
        warnings: global_warnings,
    }
}

// ── VCF scanning ──────────────────────────────────────────────────────────────

fn scan_vcf(path: &Path) -> Result<VcfScan, String> {
    let mmap_result = try_mmap(path);
    let data_owned: Vec<u8>;
    let data: &[u8] = match mmap_result {
        Ok(ref m) => m,
        Err(_) => {
            data_owned = std::fs::read(path).map_err(|e| e.to_string())?;
            &data_owned
        }
    };

    let mut chrom_set: ahash::AHashSet<String> = ahash::AHashSet::new();
    let mut per_chrom: HashMap<String, u64> = HashMap::new();
    let mut transitions: u64 = 0;
    let mut transversions: u64 = 0;
    let mut ct_transitions: u64 = 0; // C>T (or G>A)
    let mut cg_transversions: u64 = 0; // C>G (or G>C)
    let mut ta_transversions: u64 = 0; // T>A (or A>T)
    let mut total_variants: u64 = 0;
    let mut warnings: Vec<String> = Vec::new();

    let mut data_lines = 0usize;
    for line in ByteLines::new(data) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        if data_lines >= MAX_GENOMIC_LINES {
            break;
        }
        data_lines += 1;

        // Parse tab-delimited fields: CHROM POS ID REF ALT ...
        let mut fields = TabFields::new(line);
        let chrom = match fields.next() {
            Some(c) => c,
            None => continue,
        };
        // skip POS and ID
        fields.next(); // POS
        fields.next(); // ID
        let ref_allele = match fields.next() {
            Some(r) => r,
            None => continue,
        };
        let alt_allele = match fields.next() {
            Some(a) => a,
            None => continue,
        };

        // Ignore multi-alt / indels (only single-char ref and alt)
        if ref_allele.len() == 1 && alt_allele.len() == 1 {
            total_variants += 1;
            let r = ref_allele[0].to_ascii_uppercase();
            let a = alt_allele[0].to_ascii_uppercase();

            if let Some(titv) = classify_titv(r, a) {
                if titv {
                    transitions += 1;
                    // C>T signature
                    if (r == b'C' && a == b'T') || (r == b'G' && a == b'A') {
                        ct_transitions += 1;
                    }
                } else {
                    transversions += 1;
                    // C>G signature
                    if (r == b'C' && a == b'G') || (r == b'G' && a == b'C') {
                        cg_transversions += 1;
                    }
                    // T>A signature
                    if (r == b'T' && a == b'A') || (r == b'A' && a == b'T') {
                        ta_transversions += 1;
                    }
                }
            }
        }

        let chrom_str = match std::str::from_utf8(chrom) {
            Ok(s) => s.to_string(),
            Err(_) => continue,
        };
        chrom_set.insert(chrom_str.clone());
        *per_chrom.entry(chrom_str).or_insert(0) += 1;
    }

    let titv_ratio = transitions as f64 / transversions.max(1) as f64;

    // Somatic burden estimate: total SNVs / ~3000 Mb (human genome size)
    let somatic_burden_per_mb = total_variants as f64 / 3000.0;

    // Mutation signature hints
    let mutation_signature_hint = {
        let mut hint: Option<String> = None;
        if transitions > 0 {
            let ct_frac = ct_transitions as f64 / transitions as f64;
            if ct_frac > 0.60 {
                hint = Some("C>T predominant (aging/UV signature SBS1/SBS7)".into());
            }
        }
        if hint.is_none() && transversions > 0 {
            let cg_frac = cg_transversions as f64 / transversions as f64;
            if cg_frac > 0.20 {
                hint = Some("C>G predominant (tobacco/SBS4 signature)".into());
            }
        }
        if hint.is_none() && transversions > 0 {
            let ta_frac = ta_transversions as f64 / transversions as f64;
            if ta_frac > 0.15 {
                hint = Some("T>A predominant (SBS22/aflatoxin)".into());
            }
        }
        hint
    };

    // Species detection
    let species = detect_species_from_chroms(&chrom_set);

    // WGS vs WES via CV of per-chrom counts
    let assay = if per_chrom.len() >= 5 {
        let cv = coefficient_of_variation(&per_chrom);
        if cv > 1.5 {
            warnings.push(format!(
                "VCF: high per-chrom CV ({:.2}) suggests WES capture rather than WGS",
                cv
            ));
            DetectedAssay::WholeExomeSeq
        } else {
            DetectedAssay::WholeGenomeSeq
        }
    } else {
        warnings.push(format!(
            "VCF: only {} chromosomes in scanned lines — cannot distinguish WGS/WES",
            per_chrom.len()
        ));
        DetectedAssay::Unknown
    };

    // Ti/Tv range warnings
    if titv_ratio < 1.5 {
        warnings.push(format!(
            "Ti/Tv ratio {:.2} is unusually low (expected ≥1.8 for human WGS)",
            titv_ratio
        ));
    } else if titv_ratio > 3.0 {
        warnings.push(format!(
            "Ti/Tv ratio {:.2} is unusually high (expected ≤3.0)",
            titv_ratio
        ));
    }

    Ok(VcfScan {
        species,
        assay,
        titv_ratio,
        somatic_burden_per_mb,
        mutation_signature_hint,
        warnings,
    })
}

/// Returns Some(true)=transition, Some(false)=transversion, None=skip.
#[inline]
fn classify_titv(r: u8, a: u8) -> Option<bool> {
    match (r, a) {
        // Transitions: purines or pyrimidines exchanged
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => Some(true),
        // Transversions: purine <-> pyrimidine
        (b'A', b'C')
        | (b'C', b'A')
        | (b'A', b'T')
        | (b'T', b'A')
        | (b'G', b'C')
        | (b'C', b'G')
        | (b'G', b'T')
        | (b'T', b'G') => Some(false),
        _ => None,
    }
}

// ── BED scanning ──────────────────────────────────────────────────────────────

fn scan_bed(path: &Path) -> Result<BedScan, String> {
    let mmap_result = try_mmap(path);
    let data_owned: Vec<u8>;
    let data: &[u8] = match mmap_result {
        Ok(ref m) => m,
        Err(_) => {
            data_owned = std::fs::read(path).map_err(|e| e.to_string())?;
            &data_owned
        }
    };

    let mut chrom_set: ahash::AHashSet<String> = ahash::AHashSet::new();
    let mut per_chrom: HashMap<String, u64> = HashMap::new();
    let mut total_sites: u64 = 0;
    let mut warnings: Vec<String> = Vec::new();

    let mut data_lines = 0usize;
    for line in ByteLines::new(data) {
        if line.is_empty() || line[0] == b'#' || line.starts_with(b"track") || line.starts_with(b"browser") {
            continue;
        }
        if data_lines >= MAX_GENOMIC_LINES {
            break;
        }
        data_lines += 1;

        let mut fields = TabFields::new(line);
        let chrom = match fields.next() {
            Some(c) => c,
            None => continue,
        };

        let chrom_str = match std::str::from_utf8(chrom) {
            Ok(s) => s.to_string(),
            Err(_) => continue,
        };
        chrom_set.insert(chrom_str.clone());
        *per_chrom.entry(chrom_str).or_insert(0) += 1;
        total_sites += 1;
    }

    let species = detect_species_from_chroms(&chrom_set);

    // WGBS vs RRBS: use CV and total site count
    let assay = if per_chrom.len() >= 5 {
        let cv = coefficient_of_variation(&per_chrom);
        if total_sites < 2000 || cv > 1.5 {
            warnings.push(format!(
                "BED: RRBS pattern detected (sites={}, CV={:.2})",
                total_sites, cv
            ));
            DetectedAssay::ReducedRepresentationBisulfite
        } else if cv < 0.5 {
            DetectedAssay::WholeGenomeBisulfite
        } else {
            // Ambiguous — lean toward WGBS unless total sites is very low
            if total_sites < 2_000_000 {
                warnings.push(format!(
                    "BED: moderate CV ({:.2}) — possible RRBS; site count={}",
                    cv, total_sites
                ));
                DetectedAssay::ReducedRepresentationBisulfite
            } else {
                DetectedAssay::WholeGenomeBisulfite
            }
        }
    } else {
        warnings.push(format!(
            "BED: only {} chromosomes in scanned lines — cannot distinguish WGBS/RRBS",
            per_chrom.len()
        ));
        DetectedAssay::Unknown
    };

    Ok(BedScan {
        species,
        assay,
        warnings,
    })
}

// ── TSV scanning ──────────────────────────────────────────────────────────────

fn scan_tsv(path: &Path) -> Result<TsvScan, String> {
    let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
    let reader = BufReader::new(file);

    let mut n_samples = 0usize;
    for (i, line) in reader.lines().enumerate() {
        if i >= MAX_TSV_LINES {
            break;
        }
        let line = line.map_err(|e| e.to_string())?;
        if i == 0 {
            // Header row: count columns (all after gene/first column = samples)
            let ncols = line.split('\t').count();
            n_samples = ncols.saturating_sub(1);
        }
    }

    Ok(TsvScan { n_samples })
}

// ── Species detection ─────────────────────────────────────────────────────────

fn detect_species_from_chroms(chroms: &ahash::AHashSet<String>) -> DetectedSpecies {
    if chroms.is_empty() {
        return DetectedSpecies::Unknown;
    }

    // Plant: chromosome names with capital "Chr" prefix (Arabidopsis Chr1–Chr5)
    let has_capital_chr = chroms.iter().any(|c| c.starts_with("Chr"));
    if has_capital_chr {
        return DetectedSpecies::Plant;
    }

    // Determine max numeric chromosome number from names like "chr1", "chrX", "1", "22"
    let max_num = chroms.iter().filter_map(|c| {
        let stripped = c
            .trim_start_matches("chr")
            .trim_start_matches("Chr");
        stripped.parse::<u32>().ok()
    }).max();

    let has_lowercase_chr_prefix = chroms.iter().any(|c| c.starts_with("chr") && !c.starts_with("Chr"));
    let has_no_prefix_numeric = chroms.iter().any(|c| c.parse::<u32>().is_ok());

    match max_num {
        Some(n) if n > 19 => {
            // Only human has chromosomes 20, 21, 22
            DetectedSpecies::Human
        }
        Some(n) if n <= 19 && has_lowercase_chr_prefix => {
            // Could be mouse (chr1–chr19) or human (might not have scanned chr20+ yet)
            // Mouse has exactly chr1–chr19 + chrX + chrY; human also has chr1–chr19
            // Conservative: if we see exactly ≤19 numbered chroms with "chr" prefix,
            // treat as ambiguous but lean mouse if no "chr20"+ present
            let _ = n;
            DetectedSpecies::Mouse
        }
        Some(_) if has_no_prefix_numeric => {
            // No prefix, numeric: human standard (1..22, X, Y)
            DetectedSpecies::Human
        }
        _ => DetectedSpecies::Unknown,
    }
}

// ── Statistical helpers ───────────────────────────────────────────────────────

fn coefficient_of_variation(counts: &HashMap<String, u64>) -> f64 {
    if counts.len() < 2 {
        return 0.0;
    }
    let vals: Vec<f64> = counts.values().map(|&v| v as f64).collect();
    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;
    if mean == 0.0 {
        return 0.0;
    }
    let variance = vals.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n;
    variance.sqrt() / mean
}

// ── Memory-map helper ─────────────────────────────────────────────────────────

fn try_mmap(path: &Path) -> Result<memmap2::Mmap, String> {
    let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
    // SAFETY: we treat the mapped memory as read-only bytes
    unsafe { memmap2::Mmap::map(&file).map_err(|e| e.to_string()) }
}

// ── Concordance ───────────────────────────────────────────────────────────────

fn build_concordance(
    vcf: Option<&VcfScan>,
    bed: Option<&BedScan>,
    titv_ratio: Option<f64>,
    warnings: &mut Vec<String>,
) -> ConcordanceCheck {
    let mut concordance_warnings: Vec<String> = Vec::new();

    // Species concordance
    let species_concordant = match (vcf, bed) {
        (Some(v), Some(b)) => {
            let ok = v.species == b.species || v.species == DetectedSpecies::Unknown || b.species == DetectedSpecies::Unknown;
            if !ok {
                let msg = format!(
                    "Species mismatch: VCF suggests {:?} but BED suggests {:?}",
                    v.species, b.species
                );
                concordance_warnings.push(msg.clone());
                warnings.push(msg);
            }
            ok
        }
        _ => true,
    };

    // Ti/Tv expression concordance
    let titv_expression_concordant = match titv_ratio {
        Some(r) if r < 1.0 => {
            let msg = format!(
                "Ti/Tv ratio {:.2} < 1.0 — possible RNA editing artifact or data quality issue",
                r
            );
            concordance_warnings.push(msg.clone());
            warnings.push(msg);
            false
        }
        _ => true,
    };

    ConcordanceCheck {
        species_concordant,
        titv_expression_concordant,
        warnings: concordance_warnings,
    }
}

// ── Preset suggestion ─────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn suggest_preset(
    species: &DetectedSpecies,
    genomics_assay: Option<&DetectedAssay>,
    epigenomics_assay: Option<&DetectedAssay>,
    somatic_burden_per_mb: Option<f64>,
    titv_ratio: Option<f64>,
    tsv_scan: Option<&TsvScan>,
    vcf_path: Option<&Path>,
    bed_path: Option<&Path>,
    atac_path: Option<&Path>,
) -> Option<PresetSuggestion> {
    let mut candidates: Vec<PresetSuggestion> = Vec::new();

    // ATAC: strongest signal — explicit atac path given
    if atac_path.is_some() {
        candidates.push(PresetSuggestion {
            preset: "atac".into(),
            confidence: 0.90,
            reasons: vec!["ATAC-seq narrowPeak file provided".into()],
        });
    }

    // Plant
    if *species == DetectedSpecies::Plant {
        candidates.push(PresetSuggestion {
            preset: "plant".into(),
            confidence: 0.95,
            reasons: vec!["Plant chromosome names detected (Chr1–Chr5 pattern)".into()],
        });
    }

    // WGBS
    if let Some(DetectedAssay::WholeGenomeBisulfite) = epigenomics_assay {
        candidates.push(PresetSuggestion {
            preset: "wgbs".into(),
            confidence: 0.90,
            reasons: vec!["WGBS pattern detected in BED (uniform per-chrom distribution)".into()],
        });
    }

    // Clinical: human + WES
    if *species == DetectedSpecies::Human {
        if let Some(DetectedAssay::WholeExomeSeq) = genomics_assay {
            candidates.push(PresetSuggestion {
                preset: "clinical".into(),
                confidence: 0.80,
                reasons: vec![
                    "Human chromosomes detected".into(),
                    "WES capture pattern detected (high per-chrom variant CV)".into(),
                ],
            });
        }
    }

    // Cancer: somatic burden OR aberrant Ti/Tv
    {
        let mut cancer_reasons: Vec<String> = Vec::new();
        let mut cancer_signals = 0u32;

        if let Some(burden) = somatic_burden_per_mb {
            if burden > 10.0 {
                cancer_reasons.push(format!(
                    "High somatic burden estimate {:.1} variants/Mb (>10 suggests cancer)",
                    burden
                ));
                cancer_signals += 2;
            }
        }
        if let Some(titv) = titv_ratio {
            if !(1.8..=2.5).contains(&titv) {
                cancer_reasons.push(format!(
                    "Ti/Tv ratio {:.2} outside normal range [1.8, 2.5]",
                    titv
                ));
                cancer_signals += 1;
            }
        }

        if cancer_signals > 0 {
            let confidence = match cancer_signals {
                1 => 0.55,
                2 => 0.75,
                _ => 0.85,
            };
            candidates.push(PresetSuggestion {
                preset: "cancer".into(),
                confidence,
                reasons: cancer_reasons,
            });
        }
    }

    // RNA-seq: only TSV provided, no VCF and no BED
    if vcf_path.is_none() && bed_path.is_none() && tsv_scan.is_some() {
        candidates.push(PresetSuggestion {
            preset: "rna-seq".into(),
            confidence: 0.85,
            reasons: vec!["Only transcriptomics input provided — no genomics/epigenomics".into()],
        });
    }

    // Pick best candidate by confidence; require ≥0.5 to suggest
    candidates
        .into_iter()
        .filter(|c| c.confidence >= 0.5)
        .max_by(|a, b| a.confidence.partial_cmp(&b.confidence).unwrap_or(std::cmp::Ordering::Equal))
}
