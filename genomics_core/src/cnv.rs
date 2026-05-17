use ahash::AHashMap;
use serde::{Deserialize, Serialize};

use biomics_core::parse::{info_value_bytes, parse_f32, ByteLines, TabFields};

/// Copy-number state classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CnvClass {
    /// CN = 0: homozygous deletion.
    HomozygousDeletion,
    /// CN = 1: heterozygous deletion (loss of one copy).
    HeterozygousDeletion,
    /// CN = 2: diploid reference state.
    Diploid,
    /// CN 3–4: low-level amplification.
    LowAmplification,
    /// CN >= 5: high-level amplification.
    HighAmplification,
}

impl CnvClass {
    pub fn from_cn(cn: f32) -> Self {
        if cn < 0.5 {
            CnvClass::HomozygousDeletion
        } else if cn < 1.5 {
            CnvClass::HeterozygousDeletion
        } else if cn < 2.5 {
            CnvClass::Diploid
        } else if cn < 4.5 {
            CnvClass::LowAmplification
        } else {
            CnvClass::HighAmplification
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            CnvClass::HomozygousDeletion => "homdel",
            CnvClass::HeterozygousDeletion => "hetdel",
            CnvClass::Diploid => "diploid",
            CnvClass::LowAmplification => "lowamp",
            CnvClass::HighAmplification => "highamp",
        }
    }
}

/// A copy-number variant extracted from VCF CN/CNA/CNVTYPE INFO fields.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CnvRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub copy_number: f32,
    pub class: CnvClass,
    /// Gene(s) overlapping this CNV segment, if annotated in INFO.
    pub gene: Option<String>,
    /// Log2 ratio (log2(CN/2)), if available from LOG2 INFO field.
    pub log2_ratio: Option<f32>,
}

impl CnvRecord {
    pub fn length(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }
}

/// Per-chromosome CNV burden summary.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ChromCnvStats {
    pub total_segments: u64,
    pub deleted_bp: u64,
    pub amplified_bp: u64,
    pub diploid_bp: u64,
    pub mean_copy_number: f64,
}

/// Copy-number summary for the whole sample.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CnvSummary {
    pub total_segments: u64,
    pub homdel_count: u64,
    pub hetdel_count: u64,
    pub diploid_count: u64,
    pub lowamp_count: u64,
    pub highamp_count: u64,
    pub total_deleted_bp: u64,
    pub total_amplified_bp: u64,
    pub per_chrom: AHashMap<String, ChromCnvStats>,
    /// High-priority amplifications: CN >= 5, includes gene name.
    pub driver_amplifications: Vec<CnvRecord>,
    /// High-priority deletions: homozygous, includes gene name.
    pub driver_deletions: Vec<CnvRecord>,
    /// Estimated ploidy (weighted mean CN across all segments).
    pub estimated_ploidy: f64,
    /// Fraction of the genome altered (deleted + amplified).
    pub fraction_genome_altered: f64,
}

/// Extract CNV records from VCF INFO fields.
///
/// Handles the following VCF CNV representations:
/// 1. SVTYPE=CNV/DEL/DUP with CN= INFO field
/// 2. CNVTYPE= with CN= (Sequenza / PURPLE style)
/// 3. Standard FORMAT/CN or FORMAT/TCN fields
/// 4. INFO/CNA with numeric copy number
///
/// Lines without any CN information are skipped silently.
pub fn parse_cnv_from_vcf_line(line: &[u8]) -> Option<CnvRecord> {
    let mut fields = TabFields::new(line);

    let chrom_b = fields.next()?;
    let start = parse_f32(fields.next()?)? as u64;
    let _ = fields.next()?; // ID
    let _ = fields.next()?; // REF
    let _ = fields.next()?; // ALT
    let _ = fields.next()?; // QUAL
    let _ = fields.next()?; // FILTER
    let info = fields.next().unwrap_or(b".");

    // Require either SVTYPE=CNV/DEL/DUP or a CN= value present
    let svtype = info_value_bytes(info, b"SVTYPE");
    let has_cnv_svtype =
        svtype.is_some_and(|sv| matches!(sv, b"CNV" | b"DEL" | b"DUP" | b"GAIN" | b"LOSS"));

    // Extract CN value: try CN=, CNA=, TCN=, CNVTYPE= numeric
    let cn_opt = info_value_bytes(info, b"CN")
        .or_else(|| info_value_bytes(info, b"CNA"))
        .or_else(|| info_value_bytes(info, b"TCN"))
        .and_then(parse_f32);

    if cn_opt.is_none() && !has_cnv_svtype {
        return None;
    }

    let inferred_cn = match svtype {
        Some(b"DUP") | Some(b"GAIN") => 3.0_f32,
        Some(b"DEL") | Some(b"LOSS") => 1.0,
        _ => 2.0,
    };
    let copy_number = cn_opt.unwrap_or(inferred_cn);

    // END coordinate
    let end = info_value_bytes(info, b"END")
        .and_then(biomics_core::parse::parse_u64)
        .unwrap_or(start + 1);

    let log2_ratio = info_value_bytes(info, b"LOG2")
        .or_else(|| info_value_bytes(info, b"RATIO"))
        .and_then(parse_f32);

    let gene = info_value_bytes(info, b"GENE")
        .or_else(|| info_value_bytes(info, b"Gene"))
        .and_then(|v| std::str::from_utf8(v).ok().map(|s| s.to_string()));

    let chrom = std::str::from_utf8(chrom_b).ok()?.to_string();

    Some(CnvRecord {
        chrom,
        start,
        end,
        copy_number,
        class: CnvClass::from_cn(copy_number),
        gene,
        log2_ratio,
    })
}

/// Parse all CNV records from a VCF file, ignoring non-CNV lines.
pub fn parse_cnv_vcf(path: &std::path::Path) -> anyhow::Result<Vec<CnvRecord>> {
    use anyhow::Context;
    use memmap2::Mmap;

    let file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open CNV VCF '{}'", path.display()))?;
    let mmap = unsafe { Mmap::map(&file) }
        .with_context(|| format!("Cannot mmap CNV VCF '{}'", path.display()))?;
    let _ = mmap.advise(memmap2::Advice::Sequential);

    let data = mmap.as_ref();
    let mut records = Vec::with_capacity(data.len() / 300);

    for line in ByteLines::new(data) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        if let Some(rec) = parse_cnv_from_vcf_line(line) {
            records.push(rec);
        }
    }

    log::info!(
        "Parsed {} CNV segments from '{}'",
        records.len(),
        path.display()
    );
    Ok(records)
}

/// Compute a `CnvSummary` from a slice of `CnvRecord`s.
pub fn summarize_cnv(records: &[CnvRecord]) -> CnvSummary {
    let mut summary = CnvSummary::default();
    let mut total_bp = 0u64;
    let mut weighted_cn_sum = 0.0f64;

    let mut per_chrom: AHashMap<String, ChromCnvStats> = AHashMap::new();

    for rec in records {
        summary.total_segments += 1;
        let bp = rec.length();
        total_bp += bp;
        weighted_cn_sum += rec.copy_number as f64 * bp as f64;

        let stats = per_chrom.entry(rec.chrom.clone()).or_default();
        stats.total_segments += 1;
        stats.mean_copy_number += rec.copy_number as f64;

        match rec.class {
            CnvClass::HomozygousDeletion => {
                summary.homdel_count += 1;
                summary.total_deleted_bp += bp;
                stats.deleted_bp += bp;
                if rec.gene.is_some() {
                    summary.driver_deletions.push(rec.clone());
                }
            }
            CnvClass::HeterozygousDeletion => {
                summary.hetdel_count += 1;
                summary.total_deleted_bp += bp;
                stats.deleted_bp += bp;
            }
            CnvClass::Diploid => {
                summary.diploid_count += 1;
                stats.diploid_bp += bp;
            }
            CnvClass::LowAmplification => {
                summary.lowamp_count += 1;
                summary.total_amplified_bp += bp;
                stats.amplified_bp += bp;
            }
            CnvClass::HighAmplification => {
                summary.highamp_count += 1;
                summary.total_amplified_bp += bp;
                stats.amplified_bp += bp;
                if rec.gene.is_some() {
                    summary.driver_amplifications.push(rec.clone());
                }
            }
        }
    }

    // Finalize per-chrom means
    for stats in per_chrom.values_mut() {
        if stats.total_segments > 0 {
            stats.mean_copy_number /= stats.total_segments as f64;
        }
    }

    summary.per_chrom = per_chrom;
    summary.estimated_ploidy = if total_bp > 0 {
        weighted_cn_sum / total_bp as f64
    } else {
        2.0
    };
    summary.fraction_genome_altered = if total_bp > 0 {
        (summary.total_deleted_bp + summary.total_amplified_bp) as f64 / total_bp as f64
    } else {
        0.0
    };

    // Keep top driver events
    summary.driver_amplifications.sort_unstable_by(|a, b| {
        b.copy_number
            .partial_cmp(&a.copy_number)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    summary.driver_amplifications.truncate(50);
    summary.driver_deletions.truncate(50);

    summary
}
