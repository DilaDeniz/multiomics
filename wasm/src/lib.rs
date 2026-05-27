//! WebAssembly bindings for the multiomics analysis engine.
//!
//! All parsing and accumulation is implemented inline here so that this crate
//! has no C-based transitive dependencies and compiles cleanly for the
//! `wasm32-unknown-unknown` target.
//!
//! The core crates (genomics_core, transcriptomics_core, epigenomics_core,
//! integration_layer) are intentionally NOT listed as dependencies here because
//! they pull in `needletail`/`xz2`/`lzma-sys` (C code) which cannot be
//! cross-compiled for wasm32.  All needed logic is re-implemented against the
//! pure-Rust `biomics_core` primitives only.

use wasm_bindgen::prelude::*;

use ahash::AHashMap;
use biomics_core::parse::{
    info_value_bytes, nth_pipe_field, parse_f32, parse_f64, parse_u64, trim_bytes, ByteLines,
    TabFields,
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// Types
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
enum TiTvClass {
    Transition,
    Transversion,
    Indel,
    Other,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct VariantRecord {
    chrom: String,
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    qual: f32,
    titv: TiTvClass,
    af: Option<f32>,
    gene: Option<String>,
}

#[derive(Debug, Clone)]
struct GeneRecord {
    gene_id: String,
    samples: Vec<f64>,
}

#[derive(Debug, Clone)]
struct MethylationRecord {
    chrom: String,
    start: u64,
    end: u64,
    methylation: f64,
    gene: Option<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Output summary types
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct ChromDensity {
    total: u64,
    snps: u64,
    indels: u64,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct GenomicsSummary {
    total_variants: u64,
    snp_count: u64,
    indel_count: u64,
    titv_ratio: f64,
    per_chrom: HashMap<String, ChromDensity>,
    af_histogram: Vec<u64>,
    unique_positions: u64,
    high_impact_genes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GeneStats {
    mean: f64,
    std: f64,
    max: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DiffExprResult {
    gene_id: String,
    log2_fold_change: f64,
    mean_s1: f64,
    mean_s2: f64,
    p_value: f64,
    padj: f64,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct TranscriptomicsSummary {
    total_genes: u64,
    expressed_genes: u64,
    low_expression_genes: Vec<String>,
    gene_stats: HashMap<String, GeneStats>,
    top_expressed: Vec<(String, f64)>,
    diff_expr: Option<Vec<DiffExprResult>>,
    sample_count: usize,
    sample_names: Vec<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct ChromMethylation {
    total_sites: u64,
    mean_methylation: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
enum RegionKind {
    Hypermethylated,
    Hypomethylated,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct MethylationRegion {
    chrom: String,
    start: u64,
    end: u64,
    mean_methylation: f64,
    kind: RegionKind,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct EpigenomicsSummary {
    total_sites: u64,
    global_methylation_pct: f64,
    per_chrom: HashMap<String, ChromMethylation>,
    hypermethylated: Vec<MethylationRegion>,
    hypomethylated: Vec<MethylationRegion>,
    gene_methylation: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct IntegrationSummary {
    high_impact_genes_in_genomics: Vec<String>,
    expressed_genes: u64,
    global_methylation_pct: f64,
    titv_ratio: f64,
    total_variants: u64,
    total_transcriptome_genes: u64,
    total_methylation_sites: u64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Parsers
// ─────────────────────────────────────────────────────────────────────────────

#[inline(always)]
fn classify_titv(ref_base: u8, alt_base: u8) -> TiTvClass {
    match (ref_base | 0x20, alt_base | 0x20) {
        (b'a', b'g') | (b'g', b'a') | (b'c', b't') | (b't', b'c') => TiTvClass::Transition,
        (b'a', b'c')
        | (b'a', b't')
        | (b'g', b'c')
        | (b'g', b't')
        | (b'c', b'a')
        | (b'c', b'g')
        | (b't', b'a')
        | (b't', b'g') => TiTvClass::Transversion,
        _ => TiTvClass::Other,
    }
}

fn parse_vcf_line(line: &[u8]) -> Option<VariantRecord> {
    let mut fields = TabFields::new(line);
    let chrom_b = fields.next()?;
    let pos = parse_u64(fields.next()?)?;
    let _ = fields.next()?; // ID
    let ref_b = fields.next()?;
    let alt_b = fields.next()?;
    let qual_b = fields.next()?;
    let _ = fields.next()?; // FILTER
    let info = fields.next().unwrap_or(b".");

    let chrom = std::str::from_utf8(chrom_b).ok()?.to_string();
    let ref_allele = std::str::from_utf8(ref_b).ok()?.to_string();
    let alt_allele = std::str::from_utf8(alt_b).ok()?.to_string();

    let qual = if qual_b == b"." {
        0.0f32
    } else {
        parse_f32(qual_b).unwrap_or(0.0)
    };

    let titv = if ref_b.len() == 1 && alt_b.len() == 1 {
        classify_titv(ref_b[0], alt_b[0])
    } else if ref_b.len() == 1 {
        let first_alt = memchr::memchr(b',', alt_b).map_or(alt_b, |n| &alt_b[..n]);
        if first_alt.len() == 1 {
            classify_titv(ref_b[0], first_alt[0])
        } else if ref_b.len() != first_alt.len() {
            TiTvClass::Indel
        } else {
            TiTvClass::Other
        }
    } else if ref_b.len() != alt_b.len() {
        TiTvClass::Indel
    } else {
        TiTvClass::Other
    };

    let af = info_value_bytes(info, b"AF")
        .or_else(|| info_value_bytes(info, b"AF1"))
        .and_then(|v| {
            let first = memchr::memchr(b',', v).map_or(v, |n| &v[..n]);
            parse_f32(first)
        });

    let gene = info_value_bytes(info, b"GENE")
        .or_else(|| info_value_bytes(info, b"Gene"))
        .and_then(|v| std::str::from_utf8(v).ok().map(|s| s.to_string()))
        .or_else(|| {
            info_value_bytes(info, b"ANN").and_then(|ann| {
                nth_pipe_field(ann, 3)
                    .and_then(|g| std::str::from_utf8(g).ok())
                    .filter(|s| !s.is_empty())
                    .map(|s| s.to_string())
            })
        });

    Some(VariantRecord {
        chrom,
        pos,
        ref_allele,
        alt_allele,
        qual,
        titv,
        af,
        gene,
    })
}

fn parse_vcf_bytes(data: &[u8]) -> Vec<VariantRecord> {
    let mut records = Vec::new();
    for line in ByteLines::new(data) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        if let Some(rec) = parse_vcf_line(line) {
            records.push(rec);
        }
    }
    records
}

fn parse_tsv_bytes(data: &[u8]) -> (Vec<GeneRecord>, Vec<String>) {
    let mut lines = ByteLines::new(data);
    let header = match lines.next() {
        Some(h) => h,
        None => return (Vec::new(), Vec::new()),
    };
    let header_fields: Vec<&[u8]> = TabFields::new(header).collect();
    if header_fields.len() < 2 {
        return (Vec::new(), Vec::new());
    }
    let sample_names: Vec<String> = header_fields[1..]
        .iter()
        .map(|f| String::from_utf8_lossy(trim_bytes(f)).into_owned())
        .collect();
    let n_samples = sample_names.len();

    let mut records = Vec::new();
    for line in lines {
        let line = trim_bytes(line);
        if line.is_empty() {
            continue;
        }
        let mut fields = TabFields::new(line);
        let gene_id_bytes = match fields.next() {
            Some(b) if !b.is_empty() => b,
            _ => continue,
        };
        let gene_id = String::from_utf8_lossy(trim_bytes(gene_id_bytes)).into_owned();
        let mut samples = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            let val = fields
                .next()
                .and_then(|b| parse_f64(trim_bytes(b)))
                .unwrap_or(0.0);
            samples.push(val);
        }
        records.push(GeneRecord { gene_id, samples });
    }
    (records, sample_names)
}

fn parse_bed_line(line: &[u8]) -> Option<MethylationRecord> {
    let mut cols = TabFields::new(line);
    let chrom = std::str::from_utf8(cols.next()?).ok()?.to_string();
    let start = parse_u64(cols.next()?)?;
    let end = parse_u64(cols.next()?)?;
    let col4 = cols.next();
    let col5 = cols.next();

    let gene: Option<String> = col4.and_then(|b| {
        let s = std::str::from_utf8(trim_bytes(b)).ok()?;
        if col5.is_some() {
            if s.parse::<f64>().is_err() && !s.is_empty() {
                Some(s.to_string())
            } else {
                None
            }
        } else {
            None
        }
    });

    let methylation = if let Some(score_b) = col5 {
        let score: f64 = parse_f64(trim_bytes(score_b))?;
        score / 10.0
    } else if let Some(val_b) = col4 {
        let val: f64 = parse_f64(trim_bytes(val_b))?;
        if val > 1.0 {
            val.min(100.0)
        } else {
            val * 100.0
        }
    } else {
        return None;
    };

    Some(MethylationRecord {
        chrom,
        start,
        end,
        methylation: methylation.clamp(0.0, 100.0),
        gene,
    })
}

fn parse_bed_bytes(data: &[u8]) -> Vec<MethylationRecord> {
    let mut records = Vec::new();
    for line in ByteLines::new(data) {
        let line = trim_bytes(line);
        if line.is_empty()
            || line[0] == b'#'
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }
        if let Some(rec) = parse_bed_line(line) {
            records.push(rec);
        }
    }
    records
}

// ─────────────────────────────────────────────────────────────────────────────
// Accumulators (inline, no C deps)
// ─────────────────────────────────────────────────────────────────────────────

fn fold_genomics(records: &[VariantRecord]) -> GenomicsSummary {
    let mut total = 0u64;
    let mut snps = 0u64;
    let mut indels = 0u64;
    let mut transitions = 0u64;
    let mut transversions = 0u64;
    let mut per_chrom: AHashMap<String, ChromDensity> = AHashMap::new();
    let mut af_histogram = [0u64; 20];
    let mut unique_positions: ahash::AHashSet<u64> = ahash::AHashSet::new();
    let mut high_impact_genes: Vec<String> = Vec::new();

    for r in records {
        total += 1;
        let entry = per_chrom.entry(r.chrom.clone()).or_default();
        entry.total += 1;

        match r.titv {
            TiTvClass::Transition => {
                snps += 1;
                transitions += 1;
                entry.snps += 1;
            }
            TiTvClass::Transversion => {
                snps += 1;
                transversions += 1;
                entry.snps += 1;
            }
            TiTvClass::Indel => {
                indels += 1;
                entry.indels += 1;
            }
            TiTvClass::Other => {}
        }

        if r.qual > 30.0 {
            if let Some(ref gene) = r.gene {
                high_impact_genes.push(gene.clone());
            }
        }

        if let Some(af) = r.af {
            let bin = ((af as f64).clamp(0.0, 0.9999) * 20.0) as usize;
            af_histogram[bin] += 1;
        }

        // Compact position key: FNV-1a hash of (chrom, pos)
        let mut h: u64 = 14_695_981_039_346_656_037;
        for b in r.chrom.bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(1_099_511_628_211);
        }
        h ^= r.pos;
        h = h.wrapping_mul(1_099_511_628_211);
        unique_positions.insert(h ^ (h >> 17));
    }

    let titv_ratio = if transversions == 0 {
        0.0
    } else {
        transitions as f64 / transversions as f64
    };

    high_impact_genes.sort_unstable();
    high_impact_genes.dedup();

    GenomicsSummary {
        total_variants: total,
        snp_count: snps,
        indel_count: indels,
        titv_ratio,
        per_chrom: per_chrom.into_iter().collect(),
        af_histogram: af_histogram.to_vec(),
        unique_positions: unique_positions.len() as u64,
        high_impact_genes,
    }
}

fn fold_transcriptomics(
    records: &[GeneRecord],
    sample_names: Vec<String>,
) -> TranscriptomicsSummary {
    let n_samples = sample_names.len();
    let mut gene_sums: AHashMap<String, Vec<f64>> = AHashMap::new();
    let mut gene_sq_sums: AHashMap<String, Vec<f64>> = AHashMap::new();
    let mut gene_counts: AHashMap<String, u64> = AHashMap::new();

    for r in records {
        let sums = gene_sums
            .entry(r.gene_id.clone())
            .or_insert_with(|| vec![0.0; r.samples.len()]);
        let sq_sums = gene_sq_sums
            .entry(r.gene_id.clone())
            .or_insert_with(|| vec![0.0; r.samples.len()]);
        let count = gene_counts.entry(r.gene_id.clone()).or_insert(0);
        for (i, &val) in r.samples.iter().enumerate() {
            if i < sums.len() {
                sums[i] += val;
                sq_sums[i] += val * val;
            }
        }
        *count += 1;
    }

    let mut gene_stats: HashMap<String, GeneStats> = HashMap::new();
    let mut expressed_genes = 0u64;
    let mut low_expression_genes = Vec::new();
    let mut mean_by_gene: Vec<(String, f64)> = Vec::new();

    for (gene, sums) in &gene_sums {
        let count = *gene_counts.get(gene.as_str()).unwrap_or(&1) as f64;
        let sq_sums_v = gene_sq_sums.get(gene.as_str());
        let means: Vec<f64> = sums.iter().map(|&s| s / count).collect();
        let overall_mean = if means.is_empty() {
            0.0
        } else {
            means.iter().sum::<f64>() / means.len() as f64
        };
        let std = if let Some(sq) = sq_sums_v {
            let var = sq
                .iter()
                .zip(sums.iter())
                .map(|(&sq_s, &s)| {
                    let m = s / count;
                    (sq_s / count - m * m).max(0.0)
                })
                .sum::<f64>()
                / sq.len().max(1) as f64;
            var.sqrt()
        } else {
            0.0
        };
        let max = means.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        gene_stats.insert(gene.clone(), GeneStats { mean: overall_mean, std, max });
        if overall_mean >= 1.0 {
            expressed_genes += 1;
        } else {
            low_expression_genes.push(gene.clone());
        }
        mean_by_gene.push((gene.clone(), overall_mean));
    }

    mean_by_gene.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let top_expressed: Vec<(String, f64)> = mean_by_gene.into_iter().take(100).collect();

    // Simple differential expression: compare first half of samples vs second half.
    let diff_expr = if n_samples >= 2 {
        let mid = n_samples / 2;
        let mut de_results: Vec<DiffExprResult> = records
            .iter()
            .map(|r| {
                let group1: Vec<f64> = r.samples[..mid.min(r.samples.len())].to_vec();
                let group2: Vec<f64> = r.samples[mid.min(r.samples.len())..].to_vec();
                let mean_s1 = if group1.is_empty() {
                    0.0
                } else {
                    group1.iter().sum::<f64>() / group1.len() as f64
                };
                let mean_s2 = if group2.is_empty() {
                    0.0
                } else {
                    group2.iter().sum::<f64>() / group2.len() as f64
                };
                let lfc = if mean_s1 > 0.0 {
                    (mean_s2 / mean_s1).max(f64::MIN_POSITIVE).log2()
                } else if mean_s2 > 0.0 {
                    f64::INFINITY
                } else {
                    0.0
                };
                DiffExprResult {
                    gene_id: r.gene_id.clone(),
                    log2_fold_change: lfc,
                    mean_s1,
                    mean_s2,
                    p_value: f64::NAN,
                    padj: f64::NAN,
                }
            })
            .collect();
        de_results.sort_unstable_by(|a, b| {
            b.log2_fold_change
                .abs()
                .partial_cmp(&a.log2_fold_change.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        Some(de_results)
    } else {
        None
    };

    TranscriptomicsSummary {
        total_genes: gene_sums.len() as u64,
        expressed_genes,
        low_expression_genes,
        gene_stats,
        top_expressed,
        diff_expr,
        sample_count: n_samples,
        sample_names,
    }
}

fn fold_epigenomics(records: &[MethylationRecord]) -> EpigenomicsSummary {
    let mut total_sites = 0u64;
    let mut sum_methylation = 0.0f64;
    let mut per_chrom_sums: AHashMap<String, (u64, f64)> = AHashMap::new();
    let mut chrom_sites: AHashMap<String, Vec<(u64, u64, f64)>> = AHashMap::new();
    let mut gene_sites: AHashMap<String, Vec<f64>> = AHashMap::new();

    for r in records {
        total_sites += 1;
        sum_methylation += r.methylation;
        let entry = per_chrom_sums.entry(r.chrom.clone()).or_insert((0, 0.0));
        entry.0 += 1;
        entry.1 += r.methylation;
        chrom_sites
            .entry(r.chrom.clone())
            .or_default()
            .push((r.start, r.end, r.methylation));
        if let Some(ref gene) = r.gene {
            gene_sites.entry(gene.clone()).or_default().push(r.methylation);
        }
    }

    let global_methylation_pct = if total_sites == 0 {
        0.0
    } else {
        sum_methylation / total_sites as f64
    };

    let per_chrom: HashMap<String, ChromMethylation> = per_chrom_sums
        .into_iter()
        .map(|(chrom, (n, sum))| {
            (
                chrom,
                ChromMethylation {
                    total_sites: n,
                    mean_methylation: if n == 0 { 0.0 } else { sum / n as f64 },
                },
            )
        })
        .collect();

    // Sliding-window region detection (window=5 sites)
    let mut hypermethylated = Vec::new();
    let mut hypomethylated = Vec::new();
    for (chrom, mut sites) in chrom_sites {
        sites.sort_unstable_by_key(|s| s.0);
        if sites.len() < 5 {
            continue;
        }
        for window in sites.windows(5) {
            let start = window[0].0;
            let end = window[4].1;
            let mean = window.iter().map(|s| s.2).sum::<f64>() / 5.0;
            if mean > 80.0 {
                hypermethylated.push(MethylationRegion {
                    chrom: chrom.clone(),
                    start,
                    end,
                    mean_methylation: mean,
                    kind: RegionKind::Hypermethylated,
                });
            } else if mean < 20.0 {
                hypomethylated.push(MethylationRegion {
                    chrom: chrom.clone(),
                    start,
                    end,
                    mean_methylation: mean,
                    kind: RegionKind::Hypomethylated,
                });
            }
        }
    }

    let gene_methylation: HashMap<String, f64> = gene_sites
        .into_iter()
        .filter_map(|(gene, vals)| {
            if vals.is_empty() {
                None
            } else {
                Some((gene, vals.iter().sum::<f64>() / vals.len() as f64))
            }
        })
        .collect();

    EpigenomicsSummary {
        total_sites,
        global_methylation_pct,
        per_chrom,
        hypermethylated,
        hypomethylated,
        gene_methylation,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn error_json(msg: &str) -> String {
    format!(
        "{{\"error\":{}}}",
        serde_json::Value::String(msg.to_string())
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Public wasm-bindgen exports
// ─────────────────────────────────────────────────────────────────────────────

/// Returns the crate version string.
#[wasm_bindgen]
pub fn get_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// Parse and analyse a VCF file from raw bytes.
/// Returns a JSON string of `GenomicsSummary`, or `{"error":"..."}` on failure.
#[wasm_bindgen]
pub fn analyze_genomics(data: &[u8]) -> String {
    let records = parse_vcf_bytes(data);
    let summary = fold_genomics(&records);
    serde_json::to_string(&summary).unwrap_or_else(|e| error_json(&e.to_string()))
}

/// Parse and analyse an expression-matrix TSV from raw bytes.
/// Returns a JSON string of `TranscriptomicsSummary`, or `{"error":"..."}` on failure.
#[wasm_bindgen]
pub fn analyze_transcriptomics(data: &[u8]) -> String {
    let (records, sample_names) = parse_tsv_bytes(data);
    let summary = fold_transcriptomics(&records, sample_names);
    serde_json::to_string(&summary).unwrap_or_else(|e| error_json(&e.to_string()))
}

/// Parse and analyse a BED methylation file from raw bytes.
/// Returns a JSON string of `EpigenomicsSummary`, or `{"error":"..."}` on failure.
#[wasm_bindgen]
pub fn analyze_epigenomics(data: &[u8]) -> String {
    let records = parse_bed_bytes(data);
    let summary = fold_epigenomics(&records);
    serde_json::to_string(&summary).unwrap_or_else(|e| error_json(&e.to_string()))
}

/// Run all three modalities and a lightweight integration summary from raw bytes.
/// Returns a JSON string of `{ genomics, transcriptomics, epigenomics, integration }`,
/// or `{"error":"..."}` on failure.
#[wasm_bindgen]
pub fn analyze_all(vcf: &[u8], tsv: &[u8], bed: &[u8]) -> String {
    let vcf_records = parse_vcf_bytes(vcf);
    let (tsv_records, sample_names) = parse_tsv_bytes(tsv);
    let bed_records = parse_bed_bytes(bed);

    let genomics = fold_genomics(&vcf_records);
    let transcriptomics = fold_transcriptomics(&tsv_records, sample_names);
    let epigenomics = fold_epigenomics(&bed_records);

    let integration = IntegrationSummary {
        high_impact_genes_in_genomics: genomics.high_impact_genes.clone(),
        expressed_genes: transcriptomics.expressed_genes,
        global_methylation_pct: epigenomics.global_methylation_pct,
        titv_ratio: genomics.titv_ratio,
        total_variants: genomics.total_variants,
        total_transcriptome_genes: transcriptomics.total_genes,
        total_methylation_sites: epigenomics.total_sites,
    };

    let combined = serde_json::json!({
        "genomics":       genomics,
        "transcriptomics": transcriptomics,
        "epigenomics":    epigenomics,
        "integration":    integration,
    });

    serde_json::to_string(&combined).unwrap_or_else(|e| error_json(&e.to_string()))
}
