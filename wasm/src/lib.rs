use wasm_bindgen::prelude::*;

// ── In-memory parse helpers ──────────────────────────────────────────────────
// The existing parse_vcf / parse_tsv / parse_bed functions use memmap2 which
// is not available in wasm32. We re-implement the same byte-level logic here,
// operating directly on the &[u8] slice that wasm-bindgen hands us.

use biomics_core::parse::{
    info_value_bytes, nth_pipe_field, parse_f32, parse_f64, parse_u64, trim_bytes, ByteLines,
    TabFields,
};
use genomics_core::types::{TiTvClass, VariantRecord};
use transcriptomics_core::types::GeneRecord;
use epigenomics_core::types::MethylationRecord;

// ── VCF in-memory parser ─────────────────────────────────────────────────────

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

// ── TSV in-memory parser ─────────────────────────────────────────────────────

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

// ── BED in-memory parser ─────────────────────────────────────────────────────

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

// ── Sequential fold (wasm32 cannot use rayon's thread pool) ─────────────────

#[cfg(target_arch = "wasm32")]
fn sequential_fold<A>(records: &[A::Record]) -> anyhow::Result<A::Summary>
where
    A: biomics_core::accum::BatchAccum,
{
    let mut accum = A::default();
    for r in records {
        let _ = accum.process(r);
    }
    accum.finalize()
}

// ── Dispatch: parallel on native, sequential on wasm32 ──────────────────────

fn fold_genomics(
    records: &[VariantRecord],
) -> anyhow::Result<genomics_core::types::GenomicsSummary> {
    #[cfg(not(target_arch = "wasm32"))]
    {
        biomics_core::parallel_fold::<genomics_core::GenomicsAccum>(records, "genomics", None)
    }
    #[cfg(target_arch = "wasm32")]
    {
        sequential_fold::<genomics_core::GenomicsAccum>(records)
    }
}

fn fold_transcriptomics(
    records: &[GeneRecord],
) -> anyhow::Result<transcriptomics_core::types::TranscriptomicsSummary> {
    #[cfg(not(target_arch = "wasm32"))]
    {
        biomics_core::parallel_fold::<transcriptomics_core::TranscriptomicsAccum>(
            records,
            "transcriptomics",
            None,
        )
    }
    #[cfg(target_arch = "wasm32")]
    {
        sequential_fold::<transcriptomics_core::TranscriptomicsAccum>(records)
    }
}

fn fold_epigenomics(
    records: &[MethylationRecord],
) -> anyhow::Result<epigenomics_core::types::EpigenomicsSummary> {
    #[cfg(not(target_arch = "wasm32"))]
    {
        biomics_core::parallel_fold::<epigenomics_core::EpigenomicsAccum>(
            records,
            "epigenomics",
            None,
        )
    }
    #[cfg(target_arch = "wasm32")]
    {
        sequential_fold::<epigenomics_core::EpigenomicsAccum>(records)
    }
}

// ── Helper: JSON error object ─────────────────────────────────────────────────

fn error_json(msg: &str) -> String {
    format!("{{\"error\":{}}}", serde_json::Value::String(msg.to_string()))
}

// ── Public wasm-bindgen exports ───────────────────────────────────────────────

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
    match fold_genomics(&records) {
        Ok(summary) => serde_json::to_string(&summary).unwrap_or_default(),
        Err(e) => error_json(&e.to_string()),
    }
}

/// Parse and analyse an expression-matrix TSV from raw bytes.
/// Returns a JSON string of `TranscriptomicsSummary`, or `{"error":"..."}` on failure.
#[wasm_bindgen]
pub fn analyze_transcriptomics(data: &[u8]) -> String {
    let (records, sample_names) = parse_tsv_bytes(data);
    let n_samples = sample_names.len();
    match fold_transcriptomics(&records) {
        Ok(mut summary) => {
            summary.sample_names = sample_names;
            summary.sample_count = n_samples;
            // Run differential expression if there are enough samples.
            if n_samples >= 2 {
                let de = transcriptomics_core::differential_expression(&records);
                summary.diff_expr = Some(de);
            }
            serde_json::to_string(&summary).unwrap_or_default()
        }
        Err(e) => error_json(&e.to_string()),
    }
}

/// Parse and analyse a BED methylation file from raw bytes.
/// Returns a JSON string of `EpigenomicsSummary`, or `{"error":"..."}` on failure.
#[wasm_bindgen]
pub fn analyze_epigenomics(data: &[u8]) -> String {
    let records = parse_bed_bytes(data);
    match fold_epigenomics(&records) {
        Ok(summary) => serde_json::to_string(&summary).unwrap_or_default(),
        Err(e) => error_json(&e.to_string()),
    }
}

/// Run all three modalities and the integration layer from raw bytes.
/// Returns a JSON string of `{ genomics, transcriptomics, epigenomics, integration }`,
/// or `{"error":"..."}` on failure.
#[wasm_bindgen]
pub fn analyze_all(vcf: &[u8], tsv: &[u8], bed: &[u8]) -> String {
    // Parse
    let vcf_records = parse_vcf_bytes(vcf);
    let (tsv_records, sample_names) = parse_tsv_bytes(tsv);
    let bed_records = parse_bed_bytes(bed);

    // Fold each modality
    let genomics_summary = match fold_genomics(&vcf_records) {
        Ok(s) => s,
        Err(e) => return error_json(&format!("genomics: {}", e)),
    };

    let n_samples = sample_names.len();
    let transcr_summary = match fold_transcriptomics(&tsv_records) {
        Ok(mut s) => {
            s.sample_names = sample_names;
            s.sample_count = n_samples;
            if n_samples >= 2 {
                let de = transcriptomics_core::differential_expression(&tsv_records);
                s.diff_expr = Some(de);
            }
            s
        }
        Err(e) => return error_json(&format!("transcriptomics: {}", e)),
    };

    let epigen_summary = match fold_epigenomics(&bed_records) {
        Ok(s) => s,
        Err(e) => return error_json(&format!("epigenomics: {}", e)),
    };

    // Integration (skip heavy ML on wasm32 to keep bundle small)
    #[cfg(target_arch = "wasm32")]
    let skip_ml = true;
    #[cfg(not(target_arch = "wasm32"))]
    let skip_ml = false;

    let integration_summary =
        match integration_layer::run_integration(&genomics_summary, &transcr_summary, &epigen_summary, skip_ml) {
            Ok(s) => s,
            Err(e) => return error_json(&format!("integration: {}", e)),
        };

    // Bundle all four results into one object
    let combined = serde_json::json!({
        "genomics":       genomics_summary,
        "transcriptomics": transcr_summary,
        "epigenomics":    epigen_summary,
        "integration":    integration_summary,
    });

    serde_json::to_string(&combined).unwrap_or_default()
}
