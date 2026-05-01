use std::path::Path;

use anyhow::{Context, Result};
use memmap2::Mmap;

use biomics_core::parse::{
    info_value_bytes, nth_pipe_field, parse_f32, parse_u64, ByteLines, TabFields,
};
use crate::types::{TiTvClass, VariantRecord};

/// Classify a single-nucleotide substitution as transition or transversion.
///
/// Transitions: A↔G (purine-purine) and C↔T (pyrimidine-pyrimidine).
#[inline(always)]
pub fn classify_titv(ref_base: u8, alt_base: u8) -> TiTvClass {
    match (ref_base | 0x20, alt_base | 0x20) {
        // lowercase-normalise then match transition pairs
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

/// Parse one VCF data line (byte slice, no leading `#`) into a `VariantRecord`.
///
/// All number parsing uses the fast-float / manual-u64 paths from `biomics_core::parse`;
/// no intermediate `String` allocations occur until we copy the final field values.
#[inline]
fn parse_vcf_line(line: &[u8]) -> Option<VariantRecord> {
    let mut fields = TabFields::new(line);

    let chrom_b = fields.next()?;
    let pos = parse_u64(fields.next()?)?;
    let _ = fields.next()?; // ID column — skip
    let ref_b = fields.next()?;
    let alt_b = fields.next()?;
    let qual_b = fields.next()?;
    let _ = fields.next()?; // FILTER — skip
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
        // multi-allelic or long ALT: check first allele for Ti/Tv
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

    // Parse AF: try AF= first, then AF1=
    let af = info_value_bytes(info, b"AF")
        .or_else(|| info_value_bytes(info, b"AF1"))
        .and_then(|v| {
            let first = memchr::memchr(b',', v).map_or(v, |n| &v[..n]);
            parse_f32(first)
        });

    // Parse GENE: try GENE=, Gene=, or ANN= field 3
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

    Some(VariantRecord { chrom, pos, ref_allele, alt_allele, qual, titv, af, gene })
}

/// Parse a VCF file into a `Vec<VariantRecord>` using memory-mapped I/O and
/// zero-alloc byte-level line/field parsing.
///
/// - `madvise(SEQUENTIAL)` hints the kernel to prefetch pages ahead.
/// - `Vec::with_capacity` pre-allocates based on file-size estimate (≈200 bytes/variant).
/// - No intermediate `String` per line — fields are borrowed from the mmap slice.
pub fn parse_vcf(path: &Path) -> Result<Vec<VariantRecord>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open VCF file '{}'", path.display()))?;

    // SAFETY: the file is not modified while this process holds the Mmap.
    let mmap = unsafe { Mmap::map(&file) }
        .with_context(|| format!("Cannot mmap VCF file '{}'", path.display()))?;

    // Advise sequential access: kernel will read-ahead pages before we need them.
    let _ = mmap.advise(memmap2::Advice::Sequential);

    let data = mmap.as_ref();
    let mut records = Vec::with_capacity(data.len() / 200);

    for line in ByteLines::new(data) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        match parse_vcf_line(line) {
            Some(record) => records.push(record),
            None => log::warn!(
                "Skipping unparseable VCF line: {}",
                String::from_utf8_lossy(&line[..line.len().min(80)])
            ),
        }
    }

    log::info!("Parsed {} variant records from '{}'", records.len(), path.display());
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_titv() {
        assert_eq!(classify_titv(b'A', b'G'), TiTvClass::Transition);
        assert_eq!(classify_titv(b'C', b'T'), TiTvClass::Transition);
        assert_eq!(classify_titv(b'A', b'C'), TiTvClass::Transversion);
        assert_eq!(classify_titv(b'G', b'T'), TiTvClass::Transversion);
    }

    #[test]
    fn test_parse_vcf_line() {
        let line = b"chr1\t100\t.\tA\tG\t50.0\tPASS\tAF=0.35;GENE=BRCA1";
        let rec = parse_vcf_line(line).unwrap();
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.pos, 100);
        assert!((rec.qual - 50.0).abs() < 1e-4);
        assert_eq!(rec.titv, TiTvClass::Transition);
        assert!((rec.af.unwrap() - 0.35).abs() < 1e-4);
        assert_eq!(rec.gene, Some("BRCA1".to_string()));
    }
}
