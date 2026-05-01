use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use memmap2::Mmap;

use crate::types::{TiTvClass, VariantRecord};

/// Classify a single-nucleotide substitution as transition or transversion.
///
/// Transitions: A↔G (purine-purine) and C↔T (pyrimidine-pyrimidine).
/// Everything else that is a single base substitution is a transversion.
pub fn classify_titv(ref_base: u8, alt_base: u8) -> TiTvClass {
    match (ref_base.to_ascii_uppercase(), alt_base.to_ascii_uppercase()) {
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => TiTvClass::Transition,
        (b'A', b'C')
        | (b'A', b'T')
        | (b'G', b'C')
        | (b'G', b'T')
        | (b'C', b'A')
        | (b'C', b'G')
        | (b'T', b'A')
        | (b'T', b'G') => TiTvClass::Transversion,
        _ => TiTvClass::Other,
    }
}

/// Parse one INFO field string (e.g. `"AF=0.42;DP=30;GENE=TP53"`) and extract
/// the value for a given key.
fn info_value(info: &str, key: &str) -> Option<String> {
    for field in info.split(';') {
        if let Some(rest) = field.strip_prefix(key) {
            if rest.starts_with('=') {
                return Some(rest[1..].to_string());
            }
        }
    }
    None
}

/// Parse a single VCF data line (tab-separated, no leading `#`) into a
/// `VariantRecord`. Returns `None` for lines that cannot be parsed.
fn parse_vcf_line(line: &str) -> Option<VariantRecord> {
    let mut cols = line.splitn(9, '\t');
    let chrom = cols.next()?.to_string();
    let pos: u64 = cols.next()?.parse().ok()?;
    let _id = cols.next()?; // ID column — not used
    let ref_allele = cols.next()?.to_string();
    let alt_allele = cols.next()?.to_string();
    let qual_str = cols.next()?;
    let qual: f32 = if qual_str == "." {
        0.0
    } else {
        qual_str.parse().ok()?
    };
    let _filter = cols.next()?;
    let info = cols.next().unwrap_or(".");

    // Determine Ti/Tv class
    let titv = if ref_allele.len() == 1 && alt_allele.len() == 1 {
        let ref_b = ref_allele.as_bytes()[0];
        // ALT may be comma-separated (multi-allelic); take first
        let alt_b = alt_allele.split(',').next()?.as_bytes()[0];
        classify_titv(ref_b, alt_b)
    } else if ref_allele.len() != alt_allele.len() {
        TiTvClass::Indel
    } else {
        TiTvClass::Other
    };

    let af = info_value(info, "AF")
        .or_else(|| info_value(info, "AF1"))
        .and_then(|v| v.split(',').next()?.parse::<f32>().ok());

    let gene = info_value(info, "GENE")
        .or_else(|| info_value(info, "Gene"))
        .or_else(|| {
            // Try ANN field: ANN=<allele>|<effect>|<impact>|<gene>|...
            info_value(info, "ANN").and_then(|ann| {
                ann.split('|').nth(3).map(|g| g.to_string())
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

/// Parse a VCF file into a `Vec<VariantRecord>` using memory-mapped I/O.
///
/// Header lines (starting with `#`) are skipped. Malformed data lines are
/// logged and skipped — they do not cause a fatal error.
///
/// # Errors
/// Returns an error when the file cannot be opened or memory-mapped.
pub fn parse_vcf(path: &Path) -> Result<Vec<VariantRecord>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open VCF file '{}'", path.display()))?;

    // SAFETY: the file is not modified while this process holds the Mmap.
    let mmap = unsafe { Mmap::map(&file) }
        .with_context(|| format!("Cannot mmap VCF file '{}'", path.display()))?;

    let reader = BufReader::new(mmap.as_ref());
    let mut records = Vec::new();
    let mut line_no = 0usize;

    for line in reader.lines() {
        line_no += 1;
        let line =
            line.with_context(|| format!("Read error at line {} of '{}'", line_no, path.display()))?;

        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        match parse_vcf_line(&line) {
            Some(record) => records.push(record),
            None => {
                log::warn!(
                    "Skipping unparseable VCF line {} in '{}': {}",
                    line_no,
                    path.display(),
                    &line[..line.len().min(80)]
                );
            }
        }
    }

    log::info!(
        "Parsed {} variant records from '{}'",
        records.len(),
        path.display()
    );
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
    fn test_info_value() {
        let info = "AF=0.42;DP=30;GENE=TP53";
        assert_eq!(info_value(info, "AF"), Some("0.42".to_string()));
        assert_eq!(info_value(info, "GENE"), Some("TP53".to_string()));
        assert_eq!(info_value(info, "MISSING"), None);
    }

    #[test]
    fn test_parse_vcf_line() {
        let line = "chr1\t100\t.\tA\tG\t50.0\tPASS\tAF=0.35;GENE=BRCA1";
        let rec = parse_vcf_line(line).unwrap();
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.pos, 100);
        assert_eq!(rec.qual, 50.0);
        assert_eq!(rec.titv, TiTvClass::Transition);
        assert_eq!(rec.af, Some(0.35));
        assert_eq!(rec.gene, Some("BRCA1".to_string()));
    }
}
