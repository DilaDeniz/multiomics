use std::path::Path;

use anyhow::{Context, Result};

use crate::types::{site_flags, Impact, TiTvClass, VariantKind, VariantRecord};
use biomics_core::parse::{
    info_value_bytes, nth_pipe_field, parse_f32, parse_u64, ByteLines, TabFields,
};

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

/// Outcome of classifying one ALT allele against REF.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleClass {
    /// Not a variant: `.`, `*` (spanning deletion), `<*>`, `<NON_REF>`,
    /// mpileup's `X`/`<X>`, or an ALT identical to REF.
    NotVariant,
    /// A variant allele of the given class.
    Variant(VariantKind, TiTvClass),
}

/// Classify one ALT allele relative to REF.
///
/// Follows htslib's `bcf_set_variant_type` (used by `bcftools stats`): shared
/// leading bases are trimmed, then shared trailing bases (never past the first
/// remaining base). Two deliberate refinements: case is ignored throughout, and
/// an SNV embedded in longer alleles (`AC>AT`) is given the Ti/Tv class of the
/// base that actually changed rather than of the first base.
pub fn classify_allele(r: &[u8], a: &[u8]) -> AlleleClass {
    use AlleleClass::{NotVariant, Variant};

    if a.is_empty() || r.is_empty() || a == b"." || a == b"*" {
        return NotVariant;
    }
    if r.len() == 1 && a.len() == 1 {
        if a[0] == b'X' || r[0].eq_ignore_ascii_case(&a[0]) {
            return NotVariant;
        }
        return Variant(VariantKind::Snv, classify_titv(r[0], a[0]));
    }
    if a[0] == b'<' {
        return match a {
            b"<*>" | b"<X>" | b"<NON_REF>" => NotVariant,
            _ => Variant(VariantKind::Other, TiTvClass::Other),
        };
    }
    let is_breakend = |b: &[u8]| b.contains(&b'[') || b.contains(&b']');
    if is_breakend(a) || a[0] == b'.' || a[a.len() - 1] == b'.' {
        return Variant(VariantKind::Breakend, TiTvClass::Other);
    }

    let prefix = r
        .iter()
        .zip(a)
        .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
        .count();
    let (rr, aa) = (&r[prefix..], &a[prefix..]);
    match (rr.is_empty(), aa.is_empty()) {
        (true, true) => return NotVariant,
        (true, false) | (false, true) => return Variant(VariantKind::Indel, TiTvClass::Indel),
        (false, false) => {}
    }

    let (mut re, mut ae) = (rr.len() - 1, aa.len() - 1);
    while re > 0 && ae > 0 && rr[re].eq_ignore_ascii_case(&aa[ae]) {
        re -= 1;
        ae -= 1;
    }
    if re == 0 && ae == 0 {
        Variant(VariantKind::Snv, classify_titv(rr[0], aa[0]))
    } else if re == 0 || ae == 0 {
        // One side is down to a single base: an indel if that base still
        // matches (e.g. ACGT>AT after prefix trimming), otherwise complex.
        if rr[re].eq_ignore_ascii_case(&aa[ae]) {
            Variant(VariantKind::Indel, TiTvClass::Indel)
        } else {
            Variant(VariantKind::Other, TiTvClass::Other)
        }
    } else if re == ae {
        Variant(VariantKind::Mnp, TiTvClass::Other)
    } else {
        Variant(VariantKind::Other, TiTvClass::Other)
    }
}

/// Return the `i`-th comma-separated value of `v`.
#[inline]
fn nth_comma_field(v: &[u8], i: usize) -> Option<&[u8]> {
    v.split(|&b| b == b',').nth(i)
}

/// Layout of the VEP `CSQ` INFO field, read from its `##INFO` header line
/// (`... Format: Allele|Consequence|IMPACT|SYMBOL|...`).
#[derive(Debug, Clone, Default)]
pub struct CsqLayout {
    allele: usize,
    impact: Option<usize>,
    symbol: Option<usize>,
}

impl CsqLayout {
    /// Parse the `Format:` list of a `##INFO=<ID=CSQ,...>` header line.
    pub fn from_header(line: &[u8]) -> Option<Self> {
        let text = std::str::from_utf8(line).ok()?;
        if !text.starts_with("##INFO=<ID=CSQ,") {
            return None;
        }
        let fmt = text.split("Format: ").nth(1)?;
        let fmt = fmt.trim_end_matches(">").trim_end_matches('"');
        let fields: Vec<&str> = fmt.split('|').map(str::trim).collect();
        Some(Self {
            allele: fields.iter().position(|f| *f == "Allele")?,
            impact: fields.iter().position(|f| *f == "IMPACT"),
            symbol: fields.iter().position(|f| *f == "SYMBOL"),
        })
    }
}

/// The allele with its leading base removed (`-` if nothing remains), as VEP
/// writes indel alleles when every allele at a site shares the first base.
fn trimmed_allele(alt: &[u8]) -> Vec<u8> {
    if alt.len() > 1 {
        alt[1..].to_vec()
    } else {
        b"-".to_vec()
    }
}

/// Gene and most severe impact for one ALT allele from per-transcript
/// annotation entries (comma-separated, `|`-delimited fields). The candidate
/// spellings in `alleles` are tried in order; if no entry matches any of them,
/// all entries are used.
fn allele_annotation<'a>(
    value: &'a [u8],
    alleles: &[&[u8]],
    allele_field: usize,
    gene_field: Option<usize>,
    impact_field: Option<usize>,
) -> (Option<&'a [u8]>, Option<Impact>) {
    let entries: Vec<&[u8]> = value.split(|&b| b == b',').collect();
    let matching: Vec<&[u8]> = alleles
        .iter()
        .map(|allele| {
            entries
                .iter()
                .copied()
                .filter(|e| nth_pipe_field(e, allele_field) == Some(*allele))
                .collect::<Vec<_>>()
        })
        .find(|m| !m.is_empty())
        .unwrap_or_default();
    let chosen = if matching.is_empty() {
        entries
    } else {
        matching
    };
    let gene = gene_field.and_then(|gf| {
        chosen
            .iter()
            .find_map(|e| nth_pipe_field(e, gf).filter(|g| !g.is_empty()))
    });
    let impact = impact_field.and_then(|imf| {
        chosen
            .iter()
            .filter_map(|e| nth_pipe_field(e, imf).and_then(Impact::parse))
            .max()
    });
    (gene, impact)
}

/// Number of copies of allele `index` in a GT value such as `0/1` or `1|2`;
/// `None` if any allele is missing (`.`).
fn allele_dosage(gt: &[u8], index: usize) -> Option<u8> {
    let mut count = 0u8;
    for allele in gt.split(|&b| b == b'/' || b == b'|') {
        if allele == b"." || allele.is_empty() {
            return None;
        }
        let a: usize = std::str::from_utf8(allele).ok()?.parse().ok()?;
        if a == index {
            count = count.saturating_add(1);
        }
    }
    Some(count)
}

/// Parse one VCF data line, appending one [`VariantRecord`] per variant ALT
/// allele to `out` (multi-allelic records are split).
///
/// Returns `None` for a malformed line. Records whose ALT alleles are all
/// non-variant (`.`, `<NON_REF>`, `*` …) are valid but append nothing.
///
/// All number parsing uses the fast-float / manual-u64 paths from `biomics_core::parse`;
/// no intermediate `String` allocations occur until we copy the final field values.
#[inline]
fn parse_vcf_line(
    line: &[u8],
    csq: Option<&CsqLayout>,
    out: &mut Vec<VariantRecord>,
) -> Option<()> {
    let mut fields = TabFields::new(line);

    let chrom_b = fields.next()?;
    let pos = parse_u64(fields.next()?)?;
    let _ = fields.next()?; // ID column — skip
    let ref_b = fields.next()?;
    let alt_b = fields.next()?;
    let qual_b = fields.next()?;
    let _ = fields.next()?; // FILTER — skip
    let info = fields.next().unwrap_or(b".");
    // First sample's GT (FORMAT column, then the first sample column).
    let gt: Option<&[u8]> = match (fields.next(), fields.next()) {
        (Some(format), Some(sample)) => format
            .split(|&b| b == b':')
            .position(|k| k == b"GT")
            .and_then(|i| sample.split(|&b| b == b':').nth(i)),
        _ => None,
    };

    let chrom = std::str::from_utf8(chrom_b).ok()?;
    let ref_allele = std::str::from_utf8(ref_b).ok()?;
    std::str::from_utf8(alt_b).ok()?;

    let qual = if qual_b == b"." {
        0.0f32
    } else {
        parse_f32(qual_b).unwrap_or(0.0)
    };

    // AF (Number=A): one value per ALT allele; fall back to AF1 (bcftools call).
    let af_field = info_value_bytes(info, b"AF").or_else(|| info_value_bytes(info, b"AF1"));
    // Site-level gene annotation (GENE= / Gene=), else per-allele SnpEff ANN
    // (Allele|Annotation|Annotation_Impact|Gene_Name|...) or VEP CSQ.
    let site_gene = info_value_bytes(info, b"GENE")
        .or_else(|| info_value_bytes(info, b"Gene"))
        .and_then(|v| std::str::from_utf8(v).ok());
    let ann = info_value_bytes(info, b"ANN");
    let csq_value = csq.zip(info_value_bytes(info, b"CSQ"));

    let n_alts = alt_b.split(|&b| b == b',').count();
    let first_new = out.len();
    let mut flags = 0u8;

    for (i, alt) in alt_b.split(|&b| b == b',').enumerate() {
        let (kind, titv) = match classify_allele(ref_b, alt) {
            AlleleClass::NotVariant => continue,
            AlleleClass::Variant(kind, titv) => (kind, titv),
        };
        flags |= match kind {
            VariantKind::Snv => site_flags::SNV,
            VariantKind::Mnp => site_flags::MNP,
            VariantKind::Indel => site_flags::INDEL,
            VariantKind::Other => site_flags::OTHER,
            VariantKind::Breakend => 0,
        };

        let af = af_field
            .and_then(|v| nth_comma_field(v, i))
            .and_then(parse_f32);
        let (ann_gene, ann_impact) = ann.map_or((None, None), |a| {
            allele_annotation(a, &[alt], 0, Some(3), Some(2))
        });
        let (csq_gene, csq_impact) = csq_value.map_or((None, None), |(layout, v)| {
            // VEP writes the full ALT, or ALT without the shared first base
            // when every allele at the site starts with the same base.
            allele_annotation(
                v,
                &[alt, &trimmed_allele(alt)],
                layout.allele,
                layout.symbol,
                layout.impact,
            )
        });
        let gene = site_gene.map(str::to_string).or_else(|| {
            ann_gene
                .or(csq_gene)
                .and_then(|g| std::str::from_utf8(g).ok())
                .map(str::to_string)
        });
        let impact = ann_impact.or(csq_impact);
        let dosage = gt.and_then(|g| allele_dosage(g, i + 1));

        out.push(VariantRecord {
            chrom: chrom.to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: std::str::from_utf8(alt).ok()?.to_string(),
            qual,
            titv,
            af,
            gene,
            impact,
            dosage,
            kind,
            site_flags: 0,
        });
    }

    if let Some(first) = out.get_mut(first_new) {
        if n_alts > 1 {
            flags |= site_flags::MULTIALLELIC;
        }
        first.site_flags = flags | site_flags::FIRST;
    }
    Some(())
}

/// Parse a VCF file (plain, gzip or BGZF) into a `Vec<VariantRecord>` using
/// zero-alloc byte-level line/field parsing.
///
/// - Plain files are memory-mapped with `madvise(SEQUENTIAL)`; `.vcf.gz` files
///   are decompressed in memory (see [`biomics_core::io::open_input`]).
/// - `Vec::with_capacity` pre-allocates based on file-size estimate (≈200 bytes/variant).
/// - No intermediate `String` per line — fields are borrowed from the input slice.
pub fn parse_vcf(path: &Path) -> Result<Vec<VariantRecord>> {
    let input = biomics_core::io::open_input(path)
        .with_context(|| format!("Cannot read VCF file '{}'", path.display()))?;
    let data: &[u8] = &input;
    let mut records = Vec::with_capacity(data.len() / 200);
    let mut csq: Option<CsqLayout> = None;

    for line in ByteLines::new(data) {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'#' {
            if csq.is_none() {
                csq = CsqLayout::from_header(line);
            }
            continue;
        }
        match parse_vcf_line(line, csq.as_ref(), &mut records) {
            Some(()) => {}
            None => log::warn!(
                "Skipping unparseable VCF line: {}",
                String::from_utf8_lossy(&line[..line.len().min(80)])
            ),
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

    fn parse_one(line: &[u8]) -> Vec<VariantRecord> {
        let mut out = Vec::new();
        parse_vcf_line(line, None, &mut out).expect("line should parse");
        out
    }

    #[test]
    fn test_parse_vcf_line() {
        let recs = parse_one(b"chr1\t100\t.\tA\tG\t50.0\tPASS\tAF=0.35;GENE=BRCA1");
        assert_eq!(recs.len(), 1);
        let rec = &recs[0];
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.pos, 100);
        assert!((rec.qual - 50.0).abs() < 1e-4);
        assert_eq!(rec.titv, TiTvClass::Transition);
        assert_eq!(rec.kind, VariantKind::Snv);
        assert!((rec.af.unwrap() - 0.35).abs() < 1e-4);
        assert_eq!(rec.gene, Some("BRCA1".to_string()));
        assert_eq!(rec.site_flags, site_flags::FIRST | site_flags::SNV);
    }

    #[test]
    fn classify_allele_matches_htslib() {
        use AlleleClass::{NotVariant, Variant};
        let kind = |r: &str, a: &str| match classify_allele(r.as_bytes(), a.as_bytes()) {
            Variant(k, _) => Some(k),
            NotVariant => None,
        };
        assert_eq!(kind("A", "G"), Some(VariantKind::Snv));
        assert_eq!(kind("g", "t"), Some(VariantKind::Snv));
        assert_eq!(kind("A", "a"), None);
        assert_eq!(kind("A", "."), None);
        assert_eq!(kind("T", "*"), None);
        assert_eq!(kind("A", "<NON_REF>"), None);
        assert_eq!(kind("A", "<*>"), None);
        assert_eq!(kind("A", "<DEL>"), Some(VariantKind::Other));
        assert_eq!(kind("G", "G]chr2:5000]"), Some(VariantKind::Breakend));
        assert_eq!(kind("G", "[chr2:5000[G"), Some(VariantKind::Breakend));
        assert_eq!(kind("AC", "A"), Some(VariantKind::Indel));
        assert_eq!(kind("A", "ACGT"), Some(VariantKind::Indel));
        assert_eq!(kind("TAC", "AC"), Some(VariantKind::Indel));
        assert_eq!(kind("AC", "GT"), Some(VariantKind::Mnp));
        assert_eq!(kind("AC", "AT"), Some(VariantKind::Snv));
        assert_eq!(kind("AT", "G"), Some(VariantKind::Other));
        // An SNV inside longer alleles takes the Ti/Tv class of the changed base.
        assert_eq!(
            classify_allele(b"AC", b"AT"),
            Variant(VariantKind::Snv, TiTvClass::Transition)
        );
    }

    #[test]
    fn multiallelic_records_are_split_per_allele() {
        let recs = parse_one(b"chr1\t900\t.\tA\tG,T,<NON_REF>\t50\tPASS\tAF=0.3,0.2,0");
        assert_eq!(recs.len(), 2, "<NON_REF> must not produce a record");
        assert_eq!(recs[0].alt_allele, "G");
        assert_eq!(recs[1].alt_allele, "T");
        assert_eq!(recs[0].titv, TiTvClass::Transition);
        assert_eq!(recs[1].titv, TiTvClass::Transversion);
        assert!((recs[0].af.unwrap() - 0.3).abs() < 1e-6);
        assert!((recs[1].af.unwrap() - 0.2).abs() < 1e-6);
        assert_eq!(
            recs[0].site_flags,
            site_flags::FIRST | site_flags::SNV | site_flags::MULTIALLELIC
        );
        assert_eq!(recs[1].site_flags, 0);
    }

    #[test]
    fn mixed_site_sets_both_type_flags() {
        let recs = parse_one(b"chr1\t1100\t.\tA\tG,AT\t50\tPASS\t.");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[1].kind, VariantKind::Indel);
        let f = recs[0].site_flags;
        assert!(f & site_flags::SNV != 0 && f & site_flags::INDEL != 0);
    }

    #[test]
    fn no_alt_records_emit_nothing() {
        assert!(parse_one(b"chr1\t2000\t.\tC\t.\t50\tPASS\t.").is_empty());
        assert!(parse_one(b"chr1\t1700\t.\tA\t<NON_REF>\t.\tPASS\t.").is_empty());
        assert!(parse_one(b"chr1\t1800\t.\tT\t*\t50\tPASS\t.").is_empty());
    }

    #[test]
    fn ann_gene_is_chosen_per_allele() {
        let line = b"chr1\t5\t.\tA\tG,T\t50\tPASS\t\
ANN=G|missense_variant|MODERATE|KRAS|ENSG1|transcript|ENST1,T|stop_gained|HIGH|TP53|ENSG2|transcript|ENST2";
        let recs = parse_one(line);
        assert_eq!(recs[0].gene.as_deref(), Some("KRAS"));
        assert_eq!(recs[1].gene.as_deref(), Some("TP53"));
    }

    #[test]
    fn snpeff_impact_is_most_severe_per_allele() {
        let line = b"chr1\t5\t.\tA\tG,T\t50\tPASS\t\
ANN=G|missense_variant|MODERATE|KRAS|E1|transcript|T1,G|upstream_gene_variant|MODIFIER|KRAS|E1|transcript|T2,\
T|stop_gained|HIGH|TP53|E2|transcript|T3";
        let recs = parse_one(line);
        assert_eq!(recs[0].impact, Some(Impact::Moderate));
        assert_eq!(recs[1].impact, Some(Impact::High));
        assert_eq!(recs[1].gene.as_deref(), Some("TP53"));
    }

    #[test]
    fn vep_csq_uses_header_layout_and_vep_alleles() {
        let header = b"##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations \
from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene\">";
        let layout = CsqLayout::from_header(header).expect("CSQ header");
        // Mixed first bases (AT > A, GT): VEP keeps the alleles untrimmed.
        let line = b"chr1\t10\t.\tAT\tA,GT\t60\tPASS\t\
CSQ=A|frameshift_variant|HIGH|BRCA2|ENSG1,GT|synonymous_variant|LOW|BRCA2|ENSG1";
        let mut out = Vec::new();
        parse_vcf_line(line, Some(&layout), &mut out).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].impact, Some(Impact::High));
        assert_eq!(out[0].gene.as_deref(), Some("BRCA2"));
        assert_eq!(out[1].impact, Some(Impact::Low));
        // Shared first base (A > AT, AC): VEP writes "T" and "C".
        let shared = b"chr1\t20\t.\tA\tAT,AC\t60\tPASS\t\
CSQ=T|frameshift_variant|HIGH|ATM|ENSG2,C|missense_variant|MODERATE|ATM|ENSG2";
        let mut out2 = Vec::new();
        parse_vcf_line(shared, Some(&layout), &mut out2).unwrap();
        assert_eq!(out2[0].impact, Some(Impact::High));
        assert_eq!(out2[1].impact, Some(Impact::Moderate));
        // Without a CSQ header the field is ignored.
        let mut plain = Vec::new();
        parse_vcf_line(line, None, &mut plain).unwrap();
        assert!(plain.iter().all(|r| r.impact.is_none()));
    }

    #[test]
    fn genotype_dosage_is_per_alt_allele() {
        let recs = parse_one(b"chr1\t5\t.\tA\tG,T\t50\tPASS\t.\tGT:DP\t1|2:30");
        assert_eq!(recs[0].dosage, Some(1));
        assert_eq!(recs[1].dosage, Some(1));
        let hom = parse_one(b"chr1\t6\t.\tA\tG\t50\tPASS\t.\tDP:GT\t12:1/1");
        assert_eq!(hom[0].dosage, Some(2));
        let missing = parse_one(b"chr1\t7\t.\tA\tG\t50\tPASS\t.\tGT\t./.");
        assert_eq!(missing[0].dosage, None);
        assert_eq!(parse_one(b"chr1\t8\t.\tA\tG\t50\tPASS\t.")[0].dosage, None);
    }
}
