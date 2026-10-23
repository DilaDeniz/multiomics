//! Validate VCF summary statistics against `bcftools stats`.
//!
//! The reference outputs in `tests/golden/bcftools/` were produced by
//! `validation/bcftools_stats.sh` (bcftools 1.24) on:
//!
//! - `bcftools_edge_cases.vcf` — hand-written records covering multi-allelic
//!   SNV and SNV+indel sites, MNPs, symbolic alleles (`<DEL>`, `<NON_REF>`),
//!   spanning deletions (`*`), breakends, `ALT=.`, lower-case bases and
//!   non-PASS filters;
//! - `1kg_phase3_chr22_16.05-16.6Mb.vcf.gz` — a real BGZF-compressed slice of
//!   the 1000 Genomes phase 3 release (8 655 records, including multi-allelic
//!   `<CN…>` copy-number alleles).
//!
//! Per-type counts are per record, exactly as bcftools reports them; Ts/Tv is
//! counted per SNV ALT allele (bcftools' headline `ts/tv` column).

use std::collections::HashMap;
use std::path::PathBuf;

use genomics_core::analyze_vcf;

struct BcftoolsStats {
    sn: HashMap<String, u64>,
    ts: u64,
    tv: u64,
}

fn parse_bcftools_stats(text: &str) -> BcftoolsStats {
    let mut sn = HashMap::new();
    let (mut ts, mut tv) = (None, None);
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let cols: Vec<&str> = line.split('\t').collect();
        match cols[0] {
            "SN" => {
                let key = cols[2].trim_end_matches(':').to_string();
                sn.insert(key, cols[3].parse().unwrap());
            }
            "TSTV" => {
                ts = Some(cols[2].parse().unwrap());
                tv = Some(cols[3].parse().unwrap());
            }
            _ => {}
        }
    }
    BcftoolsStats {
        sn,
        ts: ts.expect("TSTV line"),
        tv: tv.expect("TSTV line"),
    }
}

fn data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../test_data/validation")
        .join(name)
}

fn assert_matches_bcftools(vcf: &str, golden: &str) -> genomics_core::GenomicsSummary {
    let reference = parse_bcftools_stats(golden);
    let summary = analyze_vcf(&data_path(vcf), None).expect("VCF should parse");

    let checks = [
        ("number of SNPs", summary.snp_count),
        ("number of MNPs", summary.mnp_count),
        ("number of indels", summary.indel_count),
        ("number of others", summary.other_count),
        ("number of multiallelic sites", summary.multiallelic_sites),
    ];
    for (key, ours) in checks {
        assert_eq!(ours, reference.sn[key], "{vcf}: {key}");
    }
    assert_eq!(summary.transitions, reference.ts, "{vcf}: transitions");
    assert_eq!(summary.transversions, reference.tv, "{vcf}: transversions");
    let expected_titv = reference.ts as f64 / reference.tv as f64;
    assert!(
        (summary.titv_ratio - expected_titv).abs() < 1e-12,
        "{vcf}: Ts/Tv {} vs bcftools {}",
        summary.titv_ratio,
        expected_titv
    );
    summary
}

#[test]
fn edge_cases_match_bcftools_stats() {
    let summary = assert_matches_bcftools(
        "bcftools_edge_cases.vcf",
        include_str!("../../tests/golden/bcftools/bcftools_edge_cases.stats.txt"),
    );
    // 24 records − 2 no-ALT records (`.`, `<NON_REF>`) − 1 spanning-deletion-only (`*`).
    assert_eq!(summary.total_variants, 21);
    // Indel-based cancer metrics must not count multi-allelic SNVs, symbolic
    // SVs or breakends as indels (bcftools: 5 indel records, one allele each).
    assert_eq!(summary.hrd.as_ref().unwrap().total_indels, 5);
    assert_eq!(summary.msi.as_ref().unwrap().total_indels, 5);
}

#[test]
fn thousand_genomes_slice_matches_bcftools_stats() {
    let golden =
        include_str!("../../tests/golden/bcftools/1kg_phase3_chr22_16.05-16.6Mb.stats.txt");
    let reference = parse_bcftools_stats(golden);
    let summary = assert_matches_bcftools("1kg_phase3_chr22_16.05-16.6Mb.vcf.gz", golden);
    assert_eq!(
        summary.total_variants,
        reference.sn["number of records"] - reference.sn["number of no-ALTs"]
    );
}
