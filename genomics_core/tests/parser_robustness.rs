//! Robustness tests: the VCF parser must degrade gracefully on adversarial or
//! malformed input (skipping unparseable lines) rather than panicking. This
//! matters because inputs are frequently user-supplied and untrusted.

use std::io::Write;

use genomics_core::parse_vcf;

/// Write `content` to a uniquely named temp file and return its path.
fn temp_vcf(tag: &str, content: &[u8]) -> std::path::PathBuf {
    let mut path = std::env::temp_dir();
    // Combine pid + tag for a collision-resistant, deterministic-per-test name.
    path.push(format!(
        "multiomics_robust_{}_{}.vcf",
        std::process::id(),
        tag
    ));
    let mut f = std::fs::File::create(&path).expect("create temp vcf");
    f.write_all(content).expect("write temp vcf");
    path
}

#[test]
fn malformed_lines_are_skipped_not_panicked() {
    let content = b"##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t50\tPASS\tAF=0.3\n\
this is not a valid vcf line at all\n\
chr1\tNOT_A_NUMBER\t.\tA\tG\t50\tPASS\t.\n\
chr2\t200\t.\tC\tT\t.\tPASS\t.\n\
\t\t\t\t\n\
chr3\t300\n";
    let path = temp_vcf("malformed", content);
    let variants = parse_vcf(&path).expect("parser should not error on malformed lines");
    // The two well-formed records survive; the rest are skipped.
    assert_eq!(variants.len(), 2);
    assert_eq!(variants[0].pos, 100);
    assert_eq!(variants[1].pos, 200);
    let _ = std::fs::remove_file(&path);
}

#[test]
fn header_only_file_yields_no_variants() {
    let content = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let path = temp_vcf("header_only", content);
    let variants = parse_vcf(&path).expect("header-only file should parse to empty");
    assert!(variants.is_empty());
    let _ = std::fs::remove_file(&path);
}

#[test]
fn non_utf8_bytes_do_not_panic() {
    // Invalid UTF-8 in the CHROM/REF fields must be rejected per-line, not crash.
    let mut content: Vec<u8> = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec();
    content.extend_from_slice(&[0xff, 0xfe, b'\t']); // bad chrom
    content.extend_from_slice(b"100\t.\tA\tG\t50\tPASS\t.\n");
    content.extend_from_slice(b"chr1\t200\t.\tA\tG\t50\tPASS\t.\n");
    let path = temp_vcf("non_utf8", &content);
    let variants = parse_vcf(&path).expect("non-utf8 lines skipped, not panic");
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].pos, 200);
    let _ = std::fs::remove_file(&path);
}

#[test]
fn empty_file_yields_no_variants() {
    // A zero-byte upload must not crash the memory-mapped reader.
    let path = temp_vcf("empty", b"");
    let variants = parse_vcf(&path).expect("empty file should parse to empty, not error");
    assert!(variants.is_empty());
    let _ = std::fs::remove_file(&path);
}

#[test]
fn extreme_position_values_are_handled() {
    let content = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chrX\t18446744073709551615\t.\tA\tG\t99\tPASS\t.\n";
    let path = temp_vcf("extreme_pos", content);
    let variants = parse_vcf(&path).expect("u64::MAX position should parse");
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].pos, u64::MAX);
    let _ = std::fs::remove_file(&path);
}
