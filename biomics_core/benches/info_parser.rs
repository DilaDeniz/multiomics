//! Benchmark: VCF INFO field extraction — Multiomics vs competitor patterns.
//!
//! INFO parsing is the hot path in any VCF tool. This benchmark compares how
//! different tools extract key=value pairs from a VCF INFO string:
//!
//!  A) python_str_split   — .split(';') + .split('='), approximate Python/pysam
//!  B) bcftools_scan      — sequential byte scan per key (C bcftools / htslib)
//!  C) repeated_memchr    — info_value_bytes called N times (naive Rust)
//!  D) aho_corasick_1x    — InfoMultiParser, one pass for all keys (Multiomics)
//!
//! Run:
//!   cargo bench --bench info_parser -p biomics_core

use biomics_core::parse::{info_value_bytes, InfoMultiParser};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

// ── INFO test strings ─────────────────────────────────────────────────────────

const INFO_SPARSE: &[u8] = b"DP=45;AF=0.32;MQ=60;SVTYPE=DEL;END=5000;CN=1";
const INFO_DENSE: &[u8] =
    b"SVTYPE=CNV;CN=3;CNA=3;TCN=4;END=123456;LOG2=0.58;RATIO=0.61;GENE=BRCA1;DP=88;AF=0.51";
const INFO_LONG: &[u8] = b"DP=100;MQ=59;MQ0=0;ExcessHet=3.01;FS=0.000;SOR=0.693;QD=28.5;\
      BaseQRankSum=-0.32;ClippingRankSum=0.00;ReadPosRankSum=-0.132;\
      SVTYPE=DUP;CN=4;END=987654;LOG2=1.02;GENE=MYC;VQSLOD=4.5;CULPRIT=MQ";

static KEYS: &[&str] = &[
    "SVTYPE", "CN", "CNA", "TCN", "END", "LOG2", "RATIO", "GENE", "Gene",
];

// ── Competitor A: Python/pysam-style string splitting ─────────────────────────
// Simulates: {k: v for k, v in (f.split('=', 1) for f in info.split(';'))}
// This is what pysam, pyvcf, and most Python VCF parsers do.

fn python_split_style(info: &[u8]) -> Vec<(&[u8], &[u8])> {
    let mut fields = Vec::with_capacity(16);
    let mut pos = 0usize;
    while pos < info.len() {
        let end = memchr::memchr(b';', &info[pos..])
            .map(|n| pos + n)
            .unwrap_or(info.len());
        let field = &info[pos..end];
        if let Some(eq) = memchr::memchr(b'=', field) {
            fields.push((&field[..eq], &field[eq + 1..]));
        }
        pos = end + 1;
    }
    fields
}

// ── Competitor B: bcftools/htslib sequential scan ─────────────────────────────
// Simulates the C bcftools approach: iterate fields, strcmp per target key.
// (We use Rust byte comparison as a faithful port of the C pattern.)

fn bcftools_style<'a>(info: &'a [u8], target_keys: &[&[u8]]) -> Vec<Option<&'a [u8]>> {
    let mut results = vec![None; target_keys.len()];
    let mut pos = 0usize;
    while pos < info.len() {
        let end = memchr::memchr(b';', &info[pos..])
            .map(|n| pos + n)
            .unwrap_or(info.len());
        let field = &info[pos..end];
        if let Some(eq) = memchr::memchr(b'=', field) {
            let key = &field[..eq];
            let val = &field[eq + 1..];
            for (i, &tk) in target_keys.iter().enumerate() {
                if key == tk {
                    results[i] = Some(val);
                    break;
                }
            }
        }
        pos = end + 1;
    }
    results
}

// ── Competitor C: repeated info_value_bytes ───────────────────────────────────
// What naive Rust code (and most Rust VCF libs) does — call per-key scan N times.

fn repeated_memchr_9x(info: &[u8]) -> (Option<&[u8]>, Option<&[u8]>, Option<&[u8]>) {
    let svtype = info_value_bytes(info, b"SVTYPE");
    let cn = info_value_bytes(info, b"CN");
    let cna = info_value_bytes(info, b"CNA");
    let tcn = info_value_bytes(info, b"TCN");
    let end = info_value_bytes(info, b"END");
    let log2 = info_value_bytes(info, b"LOG2");
    let ratio = info_value_bytes(info, b"RATIO");
    let gene = info_value_bytes(info, b"GENE");
    let gene_lc = info_value_bytes(info, b"Gene");
    // Collapse to prevent dead-code elim while keeping the call overhead real
    (
        svtype.or(cn).or(cna).or(tcn),
        end.or(log2).or(ratio),
        gene.or(gene_lc),
    )
}

// ── Benchmark: per INFO string ────────────────────────────────────────────────

fn bench_per_record(c: &mut Criterion) {
    let parser = InfoMultiParser::new(KEYS);
    let target_keys_bytes: Vec<&[u8]> = KEYS.iter().map(|k| k.as_bytes()).collect();

    let mut group = c.benchmark_group("info_extraction_per_record");

    for (label, info) in [
        ("sparse_44b", INFO_SPARSE),
        ("dense_90b", INFO_DENSE),
        ("long_271b", INFO_LONG),
    ] {
        group.throughput(Throughput::Bytes(info.len() as u64));

        group.bench_with_input(
            BenchmarkId::new("A_python_pysam_split", label),
            &info,
            |b, &info| b.iter(|| black_box(python_split_style(black_box(info)))),
        );

        group.bench_with_input(
            BenchmarkId::new("B_bcftools_htslib_scan", label),
            &info,
            |b, &info| b.iter(|| black_box(bcftools_style(black_box(info), &target_keys_bytes))),
        );

        group.bench_with_input(
            BenchmarkId::new("C_naive_rust_repeated_memchr", label),
            &info,
            |b, &info| b.iter(|| black_box(repeated_memchr_9x(black_box(info)))),
        );

        group.bench_with_input(
            BenchmarkId::new("D_multiomics_aho_corasick", label),
            &info,
            |b, &info| b.iter(|| black_box(parser.extract(black_box(info)))),
        );
    }

    group.finish();
}

// ── Throughput: 10 000 records ────────────────────────────────────────────────

fn bench_record_throughput(c: &mut Criterion) {
    let parser = InfoMultiParser::new(KEYS);
    let target_keys_bytes: Vec<&[u8]> = KEYS.iter().map(|k| k.as_bytes()).collect();

    // Realistic mixed workload
    let infos: Vec<&[u8]> = (0..10_000)
        .map(|i| match i % 3 {
            0 => INFO_SPARSE,
            1 => INFO_DENSE,
            _ => INFO_LONG,
        })
        .collect();

    let mut group = c.benchmark_group("info_throughput_10k_records");
    group.throughput(Throughput::Elements(infos.len() as u64));

    group.bench_function("A_python_pysam_split", |b| {
        b.iter(|| {
            infos
                .iter()
                .map(|&info| python_split_style(info).len())
                .sum::<usize>()
        })
    });

    group.bench_function("B_bcftools_htslib_scan", |b| {
        b.iter(|| {
            infos
                .iter()
                .map(|&info| {
                    bcftools_style(info, &target_keys_bytes)
                        .iter()
                        .filter(|v: &&Option<&[u8]>| v.is_some())
                        .count()
                })
                .sum::<usize>()
        })
    });

    group.bench_function("C_naive_rust_repeated_memchr", |b| {
        b.iter(|| {
            infos
                .iter()
                .map(|&info| {
                    let (a, b, g) = repeated_memchr_9x(info);
                    a.is_some() as usize + b.is_some() as usize + g.is_some() as usize
                })
                .sum::<usize>()
        })
    });

    group.bench_function("D_multiomics_aho_corasick", |b| {
        b.iter(|| {
            infos
                .iter()
                .map(|&info| parser.extract(info).iter().filter(|v| v.is_some()).count())
                .sum::<usize>()
        })
    });

    group.finish();
}

criterion_group!(benches, bench_per_record, bench_record_throughput);
criterion_main!(benches);
