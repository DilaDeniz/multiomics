//! Benchmark: HyperLogLog vs competitor cardinality estimation approaches.
//!
//! Compares four approaches for counting unique genomic positions:
//!
//!  A) python_counter_style  — HashMap<u64, u32> count map (Python Counter / R table)
//!  B) sorted_vec_dedup      — Vec::sort + dedup (common in older bioinformatics tools)
//!  C) ahashset_u64          — AHashSet<u64>, what most Rust tools use today
//!  D) hyperloglog           — Multiomics, fixed 16 KB, O(1) merge
//!
//! For WGS data (n = 10M), A and B become impractical; only C and D scale.
//! The merge benchmark is the rayon reduce cost (critical for parallel folds).
//!
//! Run:
//!   cargo bench --bench hll_vs_hashset -p biomics_core

use ahash::{AHashMap, AHashSet};
use biomics_core::HyperLogLog;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

#[inline(always)]
fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

// ── Insert throughput ─────────────────────────────────────────────────────────

fn bench_insert_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("unique_position_counting");

    for &n in &[1_000u64, 100_000, 1_000_000] {
        group.throughput(Throughput::Elements(n));

        // A: Python Counter / R table style — HashMap counting
        group.bench_with_input(
            BenchmarkId::new("A_python_counter_hashmap", n),
            &n,
            |b, &n| {
                b.iter(|| {
                    let mut map: AHashMap<u64, u32> = AHashMap::with_capacity(n as usize);
                    for i in 0..n {
                        *map.entry(splitmix64(i)).or_insert(0) += 1;
                    }
                    map.len()
                })
            },
        );

        // B: sorted Vec + dedup — common in C/awk bioinformatics pipelines
        group.bench_with_input(BenchmarkId::new("B_sorted_vec_dedup", n), &n, |b, &n| {
            b.iter(|| {
                let mut v: Vec<u64> = (0..n).map(splitmix64).collect();
                v.sort_unstable();
                v.dedup();
                v.len()
            })
        });

        // C: AHashSet<u64> — what most Rust bioinformatics tools use
        group.bench_with_input(
            BenchmarkId::new("C_ahashset_current_rust_tools", n),
            &n,
            |b, &n| {
                b.iter(|| {
                    let mut set: AHashSet<u64> = AHashSet::with_capacity(n as usize);
                    for i in 0..n {
                        set.insert(splitmix64(i));
                    }
                    set.len()
                })
            },
        );

        // D: HyperLogLog — Multiomics
        group.bench_with_input(
            BenchmarkId::new("D_hyperloglog_multiomics", n),
            &n,
            |b, &n| {
                b.iter(|| {
                    let mut hll = HyperLogLog::new();
                    for i in 0..n {
                        hll.insert_hashed(splitmix64(i));
                    }
                    hll.cardinality()
                })
            },
        );
    }

    group.finish();
}

// ── Merge cost (rayon parallel-fold reduce step) ──────────────────────────────
//
// This is the CRITICAL bottleneck for parallel VCF processing.
// Each rayon worker builds its own accumulator, then they're merged.
// With n=1M and 16 workers: the merge must handle 16 × 1M/16 = 1M elements.

fn bench_merge_cost(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_merge_cost");

    for &n in &[100_000u64, 1_000_000] {
        group.throughput(Throughput::Elements(n * 2));

        // Pre-build both halves
        let mut hll_a = HyperLogLog::new();
        let mut hll_b = HyperLogLog::new();
        let mut set_a: AHashSet<u64> = AHashSet::with_capacity(n as usize);
        let mut set_b: AHashSet<u64> = AHashSet::with_capacity(n as usize);
        let mut vec_a: Vec<u64> = Vec::with_capacity(n as usize);
        let mut vec_b: Vec<u64> = Vec::with_capacity(n as usize);

        for i in 0..n {
            hll_a.insert_hashed(splitmix64(i));
            set_a.insert(splitmix64(i));
            vec_a.push(splitmix64(i));
        }
        for i in n..2 * n {
            hll_b.insert_hashed(splitmix64(i));
            set_b.insert(splitmix64(i));
            vec_b.push(splitmix64(i));
        }
        vec_a.sort_unstable();
        vec_b.sort_unstable();

        // C: AHashSet extend — current approach before HLL
        group.bench_with_input(BenchmarkId::new("C_ahashset_extend", n), &n, |b, _| {
            b.iter(|| {
                let mut merged = set_a.clone();
                merged.extend(set_b.iter().copied());
                merged.len()
            })
        });

        // B: sorted Vec merge (merge-sort style)
        group.bench_with_input(BenchmarkId::new("B_sorted_vec_merge", n), &n, |b, _| {
            b.iter(|| {
                let mut merged = Vec::with_capacity(vec_a.len() + vec_b.len());
                merged.extend_from_slice(&vec_a);
                merged.extend_from_slice(&vec_b);
                merged.sort_unstable();
                merged.dedup();
                merged.len()
            })
        });

        // D: HyperLogLog merge — O(16384) register-max, always fast
        group.bench_with_input(BenchmarkId::new("D_hyperloglog_merge", n), &n, |b, _| {
            b.iter(|| {
                let mut merged = HyperLogLog::new();
                merged.merge(&hll_a);
                merged.merge(&hll_b);
                merged.cardinality()
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_insert_throughput, bench_merge_cost);
criterion_main!(benches);
