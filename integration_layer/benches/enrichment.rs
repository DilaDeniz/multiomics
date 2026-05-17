//! Benchmark: pathway enrichment — BioMultiOmics vs. competitor patterns.
//!
//! Compares three implementations for the "which pathways overlap my gene set?"
//! question, modelling what each major tool roughly does under the hood:
//!
//!  A) naive_per_pathway    — rebuild an AHashSet per pathway, then intersect.
//!                            Approximates clusterProfiler / GSEA Java.
//!  B) btreemap_scan        — BTreeMap sorted set + range scan.
//!                            Approximates bcftools / samtools auxdata patterns.
//!  C) inverted_index       — gene → pathway index, O(Q) lookup.
//!                            BioMultiOmics approach.
//!
//! Sweep: 50 / 500 / 5 000 pathways × 100 / 1 000 query genes.
//!
//! Run:
//!   cargo bench --bench enrichment -p integration_layer

use ahash::{AHashMap, AHashSet};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

// ── Synthetic dataset generation ─────────────────────────────────────────────

fn make_gene_id(i: usize) -> String {
    format!("GENE{i:05}")
}

/// Build `n_pathways` pathways of `genes_per_pathway` genes each, drawn from a
/// pool of `pool_size` genes. Genes are shared across pathways (realistic).
fn make_pathways(
    n_pathways: usize,
    genes_per_pathway: usize,
    pool_size: usize,
) -> Vec<(String, Vec<String>)> {
    (0..n_pathways)
        .map(|p| {
            let name = format!("PW{p:04}");
            let genes: Vec<String> = (0..genes_per_pathway)
                .map(|g| make_gene_id((p * 7 + g * 13) % pool_size))
                .collect();
            (name, genes)
        })
        .collect()
}

/// Build query gene set: `n_query` genes, ~30% expected overlap with pathways.
fn make_query(n_query: usize, pool_size: usize) -> Vec<String> {
    (0..n_query)
        .map(|i| make_gene_id((i * 17) % pool_size))
        .collect()
}

// ── Competitor A: naive per-pathway AHashSet rebuild ─────────────────────────
// Approximates clusterProfiler (R), GSEA Java, most Python enrichment tools.
// For each pathway: build a fresh AHashSet, call intersection.count().

fn naive_per_pathway(query: &[String], pathways: &[(String, Vec<String>)]) -> usize {
    let query_set: AHashSet<&str> = query.iter().map(|s| s.as_str()).collect();
    pathways
        .iter()
        .map(|(_, genes)| {
            let pw_set: AHashSet<&str> = genes.iter().map(|s| s.as_str()).collect();
            query_set.intersection(&pw_set).count()
        })
        .sum()
}

// ── Competitor B: BTreeMap sorted-set scan ────────────────────────────────────
// Approximates tools that use sorted data structures (bcftools-style, htslib).

fn btreemap_scan(query: &[String], pathways: &[(String, Vec<String>)]) -> usize {
    let query_set: std::collections::BTreeSet<&str> = query.iter().map(|s| s.as_str()).collect();
    pathways
        .iter()
        .map(|(_, genes)| {
            genes
                .iter()
                .filter(|g| query_set.contains(g.as_str()))
                .count()
        })
        .sum()
}

// ── BioMultiOmics: inverted gene index ────────────────────────────────────────
// Build gene → pathway indices map once, then O(Q) query.

fn inverted_index(query: &[String], pathways: &[(String, Vec<String>)]) -> usize {
    // Build inverted index (amortized over multiple queries in production)
    let mut index: AHashMap<&str, Vec<usize>> = AHashMap::new();
    for (i, (_, genes)) in pathways.iter().enumerate() {
        for g in genes {
            index.entry(g.as_str()).or_default().push(i);
        }
    }

    let mut counts = vec![0usize; pathways.len()];
    for gene in query {
        if let Some(idxs) = index.get(gene.as_str()) {
            for &i in idxs {
                counts[i] += 1;
            }
        }
    }
    counts.iter().sum()
}

// ── Competitor D: Vec<Vec<String>> linear scan ────────────────────────────────
// Approximates Python list-of-lists iteration (pandas / pure Python tools).

fn python_style_linear(query: &[String], pathways: &[(String, Vec<String>)]) -> usize {
    let query_set: AHashSet<&str> = query.iter().map(|s| s.as_str()).collect();
    pathways
        .iter()
        .map(|(_, genes)| {
            genes
                .iter()
                .filter(|g| query_set.contains(g.as_str()))
                .count()
        })
        .sum()
}

// ── Benchmark groups ──────────────────────────────────────────────────────────

fn bench_enrichment_scale(c: &mut Criterion) {
    let pool = 20_000usize; // realistic gene-pool size (~human genome)

    let mut group = c.benchmark_group("enrichment_vs_competitors");

    for &(n_pw, n_q) in &[(50, 100), (500, 500), (5_000, 1_000)] {
        let pathways = make_pathways(n_pw, 40, pool);
        let query = make_query(n_q, pool);
        let label = format!("{n_pw}pw_{n_q}q");

        // Total work units: pathways × genes looked up
        group.throughput(Throughput::Elements((n_pw * n_q) as u64));

        group.bench_with_input(
            BenchmarkId::new("A_naive_per_pathway_clusterProfiler", &label),
            &(&query, &pathways),
            |b, (q, pw)| b.iter(|| naive_per_pathway(q, pw)),
        );

        group.bench_with_input(
            BenchmarkId::new("B_btreemap_scan_bcftools_style", &label),
            &(&query, &pathways),
            |b, (q, pw)| b.iter(|| btreemap_scan(q, pw)),
        );

        group.bench_with_input(
            BenchmarkId::new("C_python_linear_scan", &label),
            &(&query, &pathways),
            |b, (q, pw)| b.iter(|| python_style_linear(q, pw)),
        );

        group.bench_with_input(
            BenchmarkId::new("D_bioomics_inverted_index", &label),
            &(&query, &pathways),
            |b, (q, pw)| b.iter(|| inverted_index(q, pw)),
        );
    }

    group.finish();
}

// ── At scale: 50K pathways (MSigDB-sized), 5K query genes ────────────────────

fn bench_msigdb_scale(c: &mut Criterion) {
    let pool = 25_000usize;
    let pathways = make_pathways(50_000, 50, pool);
    let query = make_query(5_000, pool);

    let mut group = c.benchmark_group("msigdb_scale_50k_pathways");
    group.sample_size(10); // few samples — this is intentionally slow for naive

    group.bench_function("A_naive_per_pathway", |b| {
        b.iter(|| naive_per_pathway(&query, &pathways))
    });

    group.bench_function("D_bioomics_inverted_index", |b| {
        b.iter(|| inverted_index(&query, &pathways))
    });

    group.finish();
}

criterion_group!(benches, bench_enrichment_scale, bench_msigdb_scale);
criterion_main!(benches);
