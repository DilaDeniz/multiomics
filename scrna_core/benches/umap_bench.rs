use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array2;
use scrna_core::umap::{compute_fuzzy_graph_from_knn, run_umap};
use scrna_core::umap_gpu::run_umap_gpu;

fn make_data(n_cells: usize, n_dims: usize, seed: u64) -> Array2<f64> {
    let mut s = seed ^ 0x9e37_79b9_7f4a_7c15u64;
    let mut rng = move || -> f64 {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    Array2::from_shape_fn((n_cells, n_dims), |_| rng())
}

/// CPU-only UMAP at sizes below and above the GPU threshold.
/// On a GPU machine these become the baseline to compare against.
fn bench_umap_cpu(c: &mut Criterion) {
    let mut group = c.benchmark_group("umap_cpu");
    group.sample_size(10);

    // Small (CPU-only range, below 5 000 threshold)
    for &n in &[500usize, 1_000, 2_000] {
        let data = make_data(n, 20, 42);
        group.bench_with_input(BenchmarkId::new("n", n), &n, |b, _| {
            b.iter(|| run_umap(&data, 15, 100, 0.1, 1.0, 0).unwrap())
        });
    }

    // Large (these trigger the GPU path when `--features gpu` is active)
    for &n in &[5_000usize, 10_000, 20_000] {
        let data = make_data(n, 20, 42);
        group.bench_with_input(BenchmarkId::new("n", n), &n, |b, _| {
            b.iter(|| run_umap(&data, 15, 100, 0.1, 1.0, 0).unwrap())
        });
    }

    group.finish();
}

/// GPU-accelerated UMAP (falls back to CPU if no adapter / n < 5 000).
/// Run with `cargo bench --features gpu` to activate the GPU path.
fn bench_umap_gpu(c: &mut Criterion) {
    let mut group = c.benchmark_group("umap_gpu");
    group.sample_size(10);

    // Below threshold — should match umap_cpu times exactly (same code path).
    for &n in &[500usize, 1_000, 2_000] {
        let data = make_data(n, 20, 42);
        group.bench_with_input(BenchmarkId::new("n", n), &n, |b, _| {
            b.iter(|| run_umap_gpu(&data, 15, 100, 0.1, 1.0, 0).unwrap())
        });
    }

    // Above threshold — GPU path kicks in when `--features gpu` is active.
    // On an RTX 4050 expect 5–10× end-to-end speedup vs umap_cpu above.
    for &n in &[5_000usize, 10_000, 20_000] {
        let data = make_data(n, 20, 42);
        group.bench_with_input(BenchmarkId::new("n", n), &n, |b, _| {
            b.iter(|| run_umap_gpu(&data, 15, 100, 0.1, 1.0, 0).unwrap())
        });
    }

    group.finish();
}

/// Micro-benchmark: fuzzy graph construction from pre-built KNN (Phase 1 only).
fn bench_fuzzy_graph(c: &mut Criterion) {
    let mut group = c.benchmark_group("fuzzy_graph");
    group.sample_size(20);

    for &n in &[1_000usize, 5_000, 10_000] {
        let knn: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|i| {
                (1..=15usize)
                    .map(|k| ((i + k) % n, k as f64 * 0.1))
                    .collect()
            })
            .collect();

        group.bench_with_input(BenchmarkId::new("n", n), &n, |b, &n| {
            b.iter(|| compute_fuzzy_graph_from_knn(&knn, n))
        });
    }

    group.finish();
}

criterion_group!(benches, bench_umap_cpu, bench_umap_gpu, bench_fuzzy_graph);
criterion_main!(benches);
