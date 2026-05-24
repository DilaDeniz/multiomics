use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array2;
use scrna_core::umap::{run_umap, compute_fuzzy_graph_from_knn};
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

fn bench_umap_cpu(c: &mut Criterion) {
    let mut group = c.benchmark_group("umap_cpu");
    // Measure wall time with a tight sample count (UMAP is slow at large n).
    group.sample_size(10);

    for &n in &[200usize, 500, 1_000, 2_000] {
        let data = make_data(n, 20, 42);
        group.bench_with_input(
            BenchmarkId::new("run_umap", n),
            &n,
            |b, _| {
                b.iter(|| {
                    run_umap(&data, 15, 100, 0.1, 1.0, 0).unwrap()
                })
            },
        );
    }
    group.finish();
}

fn bench_umap_gpu(c: &mut Criterion) {
    let mut group = c.benchmark_group("umap_gpu");
    group.sample_size(10);

    // Same sizes — on a machine without a real GPU this exercises the CPU
    // fallback for n < 5000, demonstrating identical performance to run_umap.
    // On a machine WITH a GPU, sizes >= 5000 would show the speedup.
    for &n in &[200usize, 500, 1_000, 2_000] {
        let data = make_data(n, 20, 42);
        group.bench_with_input(
            BenchmarkId::new("run_umap_gpu", n),
            &n,
            |b, _| {
                b.iter(|| {
                    run_umap_gpu(&data, 15, 100, 0.1, 1.0, 0).unwrap()
                })
            },
        );
    }
    group.finish();
}

fn bench_fuzzy_graph_from_knn(c: &mut Criterion) {
    // Micro-benchmark: isolate Phase 1 graph construction from pre-built KNN.
    let mut group = c.benchmark_group("fuzzy_graph");
    group.sample_size(20);

    for &n in &[500usize, 1_000, 2_000] {
        // Pre-build a synthetic KNN (nearest 15 for each cell).
        let knn: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|i| {
                (1..=15usize)
                    .map(|k| ((i + k) % n, k as f64 * 0.1))
                    .collect()
            })
            .collect();

        group.bench_with_input(
            BenchmarkId::new("compute_fuzzy_graph_from_knn", n),
            &n,
            |b, &n| {
                b.iter(|| compute_fuzzy_graph_from_knn(&knn, n))
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_umap_cpu,
    bench_umap_gpu,
    bench_fuzzy_graph_from_knn,
);
criterion_main!(benches);
