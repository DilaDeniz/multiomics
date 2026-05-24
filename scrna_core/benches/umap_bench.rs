use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array2;
use scrna_core::umap::run_umap;
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

/// CPU UMAP — only small sizes to keep runtime reasonable.
fn bench_umap_cpu(c: &mut Criterion) {
    let mut group = c.benchmark_group("umap_cpu");
    group.sample_size(10);

    for &n in &[500usize, 1_000, 2_000, 5_000] {
        let data = make_data(n, 20, 42);
        group.bench_with_input(BenchmarkId::new("n", n), &n, |b, _| {
            b.iter(|| run_umap(&data, 15, 100, 0.1, 1.0, 0).unwrap())
        });
    }

    group.finish();
}

/// GPU UMAP — includes large sizes where GPU kicks in (n >= 5000).
/// Run with `cargo bench --features gpu` to activate the GPU path.
fn bench_umap_gpu(c: &mut Criterion) {
    let mut group = c.benchmark_group("umap_gpu");
    group.sample_size(10);

    for &n in &[5_000usize, 10_000, 20_000, 50_000] {
        let data = make_data(n, 20, 42);
        group.bench_with_input(BenchmarkId::new("n", n), &n, |b, _| {
            b.iter(|| run_umap_gpu(&data, 15, 100, 0.1, 1.0, 0).unwrap())
        });
    }

    group.finish();
}

criterion_group!(benches, bench_umap_cpu, bench_umap_gpu);
criterion_main!(benches);
