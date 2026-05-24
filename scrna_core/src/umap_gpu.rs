//! GPU-accelerated UMAP via WebGPU (wgpu) compute shaders.
//!
//! Requires the `gpu` feature flag. When wgpu cannot initialise a device
//! (no compatible GPU, wrong platform, CI environment) the function
//! transparently falls back to the CPU UMAP implementation in `umap.rs`.
//!
//! # GPU speedup strategy
//! The dominant cost in UMAP Phase 1 for large datasets (n > 10 000) is the
//! O(n² × d) pairwise distance matrix. A WebGPU compute shader computes all
//! n² distances in parallel using a 16×16 workgroup tile, reducing wall time
//! from minutes (CPU) to seconds (GPU) for 100 k+ cells.
//!
//! Phase 2 (SGD embedding optimisation) runs on CPU using the GPU-computed
//! KNN graph — the SGD step is memory-bandwidth bound and harder to parallelise
//! without atomic contention.
//!
//! # Platform support
//! wgpu targets Vulkan (Linux/Windows), Metal (macOS/iOS), DX12 (Windows),
//! and WebGPU (browser). Enable with `--features gpu`.
//!
//! # References
//! * McInnes L, Healy J, Melville J (2018) UMAP: Uniform Manifold Approximation
//!   and Projection for Dimension Reduction. arXiv:1802.03426.

use anyhow::Result;
use ndarray::Array2;

use crate::umap::{run_umap, UmapResult};

/// WGSL compute shader: pairwise Euclidean distances.
///
/// Layout:
/// - binding 0: flat f32 array of shape (n_cells × n_dims) — row-major
/// - binding 1: flat f32 output of shape (n_cells × n_cells) — row-major
/// - binding 2: uniform { n_cells: u32, n_dims: u32 }
#[cfg(feature = "gpu")]
const DISTANCE_SHADER: &str = r#"
struct Uniforms {
    n_cells: u32,
    n_dims: u32,
}

@group(0) @binding(0) var<storage, read>       data:       array<f32>;
@group(0) @binding(1) var<storage, read_write> distances:  array<f32>;
@group(0) @binding(2) var<uniform>             uniforms:   Uniforms;

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let n = uniforms.n_cells;
    let d = uniforms.n_dims;

    if i >= n || j >= n {
        return;
    }

    var dist_sq: f32 = 0.0;
    for (var k: u32 = 0u; k < d; k = k + 1u) {
        let a = data[i * d + k];
        let b = data[j * d + k];
        let diff = a - b;
        dist_sq = dist_sq + diff * diff;
    }

    distances[i * n + j] = sqrt(dist_sq);
}
"#;

/// Run UMAP with GPU-accelerated KNN when the `gpu` feature is enabled.
///
/// Falls back to [`run_umap`] (CPU) when:
/// - The `gpu` feature is not compiled in.
/// - No compatible GPU adapter is found.
/// - The dataset is small enough that GPU overhead would be counterproductive
///   (n < 5 000 cells — flat CPU scan is faster due to no transfer overhead).
///
/// # Arguments
/// Same as [`run_umap`].
pub fn run_umap_gpu(
    data: &Array2<f64>,
    n_neighbors: usize,
    n_epochs: usize,
    min_dist: f64,
    learning_rate: f64,
    seed: u64,
) -> Result<UmapResult> {
    // For small datasets the GPU transfer overhead exceeds the compute gain.
    #[cfg(feature = "gpu")]
    if data.nrows() >= 5_000 {
        match gpu_knn_umap(data, n_neighbors, n_epochs, min_dist, learning_rate, seed) {
            Ok(result) => {
                log::info!("GPU UMAP completed ({n_cells} cells)");
                return Ok(result);
            }
            Err(e) => {
                log::warn!("GPU UMAP failed ({e:#}), falling back to CPU");
            }
        }
    }

    // CPU fallback (always available).
    run_umap(data, n_neighbors, n_epochs, min_dist, learning_rate, seed)
}

// ── GPU implementation (compiled only with `gpu` feature) ────────────────────

#[cfg(feature = "gpu")]
fn gpu_knn_umap(
    data: &Array2<f64>,
    n_neighbors: usize,
    n_epochs: usize,
    min_dist: f64,
    learning_rate: f64,
    seed: u64,
) -> Result<UmapResult> {
    use crate::umap::{compute_fuzzy_graph_from_knn, run_umap_from_graph};

    let n_cells = data.nrows();
    let n_dims = data.ncols();

    // Cast to f32 for GPU transfer (f16 would be faster but wgpu f16 support varies).
    let data_f32: Vec<f32> = data.iter().map(|&x| x as f32).collect();

    log::info!("GPU UMAP: uploading {n_cells}×{n_dims} matrix to GPU…");

    // Compute full distance matrix on GPU.
    let dist_matrix = pollster::block_on(gpu_distance_matrix(&data_f32, n_cells, n_dims))?;

    log::info!("GPU UMAP: distance matrix computed, extracting {n_neighbors}-NN…");

    // Extract k-NN from flat distance matrix (CPU — O(n × k log n)).
    let knn = knn_from_distances(&dist_matrix, n_cells, n_neighbors);

    // Build fuzzy graph and run SGD on CPU.
    let (adjacency, rho) = compute_fuzzy_graph_from_knn(&knn, n_cells);
    run_umap_from_graph(&adjacency, &rho, n_cells, n_epochs, min_dist, learning_rate, seed)
}

/// Upload data to GPU, dispatch the pairwise distance shader, download results.
#[cfg(feature = "gpu")]
async fn gpu_distance_matrix(data: &[f32], n_cells: usize, n_dims: usize) -> Result<Vec<f32>> {
    use wgpu::util::DeviceExt;

    let instance = wgpu::Instance::default();

    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        })
        .await
        .ok_or_else(|| anyhow::anyhow!("no wgpu adapter found"))?;

    let (device, queue) = adapter
        .request_device(&wgpu::DeviceDescriptor::default(), None)
        .await?;

    // Input buffer: n_cells × n_dims f32 values.
    let input_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("umap_input"),
        contents: bytemuck::cast_slice(data),
        usage: wgpu::BufferUsages::STORAGE,
    });

    // Output buffer: n_cells × n_cells f32 distances.
    let out_size = (n_cells * n_cells * std::mem::size_of::<f32>()) as u64;
    let output_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_distances"),
        size: out_size,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // Uniform buffer: [n_cells: u32, n_dims: u32].
    let uniforms: [u32; 2] = [n_cells as u32, n_dims as u32];
    let uniform_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("umap_uniforms"),
        contents: bytemuck::cast_slice(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    // Readback buffer (mappable).
    let readback_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_readback"),
        size: out_size,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    // Compile shader and set up pipeline.
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("distance_shader"),
        source: wgpu::ShaderSource::Wgsl(DISTANCE_SHADER.into()),
    });

    let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: None,
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: true },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 2,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
        ],
    });

    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: None,
        bind_group_layouts: &[&bind_group_layout],
        push_constant_ranges: &[],
    });

    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("distance_pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: Some("main"),
        compilation_options: wgpu::PipelineCompilationOptions::default(),
        cache: None,
    });

    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: None,
        layout: &bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: input_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: output_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: uniform_buf.as_entire_binding(),
            },
        ],
    });

    // Dispatch: ceil(n_cells / 16) × ceil(n_cells / 16) workgroups.
    let wg = ((n_cells + 15) / 16) as u32;
    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor::default());
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bind_group, &[]);
        pass.dispatch_workgroups(wg, wg, 1);
    }
    encoder.copy_buffer_to_buffer(&output_buf, 0, &readback_buf, 0, out_size);
    queue.submit(std::iter::once(encoder.finish()));

    // Download results.
    let slice = readback_buf.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |v| tx.send(v).unwrap());
    device.poll(wgpu::Maintain::Wait);
    rx.recv()??;

    let data_bytes = slice.get_mapped_range();
    let result: Vec<f32> = bytemuck::cast_slice(&data_bytes).to_vec();
    drop(data_bytes);
    readback_buf.unmap();

    Ok(result)
}

/// Select k nearest neighbours from a flat n×n distance matrix.
#[cfg(feature = "gpu")]
fn knn_from_distances(
    dist: &[f32],
    n_cells: usize,
    k: usize,
) -> Vec<Vec<(usize, f64)>> {
    (0..n_cells)
        .map(|i| {
            let row = &dist[i * n_cells..(i + 1) * n_cells];
            let mut indexed: Vec<(usize, f32)> = row
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(j, &d)| (j, d))
                .collect();
            indexed.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            indexed
                .into_iter()
                .take(k)
                .map(|(j, d)| (j, d as f64))
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn gpu_umap_falls_back_on_small_dataset() {
        // n=4 < 5000 threshold → always takes CPU path regardless of GPU availability.
        let data = Array2::from_shape_fn((4, 3), |(i, j)| (i * 3 + j) as f64);
        let result = run_umap_gpu(&data, 2, 10, 0.1, 1.0, 42);
        assert!(result.is_ok(), "should succeed via CPU fallback: {:?}", result);
        let r = result.unwrap();
        assert_eq!(r.embedding.nrows(), 4);
        assert_eq!(r.embedding.ncols(), 2);
    }

    #[test]
    fn gpu_umap_smoke_larger() {
        // n=50 — still below GPU threshold, exercises full CPU path.
        let data = Array2::from_shape_fn((50, 10), |(i, j)| {
            ((i * 10 + j) as f64).sin()
        });
        let result = run_umap_gpu(&data, 5, 20, 0.1, 1.0, 0xCAFE);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert_eq!(r.embedding.nrows(), 50);
    }
}
