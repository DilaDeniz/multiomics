//! GPU-accelerated UMAP via WebGPU (wgpu) compute shaders.
//!
//! Requires the `gpu` feature flag. When wgpu cannot initialise a device
//! (no compatible GPU, wrong platform, CI environment) the function
//! transparently falls back to the CPU UMAP implementation in `umap.rs`.
//!
//! # GPU speedup strategy
//! The dominant cost in UMAP Phase 1 for large datasets (n > 10 000) is the
//! O(n² × d) pairwise distance matrix. A WebGPU compute shader computes all
//! distances in parallel using a 16×16 workgroup tile, reducing wall time
//! from minutes (CPU) to seconds (GPU) for 100 k+ cells.
//!
//! # VRAM-safe tiled KNN
//! Storing the full n×n float32 distance matrix on the GPU requires n²×4 bytes
//! of VRAM — 40 GB for 100 k cells, far exceeding any consumer GPU. Instead the
//! kernel is dispatched in **row-tiles**: each tile computes distances for
//! `tile_rows` cells against all n cells, producing a `tile_rows × n` matrix
//! that fits comfortably in the VRAM budget. The k-NN for those rows is
//! extracted and the tile is discarded before the next dispatch. Peak VRAM
//! usage is O(tile_rows × n) ≈ 1.5 GB regardless of total n.
//!
//! For an RTX 4050 (6 GB VRAM) with n = 100 k cells this means ~67 tiles of
//! 1 500 rows each — still orders of magnitude faster than the CPU O(n²) scan.
//!
//! # Phase 2
//! SGD embedding optimisation runs on CPU using the GPU-computed KNN graph.
//! The SGD step is memory-bandwidth bound and harder to parallelise without
//! atomic contention.
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

/// WGSL compute shader: tiled pairwise Euclidean distances.
///
/// Computes distances from a tile of rows (`tile_offset..tile_offset+tile_rows`)
/// to ALL n cells.  Output is row-major `[tile_rows × n_cells]`.
///
/// Uniforms layout: `{ n_cells: u32, n_dims: u32, tile_offset: u32, tile_rows: u32 }`
#[cfg(feature = "gpu")]
const DISTANCE_SHADER_TILED: &str = r#"
struct Uniforms {
    n_cells:     u32,
    n_dims:      u32,
    tile_offset: u32,
    tile_rows:   u32,
}

@group(0) @binding(0) var<storage, read>       data:       array<f32>;
@group(0) @binding(1) var<storage, read_write> distances:  array<f32>;
@group(0) @binding(2) var<uniform>             uniforms:   Uniforms;

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let local_i = gid.x;   // row within this tile  (0 .. tile_rows)
    let j       = gid.y;   // column — all cells     (0 .. n_cells)
    let n = uniforms.n_cells;
    let d = uniforms.n_dims;

    if local_i >= uniforms.tile_rows || j >= n {
        return;
    }

    let global_i = uniforms.tile_offset + local_i;

    var dist_sq: f32 = 0.0;
    for (var k: u32 = 0u; k < d; k = k + 1u) {
        let a = data[global_i * d + k];
        let b = data[j * d + k];
        let diff = a - b;
        dist_sq = dist_sq + diff * diff;
    }

    distances[local_i * n + j] = sqrt(dist_sq);
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
    #[cfg(feature = "gpu")]
    if data.nrows() >= 5_000 {
        let n_cells = data.nrows();
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

    let data_f32: Vec<f32> = data.iter().map(|&x| x as f32).collect();

    log::info!("GPU UMAP: uploading {n_cells}×{n_dims} matrix, computing tiled KNN…");

    let knn = pollster::block_on(gpu_knn_tiled(&data_f32, n_cells, n_dims, n_neighbors))?;

    let (adjacency, rho) = compute_fuzzy_graph_from_knn(&knn, n_cells);
    run_umap_from_graph(&adjacency, &rho, n_cells, n_epochs, min_dist, learning_rate, seed)
}

/// GPU tiled KNN: dispatches the distance shader in row-tiles to stay within
/// VRAM budget, extracts k-NN per tile, and returns the full KNN graph.
///
/// Peak VRAM = O(tile_rows × n_cells) ≈ 1.5 GB regardless of n_cells.
#[cfg(feature = "gpu")]
async fn gpu_knn_tiled(
    data: &[f32],
    n_cells: usize,
    n_dims: usize,
    k: usize,
) -> Result<Vec<Vec<(usize, f64)>>> {
    use wgpu::util::DeviceExt;

    let instance = wgpu::Instance::default();
    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        })
        .await
        .map_err(|_| anyhow::anyhow!("no wgpu adapter found"))?;

    let (device, queue) = adapter
        .request_device(&wgpu::DeviceDescriptor::default())
        .await?;

    // Two separate GPU limits apply:
    // - max_buffer_size: how large a buffer can be allocated
    // - max_storage_buffer_binding_size: how large a storage binding range can be
    // The tile must fit within both. Use 90% of the smaller one for safety.
    let limits = device.limits();
    let max_tile_bytes = (limits.max_buffer_size as usize)
        .min(limits.max_storage_buffer_binding_size as usize)
        * 9 / 10;
    let tile_rows = (max_tile_bytes / (n_cells * std::mem::size_of::<f32>()))
        .max(1)
        .min(n_cells);
    let n_tiles = (n_cells + tile_rows - 1) / tile_rows;

    log::info!(
        "GPU UMAP: tile_rows={tile_rows} ({n_tiles} tiles, \
         {:.0} MB/tile, k={k})",
        (tile_rows * n_cells * 4) as f64 / 1e6
    );

    // ── Static GPU resources (created once, reused across all tiles) ─────────

    let input_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("umap_input"),
        contents: bytemuck::cast_slice(data),
        usage: wgpu::BufferUsages::STORAGE,
    });

    // Output + readback buffers sized for the largest tile (tile_rows rows).
    let tile_buf_bytes = (tile_rows * n_cells * std::mem::size_of::<f32>()) as u64;
    let output_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_tile_out"),
        size: tile_buf_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let readback_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_readback"),
        size: tile_buf_bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    // Uniform buffer: [n_cells, n_dims, tile_offset, tile_rows] — updated via
    // queue.write_buffer() each tile so no reallocation is needed.
    let uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_uniforms"),
        size: 16,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("umap_distance_tiled"),
        source: wgpu::ShaderSource::Wgsl(DISTANCE_SHADER_TILED.into()),
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

    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("umap_distance_pipeline"),
        layout: Some(&device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        })),
        module: &shader,
        entry_point: Some("main"),
        compilation_options: wgpu::PipelineCompilationOptions::default(),
        cache: None,
    });

    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: None,
        layout: &bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: input_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: output_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: uniform_buf.as_entire_binding() },
        ],
    });

    // ── Tile loop ─────────────────────────────────────────────────────────────

    let mut knn: Vec<Vec<(usize, f64)>> = vec![vec![]; n_cells];

    for tile_idx in 0..n_tiles {
        let row_offset = tile_idx * tile_rows;
        let this_tile = tile_rows.min(n_cells - row_offset);

        log::debug!(
            "GPU UMAP: tile {}/{} rows {}..{}",
            tile_idx + 1,
            n_tiles,
            row_offset,
            row_offset + this_tile
        );

        // Update tile-specific uniforms (16 bytes, no reallocation).
        let uniforms: [u32; 4] = [
            n_cells as u32,
            n_dims as u32,
            row_offset as u32,
            this_tile as u32,
        ];
        queue.write_buffer(&uniform_buf, 0, bytemuck::cast_slice(&uniforms));

        // Dispatch: ceil(this_tile/16) × ceil(n_cells/16) workgroups.
        let wg_x = ((this_tile + 15) / 16) as u32;
        let wg_y = ((n_cells + 15) / 16) as u32;
        let actual_bytes = (this_tile * n_cells * std::mem::size_of::<f32>()) as u64;

        let mut encoder =
            device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
        {
            let mut pass =
                encoder.begin_compute_pass(&wgpu::ComputePassDescriptor::default());
            pass.set_pipeline(&pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(wg_x, wg_y, 1);
        }
        encoder.copy_buffer_to_buffer(&output_buf, 0, &readback_buf, 0, actual_bytes);
        queue.submit(std::iter::once(encoder.finish()));

        // Map readback buffer and extract k-NN for this tile.
        let slice = readback_buf.slice(..actual_bytes);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |v| tx.send(v).unwrap());
        device.poll(wgpu::PollType::Wait { submission_index: None, timeout: None })?;
        rx.recv()??;

        {
            let mapped = slice.get_mapped_range();
            let dists: &[f32] = bytemuck::cast_slice(&mapped);

            for local_i in 0..this_tile {
                let global_i = row_offset + local_i;
                let row = &dists[local_i * n_cells..(local_i + 1) * n_cells];

                let mut indexed: Vec<(usize, f32)> = row
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| j != global_i)
                    .map(|(j, &d)| (j, d))
                    .collect();
                indexed.sort_unstable_by(|a, b| {
                    a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
                });
                knn[global_i] = indexed
                    .into_iter()
                    .take(k)
                    .map(|(j, d)| (j, d as f64))
                    .collect();
            }
        }
        readback_buf.unmap();
    }

    Ok(knn)
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
        let data = Array2::from_shape_fn((50, 10), |(i, j)| ((i * 10 + j) as f64).sin());
        let result = run_umap_gpu(&data, 5, 20, 0.1, 1.0, 0xCAFE);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert_eq!(r.embedding.nrows(), 50);
    }
}
