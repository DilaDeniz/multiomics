//! GPU-accelerated UMAP via WebGPU (wgpu) compute shaders.
//!
//! Requires the `gpu` feature flag. Falls back to CPU transparently when no
//! compatible adapter is found or the dataset is small (n < 5 000).
//!
//! # Algorithm: GPU k-selection
//! Instead of computing the full n×n distance matrix and copying it to CPU,
//! a single compute shader does **per-cell k-NN selection entirely on the GPU**:
//!
//! - One GPU thread handles one cell.
//! - Each thread scans all n cells, computing distances, and maintains a
//!   private max-heap of the k nearest neighbours found so far.
//! - After the scan the heap is sorted in-place and written to the output buffer.
//!
//! Data transferred = n × k × 8 bytes (indices + distances).
//! For n = 100 k cells, k = 15: **12 MB** — vs 40 GB for the naïve n×n approach.
//! No tile loop, no readback stalls, single dispatch.
//!
//! # Platform support
//! wgpu targets Vulkan (Linux/Windows), Metal (macOS/iOS), DX12 (Windows).
//! Enable with `--features gpu`.

use anyhow::Result;
use ndarray::Array2;

use crate::umap::{run_umap, UmapResult};

/// Maximum k supported by the GPU shader (private heap size).
/// Covers all practical UMAP neighbour counts (default 15, max recommended 50).
#[cfg(feature = "gpu")]
const MAX_K: usize = 64;

/// WGSL compute shader: per-cell k-nearest-neighbour selection.
///
/// One invocation per cell. Each thread keeps a private max-heap of size k,
/// scanning all n cells to find the k nearest. No global distance matrix
/// is written — output is n×k indices and n×k distances only.
#[cfg(feature = "gpu")]
const KNN_SELECT_SHADER: &str = r#"
const MAX_K: u32 = 64u;

struct Uniforms {
    n_cells: u32,
    n_dims:  u32,
    k:       u32,
    _pad:    u32,
}

@group(0) @binding(0) var<storage, read>       data:        array<f32>;
@group(0) @binding(1) var<storage, read_write> out_indices: array<u32>;
@group(0) @binding(2) var<storage, read_write> out_dists:   array<f32>;
@group(0) @binding(3) var<uniform>             uniforms:    Uniforms;

// Per-invocation (thread-private) heap — no shared memory needed.
var<private> hp_dist: array<f32, 64>;
var<private> hp_idx:  array<u32, 64>;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let n = uniforms.n_cells;
    let d = uniforms.n_dims;
    let k = min(uniforms.k, MAX_K);

    if i >= n { return; }

    // Initialise heap slots to +inf so any real distance wins immediately.
    for (var h = 0u; h < MAX_K; h++) {
        hp_dist[h] = 1e38f;
        hp_idx[h]  = 0u;
    }
    var heap_size: u32 = 0u;

    // Scan every other cell, maintain top-k by max-heap.
    for (var j = 0u; j < n; j++) {
        if j == i { continue; }

        var dsq: f32 = 0.0;
        for (var dim = 0u; dim < d; dim++) {
            let diff = data[i * d + dim] - data[j * d + dim];
            dsq += diff * diff;
        }
        let dist = sqrt(dsq);

        if heap_size < k {
            hp_dist[heap_size] = dist;
            hp_idx[heap_size]  = j;
            heap_size++;
        } else {
            // Find the current worst (furthest) neighbour.
            var worst = 0u;
            for (var h = 1u; h < k; h++) {
                if hp_dist[h] > hp_dist[worst] { worst = h; }
            }
            // Replace it if this cell is closer.
            if dist < hp_dist[worst] {
                hp_dist[worst] = dist;
                hp_idx[worst]  = j;
            }
        }
    }

    // Insertion-sort the heap ascending by distance before writing.
    for (var a = 1u; a < k; a++) {
        let kd = hp_dist[a];
        let ki = hp_idx[a];
        var b  = a;
        while b > 0u && hp_dist[b - 1u] > kd {
            hp_dist[b] = hp_dist[b - 1u];
            hp_idx[b]  = hp_idx[b - 1u];
            b--;
        }
        hp_dist[b] = kd;
        hp_idx[b]  = ki;
    }

    // Write k nearest neighbours for cell i.
    let base = i * k;
    for (var h = 0u; h < k; h++) {
        out_indices[base + h] = hp_idx[h];
        out_dists[base + h]   = hp_dist[h];
    }
}
"#;

// ── Public entry point ────────────────────────────────────────────────────────

/// Run UMAP with GPU-accelerated k-NN when the `gpu` feature is enabled.
///
/// Falls back to [`run_umap`] (CPU) when:
/// - The `gpu` feature is not compiled in.
/// - No compatible GPU adapter is found.
/// - n < 5 000 (GPU launch overhead exceeds the compute saving).
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

// ── GPU implementation ────────────────────────────────────────────────────────

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
    let k = n_neighbors.min(MAX_K).min(n_cells.saturating_sub(1));

    let data_f32: Vec<f32> = data.iter().map(|&x| x as f32).collect();

    log::info!("GPU UMAP: {n_cells} cells × {n_dims} dims, k={k} (k-selection shader)");

    let knn = pollster::block_on(gpu_knn_exact(&data_f32, n_cells, n_dims, k))?;

    let (adjacency, rho) = compute_fuzzy_graph_from_knn(&knn, n_cells);
    run_umap_from_graph(&adjacency, &rho, n_cells, n_epochs, min_dist, learning_rate, seed)
}

/// GPU exact k-NN via per-cell k-selection shader.
///
/// Single dispatch — no tile loop. Output: n×k indices + n×k distances.
/// Transfer cost: n × k × 8 bytes regardless of n.
#[cfg(feature = "gpu")]
async fn gpu_knn_exact(
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

    log::info!(
        "GPU: {} ({})",
        adapter.get_info().name,
        format!("{:?}", adapter.get_info().backend)
    );

    // ── Buffers ───────────────────────────────────────────────────────────────

    let input_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("umap_data"),
        contents: bytemuck::cast_slice(data),
        usage: wgpu::BufferUsages::STORAGE,
    });

    // Output: n × k u32 indices  +  n × k f32 distances
    let out_bytes = (n_cells * k * std::mem::size_of::<f32>()) as u64;

    let idx_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_knn_idx"),
        size: out_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let dist_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_knn_dist"),
        size: out_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let idx_rb = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_knn_idx_rb"),
        size: out_bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });
    let dist_rb = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("umap_knn_dist_rb"),
        size: out_bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    // Uniform: [n_cells, n_dims, k, pad] — 16-byte aligned
    let uniforms: [u32; 4] = [n_cells as u32, n_dims as u32, k as u32, 0];
    let uniform_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("umap_uniforms"),
        contents: bytemuck::cast_slice(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    // ── Pipeline ──────────────────────────────────────────────────────────────

    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("umap_kselect"),
        source: wgpu::ShaderSource::Wgsl(KNN_SELECT_SHADER.into()),
    });

    let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: None,
        entries: &[
            storage_entry(0, true),
            storage_entry(1, false),
            storage_entry(2, false),
            uniform_entry(3),
        ],
    });

    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("umap_kselect_pipeline"),
        layout: Some(&device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[Some(&bgl)],
            immediate_size: 0,
        })),
        module: &shader,
        entry_point: Some("main"),
        compilation_options: wgpu::PipelineCompilationOptions::default(),
        cache: None,
    });

    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: None,
        layout: &bgl,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: input_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: idx_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: dist_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: uniform_buf.as_entire_binding() },
        ],
    });

    // ── Dispatch (one thread per cell) ────────────────────────────────────────

    let wg_count = ((n_cells as u32) + 255) / 256;
    let mut encoder =
        device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
    {
        let mut pass =
            encoder.begin_compute_pass(&wgpu::ComputePassDescriptor::default());
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bind_group, &[]);
        pass.dispatch_workgroups(wg_count, 1, 1);
    }
    encoder.copy_buffer_to_buffer(&idx_buf, 0, &idx_rb, 0, out_bytes);
    encoder.copy_buffer_to_buffer(&dist_buf, 0, &dist_rb, 0, out_bytes);
    queue.submit(std::iter::once(encoder.finish()));

    // ── Single tiny readback (n × k × 8 bytes) ───────────────────────────────

    let idx_slice = idx_rb.slice(..);
    let dist_slice = dist_rb.slice(..);
    let (tx_i, rx_i) = std::sync::mpsc::channel();
    let (tx_d, rx_d) = std::sync::mpsc::channel();
    idx_slice.map_async(wgpu::MapMode::Read, move |v| tx_i.send(v).unwrap());
    dist_slice.map_async(wgpu::MapMode::Read, move |v| tx_d.send(v).unwrap());
    device.poll(wgpu::PollType::Wait { submission_index: None, timeout: None })?;
    rx_i.recv()??;
    rx_d.recv()??;

    let indices: Vec<u32> = bytemuck::cast_slice(&idx_slice.get_mapped_range()).to_vec();
    let dists: Vec<f32> = bytemuck::cast_slice(&dist_slice.get_mapped_range()).to_vec();
    idx_rb.unmap();
    dist_rb.unmap();

    // Convert to the format expected by compute_fuzzy_graph_from_knn.
    let knn = (0..n_cells)
        .map(|i| {
            let base = i * k;
            (0..k)
                .map(|ki| (indices[base + ki] as usize, dists[base + ki] as f64))
                .collect()
        })
        .collect();

    Ok(knn)
}

// ── Bind group layout helpers ─────────────────────────────────────────────────

#[cfg(feature = "gpu")]
fn storage_entry(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

#[cfg(feature = "gpu")]
fn uniform_entry(binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn gpu_umap_falls_back_on_small_dataset() {
        let data = Array2::from_shape_fn((4, 3), |(i, j)| (i * 3 + j) as f64);
        let result = run_umap_gpu(&data, 2, 10, 0.1, 1.0, 42);
        assert!(result.is_ok(), "should succeed via CPU fallback: {:?}", result);
        let r = result.unwrap();
        assert_eq!(r.embedding.nrows(), 4);
        assert_eq!(r.embedding.ncols(), 2);
    }

    #[test]
    fn gpu_umap_smoke_larger() {
        let data = Array2::from_shape_fn((50, 10), |(i, j)| ((i * 10 + j) as f64).sin());
        let result = run_umap_gpu(&data, 5, 20, 0.1, 1.0, 0xCAFE);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().embedding.nrows(), 50);
    }
}
