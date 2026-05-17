use crossbeam_channel::Sender;
use rayon::prelude::*;
use std::time::Instant;

use crate::accum::BatchAccum;

/// Default chunk size: 64 KiB worth of records per rayon task.
///
/// The actual chunk size used at runtime is computed dynamically in
/// `parallel_fold` by targeting `TARGET_CHUNK_BYTES` (256 KB) of record data
/// per task.  This constant is the fallback when `size_of::<A::Record>() == 0`
/// or the computed value would otherwise be 0.
pub const DEFAULT_CHUNK: usize = 65_536;

/// Target approximately 256 KB of record data per rayon task so that each
/// chunk fits comfortably in the L2 cache of a modern core.
const TARGET_CHUNK_BYTES: usize = 256 * 1024;

/// Retained for backwards compatibility with call sites that reference it.
pub const BATCH_SIZE: usize = DEFAULT_CHUNK;

/// How many records ahead to prefetch when the `prefetch` feature is enabled.
pub const PREFETCH_DISTANCE: usize = 16;

/// Progress event sent over the crossbeam channel to the TUI or orchestration thread.
#[derive(Debug, Clone)]
pub struct ProgressEvent {
    pub modality: &'static str,
    pub records_processed: u64,
    pub total_records: u64,
    pub records_per_sec: f64,
}

impl ProgressEvent {
    #[inline]
    pub fn fraction(&self) -> f64 {
        if self.total_records == 0 {
            1.0
        } else {
            (self.records_processed as f64 / self.total_records as f64).min(1.0)
        }
    }
}

/// Compute the chunk size for a given record type, targeting `TARGET_CHUNK_BYTES`.
///
/// Falls back to `DEFAULT_CHUNK` for ZSTs or when the size would overflow.
#[inline]
fn chunk_size_for<T>() -> usize {
    let record_size = std::mem::size_of::<T>();
    if record_size == 0 {
        DEFAULT_CHUNK
    } else {
        (TARGET_CHUNK_BYTES / record_size).max(1)
    }
}

/// Execute a lock-free parallel fold over a pre-collected slice of records.
///
/// Two code paths are compiled:
/// - **With progress** (`progress_tx = Some(…)`): one `Sender::send` per chunk.
/// - **Without progress** (`progress_tx = None`): zero branch overhead per chunk.
///
/// Both paths share the same rayon `par_chunks` + `reduce` skeleton; the
/// `if let Some(tx)` branch is fully eliminated at compile time via monomorphism
/// because callers almost always pass a concrete `Option<&Sender<…>>` literal.
///
/// Chunk size is derived dynamically from `TARGET_CHUNK_BYTES` (256 KB) to keep
/// each task's working set inside the L2 cache.  Override by calling
/// `parallel_fold_with_chunk` directly.
pub fn parallel_fold<A>(
    records: &[A::Record],
    modality: &'static str,
    progress_tx: Option<&Sender<ProgressEvent>>,
) -> anyhow::Result<A::Summary>
where
    A: BatchAccum,
{
    let chunk_size = chunk_size_for::<A::Record>();
    if let Some(tx) = progress_tx {
        parallel_fold_with_progress::<A>(records, modality, tx, chunk_size)
    } else {
        parallel_fold_bare::<A>(records, chunk_size)
    }
}

/// Inner fold without any progress reporting — zero branch cost per chunk.
#[inline]
fn parallel_fold_bare<A>(records: &[A::Record], chunk_size: usize) -> anyhow::Result<A::Summary>
where
    A: BatchAccum,
{
    let merged: A = records
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut acc = A::default();

            #[cfg(feature = "prefetch")]
            {
                // Touch each record PREFETCH_DISTANCE ahead to warm the cache
                // line before the main loop reaches it.  `black_box` prevents
                // the compiler from eliding the load as dead code.
                for (i, record) in chunk.iter().enumerate() {
                    if let Some(ahead) = chunk.get(i + PREFETCH_DISTANCE) {
                        let _ = std::hint::black_box(ahead as *const _ as usize);
                    }
                    if let Err(e) = acc.process(record) {
                        log::warn!("skipping record: {}", e);
                    }
                }
            }

            #[cfg(not(feature = "prefetch"))]
            for record in chunk {
                if let Err(e) = acc.process(record) {
                    log::warn!("skipping record: {}", e);
                }
            }

            acc
        })
        .reduce(A::default, |mut left, right| {
            left.merge(right);
            left
        });
    merged.finalize()
}

/// Inner fold with per-chunk progress events.
#[inline]
fn parallel_fold_with_progress<A>(
    records: &[A::Record],
    modality: &'static str,
    tx: &Sender<ProgressEvent>,
    chunk_size: usize,
) -> anyhow::Result<A::Summary>
where
    A: BatchAccum,
{
    let total = records.len() as u64;
    let start = Instant::now();

    let merged: A = records
        .par_chunks(chunk_size)
        .enumerate()
        .map(|(chunk_idx, chunk)| {
            let mut acc = A::default();

            #[cfg(feature = "prefetch")]
            {
                for (i, record) in chunk.iter().enumerate() {
                    if let Some(ahead) = chunk.get(i + PREFETCH_DISTANCE) {
                        let _ = std::hint::black_box(ahead as *const _ as usize);
                    }
                    if let Err(e) = acc.process(record) {
                        log::warn!(
                            "[{}] skipping record in chunk {}: {}",
                            modality,
                            chunk_idx,
                            e
                        );
                    }
                }
            }

            #[cfg(not(feature = "prefetch"))]
            for record in chunk {
                if let Err(e) = acc.process(record) {
                    log::warn!(
                        "[{}] skipping record in chunk {}: {}",
                        modality,
                        chunk_idx,
                        e
                    );
                }
            }

            let processed = (((chunk_idx + 1) * chunk_size) as u64).min(total);
            let elapsed = start.elapsed().as_secs_f64().max(1e-9);
            let _ = tx.send(ProgressEvent {
                modality,
                records_processed: processed,
                total_records: total,
                records_per_sec: processed as f64 / elapsed,
            });
            acc
        })
        .reduce(A::default, |mut left, right| {
            left.merge(right);
            left
        });

    merged.finalize()
}

/// Streaming variant for sources too large to pre-collect.
///
/// Records arrive on `record_rx`; `thread_count` workers drain it in parallel,
/// each building a local accumulator. Workers merge at the end via a bounded
/// crossbeam channel — still lock-free during the processing phase.
pub fn streaming_fold<A>(
    record_rx: crossbeam_channel::Receiver<A::Record>,
    thread_count: usize,
    modality: &'static str,
    progress_tx: Option<Sender<ProgressEvent>>,
) -> anyhow::Result<A::Summary>
where
    A: BatchAccum + 'static,
    A::Record: Send + 'static,
{
    let (result_tx, result_rx) = crossbeam_channel::bounded::<A>(thread_count);

    rayon::scope(|s| {
        for _ in 0..thread_count {
            let rx = record_rx.clone();
            let tx = result_tx.clone();
            let ptx = progress_tx.clone();
            s.spawn(move |_| {
                let mut acc = A::default();
                let mut count = 0u64;
                let start = Instant::now();
                for record in rx.iter() {
                    if let Err(e) = acc.process(&record) {
                        log::warn!("[{}] skipping record: {}", modality, e);
                    }
                    count += 1;
                    if count.is_multiple_of(DEFAULT_CHUNK as u64) {
                        if let Some(ref p) = ptx {
                            let elapsed = start.elapsed().as_secs_f64().max(1e-9);
                            let _ = p.send(ProgressEvent {
                                modality,
                                records_processed: count,
                                total_records: 0,
                                records_per_sec: count as f64 / elapsed,
                            });
                        }
                    }
                }
                let _ = tx.send(acc);
            });
        }
        drop(result_tx);
    });

    let merged = result_rx
        .iter()
        .reduce(|mut a, b| {
            a.merge(b);
            a
        })
        .unwrap_or_default();

    merged.finalize()
}
