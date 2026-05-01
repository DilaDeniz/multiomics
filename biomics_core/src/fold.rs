use crossbeam_channel::Sender;
use rayon::prelude::*;
use std::time::Instant;

use crate::accum::BatchAccum;

/// Number of records processed per rayon task.
pub const BATCH_SIZE: usize = 50_000;

/// Progress event sent over the crossbeam channel to the TUI or orchestration thread.
///
/// Delivered approximately once per `BATCH_SIZE` records processed.
#[derive(Debug, Clone)]
pub struct ProgressEvent {
    /// Human-readable modality name, e.g. "genomics".
    pub modality: &'static str,
    /// Cumulative records processed so far (approximate in parallel context).
    pub records_processed: u64,
    /// Total records in the input (known at fold-start).
    pub total_records: u64,
    /// Instantaneous throughput estimate (records / second).
    pub records_per_sec: f64,
}

impl ProgressEvent {
    /// Returns completion fraction in [0.0, 1.0].
    pub fn fraction(&self) -> f64 {
        if self.total_records == 0 {
            1.0
        } else {
            (self.records_processed as f64 / self.total_records as f64).min(1.0)
        }
    }
}

/// Execute a lock-free parallel fold over a pre-collected slice of records.
///
/// ## Algorithm
/// 1. `par_chunks(BATCH_SIZE)` partitions the slice into independent chunks
///    that rayon distributes across the thread pool.
/// 2. Each chunk is folded locally into a fresh `A` (no shared state).
/// 3. Per-chunk accumulators are reduced pairwise via `A::merge` — still
///    no shared state.
/// 4. One final accumulator emerges; `finalize()` converts it to a `Summary`.
///
/// `progress_tx` receives one event per completed chunk. Pass `None` to
/// suppress progress reporting.
pub fn parallel_fold<A>(
    records: &[A::Record],
    modality: &'static str,
    progress_tx: Option<&Sender<ProgressEvent>>,
) -> anyhow::Result<A::Summary>
where
    A: BatchAccum,
{
    let total = records.len() as u64;
    let start = Instant::now();

    let merged: A = records
        .par_chunks(BATCH_SIZE)
        .enumerate()
        .map(|(chunk_idx, chunk)| {
            let mut acc = A::default();
            for record in chunk {
                if let Err(e) = acc.process(record) {
                    log::warn!("[{}] skipping record in chunk {}: {}", modality, chunk_idx, e);
                }
            }
            if let Some(tx) = progress_tx {
                let processed =
                    (((chunk_idx + 1) * BATCH_SIZE) as u64).min(total);
                let elapsed = start.elapsed().as_secs_f64().max(1e-9);
                let _ = tx.send(ProgressEvent {
                    modality,
                    records_processed: processed,
                    total_records: total,
                    records_per_sec: processed as f64 / elapsed,
                });
            }
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
/// Records arrive on `record_rx`; `thread_count` workers drain it, each building
/// a local accumulator, which are merged after all workers finish. The crossbeam
/// multi-consumer channel provides lock-free work distribution.
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
                    if count % BATCH_SIZE as u64 == 0 {
                        if let Some(ref p) = ptx {
                            let elapsed = start.elapsed().as_secs_f64().max(1e-9);
                            let _ = p.send(ProgressEvent {
                                modality,
                                records_processed: count,
                                total_records: 0, // unknown in streaming mode
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
