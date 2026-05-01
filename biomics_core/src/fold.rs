use crossbeam_channel::Sender;
use rayon::prelude::*;
use std::time::Instant;

use crate::accum::BatchAccum;

/// Number of records processed per rayon task.
///
/// 64 K is the sweet spot: large enough to amortize rayon scheduling overhead,
/// small enough to keep L2 cache warm on a single core.
pub const BATCH_SIZE: usize = 64_000;

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

/// Execute a lock-free parallel fold over a pre-collected slice of records.
///
/// Two code paths are compiled:
/// - **With progress** (`progress_tx = Some(…)`): one `Sender::send` per chunk.
/// - **Without progress** (`progress_tx = None`): zero branch overhead per chunk.
///
/// Both paths share the same rayon `par_chunks` + `reduce` skeleton; the
/// `if let Some(tx)` branch is fully eliminated at compile time via monomorphism
/// because callers almost always pass a concrete `Option<&Sender<…>>` literal.
pub fn parallel_fold<A>(
    records: &[A::Record],
    modality: &'static str,
    progress_tx: Option<&Sender<ProgressEvent>>,
) -> anyhow::Result<A::Summary>
where
    A: BatchAccum,
{
    if let Some(tx) = progress_tx {
        parallel_fold_with_progress::<A>(records, modality, tx)
    } else {
        parallel_fold_bare::<A>(records)
    }
}

/// Inner fold without any progress reporting — zero branch cost per chunk.
#[inline]
fn parallel_fold_bare<A>(records: &[A::Record]) -> anyhow::Result<A::Summary>
where
    A: BatchAccum,
{
    let merged: A = records
        .par_chunks(BATCH_SIZE)
        .map(|chunk| {
            let mut acc = A::default();
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
            let processed = (((chunk_idx + 1) * BATCH_SIZE) as u64).min(total);
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
                    if count.is_multiple_of(BATCH_SIZE as u64) {
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
