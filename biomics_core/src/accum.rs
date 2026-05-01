/// A modality-specific accumulator that collects statistics over a batch of records
/// and can be merged with another accumulator in a lock-free parallel fold.
///
/// Implementors must be `Send + Default` so rayon can distribute fresh accumulators
/// to worker threads and fold them together without shared mutable state.
pub trait BatchAccum: Send + Default {
    /// The raw record type fed into this accumulator.
    type Record: Send + Sync;

    /// The finalized, serializable result produced after all merges.
    type Summary: Send + serde::Serialize;

    /// Process a single record, updating internal counters in-place.
    ///
    /// Returning `Err` causes the record to be skipped with a warning logged;
    /// it does not abort the analysis.
    fn process(&mut self, record: &Self::Record) -> anyhow::Result<()>;

    /// Merge `other` into `self`, consuming `other`.
    ///
    /// Must be commutative and associative — rayon may call this in any order.
    fn merge(&mut self, other: Self);

    /// Finalize internal state into a summary, consuming `self`.
    ///
    /// Called exactly once after the final merge.
    fn finalize(self) -> anyhow::Result<Self::Summary>;
}
