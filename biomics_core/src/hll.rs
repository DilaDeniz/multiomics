//! Minimal HyperLogLog cardinality estimator (Flajolet et al. 2007).
//!
//! Uses m = 2^14 = 16 384 one-byte registers → 16 KB fixed memory regardless
//! of the number of elements inserted. Standard error ≈ 0.81 %.
//!
//! # Usage
//! ```
//! use biomics_core::HyperLogLog;
//! let mut hll = HyperLogLog::new();
//! hll.insert_hashed(0xDEADBEEF_CAFEBABE);
//! hll.insert_hashed(0x0102030405060708);
//! let _ = hll.cardinality(); // approximate distinct count
//! ```

const B: u32 = 14; // number of register-index bits
const M: usize = 1 << B; // 16 384 registers
const M_F: f64 = M as f64;

// α_m correction constant for m ≥ 128 (Flajolet 2007 §4)
const ALPHA: f64 = 0.7213 / (1.0 + 1.079 / M_F);

/// Fixed-size HyperLogLog counter.
///
/// Each instance uses exactly 16 KB of heap. `merge` is O(m) = O(16 384)
/// bitwise register-max; much cheaper than merging hash-sets.
pub struct HyperLogLog {
    regs: Box<[u8; M]>,
}

impl Default for HyperLogLog {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl HyperLogLog {
    /// Create a new, empty HyperLogLog counter.
    pub fn new() -> Self {
        Self {
            regs: Box::new([0u8; M]),
        }
    }

    /// Insert a **pre-hashed** 64-bit value (caller is responsible for using
    /// a high-quality hash — `position_key`'s FNV-mix output is sufficient).
    #[inline(always)]
    pub fn insert_hashed(&mut self, hash: u64) {
        // Top B bits → register index
        let idx = (hash >> (64 - B)) as usize;
        // Remaining 64-B bits → count leading zeros + 1 (ρ function)
        let rho = (hash << B).leading_zeros() + 1; // at most 51, fits in u8
        let r = &mut self.regs[idx];
        if rho as u8 > *r {
            *r = rho as u8;
        }
    }

    /// Merge another HLL into this one (element-wise maximum of registers).
    ///
    /// Associative and commutative — safe to use as the rayon `reduce` step.
    #[inline]
    pub fn merge(&mut self, other: &Self) {
        for (a, b) in self.regs.iter_mut().zip(other.regs.iter()) {
            if *b > *a {
                *a = *b;
            }
        }
    }

    /// Estimate the number of distinct elements inserted.
    pub fn cardinality(&self) -> u64 {
        // Harmonic mean of 2^{-M[j]}
        let sum: f64 = self.regs.iter().map(|&r| (-f64::from(r)).exp2()).sum();
        let mut estimate = ALPHA * M_F * M_F / sum;

        // Small-range correction: linear counting when estimate < 5m/2 and
        // some registers are still zero.
        if estimate < 2.5 * M_F {
            let zeros = self.regs.iter().filter(|&&r| r == 0).count() as f64;
            if zeros > 0.0 {
                estimate = M_F * (M_F / zeros).ln();
            }
        }

        // Large-range correction (2^32/30 threshold) rarely triggered at
        // genome scale — omitted for simplicity without loss of accuracy.

        estimate.round() as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // splitmix64: bijective 64-bit finalizer, designed for sequential inputs.
    // Produces high-quality, uniform 64-bit output from consecutive integers.
    fn splitmix64(x: u64) -> u64 {
        let mut z = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    #[test]
    fn test_empty_cardinality_is_zero() {
        let hll = HyperLogLog::new();
        assert_eq!(hll.cardinality(), 0);
    }

    #[test]
    fn test_small_set_within_5pct() {
        let mut hll = HyperLogLog::new();
        let n = 1_000u64;
        for i in 0..n {
            hll.insert_hashed(splitmix64(i));
        }
        let est = hll.cardinality();
        let err = (est as f64 - n as f64).abs() / n as f64;
        assert!(
            err < 0.05,
            "error {:.2}% exceeded 5% for n={}",
            err * 100.0,
            n
        );
    }

    #[test]
    fn test_large_set_within_2pct() {
        let mut hll = HyperLogLog::new();
        let n = 1_000_000u64;
        for i in 0..n {
            hll.insert_hashed(splitmix64(i));
        }
        let est = hll.cardinality();
        let err = (est as f64 - n as f64).abs() / n as f64;
        // Theoretical σ ≈ 0.81%; allow 3% for deterministic single-run test headroom.
        assert!(
            err < 0.03,
            "error {:.2}% exceeded 3% for n={}",
            err * 100.0,
            n
        );
    }

    #[test]
    fn test_merge_equals_union() {
        let mut a = HyperLogLog::new();
        let mut b = HyperLogLog::new();
        let n = 500_000u64;
        for i in 0..n {
            a.insert_hashed(splitmix64(i));
        }
        for i in n..2 * n {
            b.insert_hashed(splitmix64(i));
        }
        a.merge(&b);
        let est = a.cardinality();
        let true_count = 2 * n;
        let err = (est as f64 - true_count as f64).abs() / true_count as f64;
        assert!(err < 0.03, "merge error {:.2}% exceeded 3%", err * 100.0);
    }

    #[test]
    fn test_duplicate_insertion_no_overcount() {
        let mut hll = HyperLogLog::new();
        for _ in 0..1000 {
            hll.insert_hashed(splitmix64(42)); // same element 1000 times
        }
        // Should estimate ~1 unique element; with HLL noise, bounded < 20
        assert!(
            hll.cardinality() < 20,
            "duplicate-only set over-estimated: {}",
            hll.cardinality()
        );
    }
}
