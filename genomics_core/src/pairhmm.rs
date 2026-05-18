//! Pair-HMM forward algorithm for P(read | haplotype).
//!
//! Computes the log-probability of a sequenced read given a candidate
//! haplotype, accounting for base substitution errors (from Phred scores),
//! insertion errors, and deletion errors via affine gap penalties.
//!
//! Three states: M (match/mismatch), I (insertion in read), D (deletion in read).
//!
//! # References
//! * Durbin et al. (1998) "Biological Sequence Analysis." Cambridge UP.
//! * Li & Durbin (2009) "Fast and accurate short read alignment with
//!   Burrows-Wheeler Aligner." Bioinformatics 25(14):1754–1760.
//! * McKenna et al. (2010) GATK. Genome Research 20:1297–1303.

// ── Model constants ───────────────────────────────────────────────────────────

/// Probability of opening a gap (either insertion or deletion).
const GAP_OPEN: f64 = 1e-4;
/// Probability of extending an open gap by one more base.
const GAP_EXTEND: f64 = 0.1;

// Precomputed log transition probabilities.
const LN_MM: f64 = {
    // ln(1 - 2 * GAP_OPEN) — computed at compile time is tricky; use lazy below.
    // We just store the constant values as statics.
    0.0 // placeholder — replaced by lazy statics below
};

use std::sync::LazyLock;

/// Precomputed transition log-probabilities.
static LN_TRANS: LazyLock<[f64; 6]> = LazyLock::new(|| {
    [
        (1.0 - 2.0 * GAP_OPEN).ln(), // M → M
        GAP_OPEN.ln(),                // M → I
        GAP_OPEN.ln(),                // M → D  (== M → I by symmetry)
        (1.0 - GAP_EXTEND).ln(),      // I → M  (also D → M)
        GAP_EXTEND.ln(),              // I → I  (also D → D)
        0.25_f64.ln(),                // insertion emission (uniform)
    ]
});

// Index aliases into LN_TRANS.
const LN_MM_IDX: usize = 0;
const LN_MI_IDX: usize = 1;
const LN_MD_IDX: usize = 2;
const LN_XM_IDX: usize = 3; // I→M or D→M
const LN_XX_IDX: usize = 4; // I→I or D→D
const LN_INS_EMIT: usize = 5;

// ── Public API ────────────────────────────────────────────────────────────────

/// Compute `ln P(read | haplotype)` using the pair-HMM forward algorithm.
///
/// Returns the log-probability in natural-log space.  The DP tables are
/// allocated as flat `Vec<f64>` with shape `[state][hap_len+1][read_len+1]`
/// stored in row-major order.
///
/// # Arguments
/// * `read_seq`  – Nucleotide bases of the read (ASCII uppercase).
/// * `read_qual` – Phred-scaled base quality scores (same length as `read_seq`).
/// * `hap`       – Candidate haplotype sequence.
pub fn pair_hmm_log_prob(read_seq: &[u8], read_qual: &[u8], hap: &[u8]) -> f64 {
    let t = &*LN_TRANS;
    let _ = LN_MM; // suppress unused warning on the placeholder

    let n = read_seq.len(); // read length
    let m = hap.len();      // haplotype length

    if n == 0 || m == 0 {
        return f64::NEG_INFINITY;
    }

    let sz = (m + 1) * (n + 1);
    // Three DP layers: M, I, D.
    let mut fm = vec![f64::NEG_INFINITY; sz];
    let mut fi = vec![f64::NEG_INFINITY; sz];
    let mut fd = vec![f64::NEG_INFINITY; sz];

    let idx = |i: usize, j: usize| i * (n + 1) + j;

    // Initialise: start in M at (0,0) with prob ≈ 1.
    fm[idx(0, 0)] = t[LN_MM_IDX]; // accounts for first M→M self transition
    // Allow leading deletions (consume haplotype without read bases).
    fd[idx(1, 0)] = t[LN_MD_IDX] + t[LN_XX_IDX];
    for i in 2..=m {
        let prev = fd[idx(i - 1, 0)];
        if prev.is_finite() {
            fd[idx(i, 0)] = prev + t[LN_XX_IDX];
        }
    }

    // Fill DP.
    for i in 1..=m {
        for j in 1..=n {
            // ── M state ────────────────────────────────────────────────────
            let emit_m = ln_emit_match(read_seq[j - 1], read_qual[j - 1], hap[i - 1]);
            let from_m_m = add_ln(fm[idx(i - 1, j - 1)], t[LN_MM_IDX]);
            let from_i_m = add_ln(fi[idx(i - 1, j - 1)], t[LN_XM_IDX]);
            let from_d_m = add_ln(fd[idx(i - 1, j - 1)], t[LN_XM_IDX]);
            fm[idx(i, j)] = emit_m + log_sum_exp3(from_m_m, from_i_m, from_d_m);

            // ── I state (insertion in read: read base, no hap advance) ─────
            let from_m_i = add_ln(fm[idx(i, j - 1)], t[LN_MI_IDX]);
            let from_i_i = add_ln(fi[idx(i, j - 1)], t[LN_XX_IDX]);
            fi[idx(i, j)] = t[LN_INS_EMIT] + log_sum_exp(from_m_i, from_i_i);

            // ── D state (deletion in read: hap base consumed, no read adv) ─
            let from_m_d = add_ln(fm[idx(i - 1, j)], t[LN_MD_IDX]);
            let from_d_d = add_ln(fd[idx(i - 1, j)], t[LN_XX_IDX]);
            fd[idx(i, j)] = log_sum_exp(from_m_d, from_d_d);
        }
    }

    // Sum M and I over all haplotype positions at the final read position j=n.
    let mut total = f64::NEG_INFINITY;
    for i in 0..=m {
        total = log_sum_exp(total, fm[idx(i, n)]);
        total = log_sum_exp(total, fi[idx(i, n)]);
    }
    total
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Log-probability of emitting `read_base` in state M given `hap_base` and
/// Phred quality `qual`.
///
/// `ε = 10^(-qual/10)`.  If bases match: `ln(1 - ε)`; else: `ln(ε / 3)`.
#[inline]
fn ln_emit_match(read_base: u8, qual: u8, hap_base: u8) -> f64 {
    let eps = 10.0_f64.powf(-(qual as f64) / 10.0);
    if read_base.eq_ignore_ascii_case(&hap_base) {
        (1.0 - eps).max(f64::MIN_POSITIVE).ln()
    } else {
        (eps / 3.0).max(f64::MIN_POSITIVE).ln()
    }
}

/// `ln(exp(a) + exp(b))` computed in a numerically stable way.
///
/// If both inputs are `-∞` the result is `-∞`.
#[inline]
pub fn log_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY && b == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (big, small) = if a >= b { (a, b) } else { (b, a) };
    big + (1.0 + (small - big).exp()).ln()
}

/// Three-way log-sum-exp.
#[inline]
fn log_sum_exp3(a: f64, b: f64, c: f64) -> f64 {
    log_sum_exp(log_sum_exp(a, b), c)
}

/// Add two log-space values; if either is `-∞` return `-∞`.
#[inline]
fn add_ln(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY || b == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    } else {
        a + b
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_sum_exp_basic() {
        // ln(e^0 + e^0) = ln(2)
        let result = log_sum_exp(0.0, 0.0);
        assert!((result - 2.0_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn log_sum_exp_neginf() {
        assert_eq!(log_sum_exp(f64::NEG_INFINITY, 0.0), 0.0);
        assert_eq!(log_sum_exp(0.0, f64::NEG_INFINITY), 0.0);
        assert_eq!(
            log_sum_exp(f64::NEG_INFINITY, f64::NEG_INFINITY),
            f64::NEG_INFINITY
        );
    }

    #[test]
    fn pair_hmm_exact_read_higher_than_mismatch() {
        // A perfect-match read should have higher probability than one with
        // mismatches.
        let hap = b"ACGTACGT";
        let qual = vec![30u8; 8];
        let perfect = b"ACGTACGT";
        let mismatch = b"ACGTTCGT"; // one mismatch

        let lp_perfect = pair_hmm_log_prob(perfect, &qual, hap);
        let lp_mismatch = pair_hmm_log_prob(mismatch, &qual, hap);

        assert!(
            lp_perfect > lp_mismatch,
            "perfect match ({lp_perfect}) should exceed mismatch ({lp_mismatch})"
        );
    }

    #[test]
    fn pair_hmm_empty_inputs() {
        assert_eq!(pair_hmm_log_prob(b"", b"", b"ACGT"), f64::NEG_INFINITY);
        assert_eq!(pair_hmm_log_prob(b"ACGT", &[30u8; 4], b""), f64::NEG_INFINITY);
    }

    #[test]
    fn pair_hmm_returns_finite_for_valid_input() {
        let hap = b"ACGTACGTACGT";
        let read = b"ACGTACGT";
        let qual = vec![30u8; 8];
        let lp = pair_hmm_log_prob(read, &qual, hap);
        assert!(lp.is_finite(), "log-prob should be finite: {lp}");
        assert!(lp < 0.0, "log-prob should be negative: {lp}");
    }
}
