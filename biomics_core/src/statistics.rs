// Rigorous statistical functions used across all modality crates.
//
// Implements:
// - Lanczos ln Γ — underpins Fisher's exact test
// - Hypergeometric one-sided p-value — Fisher's exact test for pathway enrichment
// - Benjamini-Hochberg FDR correction
// - Regularised incomplete beta — underpins t-test p-values
// - Welch's two-sample t-test with Satterthwaite df
// - Spearman rank correlation

// ── Gamma function ────────────────────────────────────────────────────────────

/// Natural log of the gamma function via the Lanczos g=7 approximation.
/// Accurate to ~15 significant figures for Re(x) > 0.
pub fn ln_gamma(x: f64) -> f64 {
    // Reflection formula for x < 0.5
    if x < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - ln_gamma(1.0 - x);
    }
    const G: f64 = 7.0;
    #[allow(clippy::excessive_precision)]
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let xm1 = x - 1.0;
    // Accumulate the Lanczos series using fused multiply-add where possible.
    // Each term is ci / (xm1 + i + 1); we keep a running sum.
    let mut a = C[0];
    for (i, &ci) in C[1..].iter().enumerate() {
        a = ci.mul_add((xm1 + i as f64 + 1.0).recip(), a);
    }
    let t = xm1 + G + 0.5;
    // Use mul_add for the final linear combination to reduce rounding error.
    let log_prefix = (2.0 * std::f64::consts::PI).sqrt().ln();
    (xm1 + 0.5).mul_add(t.ln(), log_prefix + a.ln()) - t
}

/// Natural log of the binomial coefficient C(n, k).
#[inline]
pub fn ln_choose(n: u64, k: u64) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let k = k.min(n - k);
    if k == 0 {
        return 0.0;
    }
    ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0)
}

// ── Hypergeometric / Fisher's exact test ─────────────────────────────────────

/// One-sided (upper tail) hypergeometric p-value: P(X ≥ k).
///
/// Parameters follow the Fisher's exact test convention:
/// - `k`  — observed overlap between query and pathway
/// - `n`  — query set size
/// - `kk` — pathway size (successes in population)
/// - `nn` — background gene universe size
///
/// Computed in log-space to avoid underflow on large gene sets.
pub fn hypergeometric_pvalue(k: usize, n: usize, kk: usize, nn: usize) -> f64 {
    let (k, n, kk, nn) = (k as u64, n as u64, kk as u64, nn as u64);
    let max_k = n.min(kk);
    if k > max_k {
        return 0.0;
    }
    if k == 0 {
        return 1.0;
    }
    let log_denom = ln_choose(nn, n);
    let p: f64 = (k..=max_k)
        .filter(|&i| n >= i && nn >= kk && nn - kk >= n - i) // feasibility guard
        .map(|i| (ln_choose(kk, i) + ln_choose(nn - kk, n - i) - log_denom).exp())
        .sum();
    p.clamp(0.0, 1.0)
}

// ── Benjamini-Hochberg FDR correction ────────────────────────────────────────

/// Adjust a slice of p-values using the Benjamini-Hochberg (1995) procedure.
///
/// Returns adjusted p-values (q-values) in the same order as the input.
/// The output is guaranteed to be monotone and ≤ 1.0.
pub fn benjamini_hochberg(pvalues: &[f64]) -> Vec<f64> {
    let m = pvalues.len();
    if m == 0 {
        return Vec::new();
    }

    // Sort indices by ascending p-value
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_unstable_by(|&a, &b| {
        pvalues[a]
            .partial_cmp(&pvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut padj = vec![1.0f64; m];
    let mut running_min = 1.0f64;

    // Traverse from largest rank downward, taking the running minimum
    for rank in (0..m).rev() {
        let i = order[rank];
        let adjusted = (pvalues[i] * m as f64 / (rank + 1) as f64).min(1.0);
        running_min = running_min.min(adjusted);
        padj[i] = running_min;
    }

    padj
}

// ── t-distribution p-value ────────────────────────────────────────────────────

/// Two-tailed p-value for a t-statistic with `df` degrees of freedom.
///
/// Uses the relationship: P(|T| > t) = I(df/(df+t²); df/2, 1/2)
/// where I is the regularised incomplete beta function.
pub fn t_pvalue(t: f64, df: f64) -> f64 {
    if df <= 0.0 {
        return 1.0;
    }
    let t2 = t * t;
    let x = df / (df + t2);
    regularised_incomplete_beta(x, df / 2.0, 0.5)
}

/// Regularised incomplete beta function I(x; a, b) via Lentz continued fractions.
///
/// Switches tails at x = (a+1)/(a+b+2) for better convergence.
fn regularised_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }

    let threshold = (a + 1.0) / (a + b + 2.0);
    let (xx, aa, bb, flipped) = if x < threshold {
        (x, a, b, false)
    } else {
        (1.0 - x, b, a, true)
    };

    let ln_prefix =
        aa * xx.ln() + bb * (1.0 - xx).ln() + ln_gamma(aa + bb) - ln_gamma(aa) - ln_gamma(bb);

    let cf = betacf(xx, aa, bb);
    let ibeta = (ln_prefix.exp() * cf / aa).clamp(0.0, 1.0);

    if flipped {
        1.0 - ibeta
    } else {
        ibeta
    }
}

/// Lentz's continued fraction for the incomplete beta function.
fn betacf(x: f64, a: f64, b: f64) -> f64 {
    const MAX_ITER: usize = 200;
    const EPS: f64 = 3.0e-7;
    const FPMIN: f64 = 1.0e-300;

    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0f64;
    let mut d = {
        let v = 1.0 - qab * x / qap;
        if v.abs() < FPMIN { FPMIN } else { v }.recip()
    };
    let mut h = d;

    for m in 1..=MAX_ITER {
        let mf = m as f64;
        let m2 = 2.0 * mf;

        // Even step
        let aa = mf * (b - mf) * x / ((qam + m2) * (a + m2));
        d = {
            let v = 1.0 + aa * d;
            if v.abs() < FPMIN { FPMIN } else { v }.recip()
        };
        c = {
            let v = 1.0 + aa / c;
            if v.abs() < FPMIN {
                FPMIN
            } else {
                v
            }
        };
        h *= d * c;

        // Odd step
        let aa = -(a + mf) * (qab + mf) * x / ((a + m2) * (qap + m2));
        d = {
            let v = 1.0 + aa * d;
            if v.abs() < FPMIN { FPMIN } else { v }.recip()
        };
        c = {
            let v = 1.0 + aa / c;
            if v.abs() < FPMIN {
                FPMIN
            } else {
                v
            }
        };
        let delta = d * c;
        h *= delta;

        if (delta - 1.0).abs() < EPS {
            break;
        }
    }

    h
}

// ── Welch's t-test ───────────────────────────────────────────────────────────

/// Two-sample Welch's t-test.
///
/// Returns `(t_statistic, p_value)` where the p-value is two-tailed.
/// Returns `None` when either group has fewer than 2 observations or zero variance.
pub fn welch_t_test(group1: &[f64], group2: &[f64]) -> Option<(f64, f64)> {
    let n1 = group1.len();
    let n2 = group2.len();
    if n1 < 2 || n2 < 2 {
        return None;
    }

    let mean1 = group1.iter().sum::<f64>() / n1 as f64;
    let mean2 = group2.iter().sum::<f64>() / n2 as f64;

    let var1 = group1.iter().map(|&x| (x - mean1).powi(2)).sum::<f64>() / (n1 - 1) as f64;
    let var2 = group2.iter().map(|&x| (x - mean2).powi(2)).sum::<f64>() / (n2 - 1) as f64;

    let se_sq = var1 / n1 as f64 + var2 / n2 as f64;
    if se_sq < 1e-30 {
        return None;
    }
    let se = se_sq.sqrt();
    let t = (mean1 - mean2) / se;

    // Welch-Satterthwaite degrees of freedom
    let df_num = se_sq.powi(2);
    let df_den =
        (var1 / n1 as f64).powi(2) / (n1 - 1) as f64 + (var2 / n2 as f64).powi(2) / (n2 - 1) as f64;
    let df = if df_den < 1e-30 {
        (n1 + n2 - 2) as f64
    } else {
        df_num / df_den
    };

    Some((t, t_pvalue(t, df)))
}

// ── Spearman rank correlation ─────────────────────────────────────────────────

/// Spearman rank correlation coefficient between two equal-length slices.
///
/// Handles ties by mid-rank averaging. Returns 0.0 for constant vectors.
pub fn spearman_r(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len();
    if n < 3 || n != b.len() {
        return 0.0;
    }
    let ra = rank_array(a);
    let rb = rank_array(b);
    pearson_of_ranks(&ra, &rb)
}

fn rank_array(v: &[f64]) -> Vec<f64> {
    let n = v.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_unstable_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut ranks = vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let val = v[idx[i]];
        let mut j = i + 1;
        while j < n && v[idx[j]] == val {
            j += 1;
        }
        // 1-indexed average rank for tied block
        let avg = (i + j + 1) as f64 / 2.0;
        for k in i..j {
            ranks[idx[k]] = avg;
        }
        i = j;
    }
    ranks
}

fn pearson_of_ranks(ra: &[f64], rb: &[f64]) -> f64 {
    let n = ra.len() as f64;
    let ma = ra.iter().sum::<f64>() / n;
    let mb = rb.iter().sum::<f64>() / n;
    let cov: f64 = ra
        .iter()
        .zip(rb.iter())
        .map(|(x, y)| (x - ma) * (y - mb))
        .sum();
    let sa = ra.iter().map(|x| (x - ma).powi(2)).sum::<f64>().sqrt();
    let sb = rb.iter().map(|y| (y - mb).powi(2)).sum::<f64>().sqrt();
    if sa < 1e-12 || sb < 1e-12 {
        0.0
    } else {
        (cov / (sa * sb)).clamp(-1.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ln_gamma() {
        // Γ(1) = 1, Γ(2) = 1, Γ(0.5) = √π
        assert!((ln_gamma(1.0)).abs() < 1e-10);
        assert!((ln_gamma(2.0)).abs() < 1e-10);
        let ln_sqrt_pi = (std::f64::consts::PI.sqrt()).ln();
        assert!((ln_gamma(0.5) - ln_sqrt_pi).abs() < 1e-10);
    }

    #[test]
    fn test_ln_choose() {
        // C(10,3) = 120
        assert!((ln_choose(10, 3) - 120.0_f64.ln()).abs() < 1e-10);
        assert_eq!(ln_choose(5, 0), 0.0);
        assert_eq!(ln_choose(5, 5), 0.0);
    }

    #[test]
    fn test_hypergeometric_pvalue() {
        // P(X >= 0) = 1.0 by definition
        let p_all = hypergeometric_pvalue(0, 10, 5, 100);
        assert!((p_all - 1.0).abs() < 1e-10);
        // k > max possible overlap → p = 0
        assert_eq!(hypergeometric_pvalue(11, 10, 10, 20), 0.0);
        // High overlap: P(X >= 8 | n=10, K=10, N=20) ≈ 0.0115 — significant
        // C(10,8)*C(10,2)/C(20,10) + ... ≈ 0.0115
        let p_sig = hypergeometric_pvalue(8, 10, 10, 20);
        assert!(
            p_sig < 0.05,
            "Expected significant enrichment, got p={p_sig}"
        );
        assert!(p_sig > 0.0);
    }

    #[test]
    fn test_benjamini_hochberg() {
        let pvals = vec![0.01, 0.04, 0.03, 0.20];
        let padj = benjamini_hochberg(&pvals);
        // BH: sorted [0.01, 0.03, 0.04, 0.20]; adjusted = [0.04, 0.06, 0.08, 0.20], then cummin from right
        assert!(padj[0] <= padj[1]); // monotone
        assert!(padj.iter().all(|&p| p <= 1.0));
    }

    #[test]
    fn test_welch_t_test() {
        // Two clearly different groups
        let g1 = vec![1.0, 1.1, 0.9, 1.05, 0.95];
        let g2 = vec![5.0, 5.1, 4.9, 5.05, 4.95];
        let (t, p) = welch_t_test(&g1, &g2).unwrap();
        assert!(t.abs() > 5.0); // large t
        assert!(p < 0.001); // highly significant
    }

    #[test]
    fn test_spearman_r() {
        // Perfect monotone → rs = 1
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let b = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((spearman_r(&a, &b) - 1.0).abs() < 1e-10);
        // Reversed → rs = -1
        let c = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        assert!((spearman_r(&a, &c) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_t_pvalue_symmetry() {
        // t=0 → p=1
        assert!((t_pvalue(0.0, 10.0) - 1.0).abs() < 1e-6);
        // large |t| → p ≈ 0
        assert!(t_pvalue(10.0, 100.0) < 1e-10);
        // symmetric
        assert!((t_pvalue(2.0, 30.0) - t_pvalue(-2.0, 30.0)).abs() < 1e-10);
    }
}
