//! Outlier detection by Cook's distance and independent filtering of p-values
//! by mean count (Bourgon, Gentleman & Huber 2010), following DESeq2's
//! `results()` defaults: Cook's cutoff F(0.99; p, m − p) and a filter
//! threshold chosen from 50 quantiles of the mean count at α = 0.1.

use biomics_core::statistics::{benjamini_hochberg, f_quantile};

use super::dispersion::trimmed_mean;

/// Samples grouped by identical design rows ("cells").
pub(crate) struct DesignCells {
    /// Cell index for each sample.
    pub cell_of: Vec<usize>,
    /// Number of samples in each cell.
    pub sizes: Vec<usize>,
}

impl DesignCells {
    pub(crate) fn new(x: &[Vec<f64>]) -> Self {
        let mut keys: Vec<Vec<u64>> = Vec::new();
        let mut cell_of = Vec::with_capacity(x.len());
        for row in x {
            let key: Vec<u64> = row.iter().map(|v| v.to_bits()).collect();
            let idx = keys.iter().position(|k| *k == key).unwrap_or_else(|| {
                keys.push(key);
                keys.len() - 1
            });
            cell_of.push(idx);
        }
        let mut sizes = vec![0; keys.len()];
        for &c in &cell_of {
            sizes[c] += 1;
        }
        Self { cell_of, sizes }
    }

    fn in_large_cell(&self, sample: usize) -> bool {
        self.sizes[self.cell_of[sample]] >= 3
    }
}

/// Trim proportion and variance scale for a cell of `n` replicates.
fn trim_and_scale(n: usize) -> (f64, f64) {
    if (n as f64) <= 3.5 {
        (1.0 / 3.0, 2.04)
    } else if (n as f64) <= 23.5 {
        (1.0 / 4.0, 1.86)
    } else {
        (1.0 / 8.0, 1.51)
    }
}

/// Outlier-robust method-of-moments dispersion used for Cook's distances:
/// trimmed within-cell variances (cells with ≥ 3 replicates) or a trimmed
/// overall variance, floored at 0.04.
fn robust_dispersion(norm: &[f64], cells: &DesignCells) -> f64 {
    let large: Vec<usize> = (0..cells.sizes.len())
        .filter(|&c| cells.sizes[c] >= 3)
        .collect();
    let variance = if large.is_empty() {
        let center = trimmed_mean(norm, 1.0 / 8.0);
        let sq: Vec<f64> = norm.iter().map(|v| (v - center).powi(2)).collect();
        1.51 * trimmed_mean(&sq, 1.0 / 8.0)
    } else {
        large
            .iter()
            .map(|&c| {
                let vals: Vec<f64> = norm
                    .iter()
                    .zip(&cells.cell_of)
                    .filter(|(_, &cc)| cc == c)
                    .map(|(&v, _)| v)
                    .collect();
                let (trim, scale) = trim_and_scale(vals.len());
                let center = trimmed_mean(&vals, trim);
                let sq: Vec<f64> = vals.iter().map(|v| (v - center).powi(2)).collect();
                scale * trimmed_mean(&sq, trim)
            })
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let mean = norm.iter().sum::<f64>() / norm.len() as f64;
    ((variance - mean) / (mean * mean)).max(0.04)
}

/// Cook's distances for one gene and the maximum over samples in cells with at
/// least three replicates (`None` when no such cell exists or m ≤ p).
pub(crate) fn cooks_distances(
    counts: &[f64],
    norm: &[f64],
    mu: &[f64],
    hat: &[f64],
    p: usize,
    cells: &DesignCells,
) -> (Vec<f64>, Option<f64>) {
    let alpha = robust_dispersion(norm, cells);
    let cooks: Vec<f64> = counts
        .iter()
        .zip(mu)
        .zip(hat)
        .map(|((&y, &m), &h)| {
            let pearson_sq = (y - m).powi(2) / (m + alpha * m * m);
            pearson_sq / p as f64 * h / (1.0 - h).powi(2)
        })
        .collect();
    let m = counts.len();
    let max = if m > p {
        (0..m)
            .filter(|&j| cells.in_large_cell(j))
            .map(|j| cooks[j])
            .reduce(f64::max)
    } else {
        None
    };
    (cooks, max)
}

/// Default Cook's distance cutoff F(0.99; p, m − p).
pub(crate) fn cooks_cutoff(m: usize, p: usize) -> f64 {
    f_quantile(0.99, p as f64, (m - p) as f64)
}

/// R's default (type 7) sample quantile of `sorted` (ascending) at `prob`.
pub(crate) fn quantile_type7(sorted: &[f64], prob: f64) -> f64 {
    let n = sorted.len();
    let index = (n - 1) as f64 * prob;
    let lo = index.floor() as usize;
    let hi = index.ceil() as usize;
    let q = sorted[lo];
    let h = index - lo as f64;
    if h > 0.0 && sorted[hi] != q {
        (1.0 - h) * q + h * sorted[hi]
    } else {
        q
    }
}

/// Benjamini–Hochberg adjustment ignoring NaN entries (R's `p.adjust` with NA).
pub(crate) fn bh_ignoring_nan(p: &[f64]) -> Vec<f64> {
    let idx: Vec<usize> = (0..p.len()).filter(|&i| !p[i].is_nan()).collect();
    let adj = benjamini_hochberg(&idx.iter().map(|&i| p[i]).collect::<Vec<_>>());
    let mut out = vec![f64::NAN; p.len()];
    for (&i, a) in idx.iter().zip(adj) {
        out[i] = a;
    }
    out
}

/// Locally weighted scatterplot smoothing (Cleveland 1979): local linear fits
/// with tricube weights over the nearest ⌊f·n⌋ points, followed by `iter`
/// robustness iterations with bisquare weights on residuals scaled by six
/// times their median. Points within `delta` of the last fitted abscissa are
/// linearly interpolated. `x` must be sorted ascending.
pub(crate) fn lowess(x: &[f64], y: &[f64], f: f64, iter: usize, delta: f64) -> Vec<f64> {
    let n = x.len();
    if n < 2 {
        return y.to_vec();
    }
    let span = ((f * n as f64 + 1e-7) as usize).clamp(2, n);
    let range = x[n - 1] - x[0];
    let mut fitted = vec![0.0; n];
    let mut robustness = vec![1.0; n];
    let mut weights = vec![0.0; n];

    for it in 0..=iter {
        let robust = it > 0;
        let (mut left, mut right) = (0usize, span - 1);
        let mut last: Option<usize> = None;
        let mut i = 0usize;
        loop {
            // Slide the window right while that brings its edge closer to x[i].
            while right + 1 < n && x[i] - x[left] > x[right + 1] - x[i] {
                left += 1;
                right += 1;
            }
            fitted[i] = local_linear_fit(
                x,
                y,
                i,
                left,
                right,
                range,
                robust,
                &robustness,
                &mut weights,
            )
            .unwrap_or(y[i]);

            if let Some(prev) = last {
                if prev + 1 < i {
                    let width = x[i] - x[prev];
                    for j in (prev + 1)..i {
                        let t = (x[j] - x[prev]) / width;
                        fitted[j] = t * fitted[i] + (1.0 - t) * fitted[prev];
                    }
                }
            }
            let mut current = i;
            let cut = x[current] + delta;
            let mut next = current + 1;
            while next < n && x[next] <= cut {
                if x[next] == x[current] {
                    fitted[next] = fitted[current];
                    current = next;
                }
                next += 1;
            }
            last = Some(current);
            if current >= n - 1 {
                break;
            }
            i = (current + 1).max(next - 1);
        }

        if it == iter {
            break;
        }
        let residuals: Vec<f64> = y.iter().zip(&fitted).map(|(a, b)| a - b).collect();
        let mean_abs = residuals.iter().map(|r| r.abs()).sum::<f64>() / n as f64;
        let mut abs_sorted: Vec<f64> = residuals.iter().map(|r| r.abs()).collect();
        abs_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let scale = 6.0 * super::dispersion::median_of_sorted(&abs_sorted);
        if scale < 1e-7 * mean_abs {
            break;
        }
        for (w, r) in robustness.iter_mut().zip(&residuals) {
            let u = r.abs();
            *w = if u <= 0.001 * scale {
                1.0
            } else if u <= 0.999 * scale {
                (1.0 - (u / scale).powi(2)).powi(2)
            } else {
                0.0
            };
        }
    }
    fitted
}

/// Weighted local linear fit at x[i] over the window [left, right] (extended
/// to the right over ties). Returns `None` when every weight is zero.
#[allow(clippy::too_many_arguments)]
fn local_linear_fit(
    x: &[f64],
    y: &[f64],
    i: usize,
    left: usize,
    right: usize,
    range: f64,
    robust: bool,
    robustness: &[f64],
    weights: &mut [f64],
) -> Option<f64> {
    let n = x.len();
    let xs = x[i];
    let h = (xs - x[left]).max(x[right] - xs);
    let (h_hi, h_lo) = (0.999 * h, 0.001 * h);
    let mut total = 0.0;
    let mut j = left;
    while j < n {
        weights[j] = 0.0;
        let r = (x[j] - xs).abs();
        if r <= h_hi {
            let mut w = if r <= h_lo {
                1.0
            } else {
                (1.0 - (r / h).powi(3)).powi(3)
            };
            if robust {
                w *= robustness[j];
            }
            weights[j] = w;
            total += w;
        } else if x[j] > xs {
            break;
        }
        j += 1;
    }
    let end = j; // exclusive
    if total <= 0.0 {
        return None;
    }
    for w in &mut weights[left..end] {
        *w /= total;
    }
    if h > 0.0 {
        let center: f64 = (left..end).map(|k| weights[k] * x[k]).sum();
        let spread: f64 = (left..end)
            .map(|k| weights[k] * (x[k] - center).powi(2))
            .sum();
        if spread.sqrt() > 0.001 * range {
            let slope = (xs - center) / spread;
            for k in left..end {
                weights[k] *= slope * (x[k] - center) + 1.0;
            }
        }
    }
    Some((left..end).map(|k| weights[k] * y[k]).sum())
}

/// Outcome of independent filtering.
pub(crate) struct Filtered {
    /// Adjusted p-values; NaN for genes removed by the filter or without a p-value.
    pub padj: Vec<f64>,
    /// Mean-count threshold below which genes were filtered.
    pub threshold: f64,
}

/// Choose the mean-count filter that maximises discoveries at `alpha` (unless
/// that maximum is ≤ 10 rejections) and adjust the surviving p-values.
pub(crate) fn independent_filtering(base_mean: &[f64], pvalue: &[f64], alpha: f64) -> Filtered {
    let n = base_mean.len();
    let lower = base_mean.iter().filter(|&&b| b == 0.0).count() as f64 / n as f64;
    let upper = if lower < 0.95 { 0.95 } else { 1.0 };
    const STEPS: usize = 50;
    let theta: Vec<f64> = (0..STEPS)
        .map(|k| {
            if k == STEPS - 1 {
                upper
            } else {
                lower + k as f64 * (upper - lower) / (STEPS - 1) as f64
            }
        })
        .collect();
    let mut sorted = base_mean.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let cutoffs: Vec<f64> = theta.iter().map(|&t| quantile_type7(&sorted, t)).collect();

    let adjusted: Vec<Vec<f64>> = cutoffs
        .iter()
        .map(|&cut| {
            let masked: Vec<f64> = base_mean
                .iter()
                .zip(pvalue)
                .map(|(&b, &p)| if b >= cut { p } else { f64::NAN })
                .collect();
            bh_ignoring_nan(&masked)
        })
        .collect();
    let rejections: Vec<f64> = adjusted
        .iter()
        .map(|padj| padj.iter().filter(|&&q| q < alpha).count() as f64)
        .collect();

    let smooth = lowess(&theta, &rejections, 1.0 / 5.0, 3, 0.01 * (upper - lower));
    let max_rej = rejections.iter().cloned().fold(0.0, f64::max);
    let j = if max_rej <= 10.0 {
        0
    } else {
        let resid: Vec<f64> = if rejections.iter().all(|&r| r == 0.0) {
            vec![0.0]
        } else {
            rejections
                .iter()
                .zip(&smooth)
                .filter(|(&r, _)| r > 0.0)
                .map(|(r, s)| r - s)
                .collect()
        };
        let rmse = (resid.iter().map(|r| r * r).sum::<f64>() / resid.len() as f64).sqrt();
        let thresh = smooth.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - rmse;
        rejections.iter().position(|&r| r > thresh).unwrap_or(0)
    };
    Filtered {
        padj: adjusted[j].clone(),
        threshold: cutoffs[j],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_type7_matches_r() {
        // R: quantile(c(1, 2, 3, 4, 10), c(0, 0.3, 0.5, 0.95, 1))
        let v = [1.0, 2.0, 3.0, 4.0, 10.0];
        let expected = [1.0, 2.2, 3.0, 8.8, 10.0];
        for (p, e) in [0.0, 0.3, 0.5, 0.95, 1.0].iter().zip(expected) {
            assert!((quantile_type7(&v, *p) - e).abs() < 1e-12);
        }
    }

    #[test]
    #[allow(clippy::excessive_precision)] // reference values copied verbatim from R
    fn lowess_matches_r() {
        // R 4.2: x <- 0:29; y <- 2*x + 1; y[16] <- y[16] + 100; lowess(x, y, f = 0.3)
        // (the median absolute residual is 0, so no robustness iteration runs)
        let x: Vec<f64> = (0..30).map(f64::from).collect();
        let mut y: Vec<f64> = x.iter().map(|v| 2.0 * v + 1.0).collect();
        y[15] += 100.0;
        let fit = lowess(&x, &y, 0.3, 3, 0.01 * 29.0);
        assert!((fit[15] - 52.579_616_324_326_5).abs() < 1e-9, "{}", fit[15]);

        // R 4.2: set.seed(1); x <- sort(runif(50)); y <- sin(6*x) + rnorm(50, 0, 0.3);
        //        y[c(10, 30)] <- y[c(10, 30)] + 4; lowess(x, y, f = 0.3)$y
        let x: [f64; 50] = [
            0.013390333158895,
            0.023331202333793,
            0.061786270467564,
            0.107943625887856,
            0.125555095961317,
            0.176556752528995,
            0.186217601411045,
            0.201681931037456,
            0.205974574899301,
            0.212142521282658,
            0.265508663142100,
            0.267220668727532,
            0.340348996687680,
            0.372123899636790,
            0.380035179434344,
            0.382387957070023,
            0.384103718213737,
            0.386114092543721,
            0.411274429643527,
            0.477230065036565,
            0.482080115471035,
            0.493541307048872,
            0.497699242085218,
            0.529719580197707,
            0.553036311641335,
            0.572853363351896,
            0.599565825425088,
            0.629114043898880,
            0.647060193819925,
            0.651673766085878,
            0.660797792486846,
            0.668466738192365,
            0.687022846657783,
            0.692731556482613,
            0.717618508264422,
            0.723710946040228,
            0.732313738670200,
            0.769841419998556,
            0.777445221319795,
            0.782932762289420,
            0.789356231689453,
            0.794239860726520,
            0.820946294115856,
            0.827373318606988,
            0.869690845720470,
            0.898389684967697,
            0.908207789994776,
            0.934705231105909,
            0.944675268605351,
            0.991906094830483,
        ];
        let y: [f64; 50] = [
            0.063416972570560,
            0.092791801722235,
            -0.078941317009757,
            0.459878297804805,
            0.809454388178614,
            1.279636754404694,
            0.868086959051133,
            1.051949812811377,
            0.928285603758709,
            4.542825014241849,
            0.875253984208004,
            0.881184035365574,
            0.873185829292365,
            1.118804564325278,
            0.987695956467360,
            0.700115405930268,
            0.666609019005927,
            0.943574431246072,
            0.791073123423191,
            0.068010385487778,
            0.034294882478779,
            0.288743391836629,
            0.385332403702852,
            -0.070420436589693,
            0.088624017404556,
            -0.171812784041939,
            -0.623790764155809,
            -0.489304136765996,
            -1.013664161254204,
            3.734885510243517,
            -0.139201210461018,
            -0.873984292817562,
            -1.144040893048067,
            -0.678457373094057,
            -0.958956348583249,
            -0.211797421686153,
            -0.961476150086406,
            -0.788725132349742,
            -0.990461081196485,
            -1.222872557012490,
            -0.943080329915960,
            -1.540080757643511,
            -0.537673588606777,
            -0.922476833066605,
            -0.223024919320555,
            -0.636207797027835,
            -0.953567613551200,
            -0.441643435466070,
            -0.857297005595073,
            -0.701786932529723,
        ];
        let expected: [f64; 50] = [
            0.047592571567251,
            0.095892302977744,
            0.278814611014958,
            0.487513513700286,
            0.564239370640114,
            0.786604624376819,
            0.829683274114142,
            0.884044208970157,
            0.894529103155691,
            0.909265190891044,
            0.923556619913958,
            0.922067237757995,
            0.873458983566495,
            0.813450724562736,
            0.778719774433783,
            0.771327379927117,
            0.765936482350906,
            0.759619914253713,
            0.664861774932387,
            0.301067341061851,
            0.274807694295742,
            0.204044553246133,
            0.176701939136162,
            -0.004385532907839,
            -0.133846237079294,
            -0.236943789989285,
            -0.389207226593852,
            -0.546558802678391,
            -0.646278425346227,
            -0.669195800700683,
            -0.713571498822129,
            -0.740282944173038,
            -0.787751672393930,
            -0.801632283703431,
            -0.852007868291888,
            -0.869317462303104,
            -0.887722739852492,
            -0.920421344393620,
            -0.912105699707683,
            -0.908583367746471,
            -0.895852347958368,
            -0.891836419141638,
            -0.831260606030075,
            -0.815775551575063,
            -0.760002216927058,
            -0.728961227075161,
            -0.718445703704725,
            -0.691946245367474,
            -0.682816080156563,
            -0.647320173721880,
        ];
        let delta = 0.01 * (x[49] - x[0]);
        let got = lowess(&x, &y, 0.3, 3, delta);
        for (i, (g, e)) in got.iter().zip(expected).enumerate() {
            assert!((g - e).abs() < 1e-9, "i={i}: {g} vs {e}");
        }
    }

    #[test]
    fn bh_skips_nan() {
        let adj = bh_ignoring_nan(&[0.01, f64::NAN, 0.04]);
        assert!(adj[1].is_nan());
        assert!((adj[0] - 0.02).abs() < 1e-12 && (adj[2] - 0.04).abs() < 1e-12);
    }
}
