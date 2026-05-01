/// Compute arithmetic mean of a slice. Returns `None` for empty input.
pub fn mean(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    Some(values.iter().sum::<f64>() / values.len() as f64)
}

/// Compute population standard deviation of a slice. Returns `None` for empty input.
pub fn std_dev(values: &[f64]) -> Option<f64> {
    let m = mean(values)?;
    let variance = values.iter().map(|v| (v - m).powi(2)).sum::<f64>() / values.len() as f64;
    Some(variance.sqrt())
}

/// Compute the p-th percentile (0–100) using linear interpolation.
/// Returns `None` for empty input.
pub fn percentile(sorted: &[f64], p: f64) -> Option<f64> {
    if sorted.is_empty() {
        return None;
    }
    if sorted.len() == 1 {
        return Some(sorted[0]);
    }
    let rank = p / 100.0 * (sorted.len() - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = (lo + 1).min(sorted.len() - 1);
    let frac = rank - lo as f64;
    Some(sorted[lo] * (1.0 - frac) + sorted[hi] * frac)
}

/// Map a value into a histogram bin index.
/// Returns the bin index in `[0, n_bins)` for `value` in `[min, max)`.
pub fn histogram_bin(value: f64, min: f64, max: f64, n_bins: usize) -> usize {
    if value >= max {
        return n_bins - 1;
    }
    if value <= min {
        return 0;
    }
    let frac = (value - min) / (max - min);
    ((frac * n_bins as f64) as usize).min(n_bins - 1)
}

/// Compute log2 fold change between two expression values.
/// Returns `None` when either value is zero (undefined).
pub fn log2_fold_change(base: f64, treatment: f64) -> Option<f64> {
    if base <= 0.0 || treatment <= 0.0 {
        return None;
    }
    Some((treatment / base).log2())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mean() {
        assert_eq!(mean(&[1.0, 2.0, 3.0]), Some(2.0));
        assert_eq!(mean(&[]), None);
    }

    #[test]
    fn test_std_dev() {
        let sd = std_dev(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]).unwrap();
        assert!((sd - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_percentile() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(percentile(&data, 0.0), Some(1.0));
        assert_eq!(percentile(&data, 100.0), Some(5.0));
        assert!((percentile(&data, 50.0).unwrap() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_histogram_bin() {
        assert_eq!(histogram_bin(0.0, 0.0, 1.0, 10), 0);
        assert_eq!(histogram_bin(0.95, 0.0, 1.0, 10), 9);
        assert_eq!(histogram_bin(1.0, 0.0, 1.0, 10), 9);
    }
}
