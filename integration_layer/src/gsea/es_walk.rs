//! Weighted Kolmogorov-Smirnov enrichment score walk.
//!
//! Unchanged from Subramanian et al. 2005, PNAS 102(43):15545-15550.
//! <https://doi.org/10.1073/pnas.0506580102>

/// Compute the enrichment score and the full running-sum trace for one gene set.
///
/// Returns `(es, running_sum, peak_index)` where `peak_index` is the position of
/// the maximum absolute deviation in the running sum.
///
/// `is_hit[i]` is `true` iff ranked gene `i` belongs to the pathway.
/// `n_hits` is the number of hits (pre-computed for efficiency).
pub fn compute_es(is_hit: &[bool], n_hits: usize) -> (f64, Vec<f64>, usize) {
    let n = is_hit.len();
    let n_miss = n - n_hits;

    if n_hits == 0 || n_miss == 0 {
        let rs = vec![0.0; n];
        return (0.0, rs, 0);
    }

    let hit_inc = ((n_miss as f64) / (n_hits as f64)).sqrt();
    let miss_dec = ((n_hits as f64) / (n_miss as f64)).sqrt();

    let mut running = 0.0_f64;
    let mut peak_val = 0.0_f64;
    let mut peak_idx = 0_usize;
    let mut rs = Vec::with_capacity(n);

    for (i, &hit) in is_hit.iter().enumerate() {
        if hit {
            running += hit_inc;
        } else {
            running -= miss_dec;
        }
        rs.push(running);

        if running.abs() > peak_val.abs() {
            peak_val = running;
            peak_idx = i;
        }
    }

    (peak_val, rs, peak_idx)
}

/// Sample `k` unique indices from `[0, n)` using partial Fisher-Yates driven by
/// xorshift64.  The same sequence is always produced for the same `seed`.
pub fn xorshift_sample(n: usize, k: usize, seed: u64) -> Vec<usize> {
    let mut state = if seed == 0 { 1 } else { seed };

    let xorshift = |s: &mut u64| -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    };

    let mut indices: Vec<usize> = (0..n).collect();
    for i in 0..k {
        let rand_val = xorshift(&mut state);
        let j = i + (rand_val as usize % (n - i));
        indices.swap(i, j);
    }
    indices[..k].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_top_pathway_has_positive_es() {
        let n = 100_usize;
        let k = 10_usize;
        let is_hit: Vec<bool> = (0..n).map(|i| i < k).collect();
        let (es, rs, _) = compute_es(&is_hit, k);
        assert!(es > 0.5, "expected ES > 0.5, got {es}");
        assert_eq!(rs.len(), n);
    }

    #[test]
    fn all_miss_pathway_returns_zero() {
        let is_hit = vec![false; 50];
        let (es, _, _) = compute_es(&is_hit, 0);
        assert_eq!(es, 0.0);
    }
}
