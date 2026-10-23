//! Small dense linear algebra for the p×p normal-equation systems of the
//! negative-binomial GLM (p = number of model coefficients, typically ≤ 10).

// Row/column indices read more clearly than iterator chains in matrix kernels.
#![allow(clippy::needless_range_loop)]

/// Row-major dense matrix.
pub(crate) type Matrix = Vec<Vec<f64>>;

/// Xᵀ diag(w) X for an m×p design `x`.
pub(crate) fn xtwx(x: &[Vec<f64>], w: &[f64]) -> Matrix {
    let p = x.first().map_or(0, Vec::len);
    let mut out = vec![vec![0.0; p]; p];
    for (row, &wj) in x.iter().zip(w) {
        for a in 0..p {
            let wa = wj * row[a];
            for b in 0..=a {
                out[a][b] += wa * row[b];
            }
        }
    }
    for a in 0..p {
        for b in (a + 1)..p {
            out[a][b] = out[b][a];
        }
    }
    out
}

/// Xᵀ diag(w) z for an m×p design `x`.
pub(crate) fn xtwz(x: &[Vec<f64>], w: &[f64], z: &[f64]) -> Vec<f64> {
    let p = x.first().map_or(0, Vec::len);
    let mut out = vec![0.0; p];
    for ((row, &wj), &zj) in x.iter().zip(w).zip(z) {
        for (o, &xa) in out.iter_mut().zip(row) {
            *o += wj * xa * zj;
        }
    }
    out
}

/// Lower-triangular Cholesky factor L with L Lᵀ = A, or `None` when A is not
/// numerically positive definite.
pub(crate) fn cholesky(a: &Matrix) -> Option<Matrix> {
    let p = a.len();
    let mut l = vec![vec![0.0; p]; p];
    for i in 0..p {
        for j in 0..=i {
            let s: f64 = (0..j).map(|k| l[i][k] * l[j][k]).sum();
            if i == j {
                let d = a[i][i] - s;
                if d <= 0.0 || !d.is_finite() {
                    return None;
                }
                l[i][i] = d.sqrt();
            } else {
                l[i][j] = (a[i][j] - s) / l[j][j];
            }
        }
    }
    Some(l)
}

/// Solve L Lᵀ x = b.
pub(crate) fn chol_solve(l: &Matrix, b: &[f64]) -> Vec<f64> {
    let p = l.len();
    let mut y = vec![0.0; p];
    for i in 0..p {
        let s: f64 = (0..i).map(|k| l[i][k] * y[k]).sum();
        y[i] = (b[i] - s) / l[i][i];
    }
    let mut x = vec![0.0; p];
    for i in (0..p).rev() {
        let s: f64 = ((i + 1)..p).map(|k| l[k][i] * x[k]).sum();
        x[i] = (y[i] - s) / l[i][i];
    }
    x
}

/// (L Lᵀ)⁻¹.
pub(crate) fn chol_inverse(l: &Matrix) -> Matrix {
    let p = l.len();
    let mut inv = vec![vec![0.0; p]; p];
    for c in 0..p {
        let mut e = vec![0.0; p];
        e[c] = 1.0;
        let col = chol_solve(l, &e);
        for r in 0..p {
            inv[r][c] = col[r];
        }
    }
    inv
}

/// ln det(L Lᵀ).
pub(crate) fn chol_logdet(l: &Matrix) -> f64 {
    2.0 * l
        .iter()
        .enumerate()
        .map(|(i, row)| row[i].ln())
        .sum::<f64>()
}

/// A · B.
pub(crate) fn mat_mul(a: &Matrix, b: &Matrix) -> Matrix {
    let n = a.len();
    let k = b.len();
    let m = b.first().map_or(0, Vec::len);
    let mut out = vec![vec![0.0; m]; n];
    for i in 0..n {
        for t in 0..k {
            let ait = a[i][t];
            for j in 0..m {
                out[i][j] += ait * b[t][j];
            }
        }
    }
    out
}

/// tr(A · B) for square matrices of equal size.
pub(crate) fn trace_of_product(a: &Matrix, b: &Matrix) -> f64 {
    let p = a.len();
    (0..p)
        .map(|i| (0..p).map(|j| a[i][j] * b[j][i]).sum::<f64>())
        .sum()
}

/// Ordinary least-squares operators for a full-rank design X:
/// the coefficient map (XᵀX)⁻¹Xᵀ (p×m) and the hat matrix X(XᵀX)⁻¹Xᵀ (m×m).
pub(crate) struct OlsOperators {
    pub coef: Matrix,
    pub hat: Matrix,
}

impl OlsOperators {
    pub(crate) fn new(x: &[Vec<f64>]) -> Option<Self> {
        let m = x.len();
        let ones = vec![1.0; m];
        let l = cholesky(&xtwx(x, &ones))?;
        let xtx_inv = chol_inverse(&l);
        let p = xtx_inv.len();
        let coef: Matrix = (0..p)
            .map(|a| {
                (0..m)
                    .map(|j| (0..p).map(|b| xtx_inv[a][b] * x[j][b]).sum())
                    .collect()
            })
            .collect();
        let hat: Matrix = (0..m)
            .map(|i| {
                (0..m)
                    .map(|j| (0..p).map(|a| x[i][a] * coef[a][j]).sum())
                    .collect()
            })
            .collect();
        Some(Self { coef, hat })
    }

    /// Least-squares coefficients for response `y`.
    pub(crate) fn coefficients(&self, y: &[f64]) -> Vec<f64> {
        self.coef
            .iter()
            .map(|row| row.iter().zip(y).map(|(a, b)| a * b).sum())
            .collect()
    }

    /// Fitted values X β̂ for response `y`.
    pub(crate) fn fitted(&self, y: &[f64]) -> Vec<f64> {
        self.hat
            .iter()
            .map(|row| row.iter().zip(y).map(|(a, b)| a * b).sum())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cholesky_round_trip() {
        let a = vec![vec![4.0, 2.0], vec![2.0, 3.0]];
        let l = cholesky(&a).unwrap();
        let x = chol_solve(&l, &[2.0, 1.0]);
        assert!((4.0 * x[0] + 2.0 * x[1] - 2.0).abs() < 1e-12);
        assert!((2.0 * x[0] + 3.0 * x[1] - 1.0).abs() < 1e-12);
        assert!((chol_logdet(&l) - 8.0f64.ln()).abs() < 1e-12);
        let inv = chol_inverse(&l);
        let id = mat_mul(&a, &inv);
        assert!((id[0][0] - 1.0).abs() < 1e-12 && id[0][1].abs() < 1e-12);
    }

    #[test]
    fn ols_two_group_fitted_values_are_group_means() {
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let ols = OlsOperators::new(&x).unwrap();
        let fitted = ols.fitted(&[1.0, 3.0, 10.0, 20.0]);
        for (f, e) in fitted.iter().zip([2.0, 2.0, 15.0, 15.0]) {
            assert!((f - e).abs() < 1e-12);
        }
        let coef = ols.coefficients(&[1.0, 3.0, 10.0, 20.0]);
        assert!((coef[0] - 2.0).abs() < 1e-12 && (coef[1] - 13.0).abs() < 1e-12);
    }
}
