//! Validate `biomics_core::statistics` against SciPy / statsmodels.
//!
//! Reference values in `tests/golden/scipy/statistics.json` are produced by
//! `validation/scipy_reference.py` (SciPy 1.18, statsmodels 0.15).

use serde_json::Value;

use biomics_core::statistics::{
    benjamini_hochberg, chi2_sf, digamma, f_quantile, hypergeometric_pvalue, ln_gamma, normal_sf,
    spearman_r, t_pvalue, trigamma, welch_t_test,
};

fn golden() -> Value {
    serde_json::from_str(include_str!("../../tests/golden/scipy/statistics.json")).unwrap()
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

fn rows(g: &Value, key: &str) -> Vec<Vec<f64>> {
    g[key]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| r.as_array().unwrap().iter().map(f).collect())
        .collect()
}

fn assert_rel(ours: f64, reference: f64, rtol: f64, what: &str) {
    let err = if reference == 0.0 {
        ours.abs()
    } else {
        ((ours - reference) / reference).abs()
    };
    assert!(
        err <= rtol,
        "{what}: ours {ours} vs SciPy {reference} (rel err {err:.2e})"
    );
}

#[test]
fn t_distribution_p_values() {
    for r in rows(&golden(), "t_pvalue") {
        let (t, df, p) = (r[0], r[1], r[2]);
        assert_rel(t_pvalue(t, df), p, 1e-9, &format!("t={t}, df={df}"));
    }
}

#[test]
fn hypergeometric_upper_tail() {
    for r in rows(&golden(), "hypergeometric") {
        let (k, n, kk, nn) = (r[0] as usize, r[1] as usize, r[2] as usize, r[3] as usize);
        assert_rel(
            hypergeometric_pvalue(k, n, kk, nn),
            r[4],
            1e-9,
            &format!("P(X>={k}; n={n}, K={kk}, N={nn})"),
        );
    }
}

#[test]
fn log_gamma() {
    for r in rows(&golden(), "ln_gamma") {
        let ours = ln_gamma(r[0]);
        assert!(
            (ours - r[1]).abs() <= 1e-12 * r[1].abs().max(1.0),
            "lnΓ({})",
            r[0]
        );
    }
}

#[test]
fn benjamini_hochberg_matches_statsmodels() {
    let g = golden();
    let p: Vec<f64> = g["benjamini_hochberg"]["p"]
        .as_array()
        .unwrap()
        .iter()
        .map(f)
        .collect();
    let expected: Vec<f64> = g["benjamini_hochberg"]["padj"]
        .as_array()
        .unwrap()
        .iter()
        .map(f)
        .collect();
    for (ours, e) in benjamini_hochberg(&p).iter().zip(expected) {
        assert!((ours - e).abs() < 1e-15);
    }
}

#[test]
fn welch_t_test_matches_scipy() {
    for case in golden()["welch"].as_array().unwrap() {
        let a: Vec<f64> = case["a"].as_array().unwrap().iter().map(f).collect();
        let b: Vec<f64> = case["b"].as_array().unwrap().iter().map(f).collect();
        let (t, p) = welch_t_test(&a, &b).unwrap();
        assert_rel(t, f(&case["t"]), 1e-12, "Welch t");
        assert_rel(p, f(&case["p"]), 1e-9, "Welch p");
    }
}

#[test]
fn spearman_matches_scipy() {
    let g = golden();
    let x: Vec<f64> = g["spearman"]["x"]
        .as_array()
        .unwrap()
        .iter()
        .map(f)
        .collect();
    let y: Vec<f64> = g["spearman"]["y"]
        .as_array()
        .unwrap()
        .iter()
        .map(f)
        .collect();
    assert_rel(
        spearman_r(&x, &y),
        f(&g["spearman"]["rho"]),
        1e-12,
        "Spearman rho",
    );
}

#[test]
fn normal_tail() {
    for r in rows(&golden(), "normal_sf") {
        if r[1] > 1e-300 {
            assert_rel(normal_sf(r[0]), r[1], 1e-12, &format!("P(Z>{})", r[0]));
        }
    }
}

#[test]
fn polygamma() {
    let g = golden();
    for r in rows(&g, "digamma") {
        assert_rel(digamma(r[0]), r[1], 1e-13, &format!("ψ({})", r[0]));
    }
    for r in rows(&g, "trigamma") {
        assert_rel(trigamma(r[0]), r[1], 1e-13, &format!("ψ₁({})", r[0]));
    }
}

#[test]
fn chi_squared_tail() {
    for r in rows(&golden(), "chi2_sf") {
        assert_rel(
            chi2_sf(r[0], r[1]),
            r[2],
            1e-11,
            &format!("χ²({}) > {}", r[1], r[0]),
        );
    }
}

#[test]
fn f_distribution_quantile() {
    for r in rows(&golden(), "f_quantile") {
        assert_rel(
            f_quantile(r[0], r[1], r[2]),
            r[3],
            1e-10,
            &format!("F⁻¹({}; {}, {})", r[0], r[1], r[2]),
        );
    }
}
