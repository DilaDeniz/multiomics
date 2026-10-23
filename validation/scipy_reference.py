#!/usr/bin/env python3
"""Regenerate SciPy/statsmodels reference values for
biomics_core/tests/scipy_reference.rs.

Usage (from the repository root, with validation/requirements.txt installed):
    python validation/scipy_reference.py
"""
import json
import pathlib

import numpy as np
import scipy
import statsmodels
from scipy import special, stats
from statsmodels.stats.multitest import multipletests

root = pathlib.Path(__file__).resolve().parent.parent
out = root / "tests" / "golden" / "scipy" / "statistics.json"

ref = {
    "scipy_version": scipy.__version__,
    "statsmodels_version": statsmodels.__version__,
}

# Two-sided Student t p-values: 2 * t.sf(|t|, df)
ref["t_pvalue"] = [
    [t, df, float(2 * stats.t.sf(abs(t), df))]
    for df in [1, 2, 2.5, 3, 4.7, 10, 30, 100, 1000, 5000, 20000, 100000, 400000]
    for t in [0.0, 0.1, 0.5, 1.0, 1.96, 2.5, 3.0, 5.0, 8.0, 12.0, 20.0]
]

# Upper-tail hypergeometric P(X >= k): hypergeom.sf(k - 1, N, K, n)
ref["hypergeometric"] = [
    [k, n, K, N, float(stats.hypergeom.sf(k - 1, N, K, n))]
    for (k, n, K, N) in [
        (9, 14, 17, 20000), (8, 10, 10, 20), (1, 5, 10, 100), (3, 50, 200, 20000),
        (50, 300, 400, 20000), (5, 500, 600, 1000), (100, 1000, 1500, 20000), (2, 2, 2, 4),
        (4, 14, 13, 290), (1, 3, 400, 500),
    ]
]

ref["ln_gamma"] = [[x, float(special.gammaln(x))]
                   for x in [0.1, 0.5, 1.0, 1.5, 2.0, 3.3, 10.0, 100.5, 1000.0, 20001.0, 1e6]]

p = [0.01, 0.04, 0.03, 0.20, 0.5, 0.001, 0.049, 0.9, 0.04, 0.04]
ref["benjamini_hochberg"] = {"p": p, "padj": multipletests(p, method="fdr_bh")[1].tolist()}

welch_cases = [
    ([120.3, 88.5], [230.1, 195.4]),
    ([1.0, 1.1, 0.9, 1.05, 0.95, 3.0], [5.0, 5.1, 4.9, 5.05, 4.95, 1.0, 2.0]),
    ([2.1, 2.5, 2.2, 2.4], [2.3, 2.6, 2.0, 2.8, 2.2]),
]
ref["welch"] = []
for a, b in welch_cases:
    res = stats.ttest_ind(a, b, equal_var=False)
    ref["welch"].append({"a": a, "b": b, "t": float(res.statistic), "p": float(res.pvalue)})

x = [1.0, 2.0, 2.0, 3.0, 5.0, 4.0, 4.0]
y = [2.0, 1.0, 4.0, 3.0, 7.0, 7.0, 6.0]
ref["spearman"] = {"x": x, "y": y, "rho": float(stats.spearmanr(x, y).statistic)}

ref["normal_sf"] = [[z, float(stats.norm.sf(z))] for z in np.arange(-38.0, 38.01, 0.25).round(2)]
xs = [1e-8, 1e-4, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.7, 9.99, 10.0, 50.0, 123.4, 1e5, 1e8]
ref["digamma"] = [[v, float(special.digamma(v))] for v in xs]
ref["trigamma"] = [[v, float(special.polygamma(1, v))] for v in xs]
ref["chi2_sf"] = [
    [v, df, float(stats.chi2.sf(v, df))]
    for df in [0.5, 1, 2, 3, 7.5, 30, 1000]
    for v in [1e-6, 0.5, 1.0, 3.84, 10.0, 50.0, 300.0, 1500.0]
    if stats.chi2.sf(v, df) > 1e-300
]
ref["f_quantile"] = [
    [q, d1, d2, float(stats.f.ppf(q, d1, d2))]
    for (d1, d2) in [(1, 1), (2, 4), (3, 7), (4, 20), (1, 100), (5, 2)]
    for q in [0.5, 0.9, 0.99, 0.999]
]

out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(ref, indent=1) + "\n")
print(f"wrote {out}")
