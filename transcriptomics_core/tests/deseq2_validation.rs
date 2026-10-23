//! Validate the DESeq2 implementation against R/DESeq2.
//!
//! Reference values in `tests/golden/deseq2/*.json` were produced by
//! `validation/deseq2_reference.R` (R 4.2, DESeq2 1.38) from count matrices in
//! `test_data/validation/`:
//!
//! - `sim_1000x6` — DESeq2's own simulator (`makeExampleDESeqDataSet`), 1 000
//!   genes, 3 vs 3, including all-zero genes, low counts and independent
//!   filtering; `~ condition`, Wald test.
//! - `sim_1000x6_lrt` — same data, likelihood-ratio test against `~ 1`.
//! - `sim_1000x8_batch` — 1 000 genes × 8 samples, `~ batch + condition`, so
//!   gene-wise fitted means come from the NB GLM rather than group means.
//! - `small_5x6` — five genes; DESeq2 falls back to a local-regression trend
//!   there, so only size factors and mean counts are compared.
//!
//! Observed agreement is ~1e-12 relative; the tolerance leaves room for libm
//! differences across platforms. The exception is a gene whose IRLS fit hits
//! DESeq2's 100-iteration cap (`beta_iter == 100`): R then refines it with
//! `optim(method = "L-BFGS-B")`, whose default tolerance only pins a very flat
//! likelihood down to ~1e-3, so such genes are compared at that precision.

use std::path::PathBuf;

use serde_json::Value;
use transcriptomics_core::deseq2::{
    design_matrix_two_group, design_matrix_with_batch, lrt_test, run_deseq2, run_deseq2_design,
    DeseqAnalysis, TrendKind,
};

const RTOL: f64 = 1e-8;
/// Tolerance for genes R fitted with L-BFGS-B instead of IRLS.
const OPTIM_RTOL: f64 = 1e-2;

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join(rel)
}

fn load_golden(name: &str) -> Value {
    let path = repo_path(&format!("tests/golden/deseq2/{name}.json"));
    serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap()
}

/// Read a genes × samples count TSV written by the R script.
fn load_counts(name: &str) -> (Vec<String>, Vec<Vec<f64>>, Vec<String>) {
    let text = std::fs::read_to_string(repo_path(&format!("test_data/validation/{name}"))).unwrap();
    let mut lines = text.lines();
    let samples: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .skip(1)
        .map(String::from)
        .collect();
    let (mut ids, mut counts) = (Vec::new(), Vec::new());
    for line in lines {
        let mut fields = line.split('\t');
        ids.push(fields.next().unwrap().to_string());
        counts.push(fields.map(|v| v.parse().unwrap()).collect());
    }
    (ids, counts, samples)
}

/// `None` for R's NA (JSON null).
fn num(v: &Value) -> Option<f64> {
    v.as_f64()
}

fn assert_close(ours: f64, reference: Option<f64>, what: &str) {
    assert_close_tol(ours, reference, what, RTOL);
}

fn assert_close_tol(ours: f64, reference: Option<f64>, what: &str, rtol: f64) {
    match reference {
        None => assert!(ours.is_nan(), "{what}: R reports NA, got {ours}"),
        Some(r) => {
            let tol = rtol * r.abs().max(1e-300);
            assert!(
                (ours - r).abs() <= tol,
                "{what}: ours {ours} vs R {r} (rel diff {:.2e})",
                (ours - r).abs() / r.abs()
            );
        }
    }
}

fn assert_size_factors(analysis: &DeseqAnalysis, golden: &Value) {
    let reference = golden["size_factors"].as_array().unwrap();
    for (j, (ours, r)) in analysis
        .size_factors
        .factors
        .iter()
        .zip(reference)
        .enumerate()
    {
        assert_close(*ours, num(r), &format!("size factor {j}"));
    }
}

fn assert_wald_matches(analysis: &DeseqAnalysis, golden: &Value) {
    assert_size_factors(analysis, golden);
    let coef = &golden["disp_fit_coefficients"];
    assert_eq!(analysis.dispersion_trend.kind, TrendKind::Parametric);
    assert_close(
        analysis.dispersion_trend.asympt_disp,
        num(&coef["asympt_disp"]),
        "asymptDisp",
    );
    assert_close(
        analysis.dispersion_trend.extra_pois,
        num(&coef["extra_pois"]),
        "extraPois",
    );
    assert_close(
        analysis.disp_prior_var,
        num(&golden["disp_prior_var"]),
        "dispPriorVar",
    );

    let genes = golden["genes"].as_array().unwrap();
    assert_eq!(genes.len(), analysis.results.len());
    for (ours, r) in analysis.results.iter().zip(genes) {
        let id = r["gene"].as_str().unwrap();
        assert_eq!(ours.gene_id, id);
        assert_close(
            ours.base_mean,
            num(&r["base_mean"]),
            &format!("{id} baseMean"),
        );
        // DESeq2 treats gene-wise estimates below 100 × minDisp (1e-6) as
        // "at the minimum" and excludes them from the trend fit; there the
        // value is dominated by floating-point cancellation, so only require
        // both implementations to land in that region.
        match num(&r["disp_gene_est"]) {
            Some(rd) if rd < 1e-6 => assert!(
                ours.disp_gene_est < 1e-6,
                "{id} dispGeneEst: ours {} vs R {rd}",
                ours.disp_gene_est
            ),
            rd => assert_close(ours.disp_gene_est, rd, &format!("{id} dispGeneEst")),
        }
        assert_close(ours.disp_fit, num(&r["disp_fit"]), &format!("{id} dispFit"));
        assert_close(
            ours.dispersion,
            num(&r["dispersion"]),
            &format!("{id} dispersion"),
        );
        let rtol = if num(&r["beta_iter"]) == Some(100.0) {
            OPTIM_RTOL
        } else {
            RTOL
        };
        let check =
            |o: f64, key: &str| assert_close_tol(o, num(&r[key]), &format!("{id} {key}"), rtol);
        check(ours.log2_fold_change, "log2_fold_change");
        check(ours.lfc_se, "lfc_se");
        check(ours.stat, "stat");
        check(ours.p_value, "pvalue");
        check(ours.padj, "padj");
    }
}

#[test]
fn wald_test_matches_r_deseq2() {
    let golden = load_golden("sim_1000x6");
    let (ids, counts, samples) = load_counts("deseq2_sim_1000x6.counts.tsv");
    let analysis = run_deseq2(&counts, &ids, &[0, 1, 2], &[3, 4, 5], &samples).unwrap();
    assert_wald_matches(&analysis, &golden);
    let significant = analysis.results.iter().filter(|r| r.padj < 0.1).count();
    assert_eq!(significant, 153, "DESeq2 calls 153 genes at padj < 0.1");
}

#[test]
fn batch_design_matches_r_deseq2() {
    let golden = load_golden("sim_1000x8_batch");
    let (ids, counts, samples) = load_counts("deseq2_sim_1000x8.counts.tsv");
    let group = [0, 0, 0, 0, 1, 1, 1, 1];
    let batch: Vec<u32> = golden["batch"]
        .as_array()
        .unwrap()
        .iter()
        .map(|b| b.as_str().unwrap().parse::<u32>().unwrap() - 1)
        .collect();
    // Columns: intercept, condition, batch2 — test the condition coefficient.
    let design = design_matrix_with_batch(&group, &batch);
    let analysis = run_deseq2_design(&counts, &ids, &design, &[0.0, 1.0, 0.0], &samples).unwrap();
    assert_wald_matches(&analysis, &golden);
}

#[test]
fn likelihood_ratio_test_matches_r_deseq2() {
    let golden = load_golden("sim_1000x6_lrt");
    let (ids, counts, samples) = load_counts("deseq2_sim_1000x6.counts.tsv");
    let full = design_matrix_two_group(&[0, 1, 2], &[3, 4, 5], 6);
    let reduced: Vec<Vec<f64>> = vec![vec![1.0]; 6];
    let results = lrt_test(&counts, &ids, &full, &reduced, &samples).unwrap();
    for (ours, r) in results.iter().zip(golden["genes"].as_array().unwrap()) {
        let id = r["gene"].as_str().unwrap();
        assert_eq!(ours.gene_id, id);
        assert_eq!(ours.df, 1);
        // The statistic is a difference of log-likelihoods: compare it on an
        // absolute scale (relative error is meaningless for values near 0).
        match num(&r["stat"]) {
            None => assert!(ours.stat.is_nan(), "{id}: R reports NA"),
            Some(r_stat) => assert!(
                (ours.stat - r_stat).abs() <= RTOL * r_stat.abs().max(1.0),
                "{id} LRT stat: ours {} vs R {r_stat}",
                ours.stat
            ),
        }
        assert_close(ours.p_value, num(&r["pvalue"]), &format!("{id} LRT pvalue"));
        assert_close(ours.padj, num(&r["padj"]), &format!("{id} LRT padj"));
    }
}

#[test]
fn small_matrix_size_factors_match_r_deseq2() {
    let golden = load_golden("small_5x6");
    let counts: Vec<Vec<f64>> = vec![
        vec![100.0, 120.0, 95.0, 200.0, 210.0, 195.0],
        vec![500.0, 480.0, 510.0, 1000.0, 990.0, 1010.0],
        vec![10.0, 8.0, 12.0, 10.0, 9.0, 11.0],
        vec![300.0, 320.0, 280.0, 150.0, 160.0, 140.0],
        vec![50.0, 55.0, 45.0, 100.0, 95.0, 105.0],
    ];
    let ids: Vec<String> = ["A", "B", "C", "D", "E"]
        .iter()
        .map(|g| format!("GENE_{g}"))
        .collect();
    let samples: Vec<String> = (1..=6).map(|i| format!("S{i}")).collect();
    let analysis = run_deseq2(&counts, &ids, &[0, 1, 2], &[3, 4, 5], &samples).unwrap();
    assert_size_factors(&analysis, &golden);
    for (ours, r) in analysis
        .results
        .iter()
        .zip(golden["genes"].as_array().unwrap())
    {
        assert_close(ours.base_mean, num(&r["base_mean"]), "baseMean");
    }
    // GENE_D is the only gene DESeq2 calls clearly down-regulated.
    let d = &analysis.results[3];
    assert!(d.log2_fold_change < -1.5 && d.p_value < 1e-20);
}
