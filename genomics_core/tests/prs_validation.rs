//! Validate polygenic risk scores against PLINK 2 `--score`.
//!
//! `test_data/validation/prs_1kg_*.vcf` hold the genotypes of 1000 Genomes
//! samples HG00096 and NA12878 at every locus of `data/prs_catalog.tsv`
//! (GRCh38 release, including homozygous-reference calls). The reference sums
//! in `tests/golden/plink2/` come from `validation/plink2_prs.sh`, which scores
//! the same weights (ln OR of the risk allele) with PLINK 2. PLINK prints six
//! significant digits, hence the tolerance.

use std::path::PathBuf;

use genomics_core::analyze_vcf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join(rel)
}

fn check_sample(sample: &str) {
    let sscore = std::fs::read_to_string(repo(&format!(
        "tests/golden/plink2/prs_1kg_{sample}.sscore"
    )))
    .unwrap();
    let mut lines = sscore.lines().filter(|l| !l.starts_with("# "));
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let values: Vec<&str> = lines.next().unwrap().split('\t').collect();

    let summary = analyze_vcf(
        &repo(&format!("test_data/validation/prs_1kg_{sample}.vcf")),
        None,
    )
    .unwrap();
    assert_eq!(summary.prs_scores.len(), header.len() - 1);
    for (col, value) in header.iter().zip(&values).skip(1) {
        let disease = col.trim_end_matches("_SUM").replace('_', " ");
        let plink: f64 = value.parse().unwrap();
        let ours = summary
            .prs_scores
            .iter()
            .find(|p| p.disease == disease)
            .unwrap_or_else(|| panic!("{disease} missing"));
        assert!(
            (ours.raw_score - plink).abs() <= 5e-6 * plink.abs().max(1.0),
            "{sample} {disease}: ours {} vs PLINK 2 {plink}",
            ours.raw_score
        );
        assert_eq!(
            ours.n_variants_matched, ours.n_variants_total,
            "{sample} {disease}"
        );
        assert!(ours.z_score.is_finite());
    }
}

#[test]
fn prs_matches_plink2_for_hg00096() {
    check_sample("HG00096");
}

#[test]
fn prs_matches_plink2_for_na12878() {
    check_sample("NA12878");
}
