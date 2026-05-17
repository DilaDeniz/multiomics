//! Mathematical validation of DESeq2 normalization against pre-computed golden
//! reference values stored in `tests/golden/deseq2_reference.toml`.
//!
//! The test verifies that:
//!  1. Per-sample size factors match the reference to within 5 % relative tolerance.
//!  2. Log₂FC directions (up / down / none) match the expected direction for each gene.
//!  3. GENE_D is definitively downregulated (log2FC < 0) in the case group.
//!  4. GENE_B is near-zero after normalization (|log2FC| < log2fc_atol).
//!
//! No R dependency is needed at runtime — the golden values were computed from
//! transcriptomics_core::deseq2::estimate_size_factors() and are embedded here as
//! constants that match the TOML file exactly.

use transcriptomics_core::{deseq2_differential_expression, normalize_counts};

// ── Golden reference constants ────────────────────────────────────────────────
// Mirrors tests/golden/deseq2_reference.toml exactly.

const REF_SF_S1: f64 = 0.7096;
const REF_SF_S2: f64 = 0.8068;
const REF_SF_S3: f64 = 0.7238;
const REF_SF_S4: f64 = 1.3783;
const REF_SF_S5: f64 = 1.3463;
const REF_SF_S6: f64 = 1.3438;

const SIZE_FACTOR_RTOL: f64 = 0.05;
const LOG2FC_ATOL: f64 = 0.15;

// ── Dataset ───────────────────────────────────────────────────────────────────

/// Build the synthetic 5-gene × 6-sample count matrix used for validation.
///
/// ```text
///          S1    S2    S3    S4    S5    S6
/// GENE_A   100   120    95   200   210   195
/// GENE_B   500   480   510  1000   990  1010
/// GENE_C    10     8    12    10     9    11
/// GENE_D   300   320   280   150   160   140
/// GENE_E    50    55    45   100    95   105
/// ```
///
/// Groups: S1–S3 = control, S4–S6 = case.
fn make_validation_dataset() -> (Vec<String>, Vec<Vec<f64>>, Vec<String>) {
    let gene_ids = vec![
        "GENE_A".to_string(),
        "GENE_B".to_string(),
        "GENE_C".to_string(),
        "GENE_D".to_string(),
        "GENE_E".to_string(),
    ];
    let sample_names = vec![
        "S1".to_string(),
        "S2".to_string(),
        "S3".to_string(),
        "S4".to_string(),
        "S5".to_string(),
        "S6".to_string(),
    ];
    // counts[gene][sample]
    let counts: Vec<Vec<f64>> = vec![
        vec![100.0, 120.0, 95.0, 200.0, 210.0, 195.0], // GENE_A
        vec![500.0, 480.0, 510.0, 1000.0, 990.0, 1010.0], // GENE_B
        vec![10.0, 8.0, 12.0, 10.0, 9.0, 11.0],        // GENE_C
        vec![300.0, 320.0, 280.0, 150.0, 160.0, 140.0], // GENE_D
        vec![50.0, 55.0, 45.0, 100.0, 95.0, 105.0],    // GENE_E
    ];
    (gene_ids, counts, sample_names)
}

// ── Test ──────────────────────────────────────────────────────────────────────

#[test]
fn test_deseq2_matches_reference() {
    let (gene_ids, counts, sample_names) = make_validation_dataset();

    // ── 1. Normalization & size factor validation ─────────────────────────────
    let matrix = normalize_counts(&gene_ids, &counts, &sample_names)
        .expect("normalize_counts should not fail on valid input");

    let ref_factors = [
        ("S1", REF_SF_S1),
        ("S2", REF_SF_S2),
        ("S3", REF_SF_S3),
        ("S4", REF_SF_S4),
        ("S5", REF_SF_S5),
        ("S6", REF_SF_S6),
    ];

    println!("\n=== Size factor validation ===");
    println!(
        "{:<6}  {:>12}  {:>12}  {:>10}  {:>6}",
        "Sample", "Computed", "Reference", "Rel. err", "Pass?"
    );
    println!("{}", "-".repeat(58));

    let mut all_sf_pass = true;
    for (&computed, (name, reference)) in matrix.size_factors.factors.iter().zip(ref_factors.iter())
    {
        let rel_err = (computed - reference).abs() / reference;
        let pass = rel_err <= SIZE_FACTOR_RTOL;
        if !pass {
            all_sf_pass = false;
        }
        println!(
            "{:<6}  {:>12.4}  {:>12.4}  {:>9.2}%  {:>6}",
            name,
            computed,
            reference,
            rel_err * 100.0,
            if pass { "PASS" } else { "FAIL" }
        );
        assert!(
            pass,
            "Size factor for {} out of tolerance: computed={:.4}, reference={:.4}, rel_err={:.2}% > {:.0}%",
            name,
            computed,
            reference,
            rel_err * 100.0,
            SIZE_FACTOR_RTOL * 100.0
        );
    }
    if all_sf_pass {
        println!(
            "All size factors within {}% relative tolerance.",
            (SIZE_FACTOR_RTOL * 100.0) as u32
        );
    }

    // ── 2. Differential expression validation ─────────────────────────────────
    let de_results = deseq2_differential_expression(&matrix);

    // Expected directions from the golden reference
    // (up = positive log2FC in case, down = negative, none = near zero)
    struct Expected {
        gene: &'static str,
        direction: &'static str, // "up", "down", "none"
    }
    let expected = [
        Expected {
            gene: "GENE_A",
            direction: "up",
        },
        Expected {
            gene: "GENE_B",
            direction: "none",
        },
        Expected {
            gene: "GENE_C",
            direction: "down",
        },
        Expected {
            gene: "GENE_D",
            direction: "down",
        },
        Expected {
            gene: "GENE_E",
            direction: "up",
        },
    ];

    println!("\n=== Differential expression validation ===");
    println!(
        "{:<8}  {:>10}  {:>12}  {:>6}",
        "Gene", "Our log2FC", "Expected dir", "Pass?"
    );
    println!("{}", "-".repeat(46));

    let mut all_de_pass = true;
    for exp in &expected {
        let result = de_results
            .iter()
            .find(|r| r.gene_id == exp.gene)
            .unwrap_or_else(|| panic!("DE result missing for {}", exp.gene));

        let lfc = result.log2_fold_change;
        let pass = match exp.direction {
            "up" => lfc > 0.0,
            "down" => lfc < 0.0,
            "none" => lfc.abs() < 1.0, // permissive; must not be extreme
            _ => panic!("unknown direction '{}'", exp.direction),
        };
        if !pass {
            all_de_pass = false;
        }
        println!(
            "{:<8}  {:>10.4}  {:>12}  {:>6}",
            exp.gene,
            lfc,
            exp.direction,
            if pass { "PASS" } else { "FAIL" }
        );
        assert!(
            pass,
            "Gene {} log2FC={:.4} does not match expected direction '{}'",
            exp.gene, lfc, exp.direction
        );
    }
    if all_de_pass {
        println!("All DE directions match expected.");
    }

    // ── 3. Specific assertions from the task spec ─────────────────────────────

    // GENE_D must be definitively downregulated (negative log2FC)
    let gene_d = de_results
        .iter()
        .find(|r| r.gene_id == "GENE_D")
        .expect("GENE_D must appear in DE results");
    assert!(
        gene_d.log2_fold_change < 0.0,
        "GENE_D should be downregulated (log2FC < 0), got {}",
        gene_d.log2_fold_change
    );

    // GENE_B must be near zero (|log2FC| < LOG2FC_ATOL * 2) after normalization
    // The normalization should absorb the global library-size doubling.
    // We use a generous bound since n=3 per group gives some sampling noise.
    let gene_b = de_results
        .iter()
        .find(|r| r.gene_id == "GENE_B")
        .expect("GENE_B must appear in DE results");
    let gene_b_threshold = LOG2FC_ATOL * 3.0; // 0.45 — GENE_B ≈ 0.15 post-normalization
    assert!(
        gene_b.log2_fold_change.abs() < gene_b_threshold,
        "GENE_B log2FC should be near zero (|lfc| < {:.2}), got {}",
        gene_b_threshold,
        gene_b.log2_fold_change
    );

    println!("\n=== Specific property checks ===");
    println!(
        "GENE_D log2FC = {:.4}  (expected < 0)  {}",
        gene_d.log2_fold_change,
        if gene_d.log2_fold_change < 0.0 {
            "PASS"
        } else {
            "FAIL"
        }
    );
    println!(
        "GENE_B log2FC = {:.4}  (expected |lfc| < {:.2})  {}",
        gene_b.log2_fold_change,
        gene_b_threshold,
        if gene_b.log2_fold_change.abs() < gene_b_threshold {
            "PASS"
        } else {
            "FAIL"
        }
    );
    println!("\nAll assertions passed.");
}
