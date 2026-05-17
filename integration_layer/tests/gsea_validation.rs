//! Mathematical validation of GSEA enrichment-score properties.
//!
//! These tests verify the theoretical bounds established by Subramanian et al. 2005
//! (PNAS 102(43):15545–15550). They are not R-specific — the properties are provable
//! from the definition of the KS enrichment statistic and do not depend on any
//! external reference package.
//!
//! Properties verified:
//!  1. A perfectly top-loaded pathway (all hits at the top of a 200-gene ranked list)
//!     produces ES > 0.8, NES > 1.5, and p < 0.05.
//!  2. A perfectly bottom-loaded pathway produces ES < -0.8.
//!  3. A randomly distributed pathway has |ES| strictly less than the top-loaded ES.
//!  4. The top-loaded NES exceeds 1.5 (strongly enriched signal).

use integration_layer::gsea_preranked;

// Number of permutations for null distribution. 500 gives stable empirical p-values
// without making the test slow in CI.
const N_PERM: usize = 500;

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Build a ranked gene list of length `n` with linearly decreasing metric scores.
/// Gene identifiers are "GENE0000", "GENE0001", …, "GENE0199".
fn make_ranked_list(n: usize) -> Vec<(String, f64)> {
    (0..n)
        .map(|i| {
            let id = format!("GENE{i:04}");
            let metric = (n as f64) - (i as f64); // highest score at rank 0
            (id, metric)
        })
        .collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[test]
fn test_gsea_matches_known_properties() {
    let n = 200_usize;
    let k = 20_usize; // pathway size — top k genes

    let ranked = make_ranked_list(n);

    // ── Pathway definitions ───────────────────────────────────────────────────

    // Top-loaded: genes at ranks 0..k (highest scores) — should be strongly enriched.
    let top_genes: Vec<&str> = ranked[..k].iter().map(|(g, _)| g.as_str()).collect();

    // Bottom-loaded: genes at ranks (n-k)..n (lowest scores) — should be strongly depleted.
    let bottom_genes: Vec<&str> = ranked[(n - k)..].iter().map(|(g, _)| g.as_str()).collect();

    // Random/uniform pathway: genes evenly spaced through the ranked list (every 10th gene).
    // This approximates a null pathway with no systematic enrichment.
    let random_genes: Vec<&str> = ranked
        .iter()
        .step_by(n / k)
        .map(|(g, _)| g.as_str())
        .collect();

    // ── Run GSEA for each pathway separately so results are independent ───────

    let top_pathways = vec![("PW_TOP", "Top-loaded pathway", top_genes.as_slice())];
    let bottom_pathways = vec![("PW_BOT", "Bottom-loaded pathway", bottom_genes.as_slice())];
    let random_pathways = vec![("PW_RND", "Random/uniform pathway", random_genes.as_slice())];

    let top_results = gsea_preranked(&ranked, &top_pathways, 5, 500, N_PERM);
    let bottom_results = gsea_preranked(&ranked, &bottom_pathways, 5, 500, N_PERM);
    let random_results = gsea_preranked(&ranked, &random_pathways, 5, 500, N_PERM);

    assert_eq!(
        top_results.len(),
        1,
        "top-loaded pathway must return one result"
    );
    assert_eq!(
        bottom_results.len(),
        1,
        "bottom-loaded pathway must return one result"
    );
    assert_eq!(
        random_results.len(),
        1,
        "random pathway must return one result"
    );

    let top = &top_results[0];
    let bottom = &bottom_results[0];
    let random = &random_results[0];

    // ── Print summary table ───────────────────────────────────────────────────
    println!("\n=== GSEA property validation ===");
    println!(
        "{:<20}  {:>8}  {:>8}  {:>8}  {:>6}",
        "Pathway", "ES", "NES", "p-value", "Genes"
    );
    println!("{}", "-".repeat(60));
    for r in [top, bottom, random] {
        println!(
            "{:<20}  {:>8.4}  {:>8.4}  {:>8.4}  {:>6}",
            r.pathway_name, r.es, r.nes, r.p_value, r.n_genes_pathway
        );
    }
    println!();

    // ── Property 1: top-loaded pathway is strongly enriched ──────────────────
    println!("Property 1: top-loaded ES > 0.8");
    assert!(
        top.es > 0.8,
        "Top-loaded pathway ES should be > 0.8 (Subramanian 2005 Table 1), got ES={}",
        top.es
    );
    println!("  ES = {:.4}  PASS", top.es);

    println!("Property 1b: top-loaded p < 0.05");
    assert!(
        top.p_value < 0.05,
        "Top-loaded pathway p-value should be < 0.05, got p={}",
        top.p_value
    );
    println!("  p = {:.4}  PASS", top.p_value);

    // ── Property 2: bottom-loaded pathway is strongly depleted ───────────────
    println!("Property 2: bottom-loaded ES < -0.8");
    assert!(
        bottom.es < -0.8,
        "Bottom-loaded pathway ES should be < -0.8, got ES={}",
        bottom.es
    );
    println!("  ES = {:.4}  PASS", bottom.es);

    // ── Property 3: random pathway |ES| < top-loaded ES ─────────────────────
    println!("Property 3: random |ES| < top-loaded ES");
    assert!(
        random.es.abs() < top.es,
        "Random pathway |ES| ({}) should be less than top-loaded ES ({})",
        random.es.abs(),
        top.es
    );
    println!(
        "  |random ES| = {:.4}  <  top ES = {:.4}  PASS",
        random.es.abs(),
        top.es
    );

    // ── Property 4: NES > 1.5 for a strongly concentrated pathway ────────────
    println!("Property 4: top-loaded NES > 1.5");
    assert!(
        top.nes > 1.5,
        "Top-loaded NES should exceed 1.5 (strong enrichment), got NES={}",
        top.nes
    );
    println!("  NES = {:.4}  PASS", top.nes);

    // ── Sanity: running-sum length matches ranked-list length ─────────────────
    assert_eq!(
        top.running_sum.len(),
        n,
        "running_sum length should equal ranked list length"
    );

    // ── Sanity: leading-edge genes are non-empty for significant pathway ───────
    assert!(
        !top.leading_edge.is_empty(),
        "leading edge should be non-empty for the top-loaded pathway"
    );

    println!("\nAll GSEA property assertions passed.");
}
