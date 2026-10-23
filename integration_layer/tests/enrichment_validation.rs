//! Validate GMT over-representation analysis against gseapy.
//!
//! `tests/golden/gseapy/enrichment.json` is produced by
//! `validation/gseapy_reference.py` (gseapy `enrich` with the union of all GMT
//! genes as background) from `test_data/validation/enrichment_sets.gmt` and
//! `enrichment_query.txt`. The data exercise query genes absent from every
//! set (dropped from the query size) and a set that lists genes twice
//! (counted once).

use std::path::PathBuf;

use integration_layer::{gmt_enrichment_analysis, parse_gmt};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../test_data/validation")
        .join(name)
}

#[test]
fn gmt_enrichment_matches_gseapy() {
    let golden: serde_json::Value =
        serde_json::from_str(include_str!("../../tests/golden/gseapy/enrichment.json")).unwrap();
    let pathways = parse_gmt(&data("enrichment_sets.gmt")).unwrap();
    let query: Vec<String> = std::fs::read_to_string(data("enrichment_query.txt"))
        .unwrap()
        .split_whitespace()
        .map(String::from)
        .collect();

    let results = gmt_enrichment_analysis(&query, &pathways, 1);
    let terms = golden["terms"].as_array().unwrap();
    assert_eq!(results.len(), terms.len(), "number of reported terms");

    let query_in_bg = golden["query_in_background"].as_u64().unwrap() as usize;
    for t in terms {
        let name = t["term"].as_str().unwrap();
        let ours = results
            .iter()
            .find(|r| r.pathway_name == name)
            .unwrap_or_else(|| panic!("{name} missing"));
        assert_eq!(
            ours.overlap as u64,
            t["overlap"].as_u64().unwrap(),
            "{name} overlap"
        );
        assert_eq!(
            ours.pathway_size as u64,
            t["set_size"].as_u64().unwrap(),
            "{name} size"
        );
        assert_eq!(ours.query_size, query_in_bg, "{name} query size");
        for (got, key) in [(ours.p_value, "p_value"), (ours.padj, "padj")] {
            let expected = t[key].as_f64().unwrap();
            assert!(
                ((got - expected) / expected).abs() < 1e-9,
                "{name} {key}: ours {got} vs gseapy {expected}"
            );
        }
    }
}
