//! Parser for Gene Matrix Transposed (.gmt) files and GMT-based pathway enrichment.
//!
//! GMT is a tab-delimited format used by MSigDB and other gene set databases:
//! ```text
//! PATHWAY_NAME\tDESCRIPTION\tGENE1\tGENE2\t...
//! ```
//! Each line is one gene set. Blank lines and comment lines (starting with `#`)
//! are silently skipped. Gene names are stored uppercased for consistent matching.
//!
//! The enrichment method uses Fisher's exact test (one-sided hypergeometric)
//! identical to [`crate::pathway::enrichment_analysis`], enabling direct
//! comparison of results across KEGG built-in pathways and custom GMT files.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use ahash::AHashSet;
use anyhow::Context;
use biomics_core::statistics::{benjamini_hochberg, hypergeometric_pvalue};
use serde::{Deserialize, Serialize};

use crate::pathway::EnrichmentResult;

// ── Public types ──────────────────────────────────────────────────────────────

/// One gene set parsed from a .gmt file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GmtPathway {
    /// Gene set name (column 1 of the GMT line).
    pub name: String,
    /// Descriptive text or URL (column 2 of the GMT line).
    pub description: String,
    /// Gene identifiers, uppercased for case-insensitive matching.
    pub genes: Vec<String>,
}

// ── Parsing ───────────────────────────────────────────────────────────────────

/// Parse a `.gmt` file into a `Vec<GmtPathway>`.
///
/// # Format
/// ```text
/// PATHWAY_NAME\tDESCRIPTION\tGENE1\tGENE2\t...
/// ```
/// - Column 1: pathway name (must be non-empty).
/// - Column 2: description or URL (may be empty — treated as an empty string).
/// - Columns 3+: gene identifiers (uppercased).
///
/// Lines beginning with `#` and blank lines are skipped.
///
/// # Errors
/// Returns an error if the file cannot be opened or if a non-comment, non-blank
/// line has fewer than 2 tab-separated columns.
pub fn parse_gmt(path: &Path) -> anyhow::Result<Vec<GmtPathway>> {
    let file =
        File::open(path).with_context(|| format!("parse_gmt: cannot open {}", path.display()))?;
    let reader = BufReader::new(file);
    parse_gmt_reader(reader)
}

/// Internal parser that works on any `BufRead` (enables in-memory testing).
pub(crate) fn parse_gmt_reader<R: BufRead>(reader: R) -> anyhow::Result<Vec<GmtPathway>> {
    let mut pathways = Vec::new();

    for (line_no, line_result) in reader.lines().enumerate() {
        let line =
            line_result.with_context(|| format!("parse_gmt: I/O error at line {line_no}"))?;

        // Skip blank lines and comment lines
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let mut fields = trimmed.splitn(3, '\t');

        let name = fields
            .next()
            .filter(|s| !s.is_empty())
            .ok_or_else(|| anyhow::anyhow!("parse_gmt: line {line_no} has no name field"))?
            .to_string();

        // Description may be empty (two adjacent tabs or missing)
        let description = fields.next().unwrap_or("").to_string();

        // Remaining genes (column 3 onward) — handle inline tab splitting
        let genes: Vec<String> = match fields.next() {
            None => Vec::new(),
            Some(rest) => rest
                .split('\t')
                .map(|g| g.trim())
                .filter(|g| !g.is_empty())
                .map(|g| g.to_uppercase())
                .collect(),
        };

        pathways.push(GmtPathway {
            name,
            description,
            genes,
        });
    }

    Ok(pathways)
}

// ── Format conversion ─────────────────────────────────────────────────────────

/// Convert parsed GMT pathways into the owned-String tuple format expected by
/// `gsea_preranked` and similar functions.
///
/// Returns `Vec<(id, name, genes)>` where `id == name` (GMT files do not have
/// a separate ID column) and genes are already uppercased.
pub fn gmt_to_query_format(pathways: &[GmtPathway]) -> Vec<(String, String, Vec<String>)> {
    pathways
        .iter()
        .map(|p| (p.name.clone(), p.name.clone(), p.genes.clone()))
        .collect()
}

// ── Enrichment analysis ───────────────────────────────────────────────────────

/// Run Fisher's exact test enrichment for custom GMT pathways.
///
/// Mirrors [`crate::pathway::enrichment_analysis`] but works with
/// [`GmtPathway`] entries instead of the built-in KEGG table.
///
/// ## Method
/// - Background universe: union of all genes appearing in *any* pathway in
///   `gmt_pathways`.
/// - Test: one-sided (upper tail) hypergeometric p-value.
/// - FDR correction: Benjamini-Hochberg across all results with
///   `overlap >= min_overlap`.
/// - Jaccard-like score: `overlap / sqrt(pathway_size × query_size)`.
///
/// Results are sorted by `padj` ascending.
///
/// # Arguments
/// - `query_genes`: the foreground gene set (e.g. DE genes). Case-insensitive.
/// - `gmt_pathways`: parsed GMT pathways.
/// - `min_overlap`: pathways with fewer overlapping genes are excluded.
pub fn gmt_enrichment_analysis(
    query_genes: &[String],
    gmt_pathways: &[GmtPathway],
    min_overlap: usize,
) -> Vec<EnrichmentResult> {
    if query_genes.is_empty() || gmt_pathways.is_empty() {
        return Vec::new();
    }

    let query_set: AHashSet<String> = query_genes.iter().map(|g| g.to_uppercase()).collect();

    // Build background universe: union of all genes in all GMT pathways
    let background: AHashSet<String> = gmt_pathways
        .iter()
        .flat_map(|p| p.genes.iter().cloned())
        .collect();
    let bg_size = background.len();

    if bg_size == 0 {
        return Vec::new();
    }

    let mut results: Vec<EnrichmentResult> = gmt_pathways
        .iter()
        .filter_map(|pathway| {
            // Pathway genes are already uppercased by the parser
            let pathway_set: AHashSet<&str> = pathway.genes.iter().map(|g| g.as_str()).collect();
            let overlap = query_set
                .iter()
                .filter(|g| pathway_set.contains(g.as_str()))
                .count();

            if overlap < min_overlap {
                return None;
            }

            let pathway_size = pathway.genes.len();
            let query_size = query_set.len();

            let score = overlap as f64 / ((pathway_size as f64) * (query_size as f64)).sqrt();

            let p_value = hypergeometric_pvalue(overlap, query_size, pathway_size, bg_size);

            Some(EnrichmentResult {
                pathway_id: pathway.name.clone(),
                pathway_name: pathway.name.clone(),
                overlap,
                pathway_size,
                query_size,
                score,
                p_value,
                padj: f64::NAN, // filled below
            })
        })
        .collect();

    if results.is_empty() {
        return Vec::new();
    }

    // BH FDR correction
    let pvals: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    let padj_vals = benjamini_hochberg(&pvals);
    for (r, padj) in results.iter_mut().zip(padj_vals) {
        r.padj = padj;
    }

    // Sort by padj ascending
    results.sort_unstable_by(|a, b| {
        a.padj
            .partial_cmp(&b.padj)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SAMPLE_GMT: &str = "\
# This is a comment line — should be skipped

HALLMARK_HYPOXIA\thttps://www.gsea-msigdb.org/HALLMARK_HYPOXIA\tVEGFA\tSLC2A1\tPGK1\tENO1\tLDHA\tPKM\tHIF1A
HALLMARK_GLYCOLYSIS\tGlycolysis gene set\tHK2\tPFK1\tALDOA\tGAPDH\tPGK1\tENO1\tPKM\tLDHA\tG6PD
CELL_CYCLE_CORE\t\tCDK1\tCDK2\tCCNA2\tCCNB1\tMCM2\tPCNA
";

    #[test]
    fn test_parse_gmt_in_memory() {
        let cursor = Cursor::new(SAMPLE_GMT);
        let pathways = parse_gmt_reader(cursor).expect("parse should succeed");

        // 3 non-comment, non-blank lines → 3 pathways
        assert_eq!(
            pathways.len(),
            3,
            "expected 3 pathways, got {}",
            pathways.len()
        );

        let hypoxia = &pathways[0];
        assert_eq!(hypoxia.name, "HALLMARK_HYPOXIA");
        assert_eq!(
            hypoxia.genes.len(),
            7,
            "HALLMARK_HYPOXIA should have 7 genes"
        );
        // Genes must be uppercased
        assert!(hypoxia.genes.contains(&"VEGFA".to_string()));
        assert!(hypoxia.genes.contains(&"HIF1A".to_string()));

        let glycolysis = &pathways[1];
        assert_eq!(glycolysis.name, "HALLMARK_GLYCOLYSIS");
        assert_eq!(
            glycolysis.genes.len(),
            9,
            "HALLMARK_GLYCOLYSIS should have 9 genes"
        );

        // Third pathway has an empty description field
        let cc = &pathways[2];
        assert_eq!(cc.name, "CELL_CYCLE_CORE");
        assert_eq!(
            cc.description, "",
            "empty description should parse as empty string"
        );
        assert_eq!(cc.genes.len(), 6);
    }

    #[test]
    fn test_gmt_enrichment_significant() {
        let cursor = Cursor::new(SAMPLE_GMT);
        let pathways = parse_gmt_reader(cursor).unwrap();

        // Query contains all glycolysis genes + a few extra → should strongly enrich HALLMARK_GLYCOLYSIS
        let query = vec![
            "HK2".to_string(),
            "PFK1".to_string(),
            "ALDOA".to_string(),
            "GAPDH".to_string(),
            "PGK1".to_string(),
            "ENO1".to_string(),
            "PKM".to_string(),
            "LDHA".to_string(),
            "G6PD".to_string(),
            "FAKE_GENE_1".to_string(),
            "FAKE_GENE_2".to_string(),
        ];

        let results = gmt_enrichment_analysis(&query, &pathways, 1);

        assert!(!results.is_empty(), "should return enrichment results");

        // The top result should be HALLMARK_GLYCOLYSIS with high overlap
        let top = &results[0];
        assert!(
            top.overlap >= 8,
            "expected >= 8 overlapping glycolysis genes, got {}",
            top.overlap
        );
        assert!(
            top.p_value < 0.05,
            "expected significant enrichment for glycolysis, p={}",
            top.p_value
        );

        // padj should be filled and not NaN
        for r in &results {
            assert!(
                !r.padj.is_nan(),
                "padj should not be NaN for {}",
                r.pathway_name
            );
        }
    }

    #[test]
    fn test_gmt_enrichment_min_overlap_filter() {
        let cursor = Cursor::new(SAMPLE_GMT);
        let pathways = parse_gmt_reader(cursor).unwrap();

        // Query overlaps CELL_CYCLE_CORE with only 1 gene
        let query = vec!["CDK1".to_string(), "UNRELATED_GENE".to_string()];

        // With min_overlap=2 the cell cycle pathway should be excluded
        let results = gmt_enrichment_analysis(&query, &pathways, 2);
        for r in &results {
            assert!(
                r.overlap >= 2,
                "result {} has overlap {} < min_overlap 2",
                r.pathway_name,
                r.overlap
            );
        }
    }

    #[test]
    fn test_gmt_to_query_format() {
        let cursor = Cursor::new(SAMPLE_GMT);
        let pathways = parse_gmt_reader(cursor).unwrap();
        let converted = gmt_to_query_format(&pathways);

        assert_eq!(converted.len(), pathways.len());
        for (i, (id, name, genes)) in converted.iter().enumerate() {
            assert_eq!(id, &pathways[i].name);
            assert_eq!(name, &pathways[i].name);
            assert_eq!(genes, &pathways[i].genes);
        }
    }
}
