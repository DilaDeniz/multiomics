//! Python bindings for Multiomics via PyO3.
//!
//! Every function returns a JSON-serialised string so callers can use any
//! Python JSON parser and avoid round-trip overhead on the Rust side.
//!
//! # Build
//! ```sh
//! pip install maturin
//! maturin develop          # editable install for development
//! maturin build --release  # produce a wheel
//! ```
//!
//! # Usage
//! ```python
//! import multiomics_py as bmo
//! import json
//!
//! genomics = json.loads(bmo.analyze_vcf("variants.vcf"))
//! transcr  = json.loads(bmo.analyze_tsv("expr.tsv"))
//! epigen   = json.loads(bmo.analyze_bed("meth.bed"))
//! result   = json.loads(bmo.run_full_pipeline("variants.vcf", "expr.tsv", "meth.bed", False))
//! ```

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::path::Path;

// ── helpers ──────────────────────────────────────────────────────────────────

/// Convert an [`anyhow::Error`] to a Python `RuntimeError`.
#[inline]
fn to_py_err(e: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(format!("{:#}", e))
}

/// Serialize `value` to a JSON string, mapping any error to `PyRuntimeError`.
#[inline]
fn to_json<T: serde::Serialize>(value: &T) -> PyResult<String> {
    serde_json::to_string(value).map_err(|e| to_py_err(anyhow::anyhow!(e)))
}

// ── individual modality functions ─────────────────────────────────────────────

/// Parse a VCF file and return a JSON-encoded `GenomicsSummary`.
///
/// # Arguments
/// * `path` – path to the VCF file (plain or gzip-compressed)
///
/// # Returns
/// JSON string — deserialise with `json.loads()` in Python.
///
/// # Errors
/// Raises `RuntimeError` if the file cannot be opened or parsed.
#[pyfunction]
fn analyze_vcf(path: &str) -> PyResult<String> {
    let summary = genomics_core::analyze_vcf(Path::new(path), None).map_err(to_py_err)?;
    to_json(&summary)
}

/// Parse a TSV expression matrix and return a JSON-encoded `TranscriptomicsSummary`.
///
/// The TSV must have a header row with sample names followed by rows of the
/// form `gene_id\t<tpm1>\t<tpm2>\t...`.
///
/// # Returns
/// JSON string — deserialise with `json.loads()` in Python.
///
/// # Errors
/// Raises `RuntimeError` if the file cannot be opened or parsed.
#[pyfunction]
fn analyze_tsv(path: &str) -> PyResult<String> {
    let summary = transcriptomics_core::analyze_tsv(Path::new(path), None).map_err(to_py_err)?;
    to_json(&summary)
}

/// Parse a BED methylation file and return a JSON-encoded `EpigenomicsSummary`.
///
/// Expected format: `chrom start end name score strand` where `score` is the
/// methylation percentage on the 0–1000 scale (ENCODE bisulfite BED).
///
/// # Returns
/// JSON string — deserialise with `json.loads()` in Python.
///
/// # Errors
/// Raises `RuntimeError` if the file cannot be opened or parsed.
#[pyfunction]
fn analyze_bed(path: &str) -> PyResult<String> {
    let summary = epigenomics_core::analyze_bed(Path::new(path), None).map_err(to_py_err)?;
    to_json(&summary)
}

// ── pathway enrichment ────────────────────────────────────────────────────────

/// Run hypergeometric pathway enrichment against the built-in KEGG pathway table.
///
/// # Arguments
/// * `genes`       – list of gene symbols (case-insensitive)
/// * `min_overlap` – minimum number of query genes that must overlap a pathway
///                   for it to appear in the results (default 2)
///
/// # Returns
/// JSON-encoded `Vec<EnrichmentResult>` sorted by BH-adjusted p-value ascending.
///
/// # Errors
/// Raises `RuntimeError` on serialization failure (extremely unlikely).
#[pyfunction]
#[pyo3(signature = (genes, min_overlap = 2))]
fn enrich_pathways(genes: Vec<String>, min_overlap: usize) -> PyResult<String> {
    let results = integration_layer::enrichment_analysis(&genes, min_overlap);
    to_json(&results)
}

// ── GSEA pre-ranked ───────────────────────────────────────────────────────────

/// Run a lightweight pre-ranked GSEA against the built-in KEGG pathway table.
///
/// This is a simplified leading-edge enrichment score (Kolmogorov-Smirnov
/// style) rather than the full Subramanian 2005 algorithm. It is suitable for
/// exploratory analysis; use a dedicated tool (GSEApy, fgsea) for publication.
///
/// # Arguments
/// * `ranked` – list of `(gene_symbol, metric)` tuples sorted **descending**
///              by metric (e.g. log₂FC × −log₁₀(p))
/// * `min_size` – minimum pathway gene count to test (default 5)
/// * `n_perm`   – number of permutations for p-value estimation (default 1000)
///
/// # Returns
/// JSON-encoded list of `GseaResult` objects sorted by NES descending.
///
/// # Errors
/// Raises `RuntimeError` if the ranked list is empty or serialization fails.
#[pyfunction]
#[pyo3(signature = (ranked, min_size = 5, n_perm = 1000))]
fn gsea_preranked(ranked: Vec<(String, f64)>, min_size: usize, n_perm: usize) -> PyResult<String> {
    if ranked.is_empty() {
        return Err(PyRuntimeError::new_err("ranked list is empty"));
    }

    let results = run_gsea_preranked(&ranked, min_size, n_perm);
    to_json(&results)
}

// ── full pipeline ─────────────────────────────────────────────────────────────

/// Run the complete multi-omics pipeline and return a JSON-encoded
/// `IntegrationSummary`.
///
/// This function parses all three input files concurrently (VCF, TSV, BED),
/// runs per-modality analyses, and then runs cross-modality integration
/// including PCA, Pearson correlation, pathway enrichment, and biological
/// insight derivation.
///
/// # Arguments
/// * `vcf`   – path to the VCF variant file
/// * `tsv`   – path to the expression matrix TSV
/// * `bed`   – path to the BED methylation file
/// * `no_ml` – when `true`, skip PCA and correlation (returns identity matrix)
///
/// # Returns
/// JSON-encoded `IntegrationSummary`.
///
/// # Errors
/// Raises `RuntimeError` if any input file cannot be parsed or integration
/// analysis fails.
#[pyfunction]
#[pyo3(signature = (vcf, tsv, bed, no_ml = false))]
fn run_full_pipeline(vcf: &str, tsv: &str, bed: &str, no_ml: bool) -> PyResult<String> {
    let vcf_path = Path::new(vcf);
    let tsv_path = Path::new(tsv);
    let bed_path = Path::new(bed);

    // Run the three modality analyses concurrently.
    let (g_res, t_res, e_res) = std::thread::scope(|s| {
        let gh = s.spawn(|| genomics_core::analyze_vcf(vcf_path, None));
        let th = s.spawn(|| transcriptomics_core::analyze_tsv(tsv_path, None));
        let eh = s.spawn(|| epigenomics_core::analyze_bed(bed_path, None));
        (gh.join(), th.join(), eh.join())
    });

    let genomics = g_res
        .map_err(|_| PyRuntimeError::new_err("genomics thread panicked"))?
        .map_err(to_py_err)?;
    let transcr = t_res
        .map_err(|_| PyRuntimeError::new_err("transcriptomics thread panicked"))?
        .map_err(to_py_err)?;
    let epigen = e_res
        .map_err(|_| PyRuntimeError::new_err("epigenomics thread panicked"))?
        .map_err(to_py_err)?;

    let integration = integration_layer::run_integration(&genomics, &transcr, &epigen, no_ml)
        .map_err(to_py_err)?;

    to_json(&integration)
}

/// Deserialize pre-computed JSON modality summaries and run integration.
///
/// This is useful when you already have summaries from previous runs and only
/// want to re-run the integration layer.
///
/// # Arguments
/// * `genomics_json`  – JSON string produced by `analyze_vcf`
/// * `transcr_json`   – JSON string produced by `analyze_tsv`
/// * `epigen_json`    – JSON string produced by `analyze_bed`
///
/// # Returns
/// JSON-encoded `IntegrationSummary`.
///
/// # Errors
/// Raises `RuntimeError` if any JSON is malformed or integration fails.
#[pyfunction]
fn run_integration_from_json(
    genomics_json: &str,
    transcr_json: &str,
    epigen_json: &str,
) -> PyResult<String> {
    let genomics: genomics_core::GenomicsSummary =
        serde_json::from_str(genomics_json).map_err(|e| to_py_err(anyhow::anyhow!(e)))?;
    let transcr: transcriptomics_core::TranscriptomicsSummary =
        serde_json::from_str(transcr_json).map_err(|e| to_py_err(anyhow::anyhow!(e)))?;
    let epigen: epigenomics_core::EpigenomicsSummary =
        serde_json::from_str(epigen_json).map_err(|e| to_py_err(anyhow::anyhow!(e)))?;

    let integration = integration_layer::run_integration(&genomics, &transcr, &epigen, false)
        .map_err(to_py_err)?;
    to_json(&integration)
}

// ── GSEA implementation ───────────────────────────────────────────────────────

/// A single GSEA result for one pathway.
#[derive(serde::Serialize)]
struct GseaResult {
    pathway_id: String,
    pathway_name: String,
    /// Normalized enrichment score (positive = enriched in top of ranked list).
    nes: f64,
    /// Estimated p-value from permutation test.
    p_value: f64,
    pathway_size: usize,
    /// Gene symbols at the leading edge.
    leading_edge: Vec<String>,
}

/// Compute a KS-style enrichment score for `pathway_genes` against a ranked
/// gene list.  Returns `(ES, leading_edge_genes)`.
fn enrichment_score(
    ranked: &[(String, f64)],
    pathway_set: &std::collections::HashSet<String>,
    n_total: usize,
    n_hit: usize,
) -> (f64, Vec<String>) {
    if n_hit == 0 {
        return (0.0, Vec::new());
    }

    let n_miss = n_total - n_hit;
    let hit_weight = if n_hit == 0 { 1.0 } else { 1.0 / n_hit as f64 };
    let miss_weight = if n_miss == 0 {
        0.0
    } else {
        1.0 / n_miss as f64
    };

    let mut running = 0.0_f64;
    let mut max_dev = 0.0_f64;
    let mut max_pos = 0_usize;

    for (i, (gene, _)) in ranked.iter().enumerate() {
        if pathway_set.contains(gene.to_uppercase().as_str()) {
            running += hit_weight;
        } else {
            running -= miss_weight;
        }
        if running.abs() > max_dev.abs() {
            max_dev = running;
            max_pos = i;
        }
    }

    // Collect leading edge: genes before (and including) the peak position that
    // are in the pathway set.
    let leading_edge: Vec<String> = ranked[..=max_pos]
        .iter()
        .filter(|(g, _)| pathway_set.contains(g.to_uppercase().as_str()))
        .map(|(g, _)| g.clone())
        .collect();

    (max_dev, leading_edge)
}

fn run_gsea_preranked(ranked: &[(String, f64)], min_size: usize, n_perm: usize) -> Vec<GseaResult> {
    use integration_layer::KEGG_PATHWAYS;

    let n_total = ranked.len();

    // Build a lookup of uppercase gene -> rank-metric for fast membership test.
    let ranked_upper: Vec<(String, f64)> =
        ranked.iter().map(|(g, m)| (g.to_uppercase(), *m)).collect();

    let mut results: Vec<GseaResult> = KEGG_PATHWAYS
        .iter()
        .filter_map(|pw| {
            let pathway_set: std::collections::HashSet<String> =
                pw.genes.iter().map(|g| g.to_uppercase()).collect();

            let n_hit = ranked_upper
                .iter()
                .filter(|(g, _)| pathway_set.contains(g.as_str()))
                .count();

            if n_hit < min_size {
                return None;
            }

            let (es, leading_edge) = enrichment_score(&ranked_upper, &pathway_set, n_total, n_hit);

            // Permutation test: shuffle gene labels, recompute ES.
            let mut null_es: Vec<f64> = Vec::with_capacity(n_perm);
            // Use a simple LCG RNG to avoid pulling in rand as a dependency.
            let mut lcg: u64 = 6_364_136_223_846_793_005;
            let mut shuffled: Vec<(String, f64)> = ranked_upper.clone();

            for _ in 0..n_perm {
                // Fisher-Yates shuffle with LCG.
                let len = shuffled.len();
                for i in (1..len).rev() {
                    lcg = lcg
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    let j = (lcg >> 33) as usize % (i + 1);
                    shuffled.swap(i, j);
                }
                let (null, _) = enrichment_score(&shuffled, &pathway_set, n_total, n_hit);
                null_es.push(null);
            }

            // NES: normalize by mean of null ES with the same sign.
            let same_sign: Vec<f64> = null_es
                .iter()
                .copied()
                .filter(|&v| if es >= 0.0 { v >= 0.0 } else { v < 0.0 })
                .collect();

            let mean_null = if same_sign.is_empty() {
                1.0
            } else {
                same_sign.iter().sum::<f64>() / same_sign.len() as f64
            };

            let nes = if mean_null.abs() < f64::EPSILON {
                es
            } else {
                es / mean_null.abs()
            };

            // p-value: fraction of permutation ES at least as extreme as observed ES.
            let extreme_count = null_es
                .iter()
                .filter(|&&v| if es >= 0.0 { v >= es } else { v <= es })
                .count();
            let p_value = (extreme_count as f64 + 1.0) / (n_perm as f64 + 1.0);

            Some(GseaResult {
                pathway_id: pw.id.to_string(),
                pathway_name: pw.name.to_string(),
                nes,
                p_value,
                pathway_size: pw.genes.len(),
                leading_edge,
            })
        })
        .collect();

    // Sort by NES descending (most enriched first).
    results.sort_unstable_by(|a, b| {
        b.nes
            .partial_cmp(&a.nes)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

// ── module definition ─────────────────────────────────────────────────────────

/// Python module `multiomics_py._core`.
///
/// Do not import this directly from Python — use `multiomics_py` (the wrapper
/// package) which parses JSON results into native dicts automatically.
#[pymodule]
fn multiomics_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Expose the crate version as `multiomics_py.__version__`.
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    m.add_function(wrap_pyfunction!(analyze_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_tsv, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_bed, m)?)?;
    m.add_function(wrap_pyfunction!(enrich_pathways, m)?)?;
    m.add_function(wrap_pyfunction!(gsea_preranked, m)?)?;
    m.add_function(wrap_pyfunction!(run_full_pipeline, m)?)?;
    m.add_function(wrap_pyfunction!(run_integration_from_json, m)?)?;

    Ok(())
}
