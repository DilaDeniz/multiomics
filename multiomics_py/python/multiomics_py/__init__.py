"""multiomics_py — Python bindings for Multiomics.

This thin wrapper imports the compiled Rust extension ``multiomics_py._core`` and
re-exports every public function with automatic JSON decoding so that callers
always receive plain Python ``dict`` / ``list`` objects.

Example usage
-------------
>>> import multiomics_py as bmo
>>> genomics = bmo.analyze_vcf("variants.vcf")
>>> transcr  = bmo.analyze_tsv("expr.tsv")
>>> epigen   = bmo.analyze_bed("meth.bed")
>>> result   = bmo.run_full_pipeline("variants.vcf", "expr.tsv", "meth.bed")
>>> enriched = bmo.enrich_pathways(["TP53", "KRAS", "BRCA1"], min_overlap=2)
"""

from __future__ import annotations

import json
from typing import Any

from multiomics_py import _core  # type: ignore[import]

__version__: str = _core.__version__
__all__ = [
    "__version__",
    "analyze_vcf",
    "analyze_tsv",
    "analyze_bed",
    "enrich_pathways",
    "gsea_preranked",
    "run_full_pipeline",
    "run_integration_from_json",
]


def analyze_vcf(path: str) -> dict[str, Any]:
    """Parse a VCF file and return a ``GenomicsSummary`` dict.

    Parameters
    ----------
    path:
        Filesystem path to the VCF file (plain or gzip-compressed).

    Returns
    -------
    dict
        Keys mirror the ``GenomicsSummary`` Rust struct:
        ``total_variants``, ``snp_count``, ``indel_count``, ``titv_ratio``,
        ``per_chrom``, ``high_impact``, ``af_histogram``, ``unique_positions``,
        ``high_impact_genes``.

    Raises
    ------
    RuntimeError
        If the file cannot be opened or parsed.
    """
    return json.loads(_core.analyze_vcf(path))


def analyze_tsv(path: str) -> dict[str, Any]:
    """Parse a gene-expression TSV file and return a ``TranscriptomicsSummary`` dict.

    The TSV must have a header row with sample names followed by rows of the
    form ``gene_id\\t<tpm1>\\t<tpm2>\\t...``.

    Parameters
    ----------
    path:
        Filesystem path to the expression matrix TSV.

    Returns
    -------
    dict
        Keys mirror the ``TranscriptomicsSummary`` Rust struct:
        ``total_genes``, ``expressed_genes``, ``low_expression_genes``,
        ``gene_stats``, ``top_100_expressed``, ``diff_expr``,
        ``sample_count``, ``sample_names``.

    Raises
    ------
    RuntimeError
        If the file cannot be opened or parsed.
    """
    return json.loads(_core.analyze_tsv(path))


def analyze_bed(path: str) -> dict[str, Any]:
    """Parse a BED methylation file and return an ``EpigenomicsSummary`` dict.

    Expected ENCODE bisulfite BED format:
    ``chrom start end name score strand``
    where *score* is the methylation percentage on the 0–1000 scale.

    Parameters
    ----------
    path:
        Filesystem path to the BED file.

    Returns
    -------
    dict
        Keys mirror the ``EpigenomicsSummary`` Rust struct:
        ``total_sites``, ``global_methylation_pct``, ``per_chrom``,
        ``cpg_islands``, ``hypermethylated``, ``hypomethylated``.

    Raises
    ------
    RuntimeError
        If the file cannot be opened or parsed.
    """
    return json.loads(_core.analyze_bed(path))


def enrich_pathways(
    genes: list[str],
    min_overlap: int = 2,
) -> list[dict[str, Any]]:
    """Run hypergeometric pathway enrichment against the built-in KEGG table.

    Parameters
    ----------
    genes:
        List of HGNC gene symbols (case-insensitive).
    min_overlap:
        Minimum number of query genes that must overlap a pathway for it to
        appear in the results. Default: 2.

    Returns
    -------
    list[dict]
        Each element contains:
        ``pathway_id``, ``pathway_name``, ``overlap``, ``pathway_size``,
        ``query_size``, ``score``, ``p_value``, ``padj``.
        Sorted by BH-adjusted p-value ascending (most significant first).

    Raises
    ------
    RuntimeError
        On serialisation failure (extremely unlikely).
    """
    return json.loads(_core.enrich_pathways(genes, min_overlap))


def gsea_preranked(
    ranked: list[tuple[str, float]],
    min_size: int = 5,
    n_perm: int = 1000,
) -> list[dict[str, Any]]:
    """Run a lightweight pre-ranked GSEA against the built-in KEGG pathway table.

    This is a KS-style enrichment score (Kolmogorov-Smirnov inspired) with a
    permutation-based p-value, suitable for exploratory analysis. For
    publication-quality results use a dedicated tool (GSEApy, fgsea).

    Parameters
    ----------
    ranked:
        List of ``(gene_symbol, metric)`` tuples sorted **descending** by
        metric (e.g. log₂FC × −log₁₀(p)).
    min_size:
        Minimum number of ranked genes that must be present in a pathway for
        it to be tested. Default: 5.
    n_perm:
        Number of permutations for p-value estimation. Default: 1000.

    Returns
    -------
    list[dict]
        Each element contains:
        ``pathway_id``, ``pathway_name``, ``nes``, ``p_value``,
        ``pathway_size``, ``leading_edge``.
        Sorted by NES descending.

    Raises
    ------
    RuntimeError
        If ``ranked`` is empty or serialisation fails.
    """
    return json.loads(_core.gsea_preranked(ranked, min_size, n_perm))


def run_full_pipeline(
    vcf: str,
    tsv: str,
    bed: str,
    no_ml: bool = False,
) -> dict[str, Any]:
    """Run the complete multi-omics integration pipeline.

    Parses all three input files concurrently, runs per-modality analyses, and
    then performs cross-modality integration: Pearson correlation, PCA, pathway
    enrichment, and biological insight derivation.

    Parameters
    ----------
    vcf:
        Path to the VCF variant file.
    tsv:
        Path to the gene-expression matrix TSV.
    bed:
        Path to the BED methylation file.
    no_ml:
        When ``True`` skip PCA and correlation (returns identity matrix).
        Default: ``False``.

    Returns
    -------
    dict
        Keys mirror the ``IntegrationSummary`` Rust struct:
        ``correlation_matrix``, ``pca``, ``top_pathways``, ``insights``.

    Raises
    ------
    RuntimeError
        If any input file cannot be parsed or integration analysis fails.
    """
    return json.loads(_core.run_full_pipeline(vcf, tsv, bed, no_ml))


def run_integration_from_json(
    genomics_json: str,
    transcr_json: str,
    epigen_json: str,
) -> dict[str, Any]:
    """Re-run integration using pre-computed JSON modality summaries.

    Useful when you already have summaries from previous runs and only want to
    re-run the integration layer without re-parsing raw files.

    Parameters
    ----------
    genomics_json:
        JSON string produced by :func:`analyze_vcf` (or its ``_core`` counterpart).
    transcr_json:
        JSON string produced by :func:`analyze_tsv`.
    epigen_json:
        JSON string produced by :func:`analyze_bed`.

    Returns
    -------
    dict
        Same structure as :func:`run_full_pipeline`.

    Raises
    ------
    RuntimeError
        If any JSON string is malformed or integration fails.
    """
    return json.loads(_core.run_integration_from_json(genomics_json, transcr_json, epigen_json))
