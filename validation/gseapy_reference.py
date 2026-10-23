#!/usr/bin/env python3
"""Regenerate the gseapy over-representation reference used by
integration_layer/tests/enrichment_validation.rs.

The background is the union of all genes in the GMT, matching
integration_layer::pathway::over_representation.

Usage (from the repository root): python validation/gseapy_reference.py
"""
import json
import pathlib

import gseapy

root = pathlib.Path(__file__).resolve().parent.parent
gmt = root / "test_data" / "validation" / "enrichment_sets.gmt"
query = [g.strip() for g in (root / "test_data" / "validation" / "enrichment_query.txt").read_text().split() if g.strip()]

gene_sets = {}
for line in gmt.read_text().splitlines():
    fields = line.split("\t")
    gene_sets[fields[0]] = [g for g in fields[2:] if g]
background = set(g for genes in gene_sets.values() for g in genes)

enr = gseapy.enrich(gene_list=query, gene_sets=gene_sets, background=background, outdir=None, no_plot=True)
res = enr.results
terms = []
for _, row in res.iterrows():
    overlap, size = row["Overlap"].split("/")
    terms.append({
        "term": row["Term"],
        "overlap": int(overlap),
        "set_size": int(size),
        "p_value": float(row["P-value"]),
        "padj": float(row["Adjusted P-value"]),
    })
out = root / "tests" / "golden" / "gseapy" / "enrichment.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps({
    "gseapy_version": gseapy.__version__,
    "background_size": len(background),
    "query_in_background": len(set(query) & background),
    "terms": terms,
}, indent=1) + "\n")
print(f"wrote {out} ({len(terms)} terms)")
