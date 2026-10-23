#!/usr/bin/env python3
"""Rebuild genomics_core/data/prs_catalog.tsv from authoritative sources.

For each (rsID, disease) candidate below:
  * GRCh38 position, reference/alternate alleles and the 1000 Genomes phase 3
    European (EUR) allele frequencies come from the Ensembl REST API;
  * the risk allele, odds ratio and PubMed ID come from the NHGRI-EBI GWAS
    Catalog REST API v2: the association with an odds ratio, a resolvable
    effect allele, p <= 5e-8 and a trait matching the disease, with the
    smallest p-value. Odds ratios below 1 are inverted onto the other allele.
Candidates without such an association are dropped (and reported).

Usage (network access required): python validation/prs_catalog.py
"""
import json
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "genomics_core" / "data" / "prs_catalog.tsv"

# Candidate loci (the original hand-written table); positions are looked up.
CANDIDATES = [
    ("rs6983267", "Colorectal cancer"), ("rs16892766", "Colorectal cancer"),
    ("rs4939827", "Colorectal cancer"), ("rs10795668", "Colorectal cancer"),
    ("rs4779584", "Colorectal cancer"), ("rs3824999", "Colorectal cancer"),
    ("rs4444235", "Colorectal cancer"), ("rs9929218", "Colorectal cancer"),
    ("rs2981582", "Breast cancer"), ("rs3803662", "Breast cancer"),
    ("rs889312", "Breast cancer"), ("rs1045485", "Breast cancer"),
    ("rs13281615", "Breast cancer"), ("rs8051542", "Breast cancer"),
    ("rs8176232", "Breast cancer"), ("rs11571833", "Breast cancer"),
    ("rs8034191", "Lung cancer"), ("rs2736100", "Lung cancer"),
    ("rs3117582", "Lung cancer"), ("rs6489769", "Lung cancer"),
    ("rs1051730", "Lung cancer"),
    ("rs6983267", "Prostate cancer"), ("rs4430796", "Prostate cancer"),
    ("rs10993994", "Prostate cancer"), ("rs10486567", "Prostate cancer"),
    ("rs7931342", "Prostate cancer"), ("rs1465618", "Prostate cancer"),
    ("rs258322", "Melanoma"), ("rs4911414", "Melanoma"), ("rs3731239", "Melanoma"),
    ("rs10757257", "Melanoma"), ("rs401681", "Melanoma"),
    ("rs3814113", "Ovarian cancer"), ("rs2072590", "Ovarian cancer"),
    ("rs8176232", "Ovarian cancer"), ("rs11571833", "Ovarian cancer"),
    ("rs10088218", "Ovarian cancer"),
]

# EFO trait labels accepted for each disease (case-insensitive exact match).
TRAITS = {
    "Colorectal cancer": {"colorectal cancer", "colorectal carcinoma", "colorectal adenocarcinoma"},
    "Breast cancer": {"breast carcinoma", "breast cancer", "estrogen-receptor positive breast cancer",
                      "estrogen-receptor negative breast cancer"},
    "Lung cancer": {"lung carcinoma", "lung cancer", "lung adenocarcinoma", "non-small cell lung carcinoma",
                    "squamous cell lung carcinoma", "small cell lung carcinoma"},
    "Prostate cancer": {"prostate carcinoma", "prostate cancer", "prostate adenocarcinoma"},
    "Melanoma": {"cutaneous melanoma", "melanoma"},
    "Ovarian cancer": {"ovarian carcinoma", "ovarian cancer", "ovarian serous carcinoma",
                       "high grade ovarian serous adenocarcinoma", "epithelial ovarian cancer"},
}
P_THRESHOLD = 5e-8


def get_json(url, data=None):
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, data=json.dumps(data).encode() if data else None, headers=headers)
            return json.load(urllib.request.urlopen(req, timeout=60))
        except Exception as exc:  # network hiccups: retry with backoff
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
            print(f"retrying {url}: {exc}", file=sys.stderr)


def ensembl_variants(ids):
    return get_json("https://rest.ensembl.org/variation/homo_sapiens?pops=1", {"ids": sorted(set(ids))})


def gwas_associations(rsid):
    url = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations?" + urllib.parse.urlencode(
        {"rs_id": rsid, "size": 500})
    return get_json(url).get("_embedded", {}).get("associations", [])


def main():
    ens = ensembl_variants([r for r, _ in CANDIDATES])
    rows, dropped = [], []
    for rsid, disease in CANDIDATES:
        var = ens.get(rsid)
        if var is None:
            dropped.append((rsid, disease, "not in Ensembl"))
            continue
        maps = [m for m in var["mappings"] if m["assembly_name"] == "GRCh38" and "_" not in m["seq_region_name"]]
        if not maps:
            dropped.append((rsid, disease, "no primary GRCh38 mapping"))
            continue
        m = maps[0]
        alleles = m["allele_string"].split("/")
        ref = alleles[0]
        eur = {p["allele"]: p["frequency"] for p in var.get("populations", [])
               if p["population"] == "1000GENOMES:phase_3:EUR"}
        # Use the two alleles observed in EUR (dbSNP often lists rare extras).
        pair = [a for a in alleles if eur.get(a, 0.0) > 0.0]
        if len(pair) != 2 or ref not in pair:
            dropped.append((rsid, disease, f"not biallelic with REF in 1000G EUR ({m['allele_string']}, {eur})"))
            continue
        alt = pair[1] if pair[0] == ref else pair[0]

        best = None
        for a in gwas_associations(rsid):
            traits = {t["efo_trait"].lower() for t in a.get("efo_traits", [])}
            if not traits & TRAITS[disease]:
                continue
            p = a.get("p_value")
            match = re.match(r"\s*([0-9]*\.?[0-9]+)", str(a.get("or_value") or ""))
            if p is None or p > P_THRESHOLD or not match or float(match.group(1)) <= 0:
                continue
            odds = float(match.group(1))
            effect = [s["effect_allele"] for s in a.get("snp_allele", []) if s.get("rs_id") == rsid]
            if len(effect) != 1 or effect[0] not in pair:
                continue
            if best is None or p < best[0]:
                genes = ",".join(a.get("mapped_genes", [])[:3]) or "."
                best = (p, odds, effect[0], a.get("pubmed_id", ""), sorted(traits & TRAITS[disease])[0], genes)
        if best is None:
            dropped.append((rsid, disease, "no qualifying GWAS Catalog association with an OR"))
            continue
        p, odds, effect, pmid, trait, genes = best
        risk = effect
        if odds < 1.0:
            risk, odds = (alt if effect == ref else ref), 1.0 / odds
        rows.append({
            "rsid": rsid, "chrom": m["seq_region_name"], "pos": m["start"], "ref": ref,
            "alt": alt, "risk_allele": risk,
            "risk_is_ref": int(risk == ref), "odds_ratio": round(odds, 4),
            "eur_risk_freq": round(eur[risk], 4), "disease": disease, "genes": genes, "trait": trait,
            "p_value": p, "pmid": pmid,
        })
    header = ["rsid", "chrom", "pos", "ref", "alt", "risk_allele", "risk_is_ref", "odds_ratio",
              "eur_risk_freq", "disease", "genes", "trait", "p_value", "pmid"]
    with OUT.open("w") as fh:
        fh.write("# Cancer GWAS loci for genomics_core::prs. Generated by validation/prs_catalog.py\n")
        fh.write("# Positions/alleles/EUR frequencies: Ensembl REST (GRCh38, 1000 Genomes phase 3 EUR).\n")
        fh.write("# Risk allele, odds ratio, p-value, PMID: NHGRI-EBI GWAS Catalog REST API v2.\n")
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[h]) for h in header) + "\n")
    print(f"wrote {OUT}: {len(rows)} loci kept")
    for d in dropped:
        print("dropped:", *d)


if __name__ == "__main__":
    main()
