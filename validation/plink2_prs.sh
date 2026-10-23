#!/usr/bin/env bash
# Regenerate the PLINK 2 polygenic-score reference used by
# genomics_core/tests/prs_validation.rs.
#
# Inputs: test_data/validation/prs_1kg_{HG00096,NA12878}.vcf — genotypes of two
# 1000 Genomes samples at the loci of genomics_core/data/prs_catalog.tsv, from
# the GRCh38 20190312 biallelic SNV/indel release, fetched with
#   bcftools view -R <catalog sites> -s HG00096,NA12878 \
#     .../ALL.chrN.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
# Requires python3 and plink2 on PATH (see validation/environment.yml).
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out="$root/tests/golden/plink2"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$out"

# One weight column per disease: ln(OR) of the risk allele, 0 elsewhere.
python3 - "$root/genomics_core/data/prs_catalog.tsv" "$tmp/scores.tsv" <<'PY'
import csv, math, sys
rows = [r for r in csv.DictReader((l for l in open(sys.argv[1]) if not l.startswith("#")), delimiter="\t")]
diseases = sorted({r["disease"] for r in rows})
weights = {}
for r in rows:
    key = (f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}", r["risk_allele"])
    weights.setdefault(key, {})[r["disease"]] = math.log(float(r["odds_ratio"]))
with open(sys.argv[2], "w") as f:
    f.write("ID\tA1\t" + "\t".join(d.replace(" ", "_") for d in diseases) + "\n")
    for (vid, a1), w in sorted(weights.items()):
        f.write(f"{vid}\t{a1}\t" + "\t".join(repr(w.get(d, 0.0)) for d in diseases) + "\n")
PY

for sample in HG00096 NA12878; do
    plink2 --vcf "$root/test_data/validation/prs_1kg_$sample.vcf" \
        --set-all-var-ids '@:#:$r:$a' \
        --score "$tmp/scores.tsv" 1 2 header-read cols=scoresums no-mean-imputation \
        --score-col-nums 3-8 --out "$tmp/$sample" > /dev/null
    { echo "# $(plink2 --version)"; cat "$tmp/$sample.sscore"; } > "$out/prs_1kg_$sample.sscore"
    echo "wrote $out/prs_1kg_$sample.sscore"
done
