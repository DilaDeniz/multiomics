#!/usr/bin/env bash
# Regenerate the bcftools reference statistics used by
# genomics_core/tests/bcftools_validation.rs.
#
# Requires bcftools on PATH (see validation/environment.yml).
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out="$root/tests/golden/bcftools"
mkdir -p "$out"

for vcf in \
    "$root/test_data/validation/bcftools_edge_cases.vcf" \
    "$root/test_data/validation/1kg_phase3_chr22_16.05-16.6Mb.vcf.gz"; do
    name="$(basename "$vcf")"
    name="${name%.gz}"
    name="${name%.vcf}"
    {
        echo "# $(bcftools --version | head -n1)"
        echo "# bcftools stats test_data/validation/$(basename "$vcf") | grep -E '^(SN|TSTV)'"
        bcftools stats "$vcf" | grep -E '^(SN|TSTV)'
    } > "$out/$name.stats.txt"
    echo "wrote $out/$name.stats.txt"
done
