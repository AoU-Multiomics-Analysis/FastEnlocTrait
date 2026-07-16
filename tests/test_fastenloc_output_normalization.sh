#!/usr/bin/env bash

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

check_fixture() {
    local kind="$1"
    local input="$tmp/fixture.enloc.${kind}.out"
    local output="$tmp/${kind}.output"
    local expected="$2"

    case "$kind" in
        sig)
            printf 'Signal\tNum_SNP\tCPIP_qtl\tCPIP_gwas_marginal\tCPIP_gwas_qtl_prior\tRCP\tLCP\n' > "$input"
            printf 'edQTL1:edQTL1_L1(@)Trait;study_chr1.10.20_L1=8e-1[9e-1:2]  2  9e-1  8e-1  8.5e-1  7.5e-1\t7.6e-1\n' >> "$input"
            ;;
        snp)
            printf 'Signal\tSNP\tPIP_qtl\tPIP_gwas_marginal\tPIP_gwas_qtl_prior\tSCP\n' > "$input"
            printf 'edQTL1:edQTL1_L1(@)Trait;study_chr1.10.20_L1=8e-1[9e-1:2]  chr1_10_A_G  9e-1  8e-1  8.5e-1  7e-1\n' >> "$input"
            ;;
        gene)
            printf 'Gene\t\tGRCP\tGLCP\n' > "$input"
            printf 'edQTL1\t\t8e-1\t9e-1\n' >> "$input"
            ;;
        mi)
            printf 'a0\ta1\tp_eqtl\tp_gwas\n' > "$input"
            printf '%s\n' '-12.0  6.0     1e-3  2e-4' >> "$input"
            ;;
        enrich)
            printf 'Intercept   -13.0   -\n' > "$input"
            printf 'Enrichment (no shrinkage)  6.5  0.16\n' >> "$input"
            ;;
    esac

    bash "$repo_root/scripts/normalize_fastenloc_output.sh" "$input" "$output"
    awk -F'\t' -v expected="$expected" '
        NF != expected { exit 1 }
        END { if (NR != 2) exit 1 }
    ' "$output"
}

check_fixture sig 7
check_fixture snp 6
check_fixture gene 3
check_fixture mi 4
check_fixture enrich 3

echo "fastENLOC output normalization tests passed"
