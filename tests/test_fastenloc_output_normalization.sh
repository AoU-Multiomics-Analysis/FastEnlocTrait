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
        END { if (NR != expected_rows) exit 1 }
    ' expected_rows="$([ "$kind" = enrich ] && printf 3 || printf 2)" "$output"

    if [ "$kind" = enrich ]; then
        expected_output=$(printf 'term\testimate\tstandard_error\nIntercept\t-13.0\t-\nEnrichment (no shrinkage)\t6.5\t0.16')
        if [ "$(cat "$output")" != "$expected_output" ]; then
            echo "Unexpected normalized enrichment output" >&2
            diff -u <(printf '%s\n' "$expected_output") "$output" >&2 || true
            exit 1
        fi
    fi
}

check_fixture sig 7
check_fixture snp 6
check_fixture gene 3
check_fixture mi 4
check_fixture enrich 3

# Reproduce the RunColocShard failure mode: enrichment estimates differ across
# QTL layers, but the synthesized headers must remain identical so that the
# per-GWAS aggregation succeeds.
printf 'Intercept   -13.0   -\nEnrichment (no shrinkage)  6.5  0.16\n' \
    > "$tmp/eqtl.enloc.enrich.out"
printf 'Intercept   -11.5   -\nEnrichment (no shrinkage)  4.2  0.21\n' \
    > "$tmp/sqtl.enloc.enrich.out"

for label in eqtl sqtl; do
    raw="$tmp/${label}.enloc.enrich.out"
    normalized="$tmp/${label}.normalized"
    pair="$tmp/${label}.pair"
    bash "$repo_root/scripts/normalize_fastenloc_output.sh" "$raw" "$normalized"
    header=$(head -n 1 "$normalized")
    {
        printf 'trait\t%s\n' "$header"
        tail -n +2 "$normalized" | awk 'BEGIN{OFS="\t"}{print "Celiac disease",$0}'
    } > "$pair"
done

expected_header=$(head -n 1 "$tmp/eqtl.pair")
if [ "$(head -n 1 "$tmp/sqtl.pair")" != "$expected_header" ]; then
    echo "Enrichment headers still vary across QTL layers" >&2
    exit 1
fi

{
    printf 'qtl_label\t%s\n' "$expected_header"
    tail -n +2 "$tmp/eqtl.pair" | awk 'BEGIN{OFS="\t"}{print "eQTL",$0}'
    tail -n +2 "$tmp/sqtl.pair" | awk 'BEGIN{OFS="\t"}{print "sQTL",$0}'
} > "$tmp/per_gwas.combined.enloc.enrich.out"

awk -F'\t' '
    NR == 1 && $0 != "qtl_label\ttrait\tterm\testimate\tstandard_error" { exit 1 }
    NF != 5 { exit 1 }
    $1 == "eQTL" { eqtl_rows++ }
    $1 == "sQTL" { sqtl_rows++ }
    END { if (NR != 5 || eqtl_rows != 2 || sqtl_rows != 2) exit 1 }
' "$tmp/per_gwas.combined.enloc.enrich.out"

echo "fastENLOC output normalization tests passed"
