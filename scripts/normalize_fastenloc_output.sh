#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 INPUT OUTPUT" >&2
    exit 2
fi

input_name="$1"
output_name="$2"

# fastENLOC data rows are fixed-width whitespace while its headers can contain
# tabs (including empty spacer fields). The enrichment output is the exception:
# it has no header, and its row labels can contain spaces. Give it a stable
# header and preserve each label by reading the two metric fields from the
# right.
case "$input_name" in
    *.enloc.enrich.out)
        awk 'BEGIN{OFS="\t"; print "term", "estimate", "standard_error"}
            NF {
                if (NF < 3) {
                    printf "Malformed fastENLOC enrichment row at line %d: expected at least 3 whitespace fields\n", NR > "/dev/stderr"
                    exit 1
                }
                label=$1
                for (i=2; i<=NF-2; i++) label=label " " $i
                print label,$(NF-1),$NF
            }
        ' "$input_name" > "$output_name"
        ;;
    *)
        awk 'BEGIN{OFS="\t"} NF {$1=$1; print}' "$input_name" > "$output_name"
        ;;
esac

expected_fields=$(awk -F'\t' 'NR == 1 { print NF }' "$output_name")
if [ -z "$expected_fields" ]; then
    echo "Empty fastENLOC output: $input_name" >&2
    exit 1
fi

awk -F'\t' -v expected="$expected_fields" -v path="$input_name" '
    NF != expected {
        printf "Malformed fastENLOC row in %s at line %d: expected %d TSV fields, found %d\n", path, NR, expected, NF > "/dev/stderr"
        exit 1
    }
' "$output_name"
