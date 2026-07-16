version 1.0

task AggregateFiles {

    input {
        Array[File] files
        String output_name
    }
    command <<<
        set -euo pipefail
        files_file="~{write_lines(files)}"
        first_file=$(head -n 1 "$files_file")
        expected_header=$(head -n 1 "$first_file")
        expected_fields=$(awk -F'\t' 'NR == 1 { print NF }' "$first_file")
        printf '%s\n' "$expected_header" > ~{output_name}
        while read -r f; do
          if [ "$(head -n 1 "$f")" != "$expected_header" ]; then
            echo "Header mismatch while aggregating $f into ~{output_name}" >&2
            exit 1
          fi
          awk -F'\t' -v expected="$expected_fields" -v path="$f" '
            NR > 1 && NF != expected {
              printf "Malformed TSV row in %s at line %d: expected %d fields, found %d\n", path, NR, expected, NF > "/dev/stderr"
              exit 1
            }
          ' "$f"
          tail -n +2 "$f" >> ~{output_name}
        done < "$files_file"
    >>>

    output {
        File combined = "~{output_name}"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "64G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}

task AggregateFilesWithQTLLabel {

    input {
        Array[File] files
        Array[String] qtl_labels
        String output_name
    }
    command <<<
        set -euo pipefail
        files_file="~{write_lines(files)}"
        labels_file="~{write_lines(qtl_labels)}"
        first_file=$(head -n 1 "$files_file")
        expected_header=$(head -n 1 "$first_file")
        expected_fields=$(awk -F'\t' 'NR == 1 { print NF }' "$first_file")

        {
            printf "qtl_label\t"
            printf '%s\n' "$expected_header"
        } > ~{output_name}

        exec 3< "$files_file"
        exec 4< "$labels_file"
        while read -r f <&3 && read -r label <&4; do
            if [ "$(head -n 1 "$f")" != "$expected_header" ]; then
                echo "Header mismatch while aggregating $f into ~{output_name}" >&2
                exit 1
            fi
            awk -F'\t' -v expected="$expected_fields" -v path="$f" '
                NR > 1 && NF != expected {
                    printf "Malformed TSV row in %s at line %d: expected %d fields, found %d\n", path, NR, expected, NF > "/dev/stderr"
                    exit 1
                }
            ' "$f"
            tail -n +2 "$f" | awk -v label="$label" 'BEGIN{OFS="\t"}{print label,$0}'
        done >> ~{output_name}
    >>>

    output {
        File combined = "~{output_name}"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "64G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}

task AggregateFilesWithGWASMetadata {

    input {
        Array[File] files
        Array[String] study_ids
        Array[String] traits
        Array[String] trait_categories
        Array[String] n_variants
        Array[String] n_credible_sets
        String output_name
    }
    command <<<
        set -euo pipefail
        files_file="~{write_lines(files)}"
        study_ids_file="~{write_lines(study_ids)}"
        traits_file="~{write_lines(traits)}"
        trait_categories_file="~{write_lines(trait_categories)}"
        n_variants_file="~{write_lines(n_variants)}"
        n_credible_sets_file="~{write_lines(n_credible_sets)}"
        first_file=$(head -n 1 "$files_file")
        expected_header=$(head -n 1 "$first_file")
        expected_fields=$(awk -F'\t' 'NR == 1 { print NF }' "$first_file")

        expected=$(wc -l < "$files_file" | tr -d ' ')
        for metadata_file in "$study_ids_file" "$traits_file" "$trait_categories_file" "$n_variants_file" "$n_credible_sets_file"; do
            observed=$(wc -l < "$metadata_file" | tr -d ' ')
            if [ "$observed" -ne "$expected" ]; then
                echo "Metadata array length ${observed} does not match files length ${expected}." >&2
                exit 1
            fi
        done

        {
            printf "study_id\tgwas_trait\ttrait_category\tn_variants\tn_credible_sets\t"
            printf '%s\n' "$expected_header"
        } > ~{output_name}

        exec 3< "$files_file"
        exec 4< "$study_ids_file"
        exec 5< "$traits_file"
        exec 6< "$trait_categories_file"
        exec 7< "$n_variants_file"
        exec 8< "$n_credible_sets_file"
        while read -r f <&3 && \
              read -r study_id <&4 && \
              read -r trait <&5 && \
              read -r trait_category <&6 && \
              read -r n_variant <&7 && \
              read -r n_credible_set <&8; do
            if [ "$(head -n 1 "$f")" != "$expected_header" ]; then
                echo "Header mismatch while aggregating $f into ~{output_name}" >&2
                exit 1
            fi
            awk -F'\t' -v expected="$expected_fields" -v path="$f" '
                NR > 1 && NF != expected {
                    printf "Malformed TSV row in %s at line %d: expected %d fields, found %d\n", path, NR, expected, NF > "/dev/stderr"
                    exit 1
                }
            ' "$f"
            tail -n +2 "$f" | awk \
                -v study_id="$study_id" \
                -v trait="$trait" \
                -v trait_category="$trait_category" \
                -v n_variant="$n_variant" \
                -v n_credible_set="$n_credible_set" \
                'BEGIN{OFS="\t"}{print study_id,trait,trait_category,n_variant,n_credible_set,$0}'
        done >> ~{output_name}
    >>>

    output {
        File combined = "~{output_name}"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "64G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}

task AggregateGzTsvFiles {

    input {
        Array[File] files
        String output_name
    }
    command <<<
        set -euo pipefail
        Rscript ~/aggregate_gz_tsv.R \
          --files "~{write_lines(files)}" \
          --out "~{output_name}"
        gzip -t "~{output_name}"
    >>>

    output {
        File combined = "~{output_name}"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "64G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}
