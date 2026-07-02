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
        head -n 1 "$first_file" > ~{output_name}
        while read -r f; do
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

        {
            printf "qtl_label\t"
            head -n 1 "$first_file"
        } > ~{output_name}

        exec 3< "$files_file"
        exec 4< "$labels_file"
        while read -r f <&3 && read -r label <&4; do
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
            head -n 1 "$first_file"
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
        Rscript -e '
        files <- readLines("~{write_lines(files)}")
        out <- gzfile("~{output_name}", "wt")
        on.exit(close(out))
        wrote_header <- FALSE
        for (path in files) {
          con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
          header <- readLines(con, n = 1)
          if (length(header) == 0) {
            close(con)
            next
          }
          if (!wrote_header) {
            writeLines(header, out)
            wrote_header <- TRUE
          }
          repeat {
            chunk <- readLines(con, n = 100000)
            if (length(chunk) == 0) break
            writeLines(chunk, out)
          }
          close(con)
        }
        if (!wrote_header) stop("No input rows found while aggregating gzipped TSV files")
        '
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
