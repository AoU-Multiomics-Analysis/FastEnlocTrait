version 1.0

task CreateGWASShards {
    input {
        Int gwas_count
        Int gwas_units_per_shard = 10
    }

    command <<<
        set -euo pipefail
        Rscript -e '
        gwas_count <- as.integer("~{gwas_count}")
        shard_size <- as.integer("~{gwas_units_per_shard}")
        if (is.na(gwas_count) || gwas_count < 1) {
          stop("gwas_count must be at least 1")
        }
        if (is.na(shard_size) || shard_size < 1) {
          stop("gwas_units_per_shard must be at least 1")
        }
        dir.create("shards", showWarnings = FALSE)
        starts <- seq.int(0L, gwas_count - 1L, by = shard_size)
        for (i in seq_along(starts)) {
          idx <- starts[[i]]:min(gwas_count - 1L, starts[[i]] + shard_size - 1L)
          writeLines(as.character(idx), sprintf("shards/shard_%06d.indices.txt", i - 1L))
        }
        '
    >>>

    output {
        Array[File] shard_index_files = glob("shards/shard_*.indices.txt")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "1G"
        cpu: 1
    }
}

task RunColocShard {
    input {
        File shard_indices
        Array[File] GWASData
        Array[File] QTLData
        Array[String] study_ids
        Array[String] traits
        Array[String] n_variants
        Array[String] qtl_labels
        Float min_clpp = 0.01
        String clpp_output_prefix = "clpp"
    }

    command <<<
        set -euo pipefail

        mkdir -p work pair_outputs per_gwas_raw

        gwas_files="~{write_lines(GWASData)}"
        qtl_files="~{write_lines(QTLData)}"
        study_ids_file="~{write_lines(study_ids)}"
        traits_file="~{write_lines(traits)}"
        n_variants_file="~{write_lines(n_variants)}"
        qtl_labels_file="~{write_lines(qtl_labels)}"

        : > pair_gwas_indices.txt
        : > pair_qtl_indices.txt
        : > pair_study_ids.txt
        : > pair_qtl_labels.txt

        line_at() {
            sed -n "$(($1 + 1))p" "$2"
        }

        json_int_array() {
            awk 'BEGIN { printf "[" } { printf "%s%s", sep, $0; sep="," } END { print "]" }' "$1" > "$2"
        }

        normalize_fastenloc_output() {
            bash ~/normalize_fastenloc_output.sh "$1" "$2"
        }

        aggregate_with_qtl_label() {
            local files_file="$1"
            local labels_file="$2"
            local output_name="$3"
            local first_file
            local expected_header
            local expected_fields
            first_file=$(head -n 1 "$files_file")
            expected_header=$(head -n 1 "$first_file")
            expected_fields=$(awk -F'\t' 'NR == 1 { print NF }' "$first_file")
            {
                printf "qtl_label\t"
                printf '%s\n' "$expected_header"
            } > "$output_name"

            exec 3< "$files_file"
            exec 4< "$labels_file"
            while read -r f <&3 && read -r label <&4; do
                if [ "$(head -n 1 "$f")" != "$expected_header" ]; then
                    echo "Header mismatch while aggregating $f into $output_name" >&2
                    exit 1
                fi
                awk -F'\t' -v expected="$expected_fields" -v path="$f" '
                    NR > 1 && NF != expected {
                        printf "Malformed TSV row in %s at line %d: expected %d fields, found %d\n", path, NR, expected, NF > "/dev/stderr"
                        exit 1
                    }
                ' "$f"
                tail -n +2 "$f" | awk -v label="$label" 'BEGIN{OFS="\t"}{print label,$0}'
            done >> "$output_name"
            exec 3<&-
            exec 4<&-
        }

        qtl_count=$(wc -l < "$qtl_files" | tr -d ' ')
        if [ "$qtl_count" -lt 1 ]; then
            echo "At least one QTL input is required." >&2
            exit 1
        fi

        while read -r gwas_index; do
            [ -n "$gwas_index" ] || continue
            gwas_ord=$(printf "%06d" "$gwas_index")
            gwas_data=$(line_at "$gwas_index" "$gwas_files")
            study_id=$(line_at "$gwas_index" "$study_ids_file")
            trait=$(line_at "$gwas_index" "$traits_file")
            number_variants=$(line_at "$gwas_index" "$n_variants_file")

            gwas_raw="work/${gwas_ord}.${study_id}.gwas.raw.txt"
            gwas_prepped="work/${gwas_ord}.${study_id}.gwas.fastenloc.txt"
            if [[ "$gwas_data" == *.gz ]]; then
                gzip -dc "$gwas_data" > "$gwas_raw"
            else
                cp "$gwas_data" "$gwas_raw"
            fi
            awk -F'\t' '
            NR == 1 && ($6 == "annotation" || $6 == "locus_string") { next }
            { print }
            ' "$gwas_raw" > "$gwas_prepped"

            gene_files="work/${gwas_ord}.${study_id}.gene.files"
            enrich_files="work/${gwas_ord}.${study_id}.enrich.files"
            mi_files="work/${gwas_ord}.${study_id}.mi.files"
            sig_files="work/${gwas_ord}.${study_id}.sig.files"
            snp_files="work/${gwas_ord}.${study_id}.snp.files"
            clpp_files="work/${gwas_ord}.${study_id}.clpp.files"
            labels_for_gwas="work/${gwas_ord}.${study_id}.qtl_labels"
            : > "$gene_files"; : > "$enrich_files"; : > "$mi_files"; : > "$sig_files"; : > "$snp_files"; : > "$clpp_files"
            : > "$labels_for_gwas"

            qtl_index=0
            while [ "$qtl_index" -lt "$qtl_count" ]; do
                qtl_ord=$(printf "%06d" "$qtl_index")
                qtl_file=$(line_at "$qtl_index" "$qtl_files")
                qtl_label=$(line_at "$qtl_index" "$qtl_labels_file")
                pair_prefix="pair_outputs/${gwas_ord}.${qtl_ord}.${study_id}.${qtl_label}.combined"
                clpp_output="pair_outputs/${gwas_ord}.${qtl_ord}.${study_id}.${qtl_label}.~{clpp_output_prefix}.combined.tsv"

                fastenloc \
                  -eqtl "$qtl_file" \
                  -gwas "$gwas_prepped" \
                  -total_variants "$number_variants" \
                  -prefix "$pair_prefix"

                for out in "$pair_prefix".enloc.*.out; do
                    normalize_fastenloc_output "$out" "${out}.normalized"
                    header=$(head -n 1 "${out}.normalized")
                    {
                        printf "trait\t%s\n" "$header"
                        tail -n +2 "${out}.normalized" | awk -v trait="$trait" 'BEGIN{OFS="\t"}{print trait,$0}'
                    } > "${out}.tmp"
                    mv "${out}.tmp" "$out"
                    rm "${out}.normalized"
                done

                Rscript ~/clpp_fastenloc.R \
                  --gwas "$gwas_data" \
                  --qtl "$qtl_file" \
                  --out "$clpp_output" \
                  --min_clpp ~{min_clpp}

                echo "${pair_prefix}.enloc.gene.out" >> "$gene_files"
                echo "${pair_prefix}.enloc.enrich.out" >> "$enrich_files"
                echo "${pair_prefix}.enloc.mi.out" >> "$mi_files"
                echo "${pair_prefix}.enloc.sig.out" >> "$sig_files"
                echo "${pair_prefix}.enloc.snp.out" >> "$snp_files"
                echo "$clpp_output" >> "$clpp_files"
                echo "$qtl_label" >> "$labels_for_gwas"

                echo "$gwas_index" >> pair_gwas_indices.txt
                echo "$qtl_index" >> pair_qtl_indices.txt
                echo "$study_id" >> pair_study_ids.txt
                echo "$qtl_label" >> pair_qtl_labels.txt

                qtl_index=$((qtl_index + 1))
            done

            aggregate_with_qtl_label "$gene_files" "$labels_for_gwas" "per_gwas_raw/${gwas_ord}.${study_id}.combined.enloc.gene.out"
            aggregate_with_qtl_label "$enrich_files" "$labels_for_gwas" "per_gwas_raw/${gwas_ord}.${study_id}.combined.enloc.enrich.out"
            aggregate_with_qtl_label "$mi_files" "$labels_for_gwas" "per_gwas_raw/${gwas_ord}.${study_id}.combined.enloc.mi.out"
            aggregate_with_qtl_label "$sig_files" "$labels_for_gwas" "per_gwas_raw/${gwas_ord}.${study_id}.combined.enloc.sig.out"
            aggregate_with_qtl_label "$snp_files" "$labels_for_gwas" "per_gwas_raw/${gwas_ord}.${study_id}.combined.enloc.snp.out"
            aggregate_with_qtl_label "$clpp_files" "$labels_for_gwas" "per_gwas_raw/${gwas_ord}.${study_id}.~{clpp_output_prefix}.combined.tsv"
        done < "~{shard_indices}"

        json_int_array pair_gwas_indices.txt pair_gwas_indices.json
        json_int_array pair_qtl_indices.txt pair_qtl_indices.json
    >>>

    output {
        Array[Int] pair_gwas_indices = read_json("pair_gwas_indices.json")
        Array[Int] pair_qtl_indices = read_json("pair_qtl_indices.json")
        Array[String] pair_gwas_indices_text = read_lines("pair_gwas_indices.txt")
        Array[String] pair_study_ids = read_lines("pair_study_ids.txt")
        Array[String] pair_qtl_labels = read_lines("pair_qtl_labels.txt")
        Array[File] pair_gene_outputs = glob("pair_outputs/*.combined.enloc.gene.out")
        Array[File] pair_enrich_outputs = glob("pair_outputs/*.combined.enloc.enrich.out")
        Array[File] pair_mi_outputs = glob("pair_outputs/*.combined.enloc.mi.out")
        Array[File] pair_sig_outputs = glob("pair_outputs/*.combined.enloc.sig.out")
        Array[File] pair_snp_outputs = glob("pair_outputs/*.combined.enloc.snp.out")
        Array[File] pair_clpp_outputs = glob("pair_outputs/*.combined.tsv")
        Array[File] per_gwas_combined_gene_outputs = glob("per_gwas_raw/*.combined.enloc.gene.out")
        Array[File] per_gwas_combined_enrich_outputs = glob("per_gwas_raw/*.combined.enloc.enrich.out")
        Array[File] per_gwas_combined_mi_outputs = glob("per_gwas_raw/*.combined.enloc.mi.out")
        Array[File] per_gwas_combined_sig_outputs = glob("per_gwas_raw/*.combined.enloc.sig.out")
        Array[File] per_gwas_combined_snp_outputs = glob("per_gwas_raw/*.combined.enloc.snp.out")
        Array[File] per_gwas_combined_clpp_outputs = glob("per_gwas_raw/*.combined.tsv")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "64G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}

task AggregateHarmonizedByGWAS {
    input {
        Array[File] signal_files
        Array[File] credible_set_files
        Array[File] gene_files
        Array[String] pair_gwas_indices
        Array[String] study_ids
        String output_prefix = "harmonized_coloc"
    }

    command <<<
        set -euo pipefail
        mkdir -p per_gwas_harmonized lists

        signal_files="~{write_lines(signal_files)}"
        cs_files="~{write_lines(credible_set_files)}"
        gene_files="~{write_lines(gene_files)}"
        indices_file="~{write_lines(pair_gwas_indices)}"
        study_ids_file="~{write_lines(study_ids)}"

        pair_count=$(wc -l < "$indices_file" | tr -d ' ')
        for files_file in "$signal_files" "$cs_files" "$gene_files"; do
            observed=$(wc -l < "$files_file" | tr -d ' ')
            if [ "$observed" -ne "$pair_count" ]; then
                echo "Harmonized file count ${observed} does not match pair index count ${pair_count}." >&2
                exit 1
            fi
        done

        aggregate_gz_tsv() {
            local files_file="$1"
            local output_name="$2"
            local tmp_file
            local current_file
            local wrote_header=0
            tmp_file="${output_name%.gz}"
            current_file="${tmp_file}.current"
            : > "$tmp_file"
            while read -r f; do
                gzip -dc "$f" > "$current_file"
                if [ "$wrote_header" -eq 0 ]; then
                    head -n 1 "$current_file" > "$tmp_file"
                    wrote_header=1
                fi
                tail -n +2 "$current_file" >> "$tmp_file"
                rm -f "$current_file"
            done < "$files_file"
            if [ "$wrote_header" -eq 0 ]; then
                echo "No input files found for $output_name" >&2
                exit 1
            fi
            gzip -c "$tmp_file" > "$output_name"
            rm -f "$tmp_file" "$current_file"
        }

        gwas_count=$(wc -l < "$study_ids_file" | tr -d ' ')
        gwas_index=0
        while [ "$gwas_index" -lt "$gwas_count" ]; do
            gwas_ord=$(printf "%06d" "$gwas_index")
            study_id=$(sed -n "$((gwas_index + 1))p" "$study_ids_file")
            signal_list="lists/${gwas_ord}.signal.files"
            cs_list="lists/${gwas_ord}.cs.files"
            gene_list="lists/${gwas_ord}.gene.files"
            : > "$signal_list"; : > "$cs_list"; : > "$gene_list"

            paste "$indices_file" "$signal_files" "$cs_files" "$gene_files" | \
              awk -v idx="$gwas_index" -v s="$signal_list" -v c="$cs_list" -v g="$gene_list" \
                'BEGIN{FS=OFS="\t"} $1 == idx { print $2 >> s; print $3 >> c; print $4 >> g }'

            if [ ! -s "$signal_list" ]; then
                echo "No harmonized files found for GWAS index ${gwas_index}." >&2
                exit 1
            fi

            aggregate_gz_tsv "$signal_list" "per_gwas_harmonized/${gwas_ord}.${study_id}.~{output_prefix}.signal.tsv.gz"
            aggregate_gz_tsv "$cs_list" "per_gwas_harmonized/${gwas_ord}.${study_id}.~{output_prefix}.cs.tsv.gz"
            aggregate_gz_tsv "$gene_list" "per_gwas_harmonized/${gwas_ord}.${study_id}.~{output_prefix}.gene.tsv.gz"

            gwas_index=$((gwas_index + 1))
        done
    >>>

    output {
        Array[File] signal_outputs = glob("per_gwas_harmonized/*.signal.tsv.gz")
        Array[File] credible_set_outputs = glob("per_gwas_harmonized/*.cs.tsv.gz")
        Array[File] gene_outputs = glob("per_gwas_harmonized/*.gene.tsv.gz")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "32G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}
