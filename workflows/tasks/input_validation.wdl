version 1.0

task ValidateGWASManifest {
    input {
        File gwas_manifest
    }

    command <<<
        set -euo pipefail
        # shellcheck disable=SC2016
        Rscript -e '
        manifest <- read.delim(
          "~{gwas_manifest}",
          header = TRUE,
          sep = "\t",
          stringsAsFactors = FALSE,
          check.names = FALSE,
          comment.char = "",
          quote = ""
        )

        names(manifest)[names(manifest) == "#study_id"] <- "study_id"
        required <- c("study_id", "trait", "n_variants", "gwas_path", "trait_category", "n_credible_sets")
        missing <- setdiff(required, names(manifest))
        if (length(missing) > 0) {
          stop("GWASManifest is missing required column(s): ", paste(missing, collapse = ", "))
        }

        manifest <- manifest[, required]
        row_count <- nrow(manifest)
        if (row_count == 0) stop("GWASManifest must contain at least one data row")

        normalize_label <- function(x) {
          tolower(trimws(gsub("[[:space:]]+", " ", x)))
        }
        manifest$trait <- normalize_label(manifest$trait)
        manifest$trait_category <- normalize_label(manifest$trait_category)

        fail_if <- function(condition, message) {
          if (any(condition, na.rm = TRUE)) stop(message)
        }

        fail_if(is.na(manifest$study_id) | manifest$study_id == "", "study_id cannot be empty")
        fail_if(!grepl("^[A-Za-z0-9._-]+$", manifest$study_id), "study_id values must match [A-Za-z0-9._-]+")
        if (any(duplicated(manifest$study_id))) {
          dup <- unique(manifest$study_id[duplicated(manifest$study_id)])
          stop("study_id values must be unique. Duplicate value(s): ", paste(dup, collapse = ", "))
        }

        fail_if(is.na(manifest$trait) | manifest$trait == "", "trait cannot be empty")
        fail_if(is.na(manifest$gwas_path) | manifest$gwas_path == "", "gwas_path cannot be empty")
        fail_if(is.na(manifest$trait_category) | manifest$trait_category == "", "trait_category cannot be empty")
        fail_if(!grepl("^[0-9]+$", manifest$n_variants), "n_variants must be a positive integer")
        fail_if(as.integer(manifest$n_variants) < 1, "n_variants must be a positive integer")
        fail_if(!grepl("^[0-9]+$", manifest$n_credible_sets), "n_credible_sets must be a positive integer")
        fail_if(as.integer(manifest$n_credible_sets) < 1, "n_credible_sets must be a positive integer")

        manifest$n_variants <- as.integer(manifest$n_variants)
        manifest$n_credible_sets <- as.integer(manifest$n_credible_sets)

        write.table(manifest, "gwas_manifest.normalized.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
        writeLines(manifest$study_id, "study_ids.txt")
        writeLines(manifest$trait, "traits.txt")
        writeLines(as.character(manifest$n_variants), "n_variants.txt")
        writeLines(manifest$gwas_path, "gwas_paths.txt")
        writeLines(manifest$trait_category, "trait_categories.txt")
        writeLines(as.character(manifest$n_credible_sets), "n_credible_sets.txt")
        writeLines(paste0("[", paste(manifest$n_variants, collapse = ","), "]"), "n_variants.json")
        writeLines(paste0("[", paste(manifest$n_credible_sets, collapse = ","), "]"), "n_credible_sets.json")
        '
    >>>

    output {
        File normalized_manifest = "gwas_manifest.normalized.tsv"
        Array[String] study_ids = read_lines("study_ids.txt")
        Array[String] traits = read_lines("traits.txt")
        Array[String] n_variants_text = read_lines("n_variants.txt")
        Array[Int] n_variants = read_json("n_variants.json")
        Array[String] gwas_paths = read_lines("gwas_paths.txt")
        Array[String] trait_categories = read_lines("trait_categories.txt")
        Array[String] n_credible_sets_text = read_lines("n_credible_sets.txt")
        Array[Int] n_credible_sets = read_json("n_credible_sets.json")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "1G"
        cpu: 1
    }
}

task ValidateQTLInputs {
    input {
        Int qtl_data_count
        Array[String] qtl_labels
    }

    command <<<
        set -euo pipefail
        labels_file="~{write_lines(qtl_labels)}"
        label_count=$(wc -l < "$labels_file" | tr -d ' ')

        if [ "~{qtl_data_count}" -eq 0 ]; then
            echo "At least one QTL input is required." >&2
            exit 1
        fi

        if [ "~{qtl_data_count}" -ne "$label_count" ]; then
            echo "QTLData and QTLLabels must have the same length: got ~{qtl_data_count} QTL files and ${label_count} labels." >&2
            exit 1
        fi

        if awk 'length($0) == 0 { exit 1 }' "$labels_file"; then
            :
        else
            echo "QTLLabels cannot contain empty labels." >&2
            exit 1
        fi

        if awk '$0 !~ /^[A-Za-z0-9._-]+$/ { exit 1 }' "$labels_file"; then
            :
        else
            echo "QTLLabels must be filename-safe and match [A-Za-z0-9._-]+." >&2
            exit 1
        fi

        duplicates=$(sort "$labels_file" | uniq -d)
        if [ -n "$duplicates" ]; then
            echo "QTLLabels must be unique. Duplicate labels:" >&2
            echo "$duplicates" >&2
            exit 1
        fi

        cp "$labels_file" qtl_labels.txt
    >>>

    output {
        Array[String] labels = read_lines("qtl_labels.txt")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "1G"
        cpu: 1
    }
}
