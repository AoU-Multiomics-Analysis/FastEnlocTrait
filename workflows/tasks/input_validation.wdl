version 1.0

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
