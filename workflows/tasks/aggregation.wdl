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
