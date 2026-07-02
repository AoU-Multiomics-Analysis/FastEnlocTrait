version 1.0

task LocalizeGWASData {
    input {
        String gwas_path
        String study_id
    }

    command <<<
        set -euo pipefail
        out="~{study_id}.gwas.fastenloc.tsv.gz"

        if [[ "~{gwas_path}" == gs://* ]]; then
            if [[ "~{gwas_path}" == *.gz ]]; then
                gsutil cp "~{gwas_path}" "$out"
            else
                gsutil cat "~{gwas_path}" | gzip -c > "$out"
            fi
        elif [[ "~{gwas_path}" == http://* || "~{gwas_path}" == https://* ]]; then
            if [[ "~{gwas_path}" == *.gz ]]; then
                curl -L "~{gwas_path}" -o "$out"
            else
                curl -L "~{gwas_path}" | gzip -c > "$out"
            fi
        else
            if [ ! -f "~{gwas_path}" ]; then
                echo "GWAS path for ~{study_id} is not accessible inside the task: ~{gwas_path}" >&2
                exit 1
            fi

            if [[ "~{gwas_path}" == *.gz ]]; then
                cp "~{gwas_path}" "$out"
            else
                gzip -c "~{gwas_path}" > "$out"
            fi
        fi

        if [ ! -s "$out" ]; then
            echo "Localized GWAS file is empty for ~{study_id}: ~{gwas_path}" >&2
            exit 1
        fi
    >>>

    output {
        File gwas_data = "~{study_id}.gwas.fastenloc.tsv.gz"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "8G"
        disks: "local-disk 100 SSD"
        cpu: 1
    }
}
