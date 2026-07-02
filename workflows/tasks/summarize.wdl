version 1.0

task SummarizeColoc {
    input {
        File credible_set_output
        File gene_output
        File gtf
        String coloc_rate_output_name = "coloc_rate_by_trait.tsv"
        String gene_summary_output_name = "gene_summary_by_trait.tsv"
        String gene_list_output_name = "colocalizing_genes_long.tsv"
        String gene_threshold = "any"
    }

    command <<<
        set -euo pipefail
        Rscript ~/summarize_coloc.R \
          --cs "~{credible_set_output}" \
          --gene "~{gene_output}" \
          --gtf "~{gtf}" \
          --coloc_rate_out "~{coloc_rate_output_name}" \
          --gene_summary_out "~{gene_summary_output_name}" \
          --gene_list_out "~{gene_list_output_name}" \
          --gene_threshold "~{gene_threshold}"
    >>>

    output {
        File coloc_rate_output = "~{coloc_rate_output_name}"
        File gene_summary_output = "~{gene_summary_output_name}"
        File gene_list_output = "~{gene_list_output_name}"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "32G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}
