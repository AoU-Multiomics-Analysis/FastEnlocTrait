version 1.0

task HarmonizeColoc {
    input {
        File sig_output
        File gene_enloc_output
        File clpp_output
        File gwas_data
        File consensus_map
        Float fdr_level = 0.05
        String output_prefix = "harmonized_coloc"
        String study
        String trait
        String trait_category
        Int n_variants
        Int n_credible_sets
        String layer = ""
    }

    command <<<
        set -euo pipefail
        Rscript ~/harmonize_coloc.R \
          --sig "~{sig_output}" \
          --gene "~{gene_enloc_output}" \
          --clpp "~{clpp_output}" \
          --gwas "~{gwas_data}" \
          --consensus_map "~{consensus_map}" \
          --out "~{output_prefix}.signal.tsv.gz" \
          --cs_out "~{output_prefix}.cs.tsv.gz" \
          --gene_out "~{output_prefix}.gene.tsv.gz" \
          --fdr_level ~{fdr_level} \
          --study "~{study}" \
          --trait "~{trait}" \
          --trait_category "~{trait_category}" \
          --n_variants ~{n_variants} \
          --n_credible_sets ~{n_credible_sets} \
          --layer "~{layer}"
    >>>

    output {
        File signal_output = "~{output_prefix}.signal.tsv.gz"
        File credible_set_output = "~{output_prefix}.cs.tsv.gz"
        File gene_level_output = "~{output_prefix}.gene.tsv.gz"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "64G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}
