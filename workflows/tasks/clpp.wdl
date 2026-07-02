version 1.0

task CLPPFastEnloc {
    input {
        File TraitData
        File QTLData
        Float min_clpp = 0.01
        String output_prefix = "clpp"
    }

    command <<<
        set -euo pipefail
        output_name="~{output_prefix}.pairs.tsv"
        Rscript ~/clpp_fastenloc.R \
          --gwas ~{TraitData} \
          --qtl ~{QTLData} \
          --out "$output_name" \
          --min_clpp ~{min_clpp}
    >>>

    output {
        File clpp_output = "~{output_prefix}.pairs.tsv"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "32G"
        cpu: 1
    }
}
