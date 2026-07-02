version 1.0

task SplitFastenloc {
    input {
        File GWASData
        Int traits_per_chunk = 25
    }

    command <<<
        set -euo pipefail
        Rscript ~/SplitTraitData.R \
          --input ~{GWASData} \
          --traits-per-chunk ~{traits_per_chunk}
    >>>

    output {
        File manifest = "chunk_manifest.txt"
        Array[File] chunk_files = glob("chunks/*.txt")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "32G"
        cpu: 1
    }
}
