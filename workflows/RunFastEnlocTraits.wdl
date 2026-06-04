version 1.0


task SplitFastenloc {
  input {
    File FastEnlocTraitData
    Int traits_per_chunk = 25
  }

  command <<<

    set -euo pipefail

    Rscript  /tmp/SplitTraitData.R \
      --input ~{FastEnlocTraitData} \
      --traits-per-chunk ~{traits_per_chunk}

  >>>

  output {
    File manifest = "chunk_manifest.txt"
    Array[File] chunk_files = read_lines("chunk_manifest.txt")
  }  
  
  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
    memory: "32G"
    cpu: 1
  }
}

task FastEnloc {
    input {
        File TraitData 
        File QTLData 
        Int NumberVariants
    }

    command <<<
    mkdir -p traits

    awk -F'\t' '
    {
        split($6, a, ";")
        trait = a[1]
        outfile = "traits/" trait ".txt"
        print >> outfile
        close(outfile)
    }

    ' ~{TraitData}'
    for f in traits/*.txt; do

        trait=$(basename "$f" .txt)

        /home/jupyter/fastenloc-3.2/src/fastenloc \
            -eqtl ~{QTLData} \
            -gwas "$f" \
            -total_variants ~{NumberVariants} \
            -prefix "$trait"

    done

    >>>
  output {
    Array[File] gene_outputs = glob("results/*.enloc.gene.out")
    Array[File] enrich_outputs = glob("results/*.enloc.enrich.out")
    Array[File] mi_outputs = glob("results/*.enloc.mi.out")
    Array[File] sig_outputs = glob("results/*.enloc.sig.out")
    Array[File] snp_outputs = glob("results/*.enloc.snp.out")
    Array[File] all_outputs = glob("results/*.enloc.*.out")
  }
}


workflow RunFastenloc {
  input {
    File FastEnlocTraitData
    File QTLData
    Int  NumberVariants
  }

  call SplitFastenloc {
    input:
      FastEnlocTraitData = FastEnlocTraitData
  }

  scatter (chunk in SplitFastenloc.chunk_files) {

    call FastEnloc {
      input:
        TraitData = chunk,
        QTLData = QTLData,
        NumberVariants = NumberVariants
    }

  }

  output {
    Array[File] trait_chunks = SplitFastenloc.chunk_files
    Array[Array[File]] gene_outputs = FastEnloc.gene_outputs
    Array[Array[File]] enrich_outputs = FastEnloc.enrich_outputs
    Array[Array[File]] mi_outputs = FastEnloc.mi_outputs
    Array[Array[File]] sig_outputs = FastEnloc.sig_outputs
    Array[Array[File]] snp_outputs = FastEnloc.snp_outputs
    Array[Array[File]] all_outputs = FastEnloc.all_outputs
  }

}
