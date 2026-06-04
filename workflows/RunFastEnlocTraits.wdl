version 1.0


task SplitFastenloc {
    input {
        File FastEnlocTraitData
        Int traits_per_chunk = 25
    }

    command <<<
        set -euo pipefail
        Rscript ~/SplitTraitData.R \
          --input ~{FastEnlocTraitData} \
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

task FastEnloc {
    input {
        File TraitData 
        File QTLData 
        Int NumberVariants
    }

    command <<<
    mkdir -p traits results

    awk -F'\t' "
    NR == 1 && (\$6 == \"annotation\" || \$6 == \"locus_string\") { next }
    {
        split(\$6, a, \";\")
        trait = a[1]

        outfile = \"traits/\" trait \".txt\"

        print >> outfile

        close(outfile)
    }
    " ~{TraitData}      
    for f in traits/*.txt; do
        trait=$(basename "$f" .txt)

        fastenloc \
          -eqtl ~{QTLData} \
          -gwas "$f" \
          -total_variants  ~{NumberVariants} \
          -prefix "${trait}"
        
        # add a column for the trait thats analyzed to each enloc output
        for out in ${trait}.enloc.*.out; do
            header=$(head -n1 "$out")
            {
                echo -e "trait\t${header}"
                tail -n +2 "$out" | \
                    awk -v trait="$trait" 'BEGIN{OFS="\t"}{print trait,$0}'
            } > tmp
            mv tmp "$out"
        done
      done
    >>>    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "8G"
        cpu: 1
    }

    output {
        Array[File] gene_outputs = glob("*.enloc.gene.out")
        Array[File] enrich_outputs = glob("*.enloc.enrich.out")
        Array[File] mi_outputs = glob("*.enloc.mi.out")
        Array[File] sig_outputs = glob("*.enloc.sig.out")
        Array[File] snp_outputs = glob("*.enloc.snp.out")
        Array[File] all_outputs = glob("*.enloc.*.out")
    }
}

task AggregateFiles {

    input {
        Array[File] files
        String output_name
    }
    command <<<
        set -euo pipefail
        first_file=$(head -n 1 ~{write_lines(files)})
        head -n 1 "$first_file" > ~{output_name}
        while read f; do
          tail -n +2 "$f" >> ~{output_name}
        done < ~{write_lines(files)}
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
    call AggregateFiles as AggregateGene {
      input:
        files = flatten(FastEnloc.gene_outputs),
        output_name = "combined.enloc.gene.out"
    }

    call AggregateFiles as AggregateEnrich {
      input:
        files = flatten(FastEnloc.enrich_outputs),
        output_name = "combined.enloc.enrich.out"
    }

    call AggregateFiles as AggregateMI {
      input:
        files = flatten(FastEnloc.mi_outputs),
        output_name = "combined.enloc.mi.out"
    }

    call AggregateFiles as AggregateSig {
      input:
        files = flatten(FastEnloc.sig_outputs),
        output_name = "combined.enloc.sig.out"
    }

    call AggregateFiles as AggregateSNP {
      input:
        files = flatten(FastEnloc.snp_outputs),
        output_name = "combined.enloc.snp.out"
    }

    output {
      File combined_gene_out = AggregateGene.combined
      File combined_enrich_out = AggregateEnrich.combined
      File combined_mi_out = AggregateMI.combined
      File combined_sig_out = AggregateSig.combined
      File combined_snp_out = AggregateSNP.combined
    }
}
