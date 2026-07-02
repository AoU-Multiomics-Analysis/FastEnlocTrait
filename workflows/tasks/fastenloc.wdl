version 1.0

task FastEnloc {
    input {
        File GWASData
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
    " ~{GWASData}
    for f in "traits"/*.txt; do
        trait=$(basename "$f" .txt)

        fastenloc \
          -eqtl ~{QTLData} \
          -gwas "$f" \
          -total_variants  ~{NumberVariants} \
          -prefix "${trait}"

        # add a column for the trait thats analyzed to each enloc output
        for out in "${trait}".enloc.*.out; do
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
