version 1.0

task FastEnloc {
    input {
        File GWASData
        File QTLData
        Int NumberVariants
        String trait
        String output_prefix
    }

    command <<<
    set -euo pipefail

    if [[ "~{GWASData}" == *.gz ]]; then
        gzip -dc "~{GWASData}" > gwas.raw.txt
    else
        cp "~{GWASData}" gwas.raw.txt
    fi

    awk -F'\t' '
    NR == 1 && ($6 == "annotation" || $6 == "locus_string") { next }
    { print }
    ' gwas.raw.txt > gwas.fastenloc.txt

    fastenloc \
      -eqtl ~{QTLData} \
      -gwas gwas.fastenloc.txt \
      -total_variants ~{NumberVariants} \
      -prefix "~{output_prefix}"

    for out in "~{output_prefix}".enloc.*.out; do
        header=$(head -n1 "$out")
        {
            printf "trait\t%s\n" "$header"
            tail -n +2 "$out" | \
                awk -v trait="~{trait}" 'BEGIN{OFS="\t"}{print trait,$0}'
        } > tmp
        mv tmp "$out"
    done
    >>>
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "8G"
        cpu: 1
    }

    output {
        File gene_output = "~{output_prefix}.enloc.gene.out"
        File enrich_output = "~{output_prefix}.enloc.enrich.out"
        File mi_output = "~{output_prefix}.enloc.mi.out"
        File sig_output = "~{output_prefix}.enloc.sig.out"
        File snp_output = "~{output_prefix}.enloc.snp.out"
        Array[File] all_outputs = glob("*.enloc.*.out")
    }
}
