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

task ValidateQTLInputs {
    input {
        Int qtl_data_count
        Array[String] qtl_labels
    }

    command <<<
        set -euo pipefail
        labels_file="~{write_lines(qtl_labels)}"
        label_count=$(wc -l < "$labels_file" | tr -d ' ')

        if [ "~{qtl_data_count}" -eq 0 ]; then
            echo "At least one QTL input is required." >&2
            exit 1
        fi

        if [ "~{qtl_data_count}" -ne "$label_count" ]; then
            echo "QTLData and QTLLabels must have the same length: got ~{qtl_data_count} QTL files and ${label_count} labels." >&2
            exit 1
        fi

        if awk 'length($0) == 0 { exit 1 }' "$labels_file"; then
            :
        else
            echo "QTLLabels cannot contain empty labels." >&2
            exit 1
        fi

        if awk '$0 !~ /^[A-Za-z0-9._-]+$/ { exit 1 }' "$labels_file"; then
            :
        else
            echo "QTLLabels must be filename-safe and match [A-Za-z0-9._-]+." >&2
            exit 1
        fi

        duplicates=$(sort "$labels_file" | uniq -d)
        if [ -n "$duplicates" ]; then
            echo "QTLLabels must be unique. Duplicate labels:" >&2
            echo "$duplicates" >&2
            exit 1
        fi

        cp "$labels_file" qtl_labels.txt
    >>>

    output {
        Array[String] labels = read_lines("qtl_labels.txt")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "1G"
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

task HarmonizeColoc {
    input {
        File sig_output
        File gene_enloc_output
        File clpp_output
        File gwas_data
        Float fdr_level = 0.05
        String output_prefix = "harmonized_coloc"
        String layer = ""
    }

    command <<<
        set -euo pipefail
        Rscript ~/harmonize_coloc.R \
          --sig "~{sig_output}" \
          --gene "~{gene_enloc_output}" \
          --clpp "~{clpp_output}" \
          --gwas "~{gwas_data}" \
          --out "~{output_prefix}.signal.tsv.gz" \
          --cs_out "~{output_prefix}.cs.tsv.gz" \
          --gene_out "~{output_prefix}.gene.tsv.gz" \
          --fdr_level ~{fdr_level} \
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


workflow RunFastenloc {
    input {
        File FastEnlocTraitData
        Array[File] QTLData
        Array[String] QTLLabels
        Int  NumberVariants
        Float min_clpp = 0.01
        String clpp_output_prefix = "clpp"
        Float harmonized_fdr_level = 0.05
        String harmonized_output_prefix = "harmonized_coloc"
    }

    call SplitFastenloc {
        input:
            FastEnlocTraitData = FastEnlocTraitData
    }

    call ValidateQTLInputs {
      input:
        qtl_data_count = length(QTLData),
        qtl_labels = QTLLabels
    }

    scatter (qtl_index in range(length(ValidateQTLInputs.labels))) {
        String qtl_label = ValidateQTLInputs.labels[qtl_index]
        File qtl_file = QTLData[qtl_index]

        scatter (chunk in SplitFastenloc.chunk_files) {
            call FastEnloc {
              input:
                TraitData = chunk,
                QTLData = qtl_file,
                NumberVariants = NumberVariants
            }

            call CLPPFastEnloc {
              input:
                TraitData = chunk,
                QTLData = qtl_file,
                min_clpp = min_clpp,
                output_prefix = qtl_label + "." + clpp_output_prefix
            }
        }

        call AggregateFiles as AggregateGene {
          input:
            files = flatten(FastEnloc.gene_outputs),
            output_name = qtl_label + ".combined.enloc.gene.out"
        }

        call AggregateFiles as AggregateEnrich {
          input:
            files = flatten(FastEnloc.enrich_outputs),
            output_name = qtl_label + ".combined.enloc.enrich.out"
        }

        call AggregateFiles as AggregateMI {
          input:
            files = flatten(FastEnloc.mi_outputs),
            output_name = qtl_label + ".combined.enloc.mi.out"
        }

        call AggregateFiles as AggregateSig {
          input:
            files = flatten(FastEnloc.sig_outputs),
            output_name = qtl_label + ".combined.enloc.sig.out"
        }

        call AggregateFiles as AggregateSNP {
          input:
            files = flatten(FastEnloc.snp_outputs),
            output_name = qtl_label + ".combined.enloc.snp.out"
        }

        call AggregateFiles as AggregateCLPP {
          input:
            files = CLPPFastEnloc.clpp_output,
            output_name = qtl_label + "." + clpp_output_prefix + ".combined.tsv"
        }

        call HarmonizeColoc {
          input:
            sig_output = AggregateSig.combined,
            gene_enloc_output = AggregateGene.combined,
            clpp_output = AggregateCLPP.combined,
            gwas_data = FastEnlocTraitData,
            fdr_level = harmonized_fdr_level,
            output_prefix = qtl_label + "." + harmonized_output_prefix,
            layer = qtl_label
        }
    }

    call AggregateFilesWithQTLLabel as AggregateAllGene {
      input:
        files = AggregateGene.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.gene.out"
    }

    call AggregateFilesWithQTLLabel as AggregateAllEnrich {
      input:
        files = AggregateEnrich.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.enrich.out"
    }

    call AggregateFilesWithQTLLabel as AggregateAllMI {
      input:
        files = AggregateMI.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.mi.out"
    }

    call AggregateFilesWithQTLLabel as AggregateAllSig {
      input:
        files = AggregateSig.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.sig.out"
    }

    call AggregateFilesWithQTLLabel as AggregateAllSNP {
      input:
        files = AggregateSNP.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.snp.out"
    }

    call AggregateFilesWithQTLLabel as AggregateAllCLPP {
      input:
        files = AggregateCLPP.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = clpp_output_prefix + ".combined.tsv"
    }

    call AggregateGzTsvFiles as AggregateAllHarmonizedSignal {
      input:
        files = HarmonizeColoc.signal_output,
        output_name = harmonized_output_prefix + ".signal.tsv.gz"
    }

    call AggregateGzTsvFiles as AggregateAllHarmonizedCS {
      input:
        files = HarmonizeColoc.credible_set_output,
        output_name = harmonized_output_prefix + ".cs.tsv.gz"
    }

    call AggregateGzTsvFiles as AggregateAllHarmonizedGene {
      input:
        files = HarmonizeColoc.gene_level_output,
        output_name = harmonized_output_prefix + ".gene.tsv.gz"
    }

    output {
      File combined_gene_out = AggregateAllGene.combined
      File combined_enrich_out = AggregateAllEnrich.combined
      File combined_mi_out = AggregateAllMI.combined
      File combined_sig_out = AggregateAllSig.combined
      File combined_snp_out = AggregateAllSNP.combined
      File combined_clpp_out = AggregateAllCLPP.combined
      File harmonized_signal_out = AggregateAllHarmonizedSignal.combined
      File harmonized_credible_set_out = AggregateAllHarmonizedCS.combined
      File harmonized_gene_out = AggregateAllHarmonizedGene.combined
      Array[File] per_qtl_combined_gene_out = AggregateGene.combined
      Array[File] per_qtl_combined_enrich_out = AggregateEnrich.combined
      Array[File] per_qtl_combined_mi_out = AggregateMI.combined
      Array[File] per_qtl_combined_sig_out = AggregateSig.combined
      Array[File] per_qtl_combined_snp_out = AggregateSNP.combined
      Array[File] per_qtl_combined_clpp_out = AggregateCLPP.combined
      Array[File] per_qtl_harmonized_signal_out = HarmonizeColoc.signal_output
      Array[File] per_qtl_harmonized_credible_set_out = HarmonizeColoc.credible_set_output
      Array[File] per_qtl_harmonized_gene_out = HarmonizeColoc.gene_level_output
    }
}
