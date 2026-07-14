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

task SummarizeRawColocByQTL {
    input {
        File gene_output
        File enrich_output
        File mi_output
        File sig_output
        File snp_output
        File clpp_output
        Array[String] qtl_labels
        String output_prefix = "raw_coloc_summary"
    }

    command <<<
        set -euo pipefail
        Rscript ~/summarize_raw_coloc.R \
          --gene "~{gene_output}" \
          --enrich "~{enrich_output}" \
          --mi "~{mi_output}" \
          --sig "~{sig_output}" \
          --snp "~{snp_output}" \
          --clpp "~{clpp_output}" \
          --qtl_labels "~{write_lines(qtl_labels)}" \
          --out "~{output_prefix}.all_qtl.tsv" \
          --per_qtl_dir raw_coloc_summary_by_qtl
    >>>

    output {
        File raw_coloc_summary = "~{output_prefix}.all_qtl.tsv"
        Array[File] per_qtl_raw_coloc_summary = glob("raw_coloc_summary_by_qtl/*.tsv")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "32G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}

task PlotColocSummary {
    input {
        File gene_output
        File coloc_rate_output
        String plot_data_output_name = "coloc_summary_plot_data.tsv"
        String png_output_name = "coloc_summary.png"
        String pdf_output_name = "coloc_summary.pdf"
    }

    command <<<
        set -euo pipefail
        Rscript ~/plot_coloc_summary.R \
          --gene "~{gene_output}" \
          --coloc_rate "~{coloc_rate_output}" \
          --plot_data_out "~{plot_data_output_name}" \
          --png_out "~{png_output_name}" \
          --pdf_out "~{pdf_output_name}" \
          --min_coloc_genes 1

        Rscript ~/plot_coloc_summary.R \
          --gene "~{gene_output}" \
          --coloc_rate "~{coloc_rate_output}" \
          --plot_data_out "coloc_summary_plot_data.min5.tsv" \
          --png_out "coloc_summary.min5.png" \
          --pdf_out "coloc_summary.min5.pdf" \
          --min_coloc_genes 5

        Rscript ~/plot_coloc_summary.R \
          --gene "~{gene_output}" \
          --coloc_rate "~{coloc_rate_output}" \
          --plot_data_out "coloc_summary_plot_data.min10.tsv" \
          --png_out "coloc_summary.min10.png" \
          --pdf_out "coloc_summary.min10.pdf" \
          --min_coloc_genes 10
    >>>

    output {
        File plot_data = "~{plot_data_output_name}"
        File png = "~{png_output_name}"
        File pdf = "~{pdf_output_name}"
        File plot_data_min5 = "coloc_summary_plot_data.min5.tsv"
        File png_min5 = "coloc_summary.min5.png"
        File pdf_min5 = "coloc_summary.min5.pdf"
        File plot_data_min10 = "coloc_summary_plot_data.min10.tsv"
        File png_min10 = "coloc_summary.min10.png"
        File pdf_min10 = "coloc_summary.min10.pdf"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "8G"
        disks: "local-disk 50 SSD"
        cpu: 1
    }
}

task CollectHighLevelOutputs {
    input {
        File normalized_gwas_manifest
        File localized_gwas_manifest
        File consensus_loci
        File consensus_loci_summary
        File combined_gene
        File combined_enrich
        File combined_mi
        File combined_sig
        File combined_snp
        File combined_clpp
        File raw_coloc_summary
        Array[File] per_qtl_raw_coloc_summaries
        File harmonized_signal
        File harmonized_credible_set
        File harmonized_gene
        File coloc_rate_by_trait
        File gene_summary_by_trait
        File colocalizing_genes_long
        File coloc_summary_plot_data
        File coloc_summary_png
        File coloc_summary_pdf
        File coloc_summary_plot_data_min5
        File coloc_summary_png_min5
        File coloc_summary_pdf_min5
        File coloc_summary_plot_data_min10
        File coloc_summary_png_min10
        File coloc_summary_pdf_min10
        String archive_name = "high_level_coloc_outputs.tar.gz"
    }

    command <<<
        set -euo pipefail

        outdir="high_level_coloc_outputs"
        mkdir -p \
          "$outdir/manifests" \
          "$outdir/consensus" \
          "$outdir/raw_fastenloc" \
          "$outdir/raw_summaries/per_qtl" \
          "$outdir/harmonized" \
          "$outdir/final_summaries" \
          "$outdir/figures"

        cp "~{normalized_gwas_manifest}" "$outdir/manifests/normalized_gwas_manifest.tsv"
        cp "~{localized_gwas_manifest}" "$outdir/manifests/localized_gwas_manifest.tsv"

        cp "~{consensus_loci}" "$outdir/consensus/consensus_loci.tsv.gz"
        cp "~{consensus_loci_summary}" "$outdir/consensus/consensus_loci.summary.tsv"

        cp "~{combined_gene}" "$outdir/raw_fastenloc/combined.enloc.gene.out"
        cp "~{combined_enrich}" "$outdir/raw_fastenloc/combined.enloc.enrich.out"
        cp "~{combined_mi}" "$outdir/raw_fastenloc/combined.enloc.mi.out"
        cp "~{combined_sig}" "$outdir/raw_fastenloc/combined.enloc.sig.out"
        cp "~{combined_snp}" "$outdir/raw_fastenloc/combined.enloc.snp.out"
        cp "~{combined_clpp}" "$outdir/raw_fastenloc/clpp.combined.tsv"

        cp "~{raw_coloc_summary}" "$outdir/raw_summaries/raw_coloc_summary.all_qtl.tsv"
        per_qtl_files="~{write_lines(per_qtl_raw_coloc_summaries)}"
        while read -r f; do
          [ -n "$f" ] || continue
          cp "$f" "$outdir/raw_summaries/per_qtl/$(basename "$f")"
        done < "$per_qtl_files"

        cp "~{harmonized_signal}" "$outdir/harmonized/harmonized_coloc.signal.tsv.gz"
        cp "~{harmonized_credible_set}" "$outdir/harmonized/harmonized_coloc.cs.tsv.gz"
        cp "~{harmonized_gene}" "$outdir/harmonized/harmonized_coloc.gene.tsv.gz"

        cp "~{coloc_rate_by_trait}" "$outdir/final_summaries/coloc_rate_by_trait.tsv"
        cp "~{gene_summary_by_trait}" "$outdir/final_summaries/gene_summary_by_trait.tsv"
        cp "~{colocalizing_genes_long}" "$outdir/final_summaries/colocalizing_genes_long.tsv"
        cp "~{coloc_summary_plot_data}" "$outdir/final_summaries/coloc_summary_plot_data.tsv"
        cp "~{coloc_summary_png}" "$outdir/figures/coloc_summary.png"
        cp "~{coloc_summary_pdf}" "$outdir/figures/coloc_summary.pdf"
        cp "~{coloc_summary_plot_data_min5}" "$outdir/final_summaries/coloc_summary_plot_data.min5.tsv"
        cp "~{coloc_summary_png_min5}" "$outdir/figures/coloc_summary.min5.png"
        cp "~{coloc_summary_pdf_min5}" "$outdir/figures/coloc_summary.min5.pdf"
        cp "~{coloc_summary_plot_data_min10}" "$outdir/final_summaries/coloc_summary_plot_data.min10.tsv"
        cp "~{coloc_summary_png_min10}" "$outdir/figures/coloc_summary.min10.png"
        cp "~{coloc_summary_pdf_min10}" "$outdir/figures/coloc_summary.min10.pdf"

        find "$outdir" -type f | sort > "$outdir/CONTENTS.txt"
        tar -czf "~{archive_name}" "$outdir"
    >>>

    output {
        File high_level_outputs_archive = "~{archive_name}"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "8G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}
