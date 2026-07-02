version 1.0

task MergeCredibleSets {
    input {
        Array[File] gwas_data
        Array[String] study_ids
        Array[String] traits
        Float jaccard = 0.90
        String output_prefix = "consensus_loci"
    }

    command <<<
        set -euo pipefail

        Rscript -e '
        gwas_data <- readLines("~{write_lines(gwas_data)}")
        study_ids <- readLines("~{write_lines(study_ids)}")
        traits <- readLines("~{write_lines(traits)}")

        n <- length(gwas_data)
        if (length(study_ids) != n || length(traits) != n) {
          stop("gwas_data, study_ids, and traits must have matching lengths")
        }

        manifest <- data.frame(
          study_id = study_ids,
          trait = traits,
          gwas_path = gwas_data,
          stringsAsFactors = FALSE
        )
        write.table(
          manifest,
          "localized_gwas_manifest.tsv",
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )
        '

        Rscript ~/merge_credible_sets.R \
          --manifest localized_gwas_manifest.tsv \
          --out "~{output_prefix}.tsv.gz" \
          --summary_out "~{output_prefix}.summary.tsv" \
          --jaccard ~{jaccard} \
          --group_col trait
    >>>

    output {
        File consensus_map = "~{output_prefix}.tsv.gz"
        File consensus_summary = "~{output_prefix}.summary.tsv"
        File localized_manifest = "localized_gwas_manifest.tsv"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/fastenloctrait:main"
        memory: "32G"
        disks: "local-disk 500 SSD"
        cpu: 1
    }
}
