version 1.0

import "tasks/aggregation.wdl" as aggregation
import "tasks/consensus.wdl" as consensus
import "tasks/harmonize.wdl" as harmonize
import "tasks/input_validation.wdl" as input_validation
import "tasks/localize.wdl" as localize
import "tasks/shard_coloc.wdl" as shard_coloc
import "tasks/summarize.wdl" as summarize

workflow RunFastenloc {
    input {
        File GWASManifest
        File GTF
        Array[File] QTLData
        Array[String] QTLLabels
        Float min_clpp = 0.01
        String clpp_output_prefix = "clpp"
        Float consensus_jaccard = 0.90
        String consensus_output_prefix = "consensus_loci"
        Float harmonized_fdr_level = 0.05
        String harmonized_output_prefix = "harmonized_coloc"
        String summary_gene_threshold = "any"
        String raw_coloc_summary_prefix = "raw_coloc_summary"
        String coloc_rate_output_name = "coloc_rate_by_trait.tsv"
        String gene_summary_output_name = "gene_summary_by_trait.tsv"
        String gene_list_output_name = "colocalizing_genes_long.tsv"
        String coloc_summary_plot_data_name = "coloc_summary_plot_data.tsv"
        String coloc_summary_png_name = "coloc_summary.png"
        String coloc_summary_pdf_name = "coloc_summary.pdf"
        Int coloc_summary_min_genes = 1
        String high_level_outputs_archive_name = "high_level_coloc_outputs.tar.gz"
        Int gwas_units_per_shard = 10
    }

    call input_validation.ValidateGWASManifest as ValidateGWASManifest {
      input:
        gwas_manifest = GWASManifest
    }

    call input_validation.ValidateQTLInputs as ValidateQTLInputs {
      input:
        qtl_data_count = length(QTLData),
        qtl_labels = QTLLabels
    }

    scatter (localize_index in range(length(ValidateGWASManifest.study_ids))) {
        call localize.LocalizeGWASData as LocalizeGWASData {
          input:
            gwas_path = ValidateGWASManifest.gwas_paths[localize_index],
            study_id = ValidateGWASManifest.study_ids[localize_index]
        }
    }

    call consensus.MergeCredibleSets as MergeCredibleSets {
      input:
        gwas_data = LocalizeGWASData.gwas_data,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        jaccard = consensus_jaccard,
        output_prefix = consensus_output_prefix
    }

    call shard_coloc.CreateGWASShards as CreateGWASShards {
      input:
        gwas_count = length(ValidateGWASManifest.study_ids),
        gwas_units_per_shard = gwas_units_per_shard
    }

    scatter (shard_index in range(length(CreateGWASShards.shard_index_files))) {
        call shard_coloc.RunColocShard as RunColocShard {
          input:
            shard_indices = CreateGWASShards.shard_index_files[shard_index],
            GWASData = LocalizeGWASData.gwas_data,
            QTLData = QTLData,
            study_ids = ValidateGWASManifest.study_ids,
            traits = ValidateGWASManifest.traits,
            n_variants = ValidateGWASManifest.n_variants_text,
            qtl_labels = ValidateQTLInputs.labels,
            min_clpp = min_clpp,
            clpp_output_prefix = clpp_output_prefix
        }
    }

    Array[Int] pair_gwas_indices = flatten(RunColocShard.pair_gwas_indices)
    Array[Int] pair_qtl_indices = flatten(RunColocShard.pair_qtl_indices)
    Array[String] pair_gwas_indices_text = flatten(RunColocShard.pair_gwas_indices_text)
    Array[String] pair_study_ids = flatten(RunColocShard.pair_study_ids)
    Array[String] pair_qtl_labels = flatten(RunColocShard.pair_qtl_labels)
    Array[File] per_pair_gene_outputs = flatten(RunColocShard.pair_gene_outputs)
    Array[File] per_pair_enrich_outputs = flatten(RunColocShard.pair_enrich_outputs)
    Array[File] per_pair_mi_outputs = flatten(RunColocShard.pair_mi_outputs)
    Array[File] per_pair_sig_outputs = flatten(RunColocShard.pair_sig_outputs)
    Array[File] per_pair_snp_outputs = flatten(RunColocShard.pair_snp_outputs)
    Array[File] per_pair_clpp_outputs = flatten(RunColocShard.pair_clpp_outputs)
    Array[File] per_gwas_combined_gene_files = flatten(RunColocShard.per_gwas_combined_gene_outputs)
    Array[File] per_gwas_combined_enrich_files = flatten(RunColocShard.per_gwas_combined_enrich_outputs)
    Array[File] per_gwas_combined_mi_files = flatten(RunColocShard.per_gwas_combined_mi_outputs)
    Array[File] per_gwas_combined_sig_files = flatten(RunColocShard.per_gwas_combined_sig_outputs)
    Array[File] per_gwas_combined_snp_files = flatten(RunColocShard.per_gwas_combined_snp_outputs)
    Array[File] per_gwas_combined_clpp_files = flatten(RunColocShard.per_gwas_combined_clpp_outputs)

    scatter (pair_index in range(length(pair_gwas_indices))) {
        Int pair_gwas_index = pair_gwas_indices[pair_index]
        Int pair_qtl_index = pair_qtl_indices[pair_index]

        call harmonize.HarmonizeColoc as HarmonizeColoc {
          input:
            sig_output = per_pair_sig_outputs[pair_index],
            gene_enloc_output = per_pair_gene_outputs[pair_index],
            clpp_output = per_pair_clpp_outputs[pair_index],
            gwas_data = LocalizeGWASData.gwas_data[pair_gwas_index],
            consensus_map = MergeCredibleSets.consensus_map,
            fdr_level = harmonized_fdr_level,
            output_prefix = ValidateGWASManifest.study_ids[pair_gwas_index] + "." + ValidateQTLInputs.labels[pair_qtl_index] + "." + harmonized_output_prefix,
            study = ValidateGWASManifest.study_ids[pair_gwas_index],
            trait = ValidateGWASManifest.traits[pair_gwas_index],
            trait_category = ValidateGWASManifest.trait_categories[pair_gwas_index],
            n_variants = ValidateGWASManifest.n_variants[pair_gwas_index],
            n_credible_sets = ValidateGWASManifest.n_credible_sets[pair_gwas_index],
            layer = ValidateQTLInputs.labels[pair_qtl_index]
        }
    }

    call shard_coloc.AggregateHarmonizedByGWAS as AggregateHarmonizedByGWAS {
      input:
        signal_files = HarmonizeColoc.signal_output,
        credible_set_files = HarmonizeColoc.credible_set_output,
        gene_files = HarmonizeColoc.gene_level_output,
        pair_gwas_indices = pair_gwas_indices_text,
        study_ids = ValidateGWASManifest.study_ids,
        output_prefix = harmonized_output_prefix
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllGene {
      input:
        files = per_gwas_combined_gene_files,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.gene.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllEnrich {
      input:
        files = per_gwas_combined_enrich_files,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.enrich.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllMI {
      input:
        files = per_gwas_combined_mi_files,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.mi.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllSig {
      input:
        files = per_gwas_combined_sig_files,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.sig.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllSNP {
      input:
        files = per_gwas_combined_snp_files,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.snp.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllCLPP {
      input:
        files = per_gwas_combined_clpp_files,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = clpp_output_prefix + ".combined.tsv"
    }

    call summarize.SummarizeRawColocByQTL as SummarizeRawColocByQTL {
      input:
        gene_output = AggregateAllGene.combined,
        enrich_output = AggregateAllEnrich.combined,
        mi_output = AggregateAllMI.combined,
        sig_output = AggregateAllSig.combined,
        snp_output = AggregateAllSNP.combined,
        clpp_output = AggregateAllCLPP.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_prefix = raw_coloc_summary_prefix
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedSignal {
      input:
        files = AggregateHarmonizedByGWAS.signal_outputs,
        output_name = harmonized_output_prefix + ".signal.tsv.gz"
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedCS {
      input:
        files = AggregateHarmonizedByGWAS.credible_set_outputs,
        output_name = harmonized_output_prefix + ".cs.tsv.gz"
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedGene {
      input:
        files = AggregateHarmonizedByGWAS.gene_outputs,
        output_name = harmonized_output_prefix + ".gene.tsv.gz"
    }

    call summarize.SummarizeColoc as SummarizeColoc {
      input:
        credible_set_output = AggregateAllHarmonizedCS.combined,
        gene_output = AggregateAllHarmonizedGene.combined,
        gtf = GTF,
        coloc_rate_output_name = coloc_rate_output_name,
        gene_summary_output_name = gene_summary_output_name,
        gene_list_output_name = gene_list_output_name,
        gene_threshold = summary_gene_threshold
    }

    call summarize.PlotColocSummary as PlotColocSummary {
      input:
        gene_output = AggregateAllHarmonizedGene.combined,
        coloc_rate_output = SummarizeColoc.coloc_rate_output,
        plot_data_output_name = coloc_summary_plot_data_name,
        png_output_name = coloc_summary_png_name,
        pdf_output_name = coloc_summary_pdf_name,
        min_coloc_genes = coloc_summary_min_genes
    }

    call summarize.CollectHighLevelOutputs as CollectHighLevelOutputs {
      input:
        normalized_gwas_manifest = ValidateGWASManifest.normalized_manifest,
        localized_gwas_manifest = MergeCredibleSets.localized_manifest,
        consensus_loci = MergeCredibleSets.consensus_map,
        consensus_loci_summary = MergeCredibleSets.consensus_summary,
        combined_gene = AggregateAllGene.combined,
        combined_enrich = AggregateAllEnrich.combined,
        combined_mi = AggregateAllMI.combined,
        combined_sig = AggregateAllSig.combined,
        combined_snp = AggregateAllSNP.combined,
        combined_clpp = AggregateAllCLPP.combined,
        raw_coloc_summary = SummarizeRawColocByQTL.raw_coloc_summary,
        per_qtl_raw_coloc_summaries = SummarizeRawColocByQTL.per_qtl_raw_coloc_summary,
        harmonized_signal = AggregateAllHarmonizedSignal.combined,
        harmonized_credible_set = AggregateAllHarmonizedCS.combined,
        harmonized_gene = AggregateAllHarmonizedGene.combined,
        coloc_rate_by_trait = SummarizeColoc.coloc_rate_output,
        gene_summary_by_trait = SummarizeColoc.gene_summary_output,
        colocalizing_genes_long = SummarizeColoc.gene_list_output,
        coloc_summary_plot_data = PlotColocSummary.plot_data,
        coloc_summary_png = PlotColocSummary.png,
        coloc_summary_pdf = PlotColocSummary.pdf,
        archive_name = high_level_outputs_archive_name
    }

    output {
      File normalized_gwas_manifest = ValidateGWASManifest.normalized_manifest
      File localized_gwas_manifest = MergeCredibleSets.localized_manifest
      File consensus_loci_out = MergeCredibleSets.consensus_map
      File consensus_loci_summary_out = MergeCredibleSets.consensus_summary
      File combined_gene_out = AggregateAllGene.combined
      File combined_enrich_out = AggregateAllEnrich.combined
      File combined_mi_out = AggregateAllMI.combined
      File combined_sig_out = AggregateAllSig.combined
      File combined_snp_out = AggregateAllSNP.combined
      File combined_clpp_out = AggregateAllCLPP.combined
      File raw_coloc_summary_out = SummarizeRawColocByQTL.raw_coloc_summary
      Array[File] per_qtl_raw_coloc_summary_out = SummarizeRawColocByQTL.per_qtl_raw_coloc_summary
      File harmonized_signal_out = AggregateAllHarmonizedSignal.combined
      File harmonized_credible_set_out = AggregateAllHarmonizedCS.combined
      File harmonized_gene_out = AggregateAllHarmonizedGene.combined
      File coloc_rate_by_trait_out = SummarizeColoc.coloc_rate_output
      File gene_summary_by_trait_out = SummarizeColoc.gene_summary_output
      File colocalizing_genes_long_out = SummarizeColoc.gene_list_output
      File coloc_summary_plot_data_out = PlotColocSummary.plot_data
      File coloc_summary_png_out = PlotColocSummary.png
      File coloc_summary_pdf_out = PlotColocSummary.pdf
      File high_level_outputs_archive = CollectHighLevelOutputs.high_level_outputs_archive
      Array[File] per_gwas_combined_gene_out = per_gwas_combined_gene_files
      Array[File] per_gwas_combined_enrich_out = per_gwas_combined_enrich_files
      Array[File] per_gwas_combined_mi_out = per_gwas_combined_mi_files
      Array[File] per_gwas_combined_sig_out = per_gwas_combined_sig_files
      Array[File] per_gwas_combined_snp_out = per_gwas_combined_snp_files
      Array[File] per_gwas_combined_clpp_out = per_gwas_combined_clpp_files
      Array[File] per_gwas_harmonized_signal_out = AggregateHarmonizedByGWAS.signal_outputs
      Array[File] per_gwas_harmonized_credible_set_out = AggregateHarmonizedByGWAS.credible_set_outputs
      Array[File] per_gwas_harmonized_gene_out = AggregateHarmonizedByGWAS.gene_outputs
      Array[Int] per_gwas_qtl_gwas_index = pair_gwas_indices
      Array[Int] per_gwas_qtl_qtl_index = pair_qtl_indices
      Array[String] per_gwas_qtl_study_id = pair_study_ids
      Array[String] per_gwas_qtl_qtl_label = pair_qtl_labels
      Array[File] per_gwas_qtl_combined_gene_out = per_pair_gene_outputs
      Array[File] per_gwas_qtl_combined_enrich_out = per_pair_enrich_outputs
      Array[File] per_gwas_qtl_combined_mi_out = per_pair_mi_outputs
      Array[File] per_gwas_qtl_combined_sig_out = per_pair_sig_outputs
      Array[File] per_gwas_qtl_combined_snp_out = per_pair_snp_outputs
      Array[File] per_gwas_qtl_combined_clpp_out = per_pair_clpp_outputs
      Array[File] per_gwas_qtl_harmonized_signal_out = HarmonizeColoc.signal_output
      Array[File] per_gwas_qtl_harmonized_credible_set_out = HarmonizeColoc.credible_set_output
      Array[File] per_gwas_qtl_harmonized_gene_out = HarmonizeColoc.gene_level_output
    }
}
