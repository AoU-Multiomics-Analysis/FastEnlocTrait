version 1.0

import "tasks/aggregation.wdl" as aggregation
import "tasks/clpp.wdl" as clpp
import "tasks/consensus.wdl" as consensus
import "tasks/fastenloc.wdl" as fastenloc
import "tasks/harmonize.wdl" as harmonize
import "tasks/input_validation.wdl" as input_validation
import "tasks/localize.wdl" as localize
import "tasks/split.wdl" as split

workflow RunFastenloc {
    input {
        File GWASManifest
        Array[File] QTLData
        Array[String] QTLLabels
        Float min_clpp = 0.01
        String clpp_output_prefix = "clpp"
        Float consensus_jaccard = 0.90
        String consensus_output_prefix = "consensus_loci"
        Float harmonized_fdr_level = 0.05
        String harmonized_output_prefix = "harmonized_coloc"
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

    scatter (gwas_index in range(length(ValidateGWASManifest.study_ids))) {
        String study_id = ValidateGWASManifest.study_ids[gwas_index]
        String trait = ValidateGWASManifest.traits[gwas_index]
        Int number_variants = ValidateGWASManifest.n_variants[gwas_index]
        String trait_category = ValidateGWASManifest.trait_categories[gwas_index]
        Int n_credible_sets = ValidateGWASManifest.n_credible_sets[gwas_index]
        File localized_gwas_data = LocalizeGWASData.gwas_data[gwas_index]

        call split.SplitFastenloc as SplitFastenloc {
          input:
            GWASData = localized_gwas_data
        }

        scatter (qtl_index in range(length(ValidateQTLInputs.labels))) {
            String qtl_label = ValidateQTLInputs.labels[qtl_index]
            File qtl_file = QTLData[qtl_index]

            scatter (chunk in SplitFastenloc.chunk_files) {
                call fastenloc.FastEnloc as FastEnloc {
                  input:
                    GWASData = chunk,
                    QTLData = qtl_file,
                    NumberVariants = number_variants
                }

                call clpp.CLPPFastEnloc as CLPPFastEnloc {
                  input:
                    GWASData = chunk,
                    QTLData = qtl_file,
                    min_clpp = min_clpp,
                    output_prefix = study_id + "." + qtl_label + "." + clpp_output_prefix
                }
            }

            call aggregation.AggregateFiles as AggregateGene {
              input:
                files = flatten(FastEnloc.gene_outputs),
                output_name = study_id + "." + qtl_label + ".combined.enloc.gene.out"
            }

            call aggregation.AggregateFiles as AggregateEnrich {
              input:
                files = flatten(FastEnloc.enrich_outputs),
                output_name = study_id + "." + qtl_label + ".combined.enloc.enrich.out"
            }

            call aggregation.AggregateFiles as AggregateMI {
              input:
                files = flatten(FastEnloc.mi_outputs),
                output_name = study_id + "." + qtl_label + ".combined.enloc.mi.out"
            }

            call aggregation.AggregateFiles as AggregateSig {
              input:
                files = flatten(FastEnloc.sig_outputs),
                output_name = study_id + "." + qtl_label + ".combined.enloc.sig.out"
            }

            call aggregation.AggregateFiles as AggregateSNP {
              input:
                files = flatten(FastEnloc.snp_outputs),
                output_name = study_id + "." + qtl_label + ".combined.enloc.snp.out"
            }

            call aggregation.AggregateFiles as AggregateCLPP {
              input:
                files = CLPPFastEnloc.clpp_output,
                output_name = study_id + "." + qtl_label + "." + clpp_output_prefix + ".combined.tsv"
            }

            call harmonize.HarmonizeColoc as HarmonizeColoc {
              input:
                sig_output = AggregateSig.combined,
                gene_enloc_output = AggregateGene.combined,
                clpp_output = AggregateCLPP.combined,
                gwas_data = localized_gwas_data,
                consensus_map = MergeCredibleSets.consensus_map,
                fdr_level = harmonized_fdr_level,
                output_prefix = study_id + "." + qtl_label + "." + harmonized_output_prefix,
                study = study_id,
                trait = trait,
                trait_category = trait_category,
                n_variants = number_variants,
                n_credible_sets = n_credible_sets,
                layer = qtl_label
            }
        }

        call aggregation.AggregateFilesWithQTLLabel as AggregateGWASGene {
          input:
            files = AggregateGene.combined,
            qtl_labels = ValidateQTLInputs.labels,
            output_name = study_id + ".combined.enloc.gene.out"
        }

        call aggregation.AggregateFilesWithQTLLabel as AggregateGWASEnrich {
          input:
            files = AggregateEnrich.combined,
            qtl_labels = ValidateQTLInputs.labels,
            output_name = study_id + ".combined.enloc.enrich.out"
        }

        call aggregation.AggregateFilesWithQTLLabel as AggregateGWASMI {
          input:
            files = AggregateMI.combined,
            qtl_labels = ValidateQTLInputs.labels,
            output_name = study_id + ".combined.enloc.mi.out"
        }

        call aggregation.AggregateFilesWithQTLLabel as AggregateGWASSig {
          input:
            files = AggregateSig.combined,
            qtl_labels = ValidateQTLInputs.labels,
            output_name = study_id + ".combined.enloc.sig.out"
        }

        call aggregation.AggregateFilesWithQTLLabel as AggregateGWASSNP {
          input:
            files = AggregateSNP.combined,
            qtl_labels = ValidateQTLInputs.labels,
            output_name = study_id + ".combined.enloc.snp.out"
        }

        call aggregation.AggregateFilesWithQTLLabel as AggregateGWASCLPP {
          input:
            files = AggregateCLPP.combined,
            qtl_labels = ValidateQTLInputs.labels,
            output_name = study_id + "." + clpp_output_prefix + ".combined.tsv"
        }

        call aggregation.AggregateGzTsvFiles as AggregateGWASHarmonizedSignal {
          input:
            files = HarmonizeColoc.signal_output,
            output_name = study_id + "." + harmonized_output_prefix + ".signal.tsv.gz"
        }

        call aggregation.AggregateGzTsvFiles as AggregateGWASHarmonizedCS {
          input:
            files = HarmonizeColoc.credible_set_output,
            output_name = study_id + "." + harmonized_output_prefix + ".cs.tsv.gz"
        }

        call aggregation.AggregateGzTsvFiles as AggregateGWASHarmonizedGene {
          input:
            files = HarmonizeColoc.gene_level_output,
            output_name = study_id + "." + harmonized_output_prefix + ".gene.tsv.gz"
        }
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllGene {
      input:
        files = AggregateGWASGene.combined,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.gene.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllEnrich {
      input:
        files = AggregateGWASEnrich.combined,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.enrich.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllMI {
      input:
        files = AggregateGWASMI.combined,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.mi.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllSig {
      input:
        files = AggregateGWASSig.combined,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.sig.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllSNP {
      input:
        files = AggregateGWASSNP.combined,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = "combined.enloc.snp.out"
    }

    call aggregation.AggregateFilesWithGWASMetadata as AggregateAllCLPP {
      input:
        files = AggregateGWASCLPP.combined,
        study_ids = ValidateGWASManifest.study_ids,
        traits = ValidateGWASManifest.traits,
        trait_categories = ValidateGWASManifest.trait_categories,
        n_variants = ValidateGWASManifest.n_variants_text,
        n_credible_sets = ValidateGWASManifest.n_credible_sets_text,
        output_name = clpp_output_prefix + ".combined.tsv"
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedSignal {
      input:
        files = AggregateGWASHarmonizedSignal.combined,
        output_name = harmonized_output_prefix + ".signal.tsv.gz"
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedCS {
      input:
        files = AggregateGWASHarmonizedCS.combined,
        output_name = harmonized_output_prefix + ".cs.tsv.gz"
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedGene {
      input:
        files = AggregateGWASHarmonizedGene.combined,
        output_name = harmonized_output_prefix + ".gene.tsv.gz"
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
      File harmonized_signal_out = AggregateAllHarmonizedSignal.combined
      File harmonized_credible_set_out = AggregateAllHarmonizedCS.combined
      File harmonized_gene_out = AggregateAllHarmonizedGene.combined
      Array[File] per_gwas_combined_gene_out = AggregateGWASGene.combined
      Array[File] per_gwas_combined_enrich_out = AggregateGWASEnrich.combined
      Array[File] per_gwas_combined_mi_out = AggregateGWASMI.combined
      Array[File] per_gwas_combined_sig_out = AggregateGWASSig.combined
      Array[File] per_gwas_combined_snp_out = AggregateGWASSNP.combined
      Array[File] per_gwas_combined_clpp_out = AggregateGWASCLPP.combined
      Array[File] per_gwas_harmonized_signal_out = AggregateGWASHarmonizedSignal.combined
      Array[File] per_gwas_harmonized_credible_set_out = AggregateGWASHarmonizedCS.combined
      Array[File] per_gwas_harmonized_gene_out = AggregateGWASHarmonizedGene.combined
      Array[Array[File]] per_gwas_qtl_combined_gene_out = AggregateGene.combined
      Array[Array[File]] per_gwas_qtl_combined_enrich_out = AggregateEnrich.combined
      Array[Array[File]] per_gwas_qtl_combined_mi_out = AggregateMI.combined
      Array[Array[File]] per_gwas_qtl_combined_sig_out = AggregateSig.combined
      Array[Array[File]] per_gwas_qtl_combined_snp_out = AggregateSNP.combined
      Array[Array[File]] per_gwas_qtl_combined_clpp_out = AggregateCLPP.combined
      Array[Array[File]] per_gwas_qtl_harmonized_signal_out = HarmonizeColoc.signal_output
      Array[Array[File]] per_gwas_qtl_harmonized_credible_set_out = HarmonizeColoc.credible_set_output
      Array[Array[File]] per_gwas_qtl_harmonized_gene_out = HarmonizeColoc.gene_level_output
    }
}
