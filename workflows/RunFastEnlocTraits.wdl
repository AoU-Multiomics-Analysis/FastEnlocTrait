version 1.0

import "tasks/aggregation.wdl" as aggregation
import "tasks/clpp.wdl" as clpp
import "tasks/fastenloc.wdl" as fastenloc
import "tasks/harmonize.wdl" as harmonize
import "tasks/input_validation.wdl" as input_validation
import "tasks/split.wdl" as split

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

    call split.SplitFastenloc as SplitFastenloc {
        input:
            FastEnlocTraitData = FastEnlocTraitData
    }

    call input_validation.ValidateQTLInputs as ValidateQTLInputs {
      input:
        qtl_data_count = length(QTLData),
        qtl_labels = QTLLabels
    }

    scatter (qtl_index in range(length(ValidateQTLInputs.labels))) {
        String qtl_label = ValidateQTLInputs.labels[qtl_index]
        File qtl_file = QTLData[qtl_index]

        scatter (chunk in SplitFastenloc.chunk_files) {
            call fastenloc.FastEnloc as FastEnloc {
              input:
                TraitData = chunk,
                QTLData = qtl_file,
                NumberVariants = NumberVariants
            }

            call clpp.CLPPFastEnloc as CLPPFastEnloc {
              input:
                TraitData = chunk,
                QTLData = qtl_file,
                min_clpp = min_clpp,
                output_prefix = qtl_label + "." + clpp_output_prefix
            }
        }

        call aggregation.AggregateFiles as AggregateGene {
          input:
            files = flatten(FastEnloc.gene_outputs),
            output_name = qtl_label + ".combined.enloc.gene.out"
        }

        call aggregation.AggregateFiles as AggregateEnrich {
          input:
            files = flatten(FastEnloc.enrich_outputs),
            output_name = qtl_label + ".combined.enloc.enrich.out"
        }

        call aggregation.AggregateFiles as AggregateMI {
          input:
            files = flatten(FastEnloc.mi_outputs),
            output_name = qtl_label + ".combined.enloc.mi.out"
        }

        call aggregation.AggregateFiles as AggregateSig {
          input:
            files = flatten(FastEnloc.sig_outputs),
            output_name = qtl_label + ".combined.enloc.sig.out"
        }

        call aggregation.AggregateFiles as AggregateSNP {
          input:
            files = flatten(FastEnloc.snp_outputs),
            output_name = qtl_label + ".combined.enloc.snp.out"
        }

        call aggregation.AggregateFiles as AggregateCLPP {
          input:
            files = CLPPFastEnloc.clpp_output,
            output_name = qtl_label + "." + clpp_output_prefix + ".combined.tsv"
        }

        call harmonize.HarmonizeColoc as HarmonizeColoc {
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

    call aggregation.AggregateFilesWithQTLLabel as AggregateAllGene {
      input:
        files = AggregateGene.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.gene.out"
    }

    call aggregation.AggregateFilesWithQTLLabel as AggregateAllEnrich {
      input:
        files = AggregateEnrich.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.enrich.out"
    }

    call aggregation.AggregateFilesWithQTLLabel as AggregateAllMI {
      input:
        files = AggregateMI.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.mi.out"
    }

    call aggregation.AggregateFilesWithQTLLabel as AggregateAllSig {
      input:
        files = AggregateSig.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.sig.out"
    }

    call aggregation.AggregateFilesWithQTLLabel as AggregateAllSNP {
      input:
        files = AggregateSNP.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = "combined.enloc.snp.out"
    }

    call aggregation.AggregateFilesWithQTLLabel as AggregateAllCLPP {
      input:
        files = AggregateCLPP.combined,
        qtl_labels = ValidateQTLInputs.labels,
        output_name = clpp_output_prefix + ".combined.tsv"
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedSignal {
      input:
        files = HarmonizeColoc.signal_output,
        output_name = harmonized_output_prefix + ".signal.tsv.gz"
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedCS {
      input:
        files = HarmonizeColoc.credible_set_output,
        output_name = harmonized_output_prefix + ".cs.tsv.gz"
    }

    call aggregation.AggregateGzTsvFiles as AggregateAllHarmonizedGene {
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
