#!/usr/bin/env Rscript

source("scripts/summarize_raw_coloc.R")

tmp <- tempfile("fastenloc-parser-")
dir.create(tmp)

metadata_header <- paste(
  c("study_id", "gwas_trait", "trait_category", "n_variants",
    "n_credible_sets", "qtl_label", "trait"),
  collapse = "\t"
)
metadata <- paste(c("study1", "Crohn's disease", "immune", "1000", "10",
                    "edQTL", "Crohn's disease"), collapse = "\t")
signal <- "edQTL1:edQTL1_L1(@)Crohns_disease;study1_chr1.10.20_L1=8.000e-01[9.5e-01:2]"

legacy_sig <- file.path(tmp, "legacy.sig.out")
writeLines(c(
  paste(metadata_header,
        paste(c("Signal", "Num_SNP", "CPIP_qtl", "CPIP_gwas_marginal",
                "CPIP_gwas_qtl_prior", "RCP", "LCP"), collapse = "\t"),
        sep = "\t"),
  paste(metadata,
        paste(signal, "2", "9.0e-01", "8.0e-01", "8.5e-01", "7.5e-01"),
        "7.6e-01", sep = "\t")
), legacy_sig)

normalized_sig <- file.path(tmp, "normalized.sig.out")
writeLines(c(
  paste(metadata_header,
        paste(c("Signal", "Num_SNP", "CPIP_qtl", "CPIP_gwas_marginal",
                "CPIP_gwas_qtl_prior", "RCP", "LCP"), collapse = "\t"),
        sep = "\t"),
  paste(metadata, signal, "2", "9.0e-01", "8.0e-01", "8.5e-01",
        "7.5e-01", "7.6e-01", sep = "\t")
), normalized_sig)

legacy_summary <- summarize_table(legacy_sig, "fastenloc_sig", "edQTL")
normalized_summary <- summarize_table(normalized_sig, "fastenloc_sig", "edQTL")

stopifnot(
  identical(legacy_summary$n_rows, 1L),
  identical(legacy_summary$n_gwas_credible_sets, 1L),
  identical(legacy_summary$n_qtl_signals, 1L),
  identical(legacy_summary$n_RCP_ge_0.5, 1L),
  isTRUE(all.equal(legacy_summary$max_RCP, 0.75)),
  isTRUE(all.equal(legacy_summary$max_LCP, 0.76)),
  isTRUE(all.equal(
    legacy_summary[c("n_rows", "n_gwas_credible_sets", "n_qtl_signals",
                     "max_RCP", "max_LCP", "n_RCP_ge_0.5")],
    normalized_summary[c("n_rows", "n_gwas_credible_sets", "n_qtl_signals",
                         "max_RCP", "max_LCP", "n_RCP_ge_0.5")]
  ))
)

legacy_snp <- file.path(tmp, "legacy.snp.out")
writeLines(c(
  paste(metadata_header,
        paste(c("Signal", "SNP", "PIP_qtl", "PIP_gwas_marginal",
                "PIP_gwas_qtl_prior", "SCP"), collapse = "\t"),
        sep = "\t"),
  paste(metadata,
        paste(signal, "chr1_10_A_G", "9.0e-01", "8.0e-01", "8.5e-01", "7.0e-01"),
        sep = "\t")
), legacy_snp)
snp_summary <- summarize_table(legacy_snp, "fastenloc_snp", "edQTL")
stopifnot(
  identical(snp_summary$n_variants_or_snps, 1L),
  identical(snp_summary$n_gwas_credible_sets, 1L),
  identical(snp_summary$n_qtl_signals, 1L)
)

legacy_gene <- file.path(tmp, "legacy.gene.out")
writeLines(c(
  paste(metadata_header, "Gene", "", "GRCP", "GLCP", sep = "\t"),
  paste(metadata, "edQTL1", "", "8.0e-01", "9.0e-01", sep = "\t")
), legacy_gene)
gene <- read_any(legacy_gene) %>% repair_legacy_fastenloc("fastenloc_gene")
stopifnot(
  !any(str_detect(names(gene), "^V[0-9]+$")),
  identical(gene$Gene, "edQTL1"),
  isTRUE(all.equal(gene$GRCP, 0.8)),
  isTRUE(all.equal(gene$GLCP, 0.9))
)

cat("fastENLOC parser regression tests passed\n")
