#!/usr/bin/env Rscript

source("scripts/summarize_raw_coloc.R")

input <- tibble::tribble(
  ~qtl_label, ~study_id, ~gwas_trait, ~Gene, ~GRCP, ~GLCP,
  "eQTL", "study1", "trait1", "gene1", 0.2, 0.1,
  "eQTL", "study1", "trait1", "gene2", 0.3, 0.2,
  "sQTL", "study1", "trait1", "gene3", 0.4, 0.3
)

eqtl <- summarize_one(input, "fastenloc_gene", "eQTL")
sqtl <- summarize_one(input, "fastenloc_gene", "sQTL")

stopifnot(eqtl$n_rows == 2, eqtl$n_genes == 2)
stopifnot(sqtl$n_rows == 1, sqtl$n_genes == 1)

mixed <- tempfile(fileext = ".out")
writeLines(c(
  paste("study_id", "gwas_trait", "trait_category", "n_variants",
        "n_credible_sets", "qtl_label", "trait", "a0 a1 p_eqtl p_gwas", sep = "\t"),
  paste("study1", "trait one", "category", "100", "2", "eQTL", "trait one",
        "-11.2  3.4  1.2e-2  1.7e-5", sep = "\t"),
  paste("study1", "trait one", "category", "100", "2", "sQTL", "trait one",
        "-10.2  2.4  2.2e-2  2.7e-5", sep = "\t")
), mixed)
parsed <- read_fastenloc_combined(mixed, "fastenloc_mi")
stopifnot(nrow(parsed) == 2, identical(parsed$qtl_label, c("eQTL", "sQTL")))
stopifnot(identical(parsed$a0, c("-11.2", "-10.2")))

cat("summarize_raw_coloc tests passed\n")
