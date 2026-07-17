#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

source("scripts/merge_credible_sets.R")

good <- tibble(
  group = "trait",
  study_id = "S1",
  gwas_cs = c("trait;S1_chr1.100.101_L1", "trait;S1_chr2.200.200_L2"),
  variant_id = c("chr1_100_A_G", "chr2_200_C_T")
)
good_manifest <- tibble(study_id = "S1", n_credible_sets = 2L)
stopifnot(inherits(validate_membership(good, good_manifest), "tbl_df"))

bad_id <- good
bad_id$gwas_cs[[1]] <- "trait;S1_chr1.0.0_LNone"
err <- tryCatch(
  {
    validate_membership(bad_id, good_manifest)
    NULL
  },
  error = identity
)
stopifnot(inherits(err, "error"), grepl("Placeholder", conditionMessage(err)))

bad_count_manifest <- tibble(study_id = "S1", n_credible_sets = 88L)
err <- tryCatch(
  {
    validate_membership(good, bad_count_manifest)
    NULL
  },
  error = identity
)
stopifnot(
  inherits(err, "error"),
  grepl("expected 88 but observed 2", conditionMessage(err))
)

cat("merge credible-set validation tests passed\n")
