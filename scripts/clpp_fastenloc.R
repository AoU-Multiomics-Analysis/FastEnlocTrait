#!/usr/bin/env Rscript
# ============================================================================
# clpp_fastenloc.R
#
# Compute CLPP (eCAVIAR colocalization posterior; Hormozdiari et al. 2016)
# directly from FastENLOC-format input files.
#
#   CLPP(G, Q) = sum_{v in G ∩ Q}  PIP_gwas(v) * PIP_qtl(v)
#
# where G is a GWAS credible set, Q is a QTL signal, and the sum runs over
# variants shared (by exact variant id) between the two.
#
# FastENLOC column 6 formats:
#   GWAS : "Trait;study_cs_id=<pip>[<cpip>:<n>]"          (one per row)
#   QTL  : "gene:signal_id@=<pip>[...]|gene:signal@=<pip>" ('|'-delimited,
#          one or more signals; note the '@=' vs the GWAS bare '=')
#
# Usage:
#   Rscript clpp_fastenloc.R --gwas GWAS.vcf.gz --qtl QTL.vcf.gz \
#           --out pairs.tsv.gz [--min_clpp 0.01]
#
# Output is written as a gzip-compressed, tab-separated file.
#
# NOTE ON DROPPED VARIANTS / PIP NORMALIZATION -------------------------------
# This computes the STANDARD eCAVIAR CLPP: a variant present in one dataset's
# credible set but absent from the other's panel contributes 0 to the sum
# (absence = no colocalization evidence). PIPs are NOT renormalized over the
# shared-variant subset. This is intentional and conservative -- it can only
# UNDER-estimate colocalization, never inflate it. Renormalizing to conditional
# PIPs over shared variants is NOT recommended when panel overlap is low (as it
# is here: many GWAS credible sets share only their lead variant with the QTL
# panel), because it would reconstruct near-1.0 weights from little mass and
# manufacture false positives. The correct fix for panel mismatch is upstream
# variant harmonization (fine-map both on a common variant set), not rescaling.
# The `n_shared` column is emitted so low-overlap pairs can be flagged.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)   # fread() for fast file loading
  library(optparse)     # command-line option parsing
})

# ---- command-line options --------------------------------------------------
build_option_parser <- function() {
  OptionParser(
    usage = "Rscript %prog --gwas GWAS.vcf.gz --qtl QTL.vcf.gz --out pairs.tsv.gz [options]",
    description = paste(
      "Compute CLPP (eCAVIAR colocalization posterior) from FastENLOC-format",
      "GWAS and QTL files. Output is a gzip-compressed TSV of colocalizing",
      "credible-set x QTL-signal pairs."
    ),
    option_list = list(
      make_option("--gwas", type = "character", default = NULL,
                  help = "FastENLOC-format GWAS file (.vcf.gz or plain) [required]"),
      make_option("--qtl", type = "character", default = NULL,
                  help = "FastENLOC-format QTL file (.vcf.gz or plain) [required]"),
      make_option("--out", type = "character", default = NULL,
                  help = "Output path; gzipped if it ends in .gz [required]"),
      make_option("--min_clpp", type = "double", default = 0.01,
                  help = "Minimum CLPP to report [default %default]")
    )
  )
}

# ---- read a FastENLOC file into long variant x signal PIP tibble -----------
# `mode` = "gwas" (split on '='), or "qtl" (split each '|' piece on '@=')
read_fastenloc <- function(path, mode = c("qtl", "gwas")) {
  mode <- match.arg(mode)

  # fread is much faster than read_tsv on the large QTL panels; keep only the
  # two columns we need (variant id + annotation), as character. fread reads
  # .gz natively, but piping through `gzip -dc` is portable across builds.
  fread_args <- list(
    header = FALSE, sep = "\t",
    select = c(3, 6),
    col.names = c("variant_id", "annot"),
    colClasses = "character",
    showProgress = FALSE
  )
  if (grepl("\\.gz$", path)) {
    fread_args$cmd <- paste("gzip -dc", shQuote(path))
  } else {
    fread_args$input <- path
  }
  raw <- do.call(data.table::fread, fread_args) %>%
    as_tibble()

  if (mode == "qtl") {
    df <- raw %>%
      # one row per (variant, signal): split the '|'-delimited signal list
      separate_rows(annot, sep = "\\|") %>%
      mutate(
        signal = str_extract(annot, "^[^=]+") %>% str_remove("@$"),
        pip    = as.numeric(str_match(annot, "@=([0-9eE.+-]+)")[, 2])
      )
  } else {
    df <- raw %>%
      mutate(
        signal = str_extract(annot, "^[^=]+"),
        pip    = as.numeric(str_match(annot, "=([0-9eE.+-]+)\\[")[, 2])
      )
  }

  df %>%
    filter(!is.na(pip), !is.na(signal)) %>%
    transmute(variant_id, signal, pip)
}

# ---- CLPP: join on shared variants, sum PIP products per signal pair -------
compute_clpp <- function(gwas_long, qtl_long, min_clpp = 0.01) {
  # extract ENSG gene id from the qtl signal (portion before first ':')
  gene_of <- function(sig) {
    g <- str_extract(sig, "^[^:]+")
    ifelse(str_detect(g, "^ENSG"), str_extract(g, "^ENSG[0-9]+"), g)
  }

  inner_join(
    gwas_long %>% rename(gwas_cs = signal, pip_gwas = pip),
    qtl_long  %>% rename(qtl_sig = signal, pip_qtl = pip),
    by = "variant_id",
    relationship = "many-to-many"
  ) %>%
    group_by(gwas_cs, qtl_sig) %>%
    summarise(
      n_shared = n(),
      CLPP     = sum(pip_gwas * pip_qtl),
      .groups  = "drop"
    ) %>%
    filter(CLPP >= min_clpp) %>%
    mutate(gene = gene_of(qtl_sig)) %>%
    arrange(desc(CLPP)) %>%
    select(gwas_cs, qtl_sig, gene, n_shared, CLPP)
}

# ---- main ------------------------------------------------------------------
main <- function() {
  a <- parse_args(build_option_parser())
  if (is.null(a$gwas) || is.null(a$qtl) || is.null(a$out)) {
    print_help(build_option_parser())
    stop("--gwas, --qtl and --out are all required.", call. = FALSE)
  }

  message("Reading GWAS: ", a$gwas)
  gwas_long <- read_fastenloc(a$gwas, "gwas")
  message("  ", n_distinct(gwas_long$signal), " credible sets, ",
          nrow(gwas_long), " variant rows")

  message("Reading QTL:  ", a$qtl)
  qtl_long <- read_fastenloc(a$qtl, "qtl")
  message("  ", n_distinct(qtl_long$signal), " signals, ",
          nrow(qtl_long), " variant rows")

  res <- compute_clpp(gwas_long, qtl_long, a$min_clpp)
  # readr::write_tsv gzips automatically when the path ends in .gz
  write_tsv(res, a$out)
  message("Wrote ", nrow(res), " CLPP pairs (>= ", a$min_clpp, ") to ", a$out)
}

if (identical(environment(), globalenv()) &&
    length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
}
