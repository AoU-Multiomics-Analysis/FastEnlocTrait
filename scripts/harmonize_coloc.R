#!/usr/bin/env Rscript
# ============================================================================
# harmonize_coloc.R
#
# Harmonize FastENLOC and CLPP colocalization outputs into ONE tidy table,
# so downstream analysis (thresholding, coverage, cross-metric comparison)
# is a single group_by/summarise instead of re-parsing three formats.
#
# Three input sources, three native formats:
#   1. FastENLOC .enloc.sig.out  -- signal level; whitespace-delimited.
#        col1 "Signal" = "<qtl_sig>(@)<gwas_cs>=<gwas_pip>[...]" (gwas part
#        empty when the QTL signal maps to no GWAS credible set); col6 = RCP,
#        col7 = LCP.  Key: (gwas_cs, qtl_sig).  Metric: RCP (signal coloc).
#   2. FastENLOC .enloc.gene.out -- gene level; "Gene  GRCP  GLCP".
#        Key: gene.  Metric: GLCP (locus-level regional coloc per gene).
#   3. CLPP pairs (from clpp_fastenloc.R) -- signal level; tsv(.gz) or csv,
#        cols gwas_cs, qtl_sig, gene, n_shared, CLPP.  Key: (gwas_cs, qtl_sig).
#
# THREE OUTPUT GRANULARITIES (same results, different row unit):
#   --out       signal level: one row per (gwas_cs, qtl_sig, gene), every
#               metric at native resolution + FDR pass flags. For "which gene,
#               which metric" and cross-metric comparison.
#   --cs_out    credible-set level: one row per gwas_cs (multi-gene loci
#               collapsed), coloc_* flags + per-method gene lists. For counting
#               colocalizing credible sets and coverage rates. Pass --gwas so
#               non-colocalizing credible sets are included -> correct denom.
#   --gene_out  gene level: one row per gene, GRCP/GLCP native + best RCP/CLPP +
#               how many credible sets the gene colocalizes across.
# All emitted as gzip TSV when the path ends in .gz.
#
# NOTE: any_coloc = union of the STRINGENT metrics (RCP>=0.5, CLPP>=0.05,
#       RCP-FDR, GLCP-FDR). The lenient CLPP>=0.01 screening threshold is
#       reported in its own column but NOT folded into any_coloc.
#
# Usage:
#   Rscript harmonize_coloc.R --sig S.enloc.sig.out --gene S.enloc.gene.out \
#           --clpp pairs.tsv.gz --gwas GWAS.vcf.gz \
#           --out harmonized.tsv.gz \
#           --cs_out harmonized.cs.tsv.gz --gene_out harmonized.gene.tsv.gz \
#           [--study GCST... --trait "Rheumatoid arthritis" \
#            --trait_category immune --n_variants 1000000 \
#            --n_credible_sets 250 --layer eQTL \
#            --consensus_map consensus_loci.tsv.gz \
#            --fdr_level 0.05]
#
# Inputs: FastENLOC .sig.out (RCP/LCP) + .gene.out (GRCP/GLCP), and optionally
# the CLPP pairs file from clpp_fastenloc.R. Only --sig, --gene, --out required.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

# ---- options ---------------------------------------------------------------
build_option_parser <- function() {
  OptionParser(
    usage = "Rscript %prog --sig sig.out --gene gene.out --clpp pairs.tsv.gz --out harmonized.tsv.gz [options]",
    description = paste(
      "Join FastENLOC sig.out (RCP), gene.out (GLCP) and CLPP pairs into one",
      "tidy per-signal table for downstream colocalization analysis."
    ),
    option_list = list(
      make_option("--sig",  type = "character", default = NULL,
                  help = "FastENLOC .enloc.sig.out (RCP/LCP) [required]"),
      make_option("--gene", type = "character", default = NULL,
                  help = "FastENLOC .enloc.gene.out (GRCP/GLCP) [required]"),
      make_option("--clpp", type = "character", default = NULL,
                  help = "CLPP pairs file, tsv/csv (+.gz) [optional]"),
      make_option("--out",  type = "character", default = NULL,
                  help = "Output path; gzipped if it ends in .gz [required]"),
      make_option("--study", type = "character", default = NA,
                  help = "Study id to stamp on every row [optional]"),
      make_option("--trait", type = "character", default = NA,
                  help = "Trait label to stamp on every row [optional]"),
      make_option("--trait_category", type = "character", default = NA,
                  help = "Trait category to stamp on every row [optional]"),
      make_option("--n_variants", type = "integer", default = NA,
                  help = "GWAS variant count denominator to stamp on every row [optional]"),
      make_option("--n_credible_sets", type = "integer", default = NA,
                  help = "Total GWAS credible-set denominator to stamp on every row [optional]"),
      make_option("--layer", type = "character", default = NA,
                  help = "QTL layer label (eQTL/sQTL/pQTL) [optional]"),
      make_option("--fdr_level", type = "double", default = 0.05,
                  help = "Bayesian FDR control level for pass flags [default %default]"),
      make_option("--gwas", type = "character", default = NULL,
                  help = paste("FastENLOC-format GWAS file; enumerates ALL credible",
                               "sets so the CS rollup denominator is complete [optional]")),
      make_option("--cs_out", type = "character", default = NULL,
                  help = paste("If set, also write a credible-set-level rollup",
                               "(one row per gwas_cs) to this path [optional]")),
      make_option("--gene_out", type = "character", default = NULL,
                  help = paste("If set, also write a gene-level rollup",
                               "(one row per gene) to this path [optional]")),
      make_option("--consensus_map", type = "character", default = NULL,
                  help = paste("Consensus-locus map from merge_credible_sets.R",
                               "(tsv/csv +.gz). If set, the CS rollup gains",
                               "consensus_locus_id and merge metadata so",
                               "coverage can be de-duplicated across studies",
                               "[optional]"))
    )
  )
}

# ---- Bayesian FDR threshold (FastENLOC bayesian_FDR.html procedure) ---------
# lfdr = sort(1 - prob); FDR = cumsum(lfdr) / seq_along(lfdr);
# threshold = 1 - lfdr[max(which(FDR <= level))]; Inf if none pass.
# A value "passes" if it is >= threshold.
bfdr_threshold <- function(probs, level = 0.05) {
  probs <- as.numeric(probs)
  probs <- probs[!is.na(probs)]
  if (length(probs) == 0) return(Inf)
  lfdr <- sort(1 - probs)
  fdr  <- cumsum(lfdr) / seq_along(lfdr)
  ok   <- which(fdr <= level)
  if (length(ok) == 0) return(Inf)
  1 - lfdr[max(ok)]
}

# helper: fread that reads plain or .gz uniformly, returns tibble
read_any <- function(path, ...) {
  args <- list(showProgress = FALSE, ...)
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                        args$input <- path
  as_tibble(do.call(data.table::fread, args))
}

gene_of <- function(sig) {
  g <- str_extract(sig, "^[^:]+")
  ifelse(str_detect(g, "^ENSG"), str_extract(g, "^ENSG[0-9]+"), g)
}

trait_of_gwas_cs <- function(gwas_cs) {
  str_extract(gwas_cs, "^[^;]+(?=;)")
}

metadata_value <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x) || x == "") {
    NA_character_
  } else {
    as.character(x)
  }
}

metadata_int <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x) || as.character(x) == "") {
    NA_integer_
  } else {
    as.integer(x)
  }
}

field_at <- function(fields, index) {
  if (is.na(index) || length(fields) < index) {
    NA_character_
  } else {
    fields[[index]]
  }
}

top_value_by_score <- function(value, score) {
  value <- as.character(value)
  score <- suppressWarnings(as.numeric(score))
  keep <- !is.na(value) & !is.na(score)
  if (!any(keep)) return(NA_character_)
  value[keep][which.max(score[keep])]
}

# ---- 1. sig.out (RCP) ------------------------------------------------------
# whitespace-delimited; parse the Signal token into qtl_sig + gwas_cs.
read_sig <- function(path) {
  raw <- read_lines(path)
  raw <- raw[!str_detect(raw, "^(trait\\s+)?Signal\\b")] # drop header
  raw <- raw[str_detect(raw, "\\(@\\)")]                 # keep real signal rows
  tibble(line = raw) %>%
    mutate(
      fields  = str_split(str_trim(line), "\\s+"),
      signal_index = map_int(fields, function(x) {
        hit <- which(str_detect(x, "\\(@\\)"))
        if (length(hit) == 0) NA_integer_ else hit[[1]]
      })
    ) %>%
    filter(!is.na(signal_index)) %>%
    mutate(
      source_trait = map2_chr(fields, signal_index, ~ if (.y > 1) .x[[1]] else NA_character_),
      signal  = map2_chr(fields, signal_index, field_at),
      RCP     = map2_dbl(fields, signal_index, ~ suppressWarnings(as.numeric(field_at(.x, .y + 5)))),
      LCP     = map2_dbl(fields, signal_index, ~ suppressWarnings(as.numeric(field_at(.x, .y + 6)))),
      qtl_sig = str_remove(str_extract(signal, "^.*\\(@\\)"), "\\(@\\)$"),
      gwas_part = str_remove(signal, "^.*\\(@\\)"),
      gwas_cs = na_if(str_extract(gwas_part, "^[^=]+"), "")
    ) %>%
    filter(!is.na(gwas_cs)) %>%           # keep signals mapped to a GWAS cs
    transmute(source_trait = coalesce(source_trait, trait_of_gwas_cs(gwas_cs)),
              gwas_cs, qtl_sig, gene = gene_of(qtl_sig), RCP, LCP) %>%
    # Pin column types so empty tables still join cleanly downstream.
    mutate(source_trait = as.character(source_trait),
           gwas_cs = as.character(gwas_cs), qtl_sig = as.character(qtl_sig),
           gene = as.character(gene),
           RCP = as.numeric(RCP), LCP = as.numeric(LCP))
}

# ---- 2. gene.out (GRCP/GLCP) ----------------------------------------------
# gene.out is "Gene<TAB><TAB>GRCP<TAB>GLCP" -- the double tab yields an empty
# field, so split on runs of whitespace and take the token + the last two
# numeric values (GRCP, GLCP) explicitly rather than by fixed column index.
read_gene <- function(path) {
  raw <- read_lines(path)
  raw <- raw[!str_detect(raw, "^(trait\\s+)?Gene\\b")] # drop header
  raw <- raw[str_trim(raw) != ""]                  # drop blanks
  tibble(line = raw) %>%
    mutate(
      fields = str_split(str_trim(line), "\\s+"),
      source_trait = map_chr(fields, ~ if (length(.x) >= 4) .x[[1]] else NA_character_),
      gene   = map_chr(fields, ~ if (length(.x) >= 4) .x[[2]] else .x[[1]]),
      GRCP   = map_dbl(fields, ~ suppressWarnings(as.numeric(.x[length(.x) - 1]))),
      GLCP   = map_dbl(fields, ~ suppressWarnings(as.numeric(.x[length(.x)])))
    ) %>%
    filter(!is.na(GLCP)) %>%
    transmute(source_trait = as.character(source_trait),
              gene = as.character(gene_of(gene)),
              GRCP = as.numeric(GRCP), GLCP = as.numeric(GLCP)) %>%
    distinct(source_trait, gene, .keep_all = TRUE)
}

# ---- 3. CLPP pairs ---------------------------------------------------------
read_clpp <- function(path) {
  read_any(path) %>%
    transmute(source_trait = as.character(trait_of_gwas_cs(gwas_cs)),
              gwas_cs = as.character(gwas_cs),
              qtl_sig = as.character(qtl_sig),
              n_shared = as.integer(n_shared),
              CLPP = as.numeric(CLPP))
}

# ---- consensus-locus map (from merge_credible_sets.R) ----------------------
read_consensus <- function(path) {
  read_any(path) %>%
    transmute(study_id = as.character(study_id),
              gwas_cs = as.character(gwas_cs),
              consensus_locus_id = as.character(consensus_locus_id),
              n_studies_in_locus = as.integer(n_studies_in_locus),
              n_cs_in_locus = as.integer(n_cs_in_locus),
              is_merged = as.logical(is_merged)) %>%
    distinct(study_id, gwas_cs, .keep_all = TRUE)
}

# ---- harmonize -------------------------------------------------------------
harmonize <- function(sig, gene, clpp = NULL, study = NA, trait = NA,
                      trait_category = NA, n_variants = NA_integer_,
                      n_credible_sets = NA_integer_, layer = NA,
                      fdr_level = 0.05) {
  # FDR thresholds computed at each metric's NATIVE granularity over the full
  # source set: RCP over all mapped signals in sig.out; GRCP/GLCP over all
  # distinct genes in gene.out. (Not over the joined rows -- that would
  # re-weight by how many signals map to each gene.)
  rcp_thr  <- bfdr_threshold(sig$RCP,   fdr_level)
  grcp_thr <- bfdr_threshold(gene$GRCP, fdr_level)
  glcp_thr <- bfdr_threshold(gene$GLCP, fdr_level)
  message(sprintf("FDR %.0f%% thresholds  RCP=%.4g  GRCP=%.4g  GLCP=%.4g",
                  100 * fdr_level, rcp_thr, grcp_thr, glcp_thr))

  out <- sig
  if (!is.null(clpp)) {
    out <- full_join(out, clpp, by = c("gwas_cs", "qtl_sig")) %>%
      mutate(
        source_trait = coalesce(as.character(source_trait.x),
                                as.character(source_trait.y),
                                as.character(trait_of_gwas_cs(gwas_cs))),
        gene = coalesce(as.character(gene), as.character(gene_of(qtl_sig)))
      ) %>%
      select(-source_trait.x, -source_trait.y)
  } else {
    out <- out %>% mutate(n_shared = NA_integer_, CLPP = NA_real_)
  }

  gene_join_by <- if (any(!is.na(gene$source_trait)) && any(!is.na(out$source_trait))) {
    c("source_trait", "gene")
  } else {
    "gene"
  }
  gene_for_join <- if (identical(gene_join_by, "gene")) {
    select(gene, -source_trait)
  } else {
    gene
  }
  study_value <- metadata_value(study)
  trait_value <- metadata_value(trait)
  trait_category_value <- metadata_value(trait_category)
  n_variants_value <- metadata_int(n_variants)
  n_credible_sets_value <- metadata_int(n_credible_sets)
  layer_value <- metadata_value(layer)

  out %>%
    left_join(gene_for_join, by = gene_join_by) %>%
    mutate(
      study = study_value,
      trait = if (is.na(trait_value)) source_trait else trait_value,
      trait_category = trait_category_value,
      n_variants = n_variants_value,
      n_credible_sets = n_credible_sets_value,
      layer = layer_value,
      RCP_pass_FDR  = !is.na(RCP)  & RCP  >= rcp_thr,
      GRCP_pass_FDR = !is.na(GRCP) & GRCP >= grcp_thr,
      GLCP_pass_FDR = !is.na(GLCP) & GLCP >= glcp_thr
    ) %>%
    select(study, trait, trait_category, n_variants, n_credible_sets,
           layer, gwas_cs, qtl_sig, gene,
           RCP, LCP, CLPP, n_shared, GRCP, GLCP,
           RCP_pass_FDR, GRCP_pass_FDR, GLCP_pass_FDR) %>%
    arrange(desc(coalesce(RCP, 0)), desc(coalesce(CLPP, 0)))
}

# ---- enumerate all GWAS credible sets (for a complete denominator) ---------
# reads the FastENLOC GWAS file and returns every distinct credible-set id,
# including those that colocalize with nothing.
all_gwas_cs <- function(path) {
  read_any(path, header = FALSE, select = 6, col.names = "annot",
           colClasses = "character") %>%
    mutate(gwas_cs = str_extract(annot, "^[^=]+")) %>%
    filter(!is.na(gwas_cs)) %>%
    transmute(gwas_cs, trait = trait_of_gwas_cs(gwas_cs)) %>%
    distinct(gwas_cs, trait)
}

# ---- credible-set-level rollup ---------------------------------------------
# collapse the signal-level table to ONE row per credible set, so credible
# sets can be counted directly. `all_cs` (optional) supplies the full set of
# credible sets so non-colocalizing ones appear with 0/FALSE.
rollup_cs <- function(signal_tbl, all_cs = NULL, consensus = NULL,
                      study = NA, trait = NA, trait_category = NA,
                      n_variants = NA_integer_,
                      n_credible_sets = NA_integer_, layer = NA) {
  study_value <- metadata_value(study)
  trait_value <- metadata_value(trait)
  trait_category_value <- metadata_value(trait_category)
  n_variants_value <- metadata_int(n_variants)
  n_credible_sets_value <- metadata_int(n_credible_sets)
  layer_value <- metadata_value(layer)

  cs <- signal_tbl %>%
    group_by(study, trait, trait_category, n_variants, n_credible_sets,
             layer, gwas_cs) %>%
    summarise(
      n_qtl_signals   = sum(!is.na(qtl_sig)),
      n_genes         = n_distinct(gene[!is.na(gene)]),
      best_RCP        = if (all(is.na(RCP)))  NA_real_ else max(RCP,  na.rm = TRUE),
      best_CLPP       = if (all(is.na(CLPP))) NA_real_ else max(CLPP, na.rm = TRUE),
      best_GLCP       = if (all(is.na(GLCP))) NA_real_ else max(GLCP, na.rm = TRUE),
      coloc_RCP_0.5   = any(RCP  >= 0.5,  na.rm = TRUE),
      coloc_CLPP_0.05 = any(CLPP >= 0.05, na.rm = TRUE),
      coloc_CLPP_0.01 = any(CLPP >= 0.01, na.rm = TRUE),
      coloc_RCP_FDR   = any(RCP_pass_FDR,  na.rm = TRUE),
      coloc_GLCP_FDR  = any(GLCP_pass_FDR, na.rm = TRUE),
      # per-method lists of the genes this credible set colocalizes with
      genes_RCP_0.5   = paste(sort(unique(gene[RCP  >= 0.5  & !is.na(RCP)])),  collapse = ";"),
      genes_CLPP_0.05 = paste(sort(unique(gene[CLPP >= 0.05 & !is.na(CLPP)])), collapse = ";"),
      genes_CLPP_0.01 = paste(sort(unique(gene[CLPP >= 0.01 & !is.na(CLPP)])), collapse = ";"),
      genes_RCP_FDR   = paste(sort(unique(gene[RCP_pass_FDR])),  collapse = ";"),
      genes_GLCP_FDR  = paste(sort(unique(gene[GLCP_pass_FDR])), collapse = ";"),
      top_gene        = top_value_by_score(gene, RCP),
      .groups = "drop"
    )
  # add non-colocalizing credible sets (present in GWAS, absent from signals)
  if (!is.null(all_cs)) {
    missing <- anti_join(all_cs, cs, by = "gwas_cs")
    if (nrow(missing) > 0) {
      cs <- bind_rows(cs, missing %>% mutate(
        study = study_value,
        trait = if (is.na(trait_value)) trait else trait_value,
        trait_category = trait_category_value,
        n_variants = n_variants_value,
        n_credible_sets = n_credible_sets_value,
        layer = layer_value,
        n_qtl_signals = 0L, n_genes = 0L,
        best_RCP = NA_real_, best_CLPP = NA_real_, best_GLCP = NA_real_,
        coloc_RCP_0.5 = FALSE, coloc_CLPP_0.05 = FALSE, coloc_CLPP_0.01 = FALSE,
        coloc_RCP_FDR = FALSE, coloc_GLCP_FDR = FALSE,
        genes_RCP_0.5 = "", genes_CLPP_0.05 = "", genes_CLPP_0.01 = "",
        genes_RCP_FDR = "", genes_GLCP_FDR = "", top_gene = NA_character_))
    }
  }
  cs <- cs %>%
    mutate(
      # any_coloc = union of the STRINGENT metrics only. CLPP>=0.01 is a lenient
      # screening threshold and is deliberately EXCLUDED here (it is still
      # reported in its own coloc_CLPP_0.01 column).
      any_coloc = coloc_RCP_0.5 | coloc_CLPP_0.05 | coloc_RCP_FDR | coloc_GLCP_FDR
    )

  if (!is.null(consensus)) {
    cs <- cs %>%
      left_join(consensus, by = c("study" = "study_id", "gwas_cs")) %>%
      mutate(
        consensus_locus_id = coalesce(consensus_locus_id, gwas_cs),
        n_studies_in_locus = coalesce(n_studies_in_locus, 1L),
        n_cs_in_locus = coalesce(n_cs_in_locus, 1L),
        is_merged = coalesce(is_merged, FALSE)
      ) %>%
      select(study, trait, trait_category, n_variants, n_credible_sets,
             layer, gwas_cs, consensus_locus_id,
             n_studies_in_locus, n_cs_in_locus, is_merged,
             n_qtl_signals, n_genes,
             best_RCP, best_CLPP, best_GLCP,
             coloc_RCP_0.5, coloc_CLPP_0.05, coloc_CLPP_0.01,
             coloc_RCP_FDR, coloc_GLCP_FDR, any_coloc,
             genes_RCP_0.5, genes_CLPP_0.05, genes_CLPP_0.01,
             genes_RCP_FDR, genes_GLCP_FDR, top_gene)
  } else {
    cs <- cs %>%
      select(study, trait, trait_category, n_variants, n_credible_sets,
             layer, gwas_cs, n_qtl_signals, n_genes,
             best_RCP, best_CLPP, best_GLCP,
             coloc_RCP_0.5, coloc_CLPP_0.05, coloc_CLPP_0.01,
             coloc_RCP_FDR, coloc_GLCP_FDR, any_coloc,
             genes_RCP_0.5, genes_CLPP_0.05, genes_CLPP_0.01,
             genes_RCP_FDR, genes_GLCP_FDR, top_gene)
  }

  cs %>%
    arrange(desc(any_coloc), desc(replace_na(best_RCP, -1)))
}

# ---- gene-level rollup -----------------------------------------------------
# collapse the signal-level table to ONE row per gene. GRCP/GLCP are already
# gene-native (one value per gene, constant across that gene's rows); RCP/CLPP
# are aggregated as the best signal for the gene, with a count of how many
# distinct credible sets the gene colocalizes with.
rollup_gene <- function(signal_tbl, study = NA, trait = NA,
                        trait_category = NA, n_variants = NA_integer_,
                        n_credible_sets = NA_integer_, layer = NA) {
  signal_tbl <- signal_tbl %>%
    filter(!is.na(gene))
  if (nrow(signal_tbl) == 0) {
    return(tibble(
      study = character(), trait = character(), trait_category = character(),
      n_variants = integer(), n_credible_sets = integer(), layer = character(),
      gene = character(), n_gene_credible_sets = integer(),
      n_signals = integer(), best_RCP = numeric(), best_CLPP = numeric(),
      GRCP = numeric(), GLCP = numeric(), coloc_RCP_0.5 = logical(),
      coloc_CLPP_0.05 = logical(), coloc_CLPP_0.01 = logical(),
      RCP_pass_FDR = logical(), GRCP_pass_FDR = logical(),
      GLCP_pass_FDR = logical(), any_coloc = logical(),
      top_cs = character()
    ))
  }

  signal_tbl %>%
    group_by(study, trait, trait_category, n_variants, n_credible_sets,
             layer, gene) %>%
    summarise(
      n_gene_credible_sets = n_distinct(gwas_cs[!is.na(gwas_cs)]),
      n_signals       = sum(!is.na(qtl_sig)),
      best_RCP        = if (all(is.na(RCP)))  NA_real_ else max(RCP,  na.rm = TRUE),
      best_CLPP       = if (all(is.na(CLPP))) NA_real_ else max(CLPP, na.rm = TRUE),
      GRCP            = first(GRCP),   # gene-native: constant within the gene
      GLCP            = first(GLCP),
      coloc_RCP_0.5   = any(RCP  >= 0.5,  na.rm = TRUE),
      coloc_CLPP_0.05 = any(CLPP >= 0.05, na.rm = TRUE),
      coloc_CLPP_0.01 = any(CLPP >= 0.01, na.rm = TRUE),
      RCP_pass_FDR    = any(RCP_pass_FDR,  na.rm = TRUE),
      GRCP_pass_FDR   = first(GRCP_pass_FDR),
      GLCP_pass_FDR   = first(GLCP_pass_FDR),
      top_cs          = top_value_by_score(gwas_cs, RCP),
      .groups = "drop"
    ) %>%
    mutate(
      # union of stringent metrics (CLPP>=0.01 screening threshold excluded)
      any_coloc = coloc_RCP_0.5 | coloc_CLPP_0.05 | RCP_pass_FDR | GLCP_pass_FDR
    ) %>%
    select(study, trait, trait_category, n_variants, n_credible_sets,
           layer, gene, n_gene_credible_sets, n_signals,
           best_RCP, best_CLPP, GRCP, GLCP,
           coloc_RCP_0.5, coloc_CLPP_0.05, coloc_CLPP_0.01,
           RCP_pass_FDR, GRCP_pass_FDR, GLCP_pass_FDR, any_coloc, top_cs) %>%
    arrange(desc(any_coloc), desc(replace_na(GLCP, -1)))
}

# ---- main ------------------------------------------------------------------
main <- function() {
  a <- parse_args(build_option_parser())
  if (is.null(a$sig) || is.null(a$gene) || is.null(a$out)) {
    print_help(build_option_parser())
    stop("--sig, --gene and --out are required.", call. = FALSE)
  }
  sig  <- read_sig(a$sig)
  gene <- read_gene(a$gene)
  clpp <- if (!is.null(a$clpp)) read_clpp(a$clpp) else NULL
  res  <- harmonize(
    sig = sig,
    gene = gene,
    clpp = clpp,
    study = a$study,
    trait = a$trait,
    trait_category = a$trait_category,
    n_variants = a$n_variants,
    n_credible_sets = a$n_credible_sets,
    layer = a$layer,
    fdr_level = a$fdr_level
  )
  write_tsv(res, a$out)
  message("Harmonized ", nrow(res), " signal rows -> ", a$out,
          "  (", n_distinct(res$gwas_cs), " credible sets, ",
          n_distinct(res$gene), " genes)")

  # optional credible-set-level rollup
  if (!is.null(a$cs_out)) {
    all_cs <- if (!is.null(a$gwas)) all_gwas_cs(a$gwas) else NULL
    consensus <- if (!is.null(a$consensus_map)) read_consensus(a$consensus_map) else NULL
    cs <- rollup_cs(
      signal_tbl = res,
      all_cs = all_cs,
      consensus = consensus,
      study = a$study,
      trait = a$trait,
      trait_category = a$trait_category,
      n_variants = a$n_variants,
      n_credible_sets = a$n_credible_sets,
      layer = a$layer
    )
    write_tsv(cs, a$cs_out)
    n_total <- nrow(cs)
    n_any   <- sum(cs$any_coloc)
    message("Rolled up to ", n_total, " credible sets -> ", a$cs_out,
            if (is.null(all_cs)) "  (colocalizing CS only; pass --gwas for full denominator)"
            else sprintf("  (%d/%d = %.1f%% colocalize by any metric)",
                         n_any, n_total, 100 * n_any / n_total))
    if (!is.null(consensus)) {
      loci <- cs %>%
        group_by(consensus_locus_id) %>%
        summarise(coloc = any(any_coloc), .groups = "drop")
      message(sprintf("  consensus loci: %d/%d = %.1f%% colocalize",
                      sum(loci$coloc), nrow(loci), 100 * mean(loci$coloc)))
    }
  }

  # optional gene-level rollup
  if (!is.null(a$gene_out)) {
    gtbl <- rollup_gene(
      signal_tbl = res,
      study = a$study,
      trait = a$trait,
      trait_category = a$trait_category,
      n_variants = a$n_variants,
      n_credible_sets = a$n_credible_sets,
      layer = a$layer
    )
    write_tsv(gtbl, a$gene_out)
    message("Rolled up to ", nrow(gtbl), " genes -> ", a$gene_out,
            sprintf("  (%d colocalize by any metric)", sum(gtbl$any_coloc)))
  }
}

if (identical(environment(), globalenv()) &&
    length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
}
