#!/usr/bin/env Rscript
# ============================================================================
# summarize_coloc.R
#
# Cross-trait colocalization summary that operates directly on the WDL's
# aggregated harmonized outputs (harmonized_coloc.cs.tsv.gz and
# harmonized_coloc.gene.tsv.gz -- the workflow's harmonized_credible_set_out
# and harmonized_gene_out). One run summarizes every trait x every QTL layer
# PLUS the cross-layer union, into a small number of tidy analysis files.
#
# Two families of question, two output files:
#
#  A) --coloc_rate_out : per (trait, layer|union) CREDIBLE-SET coverage.
#     "What fraction of a trait's de-duplicated GWAS credible sets colocalize?"
#     The denominator is CONSENSUS LOCI (consensus_locus_id) so credible sets
#     shared across studies of the same trait are counted once. Coverage is
#     reported at five thresholds:
#         GLCP_FDR   coloc_GLCP_FDR   (FastENLOC gene-level, Bayesian FDR)
#         RCP_FDR    coloc_RCP_FDR    (FastENLOC signal-level, Bayesian FDR)
#         RCP_0.5    coloc_RCP_0.5    (FastENLOC signal RCP >= 0.5)
#         CLPP_0.01  coloc_CLPP_0.01  (eCAVIAR CLPP >= 0.01)
#         CLPP_0.05  coloc_CLPP_0.05  (eCAVIAR CLPP >= 0.05)
#     plus two composite multi-metric unions:
#         stringent_union  RCP_FDR OR CLPP_0.05   (either stricter metric)
#         lenient_union    RCP_0.5 OR CLPP_0.01   (either looser metric)
#     A consensus locus counts as colocalizing under a threshold if ANY of its
#     member credible sets (across studies) meets it.
#
#  B) --gene_summary_out : per (trait, layer|union) GENE-level description.
#     number of genes colocalizing, number per locus, number of PROTEIN-CODING
#     genes per locus, and the colocalizing gene id list. A gene counts if it
#     colocalizes under ANY of the five thresholds (see --gene_threshold to
#     restrict). Protein-coding status and gene names come from a GENCODE GTF
#     (--gtf): gene-feature lines are parsed for gene_id/gene_type/gene_name,
#     ENSG versions stripped. Genes absent from the GTF are non-coding/unknown.
#
# The "union" layer for a trait pools all QTL layers: a consensus locus is
# covered in the union if it is covered in ANY layer; a gene colocalizes in the
# union if it colocalizes in any layer.
#
# INPUT is the WDL's aggregated files -- no per-study wrangling needed. The CS
# file must carry: trait, layer, consensus_locus_id, coloc_* flags, and the
# per-method gene-list columns (genes_RCP_0.5, genes_CLPP_0.05, genes_CLPP_0.01,
# genes_RCP_FDR, genes_GLCP_FDR). The gene file is required so the WDL summary
# stage depends on both harmonized outputs; this script derives gene-per-locus
# summaries from the CS gene-list columns so gene counts use the same
# consensus-locus denominator as the coverage table.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

THRESHOLDS <- c(
  GLCP_FDR  = "coloc_GLCP_FDR",
  RCP_FDR   = "coloc_RCP_FDR",
  RCP_0.5   = "coloc_RCP_0.5",
  CLPP_0.01 = "coloc_CLPP_0.01",
  CLPP_0.05 = "coloc_CLPP_0.05"
)

COMPOSITES <- list(
  stringent_union = c("coloc_RCP_FDR", "coloc_CLPP_0.05"),
  lenient_union   = c("coloc_RCP_0.5", "coloc_CLPP_0.01")
)

GENE_LIST_COLS <- c(
  RCP_0.5   = "genes_RCP_0.5",
  RCP_FDR   = "genes_RCP_FDR",
  GLCP_FDR  = "genes_GLCP_FDR",
  CLPP_0.05 = "genes_CLPP_0.05",
  CLPP_0.01 = "genes_CLPP_0.01"
)

build_option_parser <- function() {
  OptionParser(
    usage = "Rscript %prog --cs harmonized_coloc.cs.tsv.gz --gene harmonized_coloc.gene.tsv.gz [options]",
    description = "Cross-trait colocalization summary from WDL harmonized outputs.",
    option_list = list(
      make_option("--cs", type = "character", default = NULL,
                  help = "Aggregated harmonized CREDIBLE-SET file (.tsv/.gz) [required]"),
      make_option("--gene", type = "character", default = NULL,
                  help = "Aggregated harmonized GENE file (.tsv/.gz) [required]"),
      make_option("--gtf", type = "character", default = NULL,
                  help = "GENCODE GTF (.gtf/.gz) for protein-coding status + gene names [optional]"),
      make_option("--coloc_rate_out", type = "character", default = "coloc_rate_by_trait.tsv",
                  help = "Per (trait, layer|union) CS coverage table [default %default]"),
      make_option("--gene_summary_out", type = "character", default = "gene_summary_by_trait.tsv",
                  help = "Per (trait, layer|union) gene-level summary [default %default]"),
      make_option("--gene_list_out", type = "character", default = NULL,
                  help = "Optional long table: one row per colocalizing (trait,layer|union,gene) [optional]"),
      make_option("--gene_threshold", type = "character", default = "any",
                  help = paste0("Which threshold defines a colocalizing gene for the gene summary: ",
                                "one of any/GLCP_FDR/RCP_FDR/RCP_0.5/CLPP_0.01/CLPP_0.05 [default %default]"))
    )
  )
}

read_any <- function(path, ...) {
  args <- list(showProgress = FALSE, ...)
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                        args$input <- path
  as_tibble(do.call(data.table::fread, args))
}

require_columns <- function(tbl, cols, label) {
  missing <- setdiff(cols, names(tbl))
  if (length(missing) > 0) {
    stop(label, " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  if (is.numeric(x)) return(replace_na(x != 0, FALSE))
  v <- tolower(trimws(as.character(x)))
  case_when(
    v %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    v %in% c("false", "f", "0", "no", "n", "", "na", "nan") ~ FALSE,
    TRUE ~ FALSE
  )
}

read_gtf_biotype <- function(path) {
  awk <- "awk -F '\\t' 'BEGIN{OFS=\"\\t\"} $0 !~ /^#/ && $3 == \"gene\" {print $3,$9}'"
  cmd <- if (grepl("\\.gz$", path)) {
    paste("gzip -dc", shQuote(path), "|", awk)
  } else {
    paste(awk, shQuote(path))
  }
  gtf <- as_tibble(data.table::fread(
    cmd = cmd,
    header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c("feature", "attr"),
    showProgress = FALSE
  ))
  pull_attr <- function(a, key) {
    str_match(a, paste0(key, ' "([^"]+)"'))[, 2]
  }
  tibble(
    gene_id = str_remove(pull_attr(gtf$attr, "gene_id"), "\\.\\d+$"),
    gene_type = pull_attr(gtf$attr, "gene_type"),
    gene_name = pull_attr(gtf$attr, "gene_name")
  ) %>%
    filter(!is.na(gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE)
}

split_genes <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(character(0))
  g <- unlist(str_split(x, ";"))
  unique(g[g != "" & !is.na(g)])
}

coverage_table <- function(cs, thresholds = THRESHOLDS) {
  locus_layer <- cs %>%
    group_by(trait, layer, consensus_locus_id) %>%
    summarise(across(all_of(unname(thresholds)),
                     ~ any(as.logical(.x), na.rm = TRUE)),
              .groups = "drop")

  per_layer <- locus_layer %>%
    group_by(trait, layer) %>%
    summarise(n_consensus_loci = n_distinct(consensus_locus_id),
              across(all_of(unname(thresholds)), ~ sum(.x)),
              .groups = "drop")

  union_tbl <- locus_layer %>%
    group_by(trait, consensus_locus_id) %>%
    summarise(across(all_of(unname(thresholds)), ~ any(.x)), .groups = "drop") %>%
    group_by(trait) %>%
    summarise(layer = "union",
              n_consensus_loci = n_distinct(consensus_locus_id),
              across(all_of(unname(thresholds)), ~ sum(.x)),
              .groups = "drop")

  counts <- bind_rows(per_layer, union_tbl)
  for (lab in names(thresholds)) {
    col <- thresholds[[lab]]
    counts[[paste0("n_", lab)]] <- counts[[col]]
    counts[[paste0("pct_", lab)]] <- round(100 * counts[[col]] / counts$n_consensus_loci, 2)
    counts[[col]] <- NULL
  }
  counts %>%
    arrange(trait, factor(layer, levels = c(setdiff(unique(layer), "union"), "union")))
}

gene_summary_table <- function(cs, biotype = NULL, gene_threshold = "any") {
  allowed_thresholds <- c("any", names(GENE_LIST_COLS))
  if (!gene_threshold %in% allowed_thresholds) {
    stop("--gene_threshold must be one of: ",
         paste(allowed_thresholds, collapse = ", "), call. = FALSE)
  }
  cols <- if (gene_threshold == "any") GENE_LIST_COLS else GENE_LIST_COLS[gene_threshold]
  summary_keys <- bind_rows(
    cs %>% distinct(trait, layer),
    cs %>% distinct(trait) %>% mutate(layer = "union")
  ) %>%
    distinct(trait, layer)

  long <- cs %>%
    mutate(.genes = pmap(across(all_of(unname(cols))),
                         ~ split_genes(c(...)))) %>%
    select(trait, layer, consensus_locus_id, .genes) %>%
    unnest_longer(.genes, values_to = "gene") %>%
    filter(!is.na(gene), gene != "") %>%
    distinct(trait, layer, consensus_locus_id, gene)

  pc_set <- if (!is.null(biotype)) {
    biotype$gene_id[biotype$gene_type == "protein_coding"]
  } else {
    character(0)
  }

  union_long <- long %>%
    distinct(trait, consensus_locus_id, gene) %>%
    mutate(layer = "union")

  summary_long <- bind_rows(long, union_long)

  gene_counts <- summary_long %>%
    group_by(trait, layer) %>%
    summarise(
      n_genes_colocalizing = n_distinct(gene),
      n_protein_coding_genes = n_distinct(gene[gene %in% pc_set]),
      n_loci_with_genes = n_distinct(consensus_locus_id),
      gene_ids = paste(sort(unique(gene)), collapse = ";"),
      .groups = "drop"
    )

  locus_stats <- summary_long %>%
    group_by(trait, layer, consensus_locus_id) %>%
    summarise(
      n_genes = n_distinct(gene),
      n_pc = n_distinct(gene[gene %in% pc_set]),
      .groups = "drop"
    ) %>%
    group_by(trait, layer) %>%
    summarise(
      mean_genes_per_locus = round(mean(n_genes), 2),
      max_genes_per_locus = max(n_genes),
      mean_pc_genes_per_locus = round(mean(n_pc), 2),
      max_pc_genes_per_locus = max(n_pc),
      .groups = "drop"
    )

  summary_keys %>%
    left_join(gene_counts, by = c("trait", "layer")) %>%
    left_join(locus_stats, by = c("trait", "layer")) %>%
    mutate(
      n_genes_colocalizing = replace_na(n_genes_colocalizing, 0L),
      n_protein_coding_genes = replace_na(n_protein_coding_genes, 0L),
      n_loci_with_genes = replace_na(n_loci_with_genes, 0L),
      mean_genes_per_locus = replace_na(mean_genes_per_locus, 0),
      max_genes_per_locus = replace_na(max_genes_per_locus, 0L),
      mean_pc_genes_per_locus = replace_na(mean_pc_genes_per_locus, 0),
      max_pc_genes_per_locus = replace_na(max_pc_genes_per_locus, 0L),
      gene_ids = replace_na(gene_ids, "")
    ) %>%
    select(trait, layer, n_genes_colocalizing, n_protein_coding_genes,
           n_loci_with_genes, mean_genes_per_locus, max_genes_per_locus,
           mean_pc_genes_per_locus, max_pc_genes_per_locus, gene_ids) %>%
    arrange(trait, factor(layer, levels = c(setdiff(unique(layer), "union"), "union")))
}

gene_list_long <- function(cs, biotype = NULL) {
  cols <- GENE_LIST_COLS
  base <- cs %>%
    mutate(.genes = pmap(across(all_of(unname(cols))), ~ split_genes(c(...)))) %>%
    select(trait, layer, consensus_locus_id, .genes) %>%
    unnest_longer(.genes, values_to = "gene") %>%
    filter(!is.na(gene), gene != "") %>%
    distinct(trait, layer, consensus_locus_id, gene)
  union <- base %>%
    distinct(trait, consensus_locus_id, gene) %>%
    mutate(layer = "union")
  bind_rows(base, union) %>%
    left_join(if (!is.null(biotype)) {
                biotype
              } else {
                tibble(gene_id = character(), gene_type = character(), gene_name = character())
              },
              by = c("gene" = "gene_id")) %>%
    mutate(protein_coding = replace_na(gene_type == "protein_coding", FALSE)) %>%
    arrange(trait, layer, consensus_locus_id, gene)
}

main <- function() {
  a <- parse_args(build_option_parser())
  if (is.null(a$cs) || is.null(a$gene)) {
    print_help(build_option_parser())
    stop("--cs and --gene are required.", call. = FALSE)
  }
  read_any(a$gene, nrows = 1)
  cs <- read_any(a$cs, colClasses = list(character = "consensus_locus_id"))
  require_columns(
    cs,
    c("trait", "layer", "consensus_locus_id", unname(THRESHOLDS), unname(GENE_LIST_COLS)),
    "--cs"
  )
  for (col in unname(THRESHOLDS)) {
    cs[[col]] <- as_logical_flag(cs[[col]])
  }

  thresholds <- THRESHOLDS
  for (nm in names(COMPOSITES)) {
    pair <- COMPOSITES[[nm]]
    col_name <- paste0("coloc_", nm)
    cs[[col_name]] <- replace_na(cs[[pair[1]]], FALSE) | replace_na(cs[[pair[2]]], FALSE)
    thresholds[nm] <- col_name
  }

  biotype <- if (!is.null(a$gtf)) {
    b <- read_gtf_biotype(a$gtf)
    message("Parsed ", nrow(b), " genes from GTF (",
            sum(b$gene_type == "protein_coding", na.rm = TRUE), " protein_coding)")
    b
  } else {
    NULL
  }

  cov <- coverage_table(cs, thresholds)
  write_tsv(cov, a$coloc_rate_out)
  message("Wrote coverage (", nrow(cov), " trait x layer rows) -> ", a$coloc_rate_out)

  gsum <- gene_summary_table(cs, biotype, a$gene_threshold)
  write_tsv(gsum, a$gene_summary_out)
  message("Wrote gene summary (", nrow(gsum), " rows) -> ", a$gene_summary_out)

  if (!is.null(a$gene_list_out)) {
    gl <- gene_list_long(cs, biotype)
    write_tsv(gl, a$gene_list_out)
    message("Wrote colocalizing gene list (", nrow(gl), " rows) -> ", a$gene_list_out)
  }
}

if (identical(environment(), globalenv()) &&
    length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
}
