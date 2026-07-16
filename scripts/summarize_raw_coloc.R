#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

build_option_parser <- function() {
  OptionParser(
    usage = "Rscript %prog --gene combined.enloc.gene.out --clpp clpp.combined.tsv --qtl_labels qtl_labels.txt [options]",
    description = "Summarize raw fastENLOC and CLPP outputs by QTL layer.",
    option_list = list(
      make_option("--gene", type = "character", default = NULL,
                  help = "All-GWAS combined fastENLOC gene output [required]"),
      make_option("--enrich", type = "character", default = NULL,
                  help = "All-GWAS combined fastENLOC enrich output [required]"),
      make_option("--mi", type = "character", default = NULL,
                  help = "All-GWAS combined fastENLOC mi output [required]"),
      make_option("--sig", type = "character", default = NULL,
                  help = "All-GWAS combined fastENLOC signal output [required]"),
      make_option("--snp", type = "character", default = NULL,
                  help = "All-GWAS combined fastENLOC snp output [required]"),
      make_option("--clpp", type = "character", default = NULL,
                  help = "All-GWAS combined CLPP output [required]"),
      make_option("--qtl_labels", type = "character", default = NULL,
                  help = "Text file with one QTL label per line [required]"),
      make_option("--out", type = "character", default = "raw_coloc_summary.all_qtl.tsv",
                  help = "All-QTL raw summary output TSV [default %default]"),
      make_option("--per_qtl_dir", type = "character", default = "raw_coloc_summary_by_qtl",
                  help = "Directory for per-QTL summary TSV files [default %default]")
    )
  )
}

read_any <- function(path, ...) {
  args <- list(showProgress = FALSE, fill = TRUE, sep = "\t", ...)
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                        args$input <- path
  as_tibble(do.call(data.table::fread, args))
}

read_fastenloc_combined <- function(path, output_type) {
  args <- list(
    header = FALSE,
    sep = "\t",
    skip = 1,
    fill = TRUE,
    quote = "",
    showProgress = FALSE
  )
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                         args$input <- path
  dt <- do.call(data.table::fread, args)
  if (ncol(dt) < 8) {
    stop(output_type, " input has fewer than eight tab-delimited fields: ", path,
         call. = FALSE)
  }

  metadata_names <- c(
    "study_id", "gwas_trait", "trait_category", "n_variants",
    "n_credible_sets", "qtl_label", "trait"
  )
  setnames(dt, names(dt)[seq_along(metadata_names)], metadata_names)
  payload_cols <- names(dt)[-(seq_along(metadata_names))]
  dt[, payload := trimws(do.call(paste, c(.SD, sep = " "))), .SDcols = payload_cols]
  meta <- as_tibble(dt[, c(metadata_names, "payload"), with = FALSE])

  tokens <- data.table::tstrsplit(meta$payload, "[[:space:]]+", perl = TRUE)
  token <- function(i) {
    if (length(tokens) < i) rep(NA_character_, nrow(meta)) else tokens[[i]]
  }
  parsed <- switch(
    output_type,
    fastenloc_gene = tibble(Gene = token(1), GRCP = token(2), GLCP = token(3)),
    fastenloc_mi = tibble(a0 = token(1), a1 = token(2),
                          p_eqtl = token(3), p_gwas = token(4)),
    fastenloc_sig = tibble(
      Signal = token(1), Num_SNP = token(2), CPIP_qtl = token(3),
      CPIP_gwas_marginal = token(4), CPIP_gwas_qtl_prior = token(5),
      RCP = token(6), LCP = token(7)
    ),
    fastenloc_snp = tibble(
      Signal = token(1), SNP = token(2), PIP_qtl = token(3),
      PIP_gwas_marginal = token(4), PIP_gwas_qtl_prior = token(5),
      SCP = token(6)
    ),
    fastenloc_enrich = {
      parts <- str_match(meta$payload, "^(.*)[[:space:]]+([^[:space:]]+)[[:space:]]+([^[:space:]]+)$")
      tibble(term = parts[, 2], estimate = parts[, 3], standard_error = parts[, 4])
    },
    stop("Unsupported fastENLOC output type: ", output_type, call. = FALSE)
  )
  bind_cols(meta %>% select(-payload), parsed)
}

repair_legacy_fastenloc <- function(df, output_type) {
  # Older aggregation tasks prepended metadata with tabs but left the native
  # fastENLOC payload fixed-width. fread therefore placed the entire payload
  # in Signal and shifted/dropped the remaining metrics. Accept those files so
  # historical runs can be summarized, while new runs arrive as proper TSVs.
  if (output_type == "fastenloc_sig" &&
      all(c("Signal", "Num_SNP", "RCP", "LCP") %in% names(df)) &&
      all(is.na(as_num(df$RCP))) && nrow(df) > 0) {
    fields <- str_split(str_trim(as.character(df$Signal)), "\\s+")
    lengths <- lengths(fields)
    if (any(lengths != 6L)) {
      bad <- which(lengths != 6L)[[1]]
      stop("Cannot repair legacy fastENLOC signal row ", bad,
           ": expected 6 whitespace fields in Signal payload, found ",
           lengths[[bad]], call. = FALSE)
    }
    legacy_lcp <- as_num(df$Num_SNP)
    df$Signal <- map_chr(fields, 1L)
    df$Num_SNP <- map_dbl(fields, ~ as_num(.x[[2]]))
    df$CPIP_qtl <- map_dbl(fields, ~ as_num(.x[[3]]))
    df$CPIP_gwas_marginal <- map_dbl(fields, ~ as_num(.x[[4]]))
    df$CPIP_gwas_qtl_prior <- map_dbl(fields, ~ as_num(.x[[5]]))
    df$RCP <- map_dbl(fields, ~ as_num(.x[[6]]))
    df$LCP <- legacy_lcp
  }

  if (output_type == "fastenloc_snp" &&
      all(c("Signal", "SNP", "SCP") %in% names(df)) &&
      all(is.na(df$SNP)) && nrow(df) > 0) {
    fields <- str_split(str_trim(as.character(df$Signal)), "\\s+")
    lengths <- lengths(fields)
    if (any(lengths != 6L)) {
      bad <- which(lengths != 6L)[[1]]
      stop("Cannot repair legacy fastENLOC SNP row ", bad,
           ": expected 6 whitespace fields in Signal payload, found ",
           lengths[[bad]], call. = FALSE)
    }
    df$Signal <- map_chr(fields, 1L)
    df$SNP <- map_chr(fields, 2L)
    df$PIP_qtl <- map_dbl(fields, ~ as_num(.x[[3]]))
    df$PIP_gwas_marginal <- map_dbl(fields, ~ as_num(.x[[4]]))
    df$PIP_gwas_qtl_prior <- map_dbl(fields, ~ as_num(.x[[5]]))
    df$SCP <- map_dbl(fields, ~ as_num(.x[[6]]))
  }

  if (output_type == "fastenloc_gene") {
    unnamed <- names(df)[str_detect(names(df), "^V[0-9]+$")]
    empty <- unnamed[vapply(unnamed, function(col) {
      all(is.na(df[[col]]) | as.character(df[[col]]) == "")
    }, logical(1))]
    if (length(empty) > 0) df <- df %>% select(-all_of(empty))
  }

  df
}

norm_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "_", x))
}

find_col <- function(df, candidates) {
  nm <- norm_name(names(df))
  hit <- match(norm_name(candidates), nm, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit) == 0) NA_character_ else names(df)[hit[[1]]]
}

as_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_n_distinct <- function(x) {
  x <- x[!is.na(x) & x != ""]
  dplyr::n_distinct(x)
}

safe_mean <- function(x) {
  x <- as_num(x)
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- as_num(x)
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

safe_count_ge <- function(x, threshold) {
  x <- as_num(x)
  sum(!is.na(x) & x >= threshold)
}

extract_signal_parts <- function(signal) {
  tibble(signal = as.character(signal)) %>%
    mutate(
      qtl_sig = str_remove(str_extract(signal, "^.*\\(@\\)"), "\\(@\\)$"),
      gwas_part = str_remove(signal, "^.*\\(@\\)"),
      gwas_cs = str_extract(gwas_part, "^[^=]+")
    )
}

generic_numeric_summary <- function(df) {
  metadata <- c("study_id", "gwas_trait", "trait_category", "n_variants",
                "n_credible_sets", "qtl_label", "trait")
  metric_cols <- setdiff(names(df), metadata)
  numeric_cols <- metric_cols[vapply(metric_cols, function(col) {
    any(!is.na(as_num(df[[col]])))
  }, logical(1))]
  if (length(numeric_cols) == 0) {
    return(list(
      numeric_column_count = 0L,
      numeric_columns = "",
      max_numeric_value = NA_real_,
      mean_numeric_value = NA_real_
    ))
  }
  values <- unlist(lapply(numeric_cols, function(col) as_num(df[[col]])), use.names = FALSE)
  list(
    numeric_column_count = length(numeric_cols),
    numeric_columns = paste(numeric_cols, collapse = ";"),
    max_numeric_value = safe_max(values),
    mean_numeric_value = safe_mean(values)
  )
}

count_gwas_credible_sets <- function(df, gwas_cs_col, signal_parts) {
  if (!is.na(gwas_cs_col)) return(safe_n_distinct(df[[gwas_cs_col]]))
  if (nrow(signal_parts) > 0) return(safe_n_distinct(signal_parts$gwas_cs))
  NA_integer_
}

count_qtl_signals <- function(df, qtl_sig_col, signal_parts) {
  if (!is.na(qtl_sig_col)) return(safe_n_distinct(df[[qtl_sig_col]]))
  if (nrow(signal_parts) > 0) return(safe_n_distinct(signal_parts$qtl_sig))
  NA_integer_
}

summarize_one <- function(df, output_type, qtl_label) {
  qdf <- df %>% filter(.data$qtl_label == .env$qtl_label)
  gene_col <- find_col(qdf, c("Gene", "gene"))
  signal_col <- find_col(qdf, c("Signal", "signal"))
  snp_col <- find_col(qdf, c("SNP", "snp", "variant", "variant_id", "rsid"))
  gwas_cs_col <- find_col(qdf, c("gwas_cs"))
  qtl_sig_col <- find_col(qdf, c("qtl_sig"))
  rcp_col <- find_col(qdf, c("RCP"))
  lcp_col <- find_col(qdf, c("LCP"))
  grcp_col <- find_col(qdf, c("GRCP"))
  glcp_col <- find_col(qdf, c("GLCP"))
  clpp_col <- find_col(qdf, c("CLPP"))
  n_shared_col <- find_col(qdf, c("n_shared"))
  generic <- generic_numeric_summary(qdf)

  signal_parts <- if (!is.na(signal_col) && nrow(qdf) > 0) {
    extract_signal_parts(qdf[[signal_col]])
  } else {
    tibble(qtl_sig = character(), gwas_cs = character())
  }

  tibble(
    qtl_label = qtl_label,
    output_type = output_type,
    n_rows = nrow(qdf),
    n_studies = if ("study_id" %in% names(qdf)) safe_n_distinct(qdf$study_id) else NA_integer_,
    n_traits = if ("gwas_trait" %in% names(qdf)) safe_n_distinct(qdf$gwas_trait) else NA_integer_,
    n_trait_categories = if ("trait_category" %in% names(qdf)) safe_n_distinct(qdf$trait_category) else NA_integer_,
    n_genes = if (!is.na(gene_col)) safe_n_distinct(qdf[[gene_col]]) else NA_integer_,
    n_gwas_credible_sets = count_gwas_credible_sets(qdf, gwas_cs_col, signal_parts),
    n_qtl_signals = count_qtl_signals(qdf, qtl_sig_col, signal_parts),
    n_variants_or_snps = if (!is.na(snp_col)) safe_n_distinct(qdf[[snp_col]]) else NA_integer_,
    max_RCP = if (!is.na(rcp_col)) safe_max(qdf[[rcp_col]]) else NA_real_,
    mean_RCP = if (!is.na(rcp_col)) safe_mean(qdf[[rcp_col]]) else NA_real_,
    n_RCP_ge_0.1 = if (!is.na(rcp_col)) safe_count_ge(qdf[[rcp_col]], 0.1) else NA_integer_,
    n_RCP_ge_0.5 = if (!is.na(rcp_col)) safe_count_ge(qdf[[rcp_col]], 0.5) else NA_integer_,
    n_RCP_ge_0.9 = if (!is.na(rcp_col)) safe_count_ge(qdf[[rcp_col]], 0.9) else NA_integer_,
    max_LCP = if (!is.na(lcp_col)) safe_max(qdf[[lcp_col]]) else NA_real_,
    mean_LCP = if (!is.na(lcp_col)) safe_mean(qdf[[lcp_col]]) else NA_real_,
    max_GRCP = if (!is.na(grcp_col)) safe_max(qdf[[grcp_col]]) else NA_real_,
    mean_GRCP = if (!is.na(grcp_col)) safe_mean(qdf[[grcp_col]]) else NA_real_,
    n_GRCP_ge_0.5 = if (!is.na(grcp_col)) safe_count_ge(qdf[[grcp_col]], 0.5) else NA_integer_,
    max_GLCP = if (!is.na(glcp_col)) safe_max(qdf[[glcp_col]]) else NA_real_,
    mean_GLCP = if (!is.na(glcp_col)) safe_mean(qdf[[glcp_col]]) else NA_real_,
    n_GLCP_ge_0.5 = if (!is.na(glcp_col)) safe_count_ge(qdf[[glcp_col]], 0.5) else NA_integer_,
    max_CLPP = if (!is.na(clpp_col)) safe_max(qdf[[clpp_col]]) else NA_real_,
    mean_CLPP = if (!is.na(clpp_col)) safe_mean(qdf[[clpp_col]]) else NA_real_,
    n_CLPP_ge_0.01 = if (!is.na(clpp_col)) safe_count_ge(qdf[[clpp_col]], 0.01) else NA_integer_,
    n_CLPP_ge_0.05 = if (!is.na(clpp_col)) safe_count_ge(qdf[[clpp_col]], 0.05) else NA_integer_,
    n_CLPP_ge_0.1 = if (!is.na(clpp_col)) safe_count_ge(qdf[[clpp_col]], 0.1) else NA_integer_,
    mean_n_shared = if (!is.na(n_shared_col)) safe_mean(qdf[[n_shared_col]]) else NA_real_,
    max_n_shared = if (!is.na(n_shared_col)) safe_max(qdf[[n_shared_col]]) else NA_real_,
    numeric_column_count = generic$numeric_column_count,
    numeric_columns = generic$numeric_columns,
    max_numeric_value = generic$max_numeric_value,
    mean_numeric_value = generic$mean_numeric_value
  )
}

summarize_table <- function(path, output_type, qtl_labels) {
  df <- if (startsWith(output_type, "fastenloc_")) {
    read_fastenloc_combined(path, output_type)
  } else {
    read_any(path)
  }
  if (!"qtl_label" %in% names(df)) {
    stop(output_type, " input is missing required qtl_label column: ", path, call. = FALSE)
  }
  map_dfr(qtl_labels, ~ summarize_one(df, output_type, .x))
}

main <- function() {
  a <- parse_args(build_option_parser())
  required <- c("gene", "enrich", "mi", "sig", "snp", "clpp", "qtl_labels")
  missing <- required[vapply(required, function(x) is.null(a[[x]]), logical(1))]
  if (length(missing) > 0) {
    print_help(build_option_parser())
    stop("Missing required option(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  qtl_labels <- readLines(a$qtl_labels, warn = FALSE)
  qtl_labels <- qtl_labels[qtl_labels != ""]
  if (length(qtl_labels) == 0) stop("--qtl_labels must contain at least one label", call. = FALSE)

  summary <- bind_rows(
    summarize_table(a$gene, "fastenloc_gene", qtl_labels),
    summarize_table(a$enrich, "fastenloc_enrich", qtl_labels),
    summarize_table(a$mi, "fastenloc_mi", qtl_labels),
    summarize_table(a$sig, "fastenloc_sig", qtl_labels),
    summarize_table(a$snp, "fastenloc_snp", qtl_labels),
    summarize_table(a$clpp, "clpp", qtl_labels)
  ) %>%
    arrange(factor(qtl_label, levels = qtl_labels), output_type)

  write_tsv(summary, a$out)
  dir.create(a$per_qtl_dir, showWarnings = FALSE, recursive = TRUE)
  for (label in qtl_labels) {
    write_tsv(
      summary %>% filter(qtl_label == label),
      file.path(a$per_qtl_dir, paste0(label, ".raw_coloc_summary.tsv"))
    )
  }
  message("Wrote raw colocalization summary -> ", a$out)
  message("Wrote per-QTL summaries -> ", a$per_qtl_dir)
}

if (identical(environment(), globalenv()) &&
    length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
}
