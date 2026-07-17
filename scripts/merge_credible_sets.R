#!/usr/bin/env Rscript
# ============================================================================
# merge_credible_sets.R
#
# Merge fine-mapped GWAS credible sets ACROSS studies of the same trait when
# their variant-membership Jaccard index >= a threshold, so that a locus
# fine-mapped independently in two studies is not double-counted.
#
# INPUT: a manifest (tsv/csv, +.gz ok) describing the GWAS analysis units,
# needing at least these columns:
#   study_id   unique id per row
#   trait      grouping column (override with --group_col)
#   gwas_path  path to that study's FastENLOC-format GWAS file
#
# OUTPUT (--out): one row per original credible set with its consensus-locus
# assignment:
#   trait, study_id, gwas_cs, chr, start, cs_size, consensus_locus_id,
#   n_studies_in_locus, n_cs_in_locus, is_merged
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

build_option_parser <- function() {
  OptionParser(
    usage = "Rscript %prog --manifest manifest.tsv --out consensus_loci.tsv.gz [options]",
    description = paste(
      "Merge GWAS credible sets across studies of the same trait when their",
      "variant-membership Jaccard index >= threshold."
    ),
    option_list = list(
      make_option("--manifest", type = "character", default = NULL,
                  help = "Manifest tsv/csv with study_id, trait, gwas_path [required]"),
      make_option("--out", type = "character", default = NULL,
                  help = "Per-credible-set consensus map; gzipped if .gz [required]"),
      make_option("--summary_out", type = "character", default = NULL,
                  help = "Optional per-trait merge-effect summary [optional]"),
      make_option("--jaccard", type = "double", default = 0.90,
                  help = "Jaccard threshold to merge two credible sets [default %default]"),
      make_option("--group_col", type = "character", default = "trait",
                  help = "Manifest column to group studies by [default %default]")
    )
  )
}

read_any <- function(path, ...) {
  args <- list(showProgress = FALSE, ...)
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                        args$input <- path
  as_tibble(do.call(data.table::fread, args))
}

read_membership <- function(gwas_path, study_id, group_val) {
  read_any(gwas_path, header = FALSE, select = c(3, 6),
           col.names = c("variant_id", "annot"), colClasses = "character") %>%
    mutate(gwas_cs = str_extract(annot, "^[^=]+")) %>%
    filter(!is.na(gwas_cs), gwas_cs != "") %>%
    transmute(group = group_val, study_id = study_id,
              gwas_cs, variant_id) %>%
    distinct(group, study_id, gwas_cs, variant_id)
}

validate_membership <- function(memb, manifest = NULL) {
  placeholder <- grepl(
    "(_L(None|NA|NaN)$|chr([0-9]+|X|Y|M|MT)\\.0\\.0_)",
    memb$gwas_cs,
    ignore.case = TRUE
  )
  if (any(placeholder)) {
    bad <- unique(memb$gwas_cs[placeholder])
    stop(
      "Placeholder/malformed credible-set ID(s) detected: ",
      paste(head(bad, 10), collapse = ", "),
      if (length(bad) > 10) paste0(" ... and ", length(bad) - 10, " more") else "",
      ". Rebuild the GWAS fastENLOC input with distinct study-locus IDs.",
      call. = FALSE
    )
  }
  if (!is.null(manifest) && "n_credible_sets" %in% names(manifest)) {
    expected <- manifest %>%
      transmute(study_id, expected = as.integer(n_credible_sets))
    observed <- memb %>%
      distinct(study_id, gwas_cs) %>%
      count(study_id, name = "observed")
    audit <- expected %>%
      left_join(observed, by = "study_id") %>%
      mutate(observed = replace_na(observed, 0L))
    mismatch <- audit %>% filter(expected != observed)
    if (nrow(mismatch) > 0) {
      detail <- paste0(
        mismatch$study_id, " expected ", mismatch$expected,
        " but observed ", mismatch$observed
      )
      stop(
        "Manifest credible-set count mismatch: ",
        paste(detail, collapse = "; "),
        call. = FALSE
      )
    }
  }
  invisible(memb)
}

make_uf <- function(n) {
  parent <- seq_len(n)
  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  list(
    union = function(a, b) {
      ra <- find(a)
      rb <- find(b)
      if (ra != rb) parent[rb] <<- ra
    },
    roots = function() vapply(seq_len(n), find, integer(1))
  )
}

chrom_order <- function(chr) {
  case_when(
    chr %in% as.character(1:22) ~ as.integer(chr),
    chr == "X" ~ 23L,
    chr == "Y" ~ 24L,
    chr %in% c("M", "MT") ~ 25L,
    TRUE ~ 99L
  )
}

merge_one_trait <- function(memb, jaccard = 0.90) {
  cs <- memb %>%
    distinct(study_id, gwas_cs) %>%
    mutate(cs_idx = row_number())
  sizes <- memb %>% count(study_id, gwas_cs, name = "cs_size")
  cs <- cs %>% left_join(sizes, by = c("study_id", "gwas_cs"))

  dt <- as.data.table(
    memb %>%
      left_join(cs %>% select(study_id, gwas_cs, cs_idx),
                by = c("study_id", "gwas_cs"))
  )[, .(cs_idx, variant_id, study_id)]
  edges <- data.table(a = integer(0), b = integer(0))
  if (nrow(dt) > 0) {
    j <- merge(dt, dt, by = "variant_id", allow.cartesian = TRUE)
    j <- j[cs_idx.x < cs_idx.y & study_id.x != study_id.y]
    if (nrow(j) > 0) {
      shared <- j[, .(n_shared = .N), by = .(a = cs_idx.x, b = cs_idx.y)]
      sz <- setNames(cs$cs_size, cs$cs_idx)
      shared[, jac := n_shared / (sz[as.character(a)] + sz[as.character(b)] - n_shared)]
      edges <- shared[jac >= jaccard, .(a, b)]
    }
  }

  uf <- make_uf(nrow(cs))
  if (nrow(edges) > 0) for (k in seq_len(nrow(edges))) uf$union(edges$a[k], edges$b[k])
  cs$root <- uf$roots()

  cs <- cs %>%
    mutate(
      chr = str_match(gwas_cs, "chr([0-9XYMT]+)\\.")[, 2],
      start = suppressWarnings(as.integer(str_match(gwas_cs, "chr[0-9XYMT]+\\.([0-9]+)")[, 2])),
      ord_chr = chrom_order(chr)
    )

  comp <- cs %>%
    group_by(root) %>%
    summarise(ord_chr = min(ord_chr, na.rm = TRUE),
              ord_start = min(start, na.rm = TRUE),
              n_studies_in_locus = n_distinct(study_id),
              n_cs_in_locus = n(),
              .groups = "drop") %>%
    arrange(ord_chr, ord_start) %>%
    mutate(locus_num = row_number())

  cs %>%
    left_join(comp %>% select(root, locus_num, n_studies_in_locus, n_cs_in_locus),
              by = "root") %>%
    mutate(is_merged = n_cs_in_locus > 1) %>%
    select(study_id, gwas_cs, chr, start, cs_size,
           locus_num, n_studies_in_locus, n_cs_in_locus, is_merged)
}

main <- function() {
  a <- parse_args(build_option_parser())
  if (is.null(a$manifest) || is.null(a$out)) {
    print_help(build_option_parser())
    stop("--manifest and --out are required.", call. = FALSE)
  }
  man <- read_any(a$manifest)
  need <- c("study_id", a$group_col, "gwas_path")
  miss <- setdiff(need, names(man))
  if (length(miss) > 0) stop("manifest missing columns: ", paste(miss, collapse = ", "))

  memb_all <- pmap_dfr(
    list(man$gwas_path, man$study_id, man[[a$group_col]]),
    function(p, s, g) {
      message("Reading ", s, " (", g, "): ", basename(p))
      read_membership(p, s, g)
    }
  )
  validate_membership(memb_all, man)

  result <- memb_all %>%
    group_split(group) %>%
    map_dfr(function(g_memb) {
      grp <- g_memb$group[1]
      out <- merge_one_trait(g_memb %>% select(study_id, gwas_cs, variant_id), a$jaccard)
      grp_tag <- str_replace_all(grp, "[^A-Za-z0-9]+", "_")
      out %>%
        mutate(trait = grp,
               consensus_locus_id = sprintf("%s_L%04d", grp_tag, locus_num)) %>%
        select(trait, study_id, gwas_cs, chr, start, cs_size,
               consensus_locus_id, n_studies_in_locus, n_cs_in_locus, is_merged)
    })
  write_tsv(result, a$out)

  summ <- result %>%
    group_by(trait) %>%
    summarise(
      n_studies = n_distinct(study_id),
      n_raw_cs = n(),
      n_consensus_loci = n_distinct(consensus_locus_id),
      n_merged_away = n_raw_cs - n_consensus_loci,
      .groups = "drop"
    ) %>%
    arrange(desc(n_merged_away))
  message("\n--- merge summary (Jaccard >= ", a$jaccard, ") ---")
  summ %>%
    mutate(line = sprintf("  %-32s %d studies  %4d CS -> %4d loci  (%d merged)",
                          trait, n_studies, n_raw_cs, n_consensus_loci, n_merged_away)) %>%
    pull(line) %>%
    walk(message)
  if (!is.null(a$summary_out)) {
    write_tsv(summ, a$summary_out)
    message("Wrote per-trait summary -> ", a$summary_out)
  }
  message("Wrote ", nrow(result), " credible sets (",
          n_distinct(result$consensus_locus_id), " consensus loci) -> ", a$out)
}

if (identical(environment(), globalenv()) &&
    length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
}
