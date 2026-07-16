#!/usr/bin/env Rscript

source("scripts/plot_coloc_summary.R")

gene <- tibble::tribble(
  ~study, ~gene, ~trait, ~trait_category, ~layer, ~any_coloc,
  "s1", "GENE1", "Trait A", "autoimmune", "eQTL", TRUE,
  "s1", "GENE1", "Trait A", "autoimmune", "sQTL", TRUE,
  "s2", "GENE2", "Trait A", "autoimmune", "eQTL", TRUE,
  "s2", "GENE3", "Trait A", "autoimmune", "eQTL", FALSE,
  "s3", "GENE4", "Trait B", "new_category", "eQTL", TRUE
)

rate <- tibble::tribble(
  ~trait, ~layer, ~pct_stringent_union, ~pct_lenient_union,
  "Trait A", "eQTL", 10, 15,
  "Trait A", "union", 20, 30,
  "Trait B", "union", 40, 50
)

combined <- build_plot_data(gene, rate)
stopifnot(nrow(combined) == 2)
stopifnot(combined$count[combined$trait == "Trait A"] == 2)
stopifnot(combined$pct_coloc_union[combined$trait == "Trait A"] == 20)
stopifnot(combined$pct_stringent_union[combined$trait == "Trait A"] == 20)
stopifnot(combined$pct_lenient_union[combined$trait == "Trait A"] == 30)
lenient <- build_plot_data(gene, rate, rate_type = "lenient")
stopifnot(lenient$pct_coloc_union[lenient$trait == "Trait A"] == 30)
stopifnot(all(c("pct_stringent_union", "pct_lenient_union") %in% names(lenient)))
stopifnot(all(lenient$rate_type == "lenient"))
stopifnot(!is.na(complete_category_colors(combined$trait_category)[["new_category"]]))
stopifnot(identical(as_logical_flag(c("TRUE", "0", "yes")), c(TRUE, FALSE, TRUE)))
stopifnot(identical(as_logical_flag(c(1, 0, NA)), c(TRUE, FALSE, FALSE)))
stopifnot(nrow(build_plot_data(dplyr::filter(gene, !any_coloc), rate)) == 0)
stopifnot(identical(build_plot_data(gene, rate, min_coloc_genes = 2L)$trait, "Trait A"))
stopifnot(nrow(build_plot_data(gene, rate, min_coloc_genes = 3L)) == 0)
stopifnot(inherits(try(build_plot_data(gene, rate, min_coloc_genes = 0L), silent = TRUE), "try-error"))
stopifnot(inherits(try(build_plot_data(gene, rate, rate_type = "other"), silent = TRUE), "try-error"))

png <- tempfile(fileext = ".png")
multi_trait_category <- dplyr::bind_rows(
  combined,
  dplyr::mutate(combined[combined$trait == "Trait A", ], trait = "Trait C", count = 1L)
)
ggplot2::ggsave(png, build_plot(multi_trait_category), width = 8, height = 6, dpi = 72)
stopifnot(file.exists(png), file.info(png)$size > 0)

empty_png <- tempfile(fileext = ".png")
ggplot2::ggsave(empty_png, build_empty_plot(), width = 8, height = 6, dpi = 72)
stopifnot(file.exists(empty_png), file.info(empty_png)$size > 0)

cat("plot_coloc_summary tests passed\n")
