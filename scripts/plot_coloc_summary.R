#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
  library(ggh4x)
})

GROUP_LABELS <- c(
  autoimmune = "AutoImm",
  inflammatory_biomarker = "Inflam",
  allergic_respiratory = "Allerg",
  dermatological = "Derm",
  cardiometabolic = "CardioM",
  lipid_metabolic = "Lipid",
  lipid = "Lipid",
  metabolic = "Metabolic",
  clinical_chemistry = "ClinChem",
  hematological = "Heme",
  anthropometric = "Anthro",
  renal = "Renal",
  renal_urological = "Renal/Urol",
  gastrointestinal = "GI",
  reproductive_endocrine = "ReproEnd",
  musculoskeletal = "MSK",
  cancer = "Cancer",
  neurodegenerative = "Neurodeg",
  neurological = "Neurol",
  neuropsychiatric = "NeuroPsy",
  ophthalmological = "Ophth"
)

GROUP_COLORS <- c(
  autoimmune = "#D62828",
  inflammatory_biomarker = "#F77F00",
  allergic_respiratory = "#2A9D8F",
  dermatological = "#E76F51",
  cardiometabolic = "#277DA1",
  lipid_metabolic = "#6A994E",
  lipid = "#6A994E",
  metabolic = "#80B918",
  clinical_chemistry = "#F4A261",
  hematological = "#A23B72",
  anthropometric = "#E9C46A",
  renal = "#457B9D",
  renal_urological = "#8D99AE",
  gastrointestinal = "#BC6C25",
  reproductive_endocrine = "#E5989B",
  musculoskeletal = "#7F5539",
  cancer = "#9B5DE5",
  neurodegenerative = "#4361EE",
  neurological = "#3A0CA3",
  neuropsychiatric = "#B5179E",
  ophthalmological = "#48CAE4"
)

build_option_parser <- function() {
  OptionParser(
    usage = paste(
      "Rscript %prog --gene harmonized_coloc.gene.tsv.gz",
      "--coloc_rate coloc_rate_by_trait.tsv [options]"
    ),
    description = paste(
      "Plot distinct colocalizing-gene counts and stringent consensus-locus",
      "colocalization rates by trait category."
    ),
    option_list = list(
      make_option("--gene", type = "character", default = NULL,
                  help = "Aggregated harmonized gene table (.tsv/.gz) [required]"),
      make_option("--coloc_rate", type = "character", default = NULL,
                  help = "Trait/layer coloc-rate table from summarize_coloc.R [required]"),
      make_option("--plot_data_out", type = "character", default = "coloc_summary_plot_data.tsv",
                  help = "Joined trait-level plotting data [default %default]"),
      make_option("--png_out", type = "character", default = "coloc_summary.png",
                  help = "PNG figure [default %default]"),
      make_option("--pdf_out", type = "character", default = "coloc_summary.pdf",
                  help = "PDF figure [default %default]"),
      make_option("--width", type = "double", default = 10,
                  help = "Figure width in inches [default %default]"),
      make_option("--height", type = "double", default = NA_real_,
                  help = "Figure height in inches; automatically sized when omitted"),
      make_option("--dpi", type = "integer", default = 300,
                  help = "PNG resolution [default %default]")
    )
  )
}

read_any <- function(path, ...) {
  args <- list(showProgress = FALSE, ...)
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                         args$input <- path
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
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

complete_category_colors <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  missing <- setdiff(categories, names(GROUP_COLORS))
  fallback <- if (length(missing) > 0) {
    setNames(scales::hue_pal(l = 55, c = 90)(length(missing)), missing)
  } else {
    character()
  }
  c(GROUP_COLORS, fallback)
}

blend_with_white <- function(hex, amount = 0.65) {
  x <- grDevices::col2rgb(hex) / 255
  y <- x * (1 - amount) + amount
  grDevices::rgb(y[1, ], y[2, ], y[3, ])
}

make_fill <- function(category, intensity, colors) {
  base <- unname(colors[as.character(category)])
  light <- blend_with_white(base, 0.72)
  vapply(seq_along(base), function(i) {
    palette <- grDevices::colorRampPalette(c(light[[i]], base[[i]]))(100)
    palette[pmax(1, pmin(100, round(intensity[[i]] * 99) + 1))]
  }, character(1))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  scale_y_discrete(..., labels = function(x) gsub(reg, "", x))
}

build_plot_data <- function(gene, coloc_rate) {
  require_columns(gene, c("gene", "trait", "trait_category", "any_coloc"), "--gene")
  require_columns(coloc_rate, c("trait", "layer", "pct_stringent_union"), "--coloc_rate")

  gene$any_coloc <- as_logical_flag(gene$any_coloc)
  coloc_counts <- gene %>%
    filter(any_coloc, !is.na(gene), gene != "", !is.na(trait_category), trait_category != "") %>%
    distinct(gene, trait, trait_category) %>%
    count(trait, trait_category, name = "count")

  pct_coloc <- coloc_rate %>%
    filter(layer == "union") %>%
    select(trait, pct_stringent_union)
  duplicated_traits <- pct_coloc %>% count(trait) %>% filter(n > 1) %>% pull(trait)
  if (length(duplicated_traits) > 0) {
    stop("--coloc_rate has multiple union rows for trait(s): ",
         paste(duplicated_traits, collapse = ", "), call. = FALSE)
  }

  combined <- coloc_counts %>% left_join(pct_coloc, by = "trait")
  missing_rate <- combined %>% filter(is.na(pct_stringent_union)) %>% pull(trait)
  if (length(missing_rate) > 0) {
    message("Dropped (no union stringent coloc rate): ",
            paste(unique(missing_rate), collapse = ", "))
  }
  combined %>% filter(!is.na(pct_stringent_union))
}

build_empty_plot <- function() {
  ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label = "No traits have both colocalizing genes and a union stringent coloc rate",
      color = "grey35", size = 4
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void()
}

build_plot <- function(combined) {
  colors <- complete_category_colors(combined$trait_category)
  plotdat <- combined %>%
    group_by(trait_category) %>%
    mutate(
      intensity = if_else(max(count) == min(count), 0.75,
                          (count - min(count)) / (max(count) - min(count))),
      intensity = 0.25 + 0.75 * intensity,
      fill_col = make_fill(trait_category, intensity, colors)
    ) %>%
    ungroup() %>%
    mutate(
      trait_category = forcats::fct_drop(
        forcats::fct_reorder(trait_category, count, .fun = sum)
      ),
      trait = reorder_within(trait, count, trait_category),
      base_col = unname(colors[as.character(trait_category)])
    )

  long <- bind_rows(
    plotdat %>% transmute(trait_category, trait, base_col, fill_col,
                          metric = "a_count", value = count),
    plotdat %>% transmute(trait_category, trait, base_col, fill_col,
                          metric = "b_rate", value = pct_stringent_union)
  )
  xline <- -max(plotdat$count) * 0.03
  seg <- plotdat %>% count(trait_category, name = "nbar") %>% mutate(metric = "a_count")
  vref <- tibble(metric = "b_rate", xi = median(plotdat$pct_stringent_union, na.rm = TRUE))
  metric_labels <- c(a_count = "colocalizing genes", b_rate = "CS coloc rate")
  category_labeller <- c(GROUP_LABELS, setNames(
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS)),
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS))
  ))

  ggplot(long, aes(y = trait)) +
    geom_col(data = ~ subset(.x, metric == "a_count"),
             aes(x = value, fill = fill_col), width = 0.85, show.legend = FALSE) +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(y = 0.5, yend = nbar + 0.5, x = xline, xend = xline),
                 color = "grey60", linewidth = 0.5) +
    geom_vline(data = vref, aes(xintercept = xi),
               linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_segment(data = ~ subset(.x, metric == "b_rate"),
                 aes(x = 0, xend = value, yend = trait, color = base_col), linewidth = 0.4) +
    geom_point(data = ~ subset(.x, metric == "b_rate"),
               aes(x = value, color = base_col), size = 1.6) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_y_reordered() +
    facet_grid(
      trait_category ~ metric,
      scales = "free",
      space = "free_y",
      switch = "y",
      labeller = labeller(
        trait_category = as_labeller(category_labeller, default = label_value),
        metric = as_labeller(metric_labels)
      )
    ) +
    ggh4x::facetted_pos_scales(x = list(
      scale_x_continuous(),
      scale_x_continuous(
        labels = function(x) paste0(x, "%"),
        expand = expansion(mult = c(0, 0.08))
      )
    )) +
    ggh4x::force_panelsizes(cols = c(3, 1.3)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.spacing.x = grid::unit(8, "pt"),
      panel.spacing.y = grid::unit(2, "pt"),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
      strip.text.x = element_text(size = 9)
    )
}

main <- function() {
  a <- parse_args(build_option_parser())
  if (is.null(a$gene) || is.null(a$coloc_rate)) {
    print_help(build_option_parser())
    stop("--gene and --coloc_rate are required.", call. = FALSE)
  }
  if (a$width <= 0 || (!is.na(a$height) && a$height <= 0) || a$dpi <= 0) {
    stop("--width, --height, and --dpi must be positive.", call. = FALSE)
  }

  combined <- build_plot_data(read_any(a$gene), read_any(a$coloc_rate))
  write_tsv(combined %>% arrange(trait_category, desc(count), trait), a$plot_data_out)

  height <- if (is.na(a$height)) max(6, 2 + 0.18 * nrow(combined)) else a$height
  p <- if (nrow(combined) == 0) build_empty_plot() else build_plot(combined)
  ggsave(a$png_out, p, width = a$width, height = height, units = "in", dpi = a$dpi, bg = "white")
  ggsave(a$pdf_out, p, width = a$width, height = height, units = "in", bg = "white")
  message("Wrote plot data (", nrow(combined), " traits) -> ", a$plot_data_out)
  message("Wrote summary figures -> ", a$png_out, ", ", a$pdf_out)
}

if (identical(environment(), globalenv()) &&
    length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
}
