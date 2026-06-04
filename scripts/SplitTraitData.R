#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Input fastENLOC file"
  ),
  make_option(
    c("-n", "--traits-per-chunk"),
    type = "integer",
    default = 25,
    help = "Number of traits per output chunk [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) stop("Must provide --input")
if (opt$`traits-per-chunk` < 1) stop("--traits-per-chunk must be >= 1")

outdir <- "chunks"
manifest_file <- "chunk_manifest.txt"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fastenloc <- fread(
  opt$input,
  header = FALSE,
  sep = "\t",
  col.names = c("chr", "pos", "variant", "ref", "alt", "annotation")
) %>%
  as_tibble() %>%
  mutate(
    trait = str_extract(annotation, "^[^;]+")
  )

trait_map <- fastenloc %>%
  distinct(trait) %>%
  arrange(trait) %>%
  mutate(
    chunk_id = ceiling(row_number() / opt$`traits-per-chunk`),
    chunk_file = file.path(
      outdir,
      sprintf("fastenloc_chunk_%03d.txt", chunk_id)
    )
  )

chunked <- fastenloc %>%
  left_join(trait_map, by = "trait")

chunk_files <- trait_map %>%
  distinct(chunk_id, chunk_file) %>%
  arrange(chunk_id) %>%
  pull(chunk_file)

walk(chunk_files, function(outfile) {
  chunked %>%
    filter(chunk_file == outfile) %>%
    select(chr, pos, variant, ref, alt, annotation) %>%
    fwrite(
      file = outfile,
      sep = "\t",
      col.names = FALSE,
      quote = FALSE
    )
})

write_lines(chunk_files, manifest_file)

message("Wrote ", length(chunk_files), " chunk files")
message("Manifest written to: ", manifest_file)
