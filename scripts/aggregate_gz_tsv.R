#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

parser <- OptionParser(option_list = list(
  make_option("--files", type = "character",
              help = "Text file containing one input TSV or TSV.GZ path per line"),
  make_option("--out", type = "character",
              help = "Output .tsv.gz path")
))

main <- function() {
  args <- parse_args(parser)
  if (is.null(args$files) || is.null(args$out)) {
    print_help(parser)
    stop("--files and --out are required", call. = FALSE)
  }

  files <- readLines(args$files, warn = FALSE)
  files <- files[files != ""]
  if (length(files) == 0) stop("No input files were provided", call. = FALSE)

  out <- gzfile(args$out, "wt")
  on.exit(if (!is.null(out)) close(out), add = TRUE)
  expected_header <- NULL

  for (path in files) {
    con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
    header <- readLines(con, n = 1, warn = FALSE)
    if (length(header) == 0) {
      close(con)
      next
    }
    if (is.null(expected_header)) {
      expected_header <- header
      writeLines(header, out)
    } else if (!identical(header, expected_header)) {
      close(con)
      stop("Header mismatch while aggregating ", path, call. = FALSE)
    }

    repeat {
      chunk <- readLines(con, n = 100000, warn = FALSE)
      if (length(chunk) == 0) break
      writeLines(chunk, out)
    }
    close(con)
  }

  if (is.null(expected_header)) {
    stop("No input rows found while aggregating gzipped TSV files", call. = FALSE)
  }
  close(out)
  out <- NULL
}

main()
