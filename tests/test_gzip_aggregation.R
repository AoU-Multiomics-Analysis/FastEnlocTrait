#!/usr/bin/env Rscript

tmp <- tempfile("gzip-aggregation-")
dir.create(tmp)

first <- file.path(tmp, "first.tsv")
second <- file.path(tmp, "second.tsv.gz")
files <- file.path(tmp, "files.txt")
output <- file.path(tmp, "combined.tsv.gz")

writeLines(c("a\tb", "1\t2"), first)
con <- gzfile(second, "wt")
writeLines(c("a\tb", "3\t4"), con)
close(con)
writeLines(c(first, second), files)

status <- system2(
  "Rscript",
  c("scripts/aggregate_gz_tsv.R", "--files", files, "--out", output)
)
stopifnot(status == 0L)
stopifnot(system2("gzip", c("-t", output)) == 0L)

con <- gzfile(output, "rt")
rows <- readLines(con, warn = FALSE)
close(con)
stopifnot(identical(rows, c("a\tb", "1\t2", "3\t4")))

cat("gzip aggregation regression tests passed\n")
