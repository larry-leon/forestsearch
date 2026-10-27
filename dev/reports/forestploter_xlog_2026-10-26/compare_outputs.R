# usage: Rscript compare_outputs.R <before.html> <after.html> <outdir>
a <- commandArgs(TRUE)
source("dev/reports/rerender_2026-10-26/extract_lib.R")
b <- extract_outputs(a[1]); n <- extract_outputs(a[2])
writeLines(b, file.path(a[3], "outputs_before.txt")); writeLines(n, file.path(a[3], "outputs_after.txt"))
cat("lines before:", length(b), " after:", length(n), " identical:", identical(b, n), "\n")
