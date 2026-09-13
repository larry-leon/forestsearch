# Fidelity check for the curated sections of current_status.md
# (TASK_status_curated_restore_2026-09-13 §6).
#
#   1. Each marked block of status_curated.md appears in current_status.md as one
#      contiguous run of lines, byte-identical (identical() on the raw lines).
#   2. Against a reference commit's current_status.md: the curated source ranges
#      are byte-identical to status_curated.md's blocks, apart from lines ADDED to
#      status_curated.md since the reference, which are listed.
#
# usage (from scripts_dinamr/): Rscript check_status_curated.R <ref-commit> <ranges>
#   e.g.  Rscript check_status_curated.R 0ab5d1c5 5-6,10-90,132-147,149-200
# Exit 0 = PASS, 1 = FAIL.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript check_status_curated.R <ref-commit> <ranges>")
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
qmd  <- normalizePath(file.path(here, ".."))
repo <- system2("git", c("-C", qmd, "rev-parse", "--show-toplevel"), stdout = TRUE)
rq   <- sub(paste0("^", repo, "/"), "", qmd)
rd <- function(p) readLines(p, encoding = "UTF-8", warn = FALSE)
cur <- rd(file.path(qmd, "status_curated.md")); cs <- rd(file.path(qmd, "current_status.md"))
marks <- c("preamble", "before-inventory", "after-inventory-table", "after-inventory")
pos <- match(sprintf("<!-- curated:%s -->", marks), cur)
stopifnot(!anyNA(pos), !is.unsorted(pos))
blocks <- lapply(seq_along(pos), function(i) cur[(pos[i] + 1L):(if (i < length(pos)) pos[i + 1L] - 1L else length(cur))])
names(blocks) <- marks
ok <- TRUE
find_run <- function(hay, needle) {
  for (i in seq_len(length(hay) - length(needle) + 1L))
    if (identical(hay[i:(i + length(needle) - 1L)], needle)) return(i)
  NA_integer_ }

cat("== 1. status_curated.md blocks inside current_status.md\n")
for (k in marks) {
  i <- find_run(cs, blocks[[k]])
  cat(sprintf("  %-22s %3d lines, %6d bytes : %s\n", k, length(blocks[[k]]), sum(nchar(blocks[[k]], "bytes") + 1L),
      if (is.na(i)) { ok <- FALSE; "NOT FOUND" } else sprintf("byte-identical at current_status.md lines %d-%d", i, i + length(blocks[[k]]) - 1L)))
}

cat(sprintf("== 2. against %s:%s/current_status.md, ranges %s\n", args[1], rq, args[2]))
ref <- system2("git", c("-C", repo, "show", sprintf("%s:%s/current_status.md", args[1], rq)), stdout = TRUE)
Encoding(ref) <- "UTF-8"
rg <- lapply(strsplit(strsplit(args[2], ",")[[1]], "-"), as.integer)
if (length(rg) != length(marks)) stop("need one range per block (", length(marks), ")")
for (j in seq_along(marks)) {
  r <- ref[rg[[j]][1]:rg[[j]][2]]; b <- blocks[[j]]
  if (identical(r, b)) {
    cat(sprintf("  %-22s ref lines %d-%d: byte-identical (%d lines)\n", marks[j], rg[[j]][1], rg[[j]][2], length(r)))
  } else {
    # the only permitted difference: whole lines inserted into status_curated.md
    added <- integer(0); ri <- 1L
    for (bi in seq_along(b)) { if (ri <= length(r) && identical(b[bi], r[ri])) ri <- ri + 1L else added <- c(added, bi) }
    if (ri == length(r) + 1L) {
      cat(sprintf("  %-22s ref lines %d-%d: byte-identical apart from %d added line(s):\n", marks[j], rg[[j]][1], rg[[j]][2], length(added)))
      for (a in added) cat("      + ", b[a], "\n", sep = "")
    } else { ok <- FALSE
      cat(sprintf("  %-22s ref lines %d-%d: DIFFERS (reference line %d unmatched)\n", marks[j], rg[[j]][1], rg[[j]][2], rg[[j]][1] + ri - 1L)) }
  }
}
cat(if (ok) "PASS\n" else "FAIL\n")
quit(status = if (ok) 0L else 1L)
