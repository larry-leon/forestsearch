# §3.2 gate: the two display-only ```r listings of vignettes/articles/methodology.Rmd
# against formals(forestsearch::forestsearch) (installed build = HEAD's R/ tree).
# Asserts: every named argument (other than the parallel_args elements workers /
# show_message) is a current formal; every "# default <text>" comment parses to the
# deparsed formal default (or, for "# default \"<x>\" (match.arg)", to the first element of
# the formal's c() choice vector); none of vi.grf.min, max_n_confounders, max.minutes appears.
# usage: Rscript methodology_listings_check.R <methodology.Rmd> [--report]
args   <- commandArgs(TRUE)
report <- "--report" %in% args
x  <- readLines(args[1], encoding = "UTF-8")
fm <- formals(forestsearch::forestsearch)
open  <- which(x == "```r")
close <- which(x == "```")
stopifnot(length(open) == 2L)
blocks <- lapply(open, function(o) seq(o + 1L, min(close[close > o]) - 1L))
pat  <- "^\\s*([A-Za-z_.][A-Za-z0-9_.]*)\\s*=\\s*(.*?)\\s*#\\s*default\\s+(.*?)\\s*(<-- FOCAL)?\\s*$"
skip <- c("workers", "show_message")
bad <- character(0)
for (b in seq_along(blocks)) {
  for (i in blocks[[b]]) {
    ln <- x[i]
    for (w in c("vi.grf.min", "max_n_confounders", "max.minutes"))
      if (grepl(w, ln, fixed = TRUE)) bad <- c(bad, sprintf("block %d line %d: forbidden %s", b, i, w))
    m <- regmatches(ln, regexec(pat, ln, perl = TRUE))[[1]]
    if (!length(m)) next
    nm <- m[2]; txt <- m[4]
    if (nm %in% skip) next
    if (!nm %in% names(fm)) { bad <- c(bad, sprintf("block %d line %d: %s not a formal", b, i, nm)); next }
    want <- paste(deparse(fm[[nm]]), collapse = "")
    # '# default "<first choice>" (match.arg)': the formal is a c() of character
    # choices and the comment names the first one, the effective default
    if (grepl("\\s*\\(match\\.arg\\)$", txt)) {
      txt <- sub("\\s*\\(match\\.arg\\)$", "", txt)
      d <- fm[[nm]]
      ok <- is.call(d) && identical(d[[1]], as.name("c")) && length(d) > 1L && is.character(d[[2]])
      want <- if (ok) deparse(d[[2]]) else paste0("<not a match.arg choice vector: ", want, ">")
    }
    got  <- tryCatch(paste(deparse(parse(text = txt, keep.source = FALSE)[[1]]), collapse = ""),
                     error = function(e) NA_character_)
    if (!identical(got, want))
      bad <- c(bad, sprintf("block %d line %d: %s  comment '# default %s'  formal default %s", b, i, nm, txt, want))
  }
}
cat(sprintf("checked %d blocks (lines %s)\n", length(blocks),
            paste(sapply(blocks, function(v) paste0(min(v), "-", max(v))), collapse = ", ")))
if (length(bad)) { cat("MISMATCHES:\n", paste0("  ", bad, "\n"), sep = ""); if (!report) quit(status = 1) } else cat("ALL CHECKS PASS\n")
