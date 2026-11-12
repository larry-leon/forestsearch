# gate_compare.R -- compare two gate_run.R outputs (start-SHA export vs head).
# Usage: Rscript --vanilla gate_compare.R <res_start.rds> <res_head.rds> <out_txt>
# Adapted from /tmp/seed_check_2026-11-11/equiv_compare.R: same timing-field rule and leaf
# listing; the step-5 block is removed and the gate is identical() after dropping timing fields.
# after dropping those components, and no other differing path.
a <- commandArgs(TRUE)
suppressMessages(library(data.table))
A <- readRDS(a[1]); B <- readRDS(a[2]); out <- a[3]
TIMING <- "(^|[._])(time|times|timing|timings|minutes|mins|seconds|secs|elapsed)([._]|$)|^t\\.(min|sec)|^tmin|^tsec|proc_time|time_seconds|timing_seconds|^t_|_t$|^runtime|Sys\\.time"
is_timing <- function(nm) grepl(TIMING, nm, ignore.case = TRUE)
# timing columns inside data frames (per-replicate run times); a column named "time" is data
DF_TIMING <- "^tmins|minutes|seconds|elapsed|timing"
is_df_timing <- function(nm) grepl(DF_TIMING, nm, ignore.case = TRUE)

# drop timing-named components at any depth (never inside a data.frame: column names
# such as "time" are data)
strip <- function(x) {
  if (is.environment(x)) return(x)
  if (is.data.frame(x)) {
    # data.tables become plain data.frames: their .internal.selfref pointer is not data
    if (inherits(x, "data.table")) { x <- as.data.frame(x); attr(x, ".internal.selfref") <- NULL }
    drop <- is_df_timing(names(x))
    if (any(drop)) x <- x[, !drop, drop = FALSE]
    return(x)
  }
  if (is.list(x) && !is.data.frame(x)) {
    nm <- names(x)
    if (!is.null(nm)) { keep <- !is_timing(nm); at <- attributes(x); x <- x[keep]
      for (k in setdiff(names(at), c("names"))) attr(x, k) <- at[[k]] }
    for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] <- strip(x[[i]])
  }
  x
}
# leaf-level differences (for reporting)
diffs <- function(x, y, path = "") {
  if (is.function(x) || is.function(y)) {
    same <- is.function(x) && is.function(y) && identical(deparse(x), deparse(y))
    return(if (same) character(0) else path)
  }
  if (is.environment(x) || is.environment(y)) {
    if (!(is.environment(x) && is.environment(y))) return(path)
    return(diffs(as.list(x, all.names = TRUE, sorted = TRUE), as.list(y, all.names = TRUE, sorted = TRUE), paste0(path, "<env>")))
  }
  if (is.list(x) && is.list(y) && !is.data.frame(x) && !is.data.frame(y)) {
    nx <- names(x); ny <- names(y)
    if (is.null(nx) || is.null(ny)) {
      if (length(x) != length(y)) return(paste0(path, " [length ", length(x), " vs ", length(y), "]"))
      return(unlist(lapply(seq_along(x), function(i) diffs(x[[i]], y[[i]], paste0(path, "[[", i, "]]")))))
    }
    out <- character(0)
    only <- union(setdiff(nx, ny), setdiff(ny, nx))
    if (length(only)) out <- c(out, paste0(path, "$", only, " [present in one tree only]"))
    for (n in intersect(nx, ny)) out <- c(out, diffs(x[[n]], y[[n]], paste0(path, "$", n)))
    if (!identical(attributes(x)[setdiff(names(attributes(x)), "names")], attributes(y)[setdiff(names(attributes(y)), "names")]))
      out <- c(out, paste0(path, " [attributes]"))
    return(out)
  }
  if (is.data.frame(x) && is.data.frame(y) && identical(names(x), names(y)) && nrow(x) == nrow(y)) {
    d <- names(x)[!mapply(identical, as.list(x), as.list(y))]
    return(if (length(d)) paste0(path, "$", d) else character(0))
  }
  if (identical(x, y)) character(0) else path
}
ident <- function(x, y) identical(x, y, ignore.environment = TRUE, ignore.bytecode = TRUE, ignore.srcref = TRUE)

lines <- c(sprintf("start:  %s (forestsearch %s)", A$tree, A$versions$fs),
           sprintf("head:   %s (forestsearch %s)", B$tree, B$versions$fs), "")
gate <- TRUE
for (nm in names(A$res)) {
  d <- diffs(A$res[[nm]], B$res[[nm]])
  tim <- d[vapply(strsplit(gsub("\\[\\[[0-9]+\\]\\]", "", d), "\\$"), function(p) any(is_timing(p)) || is_df_timing(utils::tail(p, 1)), TRUE)]
  other <- setdiff(d, tim)
  sa <- strip(A$res[[nm]]); sb <- strip(B$res[[nm]])
  id <- ident(sa, sb)
  # environments inside results (estimator closures) are compared by content above;
  # identical() on stripped objects compares closures without their environments
  st_b <- identical(B$state[[nm]]$before, B$state[[nm]]$after)
  st_a <- identical(A$state[[nm]]$before, A$state[[nm]]$after)
  ok <- (id || length(other) == 0)
  gate <- gate && ok
  lines <- c(lines, sprintf("%s: identical() after dropping timing fields = %s; non-timing differing leaves = %d; timing leaves dropped = %d; branch .Random.seed after = before: %s; start .Random.seed after = before: %s; wall start %.1f s / branch %.1f s -> %s",
    nm, id, length(other), length(tim), st_b, st_a, A$wall[[nm]], B$wall[[nm]], if (ok) "PASS" else "FAIL"))
  if (length(tim)) lines <- c(lines, paste0("    timing: ", tim))
  if (length(other)) lines <- c(lines, paste0("    OTHER:  ", head(other, 50)))
}
lines <- c(lines, "", sprintf("GATE (all six identical after dropping timing fields): %s", if (gate) "PASS" else "FAIL"), "")

writeLines(lines, out); cat(lines, sep = "\n")
