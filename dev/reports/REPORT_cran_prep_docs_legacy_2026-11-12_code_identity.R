# Code-identity gate for TASK_cran_prep_docs_legacy_2026-11-12 (Part 6 step 2).
#
# Exports R/ at the start SHA and at the head with `git archive` into /tmp,
# parses every file with parse(keep.source = FALSE) (so comments and roxygen
# are dropped), deparses each top-level expression and compares by assigned
# name.  Gate:
#   - names only at the start  = exactly the functions Part 1 removed;
#   - names only at the head   = none;
#   - names in both with a different deparsed body = none;
#   - any other top-level difference is listed, and is allowed only as the
#     removal of names from utils::globalVariables().
#
# Usage (from the repository root):
#   Rscript --vanilla dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12_code_identity.R \
#     <start_sha> <head_sha> <expected_removed_comma_list>

args <- commandArgs(trailingOnly = TRUE)
start_sha <- if (length(args) >= 1) args[[1]] else "d93e6621"
head_sha  <- if (length(args) >= 2) args[[2]] else "HEAD"
expected_removed <- if (length(args) >= 3) {
  strsplit(args[[3]], ",", fixed = TRUE)[[1]]
} else {
  c("create_gbsg_dgm", "simulate_from_gbsg_dgm", ".simulate_from_gbsg_dgm_",
    "get_dgm_with_output")
}

resolve <- function(ref) {
  out <- system2("git", c("rev-parse", shQuote(ref)), stdout = TRUE)
  stopifnot(length(out) == 1L, nchar(out) == 40L)
  out
}
start_full <- resolve(start_sha)
head_full  <- resolve(head_sha)

root <- "/tmp/cran_prep_docs_legacy_2026-11-12/code_identity"
export_r <- function(sha, label) {
  dir <- file.path(root, label)
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)
  tarball <- file.path(root, paste0(label, ".tar"))
  st <- system2("git", c("archive", "--format=tar", "-o", shQuote(tarball),
                         sha, "R"))
  stopifnot(st == 0L)
  utils::untar(tarball, exdir = dir)
  file.path(dir, "R")
}

top_level <- function(rdir) {
  files <- sort(list.files(rdir, pattern = "\\.[Rr]$", full.names = TRUE))
  rows <- list()
  for (f in files) {
    ex <- parse(f, keep.source = FALSE)
    for (i in seq_along(ex)) {
      e <- ex[[i]]
      txt <- paste(deparse(e, width.cutoff = 500L), collapse = "\n")
      is_assign <- is.call(e) && is.name(e[[1]]) &&
        as.character(e[[1]]) %in% c("<-", "=", "<<-") &&
        (is.name(e[[2]]) || is.character(e[[2]]))
      key <- if (is_assign) as.character(e[[2]]) else NA_character_
      body <- if (is_assign) {
        paste(deparse(e[[3]], width.cutoff = 500L), collapse = "\n")
      } else txt
      rows[[length(rows) + 1L]] <- data.frame(
        file = basename(f), index = i, name = key, body = body,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

r_start <- export_r(start_full, "start")
r_head  <- export_r(head_full, "head")
A <- top_level(r_start)
B <- top_level(r_head)

out <- character(0)
say <- function(...) out <<- c(out, paste0(...))

say("Code-identity gate, TASK_cran_prep_docs_legacy_2026-11-12, Part 6 step 2")
say("R version: ", R.version.string)
say("start: ", start_sha, " = ", start_full)
say("head:  ", head_sha, " = ", head_full)
say("R/ files: start ", length(unique(A$file)), ", head ", length(unique(B$file)))
say("top-level expressions: start ", nrow(A), ", head ", nrow(B))
say("assigned names: start ", sum(!is.na(A$name)), " (", length(unique(na.omit(A$name))),
    " distinct), head ", sum(!is.na(B$name)), " (", length(unique(na.omit(B$name))),
    " distinct)")
say("")

# --- named top-level assignments, compared as multisets by name -------------
nA <- split(A$body[!is.na(A$name)], A$name[!is.na(A$name)])
nB <- split(B$body[!is.na(B$name)], B$name[!is.na(B$name)])
only_start <- sort(setdiff(names(nA), names(nB)))
only_head  <- sort(setdiff(names(nB), names(nA)))
both <- intersect(names(nA), names(nB))
changed <- both[vapply(both, function(n) {
  !identical(sort(nA[[n]]), sort(nB[[n]]))
}, logical(1))]

where <- function(D, n) paste(unique(D$file[!is.na(D$name) & D$name == n]),
                              collapse = ", ")
say("Names only at the start (", length(only_start), "):")
for (n in only_start) say("  ", n, "   [", where(A, n), "]")
say("Names only at the head (", length(only_head), "):")
for (n in only_head) say("  ", n, "   [", where(B, n), "]")
say("Names in both with a different deparsed body (", length(changed), "):")
for (n in changed) say("  ", n, "   [start: ", where(A, n), "; head: ", where(B, n), "]")
say("Names in both with an identical deparsed body: ", length(both) - length(changed))
say("")

# --- other top-level expressions (not a named assignment) -------------------
oA <- A$body[is.na(A$name)]
oB <- B$body[is.na(B$name)]
other_only_start <- setdiff(oA, oB)
other_only_head  <- setdiff(oB, oA)
say("Other top-level expressions: start ", length(oA), ", head ", length(oB))
say("  only at the start: ", length(other_only_start))
for (t in other_only_start) say("    ", substr(gsub("\n", " ", t), 1, 300))
say("  only at the head: ", length(other_only_head))
for (t in other_only_head) say("    ", substr(gsub("\n", " ", t), 1, 300))

# An allowed difference: a globalVariables() call that only lost names.
gv_names <- function(t) {
  e <- tryCatch(str2lang(t), error = function(e) NULL)
  if (is.null(e) || !is.call(e)) return(NULL)
  fn <- paste(deparse(e[[1]]), collapse = "")
  if (!fn %in% c("utils::globalVariables", "globalVariables")) return(NULL)
  v <- tryCatch(eval(e[[2]], baseenv()), error = function(e) NULL)
  if (!is.character(v)) return(NULL)
  v
}
other_ok <- TRUE
if (length(other_only_start) || length(other_only_head)) {
  gs <- lapply(other_only_start, gv_names)
  gh <- lapply(other_only_head, gv_names)
  if (any(vapply(gs, is.null, logical(1))) || any(vapply(gh, is.null, logical(1))) ||
      length(gs) != length(gh)) {
    other_ok <- FALSE
  } else {
    # every head call must be a subset of some start call
    other_ok <- all(vapply(gh, function(h) any(vapply(gs, function(s)
      all(h %in% s), logical(1))), logical(1)))
    removed <- setdiff(unlist(gs), unlist(gh))
    say("  globalVariables() names removed: ", length(removed),
        if (length(removed)) paste0(" (", paste(removed, collapse = ", "), ")") else "")
  }
}
say("")

# --- verdicts -----------------------------------------------------------------
g1 <- setequal(only_start, expected_removed)
g2 <- length(only_head) == 0L
g3 <- length(changed) == 0L
g4 <- other_ok
say("Expected removed: ", paste(expected_removed, collapse = ", "))
say("GATE names only at the start = expected removed : ", if (g1) "PASS" else "FAIL")
say("GATE names only at the head = none               : ", if (g2) "PASS" else "FAIL")
say("GATE names in both with different body = none    : ", if (g3) "PASS" else "FAIL")
say("GATE other differences only globalVariables removal: ", if (g4) "PASS" else "FAIL")
say("OVERALL: ", if (g1 && g2 && g3 && g4) "PASS" else "FAIL")

writeLines(out)
invisible(g1 && g2 && g3 && g4)
