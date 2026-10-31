# vignette_cap_gate.R -- output-only gate on the built primary vignette (2026-10-31).
#
#   Rscript vignette_cap_gate.R <forestsearch.html from the built tarball's inst/doc> <out.txt>
#
# Checks: Table 1 and the declaration table are unchanged cell for cell; the new
# cap table carries the expected cutoffs; the new prose is present once and the
# replaced wording is gone.  The cap diagnostic is recorded against the paper's
# 0.186, not gated.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, file.exists(args[1]))
doc <- xml2::read_html(args[1])
xml2::xml_remove(xml2::xml_find_all(doc, "//script | //style"))
squash <- function(x) trimws(gsub("\\s+", " ", x))
txt <- squash(xml2::xml_text(doc))

ok <- TRUE
out <- character()
put <- function(label, pass, found) {
  ok <<- ok && pass
  out <<- c(out, sprintf("%-78s %s  %s", label, if (pass) "PASS" else "FAIL", found))
}

# ---- tables, located by caption -------------------------------------------
table_by_caption <- function(prefix) {
  tabs <- xml2::xml_find_all(doc, "//table")
  caps <- vapply(tabs, function(t) {
    cap <- xml2::xml_find_first(t, ".//caption")
    if (inherits(cap, "xml_missing")) "" else squash(xml2::xml_text(cap))
  }, character(1))
  hit <- which(startsWith(caps, prefix))
  if (length(hit) != 1L) return(NULL)
  rows <- xml2::xml_find_all(tabs[[hit]], ".//tbody/tr")
  lapply(rows, function(r) squash(xml2::xml_text(xml2::xml_find_all(r, "./th | ./td"))))
}
check_table <- function(label, prefix, expected) {
  got <- table_by_caption(prefix)
  if (is.null(got)) {
    put(label, FALSE, "table not found exactly once")
    return(invisible(NULL))
  }
  same <- length(got) == length(expected) &&
    all(mapply(identical, got, expected))
  bad <- if (same) "" else {
    k <- which(!mapply(function(a, b) identical(a, b),
                       got[seq_len(min(length(got), length(expected)))],
                       expected[seq_len(min(length(got), length(expected)))]))
    sprintf("rows %d/%d; first differing row: %s", length(got), length(expected),
            if (length(k)) paste(got[[k[1]]], collapse = " | ") else "(row count)")
  }
  put(label, same, if (same) sprintf("%d rows identical", length(got)) else bad)
}

check_table("Table 1: every cell unchanged",
            "Selection-adjusted analysis of the three identified regions",
            list(
              c("Region", "{er <= 0} & {pgr <= 26}", "{grade3 >= 1} & {pgr <= 20}", "{er <= 0}"),
              c("N (events)", "75 (41)", "113 (61)", "82 (45)"),
              c("Unadjusted, 95% CI", "2.222 (1.169, 4.221)", "1.315 (0.746, 2.317)", "1.951 (1.037, 3.673)"),
              c("Corrected", "1.276", "1.076", "1.111"),
              c("Field, lower bound", "0.634", "0.638", "0.510"),
              c("Complement N (events)", "611 (258)", "573 (238)", "604 (254)"),
              c("Unadjusted, 95% CI", "0.605 (0.465, 0.788)", "0.641 (0.489, 0.839)", "0.615 (0.472, 0.801)"),
              c("Corrected", "0.645", "0.658", "0.663"),
              c("Field-s, upper bound", "0.804", "0.842", "0.846"),
              c("Bonferroni pair", "0.509 / 0.840", "0.572 / 0.887", "0.404 / 0.874")))

check_table("Declaration table: every cell unchanged",
            "Declaration calibration at alpha = 0.10",
            list(
              c("HR 0.70", "2.578", "0.990", "0.558", "0"),
              c("HR 0.75 (primary)", "2.726", "0.994", "0.667", "0"),
              c("HR 0.80", "2.891", "0.996", "0.769", "0"),
              c("HR 0.85", "3.039", "0.998", "0.844", "0"),
              c("HR 1.00 (claim threshold)", "3.503", "1.000", "0.974", "0")))

# ---- the cap table -----------------------------------------------------------
cap <- table_by_caption("The screen at the recommended cap")
if (is.null(cap) || length(cap) != 5L || !all(lengths(cap) == 2L)) {
  put("cap table: found once, 5 rows x 2 columns", FALSE, "not found or wrong shape")
  fw_cap <- NA_real_
} else {
  put("cap table: found once, 5 rows x 2 columns", TRUE, "")
  v <- vapply(cap, `[`, character(1), 2L)
  put("cap: cutoff at p* = 0.99 is 2.432 (qnorm(0.9925))", identical(v[1], "2.432"), v[1])
  put("cap: T of the declared region is 2.438", identical(v[2], "2.438"), v[2])
  m <- suppressWarnings(as.numeric(v[3]))
  put("cap: margin is positive and below 0.010", isTRUE(m > 0 && m < 0.010), v[3])
  put("cap: uncapped calibrated cutoff is 2.726", identical(v[4], "2.726"), v[4])
  fw_cap <- suppressWarnings(as.numeric(v[5]))
  put("cap: diagnostic at the capped cutoff in (0, 0.667]",
      isTRUE(fw_cap > 0 && fw_cap <= 0.667), v[5])
}

# ---- headings and prose ------------------------------------------------------
h3 <- squash(xml2::xml_text(xml2::xml_find_all(doc, "//h3")))
n_h3 <- sum(startsWith(h3, "The recommended cap"))
put("heading 'The recommended cap' once", n_h3 == 1L, sprintf("found %d", n_h3))

count <- function(s) {
  m <- gregexpr(s, txt, fixed = TRUE)[[1]]
  if (m[1] == -1L) 0L else length(m)
}
for (w in list(
  c("(the studentized complement bound)", 1),
  c("studentized complement bound,", 0),
  c("adds a realized-rate simulation on this design", 0),
  c("on 0.432 of such trials", 1),
  c("on 0.074 at 0.99 and on 0.022 at 1.00", 1),
  c("the three rates are quoted from the", 1),
  c("marginal on either reading", 1),
  c("the narrowest of margins, so the capped screen declares the region", 1),
  c("NaN", 0),
  c("?@", 0),
  c("#> Error", 0))) {
  n <- count(w[1])
  put(sprintf("count \"%s\" == %s", w[1], w[2]), n == as.integer(w[2]), sprintf("found %d", n))
}

out <- c(out, sprintf("recorded (not gated): cap diagnostic %s; the paper reports 0.186 (%s)",
                      format(fw_cap), if (isTRUE(abs(fw_cap - 0.186) < 5e-4)) "MATCH" else "DIFFERS"))
out <- c(out, sprintf("RESULT: %s", if (ok) "OK" else "FAIL"))
writeLines(out, args[2])
cat(out, sep = "\n")
quit(status = if (ok) 0L else 1L)
