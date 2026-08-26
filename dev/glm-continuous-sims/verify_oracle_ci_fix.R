# Behavioural verification of the patched .oracle_md_on().
#
# Extracts the function from a patched .qmd and exercises it on synthetic data
# under both orientations, checking the three properties the defect violated:
#   lo < hi ; width == 2 * 1.96 * se ; names aligned with positions.
#
# Usage: Rscript verify_oracle_ci_fix_v1.R <patched-document.qmd>

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript verify_oracle_ci_fix_v1.R <doc.qmd>")
doc <- args[[1L]]
txt <- paste(readLines(doc, warn = FALSE), collapse = "\n")

m <- regmatches(txt, regexpr("\\.oracle_md_on <- function.*?\\n\\}\\n", txt))
if (!length(m)) stop("could not extract .oracle_md_on() from ", doc)
tmp <- tempfile(fileext = ".R"); writeLines(m, tmp)

treat_name <- "treat_sim"; outcome_name <- "y_sim"; harm_col <- "flag_harm"
source(tmp)

cat("Verification of .oracle_md_on() in", basename(doc), "\n")
cat("--------------------------------------------------------------\n")
set.seed(20260826)
ok <- TRUE
for (adv in c(FALSE, TRUE)) {
  adverse_outcome <- adv
  n  <- 400
  df <- data.frame(treat_sim = stats::rbinom(n, 1, 0.5), flag_harm = 1L)
  df$y_sim <- 100 - 40 * df$treat_sim + stats::rnorm(n, sd = 120)
  v <- .oracle_md_on(df, rep(TRUE, n))

  c1 <- isTRUE(unname(v["lo"]) < unname(v["hi"]))
  c2 <- isTRUE(all.equal(unname(v["hi"] - v["lo"]), unname(2 * 1.96 * v["se"])))
  c3 <- identical(names(v), c("est", "lo", "hi", "se"))
  ok <- ok && c1 && c2 && c3

  cat(sprintf("adverse_outcome = %-5s  est %9.4f  lo %9.4f  hi %9.4f  se %8.4f\n",
              adv, v["est"], v["lo"], v["hi"], v["se"]))
  cat(sprintf("  [%s] lo < hi\n",                       if (c1) "PASS" else "FAIL"))
  cat(sprintf("  [%s] width == 2 * 1.96 * se  (%.4f)\n", if (c2) "PASS" else "FAIL",
              v["hi"] - v["lo"]))
  cat(sprintf("  [%s] names aligned with positions\n\n", if (c3) "PASS" else "FAIL"))
}

g <- .oracle_md_on(data.frame(treat_sim = 1, y_sim = 1, flag_harm = 1), FALSE)
c4 <- identical(names(g), c("est", "lo", "hi", "se")) && all(is.na(g))
ok <- ok && c4
cat(sprintf("[%s] empty-region guard returns named NA in est,lo,hi,se order\n",
            if (c4) "PASS" else "FAIL"))

cat(sprintf("\n%s\n", if (ok) "ALL CHECKS PASSED" else "*** FAILURES ABOVE ***"))
if (!ok) quit(status = 1L)
