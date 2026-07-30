#!/usr/bin/env Rscript
# =============================================================================
# §5A evidence probe -- stage 1: defaults, arg capture, and the Cox/consistency
# path with debias_gate = TRUE.
# READ-ONLY with respect to the package.  Writes only under dev/terminology-work/.
# =============================================================================

suppressMessages(library(forestsearch))
source("tests/testthat/helper-synthetic-dgm.R")
dir.create("dev/terminology-work/out", showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a
hr <- function(s) cat("\n===== ", s, " =====\n", sep = "")

hr("Q1: forestsearch() formals -- is the gate on out of the box?")
f <- formals(forestsearch)
cat("debias_gate      default :", deparse(f$debias_gate), "\n")
cat("debias_gate_args default :", deparse(f$debias_gate_args), "\n")

CONF <- c("age", "stage", "sex", "noise")
d <- .make_survival_data(N = 400L, HR_harm = 3.0)
mkargs <- function(...) .fs_args_for("survival", confounders = CONF,
  extra = utils::modifyList(
    list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE), list(...)))

hr("Q3: Cox/consistency path, debias_gate = TRUE")
args3 <- mkargs(debias_gate = TRUE,
                debias_gate_args = list(draws = 300L, seed = 11L))
t0 <- proc.time()[3]
r3 <- .run_fs_capture(d, args3)
fs3 <- r3$value
cat(sprintf("elapsed: %.1f s\n", proc.time()[3] - t0))
cat("warnings:", if (!length(r3$warnings)) "(none)" else
      paste(unique(r3$warnings), collapse = " | "), "\n")
cat("sg.harm            :", paste(fs3$sg.harm, collapse = " & "), "\n")
cat("debias_gate NULL?  :", is.null(fs3$debias_gate), "\n")
cat("harm_flag_debiased :", format(fs3$harm_flag_debiased), "\n")
cat("debias_gate in args_call_all? :",
    "debias_gate" %in% names(fs3$args_call_all),
    "value =", format(fs3$args_call_all$debias_gate), "\n")

if (!is.null(fs3$debias_gate)) {
  g <- fs3$debias_gate
  cat("\n-- top-level fields --\n"); str(g, max.level = 1, give.attr = FALSE)
  cat("\n-- gate meta --\n");   str(g$gate)
  cat("\n-- naive --\n");       str(g$naive)
  cat("\n-- debiased --\n");    str(g$debiased)
  cat("\n-- complement --\n");  str(g$complement, max.level = 2)
}
saveRDS(fs3, "dev/terminology-work/out/fs3_cox_consistency.rds")
saveRDS(d,   "dev/terminology-work/out/d_survival.rds")
cat("\n[saved]\n")
