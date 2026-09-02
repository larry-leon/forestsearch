#!/usr/bin/env Rscript
# =========================================================================
# PHASE 3 MEMO ACCEPTANCE: compare_subgroup_sims() must be column-
# identical to the memo's former inline machinery on the real committed
# payloads (002f4f37), and both guards must fire.
#
# Usage: Rscript dev/accept_phase3_memo_compare.R [results_dir]
# =========================================================================
suppressMessages(library(forestsearch))
args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else
  "quarto/extreme_subgroups/fixed_random/results"

.fails <- 0L
ok <- function(cond, msg) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
  if (!isTRUE(cond)) .fails <<- .fails + 1L
  invisible(cond)
}

pl_r <- readRDS(file.path(res_dir, "extreme_sims_resample_10000_payload.rds"))
pl_f <- readRDS(file.path(res_dir, "extreme_sims_fixed_10000_payload.rds"))

# ---- transcribed former inline block (verbatim from the memo) ----------
.col_stats <- function(pl) {
  nan_to_na <- function(x) { x[is.nan(x)] <- NA_real_; x }
  data.frame(
    subgroup = colnames(pl$sim_hrs),
    N    = nan_to_na(round(colMeans(pl$sim_ns, na.rm = TRUE))),
    ub2  = nan_to_na(100 * apply(pl$sim_ubs, 2, function(x) mean(x >= 2.0, na.rm = TRUE))),
    ub3  = nan_to_na(100 * apply(pl$sim_ubs, 2, function(x) mean(x >= 3.0, na.rm = TRUE))),
    mUB  = suppressWarnings(apply(pl$sim_ubs, 2, median, na.rm = TRUE)),
    hr05 = nan_to_na(100 * apply(pl$sim_hrs, 2, function(x) mean(x < 0.5, na.rm = TRUE))),
    hr1  = nan_to_na(100 * apply(pl$sim_hrs, 2, function(x) mean(x > 1.0, na.rm = TRUE))),
    mHR  = suppressWarnings(apply(pl$sim_hrs, 2, median, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}
.st_r <- .col_stats(pl_r); .st_f <- .col_stats(pl_f)
cmp_v <- data.frame(
  subgroup = .st_r$subgroup,
  N_r    = .st_r$N,    N_f    = .st_f$N,
  ub2_r  = .st_r$ub2,  ub2_f  = .st_f$ub2,
  ub3_r  = .st_r$ub3,  ub3_f  = .st_f$ub3,
  mUB_r  = .st_r$mUB,  mUB_f  = .st_f$mUB,
  hr05_r = .st_r$hr05, hr05_f = .st_f$hr05,
  hr1_r  = .st_r$hr1,  hr1_f  = .st_f$hr1,
  mHR_r  = .st_r$mHR,  mHR_f  = .st_f$mHR,
  stringsAsFactors = FALSE
)

cmp_p <- compare_subgroup_sims(pl_r, pl_f,
                               expect_designs = c("resample", "fixed"))

ok(identical(names(cmp_v), names(cmp_p)), "column names identical")
for (j in names(cmp_v)) {
  ok(identical(cmp_v[[j]], cmp_p[[j]]),
     sprintf("column '%s' identical on the real 10k payloads", j))
}
ok(identical(attr(cmp_p, "designs"), c(x = "resample", y = "fixed")) &&
     identical(unname(attr(cmp_p, "n_sims")), c(10000L, 10000L)),
   "attributes: designs and n_sims recorded from the payloads")
e <- tryCatch(compare_subgroup_sims(pl_f, pl_r,
       expect_designs = c("resample", "fixed")),
     error = function(e) conditionMessage(e))
ok(grepl("Design labels", e), "Guard 1: swapped payloads rejected")
bad <- pl_f; colnames(bad$sim_hrs)[5] <- "Zzz"
e <- tryCatch(compare_subgroup_sims(pl_r, bad),
     error = function(e) conditionMessage(e))
ok(grepl("Subgroup names differ", e), "Guard 2: mismatched panels rejected")

cat("\n")
if (.fails == 0L) {
  cat("==== PHASE 3 MEMO ACCEPTANCE: ALL GATES PASSED ====\n")
} else {
  cat(sprintf("==== PHASE 3 MEMO ACCEPTANCE: %d GATE(S) FAILED ====\n", .fails))
  quit(status = 1L)
}
