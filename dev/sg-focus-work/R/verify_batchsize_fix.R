#!/usr/bin/env Rscript
# =============================================================================
# Verify the batch_size reproducibility fix
# =============================================================================
# The defect that motivated the sg_focus work: parallel_args$batch_size is a
# performance knob, but before the fix it changed the SELECTED SUBGROUP under
# sg_focus = "hr". Early stopping halted the scan at the first candidate
# reaching stop_threshold and then applied the focus's key to that prefix --
# and how much of the prefix had been scored depended on the batch size.
#
# This script re-checks that, standalone, against the INSTALLED package.
#
#   Rscript dev/sg-focus-work/R/verify_batchsize_fix.R
#
# Exit status 0 if every check passes, 1 otherwise, so it can gate CI.
#
# Related coverage:
#   quarto/smoke-tests/sg_focus_smoke_test.qmd  (section 4, plus 33 other checks)
#   tests/testthat/test-sg-focus-selection.R    (unit level: keys and aliases)

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
})

source("dev/sg-focus-work/R/gbsg_fixture.R")

FOCI    <- c("maxeff", "maxeffCons", "maxcons", "hr", "eff",
             "maxSG", "minSG", "hrMaxSG", "effMaxSG", "hrMinSG", "effMinSG")
BATCHES <- c(1L, 8L, 64L)

# Foci whose selection key contains a Pcons term, hence unsound to early-stop.
# forestsearch() must reset stop_threshold to NULL for each of these.
EXPECT_RESET  <- c("hr", "eff", "maxcons", "maxSG", "minSG",
                   "hrMaxSG", "effMaxSG", "hrMinSG", "effMinSG")
EXPECT_NORESET <- "maxeffCons"   # the only key with no Pcons term

fails <- 0L
say <- function(...) cat(sprintf(...))
report <- function(ok, fmt, ...) {
  if (!isTRUE(ok)) fails <<- fails + 1L
  say("  [%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", sprintf(fmt, ...))
  isTRUE(ok)
}

say("forestsearch %s | R %s.%s\n\n",
    as.character(utils::packageVersion("forestsearch")),
    R.version$major, R.version$minor)

# -----------------------------------------------------------------------------
say("1. Selected subgroup is invariant to parallel_args$batch_size\n")
# Compared via gbsg_pick(), which is the SELECTION only. It deliberately omits
# n_candidates_evaluated / n_passed / early_stop_triggered: those legitimately
# vary with batch_size under a focus that early-stops, because a larger batch
# scores more candidates before the stop condition is noticed. Including them
# would report correct behaviour as a regression.
# -----------------------------------------------------------------------------
for (f in FOCI) {
  fits  <- lapply(BATCHES, function(b) suppressWarnings(gbsg_fs(
    sg_focus = f,
    parallel_args = list(plan = "sequential", workers = 1L, batch_size = b))))
  picks <- lapply(fits, gbsg_pick)
  same  <- all(vapply(picks[-1], identical, logical(1), picks[[1]]))
  evals <- paste(vapply(fits, function(x)
    x$grp.consistency$n_candidates_evaluated, numeric(1)), collapse = "/")
  report(same, "%-11s %-28s n=%-4s evaluated=%s",
         f, picks[[1]]$sg_def, picks[[1]]$n_sel, evals)
}

# -----------------------------------------------------------------------------
say("\n2. Early stopping is enabled for maxeffCons and disabled elsewhere\n")
# -----------------------------------------------------------------------------
for (f in c(EXPECT_RESET, EXPECT_NORESET)) {
  fit  <- suppressWarnings(gbsg_fs(sg_focus = f))
  thr  <- fit$grp.consistency$stop_threshold
  want_null <- f %in% EXPECT_RESET
  ok <- if (want_null) is.null(thr) else !is.null(thr)
  report(ok, "%-11s stop_threshold reaching the engine: %-6s (expected %s)",
         f, if (is.null(thr)) "NULL" else format(thr),
         if (want_null) "NULL" else "a value")
}

# -----------------------------------------------------------------------------
say("\n3. An explicit stop_threshold warns where it is being discarded\n")
# -----------------------------------------------------------------------------
for (f in c(EXPECT_RESET, EXPECT_NORESET)) {
  w <- character(0)
  invisible(withCallingHandlers(
    gbsg_fs(sg_focus = f, stop_threshold = 0.95),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }))
  warned <- any(grepl("reset to NULL", w))
  report(warned == (f %in% EXPECT_RESET),
         "%-11s warned=%-5s (expected %s)", f, warned, f %in% EXPECT_RESET)
}

# -----------------------------------------------------------------------------
say("\n4. maxeffCons: early stopping does not change the answer\n")
# Where early stopping is sound, the prefix sort must recover the same winner
# as a full scan -- for any stop_threshold >= pconsistency.threshold.
# -----------------------------------------------------------------------------
full <- gbsg_pick(gbsg_fs(sg_focus = "maxeffCons", stop_threshold = NULL))
for (thr in list(0.90, 0.95, 0.99)) {
  p <- gbsg_pick(gbsg_fs(sg_focus = "maxeffCons", stop_threshold = thr))
  report(identical(p, full), "stop_threshold=%.2f matches the full scan (%s)",
         thr, full$sg_def)
}

# -----------------------------------------------------------------------------
say("\n5. stop_threshold tracks pconsistency.threshold (no hardcoded 0.90)\n")
# -----------------------------------------------------------------------------
for (p in c(0.90, 0.85, 0.75)) {
  fit <- gbsg_fs(sg_focus = "maxeffCons", pconsistency.threshold = p)
  got <- fit$grp.consistency$stop_threshold
  report(isTRUE(all.equal(got, p)),
         "pconsistency.threshold=%.2f -> stop_threshold=%s", p, format(got))
}

# -----------------------------------------------------------------------------
say("\n%s\n", strrep("-", 66))
if (fails == 0L) {
  say("ALL CHECKS PASSED - batch_size does not affect the selected subgroup\n")
  quit(status = 0L)
} else {
  say("%d CHECK(S) FAILED\n", fails)
  quit(status = 1L)
}
