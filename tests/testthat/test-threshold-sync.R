# ============================================================================
# test-threshold-sync.R
#
# Acceptance tests for TASK_threshold_sync_2026-09-18.md.
#
# forestsearch() resolves its two effect thresholds through branches guarded
# by user_set_threshold / user_set_consistency, and those flags are detected
# with missing() on the legacy spellings (hr.threshold / hr.consistency).
# args_call_all captures the formals, and the bootstrap
# (bootstrap_analysis_dofuture.R) and the cross-validation
# (forestsearch_cross_validation.R) replay that list with EVERY formal
# supplied -- so missing() is FALSE inside a replicate.  Before the fix that
# made a replicate resolve thresholds its parent never used, on every
# identity-scale measure.
#
# These tests execute no fit, no resample and no fold.  They evaluate
# forestsearch()'s own threshold statements -- lifted out of
# body(forestsearch) by helper-threshold-sync.R, so the tests cannot drift
# from the source -- over the argument-list matrix.
# ============================================================================

# The parent-fit resolution, frozen at d0708011 (the pre-fix tree), from
# dev/reports/baseline_threshold_sync_2026-09-18.csv.  The sync may change
# only what a REPLAY resolves; every value here must survive it untouched.
.THRESHOLD_SYNC_BASELINE_PARENT <- utils::read.csv(text = "
binary-OR,legacy,both,custom,0.40546510810816438,0.18232155679395459,1.5,1.2,log
binary-OR,legacy,both,default,0.22314355131420976,0,1.25,1,log
binary-OR,legacy,c1,custom,0.40546510810816438,0,1.5,1,log
binary-OR,legacy,c1,default,0.22314355131420976,0,1.25,1,log
binary-OR,new,both,custom,0.40546510810816438,0.18232155679395459,1.5,1.2,log
binary-OR,new,both,default,0.22314355131420976,0,1.25,1,log
binary-OR,new,c1,custom,0.40546510810816438,0,1.5,1,log
binary-OR,new,c1,default,0.22314355131420976,0,1.25,1,log
binary-OR,none,neither,-,0.22314355131420976,0,1.25,1,log
binary-RD,legacy,both,custom,0.070000000000000007,0.029999999999999999,0.070000000000000007,0.029999999999999999,identity
binary-RD,legacy,both,default,0.050000000000000003,1,0.050000000000000003,1,identity
binary-RD,legacy,c1,custom,0.070000000000000007,0,0.070000000000000007,0,identity
binary-RD,legacy,c1,default,0.050000000000000003,0,0.050000000000000003,0,identity
binary-RD,new,both,custom,0.070000000000000007,0.029999999999999999,0.070000000000000007,0.029999999999999999,identity
binary-RD,new,both,default,0.050000000000000003,1,0.050000000000000003,1,identity
binary-RD,new,c1,custom,0.070000000000000007,0,0.070000000000000007,0,identity
binary-RD,new,c1,default,0.050000000000000003,0,0.050000000000000003,0,identity
binary-RD,none,neither,-,0.050000000000000003,0,0.050000000000000003,0,identity
binary-unset,legacy,both,custom,0.40546510810816438,0.18232155679395459,1.5,1.2,log
binary-unset,legacy,both,default,0.22314355131420976,0,1.25,1,log
binary-unset,legacy,c1,custom,0.40546510810816438,0,1.5,1,log
binary-unset,legacy,c1,default,0.22314355131420976,0,1.25,1,log
binary-unset,new,both,custom,0.40546510810816438,0.18232155679395459,1.5,1.2,log
binary-unset,new,both,default,0.22314355131420976,0,1.25,1,log
binary-unset,new,c1,custom,0.40546510810816438,0,1.5,1,log
binary-unset,new,c1,default,0.22314355131420976,0,1.25,1,log
binary-unset,none,neither,-,0.22314355131420976,0,1.25,1,log
continuous-MD,legacy,both,custom,30,10,30,10,identity
continuous-MD,legacy,both,default,1.25,1,1.25,1,identity
continuous-MD,legacy,c1,custom,30,0,30,0,identity
continuous-MD,legacy,c1,default,1.25,0,1.25,0,identity
continuous-MD,new,both,custom,30,10,30,10,identity
continuous-MD,new,both,default,1.25,1,1.25,1,identity
continuous-MD,new,c1,custom,30,0,30,0,identity
continuous-MD,new,c1,default,1.25,0,1.25,0,identity
continuous-MD,none,neither,-,0,0,0,0,identity
count-IRD,legacy,both,custom,0.02,0.01,0.02,0.01,identity
count-IRD,legacy,both,default,0.01,1,0.01,1,identity
count-IRD,legacy,c1,custom,0.02,0,0.02,0,identity
count-IRD,legacy,c1,default,0.01,0,0.01,0,identity
count-IRD,new,both,custom,0.02,0.01,0.02,0.01,identity
count-IRD,new,both,default,0.01,1,0.01,1,identity
count-IRD,new,c1,custom,0.02,0,0.02,0,identity
count-IRD,new,c1,default,0.01,0,0.01,0,identity
count-IRD,none,neither,-,0.01,0,0.01,0,identity
count-IRR,legacy,both,custom,0.40546510810816438,0.18232155679395459,1.5,1.2,log
count-IRR,legacy,both,default,0.22314355131420976,0,1.25,1,log
count-IRR,legacy,c1,custom,0.40546510810816438,0,1.5,1,log
count-IRR,legacy,c1,default,0.22314355131420976,0,1.25,1,log
count-IRR,new,both,custom,0.40546510810816438,0.18232155679395459,1.5,1.2,log
count-IRR,new,both,default,0.22314355131420976,0,1.25,1,log
count-IRR,new,c1,custom,0.40546510810816438,0,1.5,1,log
count-IRR,new,c1,default,0.22314355131420976,0,1.25,1,log
count-IRR,none,neither,-,0.22314355131420976,0,1.25,1,log
survival,legacy,both,custom,0.40546510810816438,0.18232155679395459,1.5,1.2,log
survival,legacy,both,default,0.22314355131420976,0,1.25,1,log
survival,legacy,c1,custom,0.40546510810816438,0,1.5,1,log
survival,legacy,c1,default,0.22314355131420976,0,1.25,1,log
survival,new,both,custom,0.40546510810816438,0.18232155679395459,1.5,1.2,log
survival,new,both,default,0.22314355131420976,0,1.25,1,log
survival,new,c1,custom,0.40546510810816438,0,1.5,1,log
survival,new,c1,default,0.22314355131420976,0,1.25,1,log
survival,none,neither,-,0.22314355131420976,0,1.25,1,log
", header = FALSE, stringsAsFactors = FALSE, colClasses = "character",
   col.names = c("estimand", "spelling", "subset", "values",
                 "parent_screening", "parent_consistency",
                 "parent_scr_natural", "parent_con_natural", "parent_scale"))

.threshold_sync_parent_table <- function(p) {
  cols <- c("estimand", "spelling", "subset", "values", "parent_screening",
            "parent_consistency", "parent_scr_natural", "parent_con_natural",
            "parent_scale")
  out <- p[order(p$estimand, p$spelling, p$subset, p$values), cols]
  out[] <- lapply(out, function(x) trimws(as.character(x)))
  rownames(out) <- NULL
  out
}


test_that("every replicate resolves what its parent fit resolved", {
  p <- probe_threshold_sync(sync = TRUE)

  expect_identical(nrow(p), 63L)
  offenders <- p[p$violation,
                 c("estimand", "spelling", "subset", "values",
                   "parent_screening", "rep_screening",
                   "parent_consistency", "rep_consistency")]
  expect_identical(
    sum(p$violation), 0L,
    info = paste0("cells where the replay resolves something else:\n",
                  paste(utils::capture.output(print(offenders)),
                        collapse = "\n")))
})


test_that("the pre-fix tree violates on exactly the identity-scale cells", {
  # The bug, as a test: with the sync statements excluded the probe must
  # reproduce the defect the task was opened against -- 15 cells, all
  # identity-scale, none on survival or any ratio measure.  If this ever goes
  # green on its own, the probe has stopped measuring anything.
  b <- probe_threshold_sync(sync = FALSE)

  expect_identical(sum(b$violation), 15L)
  expect_setequal(unique(b$estimand[b$violation]),
                  c("binary-RD", "continuous-MD", "count-IRD"))
  expect_identical(
    sum(b$violation[b$estimand %in%
                      c("survival", "binary-unset", "binary-OR", "count-IRR")]),
    0L)
})


test_that("the sync leaves the parent fit's own resolution untouched", {
  got  <- .threshold_sync_parent_table(probe_threshold_sync(sync = TRUE))
  want <- .THRESHOLD_SYNC_BASELINE_PARENT
  want <- want[order(want$estimand, want$spelling, want$subset, want$values), ]
  rownames(want) <- NULL

  expect_identical(got, want)
})


test_that("ratio and survival replicates are untouched by the sync", {
  b <- probe_threshold_sync(sync = FALSE)
  p <- probe_threshold_sync(sync = TRUE)
  ratio <- b$estimand %in% c("survival", "binary-unset", "binary-OR",
                             "count-IRR")
  cols <- c("rep_screening", "rep_consistency", "rep_scr_natural",
            "rep_con_natural", "rep_scale", "rep_error")
  expect_identical(sum(ratio), 36L)
  expect_identical(b[ratio, cols], p[ratio, cols])
})


test_that("the replayed list carries the parent's resolved naturals", {
  # args_call_all, built the bootstrap's way, must carry NON-NULL values in
  # the two NULL-defaulted spellings, equal to what the parent resolved on the
  # natural scale.  Those spellings are is.null()-detected, so they are the
  # only ones that survive a wrapper.
  p <- probe_threshold_sync(sync = TRUE)

  expect_false(any(is.na(p$acall_effect.threshold)))
  expect_false(any(is.na(p$acall_consistency.threshold)))
  expect_identical(p$acall_effect.threshold,      p$parent_scr_natural)
  expect_identical(p$acall_consistency.threshold, p$parent_con_natural)

  # The spot check the task names: an identity-scale call that left both
  # thresholds at their defaults used to replay consistency 1.0.
  md <- p[p$estimand == "continuous-MD" & p$spelling == "none", ]
  expect_identical(nrow(md), 1L)
  expect_identical(trimws(md$acall_effect.threshold),      "0")
  expect_identical(trimws(md$acall_consistency.threshold), "0")
  expect_identical(trimws(md$rep_consistency),             "0")
})
