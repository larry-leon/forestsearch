# ============================================================================
# test-duplicate-cut-names.R
#
# A-047 (F8): get_FSdata() decides whether a GRF or DINA cut is "already
# represented" in the candidate pool by comparing NAMES EXACTLY.
#
#   - GRF and DINA merges: the variable name of a cut is the text before
#     "<=" or ">=", trimmed.  The cut is dropped only if that name is exactly
#     a categorical confounder, or if the whole cut string is exactly in the
#     current pool.
#   - Median-removal site (replace_med_grf = TRUE): a continuous confounder
#     loses its median cut only if its name is exactly the variable name of a
#     GRF cut.
#
# Every dataset below is built without random numbers, so every summary
# statistic quoted in a comment can be checked by hand.  Facts used:
#   * is.continuous(): a variable is continuous when it has at least
#     cont.cutoff (default 4) distinct values, otherwise categorical
#     (?is.continuous: "Minimum number of unique values to be considered
#     continuous").
#   * cut_type = "median" cuts each continuous confounder at
#     round(median, 2), written "name <= value"; cut_type = "default" adds
#     cuts at mean, median, Q1 and Q3 (cut_var(); type-7 quantiles,
#     rounded to 1 digit by process_conf_force_expr()).
#   * A cut whose evaluation has a single level is dropped; every cut
#     tested below splits its variable into two non-empty groups (stated
#     per dataset).
# collapse_cuts = FALSE throughout so that the returned pool reflects only
# the merge being tested (no near-redundant coarsening afterwards).
# ============================================================================

skip_on_cran()

# Dataset A: bm categorical, bmi and age continuous.
#   bm  = 0,1,0,1,...         -> 2 distinct values < 4 -> categorical
#   bmi = 25 x 20, 25 x 24, 25 x 28, 25 x 32 -> 4 distinct -> continuous
#   age = 40,50,60,70,80 repeated 20 times   -> 5 distinct -> continuous
#   "bmi <= 25.3" is TRUE for the 50 rows with bmi 20 or 24, FALSE for 50;
#   "bmi >= 25.3" is TRUE for the 50 rows with bmi 28 or 32.  Both valid.
#   Default cuts on bmi (hand): mean = (20+24+28+32)/4 = 26, median = 26,
#   Q1 = 20 + 0.75*(24-20) = 23 (position 25.75), Q3 = 28 + 0.25*(32-28) = 29
#   (position 75.25).  None of these is 25.3, and all default cuts use "<=",
#   so "bmi <= 25.3" / "bmi >= 25.3" can only enter the pool from GRF / DINA.
.f8_data_a <- function() {
  n <- 100L
  data.frame(
    time  = as.numeric(seq_len(n)),
    event = rep(c(1, 0), times = n / 2),
    bm    = rep(c(0, 1), times = n / 2),
    bmi   = rep(c(20, 24, 28, 32), each = n / 4),
    age   = rep(c(40, 50, 60, 70, 80), times = n / 5)
  )
}

# Dataset B: one continuous confounder x.
#   x = 20 x 0, 20 x 1, 20 x 2, 20 x 3, 20 x 30 -> 5 distinct -> continuous
#   median(x): sorted values 50 and 51 are both 2 -> 2, so the median cut
#   string is "x <= 2" (round(2, 2) = 2).
#   "x <= 2"    : TRUE for 60 rows (0,1,2), FALSE for 40 -> valid.
#   "x <= 25.3" : TRUE for 80 rows (0..3),  FALSE for 20 -> valid.
.f8_data_b <- function() {
  n <- 100L
  data.frame(
    time  = as.numeric(seq_len(n)),
    event = rep(c(1, 0), times = n / 2),
    x     = rep(c(0, 1, 2, 3, 30), each = n / 5)
  )
}

# Dataset C: bm and bmi both continuous.
#   bm  = 1,2,3,4 repeated 25 times -> 4 distinct -> continuous;
#         median: sorted values 50 and 51 are 2 and 3 -> 2.5 -> "bm <= 2.5"
#         (TRUE for 50 rows, FALSE for 50 -> valid).
#   bmi = as in dataset A; median 26 -> "bmi <= 26".
.f8_data_c <- function() {
  n <- 100L
  data.frame(
    time  = as.numeric(seq_len(n)),
    event = rep(c(1, 0), times = n / 2),
    bm    = rep(c(1, 2, 3, 4), times = n / 4),
    bmi   = rep(c(20, 24, 28, 32), each = n / 4)
  )
}

.f8_fsdata <- function(df, confounders.name, ...) {
  get_FSdata(
    df.analysis      = df,
    confounders.name = confounders.name,
    outcome.name     = "time",
    event.name       = "event",
    use_lasso        = FALSE,
    collapse_cuts    = FALSE,
    details          = FALSE,
    ...
  )
}


test_that("GRF cut 'bmi <= 25.3' is kept when 'bm' is a categorical confounder", {
  # Name of the cut = "bmi" (text before "<=").  The categorical confounders
  # are {"bm"}; "bmi" is not exactly "bm".  Under cut_type = "default" the
  # pool before the GRF merge holds only the categorical names {"bm"}, so the
  # whole string "bmi <= 25.3" is not in it either -> kept.  Default cuts on
  # bmi are at 26, 26, 23, 29 (header), so it appears exactly once.
  res <- .f8_fsdata(.f8_data_a(), c("bm", "bmi", "age"),
                    use_grf = TRUE, grf_cuts = "bmi <= 25.3")
  expect_identical(sum(res$confs == "bmi <= 25.3"), 1L)
})


test_that("DINA cut 'bmi >= 25.3' is kept when 'bm' is a categorical confounder", {
  # Name of the cut = "bmi" (text before ">=").  Not exactly the categorical
  # "bm"; the pool before the DINA merge is {"bm"} (no GRF), which does not
  # hold the string "bmi >= 25.3" -> kept, exactly once (all default cuts
  # use "<=").
  res <- .f8_fsdata(.f8_data_a(), c("bm", "bmi", "age"),
                    use_grf = FALSE, dina_cuts = "bmi >= 25.3")
  expect_identical(sum(res$confs == "bmi >= 25.3"), 1L)
})


test_that("GRF and DINA cuts on a categorical confounder are still dropped", {
  # Name of "bm <= 0" and of "bm >= 1" = "bm", exactly the categorical
  # confounder bm -> both dropped.  The categorical bm itself stays in the
  # pool once, as the bare name "bm" (evaluated as bm == 1: 50 / 50 rows).
  res <- .f8_fsdata(.f8_data_a(), c("bm", "bmi", "age"),
                    use_grf = TRUE, grf_cuts = "bm <= 0",
                    dina_cuts = "bm >= 1")
  expect_identical(sum(res$confs == "bm <= 0"), 0L)
  expect_identical(sum(res$confs == "bm >= 1"), 0L)
  expect_identical(sum(res$confs == "bm"), 1L)
})


test_that("GRF cut 'x <= 25.3' is kept when the pool holds the median cut 'x <= 2'", {
  # cut_type = "median", replace_med_grf = FALSE: the pool before the GRF
  # merge is the median cut {"x <= 2"} (median 2, header).  Name of
  # "x <= 25.3" = "x", not categorical (x is continuous); the whole string
  # "x <= 25.3" is not "x <= 2" -> kept.  cut_type = "median" adds no
  # default cuts, so the final pool is exactly {"x <= 2", "x <= 25.3"}.
  res <- .f8_fsdata(.f8_data_b(), "x",
                    cut_type = "median", replace_med_grf = FALSE,
                    use_grf = TRUE, grf_cuts = "x <= 25.3")
  expect_setequal(unname(res$confs), c("x <= 2", "x <= 25.3"))
  expect_length(res$confs, 2L)
})


test_that("DINA cut 'x <= 25.3' is kept when the pool holds the GRF cut 'x <= 2'", {
  # cut_type = "median", replace_med_grf = TRUE: the GRF cut "x <= 2" is on
  # x (name "x" exactly), so the median cut of x is removed and the pool
  # before the GRF merge is empty -> the GRF cut is kept: pool {"x <= 2"}.
  # DINA "x <= 25.3": name "x" not categorical; string not "x <= 2" -> kept.
  # Final pool exactly {"x <= 2", "x <= 25.3"}.
  res <- .f8_fsdata(.f8_data_b(), "x",
                    cut_type = "median", replace_med_grf = TRUE,
                    use_grf = TRUE, grf_cuts = "x <= 2",
                    dina_cuts = "x <= 25.3")
  expect_setequal(unname(res$confs), c("x <= 2", "x <= 25.3"))
  expect_length(res$confs, 2L)
})


test_that("an exact duplicate cut is still dropped (appears once)", {
  # What this block pins: the outcome "the duplicate appears once", not the
  # `grf_cuts %in% confs` / `dina_cuts %in% confs` clause of the merges by
  # itself.  Both merges end in confs <- unique(c(confs, grf_cuts_keep)) or
  # unique(c(confs, dina_cuts_keep)) (get_fsdata.R), which drops an exact
  # duplicate anyway: a copy of get_fsdata.R with that clause removed from
  # both merges still passes every check in this file.  The clause's only
  # other effect is that the duplicate stays out of grf_cuts_keep /
  # dina_cuts_keep, and so out of the `protect` list passed to
  # collapse_redundant_cuts(), which runs only with collapse_cuts = TRUE
  # (FALSE here).
  # GRF: the pool holds the median cut "x <= 2"; the GRF cut "x <= 2" is the
  # same whole string -> dropped.  Final pool exactly {"x <= 2"}.
  res_grf <- .f8_fsdata(.f8_data_b(), "x",
                        cut_type = "median", replace_med_grf = FALSE,
                        use_grf = TRUE, grf_cuts = "x <= 2")
  expect_identical(unname(res_grf$confs), "x <= 2")

  # DINA: GRF "x <= 25.3" removes the median cut of x (replace_med_grf) and
  # is kept (empty pool); the DINA cut "x <= 25.3" is the same whole string
  # -> dropped.  Final pool exactly {"x <= 25.3"}.
  res_dina <- .f8_fsdata(.f8_data_b(), "x",
                         cut_type = "median", replace_med_grf = TRUE,
                         use_grf = TRUE, grf_cuts = "x <= 25.3",
                         dina_cuts = "x <= 25.3")
  expect_identical(unname(res_dina$confs), "x <= 25.3")
})


test_that("median removal compares continuous names exactly with GRF cut names", {
  # cut_type = "median", replace_med_grf = TRUE, GRF cut "bmi <= 25.3"
  # (name "bmi").  Continuous confounders bm and bmi (dataset C).
  #   bmi: name exactly "bmi" -> median cut "bmi <= 26" removed.
  #   bm : "bm" is not exactly "bmi" -> median cut "bm <= 2.5" kept.
  # Pool before the GRF merge {"bm <= 2.5"}; "bmi <= 25.3" has name "bmi"
  # (no categorical confounders) and is not in the pool -> kept.
  # Final pool exactly {"bm <= 2.5", "bmi <= 25.3"}.
  res <- .f8_fsdata(.f8_data_c(), c("bm", "bmi"),
                    cut_type = "median", replace_med_grf = TRUE,
                    use_grf = TRUE, grf_cuts = "bmi <= 25.3")
  expect_setequal(unname(res$confs), c("bm <= 2.5", "bmi <= 25.3"))
  expect_length(res$confs, 2L)
})
