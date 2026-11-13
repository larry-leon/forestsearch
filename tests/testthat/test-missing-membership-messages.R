# F3 (TASK_cran_prep_code_fixes_2026-11-12): messages when subgroup
# membership is missing (NA).  No returned value changes at any site.
#
#   F3a  Simulation classification (.extract_fs_estimates_gen() and
#        .extract_grf_estimates_gen()).  Where NA membership made the code stop
#        with R's "missing value where TRUE/FALSE needed", it stops instead
#        with the number of subjects with no subgroup assignment out of the
#        total, saying a covariate that defines the identified subgroup has
#        missing values.  Exactly the same cases stop; no na.rm is added.
#   F3b  forestsearch(subgroup_method = "grf"): the selected cut applied to
#        every subject gives NA membership -> one warning with the count and
#        the covariates that define the subgroup.
#   F3c  forestsearch_Kfold() / forestsearch_tenfold(): NA held-out
#        assignments -> one warning per call with the count, issued after the
#        folds are collected.
#
# Every input is built by hand, and every expected count, membership vector
# and metric is worked out from that construction in the comments.  No search
# is run: forestsearch()'s GRF route is given a supplied GRF result (so the
# GRF forest is not fitted), and the CV fold refits are replaced by a mock that
# returns a fixed subgroup.

skip_on_cran()

.f3_pa <- list(plan = "sequential", workers = 1L, show_message = FALSE)

# Run `expr`, muffle and record every warning (message and whether it
# carries a call).
.f3_capture <- function(expr) {
  msgs <- character(0)
  has_call <- logical(0)
  value <- withCallingHandlers(
    expr,
    warning = function(cnd) {
      msgs <<- c(msgs, conditionMessage(cnd))
      has_call <<- c(has_call, !is.null(conditionCall(cnd)))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = msgs, has_call = has_call)
}


# =============================================================================
# F3a -- simulation classification extractors
# =============================================================================

# 40 subjects.  Rows 1-10 are truly in H (flag_harm = 1), rows 11-40 are not.
# The identified subgroup (treat.recommend = 0, sg_hat = 1) is rows 1-8 and
# 11-14: 12 members.  Hence, by counting rows:
#   tp = rows 1-8   = 8     fn = rows 9-10  = 2
#   fp = rows 11-14 = 4     tn = rows 15-40 = 26
#   sens = 8/10, spec = 26/30, ppv = 8/12, npv = 26/28,
#   size.H = 12, size.Hc = 28.
.f3a_rec <- function() c(rep(0, 8), rep(1, 2), rep(0, 4), rep(1, 26))

.f3a_data <- function() {
  set.seed(1)
  N <- 40L
  data.frame(id = seq_len(N), y = round(rexp(N), 3),
             ev = rep(c(1L, 1L, 0L, 1L), 10), tr = rep(c(0L, 1L), 20),
             x1 = rep(c(0.2, 0.8), 20),
             grp = ifelse(.f3a_rec() == 0, "1", "0"),
             flag_harm = c(rep(1L, 10), rep(0L, 30)),
             stringsAsFactors = FALSE)
}

# Shape A: NA on truly-H rows (9, 10) -> 2 of 40 with no assignment.
# Shape B: NA only on truly-Hc rows (20, 21, 22) -> 3 of 40.
.f3a_shapes <- list(A = c(9L, 10L), B = c(20L, 21L, 22L))

.f3a_fs <- function(df, fs_res = NULL, fs_full = NULL) {
  .extract_fs_estimates_gen(
    df = df, fs_res = fs_res, fs_full = fs_full, dgm = NULL,
    cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y", event_name = "ev", treat_name = "tr",
    harm_col = "flag_harm", analysis = "FS", verbose = FALSE)
}

.f3a_grf <- function(df, grf_est) {
  .extract_grf_estimates_gen(
    df = df, grf_est = grf_est, dgm = NULL,
    cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y", event_name = "ev", treat_name = "tr",
    harm_col = "flag_harm", analysis = "GRF", verbose = FALSE, debug = FALSE)
}

# The message the task specifies: the count with no assignment out of the
# total, and that a covariate defining the identified subgroup has missing
# values; the covariates are named when the definition is at hand.
.f3a_msg <- function(n_na, n, covs = NULL) {
  sprintf(paste0("%d of %d subjects have no subgroup assignment: a covariate ",
                 "that defines the identified subgroup%s has missing values."),
          n_na, n,
          if (length(covs)) sprintf(" (%s)", paste(covs, collapse = ", ")) else "")
}

test_that("F3a FS extractor, df.est route: both NA shapes stop with the new message", {
  for (s in names(.f3a_shapes)) {
    df <- .f3a_data()
    rec <- .f3a_rec()
    rec[.f3a_shapes[[s]]] <- NA
    df$x1[.f3a_shapes[[s]]] <- NA
    n_na <- length(.f3a_shapes[[s]])                 # 2 (A) or 3 (B)
    # Definition at hand (fs_full$sg.harm = "{x1 <= 0.5}"): x1 is named.
    err <- expect_error(
      .f3a_fs(df, fs_full = list(df.est = data.frame(treat.recommend = rec),
                                 sg.harm = "{x1 <= 0.5}")),
      .f3a_msg(n_na, 40L, "x1"), fixed = TRUE)
    expect_null(conditionCall(err))                  # call. = FALSE
    # No definition at hand (fs_full carries no sg.harm): count alone.
    expect_error(
      .f3a_fs(df, fs_full = list(df.est = data.frame(treat.recommend = rec))),
      .f3a_msg(n_na, 40L), fixed = TRUE)
  }
})

test_that("F3a FS extractor, M.* route: NA membership stops with the new message", {
  # fs_res M.1 = "grp.1": create_subgroup_indicator() puts grp == "1" in the
  # subgroup (rows 1-8, 11-14) and gives NA where grp is NA.  "grp.1" names
  # the column grp (the part before the dot).
  for (s in names(.f3a_shapes)) {
    df <- .f3a_data()
    df$grp[.f3a_shapes[[s]]] <- NA
    expect_error(
      .f3a_fs(df, fs_res = data.frame(M.1 = "grp.1", stringsAsFactors = FALSE)),
      .f3a_msg(length(.f3a_shapes[[s]]), 40L, "grp"), fixed = TRUE)
  }
})

test_that("F3a GRF extractor: both NA shapes stop with the new message", {
  for (s in names(.f3a_shapes)) {
    df <- .f3a_data()
    rec <- .f3a_rec()
    rec[.f3a_shapes[[s]]] <- NA
    df$x1[.f3a_shapes[[s]]] <- NA
    # sg.harm.id "x1 <= 0.5" is the definition at hand: x1 is named.
    err <- expect_error(
      .f3a_grf(df, grf_est = list(sg.harm.id = "x1 <= 0.5",
                                  data = cbind(df, treat.recommend = rec))),
      .f3a_msg(length(.f3a_shapes[[s]]), 40L, "x1"), fixed = TRUE)
    expect_null(conditionCall(err))
  }
})

test_that("F3a: a 0/1 membership completes with the hand-derived values", {
  df <- .f3a_data()
  rec <- .f3a_rec()
  want <- c(sens = 8 / 10, spec = 26 / 30, ppv = 8 / 12, npv = 26 / 28)
  outs <- list(
    `FS df.est` = .f3a_fs(df, fs_full = list(df.est = data.frame(treat.recommend = rec),
                                             sg.harm = "{x1 <= 0.5}")),
    `FS M.*`    = .f3a_fs(df, fs_res = data.frame(M.1 = "grp.1",
                                                  stringsAsFactors = FALSE)),
    GRF         = .f3a_grf(df, grf_est = list(sg.harm.id = "x1 <= 0.5",
                                              data = cbind(df, treat.recommend = rec))))
  for (nm in names(outs)) {
    out <- outs[[nm]]
    # any.H = 1: a subgroup is identified (12 members, rows 1-8 and 11-14).
    expect_identical(out$any.H, 1L, label = paste(nm, "any.H"))
    expect_identical(as.integer(out$size.H), 12L, label = paste(nm, "size.H"))
    expect_identical(as.integer(out$size.Hc), 28L, label = paste(nm, "size.Hc"))
    for (k in names(want))
      expect_equal(out[[k]], want[[k]], label = paste(nm, k))
    # sg_hat = 1 for treat.recommend 0 (rows 1-8, 11-14), else 0
    expect_identical(as.integer(attr(out, "sg_hat")), as.integer(rec == 0),
                     label = paste(nm, "sg_hat"))
  }
})

test_that("F3a: cases that complete today still complete; other stops unchanged", {
  # FS extractor without the harm column: no classification is attempted, so
  # NA membership does not stop.  sg_hat keeps the NA (rows 9, 10); size.H
  # counts the 12 members, size.Hc the 40 - 12 - 2 = 26 non-members.
  df <- .f3a_data()
  df$flag_harm <- NULL
  rec <- .f3a_rec()
  rec[c(9L, 10L)] <- NA
  out <- .f3a_fs(df, fs_full = list(df.est = data.frame(treat.recommend = rec),
                                    sg.harm = "{x1 <= 0.5}"))
  expect_identical(attr(out, "sg_hat"), as.integer(rec == 0))
  expect_identical(as.integer(out$size.H), 12L)
  expect_identical(as.integer(out$size.Hc), 26L)

  # A missing truth label (flag_harm NA) with complete membership is not a
  # membership problem: it still stops with R's own message, not the new one.
  df <- .f3a_data()
  df$flag_harm[5L] <- NA
  expect_error(
    .f3a_fs(df, fs_full = list(df.est = data.frame(treat.recommend = .f3a_rec()),
                               sg.harm = "{x1 <= 0.5}")),
    "missing value where TRUE/FALSE needed", fixed = TRUE)
})


# =============================================================================
# F3b -- forestsearch(subgroup_method = "grf"), selected cut on every subject
# =============================================================================

# 40 subjects.  x1 alternates 0.2 (odd rows) / 0.8 (even rows); x2 cycles
# 0.1, 0.1, 0.9, 0.9.  The supplied GRF definition is {x1 <= 0.5} & {x2 > 0.2},
# so the members are the odd rows with x2 = 0.9: rows 3, 7, 11, ..., 39.
# Missing values, worked out with R's three-valued logic:
#   x1 NA at rows 3, 7 (x2 = 0.9):  NA & TRUE  = NA    -> no assignment
#   x1 NA at row 1    (x2 = 0.1):   NA & FALSE = FALSE -> complement (1)
#   x2 NA at row 11   (x1 = 0.2):   TRUE & NA  = NA    -> no assignment
#   x2 NA at row 4    (x1 = 0.8):   FALSE & NA = FALSE -> complement (1)
# So 3 of 40 subjects (rows 3, 7, 11) have no assignment; members are rows
# 15, 19, 23, 27, 31, 35, 39 (treat.recommend 0); every other row is 1.
.f3b_data <- function(with_na = TRUE) {
  set.seed(7)
  N <- 40L
  df <- data.frame(id = seq_len(N), treat = rep(c(0L, 1L, 1L, 0L), 10),
                   y = rbinom(N, 1L, 0.4),
                   x1 = rep(c(0.2, 0.8), 20),
                   x2 = rep(c(0.1, 0.1, 0.9, 0.9), 10))
  if (with_na) {
    df$x1[c(1L, 3L, 7L)] <- NA
    df$x2[c(4L, 11L)] <- NA
  }
  df
}

.f3b_grf_res <- function() {
  list(sg_def = list(
    conjunctions = list(data.frame(variable = c("x1", "x2"), op = c("<=", ">"),
                                   value = c(0.5, 0.2), stringsAsFactors = FALSE)),
    labels = c("{x1 <= 0.5}", "{x2 > 0.2}"),
    definition = "x1 <= 0.5 & x2 > 0.2", is_disjunction = FALSE))
}

.f3b_fit <- function(df) {
  forestsearch(df.analysis = df, outcome_type = "binary", outcome.name = "y",
               event.name = "y", treat.name = "treat", id.name = "id",
               confounders.name = c("x1", "x2"), effect_measure = "OR",
               subgroup_method = "grf", grf_res = .f3b_grf_res(),
               parallel_args = .f3_pa, details = FALSE, quiet = TRUE,
               plot.sg = FALSE)
}

test_that("F3b: NA membership on the GRF route warns once with count and covariates", {
  run <- .f3_capture(.f3b_fit(.f3b_data()))
  # Exactly one warning, and it is this one: 3 of 40 (rows 3, 7, 11), the
  # covariates of the definition in order of appearance (x1, x2).
  expect_identical(
    run$warnings,
    paste0("3 of 40 subjects have no subgroup assignment (treat.recommend is ",
           "NA): a covariate that defines the GRF-selected subgroup (x1, x2) ",
           "has missing values."))
  expect_identical(run$has_call, FALSE)              # call. = FALSE

  # The returned membership is the hand-derived vector, NA kept.
  want <- rep(1L, 40L)
  want[c(15L, 19L, 23L, 27L, 31L, 35L, 39L)] <- 0L
  want[c(3L, 7L, 11L)] <- NA_integer_
  expect_identical(run$value$df.est$treat.recommend, want)
  expect_identical(run$value$grp.consistency$sg.harm.id,
                   as.integer(want == 0L))
  # The count in the warning equals the number of NA assignments returned.
  expect_identical(sum(is.na(run$value$df.est$treat.recommend)), 3L)
})

test_that("F3b: complete covariates give no warning", {
  run <- .f3_capture(.f3b_fit(.f3b_data(with_na = FALSE)))
  expect_identical(run$warnings, character(0))
  # Members: odd rows with x2 = 0.9 -> rows 3, 7, ..., 39 (10 rows).
  want <- rep(1L, 40L)
  want[seq(3L, 39L, by = 4L)] <- 0L
  expect_identical(run$value$df.est$treat.recommend, want)
})


# =============================================================================
# F3c -- cross-validation held-out assignments
# =============================================================================

# 40 subjects; x1 alternates 0.2 / 0.8 and is missing at rows 5, 12, 33; x2
# cycles 0.1, 0.1, 0.9, 0.9 and is complete.  The full-data fit is
# forestsearch()'s GRF route with a supplied definition on x2 alone, so it has
# no NA membership (no F3b warning).  Each fold's refit is mocked to return
# the subgroup {x1 <= 0.5}; get_dfpred() then gives a held-out subject NA
# exactly when x1 is NA.  Every subject is held out once per K-fold pass, so:
#   forestsearch_Kfold():                 3 of 40 NA (ids 5, 12, 33)
#   forestsearch_tenfold(sims = 2):  2 x 3 = 6 of 80 NA
.f3c_data <- function(with_na = TRUE) {
  set.seed(3)
  N <- 40L
  df <- data.frame(id = seq_len(N), treat = rep(c(0L, 1L), 20),
                   time = round(rexp(N, 0.1), 2),
                   event = rep(c(1L, 1L, 0L, 1L), 10),
                   x1 = rep(c(0.2, 0.8), 20),
                   x2 = rep(c(0.1, 0.1, 0.9, 0.9), 10))
  if (with_na) df$x1[c(5L, 12L, 33L)] <- NA
  df
}

.f3c_fs <- function(df) {
  sg_def <- list(conjunctions = list(data.frame(variable = "x2", op = ">",
                                                value = 0.2,
                                                stringsAsFactors = FALSE)),
                 labels = "{x2 > 0.2}", definition = "x2 > 0.2",
                 is_disjunction = FALSE)
  forestsearch(df.analysis = df, outcome.name = "time", event.name = "event",
               treat.name = "treat", id.name = "id",
               confounders.name = c("x1", "x2"), subgroup_method = "grf",
               grf_res = list(sg_def = sg_def), parallel_args = .f3_pa,
               details = FALSE, quiet = TRUE, plot.sg = FALSE)
}

.f3c_msg_kfold <- paste0(
  "3 of 40 held-out subjects have no subgroup assignment (treat.recommend is ",
  "NA): a covariate that defines their fold's subgroup has missing values.")
.f3c_msg_tenfold <- paste0(
  "6 of 80 held-out assignments across 2 simulations have no subgroup ",
  "assignment (treat.recommend is NA): a covariate that defines the fold's ",
  "subgroup has missing values.")

test_that("F3c forestsearch_Kfold: one warning after the folds, with the count", {
  fs <- .f3c_fs(.f3c_data())
  expect_false(anyNA(fs$df.est$treat.recommend))
  testthat::local_mocked_bindings(
    forestsearch = function(...) list(sg.harm = "{x1 <= 0.5}"))
  run <- .f3_capture(forestsearch_Kfold(fs, Kfolds = 4L, parallel_args = .f3_pa))
  # The fold loop runs inside suppressWarnings(), so a warning that reaches
  # the caller at all was raised after the folds were collected.
  expect_identical(run$warnings, .f3c_msg_kfold)
  expect_identical(run$has_call, FALSE)              # call. = FALSE
  cv <- run$value$resCV
  expect_identical(sort(cv$id[is.na(cv$treat.recommend)]), c(5L, 12L, 33L))
  # Held-out x1 = 0.2 -> 0, x1 = 0.8 -> 1: the non-NA assignments follow x1.
  x1 <- .f3c_data()$x1[match(cv$id, seq_len(40L))]
  ok <- !is.na(x1)
  expect_identical(cv$treat.recommend[ok], ifelse(x1[ok] <= 0.5, 0L, 1L))
})

test_that("F3c forestsearch_tenfold: one warning per call, count over all sims", {
  fs <- .f3c_fs(.f3c_data())
  testthat::local_mocked_bindings(
    forestsearch = function(...) list(sg.harm = "{x1 <= 0.5}"))
  run <- .f3_capture(forestsearch_tenfold(fs, sims = 2L, Kfolds = 4L,
                                          details = FALSE,
                                          parallel_args = .f3_pa,
                                          keep_resCV = TRUE))
  expect_identical(run$warnings, .f3c_msg_tenfold)
  expect_identical(run$has_call, FALSE)
  # The returned components are the documented ones (@return), nothing added.
  expect_identical(names(run$value),
                   c("sens_summary", "find_summary", "sens_out", "find_out",
                     "fold_summary", "resCV_all", "timing_minutes", "sims",
                     "Kfolds", "mr_replicates"))
  for (r in run$value$resCV_all)
    expect_identical(sort(r$id[is.na(r$treat.recommend)]), c(5L, 12L, 33L))
})

test_that("F3c: complete held-out assignments give no warning", {
  fs <- .f3c_fs(.f3c_data(with_na = FALSE))
  testthat::local_mocked_bindings(
    forestsearch = function(...) list(sg.harm = "{x1 <= 0.5}"))
  k <- .f3_capture(forestsearch_Kfold(fs, Kfolds = 4L, parallel_args = .f3_pa))
  expect_identical(k$warnings, character(0))
  expect_false(anyNA(k$value$resCV$treat.recommend))
  t10 <- .f3_capture(forestsearch_tenfold(fs, sims = 2L, Kfolds = 4L,
                                          details = FALSE,
                                          parallel_args = .f3_pa))
  expect_identical(t10$warnings, character(0))
})
