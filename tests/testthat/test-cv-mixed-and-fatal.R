# ============================================================================
# test-cv-mixed-and-fatal.R
#
# Tier-1 edge-case regression coverage for the CV surface, hardening the
# loud-failure and no-subgroup fixes:
#
#   1. MIXED found/not-found per-sim metric aggregation.  The per-sim
#      summaries of forestsearch_tenfold() are single-sourced through
#      forestsearch_KfoldOut() for both the identify and the no-identify
#      shape; the aggregation step (do.call(rbind, lapply(...,
#      "sens_metrics_original" / "find_metrics")) in SECTION 5) must keep
#      canonical column names/order and row integrity when the two shapes
#      are mixed in one run.
#
#   2. .collect_cv_results(): the fatal-on-any-refit-error guard must
#      re-raise the FIRST captured error verbatim regardless of its
#      position in the result list, and the end-to-end stop() message
#      must name the failing fold (and sim, for tenfold) correctly.
#
# The deterministic tests construct per-sim inputs by hand and run them
# through the REAL forestsearch_KfoldOut() / .collect_cv_results() code
# paths (no re-implementation).  The end-to-end tests run small
# forestsearch fits via the shared synthetic-DGM helpers and induce refit
# failures with a call-counting mock around forestsearch() so the error
# lands on a chosen fold/sim.
# ============================================================================


.cv_pa_t1 <- list(plan = "sequential", workers = 1L, show_message = FALSE)

# Canonical metric names, in canonical order (forestsearch_KfoldOut()).
.sens_names_t1 <- c("sens_H", "sens_Hc", "ppv_H", "ppv_Hc")
.find_names_t1 <- c("Any", "Exact", "At least 1", "Cov1", "Cov2",
                    "Cov 1 & 2", "Cov1 exact", "Cov2 exact")

# Null-DGM survival fit (no subgroup identified); premise asserted by caller.
.make_null_fs_t1 <- function(N = 160L, extra = list()) {
  df <- .make_survival_data(N = N, HR_harm = 1.0, seed = 7L)
  args <- .fs_args_for(
    "survival", confounders = c("age", "stage", "sex"),
    extra = utils::modifyList(
      list(hr.threshold = 10.0, pconsistency.threshold = 0.99), extra))
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (!is.null(fs$sg.harm)) {
    skip("null DGM unexpectedly identified a subgroup")
  }
  fs
}


# ----------------------------------------------------------------------------
# T1.1 deterministic: mixed identify / no-identify shapes through the exact
# tenfold aggregation step
# ----------------------------------------------------------------------------

test_that("mixed identify/no-identify per-sim shapes aggregate without mislabeling", {
  confs <- c("age", "sex")
  sg_analysis <- c("{age <= 50}", "{sex >= 1}")
  cv_args_min <- list(confounders.name = confs, outcome.name = "time",
                      event.name = "event", treat.name = "treat")
  base_res <- list(cv_args = cv_args_min, sg0.name = "Not recommend",
                   sg1.name = "Recommend", Kfolds = 4L)
  cvidx <- rep(1:4, each = 25L)

  # (a) identify-shaped element: folds 1, 2, 4 found a subgroup; known
  # 2x2 CV-vs-original agreement counts:
  #   rec=0&orig=0: 30, rec=0&orig=1: 10, rec=1&orig=0: 5, rec=1&orig=1: 55
  rec  <- c(rep(0, 40), rep(1, 60))
  orig <- c(rep(0, 30), rep(1, 10), rep(0, 5), rep(1, 55))
  fold_sg1 <- c("{age <= 50}", "{age <= 60}", NA, "{sex >= 1}")
  fold_sg2 <- c("{sex >= 1}", NA, NA, NA)
  res_ident <- c(base_res, list(
    resCV = data.frame(cvindex = cvidx,
                       sg1 = fold_sg1[cvidx], sg2 = fold_sg2[cvidx],
                       treat.recommend = rec, treat.recommend.original = orig,
                       stringsAsFactors = FALSE),
    sg_analysis = sg_analysis))
  out_ident <- forestsearch_KfoldOut(res_ident, outall = FALSE)

  # (b) no-identify element (primary fit found nothing: sg_analysis NULL)
  res_null <- c(base_res, list(
    resCV = data.frame(cvindex = cvidx,
                       sg1 = NA_character_, sg2 = NA_character_,
                       treat.recommend = 1, treat.recommend.original = 1,
                       stringsAsFactors = FALSE)))
  out_null <- forestsearch_KfoldOut(res_null, outall = FALSE)

  # (c) identified primary but ZERO folds found -- the shape a genuinely
  # mixed tenfold run produces for its all-ITT sims
  res_zero <- c(base_res, list(
    resCV = data.frame(cvindex = cvidx,
                       sg1 = NA_character_, sg2 = NA_character_,
                       treat.recommend = 1, treat.recommend.original = orig,
                       stringsAsFactors = FALSE),
    sg_analysis = sg_analysis))
  out_zero <- forestsearch_KfoldOut(res_zero, outall = FALSE)

  # Drive the EXACT aggregation forestsearch_tenfold() SECTION 5 uses.
  valid_results <- list(
    list(sens_metrics_original = out_null$sens_metrics_original,
         find_metrics = out_null$find_metrics, sim_id = 1L),
    list(sens_metrics_original = out_ident$sens_metrics_original,
         find_metrics = out_ident$find_metrics, sim_id = 2L),
    list(sens_metrics_original = out_zero$sens_metrics_original,
         find_metrics = out_zero$find_metrics, sim_id = 3L))
  sens_out <- do.call(rbind, lapply(valid_results, `[[`, "sens_metrics_original"))
  find_out <- do.call(rbind, lapply(valid_results, `[[`, "find_metrics"))

  # Canonical column names AND order; one row per sim
  expect_identical(dim(sens_out), c(3L, 4L))
  expect_identical(colnames(sens_out), .sens_names_t1)
  expect_identical(dim(find_out), c(3L, 8L))
  expect_identical(colnames(find_out), .find_names_t1)

  # Row integrity: each aggregated row is exactly its element's vector
  # (rbind must not silently reorder or recycle across the mixed shapes)
  for (i in seq_along(valid_results)) {
    expect_identical(unname(sens_out[i, ]),
                     unname(valid_results[[i]]$sens_metrics_original))
    expect_identical(unname(find_out[i, ]),
                     unname(valid_results[[i]]$find_metrics))
  }

  # No-identify row (row 1): find all-0, sens all-NA -- in the right row
  expect_identical(unname(find_out[1L, ]), rep(0, 8L))
  expect_true(all(is.na(sens_out[1L, ])))

  # Identify row (row 2): hand-computed values in the right columns.
  #   sens_H  = 30/35 (rec=0 among orig=0), sens_Hc = 55/65,
  #   ppv_H   = 30/40,                      ppv_Hc  = 55/60
  expect_equal(unname(sens_out[2L, ]), c(30/35, 55/65, 30/40, 55/60))
  #   Any = 3/4 folds found; Exact = 1/4 (fold 1 matched both);
  #   At least 1 = 2/4 (folds 1, 4); Cov1/Cov2 covariate-match = 2/4 each;
  #   both = 1/4; Cov1 exact = 1/4; Cov2 exact = 2/4
  expect_equal(unname(find_out[2L, ]),
               c(0.75, 0.25, 0.50, 0.50, 0.50, 0.25, 0.25, 0.50))

  # Zero-found-under-identified-primary row (row 3): find all-0; sens
  # degenerates to the single-level table branch (sens_Hc/ppv_Hc NA,
  # sens_H = 1 and ppv_H = share of orig=0 among the all-ITT frame)
  expect_identical(unname(find_out[3L, ]), rep(0, 8L))
  expect_equal(unname(sens_out[3L, ]), c(1, NA, 0.35, NA))

  # The median summaries (same apply() lines) stay canonically named
  sens_summary <- apply(sens_out, 2, median, na.rm = TRUE)
  find_summary <- apply(find_out, 2, median, na.rm = TRUE)
  expect_identical(names(sens_summary), .sens_names_t1)
  expect_identical(names(find_summary), .find_names_t1)
})


# ----------------------------------------------------------------------------
# T1.1 end-to-end: borderline DGM yielding a genuine sim-level mix
# ----------------------------------------------------------------------------

test_that("tenfold aggregates a genuine sim-level mix of found/not-found", {
  skip_on_cran()

  df <- .make_survival_data(N = 150L, HR_harm = 1.8, seed = 42L)
  args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"),
                       extra = list(pconsistency.threshold = 0.90))
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (is.null(fs$sg.harm)) {
    skip("borderline DGM did not identify a primary subgroup; retune for mixing")
  }

  sims <- 14L
  cv <- suppressWarnings(forestsearch_tenfold(
    fs.est = fs, sims = sims, Kfolds = 3L,
    details = FALSE, parallel_args = .cv_pa_t1))
  .expect_cv_return_shape(cv)

  # The mix must actually have occurred: some sims found in >= 1 fold,
  # some sims found in NO fold.  Otherwise skip rather than false-pass.
  sim_found <- tapply(cv$fold_summary$any_found, cv$fold_summary$sim, max)
  if (!(any(sim_found == 0L) && any(sim_found == 1L))) {
    skip("no sim-level mix of found/not-found occurred; retune for mixing")
  }

  # Aggregated matrices: canonical names/order, one row per sim
  expect_identical(colnames(cv$sens_out), .sens_names_t1)
  expect_identical(colnames(cv$find_out), .find_names_t1)
  expect_identical(nrow(cv$sens_out), sims)
  expect_identical(nrow(cv$find_out), sims)
  expect_identical(nrow(cv$fold_summary), sims * 3L)

  # Row mislabeling check: per-sim "Any" in find_out must equal the
  # fraction of that sim's folds with any_found == 1 in fold_summary --
  # the two outputs derive from the same per-fold subgroup labels, so a
  # row swap or misalignment between sims breaks this equality.
  any_from_folds <- as.numeric(
    tapply(cv$fold_summary$any_found, cv$fold_summary$sim, mean))
  expect_equal(unname(cv$find_out[, "Any"]), any_from_folds)

  # Not-found sims carry zero find metrics; found sims carry Any > 0
  expect_true(all(cv$find_out[sim_found == 0L, "Any"] == 0))
  expect_true(all(cv$find_out[sim_found == 1L, "Any"] > 0))

  # Median summaries well-formed
  expect_identical(names(cv$sens_summary), .sens_names_t1)
  expect_identical(names(cv$find_summary), .find_names_t1)
})


# ----------------------------------------------------------------------------
# T1.2 deterministic: .collect_cv_results() error positions
# ----------------------------------------------------------------------------

test_that(".collect_cv_results re-raises the FIRST error verbatim at any position", {
  e1 <- simpleError("first induced error")
  e2 <- simpleError("second induced error")
  ok <- data.frame(a = 1:2, b = 3:4)

  # Error at first, middle, and last position: verbatim re-raise
  expect_error(forestsearch:::.collect_cv_results(list(e1, ok, ok)),
               "first induced error", fixed = TRUE)
  expect_error(forestsearch:::.collect_cv_results(list(ok, e1, ok)),
               "first induced error", fixed = TRUE)
  expect_error(forestsearch:::.collect_cv_results(list(ok, ok, e1)),
               "first induced error", fixed = TRUE)

  # Multiple errors: the FIRST one wins
  expect_error(forestsearch:::.collect_cv_results(list(ok, e1, e2)),
               "first induced error", fixed = TRUE)
  expect_error(forestsearch:::.collect_cv_results(list(e2, ok, e1)),
               "second induced error", fixed = TRUE)

  # combine = FALSE (the tenfold guard) raises identically
  expect_error(forestsearch:::.collect_cv_results(list(ok, e1), combine = FALSE),
               "first induced error", fixed = TRUE)

  # All-success paths: combine=TRUE rbinds, combine=FALSE passes through
  expect_identical(forestsearch:::.collect_cv_results(list(ok, ok)),
                   rbind(ok, ok))
  expect_identical(forestsearch:::.collect_cv_results(list(ok, ok),
                                                      combine = FALSE),
                   list(ok, ok))
})


test_that(".collect_cv_results: malformed non-error element slips the gate (documented)", {
  # DOCUMENTED CURRENT BEHAVIOR, not an endorsement: an element that is
  # neither an error condition nor a well-formed result (e.g. a bare list
  # lacking the expected columns) passes the inherits(., "error") gate.
  #
  #   * combine = TRUE (forestsearch_Kfold): it reaches do.call(rbind, .)
  #     and fails with base R's cryptic rbind error ("names do not match
  #     previous names") -- the run halts, but without fold context.
  #   * combine = FALSE (forestsearch_tenfold's fatal guard): it passes
  #     the guard SILENTLY; downstream, tenfold's Filter(is.list, .) would
  #     retain it and the metric rbind sees a malformed element.
  #
  # If .collect_cv_results ever learns to reject malformed elements
  # loudly, update this test to assert the new (better) message.
  ok  <- data.frame(a = 1:2, b = 3:4)
  mal <- list(x = 1, y = "not a result row")

  expect_error(forestsearch:::.collect_cv_results(list(ok, mal)),
               "names do not match previous names", fixed = TRUE)

  out <- forestsearch:::.collect_cv_results(list(ok, mal), combine = FALSE)
  expect_identical(out, list(ok, mal))
})


# ----------------------------------------------------------------------------
# T1.2 end-to-end: induced refit failures name the correct fold / sim
# ----------------------------------------------------------------------------

test_that("Kfold refit failure on a non-first fold halts and names that fold", {
  skip_on_cran()

  fs <- .make_null_fs_t1()
  real_fs <- forestsearch
  calls <- 0L
  testthat::local_mocked_bindings(
    forestsearch = function(...) {
      calls <<- calls + 1L
      if (calls == 2L) stop("induced refit failure XYZZY")
      real_fs(...)
    }
  )

  err <- tryCatch(
    suppressWarnings(forestsearch_Kfold(fs.est = fs, Kfolds = 3L,
                                        parallel_args = .cv_pa_t1)),
    error = function(e) e)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)

  # The message must name the failing fold index and total, and carry the
  # verbatim forestsearch() error
  fold_idx <- as.integer(sub(".*FAILED on fold (\\d+) of (\\d+).*", "\\1", msg))
  fold_tot <- as.integer(sub(".*FAILED on fold (\\d+) of (\\d+).*", "\\2", msg))
  expect_identical(fold_idx, 2L)
  expect_identical(fold_tot, 3L)
  expect_match(msg, "induced refit failure XYZZY", fixed = TRUE)
})


test_that("tenfold refit failure names the correct sim AND fold", {
  skip_on_cran()

  fs <- .make_null_fs_t1()
  real_fs <- forestsearch
  calls <- 0L
  # sims = 2, Kfolds = 3, sequential: calls 1-3 are sim 1's folds; the
  # within-sim fold loop aborts on error, so call 5 is sim 2, fold 2.
  testthat::local_mocked_bindings(
    forestsearch = function(...) {
      calls <<- calls + 1L
      if (calls == 5L) stop("induced refit failure PLUGH")
      real_fs(...)
    }
  )

  err <- tryCatch(
    suppressWarnings(forestsearch_tenfold(fs.est = fs, sims = 2L, Kfolds = 3L,
                                          details = FALSE,
                                          parallel_args = .cv_pa_t1)),
    error = function(e) e)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)

  sim_idx  <- as.integer(sub(".*FAILED on sim (\\d+), fold (\\d+) of (\\d+).*", "\\1", msg))
  fold_idx <- as.integer(sub(".*FAILED on sim (\\d+), fold (\\d+) of (\\d+).*", "\\2", msg))
  fold_tot <- as.integer(sub(".*FAILED on sim (\\d+), fold (\\d+) of (\\d+).*", "\\3", msg))
  expect_identical(sim_idx, 2L)
  expect_identical(fold_idx, 2L)
  expect_identical(fold_tot, 3L)
  expect_match(msg, "induced refit failure PLUGH", fixed = TRUE)
})
