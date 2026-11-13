# =============================================================================
# CV summary tables use the fit's own orientation (A2, cran-prep part 5).
#
# forestsearch_KfoldOut(outall = TRUE) rebuilt its effect estimator for GLM
# outcomes (the main branch and the no-original-subgroup ITT branch) without
# the fit's adverse_outcome, so it fell back to make_effect_estimator()'s
# formal (TRUE) whatever the fit used.  For a continuous fit at the default
# orientation (adverse_outcome unset -> FALSE: higher y is better) the fit
# reports MD on -y (> 0 = treatment harm) while the CV tables showed MD on y:
# same magnitude, opposite sign.  Same rule as part 4's F12
# (tests/testthat/test-display-orientation.R): the fit's stored
# adverse_outcome (res$cv_args is the fit's args_call_all), else
# forestsearch()'s rule outcome_type %in% c("binary", "count").
#
# Invariant: the effect column ("MD (95% CI)") for a set equals the fit's own
# estimate for that set, sign included.  The arm means and their raw
# difference ("Diff") stay raw.
#
# Expected values are computed by hand from the data (no package code), as in
# test-display-orientation.R:
#   - unadjusted OLS of y on a 0/1 treatment gives the treatment coefficient
#     mean(y | T = 1) - mean(y | T = 0) and the pooled-variance standard
#     error sqrt(s2 * (1/n0 + 1/n1)), s2 = (SS0 + SS1) / (n0 + n1 - 2);
#   - orientation (?forestsearch, @param adverse_outcome): TRUE means higher
#     values indicate harm, FALSE (the continuous default) that higher values
#     indicate benefit; the package's convention is effect > 0 = treatment
#     harm, so under FALSE the fit's MD is the coefficient for -y, i.e.
#     -(mean(y | T = 1) - mean(y | T = 0)), with the same standard error;
#   - the table's limits are est -/+ 1.96 * se (the literal z of
#     analyze_subgroup_glm(), R/summary_utility_functions.R:223).
#   The hand value for the selected subgroup is also checked against the
#   fit's own reported estimate (max_sg_est), which does not depend on the
#   CV tables.
# Tolerance: the tables are built at digits = 8, so a printed value is within
# half a unit of its last digit, 5e-9, of the exact value (plus 1e-12 for the
# decimal parse).  The tie to max_sg_est uses 1e-10, as in F12: max_sg_est is
# the search's maximum candidate estimate (max(hr.out$HR), R/subgroup_search.R;
# an lm.fit of y on treatment within the candidate, with y negated first under
# FALSE by fit_glm_for_subgroup_fast()), and the hand value is the closed form
# of the same coefficient in the same subjects, so the two differ only by
# floating-point rounding.  The selected subgroup is that maximum here
# (sg_focus = "hr"); the fixture block asserts the tie as a premise.
# =============================================================================

.cvo_pa <- list(plan = "sequential", workers = 1L, show_message = FALSE)
.cvo_tol <- 5e-9 + 1e-12

# Mirror twins, as in test-display-orientation.R: `harm` has a planted region
# where treatment RAISES y by 1.5 (age <= 50 & biomarker_hi); `ben` is the
# same data with y negated, so under the continuous default (higher y =
# better) the same region is where treatment harms.
.cvo_data <- function() {
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))
  df$sex <- as.numeric(as.character(df$sex))
  ben <- df; ben$y <- -ben$y
  list(harm = df, ben = ben)
}

.cvo_fit <- function(df, ...) {
  a <- list(
    df.analysis = df,
    confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 0.5, consistency.threshold = 0,
    pconsistency.threshold = 0.5, fs.splits = 20, n.min = 30, maxk = 2,
    use_lasso = FALSE, use_grf = FALSE, sg_focus = "hr", is.RCT = TRUE,
    details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = .cvo_pa)
  a <- utils::modifyList(a, list(...))
  suppressWarnings(suppressMessages(do.call(forestsearch, a)))
}

.cvo_kfold <- function(fit) {
  suppressWarnings(forestsearch_Kfold(fs.est = fit, Kfolds = 3L,
                                      parallel_args = .cvo_pa))
}

# Fits and K-fold results, built once for the file.
.cvo_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      d <- .cvo_data()
      fit_f <- .cvo_fit(d$ben)                            # unset -> FALSE
      fit_t <- .cvo_fit(d$harm, adverse_outcome = TRUE)
      # No subgroup can reach MD >= 10 when y has sd 1, so nothing is
      # identified: forestsearch_KfoldOut() takes its ITT-only branch.
      fit_0 <- .cvo_fit(d$ben, effect.threshold = 10,
                        pconsistency.threshold = 0.99)    # unset -> FALSE
      cache <<- list(
        data = d,
        fit_f = fit_f, kf_f = .cvo_kfold(fit_f),
        fit_t = fit_t, kf_t = .cvo_kfold(fit_t),
        fit_0 = fit_0, kf_0 = .cvo_kfold(fit_0))
    }
    cache
  }
})

# Hand MD (see header): OLS coefficient and pooled SE, oriented; also the raw
# difference in arm means.
.cvo_md <- function(d, adverse) {
  y0 <- d$y[d$treat == 0]; y1 <- d$y[d$treat == 1]
  n0 <- length(y0); n1 <- length(y1)
  s2 <- (sum((y0 - mean(y0))^2) + sum((y1 - mean(y1))^2)) / (n0 + n1 - 2)
  raw <- mean(y1) - mean(y0)
  list(est = if (adverse) raw else -raw, se = sqrt(s2 * (1 / n0 + 1 / n1)),
       raw = raw)
}

# Numbers printed in a cell such as "-1.09411222 (-1.66508098, -0.52314346)".
.cvo_nums <- function(s) {
  as.numeric(regmatches(s, gregexpr("-?[0-9]+\\.[0-9]+", s))[[1]])
}

# One table row against the hand values for a set: effect cell (estimate and
# limits, sign included) oriented; Diff raw.
.cvo_check_row <- function(row, d, adverse, what) {
  ref <- .cvo_md(d, adverse)
  got <- .cvo_nums(row[["MD (95% CI)"]])
  want <- c(ref$est, ref$est - 1.96 * ref$se, ref$est + 1.96 * ref$se)
  expect_length(got, 3L)
  expect_lt(max(abs(got - want)), .cvo_tol, label = paste("MD cell", what))
  expect_identical(sign(got[1]), sign(want[1]), label = paste("MD sign", what))
  # Arm means and their raw difference are not oriented.
  expect_lt(abs(as.numeric(row[["Diff"]]) - ref$raw), .cvo_tol,
            label = paste("Diff raw", what))
}

# Every table of a KfoldOut(outall = TRUE) result with an original subgroup.
.cvo_check_tables <- function(kf, adverse) {
  out <- forestsearch_KfoldOut(kf, outall = TRUE, digits = 8)
  expect_identical(out$outcome_type, "continuous")
  expect_identical(out$effect_measure, "MD")
  d <- kf$resCV
  # ITT: every subject.
  expect_identical(out$itt_tab$Subgroup, "ITT")
  .cvo_check_row(out$itt_tab[1, ], d, adverse, "ITT")
  # Full-data subgroups: row 1 is sg0.name (treat.recommend.original == 0,
  # the selected subgroup H), row 2 sg1.name (== 1, its complement).
  so <- out$SG_tab_original
  expect_identical(so$Subgroup, c(kf$sg0.name, kf$sg1.name))
  .cvo_check_row(so[1, ], d[d$treat.recommend.original == 0, ], adverse, "FA H")
  .cvo_check_row(so[2, ], d[d$treat.recommend.original == 1, ], adverse, "FA Hc")
  # K-fold subgroups (the held-out assignments): both levels occur here.
  expect_setequal(unique(d$treat.recommend), c(0, 1))
  sk <- out$SG_tab_Kfold
  expect_identical(sk$Subgroup, c(kf$sg0.name, kf$sg1.name))
  .cvo_check_row(sk[1, ], d[d$treat.recommend == 0, ], adverse, "KfA H")
  .cvo_check_row(sk[2, ], d[d$treat.recommend == 1, ], adverse, "KfA Hc")
  # The combined table carries the same cells.
  ta <- out$tab_all
  expect_identical(rownames(ta), c("Overall", "FA_0", "KfA_0", "FA_1", "KfA_1"))
  expect_identical(ta[["MD"]], c(out$itt_tab[["MD (95% CI)"]],
                                 so[["MD (95% CI)"]][1], sk[["MD (95% CI)"]][1],
                                 so[["MD (95% CI)"]][2], sk[["MD (95% CI)"]][2]))
  invisible(out)
}


test_that("continuous fit at the default orientation: the fixture is what the tests assume", {
  skip_on_cran()
  fx <- .cvo_fixture()
  fit <- fx$fit_f
  kf <- fx$kf_f
  expect_identical(fit$args_call_all$adverse_outcome, FALSE)  # forestsearch()'s rule
  expect_identical(kf$cv_args$adverse_outcome, FALSE)         # carried into cv_args
  expect_false(is.null(fit$sg.harm))
  # The CV frame's full-data H is the fit's own H (same subjects).
  h_ids <- fit$df.est$id[fit$df.est$treat.recommend == 0]
  expect_setequal(kf$resCV$id[kf$resCV$treat.recommend.original == 0], h_ids)
  # The fit's own estimate for its subgroup is the hand MD on -y, and the
  # planted harm reads as harm (> 0).
  h <- .cvo_md(kf$resCV[kf$resCV$treat.recommend.original == 0, ], adverse = FALSE)
  expect_equal(fit$max_sg_est, h$est, tolerance = 1e-10)
  expect_gt(h$est, 0)
})

test_that("forestsearch_KfoldOut(outall = TRUE): effect cells are the fit's own estimates, sign included", {
  skip_on_cran()
  fx <- .cvo_fixture()
  out <- .cvo_check_tables(fx$kf_f, adverse = FALSE)
  # The selected subgroup's cell is the fit's max_sg_est (to the printed digit).
  got <- .cvo_nums(out$SG_tab_original[["MD (95% CI)"]][1])
  expect_lt(abs(got[1] - fx$fit_f$max_sg_est), .cvo_tol)
})

test_that("no original subgroup: the ITT-only table is the fit's own orientation", {
  skip_on_cran()
  fx <- .cvo_fixture()
  fit <- fx$fit_0
  kf <- fx$kf_0
  # Premise: nothing identified, adverse_outcome resolved to FALSE.
  expect_null(fit$sg.harm)
  expect_identical(kf$cv_args$adverse_outcome, FALSE)
  out <- forestsearch_KfoldOut(kf, outall = TRUE, digits = 8)
  expect_null(out$SG_tab_original)
  expect_null(out$SG_tab_Kfold)
  expect_identical(out$itt_tab$Subgroup, "ITT")
  # The raw ITT difference is far from 0 (both signs are distinguishable at
  # the tolerance), so the sign check is informative.
  ref <- .cvo_md(kf$resCV, adverse = FALSE)
  expect_gt(abs(ref$raw), 1e-3)
  .cvo_check_row(out$itt_tab[1, ], kf$resCV, adverse = FALSE, "ITT (no SG)")
  expect_identical(out$tab_all, out$itt_tab)
})

test_that("no stored adverse_outcome: the CV tables use forestsearch()'s rule (continuous -> FALSE)", {
  skip_on_cran()
  fx <- .cvo_fixture()
  kf <- fx$kf_f
  kf$cv_args$adverse_outcome <- NULL      # an object not built by forestsearch()
  expect_false("adverse_outcome" %in% names(kf$cv_args))
  .cvo_check_tables(kf, adverse = FALSE)
  kf0 <- fx$kf_0
  kf0$cv_args$adverse_outcome <- NULL
  out0 <- forestsearch_KfoldOut(kf0, outall = TRUE, digits = 8)
  .cvo_check_row(out0$itt_tab[1, ], kf0$resCV, adverse = FALSE, "ITT (no SG, rule)")
})

test_that("adverse_outcome = TRUE: the CV tables show the raw MD, which is the fit's own estimate", {
  skip_on_cran()
  fx <- .cvo_fixture()
  fit <- fx$fit_t
  kf <- fx$kf_t
  expect_identical(fit$args_call_all$adverse_outcome, TRUE)
  expect_false(is.null(fit$sg.harm))
  h <- .cvo_md(kf$resCV[kf$resCV$treat.recommend.original == 0, ], adverse = TRUE)
  expect_equal(fit$max_sg_est, h$est, tolerance = 1e-10)
  .cvo_check_tables(kf, adverse = TRUE)
  # Mirror twins: the FALSE fit on -y selects the same subgroup.
  expect_identical(fx$fit_f$sg.harm, fit$sg.harm)
})
