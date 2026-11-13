# =============================================================================
# sg_tables() uses the fit's own orientation (A1, cran-prep part 5).
#
# sg_tables() rebuilt its effect estimator without the fit's adverse_outcome,
# in both of its branches, so make_effect_estimator()'s formal (TRUE) applied
# whatever the fit used:
#   - forestsearch() fit (main branch): the resolved value is stored in
#     fs$args_call_all$adverse_outcome;
#   - grf.subg.harm.glm() result (class "grf_glm_result", GRF branch): the
#     value the GRF fit used is stored at the top level, fs$adverse_outcome
#     (the object has no args_call_all).
# For a continuous fit at the default orientation (forestsearch(): unset ->
# FALSE; grf.subg.harm.glm(): default FALSE) the fit's effect is MD on -y
# (> 0 = treatment harm) while the table showed MD on y: same magnitude,
# opposite sign.
#
# Invariant: the table's effect column for a set is the fit's own estimate for
# that set, sign included.  The arm means and their raw difference (Mean(C),
# Mean(T), Diff) stay raw; that is checked too.
#
# Expected values are computed by hand from the data (no package code), as in
# part 4's test-display-orientation.R:
#   - unadjusted OLS of y on a 0/1 treatment gives the treatment coefficient
#     mean(y | T = 1) - mean(y | T = 0) and the pooled-variance standard
#     error sqrt(s2 * (1/n0 + 1/n1)), s2 = (SS0 + SS1) / (n0 + n1 - 2);
#   - orientation (?forestsearch and ?grf.subg.harm.glm, @param
#     adverse_outcome): TRUE means a higher value is worse, FALSE that a
#     higher value is better; the package's convention is effect > 0 =
#     treatment harm, so under FALSE the fit's MD is the coefficient for -y,
#     i.e. -(mean(y | T = 1) - mean(y | T = 0)), with the same standard error;
#   - Diff is the raw mean(y | T = 1) - mean(y | T = 0) under either value.
# The hand value for the selected subgroup is tied to the fit's own estimate:
# forestsearch()'s max_sg_est (equal, sign included); for the GRF result, the
# sign of its harm score -cate_sg (cate_sg is a doubly robust mean, not the
# arm difference, so only the sign is compared).
# Tolerances, set before the first run: the table is printed with
# ndecimals = 6, so an estimate, an arm mean or Diff is within half a unit of
# the last digit (5e-7) of the hand value.  Those checks use 5e-7 + 1e-12;
# the 1e-12 is floating-point slack: the printed decimal is parsed back to
# the nearest double, the reference (the hand value, or max_sg_est) may
# differ from the table's unrounded value by rounding, and the difference is
# taken in doubles, each of the order of 1e-16 here; and expect_lt() is
# strict, so a cell printed at an exact half-unit tie could fail at 5e-7
# alone.  1e-12 is far below the 1e-6 unit of the last digit, so the check
# stays half a unit.  The limits are printed with
# z = 1.96 (analyze_subgroup_glm()); the hand limits use qnorm(0.975), which
# differs by 3.6e-5, i.e. below 2e-5 at these standard errors (< 0.5), so the
# limits are compared within 1e-4.
# =============================================================================

.sgo_pa <- list(plan = "sequential", workers = 1L, show_message = FALSE)
.sgo_conf <- c("age", "biomarker", "biomarker_hi", "sex")

# Mirror twins (part 4's construction): `harm` has a planted region where
# treatment RAISES y by 1.5 (age <= 50 & biomarker_hi); `ben` is the same
# data with y negated, so under the continuous default (higher y = better)
# the same region is where treatment harms.
.sgo_data <- function() {
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))
  df$sex <- as.numeric(as.character(df$sex))
  ben <- df; ben$y <- -ben$y
  list(harm = df, ben = ben)
}

.sgo_fs <- function(df, ...) {
  suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df, confounders.name = .sgo_conf,
    outcome.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 0.5, consistency.threshold = 0,
    pconsistency.threshold = 0.5, fs.splits = 20, n.min = 30, maxk = 2,
    use_lasso = FALSE, use_grf = FALSE, sg_focus = "hr", is.RCT = TRUE,
    details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = .sgo_pa,
    ...)))
}

.sgo_grf <- function(df, ...) {
  suppressWarnings(suppressMessages(grf.subg.harm.glm(
    data = df, confounders.name = .sgo_conf,
    outcome.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", n.min = 30L, ...)))
}

# One fit per orientation and per branch, built once for the file.
.sgo_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      d <- .sgo_data()
      cache <<- list(
        data  = d,
        fs_f  = .sgo_fs(d$ben),                            # unset -> FALSE
        fs_t  = .sgo_fs(d$harm, adverse_outcome = TRUE),
        grf_f = .sgo_grf(d$ben),                           # default FALSE
        grf_t = .sgo_grf(d$harm, adverse_outcome = TRUE))
    }
    cache
  }
})

# Hand MD (see header): OLS coefficient and pooled SE, oriented; raw Diff and
# arm means.
.sgo_md <- function(d, adverse) {
  y0 <- d$y[d$treat == 0]; y1 <- d$y[d$treat == 1]
  n0 <- length(y0); n1 <- length(y1)
  s2 <- (sum((y0 - mean(y0))^2) + sum((y1 - mean(y1))^2)) / (n0 + n1 - 2)
  raw <- mean(y1) - mean(y0)
  list(est = if (adverse) raw else -raw, se = sqrt(s2 * (1 / n0 + 1 / n1)),
       m0 = mean(y0), m1 = mean(y1), diff = raw)
}

# Sets: ITT (all rows), H (treat.recommend == 0), Hc (== 1) of the data the
# table uses: fs$df.est for a forestsearch() fit, fs$data for a GRF result.
.sgo_sets <- function(d) {
  list(ITT = d, H = d[d$treat.recommend == 0, ], Hc = d[d$treat.recommend == 1, ])
}

# Numbers printed in a cell such as "-1.094112 (-1.665081, -0.523143)" or
# "+0.341350".
.sgo_nums <- function(s) {
  as.numeric(regmatches(s, gregexpr("[-+]?[0-9]+\\.[0-9]+", s))[[1]])
}

# The table's display data for the three sets, against the hand values.
.sgo_check <- function(fit, d, adverse) {
  ref <- lapply(.sgo_sets(d), .sgo_md, adverse = adverse)
  tb <- sg_tables(fit, ndecimals = 6)
  td <- as.data.frame(tb$tab_estimates[["_data"]], stringsAsFactors = FALSE)
  expect_true("MD (95% CI)" %in% names(td))
  z <- stats::qnorm(0.975)
  rows <- c(ITT = "ITT", H = "Questionable", Hc = "Recommend")  # default labels
  for (k in names(rows)) {
    r <- td[td$Subgroup == rows[[k]], , drop = FALSE]
    expect_identical(nrow(r), 1L)
    got <- .sgo_nums(r[["MD (95% CI)"]])
    want <- c(ref[[k]]$est, ref[[k]]$est - z * ref[[k]]$se,
              ref[[k]]$est + z * ref[[k]]$se)
    expect_length(got, 3L)
    expect_lt(abs(got[1] - want[1]), 5e-7 + 1e-12, label = paste("MD est", k))
    expect_identical(sign(got[1]), sign(want[1]), label = paste("MD sign", k))
    expect_lt(max(abs(got[2:3] - want[2:3])), 1e-4, label = paste("MD CI", k))
    # Arm means and their raw difference are not oriented.
    raw <- c(.sgo_nums(r[["Mean(C)"]]), .sgo_nums(r[["Mean(T)"]]),
             .sgo_nums(r[["Diff"]]))
    expect_length(raw, 3L)
    expect_lt(max(abs(raw - c(ref[[k]]$m0, ref[[k]]$m1, ref[[k]]$diff))),
              5e-7 + 1e-12, label = paste("raw", k))
  }
  invisible(list(td = td, ref = ref))
}


test_that("the fixtures are what the tests assume", {
  skip_on_cran()
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  fx <- .sgo_fixture()

  # forestsearch(): resolved FALSE stored in args_call_all; the fit's own
  # estimate for its subgroup is the hand MD on -y and reads as harm.
  fit <- fx$fs_f
  expect_identical(fit$args_call_all$adverse_outcome, FALSE)  # forestsearch()'s rule
  expect_false(is.null(fit$sg.harm))
  h <- .sgo_md(.sgo_sets(fit$df.est)$H, adverse = FALSE)
  # max_sg_est is the search's largest candidate estimate (here the selected
  # subgroup's, which this expectation also asserts): the treatment
  # coefficient of an unadjusted least-squares fit (lm.fit()'s QR
  # decomposition) of the oriented outcome on an intercept and the 0/1
  # treatment over the subgroup's rows.  h$est is that coefficient by hand,
  # a difference of arm means (header).  The two are equal in exact
  # arithmetic and are computed by different floating-point operations, so
  # they can differ only by rounding (measured 2.2e-16).  1e-10, relative to
  # |h$est| as expect_equal() applies it, allows that and nothing near a real
  # difference (the opposite orientation is off by 2 * |h$est|).
  expect_equal(fit$max_sg_est, h$est, tolerance = 1e-10)
  expect_gt(h$est, 0)
  expect_identical(fx$fs_t$args_call_all$adverse_outcome, TRUE)
  expect_false(is.null(fx$fs_t$sg.harm))

  # grf.subg.harm.glm(): a "grf_glm_result" with adverse_outcome at the top
  # level and no args_call_all; a subgroup is found and its harm score
  # (-cate_sg, > 0 = treatment harm) is positive.
  g <- fx$grf_f
  expect_s3_class(g, "grf_glm_result")
  expect_identical(g$adverse_outcome, FALSE)                   # GRF default
  expect_null(g$args_call_all)
  expect_false(is.null(g$sg.harm.id))
  expect_gt(-g$cate_sg, 0)
  expect_identical(fx$grf_t$adverse_outcome, TRUE)
  expect_false(is.null(fx$grf_t$sg.harm.id))
})

test_that("sg_tables(), forestsearch() fit: the effect column is the fit's own estimate, sign included", {
  skip_on_cran()
  skip_if_not_installed("gt")
  fx <- .sgo_fixture()
  fit <- fx$fs_f
  out <- .sgo_check(fit, fit$df.est, adverse = FALSE)
  # The selected subgroup's row is max_sg_est, sign included.
  got <- .sgo_nums(out$td[["MD (95% CI)"]][out$td$Subgroup == "Questionable"])[1]
  expect_lt(abs(got - fit$max_sg_est), 5e-7 + 1e-12)
  expect_identical(sign(got), sign(fit$max_sg_est))
})

test_that("sg_tables(), grf.subg.harm.glm() result: the effect column follows the result's stored orientation", {
  skip_on_cran()
  skip_if_not_installed("gt")
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  fx <- .sgo_fixture()
  g <- fx$grf_f
  out <- .sgo_check(g, g$data, adverse = FALSE)
  # The selected subgroup's row reads as harm, as the GRF result's own score
  # (-cate_sg) does.
  got <- .sgo_nums(out$td[["MD (95% CI)"]][out$td$Subgroup == "Questionable"])[1]
  expect_identical(sign(got), sign(-g$cate_sg))
})

test_that("no stored adverse_outcome: sg_tables() uses forestsearch()'s rule (continuous -> FALSE)", {
  skip_on_cran()
  skip_if_not_installed("gt")
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  fx <- .sgo_fixture()
  fit <- fx$fs_f
  fit$args_call_all$adverse_outcome <- NULL       # an object not built by forestsearch()
  expect_false("adverse_outcome" %in% names(fit$args_call_all))
  .sgo_check(fit, fit$df.est, adverse = FALSE)
  g <- fx$grf_f
  g$adverse_outcome <- NULL                        # nor by grf.subg.harm.glm()
  expect_false("adverse_outcome" %in% names(g))
  .sgo_check(g, g$data, adverse = FALSE)
})

test_that("adverse_outcome = TRUE: sg_tables() shows the raw MD, which is the fit's own estimate", {
  skip_on_cran()
  skip_if_not_installed("gt")
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  fx <- .sgo_fixture()
  fit <- fx$fs_t
  h <- .sgo_md(.sgo_sets(fit$df.est)$H, adverse = TRUE)
  # As in the fixture block: the search's least-squares coefficient and the
  # hand difference of arm means (here on y itself) are equal in exact
  # arithmetic, so 1e-10 (relative) leaves room for rounding only.
  expect_equal(fit$max_sg_est, h$est, tolerance = 1e-10)
  .sgo_check(fit, fit$df.est, adverse = TRUE)
  .sgo_check(fx$grf_t, fx$grf_t$data, adverse = TRUE)
})

test_that("mirror twins: the FALSE fit on -y and the TRUE fit on y show the same effect column", {
  skip_on_cran()
  skip_if_not_installed("gt")
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  fx <- .sgo_fixture()
  # Same subgroup (forestsearch() orients -y under FALSE; GRF gives the forest
  # the same Y in both: -y_harm), so the same oriented estimates.
  expect_identical(fx$fs_f$sg.harm, fx$fs_t$sg.harm)
  expect_identical(fx$grf_f$sg.harm.id, fx$grf_t$sg.harm.id)
  eff <- function(fit) {
    td <- sg_tables(fit, ndecimals = 6)$tab_estimates[["_data"]]
    td[["MD (95% CI)"]]
  }
  expect_identical(eff(fx$fs_f), eff(fx$fs_t))
  expect_identical(eff(fx$grf_f), eff(fx$grf_t))
})
