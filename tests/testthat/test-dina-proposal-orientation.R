# DINA's proposal floor is applied in the admission floor's orientation
# (TASK_grf_dina_fixes_2026-09-16, P2).  The DINA fits are constructed: a
# "dina" object only needs coefficients, vcov and family here.

.fake_dina_fit <- function(slope, family = "gaussian") {
  structure(
    list(
      coefficients = c("(Intercept)" = 0, x = slope),
      vcov = diag(c(1, 0.01)),
      family = family
    ),
    class = "dina"
  )
}

.orient_df <- function() {
  data.frame(
    id = 1:20,
    x = as.numeric(1:20),
    treat = rep(0:1, 10),
    y = seq(0, 19)
  )
}

# Raw subgroup mean of tau-hat = slope * x over each candidate's set.
.raw_candidate_means <- function(cand, df, slope) {
  vapply(seq_len(nrow(cand)), function(i) {
    x <- df[[cand$v1[i]]]
    inset <- if (cand$d1[i] == "left") x <= cand$c1[i] else x >= cand$c1[i]
    mean(slope * x[inset])
  }, numeric(1))
}

test_that(".dina_tau_sign negates only the outcomes the effect estimator flips", {
  expect_identical(.dina_tau_sign("continuous", FALSE), -1)
  expect_identical(.dina_tau_sign("binary", FALSE), -1)
  expect_identical(.dina_tau_sign("continuous", TRUE), 1)
  expect_identical(.dina_tau_sign("binary", TRUE), 1)
  expect_identical(.dina_tau_sign("count", FALSE), 1)
  expect_identical(.dina_tau_sign("survival", TRUE), 1)
  expect_identical(.dina_tau_sign("survival", FALSE), 1)
})

test_that("under adverse_outcome = FALSE the proposal set lies on the harm side of m_diff", {
  df <- .orient_df()
  fit <- .fake_dina_fit(slope = -5) # harm is the negative raw effect
  raw <- dina_subgroup(fit, df, covariates = "x", m_diff = 30, n_min = 5L,
                       max_depth = 1L, sg_focus = "maxSG")
  expect_false(raw$found)
  expect_identical(raw$n_candidates_qualifying, 0L)

  sg <- dina_subgroup(fit, df, covariates = "x", m_diff = 30, n_min = 5L,
                      max_depth = 1L, sg_focus = "maxSG", tau_sign = -1)
  expect_true(sg$found)
  expect_gt(nrow(sg$candidates), 0)
  expect_true(all(sg$candidates$tau_hat >= 30))
  raw_means <- .raw_candidate_means(sg$candidates, df, slope = -5)
  expect_equal(raw_means, -sg$candidates$tau_hat)
  expect_true(all(raw_means <= -30))
  expect_equal(sg$mean_tau_hat, -mean(-5 * df$x[sg$mask]))
  expect_identical(eval(sg$call$tau_sign), -1) # recorded as the expression -1
})

test_that("under adverse_outcome = TRUE the proposal set is unchanged", {
  df <- .orient_df()
  fit <- .fake_dina_fit(slope = 5)
  sg_default <- dina_subgroup(fit, df, covariates = "x", m_diff = 30,
                              n_min = 5L, max_depth = 1L, sg_focus = "maxSG")
  sg_explicit <- dina_subgroup(fit, df, covariates = "x", m_diff = 30,
                               n_min = 5L, max_depth = 1L, sg_focus = "maxSG",
                               tau_sign = 1)
  expect_identical(sg_explicit, sg_default)
  expect_null(sg_default$call$tau_sign)
  expect_true(sg_default$found)
  raw_means <- .raw_candidate_means(sg_default$candidates, df, slope = 5)
  expect_equal(raw_means, sg_default$candidates$tau_hat)
  expect_true(all(sg_default$candidates$tau_hat >= 30))
})

test_that("tau_sign is validated", {
  df <- .orient_df()
  fit <- .fake_dina_fit(slope = 5)
  expect_error(
    dina_subgroup(fit, df, covariates = "x", m_diff = 30, tau_sign = 2),
    "tau_sign"
  )
})

test_that("forestsearch's DINA selection orients the floor by adverse_outcome", {
  df <- .orient_df()
  fit <- .fake_dina_fit(slope = -5)
  select <- function(adverse_outcome) {
    .forestsearch_dina_select(
      df = df, df.predict = NULL, df.test = NULL,
      confounders.name = "x", outcome.name = "y", event.name = NULL,
      treat.name = "treat", id.name = "id", outcome_type = "continuous",
      hr.threshold = 30, n.min = 5L, sg_focus = "maxSG",
      selection_rule = "neighborhood", effect_neighborhood = 0.10,
      dina_args = list(select_statistic = "dina", max_depth = 1L),
      dina_res = fit, seedit = 1L, details = FALSE,
      adverse_outcome = adverse_outcome
    )
  }
  harm_side <- select(adverse_outcome = FALSE)
  expect_true(harm_side$found)
  expect_true(all(harm_side$grp.consistency$out_sg$candidates$tau_hat >= 30))
  expect_false(select(adverse_outcome = TRUE)$found)
})
