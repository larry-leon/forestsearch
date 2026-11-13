# C3 (TASK_cran_prep_displays_se_docs_2026-11-12, Part C): the GRF CATE
# source note of sg_tables() names the scale of the outcome given to GRF.
#
# Before the fix, sg_tables() on a grf.subg.harm.glm() result printed
#   "**GRF CATE:** H&#x302; = <cate_sg>, H&#x302;<sup>c</sup> = <cate_sgc>
#    (risk-difference scale)"
# for every outcome type, counts included.
#
# What GRF is given for a count, from grf.subg.harm.glm()
# (R/grf_subg_harm_glm.R at the start SHA b06de4e8):
#   * Y_raw = data[[outcome.name]], the per-subject count (:401).
#   * grf_count_transform = "log" (default): Y_grf = log(Y_raw + 0.5) and
#     grf_y_transform = "log(Y + 0.5)" (:410-412).
#     grf_count_transform = "identity": Y_grf = Y_raw and
#     grf_y_transform = "identity" (:421-423).
#   * The offset never enters Y_grf: no division by exposure, no log-offset.
#     Apart from argument checks (:338-366) and a verbose message (:390-391),
#     offset.name does two things: it adds the exposure column to the
#     complete-case filter (:372-373), so rows with a missing exposure are
#     dropped before the forest is fitted, and it is recorded in the result
#     for the GLM effect estimates.
#   * adverse_outcome = TRUE negates Y_grf after grf_y_transform is set
#     (:452-457, part 4's F10), so grf_y_transform does not record the flip.
#   * cate_sg / cate_sgc are means of DR-score contrasts of the forest fitted
#     to Y_grf (:605-606), i.e. differences on the Y_grf scale.
#   * The object records grf_y_transform and offset.name on every return path
#     (roxygen @return :154-159: grf_y_transform is the "transformation
#     actually applied to Y for GRF screening").
#
# Expected source-note text (worked out from the above, not from a run of the
# fixed code).  The leading part and the "%.4f" formatting are those of the
# start SHA (R/summary_utility_functions.R :916-918); only the parenthesis
# changes, and only for outcome_type = "count":
#   count, grf_y_transform "log(Y + 0.5)" -> "(log(Y + 0.5) scale, Y = count"
#   count, grf_y_transform "identity"     -> "(count scale"
#   then ", not offset-adjusted" when offset.name is recorded (the count
#   given to GRF is not divided by exposure), and ")".
#   count, no recognised grf_y_transform  -> "(scale not recorded)".
#   binary (and every other type)         -> "(risk-difference scale)", as at
#   the start SHA.
# The sign flip of adverse_outcome = TRUE is not named: negating Y_grf does
# not change the scale, and the object does not record whether the flip was
# applied (objects made before F10 carry adverse_outcome = TRUE for a count
# without it), so the same scale is expected for both orientations.
#
# Data: the count design of test-grf-count-orientation.R (n = 1000; region z,
# Poisson means 3 / 9 / 1), and a binary analogue (control risk 0.4; treated
# 0.8 in z == 1, 0.1 in z == 0).  In both, treatment moves the outcome in
# opposite directions in the two regions, so the policy tree finds a subgroup
# and cate_sg is stored (the footnote is printed only then).  The premise
# !is.null(fit$cate_sg) is asserted before the footnote.  Seeds are fixed
# values chosen once; nothing below depends on the CATE values themselves.

skip_on_cran()

.c3_count_data <- function() {
  set.seed(20261112L)
  n <- 1000L
  z <- rep(c(0L, 1L), each = n / 2L)
  treat <- rep(c(0L, 1L), times = n / 2L)
  mu <- ifelse(treat == 0L, 3, ifelse(z == 1L, 9, 1))
  data.frame(
    id          = seq_len(n),
    treat       = treat,
    z           = z,
    x1          = stats::rnorm(n),
    x2          = stats::runif(n),
    person_time = 1,
    events      = stats::rpois(n, mu)
  )
}

.c3_binary_data <- function() {
  set.seed(20261113L)
  n <- 1000L
  z <- rep(c(0L, 1L), each = n / 2L)
  treat <- rep(c(0L, 1L), times = n / 2L)
  p <- ifelse(treat == 0L, 0.4, ifelse(z == 1L, 0.8, 0.1))
  data.frame(
    id    = seq_len(n),
    treat = treat,
    z     = z,
    x1    = stats::rnorm(n),
    x2    = stats::runif(n),
    y     = stats::rbinom(n, 1L, p)
  )
}

.c3_count_fit <- function(df, ...) {
  grf.subg.harm.glm(
    data             = df,
    confounders.name = c("z", "x1", "x2"),
    outcome.name     = "events",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "count",
    effect_measure   = "log_IRR",
    seedit           = 20261112L,
    ...
  )
}

# The source notes of Table 1, as the markdown strings sg_tables() built.
.c3_notes <- function(fit) {
  tabs <- sg_tables(fit)
  vapply(tabs$tab_estimates[["_source_notes"]], as.character, character(1))
}

# Leading part and number format as at the start SHA (:917); `scale` is the
# text inside the closing parenthesis.
.c3_expected <- function(fit, scale) {
  sprintf(
    "**GRF CATE:** H&#x302; = %.4f, H&#x302;<sup>c</sup> = %.4f (%s)",
    fit$cate_sg, fit$cate_sgc, scale
  )
}

test_that("count fit, log transform, offset recorded: note names that scale", {
  fit <- .c3_count_fit(.c3_count_data(), offset.name = "person_time")
  # Premises: what the object records, and that a CATE note is printed.
  expect_identical(fit$outcome_type, "count")
  expect_identical(fit$grf_y_transform, "log(Y + 0.5)")
  expect_identical(fit$offset.name, "person_time")
  expect_false(is.null(fit$cate_sg))

  notes <- .c3_notes(fit)
  expect_identical(
    notes,
    .c3_expected(fit, "log(Y + 0.5) scale, Y = count, not offset-adjusted")
  )
  expect_false(any(grepl("risk-difference", notes, fixed = TRUE)))
})

test_that("count fit, adverse_outcome = TRUE: same scale, sign not named", {
  fit <- .c3_count_fit(.c3_count_data(), offset.name = "person_time",
                       adverse_outcome = TRUE)
  expect_true(isTRUE(fit$adverse_outcome))
  # grf_y_transform is set before the flip (header), so it is unchanged.
  expect_identical(fit$grf_y_transform, "log(Y + 0.5)")
  expect_false(is.null(fit$cate_sg))

  expect_identical(
    .c3_notes(fit),
    .c3_expected(fit, "log(Y + 0.5) scale, Y = count, not offset-adjusted")
  )
})

test_that("count fit, identity transform: note names the count scale", {
  fit <- .c3_count_fit(.c3_count_data(), offset.name = "person_time",
                       grf_count_transform = "identity")
  expect_identical(fit$grf_y_transform, "identity")
  expect_identical(fit$offset.name, "person_time")
  expect_false(is.null(fit$cate_sg))

  expect_identical(
    .c3_notes(fit),
    .c3_expected(fit, "count scale, not offset-adjusted")
  )
})

test_that("count fit without an offset: no offset clause", {
  # grf.subg.harm.glm() warns for log_IRR without an offset (:338-345).
  expect_warning(
    fit <- .c3_count_fit(.c3_count_data(), offset.name = NULL),
    "offset.name = NULL"
  )
  expect_null(fit$offset.name)
  expect_identical(fit$grf_y_transform, "log(Y + 0.5)")
  expect_false(is.null(fit$cate_sg))

  expect_identical(
    .c3_notes(fit),
    .c3_expected(fit, "log(Y + 0.5) scale, Y = count")
  )
})

test_that("count object without a recognised grf_y_transform: says so", {
  fit <- .c3_count_fit(.c3_count_data(), offset.name = "person_time")
  # An object that does not record the transform (e.g. built by hand): the
  # note must not name a scale it cannot read from the object.
  fit$grf_y_transform <- NULL
  expect_identical(.c3_notes(fit), .c3_expected(fit, "scale not recorded"))
})

test_that("binary fit: note unchanged (risk-difference scale)", {
  fit <- grf.subg.harm.glm(
    data             = .c3_binary_data(),
    confounders.name = c("z", "x1", "x2"),
    outcome.name     = "y",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "binary",
    seedit           = 20261112L
  )
  expect_identical(fit$outcome_type, "binary")
  expect_false(is.null(fit$cate_sg))

  # Exactly the start-SHA string (:917): GRF on a 0/1 outcome estimates a
  # difference in risks.
  expect_identical(.c3_notes(fit), .c3_expected(fit, "risk-difference scale"))
})
