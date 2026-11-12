# A-064 (TASK_cran_prep_code_fixes_2026-11-12, F10): count outcomes with GRF
# are oriented.  With adverse_outcome = TRUE a higher count is worse, so the
# subgroup grf.subg.harm.glm() reports as "harm" must be where treatment
# RAISES the count.  Before the fix the count branch had no flip, so the
# reported subgroup was the region where treatment LOWERS the count, and
# adverse_outcome = TRUE and FALSE gave identical memberships.
#
# Data-generating model (made up, fixed by design, not tuned):
#   n = 1000; region z = 0/1 with 500 subjects each; treatment alternates
#   0/1 within each region, so every region x arm cell has 250 subjects;
#   x1 ~ N(0, 1) and x2 ~ U(0, 1) are pure noise covariates.
#   Y ~ Poisson(mu), person_time = 1 for everyone (log offset 0), with
#     control, both regions : mu = 3
#     treated,  z == 1      : mu = 9   (IRR 3:   treatment RAISES the count)
#     treated,  z == 0      : mu = 1   (IRR 1/3: treatment LOWERS the count)
#   So the region where treatment raises the count is exactly z == 1.
#
# Expected direction, from the roxygen of adverse_outcome ("If TRUE, a higher
# outcome value indicates a worse result ... so the 'harm' subgroup is where
# treatment increases the outcome. If FALSE (default), higher values are
# better"):
#   adverse_outcome = TRUE  -> harm subgroup = where treatment raises the
#                              count = z == 1.
#   adverse_outcome = FALSE -> higher count is better, harm = where
#                              treatment lowers the count = z == 0.
#   The two memberships therefore differ (they point at opposite regions).
#
# "Mostly" is checked as a share >= 0.8 of the identified subgroup's members
# in the expected region.  Why that holds with a wide margin (worked out by
# hand from the Poisson means above, not from any run of the code):
#   * GRF scale "log" (default): E[log(Y + 0.5)] for Poisson(mu), computed as
#     sum(dpois(y, mu) * log(y + 0.5)), is 2.197 (mu 9), 1.101 (mu 3),
#     0.169 (mu 1), with per-subject SD at most 0.71.  So the CATE on the
#     forest's scale is +1.10 in z == 1 and -0.93 in z == 0 (signs reversed
#     after the adverse flip).  A depth-1 policy tree splits on z, the only
#     effect modifier, and assigns control (the harm action) to one whole
#     region.  The DR-score contrast per subject has SD about 2 x 0.71 = 1.4
#     (W.hat = 0.5 doubles the residual).  To pull the z == 1 share of a
#     ~500-member subgroup below 0.8, a depth-2 leaf of >= 125 subjects from
#     the other region would need its mean DR contrast to move by 0.93
#     against a standard error of 1.4 / sqrt(125) = 0.13, i.e. about 7
#     standard errors.  Under the wrong orientation the share is near 0.
#   * GRF scale "identity", frontier selection: the harm-effect of a
#     candidate is mean(control DR) - mean(treated DR) on the forest's scale,
#     i.e. 9 - 3 = 6 in z == 1 and 1 - 3 = -2 in z == 0 (after the adverse
#     flip; reversed without it).  A candidate that does not condition on z
#     mixes the regions about half and half (x1, x2 independent of z), so
#     its expected harm-effect is about (6 - 2) / 2 = 2, far below the
#     effMaxSG band floor 0.9 x 6 = 5.4; a candidate on z <= 0 has about -2.
#     So the selected candidate conditions on z > 0, and all its members are
#     in z == 1 (share 1).  Under adverse_outcome = FALSE the same argument
#     with the signs reversed (z == 0 effect +2, z == 1 effect -6) puts the
#     selection on z <= 0.
#
# Seeds: set.seed() for the data and seedit for GRF are fixed values chosen
# once; the assertions are qualitative and hold by the margins above.

skip_on_cran()

.a064_count_data <- function() {
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

.a064_fit <- function(df, adverse, ...) {
  grf.subg.harm.glm(
    data             = df,
    confounders.name = c("z", "x1", "x2"),
    outcome.name     = "events",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "count",
    effect_measure   = "log_IRR",
    offset.name      = "person_time",
    adverse_outcome  = adverse,
    seedit           = 20261112L,
    ...
  )
}

# Members of the identified subgroup: treat.recommend == 0 (roxygen @return).
.a064_members <- function(res) res$data$treat.recommend == 0L

test_that("count + adverse_outcome: tree path reports where treatment raises the count", {
  df <- .a064_count_data()
  res_adv <- .a064_fit(df, adverse = TRUE)   # tree, log(Y + 0.5): defaults
  res_ben <- .a064_fit(df, adverse = FALSE)

  in_adv <- .a064_members(res_adv)
  in_ben <- .a064_members(res_ben)
  expect_false(is.null(res_adv$sg.harm.id))
  expect_false(is.null(res_ben$sg.harm.id))
  expect_gt(sum(in_adv), 0L)
  expect_gt(sum(in_ben), 0L)

  # Higher count is worse: the harm subgroup lies mostly in z == 1, the region
  # where treatment raises the count (share >= 0.8, margin in the header).
  expect_gte(mean(df$z[in_adv] == 1L), 0.8)
  # Higher count is better: harm = treatment lowers the count = z == 0.
  expect_gte(mean(df$z[in_ben] == 0L), 0.8)
  # The two orientations give different memberships.
  expect_false(identical(res_adv$data$treat.recommend,
                         res_ben$data$treat.recommend))
})

test_that("count + adverse_outcome: frontier path (identity scale) is oriented", {
  df <- .a064_count_data()
  res_adv <- .a064_fit(df, adverse = TRUE,
                       grf_selection = "frontier",
                       grf_count_transform = "identity")
  res_ben <- .a064_fit(df, adverse = FALSE,
                       grf_selection = "frontier",
                       grf_count_transform = "identity")

  in_adv <- .a064_members(res_adv)
  in_ben <- .a064_members(res_ben)
  expect_false(is.null(res_adv$sg.harm.id))
  expect_false(is.null(res_ben$sg.harm.id))
  expect_gt(sum(in_adv), 0L)
  expect_gt(sum(in_ben), 0L)

  expect_gte(mean(df$z[in_adv] == 1L), 0.8)
  expect_gte(mean(df$z[in_ben] == 0L), 0.8)
  expect_false(identical(res_adv$data$treat.recommend,
                         res_ben$data$treat.recommend))
})
