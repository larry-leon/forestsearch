# =============================================================================
# DINA orientation sign on the use_dina = TRUE screening route
# (TASK_cran_prep_code_fixes_2026-11-12, F11).
#
# Two routes call dina_subgroup():
#   1. .forestsearch_dina_select()   -- subgroup_method = "dina"
#   2. forestsearch(), SECTION 3B    -- use_dina = TRUE, selected_only = TRUE
# Route 1 passes tau_sign = .dina_tau_sign(outcome_type, adverse_outcome), so
# DINA's tau-hat is read in the harm orientation the admission floor uses.
# Route 2 must pass the same sign; otherwise, for continuous outcomes and for
# binary outcomes with adverse_outcome = FALSE, the screening floor and the
# sg_focus ranking read the raw (benefit-oriented) tau-hat and the proposed
# cut points at the benefit region.
#
# The sign is observed by wrapping dina_subgroup() (the pattern of
# test-dina-gaussian-floor.R).  Compute: small (N <= 300) sequential fits.
# =============================================================================

skip_on_cran()

# Run `f()` with dina_subgroup() wrapped; return its value and, for every
# dina_subgroup() call made while it ran, the tau_sign that call used.
# dina_subgroup()'s documented default is `tau_sign = 1` (?dina_subgroup:
# "tau_sign `1` (default) or `-1`"), so a call that does not pass the
# argument uses 1.
.ts_capture <- function(f) {
  signs <- numeric(0)
  orig <- get("dina_subgroup", envir = asNamespace("forestsearch"))
  testthat::local_mocked_bindings(
    dina_subgroup = function(...) {
      a <- list(...)
      signs <<- c(signs, if (is.null(a[["tau_sign"]])) 1 else a[["tau_sign"]])
      orig(...)
    },
    .package = "forestsearch")
  value <- suppressWarnings(suppressMessages(f()))
  list(value = value, signs = signs)
}

# forestsearch() on the use_dina = TRUE selected-cut screening route (GRF and
# LASSO off, so DINA is the only screening proposal).  `adverse_outcome = NULL`
# drops the argument (modifyList), leaving forestsearch()'s own default.
.ts_screen <- function(outcome_type, df, confounders, adverse_outcome = NULL,
                       extra = list()) {
  args <- .fs_args_for(outcome_type, confounders = confounders,
                       extra = c(list(subgroup_method = "consistency",
                                      use_dina = TRUE, use_grf = FALSE,
                                      use_lasso = FALSE,
                                      adverse_outcome = adverse_outcome),
                                 extra))
  .ts_capture(function() do.call(forestsearch, c(list(df.analysis = df), args)))
}


test_that("the screening route calls dina_subgroup() with the outcome's orientation sign", {
  # Expected signs, by hand from the documented rules (not by calling
  # .dina_tau_sign()):
  #   * ?forestsearch, adverse_outcome: NULL resolves to TRUE for binary and
  #     count, FALSE for continuous ("higher values indicate benefit").
  #   * ?dina_subgroup, tau_sign: forestsearch() "passes -1 for continuous and
  #     binary outcomes with adverse_outcome = FALSE", otherwise the default 1;
  #     the .dina_tau_sign() roxygen (R/forestsearch_helpers.R, @noRd) states
  #     that survival, count and adverse_outcome = TRUE keep 1.
  cases <- list(
    # continuous, default orientation (FALSE)            -> -1
    list(type = "continuous", ao = NULL,  expected = -1,
         df = .make_continuous_data(N = 200L, seed = 42L),
         conf = c("age", "biomarker")),
    # continuous, higher y is harm                        ->  1
    list(type = "continuous", ao = TRUE,  expected = 1,
         df = .make_continuous_data(N = 200L, seed = 42L),
         conf = c("age", "biomarker")),
    # binary (OR), event is good                          -> -1
    list(type = "binary",     ao = FALSE, expected = -1,
         df = .make_binary_data(N = 200L, seed = 42L),
         conf = c("age", "wtkg", "biomarker")),
    # binary (OR), default orientation (TRUE)             ->  1
    list(type = "binary",     ao = NULL,  expected = 1,
         df = .make_binary_data(N = 200L, seed = 42L),
         conf = c("age", "wtkg", "biomarker")),
    # count (IRR): never flipped, even with FALSE         ->  1
    list(type = "count",      ao = FALSE, expected = 1,
         df = .make_count_data(N = 200L, seed = 42L),
         conf = c("age", "biomarker")),
    # survival: never flipped                             ->  1
    list(type = "survival",   ao = NULL,  expected = 1,
         df = .make_survival_data(N = 200L, seed = 42L),
         conf = c("age", "noise"))
  )
  for (cs in cases) {
    info <- paste0(cs$type, ", adverse_outcome = ",
                   if (is.null(cs$ao)) "NULL" else cs$ao)
    run <- .ts_screen(cs$type, cs$df, cs$conf, adverse_outcome = cs$ao)
    # One screening call (subgroup_method = "consistency": route 1 is not run).
    expect_length(run$signs, 1L)
    expect_identical(run$signs, cs$expected, info = info)
  }
})

test_that("on continuous data with a known harm region, the screening cut points at it", {
  # Made-up data: treatment lowers y by 2 where x1 <= 0 and raises it by 2
  # where x1 > 0; x2 is a prognostic covariate with no interaction.
  set.seed(2026L)
  N <- 300L
  df <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                   x1 = round(rnorm(N), 3), x2 = round(rnorm(N), 3))
  df$y <- 0.5 * df$x2 + df$treat * ifelse(df$x1 <= 0, -2, 2) + rnorm(N)

  screen_cut <- function(adverse_outcome) {
    run <- .ts_screen("continuous", df, c("x1", "x2"),
                      adverse_outcome = adverse_outcome,
                      extra = list(hr.threshold = NULL, hr.consistency = NULL,
                                   effect.threshold = 1,
                                   dina_args = list(max_depth = 1L)))
    run$value$dina_cuts
  }

  # adverse_outcome = FALSE: higher y is benefit, so harm is where treatment
  # lowers y, i.e. the region x1 <= 0.  A single cut on x1 points at that
  # region iff it keeps the low side: "x1 <= c".
  cut_benefit_y <- screen_cut(adverse_outcome = FALSE)
  # max_depth = 1L: dina_subgroup() selects one single-covariate cut, so one
  # cut string.
  expect_length(cut_benefit_y, 1L)
  expect_match(cut_benefit_y, "^x1 <= ")

  # adverse_outcome = TRUE: higher y is harm, so harm is where treatment
  # raises y, i.e. x1 > 0; the cut keeps the high side: "x1 >= c".
  cut_harm_y <- screen_cut(adverse_outcome = TRUE)
  expect_length(cut_harm_y, 1L)
  expect_match(cut_harm_y, "^x1 >= ")
})
