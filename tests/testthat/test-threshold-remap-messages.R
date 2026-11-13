# =============================================================================
# test-threshold-remap-messages.R
#
# An identity-scale threshold above 1 for a risk difference (RD) or a rate
# difference (IRD) draws two warnings from forestsearch():
#   - the validator (.validate_outcome_threshold_config(), check 4), and
#   - the remap in SECTION 2B, which replaces the value by the estimand's
#     default.
# Each message must name the measure, its ratio-scale counterpart and (the
# remap) the default that applies.  The defaults are those of ?forestsearch,
# "Resolved defaults, per estimand" (man/forestsearch.Rd):
#   RD  -> c1 = 0.05      IRD -> c1 = 0.01
# The ratio-scale counterparts (?forestsearch, effect_measure): RD -> "OR",
# the binary default and the scale the RD messages already name; IRD -> "IRR",
# the ratio measure of the same rate model (IRR and IRD are listed together,
# for count and for binary outcomes).
# The remapped value itself does not change: the fit records the default.
# =============================================================================

skip_on_cran()

# Small synthetic data, built here so the test does not depend on a helper.
.trm_bin <- function() {
  set.seed(11); N <- 200L
  df <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                   x1 = round(rnorm(N), 2), x2 = round(runif(N, 0, 10), 1))
  df$y <- rbinom(N, 1L, plogis(-0.5 + 0.8 * df$treat * (df$x1 <= 0)))
  df
}

.trm_cnt <- function() {
  set.seed(12); N <- 200L
  df <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                   x1 = round(rnorm(N), 2), x2 = round(runif(N, 0, 10), 1),
                   ftime = round(runif(N, 1, 5), 2))
  df$y <- rpois(N, 0.5 * exp(0.7 * df$treat * (df$x1 <= 0)) * df$ftime)
  df
}

# One forestsearch() fit; every warning's message is collected, messages are
# muffled.  A one-factor search on two confounders keeps each fit to seconds.
.trm_fit <- function(measure, threshold) {
  base <- list(confounders.name = c("x1", "x2"), outcome.name = "y",
               treat.name = "treat", id.name = "id",
               effect_measure = measure, effect.threshold = threshold,
               n.min = 40L, use_grf = FALSE, use_lasso = FALSE, maxk = 1L,
               fs.splits = 20L, details = FALSE, quiet = TRUE,
               plot.sg = FALSE, plot.grf = FALSE,
               parallel_args = list(plan = "sequential", workers = 1L,
                                    show_message = FALSE))
  args <- switch(measure,
    RD  = c(list(df.analysis = .trm_bin(), outcome_type = "binary",
                 event.name = "y"), base),
    IRD = c(list(df.analysis = .trm_cnt(), outcome_type = "count",
                 offset.name = "ftime"), base))
  w <- character(0)
  fit <- withCallingHandlers(
    do.call(forestsearch, args),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  list(fit = fit, warnings = w)
}

# The two messages, picked by their opening words (unchanged on both trees).
.trm_validator <- function(w) {
  grep("^effect\\.threshold = [0-9.]+ but effect_measure", w, value = TRUE)
}
.trm_remap <- function(w) {
  grep("^effect\\.threshold = [0-9.]+ appears to be on a ratio scale", w,
       value = TRUE)
}

test_that("rate difference: both messages name IRD, the IRR scale and the 0.01 default", {
  r <- .trm_fit("IRD", 1.5)
  val <- .trm_validator(r$warnings)
  rem <- .trm_remap(r$warnings)
  expect_length(val, 1L)
  expect_length(rem, 1L)

  # Validator: the measure is a rate difference, its ratio counterpart IRR
  # (?forestsearch effect_measure; see the header).
  expect_match(val, "effect_measure = 'IRD'", fixed = TRUE)
  expect_match(val, "rate[- ]difference")
  expect_match(val, "effect_measure = 'IRR'", fixed = TRUE)
  expect_no_match(val, "risk difference", fixed = TRUE)
  expect_no_match(val, "'OR'", fixed = TRUE)

  # Remap: IRR scale and the IRD default 0.01 ("Resolved defaults, per
  # estimand": IRD c1 = 0.01), printed with %.2f -> "0.01".
  expect_match(rem, "effect_measure = 'IRD'", fixed = TRUE)
  expect_match(rem, "on the IRR scale", fixed = TRUE)
  expect_match(rem, "Remapping to default: 0.01.", fixed = TRUE)
  expect_no_match(rem, "OR scale", fixed = TRUE)
  expect_no_match(rem, "0.05", fixed = TRUE)
  # The percentage-point example belongs to a risk difference.
  expect_no_match(rem, "pct-point", fixed = TRUE)

  # The remapped value is unchanged: the fit records the IRD default 0.01
  # (identity scale, so threshold_config$screening is the value itself).
  expect_identical(r$fit$threshold_config$screening, 0.01)
  expect_identical(r$fit$args_call_all$effect.threshold, 0.01)
})

test_that("risk difference: both messages keep their text (RD, OR scale, 0.05 default)", {
  r <- .trm_fit("RD", 1.5)
  val <- .trm_validator(r$warnings)
  rem <- .trm_remap(r$warnings)
  expect_length(val, 1L)
  expect_length(rem, 1L)

  # The RD messages already name the right measure, scale and default; their
  # text is written out here by hand from the format strings at the start SHA
  # (R/forestsearch_main.R:302-308 and 2335-2347), with effect.threshold = 1.5
  # printed by %.2f as "1.50", the ratio counterpart "OR" (see the header) and
  # the RD default 0.05 ("Resolved defaults, per estimand": RD c1 = 0.05).
  expect_identical(val, paste0(
    "effect.threshold = 1.50 but effect_measure = 'RD' ",
    "(identity scale, range [-1, 1]).\n",
    "  A risk difference cannot exceed 1.0.\n",
    "  Did you mean effect_measure = 'OR' (ratio scale)?"))
  expect_identical(rem, paste0(
    "effect.threshold = 1.50 appears to be on a ratio scale, ",
    "but effect_measure = 'RD' uses the identity scale ",
    "(range [-1, 1]).\n",
    "  Did you mean effect.threshold = 1.50 on the OR scale?\n",
    "  Remapping to default: 0.05.\n",
    "  To use a custom RD threshold, pass a value in [-1, 1] ",
    "(e.g., effect.threshold = 0.07 for a 7 pct-point threshold)."))

  # The remapped value is unchanged: the fit records the RD default 0.05.
  expect_identical(r$fit$threshold_config$screening, 0.05)
  expect_identical(r$fit$args_call_all$effect.threshold, 0.05)
})

test_that("in-range thresholds draw neither message (RD 0.07, IRD 0.02)", {
  # Values in [-1, 1] are used as given: no validator or remap warning, and
  # the fit records the supplied value.
  for (s in list(list(m = "RD", v = 0.07), list(m = "IRD", v = 0.02))) {
    r <- .trm_fit(s$m, s$v)
    expect_length(.trm_validator(r$warnings), 0L)
    expect_length(.trm_remap(r$warnings), 0L)
    expect_identical(r$fit$threshold_config$screening, s$v, info = s$m)
  }
})
