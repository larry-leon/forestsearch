# =============================================================================
# fpr_calibration(): outcome type of the P1 approximation (audit item A-107).
#
# When neither outcome_type_approx nor fs_params$outcome_type is given, the
# approximation type is "survival", which is forestsearch()'s own default:
# forestsearch() declares outcome_type = c("survival", "binary",
# "continuous", "count") and resolves it with match.arg(), which returns the
# first element when the argument is not supplied.  The replicates that
# fpr_calibration() runs are therefore survival analyses, and the P1
# approximation must use the survival effective information.  An
# unrecognised type stops with an error that names it.
#
# The approximation type and d_eff are resolved before the replicate loop and
# do not depend on it, so the %dofuture% loop is stubbed out (the same device
# as test-plan-restored.R): each iteration reports "no subgroup found".
# =============================================================================

skip_on_cran()

# One "no subgroup found" result per iteration of the foreach loop.
.fc_fake_dofuture <- function(obj, ex) {
  rep(list(list(found = 0L)), length(obj$args[[1L]]))
}

# Survival-shaped permutation-mode inputs.  Only the treatment column is read
# before the (stubbed) replicate loop.
.fc_data <- data.frame(
  id    = seq_len(120L),
  treat = rep(0:1, 60L),
  time  = rep(c(5, 10, 15), 40L),
  event = rep(c(1L, 0L), 60L),
  x1    = rep(c(0L, 1L, 1L, 0L), 30L)
)
.fc_fs_params <- list(                # survival fs_params, no outcome_type
  outcome.name = "time", event.name = "event",
  treat.name = "treat", id.name = "id",
  confounders.name = "x1"
)

.fc_call <- function(...) {
  fpr_calibration(c1 = 1.25, c2 = 1.0, df.analysis = .fc_data,
                  nsims = 2L, n_min = 60L, quiet = TRUE, ...)
}


test_that("no type given: the approximation type is survival", {
  local_mocked_bindings(`%dofuture%` = .fc_fake_dofuture,
                        .package = "forestsearch")

  # Survival fs_params without outcome_type, prop_cens = 0.25.
  res <- .fc_call(fs_params = .fc_fs_params, prop_cens = 0.25)

  # Expected type: forestsearch()'s default outcome_type, "survival" (first
  # element of its outcome_type formal, taken by match.arg()).
  expect_identical(res$call_args$otype, "survival")
  # Expected d_eff: d_eff_survival() at n_min = 60, prop_cens = 0.25,
  # evaluated here.  By hand from its documentation, d_eff = n_H (1 - p_c)
  # = 60 * (1 - 0.25) = 45.  (The binary formula used before the fix gives
  # n_H p (1 - p) = 60 * 0.5 * 0.5 = 15 and ignores prop_cens.)
  expect_equal(res$d_eff, d_eff_survival(n_sg = 60, prop_cens = 0.25))
  expect_equal(res$d_eff, 45)
})

test_that("all defaults (fs_params = list()): the approximation type is survival", {
  local_mocked_bindings(`%dofuture%` = .fc_fake_dofuture,
                        .package = "forestsearch")

  res <- .fc_call()

  # Expected type: forestsearch()'s default, "survival" (as above).
  expect_identical(res$call_args$otype, "survival")
  # Expected d_eff: d_eff_survival() at n_min = 60 and the default
  # prop_cens = 0 (fpr_calibration()'s formals), evaluated here.  By hand,
  # 60 * (1 - 0) = 60 (binary formula before the fix: 15).
  expect_equal(res$d_eff, d_eff_survival(n_sg = 60, prop_cens = 0))
  expect_equal(res$d_eff, 60)
})

test_that("a misspelt outcome type stops with an error naming it", {
  local_mocked_bindings(`%dofuture%` = .fc_fake_dofuture,
                        .package = "forestsearch")

  # Documented values are "binary", "survival", "count", "continuous"
  # (roxygen of outcome_type_approx); "survial" is none of them, so each call
  # must stop, and the message must name the value given ("survial") and list
  # the recognised ones (so it contains "survival").
  expect_error(.fc_call(fs_params = .fc_fs_params,
                        outcome_type_approx = "survial"),
               regexp = "survial", fixed = TRUE)
  expect_error(.fc_call(fs_params = .fc_fs_params,
                        outcome_type_approx = "survial"),
               regexp = "survival", fixed = TRUE)

  # The same misspelling arriving through fs_params$outcome_type.
  fs_bad <- c(.fc_fs_params, list(outcome_type = "survial"))
  expect_error(.fc_call(fs_params = fs_bad), regexp = "survial", fixed = TRUE)
  expect_error(.fc_call(fs_params = fs_bad), regexp = "survival", fixed = TRUE)

  # An unrecognised type stops even when d_eff_approx is supplied.
  expect_error(.fc_call(fs_params = .fc_fs_params,
                        outcome_type_approx = "survial", d_eff_approx = 30),
               regexp = "survial", fixed = TRUE)
})

test_that("a type that is given is still used", {
  local_mocked_bindings(`%dofuture%` = .fc_fake_dofuture,
                        .package = "forestsearch")

  # fs_params$outcome_type = "binary", default p_event = 0.5: d_eff by hand
  # from d_eff_binary()'s documentation, n_H p (1 - p) = 60 * 0.5 * 0.5 = 15.
  res_b <- .fc_call(fs_params = c(.fc_fs_params,
                                  list(outcome_type = "binary")))
  expect_identical(res_b$call_args$otype, "binary")
  expect_equal(res_b$d_eff, 15)

  # outcome_type_approx takes precedence over fs_params$outcome_type
  # (roxygen: fs_params$outcome_type is used only when outcome_type_approx is
  # NULL).  Continuous with sigma_y = 1.5: by hand from d_eff_continuous()'s
  # documentation, n_H / sigma_y^2 = 60 / 2.25 = 80 / 3 (binary would give
  # 15, survival at the default prop_cens = 0 would give 60).
  res_c <- .fc_call(fs_params = c(.fc_fs_params,
                                  list(outcome_type = "binary")),
                    outcome_type_approx = "continuous", sigma_y = 1.5)
  expect_identical(res_c$call_args$otype, "continuous")
  expect_equal(res_c$d_eff, 80 / 3)

  # d_eff_approx, when supplied, is used as given.
  res_d <- .fc_call(fs_params = .fc_fs_params, d_eff_approx = 30)
  expect_equal(res_d$d_eff, 30)
})
