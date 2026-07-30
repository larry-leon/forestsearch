# =============================================================================
# Post-selection inference on a completed forestsearch analysis
# (multiplier resampling: the `mr_inference` surface)
# =============================================================================
# First committed against the pre-rename identifiers, so that the FB/MR
# terminology rename could be shown to be behaviour-preserving: across that
# rename only the names in this file changed, and the assertion counts did not.
#
# What is pinned here:
#   1. the surface is OFF by default, and a default run returns NULL
#   2. enabled on the Cox/consistency path it returns the documented fields
#   3. no subgroup identified -> clean skip, harm flag NA, no error
#   4. the two harm-confirmation rules ("point"/"ci") both run and can differ
#   5. the re-selection rule derived from sg_focus, for all eleven accepted foci
#
# Item 6 pins a KNOWN DEFECT as present behaviour, not as desired behaviour:
# on the default consistency path a GLM outcome silently skips the whole step.
# See the note on that test.

CONF_SURV <- c("age", "stage", "sex", "noise")

# One shared fit per configuration: forestsearch() runs ~3 s here, so the
# MR-enabled fits are built once at file scope rather than per expectation.
.mri_args <- function(...) {
  .fs_args_for("survival", confounders = CONF_SURV,
               extra = utils::modifyList(
                 list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE),
                 list(...)))
}

.mri_draws <- list(draws = 200L, seed = 11L)

.mri_data <- .make_survival_data(N = 400L, HR_harm = 3.0)


# ---------------------------------------------------------------------------
# 1. Off by default
# ---------------------------------------------------------------------------

test_that("the de-biasing step is off by default and returns NULL when off", {
  expect_identical(formals(forestsearch)$mr_inference, FALSE)
  expect_identical(formals(forestsearch)$mr_inference_args, quote(list()))

  fs <- .run_fs_capture(.mri_data, .mri_args())$value

  # The slot must EXIST (downstream consumers index it) and be NULL.
  expect_true("mr_inference" %in% names(fs))
  expect_null(fs$mr_inference)
  expect_true("mr_harm_confirmed" %in% names(fs))
  expect_true(is.na(fs$mr_harm_confirmed))

  # A default analysis is pure identification: the subgroup is still found.
  .expect_fs_identified_return(fs)
})


# ---------------------------------------------------------------------------
# 2. Enabled on the Cox/consistency path: the documented fields
# ---------------------------------------------------------------------------

.mri_fit <- local({
  r <- .run_fs_capture(.mri_data,
                       .mri_args(mr_inference = TRUE,
                                 mr_inference_args = .mri_draws))
  r
})

test_that("enabled on the Cox/consistency path it returns the documented fields", {
  fs <- .mri_fit$value
  .expect_fs_identified_return(fs)

  g <- fs$mr_inference
  expect_false(is.null(g))
  expect_true(is.list(g))

  # Top-level contract.
  expect_true(all(c("selected_index", "selected_label", "measure", "log_scale",
                    "ci_method", "naive", "debiased", "selection_bias",
                    "fixed_bias", "selection_rate", "complement", "settings",
                    "harm_flag", "n_family", "n_selected",
                    "timing_seconds") %in% names(g)))

  expect_identical(g$measure, "HR")
  expect_true(g$log_scale)
  expect_identical(g$ci_method, "ij")
  expect_true(is.finite(g$selection_bias))
  expect_true(is.finite(g$fixed_bias))
  expect_true(g$selection_rate >= 0 && g$selection_rate <= 1)
  expect_gte(g$n_family, 1L)
  expect_gte(g$n_selected, 1L)

  # Effect-scale estimates with 95% CIs (naive and de-biased).
  for (el in list(g$naive, g$debiased)) {
    expect_true(all(c("est", "lower", "upper") %in% names(el)))
    expect_true(is.finite(el$est))
    expect_lte(el$lower, el$est)
    expect_gte(el$upper, el$est)
  }

  # De-biased element carries the SE provenance the docs promise.
  expect_true(all(c("lower_1s", "se", "se_ij", "se_wald", "var_ij",
                    "ij_source", "ij_draws") %in% names(g$debiased)))
  expect_true(g$debiased$ij_source %in% c("ij", "ij_raw", "wald_fallback"))
  expect_true(is.finite(g$debiased$se) && g$debiased$se > 0)
  # Default ci_method = "ij" means the reported SE is the IJ one.
  expect_equal(g$debiased$se, g$debiased$se_ij)
  # One-sided 95% bound sits above the two-sided 95% lower bound.
  expect_gte(g$debiased$lower_1s, g$debiased$lower)

  # Harm-confirmation metadata.
  expect_true(all(c("t_confirm", "confirm_rule", "reselection",
                    "selection_rule", "multiplier", "draws")
                  %in% names(g$settings)))
  expect_identical(g$settings$confirm_rule, "point")
  expect_identical(g$settings$draws, 200L)
  expect_identical(g$settings$multiplier, "poisson")
  # Ratio measure -> the near-null default threshold is 1.
  expect_equal(g$settings$t_confirm, 1)

  # The flag is a bare logical, and the top-level mirror agrees with it.
  expect_true(is.logical(g$harm_flag))
  expect_length(g$harm_flag, 1L)
  expect_false(is.na(g$harm_flag))
  expect_identical(fs$mr_harm_confirmed, g$harm_flag)
})

test_that("the near-null default threshold is measure-dependent", {
  expect_equal(.fs_mr_confirm_null("HR"), 1)
  expect_equal(.fs_mr_confirm_null("OR"), 1)
  expect_equal(.fs_mr_confirm_null("RR"), 1)
  expect_equal(.fs_mr_confirm_null("IRR"), 1)
  expect_equal(.fs_mr_confirm_null("RD"), 0)
  expect_equal(.fs_mr_confirm_null("MD"), 0)
  expect_equal(.fs_mr_confirm_null(NULL), 1)
})

test_that("the complement is de-biased by default from forestsearch()", {
  # include_complement defaults to TRUE at every internal caller, so the
  # complement block is populated in practice even though fs_mr_inference()'s
  # own formal default is FALSE.
  cc <- .mri_fit$value$mr_inference$complement
  expect_false(is.null(cc))
  expect_true(all(c("naive", "debiased", "selection_bias", "fixed_bias", "n")
                  %in% names(cc)))
  expect_true(is.finite(cc$debiased$est))
})

test_that("enabling the de-biasing step does not alter identification", {
  # The step runs strictly after sg.harm is assigned, so it cannot feed back
  # into selection. Same data, same seed, MR on vs off.
  off <- .run_fs_capture(.mri_data, .mri_args())$value
  on  <- .mri_fit$value
  expect_identical(off$sg.harm, on$sg.harm)
  expect_identical(dim(off$df.est), dim(on$df.est))
  expect_identical(names(off), names(on))
})


# ---------------------------------------------------------------------------
# 3. No subgroup identified -> clean skip
# ---------------------------------------------------------------------------

test_that("no subgroup identified gives a clean skip, NA flag, and no error", {
  d_null <- .make_survival_data(N = 300L, HR_harm = 1.0)
  r <- .run_fs_capture(
    d_null,
    .mri_args(hr.threshold = 3.0, pconsistency.threshold = 0.99,
              mr_inference = TRUE, mr_inference_args = .mri_draws))
  fs <- r$value

  expect_null(fs$sg.harm)
  expect_null(fs$mr_inference)
  expect_true(is.na(fs$mr_harm_confirmed))
  # Clean skip: nothing about the de-biasing step is warned about.
  expect_false(any(grepl("mr_inference|harm confirmation",r$warnings,
                         ignore.case = TRUE)))
})


# ---------------------------------------------------------------------------
# 4. The two harm-confirmation rules
# ---------------------------------------------------------------------------

test_that("both harm-confirmation rules run and can disagree", {
  # "ci" is the strictly stronger requirement: it asks the one-sided 95%
  # selection-adjusted lower bound to clear the threshold, where "point" asks
  # only the point estimate to. Choose a threshold between the two so the rules
  # are forced to disagree, rather than asserting a disagreement that depends
  # on the DGM.
  g0 <- .mri_fit$value$mr_inference
  skip_if(is.null(g0), "de-biasing step did not run on the reference fit")

  est <- g0$debiased$est
  lo1 <- g0$debiased$lower_1s
  skip_if_not(is.finite(est) && is.finite(lo1) && lo1 < est,
              "reference fit has no gap between the point estimate and its bound")
  t_mid <- (est + lo1) / 2

  mk <- function(rule) {
    .run_fs_capture(
      .mri_data,
      .mri_args(mr_inference = TRUE,
                mr_inference_args = utils::modifyList(
                  .mri_draws,
                  list(confirm_rule = rule,
                       t_confirm = t_mid))))$value$mr_inference
  }
  gp <- mk("point")
  gc <- mk("ci")

  expect_false(is.null(gp))
  expect_false(is.null(gc))
  expect_identical(gp$settings$confirm_rule, "point")
  expect_identical(gc$settings$confirm_rule, "ci")
  expect_equal(gp$settings$t_confirm, t_mid)
  expect_equal(gc$settings$t_confirm, t_mid)

  # Same draws (same seed) -> identical estimates; only the decision differs.
  expect_equal(gp$debiased$est, gc$debiased$est)
  expect_true(gp$harm_flag)
  expect_false(gc$harm_flag)
})


# ---------------------------------------------------------------------------
# 5. The re-selection rule derived from sg_focus
# ---------------------------------------------------------------------------

test_that("sg_focus maps to the re-selection rule per engine, for all eleven foci", {
  # The full mapping table. The consistency engine ranks by consistency rate,
  # so the "hr" family maps to "maxcons" there and to "maxeff" for the
  # effect-ranking engines (DINA/GRF); every other focus maps identically.
  map <- data.frame(
    sg_focus    = c("hr",      "eff",     "maxcons",
                    "maxeff",  "maxeffCons",
                    "maxSG",   "minSG",
                    "hrMaxSG", "hrMinSG",
                    "effMaxSG", "effMinSG"),
    consistency = c("maxcons", "maxcons", "maxcons",
                    "maxeff",  "maxeff",
                    "maxSG",   "minSG",
                    "effMaxSG", "effMinSG",
                    "effMaxSG", "effMinSG"),
    effect      = c("maxeff",  "maxeff",  "maxeff",
                    "maxeff",  "maxeff",
                    "maxSG",   "minSG",
                    "effMaxSG", "effMinSG",
                    "effMaxSG", "effMinSG"),
    stringsAsFactors = FALSE)

  for (i in seq_len(nrow(map))) {
    expect_identical(
      .fs_mr_reselection_from_focus(map$sg_focus[i], engine = "consistency"),
      map$consistency[i],
      info = sprintf("sg_focus '%s', engine 'consistency'", map$sg_focus[i]))
    expect_identical(
      .fs_mr_reselection_from_focus(map$sg_focus[i], engine = "effect"),
      map$effect[i],
      info = sprintf("sg_focus '%s', engine 'effect'", map$sg_focus[i]))
  }

  # Every rule the mapping can produce must be one fs_mr_inference() accepts.
  accepted <- eval(formals(fs_mr_inference)$reselection)
  expect_true(all(unique(c(map$consistency, map$effect)) %in% accepted))
})

test_that("the rule reaching fs_mr_inference() is the one sg_focus implies", {
  # End-to-end, not just the mapping helper: what the search recorded.
  expect_identical(.mri_fit$value$mr_inference$settings$reselection,
                   .fs_mr_reselection_from_focus("maxSG", engine = "consistency"))

  r <- .run_fs_capture(
    .mri_data,
    .mri_args(sg_focus = "maxeffCons", mr_inference = TRUE,
              mr_inference_args = .mri_draws))$value
  skip_if(is.null(r$mr_inference), "de-biasing step did not run under maxeffCons")
  expect_identical(r$mr_inference$settings$reselection, "maxeff")
})


# ---------------------------------------------------------------------------
# 6. KNOWN DEFECT, pinned as present behaviour
# ---------------------------------------------------------------------------

test_that("GLM outcomes silently skip the step under the default consistency_method", {
  # PRESENT BEHAVIOUR, NOT DESIRED BEHAVIOUR. The consistency branch runs only
  # when `consistency_method == "resample" && !is.null(estimator_fn)` (GLM) or
  # `outcome_type == "survival" && is.null(estimator_fn)` (Cox). For a GLM
  # outcome under the default consistency_method = "split", neither holds, so
  # the user opts in, nothing runs, and nothing warns -- leaving the harm flag
  # NA, which isTRUE() renders as "harm not confirmed".
  #
  # This test exists so the rename can be shown to preserve behaviour. Fixing
  # the defect is separate work; when it is fixed, this test SHOULD fail and
  # should then be rewritten to assert the fixed behaviour.
  db <- .make_binary_data(N = 400L, OR_harm = 3.0)
  bargs <- function(cm) {
    .fs_args_for("binary",
                 confounders = c("age", "biomarker_hi", "biomarker", "sex"),
                 extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE,
                              consistency_method = cm,
                              mr_inference = TRUE,
                              mr_inference_args = .mri_draws))
  }

  r_split <- .run_fs_capture(db, bargs("split"))
  expect_false(is.null(r_split$value$sg.harm))       # a subgroup WAS found
  expect_null(r_split$value$mr_inference)             # ... yet nothing ran
  expect_true(is.na(r_split$value$mr_harm_confirmed))
  expect_false(any(grepl("mr_inference|harm confirmation",r_split$warnings,
                         ignore.case = TRUE)))       # ... and nothing warned

  # The contrast: on the resample path the same data does run the step.
  r_res <- .run_fs_capture(db, bargs("resample"))
  expect_false(is.null(r_res$value$mr_inference))
  expect_false(is.na(r_res$value$mr_harm_confirmed))
})
