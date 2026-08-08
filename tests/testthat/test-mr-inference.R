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

test_that("GLM outcomes silently skip the step under consistency_method = 'split'", {
  # PRESENT BEHAVIOUR, NOT DESIRED BEHAVIOUR. The consistency branch runs only
  # when `consistency_method == "resample" && !is.null(estimator_fn)` (GLM) or
  # `outcome_type == "survival" && is.null(estimator_fn)` (Cox). For a GLM
  # outcome under consistency_method = "split", neither holds, so the user opts
  # in, nothing runs, and nothing warns -- leaving the harm flag NA, which
  # isTRUE() renders as "harm not confirmed".  Both arms are passed explicitly
  # below, so this test does not depend on which value is the formal default.
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


# ---------------------------------------------------------------------------
# 7. mr_in_replicates: MR inside bootstrap replicates and CV folds
# ---------------------------------------------------------------------------
# Before this control existed, `mr_inference` propagated into every replicate
# and fold via args_call_all's bulk mget(), and the resulting object was never
# consumed by anything.  The default now strips it; TRUE runs it AND retains
# the result.  The two behaviours are deliberately tied to one argument so
# they can never diverge into "computed and discarded".

# Counts calls to the (internal) MR entry point while `expr` runs.  trace() on
# the namespace is visible to a sequential plan because the replicates then run
# in this same process; every test below pins plan = "sequential" for that
# reason.
.mri_count_mr <- function(expr) {
  n <- 0L
  suppressMessages(trace("fs_mr_inference", where = asNamespace("forestsearch"),
                         tracer = function() n <<- n + 1L, print = FALSE))
  on.exit(suppressMessages(untrace("fs_mr_inference",
                                   where = asNamespace("forestsearch"))),
          add = TRUE)
  force(expr)
  n
}

.mri_seq <- list(plan = "sequential", workers = 1L, show_message = FALSE)

test_that("mr_in_replicates defaults to FALSE on both resampling entry points", {
  expect_identical(formals(forestsearch_bootstrap_dofuture)$mr_in_replicates,
                   FALSE)
  expect_identical(formals(forestsearch_Kfold)$mr_in_replicates, FALSE)
})

test_that("default runs zero MR inside bootstrap replicates; top level keeps its own", {
  fs <- .run_fs_capture(.mri_data,
                        .mri_args(mr_inference = TRUE,
                                  mr_inference_args = .mri_draws))$value
  skip_if(is.null(fs$mr_inference), "MR did not run on the primary fit")

  # The top-level analysis DID get its de-biased estimate ...
  expect_false(is.null(fs$mr_inference))
  expect_false(is.na(fs$mr_harm_confirmed))

  # ... and the replicates run none at all.
  n <- .mri_count_mr(
    bc <- suppressWarnings(forestsearch_bootstrap_dofuture(
      fs.est = fs, nb_boots = 2L, parallel_args = .mri_seq)))
  expect_identical(n, 0L)
  expect_null(bc$mr_replicates)
})

test_that("default runs zero MR inside CV folds", {
  fs <- .run_fs_capture(.mri_data,
                        .mri_args(mr_inference = TRUE,
                                  mr_inference_args = .mri_draws))$value
  skip_if(is.null(fs$sg.harm), "no subgroup identified on the primary fit")

  n <- .mri_count_mr(
    cv <- suppressWarnings(forestsearch_Kfold(
      fs.est = fs, Kfolds = 2L, parallel_args = .mri_seq)))
  expect_identical(n, 0L)
  expect_null(cv$mr_replicates)
})

test_that("mr_in_replicates = TRUE retains one MR result per replicate", {
  fs <- .run_fs_capture(.mri_data,
                        .mri_args(mr_inference = TRUE,
                                  mr_inference_args = .mri_draws))$value
  skip_if(is.null(fs$mr_inference), "MR did not run on the primary fit")

  nb <- 2L
  n <- .mri_count_mr(
    bc <- suppressWarnings(suppressMessages(forestsearch_bootstrap_dofuture(
      fs.est = fs, nb_boots = nb, parallel_args = .mri_seq,
      mr_in_replicates = TRUE))))

  # MR runs once per IDENTIFYING replicate, not once per replicate: a
  # replicate that finds no subgroup has nothing to de-bias and correctly runs
  # none.  Asserting nb here would pass or fail on whether the fixture happened
  # to identify in every resample.  Tie the three counts together instead --
  # calls made, objects retained, replicates that identified -- which is the
  # invariant that actually holds and is stronger than a bare count.
  n_found <- sum(bc$results$any_found)
  expect_identical(n, as.integer(n_found))
  expect_lte(n, nb)

  expect_false(is.null(bc$mr_replicates))      # ... and was RETAINED
  expect_length(bc$mr_replicates, nb)          # one slot per replicate
  expect_true(is.list(bc$mr_replicates))

  keep <- Filter(Negate(is.null), bc$mr_replicates)
  expect_identical(length(keep), as.integer(n_found))

  # Documented structure: each element is that replicate's own MR object.
  skip_if(length(keep) == 0L, "no replicate produced an MR object")
  expect_true(all(vapply(keep, function(x) "debiased" %in% names(x), logical(1))))
})


# ---------------------------------------------------------------------------
# 8. Messaging (Part C): audible when it matters, silent under quiet
# ---------------------------------------------------------------------------

test_that("C.1 announces MR once under quiet = FALSE and never under quiet = TRUE", {
  loud <- testthat::capture_messages(
    suppressWarnings(do.call(forestsearch, c(
      list(df.analysis = .mri_data),
      .mri_args(mr_inference = TRUE, mr_inference_args = .mri_draws,
                quiet = FALSE)))))
  hits <- grep("Multiplier resampling \\(MR\\): constructing", loud)
  expect_length(hits, 1L)
  expect_match(loud[hits[1]], "200 draws")
  expect_match(loud[hits[1]], "cannot change the identified subgroup")

  quiet <- testthat::capture_messages(
    suppressWarnings(do.call(forestsearch, c(
      list(df.analysis = .mri_data),
      .mri_args(mr_inference = TRUE, mr_inference_args = .mri_draws,
                quiet = TRUE)))))
  expect_length(grep("Multiplier resampling", quiet), 0L)
})

test_that("C.2 names consistency_method = resample for a GLM outcome under split", {
  db <- .make_binary_data(N = 400L, OR_harm = 3.0)
  bargs <- function(quiet) {
    .fs_args_for("binary",
                 confounders = c("age", "biomarker_hi", "biomarker", "sex"),
                 extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = quiet,
                              consistency_method = "split",
                              mr_inference = TRUE,
                              mr_inference_args = .mri_draws))
  }

  msgs <- testthat::capture_messages(
    fs <- suppressWarnings(do.call(forestsearch,
                                   c(list(df.analysis = db), bargs(FALSE)))))
  hits <- grep("was NOT performed", msgs)
  expect_length(hits, 1L)
  expect_match(msgs[hits[1]], "consistency_method = \"resample\"", fixed = TRUE)
  expect_match(msgs[hits[1]], "NA is not evidence against harm", fixed = TRUE)

  # BEHAVIOUR UNCHANGED: the skip itself, and the NA, are exactly as before.
  expect_null(fs$mr_inference)
  expect_true(is.na(fs$mr_harm_confirmed))

  # ... and the same call says nothing under quiet = TRUE.
  qmsgs <- testthat::capture_messages(
    suppressWarnings(do.call(forestsearch,
                             c(list(df.analysis = db), bargs(TRUE)))))
  expect_length(grep("was NOT performed", qmsgs), 0L)
})

test_that("C.3 fires exactly once for the whole loop, not once per replicate", {
  fs <- .run_fs_capture(.mri_data,
                        .mri_args(mr_inference = TRUE,
                                  mr_inference_args = .mri_draws,
                                  quiet = FALSE))$value
  skip_if(is.null(fs$sg.harm), "no subgroup identified on the primary fit")

  nb <- 8L
  msgs <- testthat::capture_messages(
    suppressWarnings(forestsearch_bootstrap_dofuture(
      fs.est = fs, nb_boots = nb, parallel_args = .mri_seq,
      mr_in_replicates = TRUE)))
  expect_length(grep("will run inside each of the", msgs), 1L)

  # And nothing at all under the default.
  msgs0 <- testthat::capture_messages(
    suppressWarnings(forestsearch_bootstrap_dofuture(
      fs.est = fs, nb_boots = nb, parallel_args = .mri_seq)))
  expect_length(grep("will run inside each of the", msgs0), 0L)
})


# ---------------------------------------------------------------------------
# 9. mr_in_replicates on the repeated K-fold consumer
# ---------------------------------------------------------------------------
# forestsearch_tenfold() rebuilds the same cv_args from args_call_all, so it
# carried the same discarded-MR behaviour -- at sims * Kfolds runs rather than
# Kfolds, since each simulation runs its own full fold loop.

test_that("mr_in_replicates defaults to FALSE on forestsearch_tenfold", {
  expect_identical(formals(forestsearch_tenfold)$mr_in_replicates, FALSE)
})

test_that("default runs zero MR inside tenfold's sim x fold grid", {
  fs <- .run_fs_capture(.mri_data,
                        .mri_args(mr_inference = TRUE,
                                  mr_inference_args = .mri_draws))$value
  skip_if(is.null(fs$sg.harm), "no subgroup identified on the primary fit")

  n <- .mri_count_mr(
    tf <- suppressWarnings(forestsearch_tenfold(
      fs.est = fs, sims = 2L, Kfolds = 2L, details = FALSE,
      parallel_args = .mri_seq)))
  expect_identical(n, 0L)
  expect_null(tf$mr_replicates)
})

test_that("tenfold retains one MR result per fold, nested by simulation", {
  fs <- .run_fs_capture(.mri_data,
                        .mri_args(mr_inference = TRUE,
                                  mr_inference_args = .mri_draws))$value
  skip_if(is.null(fs$mr_inference), "MR did not run on the primary fit")

  sims <- 2L; K <- 2L
  n <- .mri_count_mr(
    tf <- suppressWarnings(suppressMessages(forestsearch_tenfold(
      fs.est = fs, sims = sims, Kfolds = K, details = FALSE,
      parallel_args = .mri_seq, mr_in_replicates = TRUE))))

  # As in the bootstrap case: MR runs once per IDENTIFYING fold, not on every
  # sim x fold cell.  A fold that finds no subgroup has nothing to de-bias.
  # fold_summary$any_found is the per-cell identification flag, so it is the
  # correct denominator; sims * K is only an upper bound.
  n_found <- sum(tf$fold_summary$any_found)
  expect_identical(n, as.integer(n_found))
  expect_lte(n, sims * K)

  expect_false(is.null(tf$mr_replicates))       # ... and was RETAINED
  expect_length(tf$mr_replicates, sims)         # one slot per simulation ...
  expect_true(all(lengths(tf$mr_replicates) == K))   # ... each with K slots
  expect_match(names(tf$mr_replicates)[1], "^sim[0-9]+$")

  flat <- unlist(tf$mr_replicates, recursive = FALSE)
  expect_identical(sum(vapply(flat, Negate(is.null), logical(1))),
                   as.integer(n_found))

  keep <- Filter(Negate(is.null), flat)
  skip_if(length(keep) == 0L, "no fold produced an MR object")
  expect_true(all(vapply(keep, function(x) "debiased" %in% names(x), logical(1))))
})

test_that("tenfold C.3 fires exactly once for the whole sim x fold grid", {
  fs <- .run_fs_capture(.mri_data,
                        .mri_args(mr_inference = TRUE,
                                  mr_inference_args = .mri_draws,
                                  quiet = FALSE))$value
  skip_if(is.null(fs$sg.harm), "no subgroup identified on the primary fit")

  msgs <- testthat::capture_messages(
    suppressWarnings(forestsearch_tenfold(
      fs.est = fs, sims = 2L, Kfolds = 2L, details = FALSE,
      parallel_args = .mri_seq, mr_in_replicates = TRUE)))
  expect_length(grep("will run inside each of the", msgs), 1L)

  msgs0 <- testthat::capture_messages(
    suppressWarnings(forestsearch_tenfold(
      fs.est = fs, sims = 2L, Kfolds = 2L, details = FALSE,
      parallel_args = .mri_seq)))
  expect_length(grep("will run inside each of the", msgs0), 0L)
})
