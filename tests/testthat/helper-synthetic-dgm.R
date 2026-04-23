# ============================================================================
# helper-synthetic-dgm.R
#
# Shared synthetic-data factories for forestsearch test files.  These build
# small (N = 100-300) deterministic datasets exercising each supported
# outcome type with a known interacting subgroup.  The goal is correctness
# of the package's error-handling and return-shape contracts, NOT
# simulation precision -- hence the small N.
#
# All factories take a seed argument (default fixed) for deterministic
# reproduction and return plain data.frames with canonical column names:
#   id, treat, y, event (where applicable), time (where applicable)
#   plus covariates depending on the factory.
# ============================================================================


#' Synthetic survival dataset with a known harm subgroup.
#'
#' Hazard ratio HR_harm > 1 within a subgroup defined by
#' age <= 50 & stage == "3".  Outcome is a Weibull time with
#' administrative censoring at tau.
#'
#' @keywords internal
.make_survival_data <- function(N = 200L, HR_harm = 2.0, tau = 60,
                                seed = 42L) {
  set.seed(seed)
  age   <- round(rnorm(N, 55, 12))
  stage <- factor(sample(c("1", "2", "3"), N, replace = TRUE),
                  levels = c("1", "2", "3"))
  sex   <- factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L))
  noise <- rnorm(N)
  treat <- rbinom(N, 1L, 0.5)

  in_harm <- age <= 50 & stage == "3"
  lambda  <- 0.05 * exp(ifelse(in_harm, log(HR_harm), 0) * treat)
  time_true <- rexp(N, rate = lambda)
  time_obs  <- pmin(time_true, tau)
  event     <- as.integer(time_true <= tau)

  data.frame(
    id     = seq_len(N),
    treat  = treat,
    time   = time_obs,
    event  = event,
    age    = age,
    stage  = stage,
    sex    = sex,
    noise  = noise,
    stringsAsFactors = FALSE
  )
}


#' Synthetic binary-outcome dataset with a known harm subgroup.
#'
#' Odds ratio OR_harm > 1 within a subgroup defined by
#' age <= 50 & biomarker_hi == "1".  Outcome is Bernoulli.
#' Covariates are a MIX of continuous and factor (the canonical
#' pattern for GLM GRF).
#'
#' @keywords internal
.make_binary_data <- function(N = 200L, OR_harm = 2.5, seed = 42L) {
  set.seed(seed)
  df <- data.frame(
    id           = seq_len(N),
    treat        = rbinom(N, 1L, 0.5),
    age          = round(rnorm(N, 55, 12)),
    wtkg         = round(rnorm(N, 75, 12), 1),
    biomarker    = rnorm(N),
    biomarker_hi = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
    sex          = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
    stringsAsFactors = FALSE
  )
  in_harm <- df$age <= 50 & df$biomarker_hi == "1"
  lin     <- df$treat * ifelse(in_harm, log(OR_harm), 0)
  df$y    <- rbinom(N, 1L, plogis(lin))
  df
}


#' Synthetic continuous-outcome dataset.
#'
#' Mean difference MD_harm > 0 within the subgroup.
#'
#' @keywords internal
.make_continuous_data <- function(N = 200L, MD_harm = 1.5, seed = 42L) {
  set.seed(seed)
  df <- data.frame(
    id           = seq_len(N),
    treat        = rbinom(N, 1L, 0.5),
    age          = round(rnorm(N, 55, 12)),
    biomarker    = rnorm(N),
    biomarker_hi = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
    sex          = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
    stringsAsFactors = FALSE
  )
  in_harm <- df$age <= 50 & df$biomarker_hi == "1"
  df$y    <- rnorm(N, mean = df$treat * ifelse(in_harm, MD_harm, 0), sd = 1)
  df
}


#' Synthetic count-outcome dataset with a follow-up time offset.
#'
#' Incidence rate ratio IRR_harm > 1 within the subgroup.
#'
#' @keywords internal
.make_count_data <- function(N = 200L, IRR_harm = 2.0, seed = 42L) {
  set.seed(seed)
  df <- data.frame(
    id           = seq_len(N),
    treat        = rbinom(N, 1L, 0.5),
    age          = round(rnorm(N, 55, 12)),
    biomarker    = rnorm(N),
    biomarker_hi = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
    sex          = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
    ftime        = round(runif(N, 1, 5), 2),
    stringsAsFactors = FALSE
  )
  in_harm <- df$age <= 50 & df$biomarker_hi == "1"
  rate    <- 0.5 * exp(df$treat * ifelse(in_harm, log(IRR_harm), 0))
  df$y    <- rpois(N, rate * df$ftime)
  df
}


#' Canonical forestsearch() argument list for a given outcome type.
#'
#' Returns a list of argument values suitable for do.call(forestsearch, ...).
#' Uses tight thresholds and small fs.splits for fast test execution.
#'
#' @keywords internal
.fs_args_for <- function(outcome_type = c("survival", "binary",
                                          "continuous", "count"),
                         confounders = NULL,
                         extra = list()) {
  outcome_type <- match.arg(outcome_type)

  base <- list(
    outcome_type           = outcome_type,
    sg_focus               = "maxSG",
    hr.threshold           = 1.25,
    hr.consistency         = 1.00,
    pconsistency.threshold = 0.80,
    max_subgroups_search   = 5L,
    use_grf                = TRUE,
    use_lasso              = FALSE,
    is.RCT                 = TRUE,
    vi.grf.min             = -0.2,
    cont.cutoff            = 4L,
    cut_type               = "default",
    seedit                 = 42L,
    maxk                   = 2L,
    n.min                  = 40L,
    d0.min                 = 5L,
    d1.min                 = 5L,
    fs.splits              = 100L,
    details                = FALSE,
    quiet                  = TRUE,
    plot.sg                = FALSE,
    plot.grf               = FALSE,
    parallel_args          = list(plan    = "sequential",
                                  workers = 1L,
                                  show_message = FALSE)
  )

  # Outcome-specific additions
  type_args <- switch(outcome_type,
    survival = list(
      outcome.name   = "time",
      event.name     = "event",
      treat.name     = "treat",
      id.name        = "id"
    ),
    binary = list(
      outcome.name   = "y",
      event.name     = "y",
      treat.name     = "treat",
      id.name        = "id",
      effect_measure = "OR"
    ),
    continuous = list(
      outcome.name   = "y",
      treat.name     = "treat",
      id.name        = "id",
      effect_measure = "MD"
    ),
    count = list(
      outcome.name   = "y",
      treat.name     = "treat",
      id.name        = "id",
      effect_measure = "IRR",
      offset.name    = "ftime"
    )
  )

  if (!is.null(confounders)) {
    base$confounders.name <- confounders
  }

  modifyList(modifyList(base, type_args), extra)
}


#' Run forestsearch() quietly and capture warnings as a character vector
#' alongside the return value.  Useful when a test needs to distinguish
#' "succeeded with warnings" from "failed silently".
#'
#' @keywords internal
.run_fs_capture <- function(df, args) {
  warnings_seen <- character(0)
  val <- withCallingHandlers(
    do.call(forestsearch, c(list(df.analysis = df), args)),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = warnings_seen)
}


#' Assertion helper: a forestsearch() return object satisfies the
#' downstream CV/bootstrap consumer contract.
#'
#' The actual contract enforced by forestsearch_tenfold() and
#' forestsearch_bootstrap_dofuture() is a PRESENCE check on the names
#' of the returned list, not a non-null check.  In particular:
#'
#'   - \code{args_call_all}: must be present AND non-null (required to
#'     re-run forestsearch() on each CV fold / bootstrap replicate)
#'   - \code{sg.harm}: slot must exist (NULL is a valid value meaning
#'     "no subgroup identified")
#'   - \code{df.est}: slot must exist (NULL when no subgroup identified,
#'     populated when a subgroup is found)
#'
#' This helper encodes that exact contract.  For tests that require a
#' subgroup to be identified (e.g. CV/bootstrap integration tests),
#' use \code{.expect_fs_identified_return()} instead.
#'
#' @keywords internal
.expect_fs_full_return <- function(fs) {
  expect_true(is.list(fs))
  expect_true(!is.null(fs$args_call_all),
              info = "args_call_all must be non-null for downstream consumers")
  expect_true("sg.harm" %in% names(fs),
              info = "sg.harm slot must exist (may be NULL)")
  expect_true("df.est" %in% names(fs),
              info = "df.est slot must exist (may be NULL)")
}


#' Assertion helper: a forestsearch() run that SHOULD have identified a
#' subgroup (and therefore should have a populated df.est).  Use this
#' for tests where the DGM is strong enough that subgroup identification
#' is an expected outcome, not an aspirational one.
#'
#' @keywords internal
.expect_fs_identified_return <- function(fs) {
  .expect_fs_full_return(fs)
  expect_true(!is.null(fs$sg.harm),
              info = "expected a subgroup to be identified")
  expect_true(!is.null(fs$df.est),
              info = "df.est must be populated when a subgroup is identified")
}


#' Assertion helper: a forestsearch_tenfold() return object has the
#' expected shape.  Used across Tier 3/4 tests that exercise the CV
#' pipeline.
#'
#' @keywords internal
.expect_cv_return_shape <- function(cv) {
  expect_true(is.list(cv))
  expect_s3_class(cv, "fs_tenfold")
  # Core slots: must all be present
  required_slots <- c("sens_summary", "find_summary", "fold_summary",
                      "sims", "Kfolds")
  missing_slots <- setdiff(required_slots, names(cv))
  expect_identical(missing_slots, character(0),
    info = sprintf("CV return missing slots: %s",
                   paste(missing_slots, collapse = ", ")))
  # fold_summary is a data.frame with required columns
  fsum <- cv$fold_summary
  expect_true(is.data.frame(fsum))
  fsum_required <- c("sim", "fold", "sg1", "sg2", "any_found",
                     "grf_cuts", "pconsistency", "training_fs_hr",
                     "n_candidates_evaluated")
  missing_cols <- setdiff(fsum_required, names(fsum))
  expect_identical(missing_cols, character(0),
    info = sprintf("fold_summary missing columns: %s",
                   paste(missing_cols, collapse = ", ")))
}


#' Assertion helper: a forestsearch() early-return object at least
#' carries args_call_all so downstream CV/bootstrap consumers can
#' validate gracefully.
#'
#' @keywords internal
.expect_fs_early_return_shape <- function(fs, expected_stage = NULL) {
  expect_true(is.list(fs))
  expect_true(!is.null(fs$args_call_all),
              info = "early-return object must carry args_call_all")
  expect_null(fs$sg.harm,
              info = "early-return object has sg.harm = NULL")
  expect_true(!is.null(fs$error_log),
              info = "early-return object has error_log")
  if (!is.null(expected_stage)) {
    expect_identical(fs$error_log$stage, expected_stage)
  }
}
