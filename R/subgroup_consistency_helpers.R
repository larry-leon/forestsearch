# =============================================================================
# Subgroup Consistency Helper Functions
# =============================================================================
#
# Helper functions for subgroup consistency evaluation that need to be
# exported for parallel worker access. When using callr or multisession
# parallelization, workers only see exported functions from installed packages.
#
# =============================================================================

# =============================================================================
# PARALLEL SETUP
# =============================================================================

#' Set up parallel processing for subgroup consistency
#'
#' Sets up parallel processing using the specified approach and number of workers.
#'
#' @param parallel_args List with \code{plan} (character), \code{workers} (integer),
#'   and \code{show_message} (logical).
#'
#' @return None. Sets up parallel backend as side effect.
#'
#' @examples
#' \donttest{
#' setup_parallel_SGcons(list(plan = "sequential", workers = 1,
#'                             show_message = FALSE))
#' future::plan(future::sequential)  # reset
#' }
#' @importFrom future plan multisession multicore sequential
#' @export
setup_parallel_SGcons <- function(parallel_args = list(
    plan = "multisession",
    workers = 4,
    show_message = TRUE
)) {
  plan_type <- parallel_args$plan
  n_workers <- parallel_args$workers
  show_message <- isTRUE(parallel_args$show_message)

  if (is.null(plan_type)) stop("parallel_args$plan must be specified.")

  allowed_plans <- c("multisession", "multicore", "callr", "sequential")

  if (!plan_type %in% allowed_plans) {
    stop("plan_type must be one of: ", paste(allowed_plans, collapse = ", "))
  }

  max_cores <- parallel::detectCores()
  if (is.null(n_workers) || !is.numeric(n_workers) || n_workers < 1) {
    n_workers <- 1
  } else {
    n_workers <- min(n_workers, max_cores)
  }

  if (plan_type == "multisession") {
    future::plan(future::multisession, workers = n_workers)
    if (show_message) message("Parallel plan: multisession with ", n_workers, " workers.")
  } else if (plan_type == "multicore") {
    future::plan(future::multicore, workers = n_workers)
    if (show_message) message("Parallel plan: multicore with ", n_workers, " workers.")
  } else if (plan_type == "callr") {
    if (!requireNamespace("future.callr", quietly = TRUE)) {
      stop("The 'future.callr' package is required for the callr approach.")
    }
    future::plan(future.callr::callr, workers = n_workers)
    if (show_message) message("Parallel plan: callr with ", n_workers, " workers.")
  } else if (plan_type == "sequential") {
    future::plan(future::sequential)
    if (show_message) message("Sequential plan")
  } else {
    stop("Unknown parallel plan: ", plan_type)
  }
}


# =============================================================================
# STATISTICAL HELPER FUNCTIONS
# =============================================================================

#' Wilson Score Confidence Interval
#'
#' Computes Wilson score confidence interval for a proportion, which has
#' better coverage properties than the normal approximation for small samples
#' and proportions near 0 or 1.
#'
#' @param x Integer. Number of successes.
#' @param n Integer. Number of trials.
#' @param conf.level Numeric. Confidence level (default 0.95).
#'
#' @return Named numeric vector with elements: estimate, lower, upper.
#'
#' @examples
#' wilson_ci(90, 100)
#' wilson_ci(5, 20, conf.level = 0.90)
#' @export
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) {
    return(c(estimate = NA_real_, lower = 0, upper = 1))
  }

  z <- qnorm(1 - (1 - conf.level) / 2)
  p_hat <- x / n

  denominator <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denominator
  margin <- (z / denominator) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))

  lower <- max(0, center - margin)
  upper <- min(1, center + margin)

  c(estimate = p_hat, lower = lower, upper = upper)
}


#' Early Stopping Decision
#'
#' Evaluates whether enough evidence exists to stop early based on
#' confidence interval for consistency proportion.
#'
#' @param n_success Integer. Number of splits meeting consistency.
#' @param n_total Integer. Total number of valid splits.
#' @param threshold Numeric. Target consistency threshold.
#' @param conf.level Numeric. Confidence level for decision (default 0.95).
#' @param min_samples Integer. Minimum samples before allowing early stop.
#'
#' @return Character. One of "continue", "pass", or "fail".
#'
#' @examples
#' early_stop_decision(95, 100, threshold = 0.90)
#' early_stop_decision(60, 100, threshold = 0.90)
#' early_stop_decision(10, 15, threshold = 0.90)  # below min_samples
#' @export
early_stop_decision <- function(n_success, n_total, threshold,
                                conf.level = 0.95, min_samples = 20) {
  if (n_total < min_samples) {
    return("continue")
  }

  ci <- wilson_ci(n_success, n_total, conf.level)

  if (ci["lower"] >= threshold) {
    return("pass")
  }

  if (ci["upper"] < threshold) {
    return("fail")
  }

  return("continue")
}


# =============================================================================
# COX MODEL HELPERS
# =============================================================================

#' Fast Cox Model HR Estimation
#'
#' Fits a minimal Cox model to estimate hazard ratio with reduced overhead.
#' Disables robust variance, model matrix storage, and other extras for speed.
#'
#' @param df data.frame or data.table with Y, Event, Treat columns.
#' @param cox_init Numeric. Initial value for coefficient (default 0).
#'   Ignored when `adjust_covariates` is supplied (the default coxph
#'   starting values are used so the init length matches the number of
#'   model coefficients).
#' @param adjust_covariates Character vector or `NULL`. Additional terms
#'   appended to the model right-hand side after `Treat`.  Terms are pasted
#'   verbatim, so `"strata(x1)"` produces a stratified Cox model.  The
#'   referenced columns must be present in `df`.  `NULL` (default) fits the
#'   unadjusted `Surv(Y, Event) ~ Treat` model.
#'
#' @return Numeric. Estimated hazard ratio for `Treat`, or NA if model fails.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' df <- data.frame(
#'   Y     = rexp(80),
#'   Event = rbinom(80, 1, 0.6),
#'   Treat = rep(0:1, 40)
#' )
#' get_split_hr_fast(df)
#' }
#' @importFrom survival coxph Surv
#' @export
get_split_hr_fast <- function(df, cox_init = 0, adjust_covariates = NULL) {
  if (nrow(df) < 2 || sum(df$Event) < 2) {
    return(NA_real_)
  }

  adj <- .fs_adjust_terms(adjust_covariates)
  if (length(adj) > 0L) {
    fmla <- stats::reformulate(c("Treat", adj), response = "survival::Surv(Y, Event)")
  } else {
    fmla <- survival::Surv(Y, Event) ~ Treat
  }

  # coxph treats `init = NULL` as a length-0 init (an error), so the init
  # argument is omitted entirely on the adjusted path; the scalar warm-start
  # only applies to the single-coefficient treatment-only model.
  fit_args <- list(
    formula = fmla, data = df,
    robust = FALSE, model = FALSE, x = FALSE, y = FALSE
  )
  if (length(adj) == 0L) fit_args$init <- cox_init

  fit <- tryCatch(
    suppressWarnings(do.call(survival::coxph, fit_args)),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NA_real_)
  cf <- fit$coefficients
  beta <- if ("Treat" %in% names(cf)) cf[["Treat"]] else cf[1]
  return(exp(beta))
}


#' Run Single Consistency Split
#'
#' Performs one random 50/50 split and evaluates whether both halves
#' meet the HR consistency threshold.
#'
#' @param df.x data.table. Subgroup data with columns Y, Event, Treat.
#' @param N.x Integer. Number of observations in subgroup.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param cox_init Numeric. Initial value for Cox model (log HR).
#' @param estimator_fn Closure or \code{NULL}. Effect-estimator closure from
#'   \code{\link{make_effect_estimator}} for GLM outcomes.  Default \code{NULL}
#'   (uses Cox model).
#' @param consistency_threshold Numeric or \code{NULL}. Threshold for GLM
#'   consistency evaluation.  Default \code{NULL} (uses \code{hr.consistency}).
#' @param adjust_covariates Character vector or \code{NULL}. Additional Cox
#'   model terms (e.g. \code{"strata(x1)"}) passed to
#'   \code{\link{get_split_hr_fast}} on the survival path.  Ignored on the
#'   GLM path.  Default \code{NULL} (unadjusted).
#'
#' @return Numeric. 1 if both splits meet threshold, 0 if not, NA if error.
#'
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(1)
#' df <- data.table(
#'   Y     = rexp(100),
#'   Event = rbinom(100, 1, 0.55),
#'   Treat = rep(0:1, 50)
#' )
#' run_single_consistency_split(df, N.x = 100, hr.consistency = 1.0)
#' }
#' @export
run_single_consistency_split <- function(df.x, N.x, hr.consistency, cox_init = 0,
                                         estimator_fn = NULL,
                                         consistency_threshold = NULL,
                                         adjust_covariates = NULL) {

  in.split1 <- tryCatch({
    sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(in.split1)) return(NA_real_)

  df.x$insplit1 <- in.split1
  df.x.split1 <- df.x[df.x$insplit1 == TRUE, ]
  df.x.split2 <- df.x[df.x$insplit1 == FALSE, ]

  if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) {
    return(NA_real_)
  }

  # --- GLM path: use estimator closure ---
  if (!is.null(estimator_fn) && !is.null(consistency_threshold)) {
    # Minimum per-arm check (replaces event-count check for non-survival)
    for (half in list(df.x.split1, df.x.split2)) {
      treat_col <- intersect(c("Treat", "treat"), names(half))
      if (length(treat_col) > 0L) {
        tc <- treat_col[1L]
        if (sum(half[[tc]] == 1L) < 3L || sum(half[[tc]] == 0L) < 3L) {
          return(NA_real_)
        }
      }
    }

    res1 <- estimator_fn(df.x.split1)
    res2 <- estimator_fn(df.x.split2)

    if (is.na(res1$estimate) || is.na(res2$estimate)) return(NA_real_)

    return(as.numeric(
      res1$estimate >= consistency_threshold &&
        res2$estimate >= consistency_threshold
    ))
  }

  # --- Survival path (existing, unchanged) ---
  if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) {
    return(NA_real_)
  }

  hr.split1 <- get_split_hr_fast(df.x.split1, cox_init, adjust_covariates = adjust_covariates)
  hr.split2 <- get_split_hr_fast(df.x.split2, cox_init, adjust_covariates = adjust_covariates)

  if (!is.na(hr.split1) && !is.na(hr.split2)) {
    as.numeric(hr.split1 > hr.consistency && hr.split2 > hr.consistency)
  } else {
    NA_real_
  }
}


# =============================================================================
# RESAMPLE-PATH HELPERS (shared by single-stage and two-stage evaluators)
# =============================================================================

#' Can the GLM resampling approximation represent this effect measure?
#'
#' The resampling approximation requires the treatment effect to be a single
#' model coefficient with a well-defined influence function. That holds for
#' OR/RR (logistic / log-binomial), RD (identity-link binomial), MD (OLS), and
#' IRR (Poisson with offset). It does not hold for IRD (a delta-method rate
#' difference) or propensity-adjusted (IPTW / G-computation) effects, for which
#' the caller falls back to literal splitting.
#'
#' @param spec List; the `glm_resample_spec` threaded from [forestsearch()].
#' @return Logical scalar.
#' @keywords internal
.glm_resample_supported <- function(spec) {
  if (is.null(spec) || is.null(spec$effect_measure)) return(FALSE)
  spec$effect_measure %in% c("OR", "RR", "RD", "MD", "IRR")
}

#' Consistency proportion via literal splitting (single-stage)
#'
#' Runs `n.splits` Bernoulli splits through [run_single_consistency_split()]
#' and returns the rounded proportion of consistent splits, or `NA_real_` when
#' too few valid splits are obtained. Survival and GLM (`estimator_fn`) paths.
#'
#' @keywords internal
.consistency_via_splits <- function(df.x, N.x, n.splits, hr.consistency,
                                    cox_init, estimator_fn,
                                    consistency_threshold, adjust_covariates,
                                    pconsistency.digits, m = NA, details = FALSE) {
  flag.consistency <- numeric(n.splits)
  for (i in seq_len(n.splits)) {
    flag.consistency[i] <- run_single_consistency_split(
      df.x, N.x, hr.consistency, cox_init,
      estimator_fn = estimator_fn,
      consistency_threshold = consistency_threshold,
      adjust_covariates = adjust_covariates
    )
  }
  n_valid <- sum(!is.na(flag.consistency))
  if (n_valid < 10) {
    if (details) cat("Subgroup ", m, ": too few valid splits (", n_valid, ")\n", sep = "")
    return(NA_real_)
  }
  tryCatch(round(mean(flag.consistency, na.rm = TRUE), pconsistency.digits),
           error = function(e) NA_real_)
}


# =============================================================================
# LABEL HELPERS
# =============================================================================

#' Convert Factor Code to Label
#'
#' Converts q-indexed codes to human-readable labels using the confs_labels
#' mapping. Supports both full format ("q1.1", "q3.0") and short format
#' ("q1", "q3"). Handles vector input via recursion.
#'
#' @param Qsg Character. Factor code in format `"q<index>.<action>"` or
#'   `"q<index>"`. For the full format, action 0 = NOT (complement),
#'   action 1 = IN (member). Short format defaults to action 1.
#'   Can also be a character vector for vectorized input.
#' @param confs_labels Character vector. Labels for each factor, indexed by
#'   factor number.
#'
#' @return Character. Human-readable label wrapped in braces, e.g.,
#'   `"{age <= 50}"` or `"!{age <= 50}"` for complement. Returns the
#'   original code if no match is found.
#'
#' @examples
#' \dontrun{
#' labels <- c("age <= 50", "tumor size <= 20", "nodes <= 3")
#' FS_labels("q1.1", labels)
#' FS_labels("q1.0", labels)
#' FS_labels("q2", labels)
#' FS_labels(c("q1", "q3"), labels)
#' }
#' @export

FS_labels <- function(Qsg, confs_labels) {
  # Handle vector input
  if (length(Qsg) > 1) {
    return(vapply(Qsg, FS_labels, character(1), confs_labels = confs_labels,
                  USE.NAMES = FALSE))
  }

  # Try full format: q<index>.<action> (e.g., "q1.1", "q3.0")
  pattern_full <- "^q(\\d+)\\.(\\d)$"
  matches <- regmatches(Qsg, regexec(pattern_full, Qsg))[[1]]

  if (length(matches) >= 3) {
    idx <- as.integer(matches[2])
    action <- matches[3]
  } else {
    # Try short format: q<index> (e.g., "q1") -- default action = 1
    pattern_short <- "^q(\\d+)$"
    matches <- regmatches(Qsg, regexec(pattern_short, Qsg))[[1]]

    if (length(matches) >= 2) {
      idx <- as.integer(matches[2])
      action <- "1"
    } else {
      return(Qsg)
    }
  }

  if (idx < 1 || idx > length(confs_labels)) return(Qsg)

  base_label <- confs_labels[idx]

  if (action == "0") {
    paste0("!{", base_label, "}")
  } else {
    paste0("{", base_label, "}")
  }
}




#' Sort Subgroups by Focus (post-consistency)
#'
#' Sorts a data.table of subgroup results according to the specified focus.
#' For \code{"hrMaxSG"} and \code{"hrMinSG"}, candidates are first
#' partitioned into a \emph{candidate inclusion set} (the "band"), then
#' ranked within the band by sample size.  The \code{selection_rule}
#' argument controls how the band is defined:
#' \itemize{
#'   \item \code{"neighborhood"} -- 1D effect-size band (default; legacy
#'     behaviour): a candidate is in-band iff its effect is within
#'     \code{effect_neighborhood} of the maximum.
#'   \item \code{"pareto"} -- 2D Pareto-non-dominated set in
#'     (effect, N) space: a candidate is in-band iff no other candidate
#'     has both a larger effect and a larger N.
#'   \item \code{"both"} -- intersection of the two: a candidate is
#'     in-band iff it satisfies the neighborhood test \emph{and} is
#'     Pareto-non-dominated.
#' }
#' See Details.
#'
#' @param result_new A data.table of subgroup results with columns
#'   \code{Pcons}, \code{hr}, \code{N}, \code{K}.
#' @param sg_focus Character. Sorting focus.  One of \code{"hr"},
#'   \code{"maxSG"}, \code{"hrMaxSG"}, \code{"minSG"}, \code{"hrMinSG"}.
#' @param selection_rule Character. Rule defining the candidate
#'   inclusion set for \code{"hrMaxSG"} / \code{"hrMinSG"}.  One of
#'   \code{"neighborhood"} (default; current behaviour),
#'   \code{"pareto"}, or \code{"both"}.  Must be \code{"neighborhood"}
#'   for other \code{sg_focus} values.
#' @param effect_neighborhood Numeric in \code{[0, 1)}.  Relative
#'   tolerance defining the effect-size neighborhood for
#'   \code{"hrMaxSG"} and \code{"hrMinSG"}.  A candidate is in the
#'   neighborhood iff its (natural-scale) effect is at least
#'   \code{(1 - effect_neighborhood) * max(effect)}.  Default
#'   \code{0.10} (i.e., within 10\% of the strongest effect).  Used by
#'   \code{selection_rule} \code{"neighborhood"} and \code{"both"};
#'   ignored by \code{"pareto"} (and rejected if non-default with
#'   \code{"pareto"}).
#' @param effect_log_scale Logical.  If \code{TRUE}, the \code{hr}
#'   column stores log-scale values (log-OR, log-IRR) and is
#'   exponentiated before applying the inclusion test.  Default
#'   \code{FALSE} (Cox stores natural-scale HR).
#'
#' @details
#' \strong{Sort keys by \code{sg_focus}:}
#' \describe{
#'   \item{\code{"hr"}}{\code{(-Pcons, -hr, K)} -- prefer high consistency
#'     and large effect.}
#'   \item{\code{"maxSG"}}{\code{(-N, -Pcons, K)} -- prefer large
#'     subgroups with high consistency.}
#'   \item{\code{"minSG"}}{\code{(N, -Pcons, K)} -- prefer small
#'     subgroups with high consistency.}
#'   \item{\code{"hrMaxSG"}}{\code{(-in_band, -N, -Pcons, -hr, K)} -- among
#'     candidates in the inclusion band defined by \code{selection_rule},
#'     prefer the largest sample size.}
#'   \item{\code{"hrMinSG"}}{\code{(-in_band, N, -Pcons, -hr, K)} -- among
#'     candidates in the inclusion band defined by \code{selection_rule},
#'     prefer the smallest sample size.}
#' }
#'
#' \code{"hrMaxSG"} and \code{"hrMinSG"} use a \emph{lexicographic with
#' candidate band} rule: effect size is primary, but a bounded set of
#' candidates near (or trading off against) the maximum is accepted to
#' optimize sample size.
#'
#' \strong{Selection rules:}
#' \describe{
#'   \item{\code{"neighborhood"}}{1D \eqn{\varepsilon}-band on the effect
#'     scale: \code{in_band = (hr >= (1 - effect_neighborhood) * max(hr))}.
#'     Setting \code{effect_neighborhood = 0} reduces this to a strict
#'     max-effect filter.  This is the legacy behaviour.}
#'   \item{\code{"pareto"}}{2D Pareto-non-dominated set: a candidate is
#'     in-band iff no other candidate has both \code{hr_j >= hr_i} and
#'     \code{N_j >= N_i} with at least one strict inequality.  Reuses
#'     the same dominance computation as \code{compute_pareto_frontier()}.
#'     Removes the \code{effect_neighborhood} tuning parameter from the
#'     selection criterion.}
#'   \item{\code{"both"}}{Intersection: in-band iff in the
#'     \eqn{\varepsilon}-band \emph{and} Pareto-non-dominated.  Strictest
#'     of the three.}
#' }
#'
#' @return A sorted data.table.  The top row is the selected subgroup
#'   under \code{sg_focus}; remaining rows are diagnostic.
#'
#' @examples
#' \donttest{
#' library(data.table)
#' dt <- data.table(Pcons = c(0.92, 0.95, 0.88, 0.90),
#'                  hr    = c(2.5,  2.4,  2.3,  2.0),
#'                  N     = c(70,   100,  150,  200),
#'                  K     = c(1, 2, 1, 2))
#' # hrMaxSG with 10% neighborhood: HR 2.5, 2.4, 2.3 are in-band;
#' # among those, N = 150 is largest -> top row is N=150, hr=2.3.
#' sort_subgroups(dt, sg_focus = "hrMaxSG", effect_neighborhood = 0.10)
#'
#' # Same data with selection_rule = "pareto": the dominance set in
#' # (hr, N) is used instead of the 1D effect band.
#' sort_subgroups(dt, sg_focus = "hrMaxSG", selection_rule = "pareto")
#' }
#' @importFrom data.table setorder
#' @export
sort_subgroups <- function(result_new, sg_focus,
                           selection_rule = "neighborhood",
                           effect_neighborhood = 0.10,
                           effect_log_scale = FALSE) {

  .validate_selection_rule(selection_rule, sg_focus, effect_neighborhood)

  # Effect maximiser, no auxiliary condition.  Distinct from "hr", which sorts
  # (-Pcons, -hr, K) and lets consistency dominate.  K breaks exact ties toward
  # the simpler rule; it cannot override the effect ordering.
  if (sg_focus == "maxeff") {
    data.table::setorder(result_new, -hr, K)
    return(result_new)
  }
  if (sg_focus == "hr") {
    data.table::setorder(result_new, -Pcons, -hr, K)
    return(result_new)
  }
  if (sg_focus == "maxSG") {
    data.table::setorder(result_new, -N, -Pcons, K)
    return(result_new)
  }
  if (sg_focus == "minSG") {
    data.table::setorder(result_new, N, -Pcons, K)
    return(result_new)
  }

  if (sg_focus %in% c("hrMaxSG", "hrMinSG")) {

    hr_vec <- as.numeric(result_new$hr)
    if (isTRUE(effect_log_scale)) hr_vec <- exp(hr_vec)

    if (all(is.na(hr_vec))) {
      data.table::setorder(result_new, -Pcons, -hr, K)
      return(result_new)
    }

    N_vec     <- as.numeric(result_new$N)
    Pcons_vec <- as.numeric(result_new$Pcons)
    K_vec     <- as.numeric(result_new$K)

    in_band <- .compute_inclusion_band(
      hr_vec              = hr_vec,
      n_vec               = N_vec,
      selection_rule      = selection_rule,
      effect_neighborhood = effect_neighborhood
    )

    ord <- if (sg_focus == "hrMaxSG") {
      order(-in_band, -N_vec, -Pcons_vec, -hr_vec, K_vec)
    } else {
      order(-in_band,  N_vec, -Pcons_vec, -hr_vec, K_vec)
    }

    return(result_new[ord, ])
  }

  stop(sprintf("Unknown sg_focus value: %s", sg_focus))
}

#' Sort Subgroups by Focus (pre-consistency)
#'
#' Sorts a data.table of candidate subgroups before consistency
#' evaluation.  Mirrors \code{\link{sort_subgroups}} but operates on the
#' pre-consistency column convention (\code{HR}, \code{n}, \code{K})
#' and omits \code{Pcons} tiebreakers (consistency not yet available).
#'
#' @param result_new A data.table of candidate subgroups with columns
#'   \code{HR}, \code{n}, \code{K}.
#' @param sg_focus Character. Sorting focus.  One of \code{"hr"},
#'   \code{"maxSG"}, \code{"hrMaxSG"}, \code{"minSG"}, \code{"hrMinSG"}.
#' @param selection_rule Character. See \code{\link{sort_subgroups}}.
#'   Default \code{"neighborhood"}.
#' @param effect_neighborhood Numeric in \code{[0, 1)}.  See
#'   \code{\link{sort_subgroups}}.  Default \code{0.10}.
#' @param effect_log_scale Logical.  See \code{\link{sort_subgroups}}.
#'   Default \code{FALSE}.
#'
#' @return A sorted data.table.
#' @importFrom data.table setorder
#' @keywords internal
sort_subgroups_preview <- function(result_new, sg_focus,
                                   selection_rule = "neighborhood",
                                   effect_neighborhood = 0.10,
                                   effect_log_scale = FALSE) {

  .validate_selection_rule(selection_rule, sg_focus, effect_neighborhood)

  if (sg_focus == "maxeff") {
    data.table::setorder(result_new, -HR, K)
    return(result_new)
  }
  if (sg_focus == "hr") {
    data.table::setorder(result_new, -HR, K)
    return(result_new)
  }
  if (sg_focus == "maxSG") {
    data.table::setorder(result_new, -n, K)
    return(result_new)
  }
  if (sg_focus == "minSG") {
    data.table::setorder(result_new, n, K)
    return(result_new)
  }

  if (sg_focus %in% c("hrMaxSG", "hrMinSG")) {

    hr_vec <- as.numeric(result_new$HR)
    if (isTRUE(effect_log_scale)) hr_vec <- exp(hr_vec)

    if (all(is.na(hr_vec))) {
      data.table::setorder(result_new, -HR, K)
      return(result_new)
    }

    n_vec <- as.numeric(result_new$n)
    K_vec <- as.numeric(result_new$K)

    in_band <- .compute_inclusion_band(
      hr_vec              = hr_vec,
      n_vec               = n_vec,
      selection_rule      = selection_rule,
      effect_neighborhood = effect_neighborhood
    )

    # Frontier-preserving preview sort:
    # Include `-on_frontier` as a lex axis between `-in_band` and `-N`
    # so the full Pareto frontier (on (effect, N)) is guaranteed
    # inclusion in the top-`stop_Kgroups` candidates that go to
    # consistency.  Without this, restrictive rules like
    # `selection_rule = "both"` can crowd low-N frontier members out
    # of the top-K by filling slots with higher-N dominated
    # candidates (see issue documented in NEWS).
    # NB: this affects only candidate FILTERING for consistency
    # evaluation -- the post-consistency winner is still chosen by
    # sort_subgroups() (unchanged), so selection semantics are
    # preserved.
    on_frontier <- as.integer(!.pareto_dominated_xy(hr_vec, n_vec))

    ord <- if (sg_focus == "hrMaxSG") {
      order(-in_band, -on_frontier, -n_vec, -hr_vec, K_vec)
    } else {
      order(-in_band, -on_frontier,  n_vec, -hr_vec, K_vec)
    }

    return(result_new[ord, ])
  }

  stop(sprintf("Unknown sg_focus value: %s", sg_focus))
}

# Internal validator for effect_neighborhood -------------------------------
.validate_effect_neighborhood <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x >= 1) {
    stop("'effect_neighborhood' must be a single numeric in [0, 1).",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Internal validator for selection_rule ------------------------------------
# Shared by sort_subgroups() and sort_subgroups_preview().
# Enforces:
#   1. selection_rule must be one of c("neighborhood", "pareto", "both").
#   2. selection_rule != "neighborhood" only valid for sg_focus in
#      {"hrMaxSG", "hrMinSG"}.
#   3. selection_rule == "pareto" forbids non-default effect_neighborhood
#      (since the rule does not consult it).

.SELECTION_RULES <- c("neighborhood", "pareto", "both")
.EFFECT_NEIGHBORHOOD_DEFAULT <- 0.10

.validate_selection_rule <- function(selection_rule, sg_focus,
                                     effect_neighborhood) {
  if (!is.character(selection_rule) || length(selection_rule) != 1L ||
      !selection_rule %in% .SELECTION_RULES) {
    stop(sprintf(
      "'selection_rule' must be one of: %s.",
      paste(shQuote(.SELECTION_RULES), collapse = ", ")),
      call. = FALSE)
  }
  if (selection_rule != "neighborhood" &&
      !sg_focus %in% c("hrMaxSG", "hrMinSG")) {
    stop(sprintf(
      "'selection_rule = \"%s\"' is only meaningful for sg_focus in {\"hrMaxSG\", \"hrMinSG\"}; got sg_focus = \"%s\".",
      selection_rule, sg_focus),
      call. = FALSE)
  }
  if (selection_rule == "pareto" &&
      !isTRUE(all.equal(effect_neighborhood,
                        .EFFECT_NEIGHBORHOOD_DEFAULT))) {
    stop(sprintf(
      "'selection_rule = \"pareto\"' does not consult 'effect_neighborhood', but a non-default value (%g) was supplied.  Use 'selection_rule = \"both\"' if you want both criteria, or leave 'effect_neighborhood' at its default (%g).",
      effect_neighborhood, .EFFECT_NEIGHBORHOOD_DEFAULT),
      call. = FALSE)
  }
  invisible(TRUE)
}

# Shared inclusion-band computation ----------------------------------------
# Used by both sort_subgroups() (post-consistency: hr / N columns) and
# sort_subgroups_preview() (pre-consistency: HR / n columns).  Vector-
# based to avoid a column-name dependency.  Returns an integer 0/1
# indicator vector (1 = in candidate band) of length length(hr_vec).
# The validator is assumed to have already enforced selection_rule
# consistency.

.compute_inclusion_band <- function(hr_vec, n_vec, selection_rule,
                                    effect_neighborhood) {
  n <- length(hr_vec)
  in_nbhd <- if (selection_rule %in% c("neighborhood", "both")) {
    .validate_effect_neighborhood(effect_neighborhood)
    hr_max   <- max(hr_vec, na.rm = TRUE)
    hr_floor <- (1 - effect_neighborhood) * hr_max
    as.integer(!is.na(hr_vec) & hr_vec >= hr_floor)
  } else {
    rep(1L, n)
  }
  in_pareto <- if (selection_rule %in% c("pareto", "both")) {
    dominated <- .pareto_dominated_xy(hr_vec, n_vec)
    as.integer(!dominated)
  } else {
    rep(1L, n)
  }
  # AND-combination: candidate is in-band only if both criteria are
  # satisfied.  For pure rules one of the two vectors is all 1s, so the
  # product reduces to the active rule's indicator.
  in_nbhd * in_pareto
}

#' Compute Pareto Frontier on (Effect, N)
#'
#' Returns the subset of candidate subgroups that are not dominated on
#' the two-objective space (effect size, sample size), where both
#' objectives are maximized.  A candidate \eqn{i} is dominated iff
#' another candidate \eqn{j} has \eqn{hr_j \ge hr_i} and
#' \eqn{N_j \ge N_i}, with at least one inequality strict.
#'
#' @param result_dt A data.table with numeric columns \code{hr} and
#'   \code{N}.  Typically the post-consistency result table.
#' @param effect_log_scale Logical.  If \code{TRUE}, exponentiate
#'   \code{hr} before comparing (so the dominance test is on the
#'   natural scale for ratio measures).  Default \code{FALSE}.
#'
#' @return A data.table of non-dominated rows, sorted by \code{hr}
#'   descending.  Returns an empty 0-row data.table if input is empty.
#'
#' @details The frontier is intended as a \strong{post-hoc reporting}
#'   artifact, not a selection criterion.  The selected subgroup
#'   (under \code{sg_focus}) may or may not appear on the frontier --
#'   in particular, \code{"hrMinSG"} may select an N-dominated point
#'   by design (preferring small subgroups).
#'
#' @keywords internal
compute_pareto_frontier <- function(result_dt, effect_log_scale = FALSE) {
  if (!data.table::is.data.table(result_dt) || nrow(result_dt) == 0L) {
    return(if (data.table::is.data.table(result_dt)) result_dt[integer(0), ]
           else data.table::data.table())
  }

  is_dominated <- pareto_dominated_flags(result_dt,
                                         effect_log_scale = effect_log_scale)
  frontier <- result_dt[!is_dominated, ]
  if (nrow(frontier) > 1L) data.table::setorder(frontier, -hr)
  frontier
}

#' Pareto Dominance Flags
#'
#' Internal helper.  For each row of \code{result_dt}, returns
#' \code{TRUE} if the row is dominated in (\code{hr}, \code{N}) space
#' (both maximized) by another row, and \code{FALSE} otherwise.  Rows
#' with \code{NA} \code{hr} or \code{N} are flagged as dominated.
#'
#' Used by both \code{compute_pareto_frontier()} (which filters by
#' the negation of this vector) and the \code{selection_rule = "pareto"}
#' branches of \code{sort_subgroups()} / \code{sort_subgroups_preview()}
#' (which use it as an inclusion indicator).
#'
#' @param result_dt Data.table of candidate subgroups with columns
#'   \code{hr} and \code{N}.
#' @param effect_log_scale Logical.  If \code{TRUE}, \code{hr} is
#'   exponentiated before comparison.
#'
#' @return Logical vector of length \code{nrow(result_dt)}.
#' @keywords internal

pareto_dominated_flags <- function(result_dt, effect_log_scale = FALSE) {
  if (!data.table::is.data.table(result_dt) || nrow(result_dt) == 0L) {
    return(logical(0))
  }
  hr_vec <- as.numeric(result_dt$hr)
  if (isTRUE(effect_log_scale)) hr_vec <- exp(hr_vec)
  N_vec <- as.numeric(result_dt$N)
  .pareto_dominated_xy(hr_vec, N_vec)
}

# Vector-based core: O(n^2) dominance loop on two maximized vectors.
# Used by pareto_dominated_flags() (post-consistency, hr / N) and
# .compute_inclusion_band() (which both consistency stages call,
# passing the appropriate effect / size vectors regardless of the
# data.table column-name convention).

.pareto_dominated_xy <- function(x_vec, y_vec) {
  stopifnot(length(x_vec) == length(y_vec))
  n <- length(x_vec)
  is_dominated <- logical(n)
  if (n == 0L) return(is_dominated)
  for (i in seq_len(n)) {
    if (is.na(x_vec[i]) || is.na(y_vec[i])) {
      is_dominated[i] <- TRUE
      next
    }
    for (j in seq_len(n)) {
      if (i == j) next
      if (is.na(x_vec[j]) || is.na(y_vec[j])) next
      if (x_vec[j] >= x_vec[i] && y_vec[j] >= y_vec[i] &&
          (x_vec[j] >  x_vec[i] || y_vec[j] >  y_vec[i])) {
        is_dominated[i] <- TRUE
        break
      }
    }
  }
  is_dominated
}

#' Extract Subgroup Information
#'
#' Extracts subgroup definition and membership from results.
#'
#' @param df Data.frame. Original analysis data.
#' @param top_result Data.table row. Top subgroup result.
#' @param index.Z Matrix. Factor indicators for all subgroups.
#' @param names.Z Character vector. Factor column names.
#' @param confs_labels Character vector. Human-readable labels.
#'
#' @return List with sg.harm, sg.harm_label, df_flag, sg.harm.id.
#'
#' @keywords internal
extract_subgroup <- function(df, top_result, index.Z, names.Z, confs_labels) {
  m <- as.integer(top_result$m)
  indexm <- as.numeric(unlist(index.Z[m, ]))
  sg.harm <- names.Z[indexm == 1]
  sg.harm_label <- vapply(sg.harm, function(q) FS_labels(q, confs_labels), character(1))

  # Create membership flag
  df_temp <- data.table::as.data.table(df)
  for (varname in sg.harm) {
    df_temp <- df_temp[df_temp[[varname]] == 1, ]
  }

  sg.harm.id <- rep(0, nrow(df))
  if (nrow(df_temp) > 0 && "id" %in% names(df_temp)) {
    sg.harm.id[df$id %in% df_temp$id] <- 1
  }

  df_flag <- data.frame(
    id = df$id,
    treat.recommend = ifelse(sg.harm.id == 1, 0, 1)
  )

  list(
    sg.harm = sg.harm,
    sg.harm_label = sg.harm_label,
    df_flag = df_flag,
    sg.harm.id = sg.harm.id
  )
}


#' Plot Subgroup Survival Curves
#'
#' Plots weighted Kaplan-Meier survival curves for a specified subgroup and its complement using the \pkg{weightedsurv} package.
#'
#' @param df.sub A data frame containing data for the subgroup of interest.
#' @param df.subC A data frame containing data for the complement subgroup.
#' @param by.risk Numeric. The risk interval for plotting (passed to \code{weightedsurv::df_counting}).
#' @param confs_labels Named character vector. Covariate label mapping (not used directly in this function, but may be used for labeling).
#' @param this.1_label Character. Label for the subgroup being plotted.
#' @param top_result Data frame row. The top subgroup result row, expected to contain a \code{Pcons} column for consistency criteria.
#'
#' @importFrom weightedsurv df_counting plot_weighted_km
#' @keywords internal

plot_subgroup <- function(df.sub, df.subC, by.risk, confs_labels, this.1_label, top_result) {
  if (requireNamespace("weightedsurv", quietly = TRUE)) {
    tte.name <- "Y"
    event.name <- "Event"
    treat.name <- "Treat"
    con.lab <- "control"
    exp.lab <- "treat"
    dfcount <- weightedsurv::df_counting(df.sub, tte.name = tte.name, event.name = event.name, treat.name = treat.name, arms = c(exp.lab, con.lab), by.risk = by.risk)
    dfcountC <- weightedsurv::df_counting(df.subC, tte.name = tte.name, event.name = event.name, treat.name = treat.name, arms = c(exp.lab, con.lab), by.risk = by.risk)
    par(mfrow = c(1, 2))
    weightedsurv::plot_weighted_km(dfcount, conf.int = TRUE, show.logrank = TRUE, put.legend.lr = "topleft", ymax = 1.05, xmed.fraction = 0.65)
    weightedsurv::plot_weighted_km(dfcountC, conf.int = TRUE, show.logrank = TRUE, put.legend.lr = "topleft", ymax = 1.05, xmed.fraction = 0.65)
    cat("*** Subgroup found:", c(this.1_label), "\n")
    cat("% consistency criteria met=", c(top_result$Pcons), "\n")
  }  else {
    message("Package 'weightedsurv' not available: skipping weighted KM plots.")
  }
}


#' Output Subgroup Consistency Results
#'
#' Returns the top subgroup(s) and recommended treatment flags.
#'
#' @param df Data.frame. Original analysis data.
#' @param result_new Data.table. Sorted subgroup results.
#' @param sg_focus Character. Sorting focus criterion.
#' @param index.Z Matrix. Subgroup factor indicators.
#' @param names.Z Character vector. Factor column names.
#' @param details Logical. Print details.
#' @param plot.sg Logical. Plot subgroup curves.
#' @param by.risk Numeric. Risk interval for plotting.
#' @param confs_labels Character vector. Human-readable labels.
#' @param is_glm Logical. If \code{TRUE}, suppresses Kaplan-Meier survival
#'   plots (which are not meaningful for GLM outcomes).  Default \code{FALSE}.
#' @param selection_rule Character. Rule defining the candidate
#'   inclusion set for \code{"hrMaxSG"} / \code{"hrMinSG"}.  See
#'   \code{\link{sort_subgroups}}.  Default \code{"neighborhood"}.
#' @param effect_neighborhood Numeric in \code{[0, 1)}.  Relative
#'   tolerance for \code{"hrMaxSG"}/\code{"hrMinSG"} selection.  See
#'   \code{\link{sort_subgroups}}.  Default \code{0.10}.
#' @param effect_log_scale Logical.  If \code{TRUE}, the \code{hr}
#'   column stores log-scale values and is exponentiated for the
#'   neighborhood test and for Pareto-frontier dominance.  Default
#'   \code{FALSE}.
#'
#' @return List with elements:
#'   \describe{
#'     \item{result}{Sorted candidate table (top row = selected subgroup).}
#'     \item{pareto_frontier}{Data.table of non-dominated candidates on
#'       (effect, N), both maximized.  See
#'       \code{\link{compute_pareto_frontier}}.  May be \code{NULL} if
#'       computation failed.}
#'     \item{sg.harm}{Factor-level cut names defining the selected subgroup.}
#'     \item{sg.harm_label}{Human-readable subgroup labels.}
#'     \item{df_flag}{Per-subject treatment-recommendation flags.}
#'     \item{sg.harm.id}{Per-subject subgroup-membership indicator.}
#'   }
#'
#' @importFrom data.table copy
#' @examples
#' \dontrun{
#' # sg_consistency_out is called internally by forestsearch().
#' # See forestsearch() for the standard entry point.
#' fs <- forestsearch(gbsg,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
#' fs$grp.consistency$out_sg
#' }
#' @export
sg_consistency_out <- function(df, result_new, sg_focus, index.Z, names.Z,
                               details = FALSE, plot.sg = FALSE,
                               by.risk = 12, confs_labels,
                               is_glm = FALSE,
                               selection_rule = "neighborhood",
                               effect_neighborhood = 0.10,
                               effect_log_scale = FALSE) {

  result_new <- sort_subgroups(result_new, sg_focus,
                               selection_rule      = selection_rule,
                               effect_neighborhood = effect_neighborhood,
                               effect_log_scale    = effect_log_scale)
  top_result <- result_new[1, ]
  subgroup_info <- extract_subgroup(df, top_result, index.Z, names.Z, confs_labels)

  # === ADD THIS PLOTTING SECTION ===


  if (details && plot.sg && !is_glm) {
    sg.harm <- subgroup_info$sg.harm

    in_harm <- Reduce(`&`, lapply(sg.harm, function(v) df[[v]] == 1L))
    df.sub  <- df[in_harm, , drop = FALSE]
    df.subC <- df[!in_harm, , drop = FALSE]

    plot_subgroup(
      df.sub = df.sub,
      df.subC = df.subC,
      by.risk = by.risk,
      confs_labels = confs_labels,
      this.1_label = subgroup_info$sg.harm_label,
      top_result = top_result
    )
  }


  # === END PLOTTING SECTION ===

  result_out <- data.table::copy(result_new)

  # Pareto frontier on (effect, N), both maximized -- post-hoc diagnostic
  pareto_frontier <- tryCatch(
    compute_pareto_frontier(result_out, effect_log_scale = effect_log_scale),
    error = function(e) {
      warning("Failed to compute Pareto frontier: ", e$message,
              call. = FALSE)
      NULL
    }
  )

  list(
    result = result_out,
    pareto_frontier = pareto_frontier,
    sg.harm = subgroup_info$sg.harm,
    sg.harm_label = subgroup_info$sg.harm_label,
    df_flag = subgroup_info$df_flag,
    sg.harm.id = subgroup_info$sg.harm.id
  )
}

# =============================================================================
# DUPLICATE REMOVAL
# =============================================================================

#' Remove Near-Duplicate Subgroups
#'
#' Removes subgroups with nearly identical statistics (HR, n, E, etc.)
#' to reduce redundancy in candidate list.
#'
#' @param hr_subgroups Data.table of subgroup results.
#' @param tolerance Numeric. Tolerance for numeric comparison (default 0.001).
#' @param details Logical. Print details about removed duplicates.
#'
#' @return Data.table with near-duplicate rows removed.
#'
#' @importFrom data.table as.data.table
#' @keywords internal
remove_near_duplicate_subgroups <- function(hr_subgroups,
                                            tolerance = 0.001,
                                            details = FALSE) {

  df <- as.data.frame(hr_subgroups)

  # Columns to check: K, n, E, d1, m1, m0, HR, L(HR), U(HR)
  cols_to_check <- 2:min(10, ncol(df))

  df_rounded <- df
  for (i in cols_to_check) {
    if (is.numeric(df_rounded[, i])) {
      df_rounded[, i] <- round(df_rounded[, i] / tolerance) * tolerance
    }
  }

  key_cols <- df_rounded[, cols_to_check, drop = FALSE]
  dup_key <- apply(key_cols, 1, function(x) paste(x, collapse = "_"))

  keep_rows <- !duplicated(dup_key)

  n_removed <- sum(!keep_rows)

  if (details && n_removed > 0) {
    cat("Removed", n_removed, "near-duplicate subgroups\n")
  }

  return(data.table::as.data.table(hr_subgroups[keep_rows, ]))
}


#' Remove Redundant Subgroups
#'
#' Removes redundant subgroups by checking for exact ties in key columns.
#'
#' @param found.hrs Data.table of found subgroups.
#'
#' @return Data.table of non-redundant subgroups.
#'
#' @importFrom data.table data.table
#' @keywords internal
remove_redundant_subgroups <- function(found.hrs) {
  found.new <- found.hrs[order(found.hrs$HR, decreasing = TRUE), ]
  f1.hrs <- found.new[1, ]
  temp <- found.new[-c(1), ]
  temp2 <- as.matrix(found.new[, c("HR", "n", "E", "K", "L(HR)", "U(HR)")])
  id_keep <- which(round(diff(temp2, lag = 1), 6) != 0)
  fkeep.hrs <- temp[id_keep, ]
  na.omit(data.table::data.table(rbind(f1.hrs, fkeep.hrs)))
}


# =============================================================================
# FIXED-SAMPLE CONSISTENCY EVALUATION
# =============================================================================

#' Evaluate Single Subgroup for Consistency (Fixed-Sample)
#'
#' Evaluates a single subgroup for consistency across random splits using
#' a fixed number of splits.
#'
#' @param m Integer. Index of the subgroup to evaluate.
#' @param index.Z Data.table or matrix. Factor indicators for all subgroups.
#' @param names.Z Character vector. Names of factor columns.
#' @param df Data.frame. Original data with Y, Event, Treat, id columns.
#' @param found.hrs Data.table. Subgroup hazard ratio results.
#' @param n.splits Integer. Number of random splits for consistency evaluation.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param pconsistency.threshold Numeric. Minimum proportion of splits meeting
#'   consistency.
#' @param pconsistency.digits Integer. Rounding digits for consistency proportion.
#' @param maxk Integer. Maximum number of factors in a subgroup.
#' @param confs_labels Character vector. Labels for confounders.
#' @param details Logical. Print details during execution.
#' @param estimator_fn Closure or \code{NULL}. Effect-estimator closure from
#'   \code{\link{make_effect_estimator}} for GLM outcomes.  Default \code{NULL}
#'   (uses Cox model).
#' @param consistency_threshold Numeric or \code{NULL}. Threshold for GLM
#'   consistency evaluation.  Default \code{NULL} (uses \code{hr.consistency}).
#' @param adjust_covariates Character vector or \code{NULL}. Additional Cox
#'   model terms (e.g. \code{"strata(x1)"}) used when scoring survival
#'   subgroups.  Referenced columns must be present in \code{df}.  Ignored
#'   on the GLM path.  Default \code{NULL} (unadjusted).
#' @param consistency_method Character. \code{"split"} (default) or
#'   \code{"resample"}; see \code{\link{subgroup.consistency}}.  When
#'   \code{"resample"} yields an unsupported/non-convergent rate the evaluator
#'   falls back to literal splitting for that subgroup.
#' @param glm_resample_spec List or \code{NULL}. GLM resampling specification
#'   threaded from \code{\link{forestsearch}}; see
#'   \code{\link{subgroup.consistency}}.  \code{NULL} on the survival path.
#'
#' @return Named numeric vector with consistency results, or NULL if criteria
#'   not met.
#'
#' @examples
#' \dontrun{
#' # evaluate_subgroup_consistency() is called internally by forestsearch().
#' # Use forestsearch() as the standard entry point:
#' fs <- forestsearch(gbsg,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
#' }
#' Build the consistency evaluation closure in a clean environment
#'
#' The parallel consistency search applies a per-candidate closure via
#' `future_lapply`.  If that closure is defined inside `subgroup.consistency()`,
#' its environment is the whole enclosing frame, which can hold non-exportable
#' references (e.g. a screening forest's `externalptr`, or a large precomputed
#' object); `future` then fails to serialize it to the workers
#' ("non-exportable reference (externalptr) in ...future.FUN").  These factory
#' functions are defined at package scope, so the closure they return captures
#' only the small, exportable arguments passed here -- nothing from the caller's
#' frame.
#'
#' @return A function of one argument `m` (the candidate index).
#' @keywords internal
#' @noRd
.make_eval_subgroup_consistency <- function(index.Z, names.Z, df, found.hrs,
                                            n.splits, hr.consistency,
                                            pconsistency.threshold,
                                            pconsistency.digits, maxk, confs_labels,
                                            estimator_fn, consistency_threshold,
                                            adjust_covariates, consistency_method,
                                            glm_resample_spec) {
  function(m) {
    evaluate_subgroup_consistency(
      m = m, index.Z = index.Z, names.Z = names.Z, df = df,
      found.hrs = found.hrs, n.splits = n.splits,
      hr.consistency = hr.consistency,
      pconsistency.threshold = pconsistency.threshold,
      pconsistency.digits = pconsistency.digits, maxk = maxk,
      confs_labels = confs_labels, details = FALSE,
      estimator_fn = estimator_fn,
      consistency_threshold = consistency_threshold,
      adjust_covariates = adjust_covariates,
      consistency_method = consistency_method,
      glm_resample_spec = glm_resample_spec
    )
  }
}

#' Two-stage variant of the clean-environment consistency evaluation factory
#' @keywords internal
#' @noRd
.make_eval_consistency_twostage <- function(index.Z, names.Z, df, found.hrs,
                                            hr.consistency, pconsistency.threshold,
                                            pconsistency.digits, maxk, confs_labels,
                                            n.splits.screen, screen.threshold,
                                            n.splits.max, batch.size, conf.level,
                                            min.valid.screen, estimator_fn,
                                            consistency_threshold, adjust_covariates,
                                            consistency_method, glm_resample_spec) {
  function(m) {
    evaluate_consistency_twostage(
      m = m, index.Z = index.Z, names.Z = names.Z, df = df,
      found.hrs = found.hrs, hr.consistency = hr.consistency,
      pconsistency.threshold = pconsistency.threshold,
      pconsistency.digits = pconsistency.digits, maxk = maxk,
      confs_labels = confs_labels, details = FALSE,
      n.splits.screen = n.splits.screen, screen.threshold = screen.threshold,
      n.splits.max = n.splits.max, batch.size = batch.size,
      conf.level = conf.level, min.valid.screen = min.valid.screen,
      estimator_fn = estimator_fn,
      consistency_threshold = consistency_threshold,
      adjust_covariates = adjust_covariates,
      consistency_method = consistency_method,
      glm_resample_spec = glm_resample_spec
    )
  }
}

#' Evaluate Subgroup Consistency (Single-Stage)
#'
#' Evaluates a single candidate subgroup for treatment-effect consistency using
#' repeated random data splits (the standard, non-two-stage algorithm).  For
#' each split the subgroup effect is scored and compared against the harm
#' threshold; the consistency rate is the proportion of splits meeting it.
#' Called internally by \code{\link{forestsearch}} (via
#' \code{\link{subgroup.consistency}}); \code{\link{evaluate_consistency_twostage}}
#' is the screened, early-stopping counterpart used when
#' \code{use_twostage = TRUE}.
#'
#' @param m Integer. Index of subgroup to evaluate.
#' @param index.Z data.table or matrix. Factor indicators for all subgroups.
#' @param names.Z Character vector. Names of factor columns.
#' @param df data.frame. Original data with Y, Event, Treat, id columns.
#' @param found.hrs data.table. Subgroup hazard ratio results.
#' @param n.splits Integer. Number of random splits used to evaluate
#'   consistency.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param pconsistency.threshold Numeric. Final consistency threshold.
#' @param pconsistency.digits Integer. Rounding digits for output.
#' @param maxk Integer. Maximum number of factors in a subgroup.
#' @param confs_labels Character vector. Labels for confounders.
#' @param details Logical. Print progress details.
#' @param estimator_fn Closure or \code{NULL}. Effect-estimator closure from
#'   \code{\link{make_effect_estimator}} for GLM outcomes.  Default \code{NULL}
#'   (uses Cox model).
#' @param consistency_threshold Numeric or \code{NULL}. Threshold for GLM
#'   consistency evaluation.  Default \code{NULL} (uses \code{hr.consistency}).
#' @param adjust_covariates Character vector or \code{NULL}. Additional Cox
#'   model terms (e.g. \code{"strata(x1)"}) used when scoring survival
#'   subgroups.  Referenced columns must be present in \code{df}.  Ignored
#'   on the GLM path.  Default \code{NULL} (unadjusted).
#' @param consistency_method Character. \code{"split"} (default) or
#'   \code{"resample"}; see \code{\link{subgroup.consistency}}.  When
#'   \code{"resample"} engages it bypasses the split evaluation (the rate is
#'   computed from a single fit); an unsupported/non-convergent GLM rate falls
#'   through to the split evaluation.
#' @param glm_resample_spec List or \code{NULL}. GLM resampling specification
#'   threaded from \code{\link{forestsearch}}; see
#'   \code{\link{subgroup.consistency}}.  \code{NULL} on the survival path.
#' @param skip_consistency Logical. When \code{TRUE}, skip the consistency
#'   estimation entirely and return the effect-carry result row (effect from
#'   \code{found.hrs}, \code{Pcons = NA}) without forming the subgroup subset or
#'   refitting.  Used by the \code{sg_focus = "maxeff"} fast path, which needs
#'   effects for the argmax but consistency for the selected winner only.
#'   Default \code{FALSE}.
#'
#' @return Named numeric vector with consistency results, or \code{NULL} if the
#'   subgroup does not meet the criteria.
#'
#' @examples
#' \dontrun{
#' # evaluate_subgroup_consistency() is called internally by forestsearch().
#' # Use forestsearch() as the entry point:
#' fs <- forestsearch(gbsg,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
#' }
#' @importFrom data.table data.table
#' @importFrom survival coxph Surv
#' @export
evaluate_subgroup_consistency <- function(
    m,
    index.Z,
    names.Z,
    df,
    found.hrs,
    n.splits,
    hr.consistency,
    pconsistency.threshold,
    pconsistency.digits = 2,
    maxk,
    confs_labels,
    details = FALSE,
    estimator_fn = NULL,
    consistency_threshold = NULL,
    adjust_covariates = NULL,
    consistency_method = "split",
    glm_resample_spec = NULL,
    skip_consistency = FALSE
) {

  # -------------------------------------------------------------------------
  # SECTION 1: VALIDATE INPUT
  # -------------------------------------------------------------------------

  if (m < 1 || m > nrow(index.Z)) {
    warning("Invalid subgroup index m=", m, ". Skipping.")
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 2: EXTRACT SUBGROUP DEFINITION
  # -------------------------------------------------------------------------

  indexm <- as.numeric(unlist(index.Z[m, ]))

  if (length(indexm) != length(names.Z)) {
    warning("Subgroup ", m, ": index length mismatch. Skipping.")
    return(NULL)
  }

  this.m <- names.Z[indexm == 1]

  if (length(this.m) == 0) {
    warning("Subgroup ", m, ": no factors selected. Skipping.")
    return(NULL)
  }

  # Convert to labels
  this.m_label <- vapply(this.m, function(q) FS_labels(q, confs_labels), character(1))

  # -------------------------------------------------------------------------
  # SECTION 2b: SKIP CONSISTENCY (Pcons = NA)
  # -------------------------------------------------------------------------
  # Effect-carry-only path used by the sg_focus = "maxeff" fast path: the argmax
  # ranks on effect alone (found.hrs$HR, already computed by the search), so
  # consistency need not be evaluated per candidate.  Build the result row from
  # found.hrs with Pcons = NA and skip the subgroup subset + Cox refits below.
  # The winner is evaluated separately with skip_consistency = FALSE.
  if (isTRUE(skip_consistency)) {
    k <- length(this.m)
    Mnames <- paste("M", seq_len(maxk), sep = ".")
    mfound <- rep("", maxk)
    mfound[seq_len(k)] <- this.m_label
    resultk <- c(NA_real_, found.hrs$HR[m], found.hrs$n[m],
                 found.hrs$E[m], found.hrs$grp[m], m, k, mfound)
    names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)
    return(resultk)
  }

  # -------------------------------------------------------------------------
  # SECTION 3: IDENTIFY SUBGROUP MEMBERS
  # -------------------------------------------------------------------------

  df.x <- data.table::as.data.table(df)

  for (varname in this.m) {
    if (!varname %in% names(df.x)) {
      warning("Subgroup ", m, ": variable '", varname, "' not in df. Skipping.")
      return(NULL)
    }
    df.x <- df.x[df.x[[varname]] == 1, ]
  }

  N.x <- nrow(df.x)

  if (N.x < 10) {
    if (details) cat("Subgroup ", m, ": too few observations (", N.x, ")\n", sep = "")
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 4: INITIALIZE EFFECT ESTIMATOR
  # -------------------------------------------------------------------------

  if (!is.null(estimator_fn)) {
    # GLM path: closure is pre-built; no Cox init needed
    cox_init <- NULL
  } else {
    # Survival path: fit initial Cox model for warm-start
    cox_init <- tryCatch({
      fit0 <- survival::coxph(
        survival::Surv(Y, Event) ~ Treat,
        data = df.x,
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
      fit0$coefficients[1]
    }, error = function(e) 0)
  }

  # -------------------------------------------------------------------------
  # SECTION 5: RUN CONSISTENCY SPLITS
  # -------------------------------------------------------------------------

  if (identical(consistency_method, "resample") && is.null(estimator_fn)) {
    # Survival resampling approximation: one subgroup fit, no repeated splits.
    rr <- consistency_resample(
      df.x, hr.consistency = hr.consistency, method = "closed",
      adjust_covariates = adjust_covariates, cox_init = cox_init
    )
    p.consistency <- round(rr$rate_closed, pconsistency.digits)
    if (is.na(p.consistency)) {
      # Inestimable closed form (e.g. non-convergent fit on a candidate):
      # fall back to literal splitting rather than dropping the subgroup.
      # Mirrors the GLM path so Cox and GLM behave identically.
      p.consistency <- .consistency_via_splits(
        df.x, N.x, n.splits, hr.consistency, cox_init, estimator_fn,
        consistency_threshold, adjust_covariates, pconsistency.digits,
        m, details)
      if (is.na(p.consistency)) return(NULL)
    }
  } else if (identical(consistency_method, "resample") &&
             !is.null(estimator_fn) && !is.null(glm_resample_spec) &&
             .glm_resample_supported(glm_resample_spec)) {
    # GLM resampling approximation. comparison_threshold is already on the
    # comparison scale (log for ratio, identity for identity), matching what
    # the splitter compares against, so no further transform is applied.
    rr <- consistency_resample(
      df.x, method = "closed",
      outcome_type         = glm_resample_spec$outcome_type,
      effect_measure       = glm_resample_spec$effect_measure,
      comparison_threshold = glm_resample_spec$comparison_threshold,
      treat.name           = glm_resample_spec$treat.name,
      outcome.name         = glm_resample_spec$outcome.name,
      offset.name          = glm_resample_spec$offset.name,
      adjust_covariates    = glm_resample_spec$adjust_covariates,
      adverse_outcome      = glm_resample_spec$adverse_outcome
    )
    p.consistency <- round(rr$rate_closed, pconsistency.digits)
    if (is.na(p.consistency)) {
      # Unsupported configuration (e.g. non-convergent identity-link binomial):
      # fall back to literal splitting rather than dropping the subgroup.
      p.consistency <- .consistency_via_splits(
        df.x, N.x, n.splits, hr.consistency, cox_init, estimator_fn,
        consistency_threshold, adjust_covariates, pconsistency.digits,
        m, details)
      if (is.na(p.consistency)) return(NULL)
    }
  } else {
    p.consistency <- .consistency_via_splits(
      df.x, N.x, n.splits, hr.consistency, cox_init, estimator_fn,
      consistency_threshold, adjust_covariates, pconsistency.digits,
      m, details)
    if (is.na(p.consistency)) return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 7: CHECK CONSISTENCY THRESHOLD
  # -------------------------------------------------------------------------

  if (isTRUE(p.consistency < pconsistency.threshold)) {
    if (details) {
      cat("*** Not met: Subgroup, % Consistency =",
          c(this.m_label, p.consistency), "\n")
    }
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 8: FORMAT AND RETURN RESULT
  # -------------------------------------------------------------------------

  k <- length(this.m)
  covsm <- rep("M", maxk)
  mindex <- seq_len(maxk)
  Mnames <- paste(covsm, mindex, sep = ".")

  mfound <- matrix(rep("", maxk))
  mfound[seq_len(k)] <- this.m_label

  resultk <- c(
    p.consistency,
    found.hrs$HR[m],
    found.hrs$n[m],
    found.hrs$E[m],
    found.hrs$grp[m],
    m,
    k,
    mfound
  )

  names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)

  if (details) {
    cat("*** Met: Subgroup, % Consistency =",
        c(this.m_label, p.consistency), "\n")
  }

  return(resultk)
}


# =============================================================================
# TWO-STAGE CONSISTENCY EVALUATION
# =============================================================================

#' Evaluate Consistency (Two-Stage Algorithm)
#'
#' Evaluates a single subgroup for consistency using a two-stage approach:
#' Stage 1 screens with fewer splits, Stage 2 uses sequential batched
#' evaluation with early stopping for efficient evaluation.
#'
#' @param m Integer. Index of subgroup to evaluate.
#' @param index.Z data.table or matrix. Factor indicators for all subgroups.
#' @param names.Z Character vector. Names of factor columns.
#' @param df data.frame. Original data with Y, Event, Treat, id columns.
#' @param found.hrs data.table. Subgroup hazard ratio results.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param pconsistency.threshold Numeric. Final consistency threshold.
#' @param pconsistency.digits Integer. Rounding digits for output.
#' @param maxk Integer. Maximum number of factors in a subgroup.
#' @param confs_labels Character vector. Labels for confounders.
#' @param details Logical. Print progress details.
#' @param n.splits.screen Integer. Number of splits for Stage 1 (default 30).
#' @param screen.threshold Numeric. Screening threshold for Stage 1
#'   (default auto-calculated).
#' @param n.splits.max Integer. Maximum total splits (default 400).
#' @param batch.size Integer. Splits per batch in Stage 2 (default 20).
#' @param conf.level Numeric. Confidence level for early stopping (default 0.95).
#' @param min.valid.screen Integer. Minimum valid splits in Stage 1 (default 10).
#' @param estimator_fn Closure or \code{NULL}. Effect-estimator closure from
#'   \code{\link{make_effect_estimator}} for GLM outcomes.  Default \code{NULL}
#'   (uses Cox model).
#' @param consistency_threshold Numeric or \code{NULL}. Threshold for GLM
#'   consistency evaluation.  Default \code{NULL} (uses \code{hr.consistency}).
#' @param adjust_covariates Character vector or \code{NULL}. Additional Cox
#'   model terms (e.g. \code{"strata(x1)"}) used when scoring survival
#'   subgroups.  Referenced columns must be present in \code{df}.  Ignored
#'   on the GLM path.  Default \code{NULL} (unadjusted).
#' @param consistency_method Character. \code{"split"} (default) or
#'   \code{"resample"}; see \code{\link{subgroup.consistency}}.  When
#'   \code{"resample"} engages it bypasses the two-stage split screening (the
#'   rate is computed from a single fit); an unsupported/non-convergent GLM
#'   rate falls through to the Stage 1/2 splitting.
#' @param glm_resample_spec List or \code{NULL}. GLM resampling specification
#'   threaded from \code{\link{forestsearch}}; see
#'   \code{\link{subgroup.consistency}}.  \code{NULL} on the survival path.
#'
#' @return Named numeric vector with consistency results, or NULL if not met.
#'
#' @examples
#' \dontrun{
#' # evaluate_consistency_twostage() is called internally by forestsearch()
#' # when use_twostage = TRUE. Use forestsearch() as the entry point:
#' fs <- forestsearch(gbsg,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", treat.name = "hormon", event.name = "status",
#'   use_twostage = TRUE)
#' }
#' @importFrom data.table data.table as.data.table
#' @importFrom survival coxph Surv
#' @export
evaluate_consistency_twostage <- function(
    m,
    index.Z,
    names.Z,
    df,
    found.hrs,
    hr.consistency,
    pconsistency.threshold,
    pconsistency.digits = 2,
    maxk,
    confs_labels,
    details = FALSE,
    n.splits.screen = 30,
    screen.threshold = NULL,
    n.splits.max = 400,
    batch.size = 20,
    conf.level = 0.95,
    min.valid.screen = 10,
    estimator_fn = NULL,
    consistency_threshold = NULL,
    adjust_covariates = NULL,
    consistency_method = "split",
    glm_resample_spec = NULL
) {

  # ===========================================================================
  # EMBEDDED HELPER FUNCTIONS (for parallel execution compatibility)
  # ===========================================================================

  .wilson_ci <- function(x, n, conf.level = 0.95) {
    if (n == 0) {
      return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    z <- qnorm(1 - (1 - conf.level) / 2)
    p_hat <- x / n
    denom <- 1 + z^2 / n
    center <- (p_hat + z^2 / (2 * n)) / denom
    margin <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
    c(estimate = p_hat, lower = max(0, center - margin), upper = min(1, center + margin))
  }

  .early_stop_decision <- function(n_success, n_total, threshold, conf.level = 0.95, min_samples = 20) {
    if (n_total < min_samples) return("continue")
    ci <- .wilson_ci(n_success, n_total, conf.level)
    if (is.na(ci["lower"]) || is.na(ci["upper"])) return("continue")
    if (ci["lower"] >= threshold) return("pass")
    if (ci["upper"] < threshold) return("fail")
    return("continue")
  }

  .get_split_hr_fast <- function(df_split, cox_initial = NULL) {
    if (nrow(df_split) < 2 || sum(df_split$Event) < 2) return(NA_real_)
    adj <- .fs_adjust_terms(adjust_covariates)
    if (length(adj) > 0L) {
      fmla <- stats::reformulate(c("Treat", adj),
                                 response = "survival::Surv(Y, Event)")
    } else {
      fmla <- survival::Surv(Y, Event) ~ Treat
    }
    fit_args <- list(
      formula = fmla, data = df_split,
      robust = FALSE, model = FALSE, x = FALSE, y = FALSE
    )
    # `init = NULL` is an error in coxph; only pass the scalar warm-start
    # for the single-coefficient treatment-only model.
    if (length(adj) == 0L) fit_args$init <- cox_initial
    fit <- tryCatch(
      suppressWarnings(do.call(survival::coxph, fit_args)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    cf <- fit$coefficients
    beta <- if ("Treat" %in% names(cf)) cf[["Treat"]] else cf[1]
    return(exp(beta))
  }

  .run_single_consistency_split <- function(df.x, N.x, hr.cons, cox_init,
                                            est_fn = NULL, cons_thresh = NULL) {
    in.split1 <- tryCatch({
      sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
    }, error = function(e) NULL)
    if (is.null(in.split1)) return(NA_real_)

    df.x$insplit1 <- in.split1
    df.x.split1 <- df.x[df.x$insplit1 == TRUE, ]
    df.x.split2 <- df.x[df.x$insplit1 == FALSE, ]

    if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) return(NA_real_)

    # GLM path
    if (!is.null(est_fn) && !is.null(cons_thresh)) {
      for (half in list(df.x.split1, df.x.split2)) {
        treat_col <- intersect(c("Treat", "treat"), names(half))
        if (length(treat_col) > 0L) {
          tc <- treat_col[1L]
          if (sum(half[[tc]] == 1L) < 3L || sum(half[[tc]] == 0L) < 3L) {
            return(NA_real_)
          }
        }
      }
      res1 <- est_fn(df.x.split1)
      res2 <- est_fn(df.x.split2)
      if (is.na(res1$estimate) || is.na(res2$estimate)) return(NA_real_)
      return(as.numeric(
        res1$estimate >= cons_thresh && res2$estimate >= cons_thresh
      ))
    }

    # Survival path (existing)
    if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) return(NA_real_)

    hr.split1 <- .get_split_hr_fast(df_split = df.x.split1, cox_initial = cox_init)
    hr.split2 <- .get_split_hr_fast(df_split = df.x.split2, cox_initial = cox_init)

    if (!is.na(hr.split1) && !is.na(hr.split2)) {
      as.numeric(hr.split1 > hr.cons && hr.split2 > hr.cons)
    } else {
      NA_real_
    }
  }


  # ---------------------------------------------------------------------------
  # Parameter initialization
  # ---------------------------------------------------------------------------

  if (is.null(screen.threshold)) {
    se_estimate <- sqrt(pconsistency.threshold * (1 - pconsistency.threshold) / n.splits.screen)
    screen.threshold <- max(0.5, pconsistency.threshold - 2.5 * se_estimate)
  }

  # ---------------------------------------------------------------------------
  # Validate and extract subgroup
  # ---------------------------------------------------------------------------

  if (m < 1 || m > nrow(index.Z)) {
    warning("Invalid subgroup index m=", m, ". Skipping.")
    return(NULL)
  }

  indexm <- as.numeric(unlist(index.Z[m, ]))

  if (length(indexm) != length(names.Z)) {
    warning("Subgroup ", m, ": index length mismatch. Skipping.")
    return(NULL)
  }

  if (!all(indexm %in% c(0, 1, NA))) {
    warning("Subgroup ", m, ": invalid index values. Skipping.")
    return(NULL)
  }

  this.m <- names.Z[indexm == 1]

  if (length(this.m) == 0) {
    warning("Subgroup ", m, ": no factors selected. Skipping.")
    return(NULL)
  }

  this.m_label <- vapply(this.m, function(q) FS_labels(q, confs_labels), character(1))

  # ---------------------------------------------------------------------------
  # Identify subgroup members
  # ---------------------------------------------------------------------------

  df.x <- data.table::as.data.table(df)

  for (varname in this.m) {
    if (!varname %in% names(df.x)) {
      warning("Subgroup ", m, ": variable '", varname, "' not in df. Skipping.")
      return(NULL)
    }
    df.x <- df.x[df.x[[varname]] == 1, ]
  }

  N.x <- nrow(df.x)

  if (N.x < 10) {
    if (details) cat("Subgroup ", m, ": too few observations (", N.x, ")\n", sep = "")
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Initialize effect estimator
  # ---------------------------------------------------------------------------

  if (!is.null(estimator_fn)) {
    cox_init <- NULL
  } else {
    cox_init <- tryCatch({
      fit0 <- survival::coxph(
        survival::Surv(Y, Event) ~ Treat,
        data = df.x,
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
      fit0$coefficients[1]
    }, error = function(e) 0)
  }

  # ---------------------------------------------------------------------------
  # Resample short-circuit: the two-stage split machinery exists only to limit
  # the number of refits via early stopping; the resampling approximation
  # returns the rate from a single fit, so Stage 1/2 are bypassed entirely.
  # Survival path always; GLM path when a resample spec is supplied for a
  # supported measure. An unsupported/non-convergent GLM configuration falls
  # through to the Stage 1/2 splitting below.
  # ---------------------------------------------------------------------------
  cox_resample <- identical(consistency_method, "resample") && is.null(estimator_fn)
  glm_resample <- identical(consistency_method, "resample") &&
                  !is.null(estimator_fn) && !is.null(glm_resample_spec) &&
                  .glm_resample_supported(glm_resample_spec)

  if (cox_resample || glm_resample) {
    rr <- if (cox_resample) {
      consistency_resample(
        df.x, hr.consistency = hr.consistency, method = "closed",
        adjust_covariates = adjust_covariates, cox_init = cox_init)
    } else {
      consistency_resample(
        df.x, method = "closed",
        outcome_type         = glm_resample_spec$outcome_type,
        effect_measure       = glm_resample_spec$effect_measure,
        comparison_threshold = glm_resample_spec$comparison_threshold,
        treat.name           = glm_resample_spec$treat.name,
        outcome.name         = glm_resample_spec$outcome.name,
        offset.name          = glm_resample_spec$offset.name,
        adjust_covariates    = glm_resample_spec$adjust_covariates,
        adverse_outcome      = glm_resample_spec$adverse_outcome)
    }
    p.consistency <- round(rr$rate_closed, pconsistency.digits)

    # Inestimable closed-form rate (NA) for EITHER family -> fall through to the
    # Stage 1/2 splitting below rather than dropping the candidate.  Cox and GLM
    # behave identically here.
    resample_unsupported <- is.na(p.consistency)
    if (!resample_unsupported) {
      if (p.consistency < pconsistency.threshold) {
        if (details) {
          cat("Subgroup ", m, ": resample Pcons=", p.consistency,
              " (threshold ", pconsistency.threshold, ")\n", sep = "")
        }
        return(NULL)
      }
      k <- length(this.m)
      covsm <- rep("M", maxk)
      mindex <- seq_len(maxk)
      Mnames <- paste(covsm, mindex, sep = ".")
      mfound <- matrix(rep("", maxk))
      mfound[seq_len(k)] <- this.m_label
      resultk <- c(p.consistency, found.hrs$HR[m], found.hrs$n[m],
                   found.hrs$E[m], found.hrs$grp[m], m, k, mfound)
      names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)
      if (details) {
        cat("Subgroup ", m, ": PASSED via resample (Pcons=", p.consistency, ")\n", sep = "")
      }
      return(resultk)
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 1: Quick screening
  # ---------------------------------------------------------------------------

  stage1_flags <- numeric(n.splits.screen)
  for (i in seq_len(n.splits.screen)) {
    stage1_flags[i] <- .run_single_consistency_split(
      df.x, N.x, hr.consistency, cox_init,
      est_fn = estimator_fn, cons_thresh = consistency_threshold
    )
  }

  n_valid_stage1 <- sum(!is.na(stage1_flags))
  n_success_stage1 <- sum(stage1_flags == 1, na.rm = TRUE)

  if (n_valid_stage1 < min.valid.screen) {
    if (details) {
      cat("Subgroup ", m, ": insufficient valid Stage 1 splits (", n_valid_stage1, ")\n", sep = "")
    }
    return(NULL)
  }

  p_stage1 <- n_success_stage1 / n_valid_stage1

  if (p_stage1 < screen.threshold) {
    if (details) {
      cat("Subgroup ", m, ": SCREENED OUT (p=", round(p_stage1, 3), ")\n", sep = "")
    }
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Stage 2: Sequential evaluation with early stopping
  # ---------------------------------------------------------------------------

  if (details) {
    cat("Subgroup ", m, ": Stage 2 sequential\n", sep = "")
  }

  all_flags <- stage1_flags
  n_total_valid <- n_valid_stage1
  n_total_success <- n_success_stage1

  n_remaining <- n.splits.max - n.splits.screen
  n_batches <- ceiling(n_remaining / batch.size)

  final_decision <- "continue"

  for (batch_num in seq_len(n_batches)) {

    decision <- .early_stop_decision(
      n_success = n_total_success,
      n_total = n_total_valid,
      threshold = pconsistency.threshold,
      conf.level = conf.level,
      min_samples = max(20, min.valid.screen)
    )

    if (decision != "continue") {
      final_decision <- decision
      if (details) {
        cat("Subgroup ", m, ": EARLY ", toupper(decision),
            " at n=", n_total_valid, "\n", sep = "")
      }
      break
    }

    batch_flags <- numeric(batch.size)
    for (i in seq_len(batch.size)) {
      batch_flags[i] <- .run_single_consistency_split(
        df.x, N.x, hr.consistency, cox_init,
        est_fn = estimator_fn, cons_thresh = consistency_threshold
      )
    }

    n_batch_valid <- sum(!is.na(batch_flags))
    n_batch_success <- sum(batch_flags == 1, na.rm = TRUE)

    n_total_valid <- n_total_valid + n_batch_valid
    n_total_success <- n_total_success + n_batch_success
  }

  # ---------------------------------------------------------------------------
  # Final evaluation
  # ---------------------------------------------------------------------------

  if (final_decision == "continue") {
    final_decision <- .early_stop_decision(
      n_success = n_total_success,
      n_total = n_total_valid,
      threshold = pconsistency.threshold,
      conf.level = conf.level,
      min_samples = 20
    )
  }

  p.consistency <- tryCatch({
    round(n_total_success / n_total_valid, pconsistency.digits)
  }, error = function(e) {
    return(NA_real_)
  })

  if (is.na(p.consistency)) {
    return(NULL)
  }

  if (final_decision == "fail" || p.consistency < pconsistency.threshold) {
    if (details) {
      cat("Subgroup ", m, ": FAILED (Pcons=", p.consistency, ")\n", sep = "")
    }
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Format and return result
  # ---------------------------------------------------------------------------

  k <- length(this.m)
  covsm <- rep("M", maxk)
  mindex <- seq_len(maxk)
  Mnames <- paste(covsm, mindex, sep = ".")

  mfound <- matrix(rep("", maxk))
  mfound[seq_len(k)] <- this.m_label

  resultk <- c(
    p.consistency,
    found.hrs$HR[m],
    found.hrs$n[m],
    found.hrs$E[m],
    found.hrs$grp[m],
    m,
    k,
    mfound
  )

  names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)

  if (details) {
    cat("Subgroup ", m, ": PASSED (Pcons=", p.consistency,
        ", n_splits=", n_total_valid, ")\n", sep = "")
  }

  return(resultk)
}

#' Collapse identical-membership candidates for sg_focus = "maxeff"
#'
#' Two candidate cut-sets can define the SAME subject-level subgroup -- e.g.
#' \code{{er <= 0}} and \code{{er <= 0} & {er <= 10}} when \code{er <= 0}
#' implies \code{er <= 10}.  For the unconditional argmax these are exact
#' duplicates; evaluating both wastes compute and lets the \code{(-hr, K)} order
#' break a tie on cut count.  This keeps ONE representative per membership: the
#' fewest cuts (smallest \code{K}), ties broken by the lexicographically
#' smallest selected-column indices (a deterministic "first index").  No
#' distinct subgroup is dropped -- a collapsed encoding always shares its
#' membership with the retained representative.
#'
#' @param found.hrs data.table/data.frame of candidates; the \code{names.Z}
#'   columns are 0/1 selection indicators.
#' @param df Data frame carrying the \code{names.Z} cut columns (subject rows).
#' @param names.Z Character vector of cut-indicator column names.
#' @return \code{found.hrs} with duplicate-membership rows removed, original
#'   row order preserved.
#' @keywords internal
#' @noRd
.maxeff_membership_dedup <- function(found.hrs, df, names.Z) {
  fh <- as.data.frame(found.hrs)
  if (nrow(fh) <= 1L) return(found.hrs)

  Zmat <- as.matrix(df[, names.Z, drop = FALSE]);  storage.mode(Zmat) <- "double"
  Smat <- as.matrix(fh[, names.Z, drop = FALSE]);  storage.mode(Smat) <- "double"
  n_subj <- nrow(Zmat)

  keys <- character(nrow(fh))   # subject-membership key
  ktie <- character(nrow(fh))   # deterministic tiebreak: sorted selected indices
  Kvec <- integer(nrow(fh))
  for (r in seq_len(nrow(fh))) {
    sel <- which(Smat[r, ] == 1)
    Kvec[r] <- length(sel)
    memb <- if (length(sel) == 0L) rep(TRUE, n_subj)
            else rowSums(Zmat[, sel, drop = FALSE]) == length(sel)
    keys[r] <- paste0(which(memb), collapse = ",")
    ktie[r] <- paste(sprintf("%05d", sort(sel)), collapse = ".")
  }

  # Order by K asc, then first-index; keep the first row seen per membership.
  ord       <- order(Kvec, ktie)
  keep_rows <- ord[!duplicated(keys[ord])]
  found.hrs[sort(keep_rows), ]
}
