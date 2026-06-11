# =============================================================================
# Helper: Sync resolved values into args_call_all
# =============================================================================

#' Sync resolved parameter values into \code{args_call_all}
#'
#' For parameters whose values are resolved late in \code{forestsearch()}
#' (after the initial argument capture in Section 1B), this helper
#' copies the current value of each named formal from the calling
#' environment back into \code{args_call_all}.  Centralising this
#' pattern avoids the historical drift where scattered
#' \code{args_call_all$X <- X} mirror lines silently fell out of sync
#' with newly-added formals.
#'
#' @param args_call_all The captured argument list to update.
#' @param env Environment to read values from (typically
#'   \code{environment()} of the caller).
#' @param names Character vector of formal names to sync.
#' @return The updated \code{args_call_all} list.
#' @noRd
.sync_args_call_all <- function(args_call_all, env, names) {
  for (n in names) {
    args_call_all[[n]] <- get(n, envir = env)
  }
  args_call_all
}


# =============================================================================
# Default Worker Count for parallel_args
# =============================================================================

#' Compute the default worker count for \code{parallel_args}
#'
#' Returns \code{floor(0.80 * N)} where \code{N} is the number of
#' physical (non-hyperthreaded) cores, with a floor of 1 worker.
#' When the environment variable \code{_R_CHECK_LIMIT_CORES_} is
#' set (as during \code{R CMD check}), the result is capped at 2
#' to comply with CRAN policy.  When \code{detectCores(logical =
#' FALSE)} returns \code{NA} (some non-x86 platforms), falls back
#' to assuming 2 cores.
#'
#' This helper is used as the default value for
#' \code{parallel_args$workers} in \code{forestsearch()} and
#' related entry points; users can override at the call site.
#'
#' @return Positive integer worker count.
#' @noRd
.default_parallel_workers <- function() {
  n_phys <- parallel::detectCores(logical = FALSE)
  if (length(n_phys) != 1L || is.na(n_phys) || n_phys < 1L) n_phys <- 2L

  cran_cap <- isTRUE(as.logical(
    Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE")))
  if (cran_cap) {
    return(min(2L, as.integer(n_phys)))
  }

  max(1L, as.integer(floor(0.80 * n_phys)))
}


# =============================================================================
# Validation Helper: Outcome / Threshold Consistency
# =============================================================================

#' Validate outcome type, effect measure, and threshold consistency
#'
#' Checks for common misconfigurations that produce silent nonsense:
#' (1) binary data declared as survival, (2) ratio-scale thresholds
#' with identity-scale measures, (3) implausibly small ratio thresholds.
#'
#' @noRd
.validate_outcome_threshold_config <- function(
    df.analysis, outcome.name, event.name,
    outcome_type, effect_measure,
    hr.threshold, hr.consistency,
    user_set_threshold, user_set_consistency) {

  Y <- df.analysis[[outcome.name]]

  # --- Check 1: Data looks binary but outcome_type is survival ---
  y_unique <- sort(unique(Y[!is.na(Y)]))
  y_is_binary <- length(y_unique) == 2 && all(y_unique %in% c(0, 1))

  if (y_is_binary && outcome_type == "survival") {
    # Could be legitimate (0/1 event indicator with separate time variable)
    # but warn if there's no time-like variation
    if (!is.null(event.name) && event.name %in% names(df.analysis)) {
      E <- df.analysis[[event.name]]
      if (identical(sort(unique(Y)), sort(unique(E)))) {
        warning(
          sprintf(paste0(
            "outcome '%s' is binary (0/1) and identical to event ",
            "column '%s'.\n",
            "  If this is a binary endpoint (not survival), set ",
            "outcome_type = 'binary'.\n",
            "  If this IS survival, the outcome should be a ",
            "time-to-event variable, not the event indicator."),
            outcome.name, event.name),
          call. = FALSE)
      }
    } else if (is.null(event.name) || !event.name %in% names(df.analysis)) {
      warning(
        sprintf(paste0(
          "outcome '%s' is binary (0/1) but outcome_type = 'survival' ",
          "and no event column found.\n",
          "  Did you mean outcome_type = 'binary'?"),
          outcome.name),
        call. = FALSE)
    }
  }

  # --- Check 2: Data looks continuous but outcome_type is binary ---
  if (outcome_type == "binary" && length(y_unique) > 10) {
    warning(
      sprintf(paste0(
        "outcome '%s' has %d unique values but outcome_type = 'binary'.\n",
        "  Binary outcomes should have exactly 2 values (0/1).\n",
        "  Did you mean outcome_type = 'continuous' or 'survival'?"),
        outcome.name, length(y_unique)),
      call. = FALSE)
  }

  # --- Check 3: Suspiciously small thresholds for ratio-scale ---
  # If outcome_type is survival (ratio scale) and user passed very
  # small thresholds, they likely intended an identity-scale measure
  is_ratio <- outcome_type == "survival" ||
    (!is.null(effect_measure) && effect_measure %in% c("OR", "RR", "IRR"))

  if (is_ratio && user_set_threshold && hr.threshold > 0 &&
      hr.threshold < 0.5) {
    warning(
      sprintf(paste0(
        "effect.threshold = %.3f is very small for a ratio-scale ",
        "measure (%s).\n",
        "  HR/OR >= %.3f means virtually all subgroups pass screening.\n",
        "  If you intended an identity-scale threshold (e.g., risk ",
        "difference), set:\n",
        "    outcome_type = 'binary', effect_measure = 'RD', ",
        "effect.threshold = %.3f"),
        hr.threshold,
        if (outcome_type == "survival") "HR" else effect_measure,
        hr.threshold, hr.threshold),
      call. = FALSE)
  }

  if (is_ratio && user_set_consistency && hr.consistency > 0 &&
      hr.consistency < 0.5) {
    warning(
      sprintf(paste0(
        "consistency.threshold = %.3f is very small for a ratio-scale ",
        "measure (%s).\n",
        "  This means virtually all splits pass consistency.\n",
        "  If you intended an identity-scale threshold, set:\n",
        "    effect_measure = 'RD', consistency.threshold = %.3f"),
        hr.consistency,
        if (outcome_type == "survival") "HR" else effect_measure,
        hr.consistency),
      call. = FALSE)
  }

  # --- Check 4: Thresholds > 1 with identity-scale measure ---
  # RD/IRD are probability/rate differences bounded in [-1, 1].
  # MD (mean difference) has no natural bounds -- skip this check for MD.
  is_bounded_identity <- !is.null(effect_measure) &&
    effect_measure %in% c("RD", "IRD")

  if (is_bounded_identity && user_set_threshold && hr.threshold > 1.0) {
    warning(
      sprintf(paste0(
        "effect.threshold = %.2f but effect_measure = '%s' ",
        "(identity scale, range [-1, 1]).\n",
        "  A risk difference cannot exceed 1.0.\n",
        "  Did you mean effect_measure = 'OR' (ratio scale)?"),
        hr.threshold, effect_measure),
      call. = FALSE)
  }

  invisible(NULL)
}


#' @title ForestSearch: Exploratory Subgroup Identification
#'
#' @description
#' Implements advanced statistical methods for exploratory subgroup
#' identification in clinical trials.  Provides tools for identifying
#' patient subgroups with differential treatment effects using machine
#' learning approaches including Generalized Random Forests ('GRF'),
#' LASSO regularization, and exhaustive combinatorial search algorithms.
#' Supports survival endpoints (Cox proportional hazards), binary
#' outcomes (log odds ratio, log relative risk, risk difference),
#' continuous outcomes (mean difference), and count / rate outcomes
#' (log incidence rate ratio via Poisson, quasi-Poisson, or
#' negative-binomial GLMs with optional person-time offset).  Features
#' bootstrap bias correction using infinitesimal jackknife methods to
#' address selection bias in post-hoc analyses.  Designed for clinical
#' researchers conducting exploratory subgroup analyses in randomized
#' controlled trials, particularly for multi-regional clinical trials
#' ('MRCT') requiring regional consistency evaluation.  Methods are
#' described in Leon et al. (2024) \doi{10.1002/sim.10163}.
#'
#' @param df.analysis Data frame. Analysis dataset with required columns.
#' @param outcome.name Character. Name of time-to-event outcome variable. Default "tte".
#' @param event.name Character. Name of event indicator (1=event, 0=censored). Default "event".
#' @param treat.name Character. Name of treatment variable (1=treatment, 0=control). Default "treat".
#' @param id.name Character. Name of subject ID variable. Default "id".
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param flag_harm.name Character. Name of true harm flag for simulations (optional).
#' @param confounders.name Character vector. Names of candidate subgroup-defining variables.
#' @param parallel_args List. Parallel processing configuration:
#'   \describe{
#'     \item{plan}{Character. One of \code{"multisession"} (default;
#'       cross-platform, reuses workers across calls), \code{"multicore"}
#'       (Linux/macOS only, fastest startup but fails in RStudio interactive
#'       sessions), \code{"callr"} (slower per-call, fully isolated R
#'       processes), or \code{"sequential"} (no parallelism).}
#'     \item{workers}{Integer. Number of parallel workers.  Default:
#'       \code{floor(0.80 * N)} where \code{N} is the number of physical
#'       (non-hyperthreaded) cores, capped to 2 when
#'       \code{_R_CHECK_LIMIT_CORES_} is set (e.g., during \code{R CMD
#'       check}).  On a 4-core laptop the default is 3 workers; on a
#'       16-core workstation it is 12.}
#'     \item{show_message}{Logical. Show parallel setup messages.}
#'   }
#'   To disable parallelism entirely, pass
#'   \code{parallel_args = list(plan = "sequential")}.
#' @param df.predict Data frame. Prediction dataset (optional).
#' @param df.test Data frame. Test dataset (optional).
#' @param is.RCT Logical. Is this a randomized controlled trial? Default TRUE.
#' @param seedit Integer. Random seed. Default 8316951.
#' @param est.scale Character. Estimation scale ("hr" or "rmst"). Default "hr".
#' @param use_lasso Logical. Use LASSO for variable selection. Default TRUE.
#' @param use_grf Logical. Use GRF for variable importance. Default TRUE.
#' @param grf_res GRF results object (optional, for reuse).
#' @param grf_cuts List. Custom GRF cut points (optional).
#' @param use_dina Logical. Generate additional screening-stage candidate
#'   cuts from a DINA per-covariate Pareto frontier, fed into
#'   \code{\link{get_FSdata}} alongside the GRF cuts. Default FALSE
#'   (opt-in); existing behaviour is unchanged when FALSE.
#' @param dina_res A fitted DINA object (optional, for reuse). When NULL
#'   and \code{use_dina = TRUE}, a DINA model is fit at the screening
#'   stage.
#' @param dina_cuts Character vector of pre-supplied DINA cut expressions
#'   (optional). When supplied, the DINA fit / frontier step is skipped
#'   and these are used directly.
#' @param dina_args Named list of DINA tuning options, all optional.
#'   Recognised keys: \code{family} (defaults to the mapping from
#'   \code{outcome_type}: survival -> cox, binary -> binomial,
#'   continuous -> gaussian, count -> poisson), \code{seed} (defaults to
#'   \code{seedit}), and the fit-only keys \code{n_folds},
#'   \code{cens_type}, \code{cens_params} (passed to \code{\link{dina}}
#'   only when set); plus the frontier keys forwarded to
#'   \code{\link{dina_frontier}}: \code{scope} (default \code{"wide"}),
#'   \code{m_diff}, \code{n_min} (defaults to the \code{n.min} of this
#'   call), \code{direction}, \code{max_per_covariate}, \code{max_subgroups},
#'   and \code{digits}; plus the screening-behavior keys \code{selected_only}
#'   (default \code{TRUE}, matching GRF's \code{return_selected_cuts_only}),
#'   \code{max_depth} (\code{2L} default -- aligning DINA with \code{maxk = 2}
#'   and \code{grf_depth = 2} so the selected cut may be an AND-conjunction of
#'   two covariates; set \code{1L} for single-covariate cuts only), and
#'   \code{grid_probs}
#'   (the per-covariate quantile grid for depth-2 pair thresholds; default
#'   interior deciles).  \code{max_depth} / \code{grid_probs} are forwarded
#'   to \code{\link{dina_subgroup}} and apply to both the selected-cut
#'   screening and \code{subgroup_method = "dina"}.  (The selector statistic
#'   formerly accepted here as \code{select_statistic} is now the top-level
#'   \code{dina_select_statistic} argument; a value supplied inside
#'   \code{dina_args} is deprecated and overridden with a warning.)
#'   When \code{selected_only = TRUE} (the default), \code{use_dina} screening
#'   contributes the cut(s) chosen by \code{\link{dina_subgroup}} -- using
#'   this call's \code{sg_focus} / \code{selection_rule} /
#'   \code{effect_neighborhood} / \code{n.min} / \code{hr.threshold}, i.e. the
#'   same cut \code{subgroup_method = "dina"} would select. Set \code{FALSE} to
#'   contribute the full frontier candidate set instead. Unknown keys raise
#'   an error.
#' @param dina_select_statistic Character, one of \code{"dina"} (default) or
#'   \code{"effect"}, used only when \code{subgroup_method = "dina"}.  The DINA
#'   analogue of \code{grf_select_statistic}.  \code{"dina"} ranks DINA's
#'   qualifying family on its native subgroup-mean tau-hat (unchanged legacy
#'   behaviour).  \code{"effect"} re-ranks the same family on the inferential
#'   effect measure (Cox HR for survival; the resolved GLM effect otherwise),
#'   scored with the same per-candidate estimator the Tier-2 de-biased gate
#'   uses, so the realized selection is the exact event the gate de-biases (and
#'   the gate's native-family restriction follows the effect rather than
#'   tau-hat).
#' @param subgroup_method Character, one of \code{"consistency"} (default),
#'   \code{"dina"}, or \code{"grf"}. \code{"consistency"} is the standard
#'   ForestSearch pipeline (GRF/LASSO screening then the consistency search).
#'   \code{"dina"} delegates subgroup \emph{selection} to
#'   \code{\link{dina_subgroup}} and bypasses GRF, LASSO and the consistency
#'   search: a DINA model is fit, a single subgroup is selected using this
#'   call's \code{sg_focus} / \code{selection_rule} /
#'   \code{effect_neighborhood} / \code{n.min} / \code{hr.threshold}, and
#'   that subgroup is returned as \code{sg.harm}. \code{dina_args} supplies
#'   only the DINA \emph{fit} tuning in this mode (frontier keys are
#'   ignored). \code{"grf"} is the GRF analogue: it delegates selection to the
#'   GRF subgroup-identification routine
#'   (\code{\link{grf.subg.harm.survival}} for survival,
#'   \code{\link{grf.subg.harm.glm}} for GLM outcomes), again bypassing LASSO,
#'   DINA and the consistency search, and returns the policy-tree harm leaf as
#'   \code{sg.harm}. GRF-selection mode is configured from this call's existing
#'   GRF arguments (\code{frac.tau}, \code{dmin.grf}, \code{grf_depth},
#'   \code{is.RCT}, \code{seedit}, and the GLM extras) -- there is no separate
#'   argument list -- so it uses the same settings GRF screening would. For
#'   both \code{"dina"} and \code{"grf"}, the existing bootstrap de-biasing
#'   applies unchanged, re-fitting and re-selecting on each replicate. The
#'   selected subgroup may be a single covariate or, where the method's depth
#'   allows (DINA \code{max_depth = 2}; GRF \code{grf_depth >= 2}), an
#'   AND-conjunction.
#' @param max_n_confounders Integer. Maximum confounders to consider. Default 1000.
#' @param grf_depth Integer. GRF tree depth. Default 2.
#' @param grf_selection Character, one of \code{"tree"} (default) or
#'   \code{"frontier"}, used only when \code{subgroup_method = "grf"}.
#'   \code{"tree"} selects via the GRF policy tree. \code{"frontier"} is an
#'   \strong{experimental} alternative that selects from the Pareto frontier of
#'   the doubly-robust scores (see \code{\link{grf.subg.harm.survival}}); the
#'   frontier rule is taken from \code{sg_focus} -- all five values map through
#'   (\code{"eff"}, \code{"effMaxSG"}, \code{"effMinSG"}, \code{"maxSG"},
#'   \code{"minSG"}), aligned with \code{dina_subgroup()} -- and the relative
#'   band from \code{effect_neighborhood}. Provided for comparison; the tree is
#'   the recommended default.
#' @param grf_select_statistic Character, one of \code{"dr"} (default) or
#'   \code{"effect"}, used only when \code{subgroup_method = "grf"} and
#'   \code{grf_selection = "frontier"}.  \code{"dr"} selects the frontier winner
#'   on GRF's native doubly-robust score (legacy behaviour).  \code{"effect"}
#'   re-ranks the same DR-candidate frontier on the inferential effect measure
#'   (Cox HR for survival; the resolved GLM effect otherwise), scored with the
#'   same per-candidate estimator the Tier-2 de-biased gate uses, so the realized
#'   selection is the exact event the gate de-biases (and the gate's
#'   native-family restriction follows the effect rather than the DR score).
#'   Ignored for \code{grf_selection = "tree"} (the policy-tree leaf is the
#'   selection, with no enumerated family to rank).
#' @param dmin.grf Numeric. Minimum events for GRF. Default 0.0.
#' @param frac.tau Numeric in (0, 1]. Multiplier on the GRF time horizon
#'   passed to \code{grf::causal_survival_forest()}. The effective
#'   horizon is
#'   \code{frac.tau * min(max(Y[treated, event]), max(Y[control, event]))},
#'   i.e. \code{frac.tau} times the smaller of the two arm-specific
#'   maximum event times.  Values < 1 truncate the horizon inward of
#'   that maximum (the conservative choice when events are sparse near
#'   the tail).  Applies to both the GRF variable-importance
#'   pre-screen inside \code{forestsearch()} and the downstream call to
#'   \code{\link{grf.subg.harm.survival}}.  Used only when
#'   \code{outcome_type = "survival"}; ignored for GLM outcomes
#'   (binary, continuous, count).  Default 0.8.
#' @param return_selected_cuts_only Logical. If TRUE (default), GRF returns only cuts from the
#'   tree depth that identified the selected subgroup meeting `dmin.grf`. If FALSE
#'   returns all cuts from all fitted trees (depths 1 through `grf_depth`).
#'   See \code{\link{grf.subg.harm.survival}} for details.
#' @param conf_force Character vector. Variables to force include (optional).
#' @param defaultcut_names Character vector. Default cut variable names (optional).
#' @param cut_type Character. Cut type ("default" or "custom"). Default "default".
#' @param exclude_cuts Character vector. Variables to exclude from cutting (optional).
#' @param replace_med_grf Logical. Replace median with GRF cuts. Default FALSE.
#' @param cont.cutoff Integer. Cutoff for continuous vs categorical. Default 4.
#' @param conf.cont_medians Named numeric vector. Median values for continuous variables (optional).
#' @param conf.cont_medians_force Named numeric vector. Forced median values (optional).
#' @param conf.cont_jcuts Named integer list (each value >= 1). Per-variable
#'   J-quantile cut override for continuous covariates: for an entry
#'   \code{X = J}, the default cut set for X is replaced by \code{J}
#'   binary cut points at the (k/(J+1))-th quantiles of X
#'   (k = 1, ..., J), defining J+1 non-overlapping intervals.  Names
#'   must be in \code{confounders.name} and must not overlap with
#'   \code{defaultcut_names} or \code{conf.cont_medians_force}.  J-cut
#'   variables bypass LASSO filtering, matching \code{defaultcut_names}
#'   semantics.  Variables not listed retain default behaviour.
#'   Forwarded to \code{\link{get_FSdata}} via \code{filter_call_args()};
#'   see there for full semantics.
#' @param collapse_cuts Logical.  If TRUE, collapse near-redundant continuous
#'   candidate cuts before the search.  Continuous covariates can generate many
#'   thresholds that are practically redundant (e.g. \code{"age <= 35.7"} and
#'   \code{"age <= 35"}); when enabled, cuts on the same variable and operator
#'   within a per-variable standard-error band are merged to a single
#'   rounded-centroid threshold, subject to a membership safety check.
#'   Forwarded to \code{\link{get_FSdata}}; see \code{\link{collapse_redundant_cuts}}
#'   for the algorithm.  Default TRUE.  Pass \code{collapse_cuts = FALSE} to
#'   recover the un-coarsened candidate pool used by earlier package versions.
#' @param collapse_cuts_args List of overrides merged onto the defaults
#'   \code{list(c = 1.0, tol = 0.05, digits = 0L)}: band multiplier \code{c}
#'   (\code{band = c * sd(x)/sqrt(n)}), membership safety tolerance \code{tol}
#'   (fraction of n when < 1, absolute count when >= 1), and representative
#'   rounding \code{digits}.  Ignored when \code{collapse_cuts = FALSE}.
#' @param n.min Integer or \code{NULL}. Minimum subgroup size. Default 60.
#'   Supplying a value (or omitting the argument) uses that fixed floor, as in
#'   prior versions. Passing \code{n.min = NULL} opts into a sample-size-adaptive
#'   floor \code{max(60, ceiling(n.min.frac * N))}, where \code{N} is the
#'   analysis sample size -- at least 60, and at least \code{n.min.frac} of
#'   \code{N}. Growing the minimum subgroup size with \code{N} curbs false
#'   discovery from very small subgroups in large samples. The resolved integer
#'   is used everywhere downstream (screening, DINA, GRF, the consistency
#'   search, and the bootstrap).
#' @param n.min.frac Numeric in (0, 1). Fraction of the analysis sample size
#'   used for the adaptive \code{n.min} floor when \code{n.min = NULL}. Default
#'   0.10. Ignored when \code{n.min} is supplied.
#' @param effect.threshold Numeric or NULL. Screening threshold for candidate
#'   subgroups.  For ratio-scale measures (OR, RR, IRR, HR): on the ratio
#'   scale (e.g., 1.5 means OR >= 1.5).  For identity-scale measures
#'   (RD, IRD, MD): on the identity scale (e.g., 0.07 means RD >= 7 pct
#'   points).  If NULL (default), falls back to \code{hr.threshold}.
#' @param consistency.threshold Numeric or NULL. Threshold for split-sample
#'   consistency.  Each random 50/50 split must produce an estimate at or
#'   above this value.  Same scale conventions as \code{effect.threshold}.
#'   If NULL (default), falls back to \code{hr.consistency}.
#' @param hr.threshold Numeric. Legacy name for \code{effect.threshold}.
#'   Retained for backward compatibility.  Default 1.25.  When
#'   \code{effect.threshold} is provided, \code{hr.threshold} is ignored.
#' @param hr.consistency Numeric. Legacy name for \code{consistency.threshold}.
#'   Retained for backward compatibility.  Default 1.0.  When
#'   \code{consistency.threshold} is provided, \code{hr.consistency} is ignored.
#' @param sg_focus Character. Subgroup selection focus. One of:
#'   \describe{
#'     \item{\code{"hr"} (or \code{"eff"})}{Pick the candidate with the
#'       highest consistency (\code{Pcons}); ties broken by largest
#'       effect.}
#'     \item{\code{"maxSG"}}{Pick the candidate with the largest sample
#'       size; ties broken by consistency.}
#'     \item{\code{"minSG"}}{Pick the candidate with the smallest sample
#'       size; ties broken by consistency.}
#'     \item{\code{"hrMaxSG"} (or \code{"effMaxSG"})}{Among candidates
#'       with effect size within \code{effect_neighborhood} of the
#'       maximum, pick the one with the \emph{largest} sample size.}
#'     \item{\code{"hrMinSG"} (or \code{"effMinSG"})}{Among candidates
#'       with effect size within \code{effect_neighborhood} of the
#'       maximum, pick the one with the \emph{smallest} sample size.}
#'   }
#'   Default \code{"hr"}.  The \code{"eff*"} forms are aliases for the
#'   \code{"hr*"} forms and read more naturally in GLM contexts
#'   (continuous MD, binary OR/RR/RD, count IRR) where there is no
#'   hazard ratio.  Both vocabularies produce identical results; pick
#'   whichever fits the outcome type.  See \code{\link{sort_subgroups}}
#'   for full sort keys and tiebreakers.
#' @param selection_rule Character. Rule defining the candidate
#'   inclusion set for \code{"hrMaxSG"} / \code{"hrMinSG"}.  One of
#'   \code{"neighborhood"} (default; current behaviour),
#'   \code{"pareto"} (2D Pareto-non-dominated set in (effect, N)),
#'   or \code{"both"} (intersection of the two).  Forwarded to
#'   \code{\link{sort_subgroups}} via \code{filter_call_args()};
#'   see there for full semantics.  Must be \code{"neighborhood"} for
#'   other \code{sg_focus} values.
#' @param effect_neighborhood Numeric in \code{[0, 1)}.  Relative
#'   tolerance for \code{"hrMaxSG"} and \code{"hrMinSG"} selection.  A
#'   candidate is in the \emph{effect-size neighborhood} iff its
#'   (natural-scale) effect is at least
#'   \code{(1 - effect_neighborhood) * max(effect)}.  Default
#'   \code{0.10} (within 10\% of the strongest effect).  Setting
#'   \code{effect_neighborhood = 0} reduces \code{"hrMaxSG"}/
#'   \code{"hrMinSG"} to a strict max-effect filter (only candidates
#'   tied at the maximum qualify).  For GLM ratio measures (OR, IRR),
#'   the neighborhood test is applied on the natural scale.  Ignored
#'   when \code{sg_focus} is \code{"hr"}, \code{"maxSG"}, or
#'   \code{"minSG"}.
#' @param stop_threshold Numeric in \code{[0, 1]} or \code{NULL}.
#'   Early stopping threshold for consistency evaluation. When a candidate
#'   subgroup's consistency probability (\code{Pcons}) meets or exceeds this
#'   threshold, evaluation stops early -- remaining candidates are skipped.
#'   Set to \code{NULL} to disable early stopping and force evaluation of
#'   all candidates up to \code{max_subgroups_search}. Default \code{0.95}.
#'
#'   \strong{Note:} Values > 1.0 are not permitted. To disable early
#'   stopping, use \code{stop_threshold = NULL}, not a value above 1.
#'
#'   Automatically reset to \code{NULL} (with a warning) when
#'   \code{sg_focus} is \code{"hrMaxSG"} or \code{"hrMinSG"}, because
#'   neighborhood-based selection requires comparing effect sizes
#'   across all candidates to determine the in-neighborhood set.
#' @param fs.splits Integer. Number of splits for consistency evaluation (or maximum
#'   splits when \code{use_twostage = TRUE}). Default 1000.
#' @param m1.threshold Numeric. Maximum median survival threshold. Default Inf.
#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.
#' @param show_candidate_summary Logical. If \code{TRUE}, prints a
#'   post-consistency summary table of all passing candidates with
#'   Frontier/InBand/Selected flags.  See
#'   \code{\link{subgroup.consistency}} for details.  Default
#'   \code{FALSE}.
#' @param d0.min Integer. Minimum per-arm filter for candidate subgroups.
#'   For \code{"survival"} and \code{"binary"} outcomes, this is the minimum
#'   number of events in the control arm (for binary, events are Y = 1).
#'   Ignored for \code{"continuous"} outcomes (only \code{n.min} applies).
#'   Default 10.
#' @param d1.min Integer. Same as \code{d0.min} for the treatment arm.
#'   Default 10.
#' @param max.minutes Numeric. \strong{Currently inert; scheduled for
#'   deprecation in v0.3.0.} Previously intended as a wall-clock time
#'   budget for the combination search, this argument is no longer
#'   enforced in the parallelized search path and has no effect on
#'   behavior. Search scope is governed by \code{maxk}, candidate-factor
#'   screening (\code{use_grf}, \code{use_lasso}, \code{conf_force}),
#'   and the number of parallel workers. The default of 3 is retained
#'   only for signature compatibility.
#' @param minp Numeric. Minimum prevalence threshold. Default 0.025.
#' @param details Logical. Print progress details. Default FALSE.
#' @param quiet Logical. If TRUE, suppress the configuration summary
#'   message printed at startup.  Useful when \code{forestsearch()} is
#'   called repeatedly inside bootstrap or cross-validation loops.
#'   Default FALSE.
#' @param maxk Integer. Maximum number of factors per subgroup. Default 2.
#' @param by.risk Integer. Risk table interval. Default 12.
#' @param plot.sg Logical. Plot subgroup survival curves. Default FALSE.
#' @param plot.grf Logical. Plot GRF results. Default FALSE.
#' @param max_subgroups_search Integer. Maximum subgroups to evaluate.
#'   Default 30.  For `sg_focus` in `"effMaxSG"` / `"effMinSG"` (equivalently
#'   `"hrMaxSG"` / `"hrMinSG"`), the subgroup that is optimal under the
#'   criterion can sit further down the ranked candidate list, so the
#'   candidate pool should be broad: keep `max_subgroups_search >= 30` (a
#'   smaller value may miss the optimal subgroup).
#' @param vi.grf.min Numeric. Minimum GRF variable importance. Default -0.2.
#' @param use_twostage Logical. Use two-stage sequential consistency algorithm for
#'   improved performance. Default FALSE for backward compatibility. When TRUE,
#'   \code{fs.splits} becomes the maximum number of splits, and early stopping
#'   is enabled. See Details.
#' @param twostage_args List. Parameters for two-stage algorithm (only used when
#'   \code{use_twostage = TRUE}):
#'   \describe{
#'     \item{n.splits.screen}{Integer. Splits for Stage 1 screening. Default 30.}
#'     \item{screen.threshold}{Numeric. Consistency threshold for Stage 1. Default
#'       is automatically calculated to provide ~2.5 SE margin.}
#'     \item{batch.size}{Integer. Splits per batch in Stage 2. Default 20.}
#'     \item{conf.level}{Numeric. Confidence level for early stopping. Default 0.95.}
#'     \item{min.valid.screen}{Integer. Minimum valid Stage 1 splits. Default 10.}
#'   }
#' @param debias_gate Logical.  When \code{TRUE}, compute a fast
#'   multiplier-bootstrap de-biased estimate of the selected subgroup and flag
#'   whether it remains consistent with harm.  Supported on the GLM
#'   resample-consistency path (\code{consistency_method = "resample"} with a
#'   GLM effect) and on the survival/Cox path (\code{outcome_type = "survival"}),
#'   where it runs as a fast standalone approximation under the default
#'   split-consistency search.  This is a Tier-2 approximation of the full
#'   bootstrap bias-correction (\code{\link{forestsearch_bootstrap_dofuture}})
#'   and does not replace it.  Default \code{FALSE}.
#' @param debias_gate_args List of optional gate controls: \code{t_gate}
#'   (effect-scale gate; near-null default \code{1} for ratio measures, \code{0}
#'   for differences -- set near the null, not at the screen), \code{gate}
#'   (\code{"point"} or \code{"ci"}), \code{reselection} (bootstrap re-selection
#'   rule, default \code{"maxcons"}), \code{draws} (default \code{2000}),
#'   \code{multiplier} (default \code{"poisson"}), \code{include_complement}
#'   (default \code{TRUE}), \code{ci_method} (\code{"ij"} default -- de-biased CI
#'   from the infinitesimal-jackknife variance, the Tier-1 analogue -- or
#'   \code{"wald"} for the subgroup robust SE), \code{seed}, and
#'   \code{family_native_neighborhood}.  The last applies
#'   only to \code{subgroup_method = "dina"} and \code{"grf"}: a value in
#'   \code{[0, 1)} restricts the gate's re-selection family to candidates within
#'   that multiplicative band of the best \emph{native} statistic (DINA
#'   subgroup-mean tau-hat; GRF doubly-robust score), so the Cox/GLM-effect
#'   re-selection mirrors the method's own selection neighbourhood rather than
#'   ranging over the full candidate family -- curbing the over-correction that
#'   can occur when the native ranking diverges from the Cox/GLM effect.  It
#'   \strong{defaults to \code{effect_neighborhood}} (matching the band DINA/GRF
#'   select within, the recommended setting); pass a value \code{>= 1} to disable
#'   the restriction and de-bias against the full enumerated family.  Note this
#'   tunes only the candidate-competition part of the bias: the gate holds the
#'   forest fixed, so for DINA/GRF the full bootstrap remains the reference
#'   regardless of this value.  Ignored by \code{"consistency"}, whose
#'   \code{maxcons} re-selection already matches its search statistic.
#' @param consistency_method Character. \code{"split"} (default) runs the
#'   literal repeated 50/50 split-and-refit consistency calculation;
#'   \code{"resample"} uses the multiplier (influence-function / \code{dfbeta})
#'   approximation (\code{\link{consistency_resample}}), returning the rate from
#'   a single subgroup fit. The \code{"resample"} path applies to the Cox
#'   (survival) outcome and to GLM outcomes whose effect is a single model
#'   coefficient (OR, RR, RD, MD, IRR); configurations it cannot represent that
#'   way (IRD, propensity-adjusted effects, or a non-convergent fit on a given
#'   candidate) fall back to \code{"split"} automatically. When
#'   \code{use_twostage = TRUE}, \code{"resample"} bypasses the two-stage
#'   split screening entirely (the rate is computed directly).
#' @param outcome_type Character. One of \code{"survival"} (default),
#'   \code{"binary"}, \code{"continuous"}, or \code{"count"}.
#' @param effect_measure Character or \code{NULL}. Effect measure for GLM
#'   outcomes.  For binary: \code{"RD"}, \code{"OR"}, \code{"RR"},
#'   \code{"IRR"}, \code{"IRD"}.  For continuous: \code{"MD"}.
#'   For count: \code{"IRR"} (default), \code{"IRD"}.
#'   Default depends on \code{outcome_type}.
#' @param offset.name Character or \code{NULL}. Name of the follow-up time
#'   column for rate-based measures (IRR, IRD).
#' @param adverse_outcome Logical or \code{NULL}. If \code{TRUE}, higher
#'   outcome values indicate harm (default for binary and count).  If \code{FALSE},
#'   higher values indicate benefit (default for continuous).
#' @param overdispersion Character. Overdispersion correction for
#'   \code{outcome_type = "count"}.  One of \code{"none"} (default,
#'   standard Poisson), \code{"quasi"} (quasi-Poisson), or \code{"negbin"}
#'   (negative binomial via \code{MASS::glm.nb}).  Ignored for other
#'   outcome types.
#' @param grf_count_transform Character. Transformation applied to the
#'   count outcome before passing to \code{grf::causal_forest()} for
#'   factor screening.  \code{"log"} (default) applies
#'   \eqn{\log(Y + 0.5)}; \code{"identity"} passes raw counts.
#'   Ignored for non-count outcomes.
#' @param tune_grf Logical. If \code{TRUE}, enables cross-validated
#'   hyperparameter tuning for GRF's \code{causal_forest()} via
#'   \code{tune.parameters = "all"}.  If \code{FALSE} (default), all
#'   GRF parameters remain at their defaults, preserving existing
#'   behavior.  Applies to both internal GRF pre-screening and
#'   standalone GRF analysis when passed through
#'   \code{run_simulation_analysis()}.
#' @param adjust_covariates Character vector or \code{NULL}. Additional terms
#'   used to adjust the Cox model that scores survival subgroups during
#'   consistency evaluation.  Terms are appended verbatim to the model
#'   right-hand side after the treatment term, so survival modelling helpers
#'   may be used directly: \code{adjust_covariates = "strata(x1)"} fits a
#'   model stratified by \code{x1}, whereas \code{adjust_covariates = "x1"}
#'   includes \code{x1} as a linear covariate.  The referenced raw columns
#'   must be present in \code{df.analysis}; they are carried into the
#'   internal scoring frame automatically.  Applies to the survival (Cox)
#'   consistency path and to GLM outcomes (the effect-estimator closure
#'   adds the terms to the outcome model, e.g. logistic for binary/OR).
#'   The candidate-search ranking is adjusted on both paths.  Mutually
#'   exclusive with \code{ps_method != "none"}.  Default \code{NULL}
#'   (treatment-only, unadjusted -- the previous behaviour).
#' @param ps_method Character or \code{NULL}. Propensity score estimation
#'   method: \code{"grf"}, \code{"lasso"}, \code{"logistic"}, or
#'   \code{"none"}.  Default: \code{"none"} for RCT, \code{"grf"} for
#'   observational.
#' @param ps_adjust_method Character. PS adjustment method: \code{"none"}
#'   (default), \code{"iptw"} (stabilized inverse probability of treatment
#'   weights), or \code{"dr_gcomp"} (IPS covariate).
#' @param ps_hat Numeric vector or \code{NULL}. User-supplied propensity
#'   scores.  Must have length equal to \code{nrow(df.analysis)}.
#'
#' @return A list of class "forestsearch" containing:
#'   \describe{
#'     \item{grp.consistency}{Consistency evaluation results including:
#'       \itemize{
#'         \item out_sg: Selected subgroup based on sg_focus
#'         \item sg_focus: Focus criterion used
#'         \item df_flag: Treatment recommendations
#'         \item sg.harm: Subgroup definition labels (character vector
#'           of cut names) -- same as the top-level \code{sg.harm}
#'         \item sg.harm.id: per-subject 0/1 membership indicator
#'           (\emph{not} a character vector of cut expressions; see
#'           the \emph{Field naming collision} section below)
#'         \item algorithm: "twostage" or "fixed"
#'         \item n_candidates_evaluated: Number evaluated
#'         \item n_passed: Number passing threshold
#'       }}
#'     \item{find.grps}{Subgroup search results}
#'     \item{confounders.candidate}{Candidate confounders considered}
#'     \item{confounders.evaluated}{Confounders after variable selection}
#'     \item{df.est}{Analysis data with treatment recommendations}
#'     \item{df.predict}{Prediction data with recommendations (if provided)}
#'     \item{df.test}{Test data with recommendations (if provided)}
#'     \item{minutes_all}{Total computation time}
#'     \item{grf_res}{GRF results object}
#'     \item{sg_focus}{Subgroup focus criterion used}
#'     \item{sg.harm}{Selected subgroup definition -- character vector
#'       of factor-level cut names.  \strong{Recommended accessor} for
#'       displaying the identified subgroup (e.g.,
#'       \code{paste(obj$sg.harm, collapse = " & ")}).  \code{NULL} if
#'       no subgroup was identified.}
#'     \item{grf_cuts}{GRF cut points used}
#'     \item{prop_maxk}{Proportion of max combinations searched}
#'     \item{max_sg_est}{Maximum subgroup HR estimate}
#'     \item{grf_plot}{GRF plot object (if plot.grf = TRUE)}
#'     \item{args_call_all}{All arguments for reproducibility}
#'     \item{debias_gate}{Tier-2 de-biased gate result, or \code{NULL} when the
#'       gate did not run.  See \code{\link{fs_debias_gate}} for fields
#'       (\code{naive}, \code{debiased}, \code{selection_bias}, \code{gate},
#'       \code{harm_flag}, \code{timing_seconds}).  The de-biased CI is
#'       approximate (subgroup robust SE) and narrower than the Tier-1
#'       infinitesimal-jackknife CI.}
#'     \item{harm_flag_debiased}{Logical.  \code{TRUE} when the selected
#'       subgroup's de-biased estimate is still consistent with harm; \code{NA}
#'       when the gate did not run.}
#'   }
#'
#' @section Field naming collision with GRF results:
#' The top-level \code{sg.harm} on this object and the
#' \code{grp.consistency$sg.harm} nested field both hold a character
#' vector of cut expressions.  The nested \code{grp.consistency$sg.harm.id}
#' field holds a \emph{per-subject 0/1 membership indicator}, not a cut
#' vector -- unlike \code{sg.harm.id} on GRF result objects, which
#' \emph{does} hold cut expressions.  See
#' \code{\link{subgroup.consistency}} and \code{\link{grf.subg.harm.glm}}
#' for the full discussion of this naming collision.  Code that must
#' handle both object types should prefer the top-level \code{sg.harm}
#' accessor on forestsearch results.
#'
#' @details
#' \strong{Algorithm Overview:}
#' \enumerate{
#'   \item \strong{Variable Selection}: GRF identifies variables with treatment
#'     effect heterogeneity; LASSO selects most predictive
#'   \item \strong{Subgroup Discovery}: Exhaustive search over factor combinations
#'     up to \code{maxk}
#'   \item \strong{Consistency Validation}: Split-sample validation ensures
#'     reproducibility
#'   \item \strong{Selection}: Choose subgroup based on \code{sg_focus} criterion
#' }
#'
#' \strong{Two-Stage Consistency Algorithm:}
#' When \code{use_twostage = TRUE}, the consistency evaluation uses an optimized
#' algorithm that can provide 3-10x speedup:
#' \itemize{
#'   \item \strong{Stage 1}: Quick screening with \code{n.splits.screen} splits
#'     eliminates clearly non-viable candidates
#'   \item \strong{Stage 2}: Sequential batched evaluation with early stopping
#'     for candidates passing Stage 1
#' }
#'
#' The two-stage algorithm is recommended for:
#' \itemize{
#'   \item Exploratory analyses with many candidate subgroups
#'   \item Large \code{fs.splits} values (>200)
#'   \item Iterative model development
#' }
#'
#' For final regulatory submissions, \code{use_twostage = FALSE} may be preferred
#' for exact reproducibility.
#'
#' @examples
#' \dontrun{
#' # Example 1: Standard analysis (backward compatible)
#' result <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",
#'   hr.threshold = 1.25,
#'   pconsistency.threshold = 0.90,
#'   fs.splits = 400,
#'   details = TRUE
#' )
#'
#' # Example 2: Fast exploratory analysis with two-stage
#' result_fast <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "maxSG",
#'   hr.threshold = 1.15,
#'   pconsistency.threshold = 0.85,
#'   fs.splits = 500,
#'   use_twostage = TRUE,
#'   details = TRUE
#' )
#'
#' # Example 3: Two-stage with custom parameters
#' result_custom <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",
#'   hr.threshold = 1.3,
#'   pconsistency.threshold = 0.95,
#'   fs.splits = 600,
#'   use_twostage = TRUE,
#'   twostage_args = list(
#'     n.splits.screen = 50,
#'     batch.size = 25,
#'     conf.level = 0.99
#'   ),
#'   parallel_args = list(plan = "multisession", workers = 4),
#'   details = TRUE
#' )
#'
#' # Example 4: Binary outcome with OR thresholds
#' result_binary <- forestsearch(
#'   df.analysis = trial_data,
#'   outcome_type = "binary",
#'   effect_measure = "OR",
#'   effect.threshold = 1.5,
#'   consistency.threshold = 1.3,
#'   pconsistency.threshold = 0.90
#' )
#'
#' # Example 5: Binary outcome with RD thresholds
#' result_rd <- forestsearch(
#'   df.analysis = trial_data,
#'   outcome_type = "binary",
#'   effect_measure = "RD",
#'   effect.threshold = 0.07,       # 7 pct-point harm
#'   consistency.threshold = 0.03   # 3 pct-point per split
#' )
#' }
#'
#' @references
#' \itemize{
#'   \item FDA Guidance for Industry: Enrichment Strategies for Clinical Trials
#'   \item Athey & Imbens (2016). Recursive partitioning for heterogeneous
#'     causal effects. PNAS.
#'   \item Wager & Athey (2018). Estimation and inference of heterogeneous
#'     treatment effects using random forests. JASA.
#' }
#'
#' @seealso
#' \code{\link{subgroup.consistency}} for consistency evaluation details
#' \code{\link{forestsearch_bootstrap_dofuture}} for bootstrap inference
#' \code{\link{forestsearch_Kfold}} for cross-validation
#'
#' Package website: \url{https://larry-leon.github.io/forestsearch/}
#'
#' Source code: \url{https://github.com/larry-leon/forestsearch}
#'
#' @importFrom survival coxph Surv
#' @importFrom grf causal_survival_forest variable_importance
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table setorder
#' @importFrom stats quantile sd median
#' @importFrom stats complete.cases
#' @importFrom future.apply future_lapply
#' @importFrom randomForest randomForest
#' @importFrom weightedsurv df_counting
#' @export
forestsearch <- function(df.analysis,
                         outcome.name = "tte",
                         event.name = "event",
                         treat.name = "treat",
                         id.name = "id",
                         potentialOutcome.name = NULL,
                         flag_harm.name = NULL,
                         confounders.name = NULL,
                         parallel_args = list(plan = "multisession",
                                              workers = .default_parallel_workers(),
                                              show_message = TRUE),
                         df.predict = NULL,
                         df.test = NULL,
                         is.RCT = TRUE,
                         seedit = 8316951,
                         est.scale = "hr",
                         use_lasso = TRUE,
                         use_grf = TRUE,
                         grf_res = NULL,
                         grf_cuts = NULL,
                         use_dina = FALSE,
                         dina_res = NULL,
                         dina_cuts = NULL,
                         dina_args = list(),
                         dina_select_statistic = c("dina", "effect"),
                         subgroup_method = c("consistency", "dina", "grf"),
                         max_n_confounders = 1000,
                         grf_depth = 2,
                         grf_selection = c("tree", "frontier"),
                         grf_select_statistic = c("dr", "effect"),
                         dmin.grf = 0.0,
                         frac.tau = 0.8,
                         return_selected_cuts_only = TRUE,
                         conf_force = NULL,
                         defaultcut_names = NULL,
                         cut_type = "default",
                         exclude_cuts = NULL,
                         replace_med_grf = FALSE,
                         cont.cutoff = 4,
                         conf.cont_medians = NULL,
                         conf.cont_medians_force = NULL,
                         conf.cont_jcuts = NULL,
                         collapse_cuts = TRUE,
                         collapse_cuts_args = list(),
                         n.min = 60,
                         n.min.frac = 0.10,
                         effect.threshold = NULL,
                         consistency.threshold = NULL,
                         hr.threshold = 1.25,
                         hr.consistency = 1.0,
                         sg_focus = "hr",
                         selection_rule = "neighborhood",
                         effect_neighborhood = 0.10,
                         fs.splits = 1000,
                         m1.threshold = Inf,
                         pconsistency.threshold = 0.90,
                         stop_threshold = 0.95,
                         show_candidate_summary = FALSE,
                         d0.min = 10,
                         d1.min = 10,
                         max.minutes = 3,
                         minp = 0.025,
                         details = FALSE,
                         quiet = FALSE,
                         maxk = 2,
                         by.risk = 12,
                         plot.sg = FALSE,
                         plot.grf = FALSE,
                         max_subgroups_search = 30,
                         vi.grf.min = -0.2,
                         # NEW: Two-stage consistency parameters
                         use_twostage = TRUE,
                         twostage_args = list(),
                         consistency_method = c("split", "resample"),
                         # NEW: GLM outcome support
                         outcome_type = c("survival", "binary", "continuous",
                                          "count"),
                         effect_measure = NULL,
                         offset.name = NULL,
                         adverse_outcome = NULL,
                         overdispersion = c("none", "quasi", "negbin"),
                         grf_count_transform = c("log", "identity"),
                         tune_grf = FALSE,
                         # Covariate adjustment for subgroup evaluation (Cox)
                         adjust_covariates = NULL,
                         # Propensity score adjustment
                         ps_method = NULL,
                         ps_adjust_method = c("none", "iptw", "dr_gcomp"),
                         ps_hat = NULL,
                         # Tier-2 de-biased gate (optional; GLM resample path)
                         debias_gate = FALSE,
                         debias_gate_args = list()) {

  # ===========================================================================
  # SECTION 1A: RESOLVE FORMAL DEFAULTS BEFORE ARGUMENT CAPTURE
  # ===========================================================================
  # Resolve match.arg / default-NULL formals up front so the captured
  # args_call_all (Section 1B) holds the values that downstream do.call
  # sites need.  This eliminates the historical pattern of scattered
  # args_call_all$X <- X mirror lines for these parameters.

  outcome_type        <- match.arg(outcome_type)
  overdispersion      <- match.arg(overdispersion)
  grf_count_transform <- match.arg(grf_count_transform)
  subgroup_method     <- match.arg(subgroup_method)
  consistency_method  <- match.arg(consistency_method)

  # DINA's selector statistic is a first-class argument, mirroring
  # `grf_select_statistic`.  It is normalized here and injected as the canonical
  # value into `dina_args`, which `.resolve_dina_args()` reads, so every
  # downstream DINA site sees a single source of truth.  A value supplied the
  # old way (inside `dina_args$select_statistic`) is deprecated: it is honoured
  # for one cycle but overridden by the top-level argument, with a warning.
  dina_select_statistic <- match.arg(dina_select_statistic)
  if (is.null(dina_args)) dina_args <- list()
  if (!is.list(dina_args)) {
    stop("`dina_args` must be a list.", call. = FALSE)
  }
  if (!is.null(dina_args[["select_statistic"]]) &&
      !identical(dina_args[["select_statistic"]], dina_select_statistic)) {
    warning("`dina_args$select_statistic` is deprecated; use the top-level ",
            "`dina_select_statistic` argument (mirrors `grf_select_statistic`)",
            ".  The top-level value is used.", call. = FALSE)
  }
  dina_args[["select_statistic"]] <- dina_select_statistic

  if (outcome_type != "survival" && is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      binary     = "RD",
      continuous = "MD",
      count      = "IRR"
    )
  }

  if (is.null(adverse_outcome)) {
    adverse_outcome <- (outcome_type %in% c("binary", "count"))
  }

  # Threshold aliases: effect.threshold / consistency.threshold are the
  # preferred names; hr.threshold / hr.consistency are retained for
  # backward compatibility.  New names take precedence when both are
  # provided.
  user_set_threshold   <- !is.null(effect.threshold)      || !missing(hr.threshold)
  user_set_consistency <- !is.null(consistency.threshold) || !missing(hr.consistency)
  if (!is.null(effect.threshold))      hr.threshold   <- effect.threshold
  if (!is.null(consistency.threshold)) hr.consistency <- consistency.threshold

  # ===========================================================================
  # SECTION 1A2: RESOLVE ADAPTIVE n.min (opt-in via n.min = NULL)
  # ===========================================================================
  # Backward compatible: the default n.min = 60 is unchanged when supplied or
  # omitted.  Passing n.min = NULL opts into a sample-size-adaptive floor
  #   n.min <- max(60, ceiling(n.min.frac * N))
  # i.e. at least 60, and at least n.min.frac of the analysis sample size N.
  # This grows the minimum subgroup size with N to curb false discovery from
  # tiny subgroups in large samples.  Resolved here (before the args_call_all
  # capture below) so the integer value propagates to screening, DINA, GRF,
  # the consistency search, and the per-resample bootstrap, and is reported
  # in their n.min logging.  Guarded on a valid, non-empty df.analysis so the
  # zero-row early-exit path is unaffected.
  if (is.null(n.min)) {
    if (!is.data.frame(df.analysis) || nrow(df.analysis) == 0L) {
      n.min <- 60L   # degenerate input; zero-row guard below handles the exit
    } else {
      if (!is.numeric(n.min.frac) || length(n.min.frac) != 1L ||
          is.na(n.min.frac) || n.min.frac <= 0 || n.min.frac >= 1) {
        stop("n.min.frac must be a single number in (0, 1).", call. = FALSE)
      }
      N_analysis <- nrow(df.analysis)
      n.min <- max(60L, as.integer(ceiling(n.min.frac * N_analysis)))
      if (n.min >= N_analysis) {
        stop(sprintf(paste0("Adaptive n.min (%d) is >= the analysis sample ",
                            "size (%d); lower n.min.frac or set n.min ",
                            "explicitly."), n.min, N_analysis), call. = FALSE)
      }
      if (!quiet) {
        message(sprintf(
          "n.min = NULL -> resolved to %d  (max(60, ceiling(%.3g * %d))).",
          n.min, n.min.frac, N_analysis))
      }
    }
  }

  # ===========================================================================
  # SECTION 1B: CAPTURE ALL ARGUMENTS FOR REPRODUCIBILITY
  # ===========================================================================
  # Capture formals only (NOT as.list(environment())): Section 1A binds
  # non-formal locals (user_set_threshold, user_set_consistency) before
  # this point.  args_call_all flows back to bootstrap / CV consumers
  # that may reconstruct the call via do.call(forestsearch, args_call_all);
  # extra entries would error as 'unused argument'.

  args_names <- names(formals())
  args_call_all <- mget(args_names, envir = environment())

  # ===========================================================================
  # SECTION 1C: DEFENSIVE EARLY EXITS (shape guards)
  # ===========================================================================
  # Reject degenerate inputs before any heavy computation.  These cases
  # produce phantom output (or segfault inside compiled dependencies such
  # as glmnet's coxnet_exp on zero-row matrices) if allowed to proceed.
  #
  # Shape mirrors the get_FSdata-failure early return below (line ~1437):
  # list(sg.harm = NULL, args_call_all = ..., error_log = list(...)).
  # Downstream consumers (CV, bootstrap) inspect $sg.harm == NULL plus
  # $args_call_all to reconstruct calls, so the same shape works here.

  # Zero-row input.  glmnet::cv.glmnet -> coxnet_exp() segfaults on a
  # zero-row design matrix on Linux (the Fortran routine dereferences a
  # null pointer); macOS happens to mask this via different allocator
  # behaviour but the bug is present on both.  Catch it here.
  if (!is.data.frame(df.analysis) || nrow(df.analysis) == 0L) {
    if (!quiet) {
      warning("df.analysis has zero rows; returning early without ",
              "running the subgroup search.", call. = FALSE)
    }
    return(list(
      sg.harm       = NULL,
      args_call_all = args_call_all,
      error_log     = list(
        stage   = "input_validation",
        message = "df.analysis has zero rows; no subgroup search performed."
      )
    ))
  }

  # ===========================================================================
  # SECTION 2: VALIDATE INPUTS
  # ===========================================================================

  # Validate selection-related arguments up front, so parameter conflicts
  # (e.g., selection_rule = "pareto" with sg_focus = "hr", an unknown
  # sg_focus value, or effect_neighborhood out of range) are reported
  # immediately rather than after data preparation, factor selection, and
  # consistency evaluation have run.  The same validators are also called
  # inside sort_subgroups() / sort_subgroups_preview() / subgroup.consistency()
  # as defense-in-depth (those functions are exported and may be invoked
  # directly).

  # 1. sg_focus must be one of the known values.  Hoisted here so that
  #    .validate_selection_rule() below sees an already-validated sg_focus.
  #
  # First normalize the GLM-natural vocabulary (effMaxSG/effMinSG/eff) to
  # the canonical internal form (hrMaxSG/hrMinSG/hr).  The two
  # vocabularies are interchangeable at the user-facing API; internals
  # are keyed on the canonical names so all downstream code is unchanged.
  # Unrecognized values pass through and are reported by the whitelist
  # check below.
  sg_focus_user <- sg_focus                       # save for error message
  sg_focus      <- .normalize_sg_focus(sg_focus)
  valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG")
  valid_sg_focus_user <- c(valid_sg_focus,
                           "eff", "effMaxSG", "effMinSG")
  if (!is.character(sg_focus) || length(sg_focus) != 1L ||
      !sg_focus %in% valid_sg_focus) {
    stop(sprintf(
      "'sg_focus' must be one of: %s.  Got: %s.",
      paste(shQuote(valid_sg_focus_user), collapse = ", "),
      if (is.character(sg_focus_user) && length(sg_focus_user) == 1L)
        shQuote(sg_focus_user) else "<invalid>"),
      call. = FALSE)
  }
  # Sync the normalized sg_focus back into args_call_all (captured at
  # line 714 from the raw formals, BEFORE normalization).  Downstream
  # consumers -- plot_pareto_combined() in particular reads
  # fs$args_call_all$sg_focus to decide whether to draw the
  # effect-neighborhood band -- expect the canonical "hr*" form here.
  # Without this line, passing sg_focus = "effMaxSG" would store
  # "effMaxSG" in args_call_all and downstream feature checks against
  # c("hrMaxSG", "hrMinSG") would silently fall through to "no band".
  args_call_all$sg_focus <- sg_focus

  # 2. selection_rule / sg_focus / effect_neighborhood compatibility.
  .validate_selection_rule(selection_rule, sg_focus, effect_neighborhood)

  # 3. effect_neighborhood range -- only meaningful for hrMaxSG / hrMinSG
  #    (matches the existing conditional pattern in sort_subgroups()).
  if (sg_focus %in% c("hrMaxSG", "hrMinSG")) {
    .validate_effect_neighborhood(effect_neighborhood)
  }

  # Validate parallel arguments
  if (length(parallel_args) > 0) {
    allowed_plans <- c("multisession", "multicore", "callr", "sequential")
    plan_type <- parallel_args$plan
    n_workers <- parallel_args$workers
    max_cores <- parallel::detectCores()

    if (is.null(plan_type)) {
      stop("parallel_args$plan must be specified.")
    }

    if (!plan_type %in% allowed_plans) {
      stop("plan_type must be one of: ", paste(allowed_plans, collapse = ", "))
    }

    if (is.null(n_workers) || !is.numeric(n_workers) || n_workers < 1) {
      n_workers <- 1
    } else {
      n_workers <- min(n_workers, max_cores)
    }

    # Hoisted from subgroup.consistency() so an invalid batch_size errors
    # before any data preparation runs.  NULL means "use default downstream";
    # the original validator (still in place as defense-in-depth) is invoked
    # only when the value is actually consumed.
    if (!is.null(parallel_args$batch_size) &&
        (!is.numeric(parallel_args$batch_size) ||
         length(parallel_args$batch_size) != 1L ||
         is.na(parallel_args$batch_size) ||
         parallel_args$batch_size < 1 ||
         parallel_args$batch_size != as.integer(parallel_args$batch_size))) {
      stop("parallel_args$batch_size must be a positive integer.",
           call. = FALSE)
    }
  }

  # Early stopping is incompatible with compound sg_focus criteria that
  # require comparing HR *and* size across all candidates.
  if (!is.null(stop_threshold) && sg_focus %in% c("hrMaxSG", "hrMinSG")) {

    # Detect whether user explicitly set stop_threshold (vs. inheriting default)
    user_explicit <- !missing(stop_threshold)

    if (user_explicit) {
      warning(
        "stop_threshold = ", stop_threshold,
        " reset to NULL for sg_focus = '", sg_focus, "'.\n",
        "Neighborhood-based selection requires evaluating all candidates ",
        "to determine the effect-size neighborhood.\n",
        "To suppress this warning, pass stop_threshold = NULL explicitly.",
        call. = FALSE
      )
    }

    stop_threshold <- NULL
    args_call_all <- .sync_args_call_all(args_call_all, environment(),
                                         "stop_threshold")
  }

  # Validate two-stage parameters
 if (!is.logical(use_twostage) || length(use_twostage) != 1) {
    stop("'use_twostage' must be a single logical value (TRUE or FALSE)")
  }

  if (use_twostage && length(twostage_args) > 0) {
    valid_ts_args <- c("n.splits.screen", "screen.threshold", "batch.size",
                       "conf.level", "min.valid.screen")
    invalid_args <- setdiff(names(twostage_args), valid_ts_args)
    if (length(invalid_args) > 0) {
      warning("Unknown twostage_args parameters ignored: ",
              paste(invalid_args, collapse = ", "))
    }
  }

  # ===========================================================================
  # SECTION 2B: GLM OUTCOME TYPE SETUP
  # ===========================================================================

  outcome_type <- match.arg(outcome_type)
  overdispersion <- match.arg(overdispersion)
  grf_count_transform <- match.arg(grf_count_transform)

  # ===========================================================================
  # SECTION 2B-i: VALIDATE OUTCOME / THRESHOLD CONSISTENCY
  # ===========================================================================
  # Catch common misconfigurations before they produce silent nonsense.

  .validate_outcome_threshold_config(
    df.analysis  = df.analysis,
    outcome.name = outcome.name,
    event.name   = event.name,
    outcome_type = outcome_type,
    effect_measure = if (exists("effect_measure")) effect_measure else NULL,
    hr.threshold = hr.threshold,
    hr.consistency = hr.consistency,
    user_set_threshold = user_set_threshold,
    user_set_consistency = user_set_consistency
  )

  # Resolve adverse_outcome default:
  #   Binary outcomes: TRUE (event = adverse, the clinical trial convention)
  #   Count outcomes:  TRUE (higher count = more events = worse)
  #   Continuous outcomes: FALSE (user must set TRUE if higher Y = worse)
  #   Survival: not used (causal_survival_forest handles sign correctly)
  if (is.null(adverse_outcome)) {
    adverse_outcome <- (outcome_type %in% c("binary", "count"))
  }

  # Build effect estimator closure for non-survival outcomes
  estimator_fn <- NULL
  consistency_threshold <- NULL
  effect_threshold <- NULL


  if (outcome_type != "survival") {

    # Resolve default effect measure
    if (is.null(effect_measure)) {
      effect_measure <- switch(outcome_type,
        binary     = "RD",
        continuous = "MD",
        count      = "IRR"
      )
    }

    # Validate offset requirement for rate-based measures
    if (effect_measure %in% c("IRR", "IRD") && is.null(offset.name)) {
      stop(
        "effect_measure = '", effect_measure,
        "' requires `offset.name` (follow-up time column).",
        call. = FALSE
      )
    }

    # Build the estimator closure
    estimator_fn <- make_effect_estimator(
      outcome_type      = outcome_type,
      treat.name        = treat.name,
      outcome.name      = outcome.name,
      event.name        = if (outcome_type == "survival") event.name else NULL,
      offset.name       = offset.name,
      effect_measure    = effect_measure,
      adverse_outcome   = adverse_outcome,
      adjust_covariates = adjust_covariates
    )

    # Resolve screening and consistency thresholds from effect_measure.
    #
    # The parameters hr.threshold and hr.consistency are named for the
    # survival (hazard ratio) case.  For GLM outcomes they serve as
    # generic threshold parameters:
    #   - Ratio-scale measures (OR, RR, IRR): thresholds are on the
    #     ratio scale (e.g., 1.5 means OR >= 1.5).  Converted to log
    #     internally.
    #   - Identity-scale measures (RD, IRD, MD): thresholds are on the
    #     identity scale (e.g., 0.07 means RD >= 7 pct points).
    #     Values > 1 are nonsensical and trigger a warning.

    effect_threshold <- hr.threshold
    consistency_threshold <- hr.consistency

    if (effect_measure %in% c("RD", "IRD")) {

      # Remap survival defaults silently for RD/IRD
      if (!user_set_threshold && hr.threshold == 1.25) {
        effect_threshold <- switch(effect_measure, RD = 0.05, IRD = 0.01)
      }
      if (!user_set_consistency && hr.consistency == 1.0) {
        consistency_threshold <- 0.0
      }

      # Detect ratio-scale values passed to identity-scale measure
      # NOTE: RD/IRD are probability/rate differences bounded in [-1, 1].
      # MD (mean difference) has no natural bounds and is NOT checked here.
      if (user_set_threshold && effect_threshold > 1.0) {
        warning(
          sprintf(
            paste0(
              "effect.threshold = %.2f appears to be on a ratio scale, ",
              "but effect_measure = '%s' uses the identity scale ",
              "(range [-1, 1]).\n",
              "  Did you mean effect.threshold = %.2f on the %s scale?\n",
              "  Remapping to default: %.2f.\n",
              "  To use a custom %s threshold, pass a value in [-1, 1] ",
              "(e.g., effect.threshold = 0.07 for a 7 pct-point threshold)."),
            effect_threshold, effect_measure,
            effect_threshold, "OR", 0.05, effect_measure),
          call. = FALSE)
        effect_threshold <- switch(effect_measure,
          RD = 0.05, IRD = 0.01)
      }

      if (consistency_threshold >= 1.0 && !user_set_consistency) {
        # Default 1.0 maps to 0.0 ("any harm direction")
        consistency_threshold <- 0.0
      } else if (consistency_threshold > 1.0) {
        warning(
          sprintf(
            paste0(
              "consistency.threshold = %.2f appears to be on a ratio scale, ",
              "but effect_measure = '%s' uses the identity scale.\n",
              "  Remapping to default: 0.0.\n",
              "  To use a custom %s consistency threshold, pass a value ",
              "in [-1, 1]."),
            consistency_threshold, effect_measure, effect_measure),
          call. = FALSE)
        consistency_threshold <- 0.0
      }

    } else if (effect_measure == "MD") {

      # MD (mean difference) is on the natural outcome scale with no
      # fixed bounds.  Thresholds of 30, 50, etc. are legitimate for
      # outcomes like CD4 counts.  Only remap survival defaults if
      # user didn't explicitly set thresholds.
      if (!user_set_threshold && hr.threshold == 1.25) {
        # Survival default 1.25 is meaningless for MD; use 0.0
        effect_threshold <- 0.0
      }
      if (!user_set_consistency && hr.consistency == 1.0) {
        consistency_threshold <- 0.0
      }

    } else {
      # Log-scale measures (OR, RR, IRR): convert to log scale.
      # Ratio measures are positive by definition; non-positive thresholds
      # would either fail log() (yielding NaN/-Inf via log(<=0)) or, under
      # the prior > 0 guard, pass through unchanged and be silently
      # interpreted downstream as already-log-scale values.  Either path
      # is a silent wrong-answer, so reject explicitly.
      if (effect_threshold <= 0) {
        stop(
          "effect.threshold = ", effect_threshold,
          " is invalid for effect_measure = '", effect_measure,
          "' (ratio scale).  Pass a positive value on the natural ",
          "ratio scale (e.g., effect.threshold = 1.25 for ",
          effect_measure, " >= 1.25).",
          call. = FALSE
        )
      }
      if (consistency_threshold <= 0) {
        stop(
          "consistency.threshold = ", consistency_threshold,
          " is invalid for effect_measure = '", effect_measure,
          "' (ratio scale).  Pass a positive value on the natural ",
          "ratio scale (e.g., consistency.threshold = 1.0 for ",
          effect_measure, " >= 1.0 per split).",
          call. = FALSE
        )
      }
      effect_threshold      <- log(effect_threshold)
      consistency_threshold <- log(consistency_threshold)
    }

    # -----------------------------------------------------------------
    # Configuration summary (always printed via message)
    # -----------------------------------------------------------------
    is_identity <- effect_measure %in% c("RD", "IRD", "MD")
    measure_labels <- c(
      RD  = "Risk Difference", OR = "Odds Ratio", RR = "Risk Ratio",
      IRR = "Incidence Rate Ratio", IRD = "Incidence Rate Difference",
      MD  = "Mean Difference"
    )
    measure_label <- measure_labels[[effect_measure]]

    if (is_identity) {
      screen_desc <- sprintf("%.4f", effect_threshold)
      consist_desc <- sprintf("%.4f", consistency_threshold)
      if (effect_measure == "RD") {
        screen_interp <- sprintf(
          "treatment increases event probability by >= %.1f pct points",
          100 * effect_threshold)
        consist_interp <- if (consistency_threshold == 0)
          "each split shows any harm direction (RD > 0)"
        else sprintf("each split shows RD >= %.1f pct points",
                      100 * consistency_threshold)
      } else {
        screen_interp <- sprintf("%s >= %s", effect_measure, screen_desc)
        consist_interp <- sprintf("%s >= %s per split",
                                   effect_measure, consist_desc)
      }
    } else {
      # Ratio scale: show on natural scale
      screen_nat <- exp(effect_threshold)
      consist_nat <- exp(consistency_threshold)
      screen_desc <- sprintf("%.2f (log: %.4f)", screen_nat, effect_threshold)
      consist_desc <- sprintf("%.2f (log: %.4f)", consist_nat, consistency_threshold)
      screen_interp <- sprintf("subgroup %s >= %.2f",
                                effect_measure, screen_nat)
      consist_interp <- sprintf("each split %s >= %.2f",
                                 effect_measure, consist_nat)
    }

    if (!quiet && subgroup_method != "dina") {
    message(sprintf(paste0(
      "\n[forestsearch] Subgroup Identification Configuration\n",
      "  Outcome type:   %s\n",
      "  Effect measure: %s (%s)\n",
      "  Scale:          %s\n",
      "\n",
      "  Screening:      %s >= %s\n",
      "  Consistency:    %s >= %s\n",
      "  Consistency rate threshold: %.0f%%\n",
      "\n",
      "  Interpretation:\n",
      "    Candidates: %s\n",
      "    Splits:     %s\n"),
      outcome_type, effect_measure, measure_label,
      if (is_identity) "identity" else "log (ratio)",
      effect_measure, screen_desc,
      effect_measure, consist_desc,
      100 * pconsistency.threshold,
      screen_interp, consist_interp))
    }

    # Store resolved config for the return object
    threshold_config <- list(
      outcome_type = outcome_type,
      effect_measure = effect_measure,
      screening = effect_threshold,
      consistency = consistency_threshold,
      screening_natural = if (is_identity) effect_threshold
                          else exp(effect_threshold),
      consistency_natural = if (is_identity) consistency_threshold
                            else exp(consistency_threshold),
      scale = if (is_identity) "identity" else "log",
      pconsistency = pconsistency.threshold,
      screening_description = screen_interp,
      consistency_description = consist_interp
    )

    # Resolve dmin.grf for GLM outcomes.
    # The survival default for dmin.grf is on the RMST scale (event count
    # threshold over a window) and would filter out all candidates for
    # binary/continuous outcomes where DR score differences are on the
    # probability/mean scale (0--1).  When the user has not supplied a
    # value, default to 0.0 (any positive harm qualifies); when the user
    # supplied a value explicitly, respect it.
    if (missing(dmin.grf)) {
      dmin.grf <- 0.0
    }
    args_call_all <- .sync_args_call_all(args_call_all, environment(),
                                         "dmin.grf")
  } else {
    # Survival path: print configuration summary
    if (!quiet && subgroup_method != "dina") {
    message(sprintf(paste0(
      "\n[forestsearch] Subgroup Identification Configuration\n",
      "  Outcome type:   survival (Cox PH)\n",
      "  Effect measure: HR (hazard ratio)\n",
      "  Scale:          log (ratio)\n",
      "\n",
      "  Screening:      HR >= %.2f  (log: %.4f)\n",
      "  Consistency:    HR >= %.2f  (log: %.4f)\n",
      "  Consistency rate threshold: %.0f%%\n",
      "\n",
      "  Interpretation:\n",
      "    Candidates: subgroup hazard ratio >= %.2f\n",
      "    Splits:     each split HR >= %.2f\n"),
      hr.threshold, log(hr.threshold),
      hr.consistency, log(max(hr.consistency, 0.001)),
      100 * pconsistency.threshold,
      hr.threshold, hr.consistency))
    }

    # Store resolved config for the return object
    threshold_config <- list(
      outcome_type = "survival",
      effect_measure = "HR",
      screening = log(hr.threshold),
      consistency = log(max(hr.consistency, 0.001)),
      screening_natural = hr.threshold,
      consistency_natural = hr.consistency,
      scale = "log",
      pconsistency = pconsistency.threshold,
      screening_description = sprintf(
        "subgroup hazard ratio >= %.2f", hr.threshold),
      consistency_description = sprintf(
        "each split HR >= %.2f", hr.consistency)
    )
  }

  # -- Search alignment diagnostic --------------------------------------
  # Only print when details = TRUE (single-data exploratory analysis).
  # Suppressed during bootstrap, cross-validation, and simulation loops
  # where quiet = TRUE or details = FALSE.
  search_diag <- interpret_search_config(
    outcome_type         = if (outcome_type == "survival") "survival"
                           else outcome_type,
    effect_measure       = if (outcome_type == "survival") "HR"
                           else effect_measure,
    adverse_outcome      = adverse_outcome,
    effect_threshold     = if (outcome_type == "survival") log(hr.threshold)
                           else effect_threshold,
    consistency_threshold = if (outcome_type == "survival")
                              log(max(hr.consistency, 0.001))
                            else consistency_threshold,
    use_lasso            = use_lasso,
    use_grf              = use_grf,
    outcome.name         = outcome.name,
    event.name           = event.name,
    treat.name           = treat.name,
    offset.name          = offset.name,
    quiet                = !details || quiet || subgroup_method == "dina"
  )

  # ===========================================================================
  # SECTION 2C: PROPENSITY SCORE ESTIMATION (optional)
  # ===========================================================================

  is_glm <- outcome_type != "survival"
  #
  # When ps_method != "none", estimate P(W=1|X) and use either:
  #   - IPTW: stabilized weights in glm (coefficient IS the effect)
  #   - DR G-comp: Bang & Robins (2005) IPS covariate + G-computation
  #
  # Defaults:
  #   RCT = TRUE  -> ps_method = "none"  (randomization ensures balance)
  #   RCT = FALSE -> ps_method = "grf"   (non-parametric, cross-fitted)
  #   ps_adjust_method = "iptw" (default)
  # ===========================================================================

  # Resolve ps_method default
  if (is.null(ps_method)) {
    ps_method <- if (is.RCT) "none" else "grf"
  }
  args_call_all <- .sync_args_call_all(args_call_all, environment(),
                                       "ps_method")

  # Resolve ps_adjust_method.  NB: the resolved value lives in
  # ps_adjust_resolved (used downstream alongside the un-overridden
  # ps_adjust_method), so this mirror remains an explicit override
  # rather than a .sync_args_call_all() call.
  ps_adjust_method <- match.arg(ps_adjust_method)
  ps_adjust_resolved <- if (ps_method == "none") "none" else ps_adjust_method
  args_call_all$ps_adjust_method <- ps_adjust_resolved

  # Regression/strata covariate adjustment (adjust_covariates) and
  # propensity-score adjustment are mutually exclusive for now: the PS branch
  # below rebuilds the GLM estimator closure around ps_adjust_method only and
  # does not also fold in regression covariates.  Guard rather than silently
  # drop the user's adjust_covariates.  (Do NOT reset adjust_covariates here;
  # it is a user-facing formal that must flow to df.fs re-attachment and the
  # consistency engine.)
  if (!is.null(adjust_covariates) && ps_method != "none") {
    stop("adjust_covariates cannot be combined with ps_method != 'none'. ",
         "Regression/strata adjustment and propensity-score adjustment are ",
         "mutually exclusive; choose one.", call. = FALSE)
  }

  if (ps_method != "none") {
    if (!is.null(ps_hat)) {
      # User-supplied propensity scores
      if (length(ps_hat) != nrow(df.analysis)) {
        stop("ps_hat must have length equal to nrow(df.analysis)", call. = FALSE)
      }
      # Reject values that would produce Inf / NaN IPTW weights downstream.
      # The internal estimate_propensity_scores() path trims; the user-
      # supplied path bypasses that protection, so guard explicitly.
      if (anyNA(ps_hat)) {
        stop("ps_hat contains NA values; remove or impute before passing.",
             call. = FALSE)
      }
      if (any(ps_hat <= 0 | ps_hat >= 1)) {
        n_bad <- sum(ps_hat <= 0 | ps_hat >= 1)
        stop("ps_hat must lie strictly in (0, 1); ", n_bad,
             " value(s) outside this range would produce Inf/NaN weights. ",
             "Trim or clip ps_hat (e.g., pmin(pmax(ps_hat, 0.01), 0.99)) ",
             "before passing.",
             call. = FALSE)
      }
      df.analysis$ps_hat <- ps_hat
      W_vec <- df.analysis[[treat.name]]
      p_treat <- mean(W_vec, na.rm = TRUE)
      df.analysis$sw <- ifelse(W_vec == 1,
                               p_treat / ps_hat,
                               (1 - p_treat) / (1 - ps_hat))
      df.analysis$ips_covar <- ifelse(W_vec == 1,
                                      1 / ps_hat,
                                      1 / (1 - ps_hat))
    } else {
      # Estimate PS using the requested method
      ps_result <- estimate_propensity_scores(
        data             = df.analysis,
        treat.name       = treat.name,
        confounders.name = confounders.name,
        method           = ps_method,
        seed             = seedit
      )
      df.analysis <- ps_result$data
      if (details) {
        cat("  PS estimation: method =", ps_result$method,
            ", adjust =", ps_adjust_resolved,
            ", trimmed =", ps_result$trimmed, "\n")
      }
    }

    # Rebuild estimator closure with PS adjustment (GLM outcomes only)
    if (is_glm) {
      estimator_fn <- make_effect_estimator(
        outcome_type     = outcome_type,
        treat.name       = treat.name,
        outcome.name     = outcome.name,
        offset.name      = offset.name,
        effect_measure   = effect_measure,
        ps_adjust_method = ps_adjust_resolved
      )
    }
  }

  args_call_all <- .sync_args_call_all(args_call_all, environment(),
                                       "adjust_covariates")

  # ===========================================================================
  # SECTION 3: INITIALIZE TIMING AND DATA
  # ===========================================================================

  t.start_all <- proc.time()[3]

  df <- df.analysis
  grf_plot <- NULL

  if (details && !is.null(conf_force)) {
    cat("Forced confounders:", paste(conf_force, collapse = ", "), "\n")
  }

  # Advisory note: for the effect-band foci, the criterion-optimal subgroup
  # can rank below a small candidate cap, so a broad pool is recommended.
  if (details &&
      .normalize_sg_focus(sg_focus) %in% c("hrMaxSG", "hrMinSG") &&
      max_subgroups_search < 30) {
    cat("Note: sg_focus = '", sg_focus, "' with max_subgroups_search = ",
        max_subgroups_search,
        ". For effMaxSG/effMinSG the optimal subgroup may rank below the cap; ",
        "increasing max_subgroups_search to >= 30 may be necessary.\n",
        sep = "")
  }

  if (details && length(conf.cont_jcuts) > 0L) {
    cat("J-quantile cuts:",
        paste0(names(conf.cont_jcuts),
               " (J=", unlist(conf.cont_jcuts), ")",
               collapse = ", "),
        "\n")
  }

  # ===========================================================================
  # SECTION 3-DINA: DINA-SELECTION MODE (subgroup_method = "dina")
  # ===========================================================================
  # Delegate subgroup *selection* to dina_subgroup(), bypassing GRF, LASSO and
  # the consistency search.  The selected subgroup is returned as sg.harm and
  # consumed by the existing (label / treat.recommend-based) estimation and
  # bootstrap machinery.  Selection criteria are inherited from this call's
  # sg_focus / selection_rule / effect_neighborhood / n.min / hr.threshold;
  # dina_args supplies only the DINA fit tuning.
  if (subgroup_method == "dina") {
    dsel <- .forestsearch_dina_select(
      df = df, df.predict = df.predict, df.test = df.test,
      confounders.name = confounders.name, outcome.name = outcome.name,
      event.name = event.name, treat.name = treat.name, id.name = id.name,
      outcome_type = outcome_type, hr.threshold = hr.threshold, n.min = n.min,
      sg_focus = sg_focus, selection_rule = selection_rule,
      effect_neighborhood = effect_neighborhood, dina_args = dina_args,
      dina_res = dina_res, seedit = seedit, details = details,
      effect_measure = effect_measure, offset.name = offset.name,
      adjust_covariates = adjust_covariates,
      adverse_outcome = if (identical(outcome_type, "survival")) TRUE
                        else adverse_outcome)

    t.min_all <- (proc.time()[3] - t.start_all) / 60

    if (details) {
      if (isTRUE(dsel$found)) {
        cat("DINA-selected subgroup:",
            paste(dsel$sg.harm, collapse = " & "), "\n")
      } else {
        cat("DINA selection: no subgroup met the criteria\n")
      }
    }

    out <- list(
      grp.consistency       = dsel$grp.consistency,
      find.grps             = NULL,
      confounders.candidate = confounders.name,
      confounders.evaluated = confounders.name,
      df.est                = dsel$df.est,
      df.predict            = dsel$df.predict,
      df.test               = dsel$df.test,
      minutes_all           = t.min_all,
      grf_res               = NULL,
      sg_focus              = sg_focus,
      sg.harm               = dsel$sg.harm,
      grf_cuts              = NULL,
      dina_res              = dsel$dina_res,
      dina_cuts             = NULL,
      prop_maxk             = NA_real_,
      max_sg_est            = NA_real_,
      grf_plot              = NULL,
      args_call_all         = args_call_all,
      consistency_algorithm = "dina",
      outcome_type          = outcome_type,
      effect_measure        = effect_measure,
      threshold_config      = if (exists("threshold_config")) threshold_config
                              else NULL,
      subgroup_method       = "dina"
    )

    # SECTION 9B (DINA): Tier-2 de-biased gate over DINA's OWN candidate family
    # (approach (i)).  Re-selection matches the DINA focus rule via sg_focus.
    # Defensive: any failure leaves out$debias_gate NULL; DINA output unaffected.
    out$debias_gate <- NULL
    if (isTRUE(debias_gate) && isTRUE(dsel$found) &&
        !is.null(dsel$grp.consistency) &&
        !is.null(dsel$grp.consistency$sg.harm.id)) {
      .dg_df   <- dsel$df.est
      # Mirror SECTION 9B: survival -> effect_measure "HR" on log scale;
      # GLM -> resolved effect_measure with the precomputed comparison-scale
      # threshold (effect_measure is NULL for survival, so do NOT test it).
      .dg_em   <- if (identical(outcome_type, "survival")) "HR" else effect_measure
      .dg_cscr <- if (identical(outcome_type, "survival")) log(hr.threshold)
                  else effect_threshold
      .dg_spec <- .fs_dg_spec(
        outcome_type, .dg_em, treat.name, outcome.name, event.name,
        offset.name, adjust_covariates,
        adverse_outcome = if (identical(outcome_type, "survival")) TRUE
                          else adverse_outcome,
        df = .dg_df)
      # Restrict the gate family to candidates competitive on the SAME statistic
      # that selected the winner, so the re-selection mirrors DINA's selection
      # neighbourhood instead of the full family.  With select_statistic =
      # "effect", DINA selected on the inferential effect, so restrict on the
      # per-candidate effect (`sel_effect`, attached by the select wrapper);
      # otherwise restrict on DINA's native tau-hat.  Default neighbourhood:
      # match effect_neighborhood (the band DINA selects within).  Set
      # debias_gate_args$family_native_neighborhood >= 1 to disable (full family).
      .dg_nbhd <- debias_gate_args$family_native_neighborhood
      if (is.null(.dg_nbhd)) .dg_nbhd <- effect_neighborhood
      .dg_stat <- if (identical(dsel$select_statistic, "effect") &&
                      "sel_effect" %in% names(dsel$candidates))
                    dsel$candidates$sel_effect else dsel$candidates$tau_hat
      .dg_tab  <- .fs_dg_restrict_native(
        dsel$candidates, .dg_stat,
        neighborhood = .dg_nbhd,
        log_scale = .dg_em %in% c("HR", "OR", "IRR"))
      .dg_fam  <- .fs_dg_family_from_table(.dg_df, .dg_tab,
                                           op_right = ">=", n_min = n.min)
      out$debias_gate <- .fs_apply_debias_gate(
        df = .dg_df, candidates = .dg_fam,
        selected_members = which(dsel$grp.consistency$sg.harm.id == 1L),
        spec = .dg_spec, c_screen = .dg_cscr, c_consistency = 0,
        p_star = pconsistency.threshold,
        effect_neighborhood = effect_neighborhood,
        reselection_default = .fs_dg_reselection_from_focus(sg_focus, engine = "effect"),
        selection_rule_default = selection_rule,
        debias_gate_args = debias_gate_args, seedit = seedit)
    }

    out$harm_flag_debiased <- if (!is.null(out$debias_gate))
      isTRUE(out$debias_gate$harm_flag) else NA
    class(out) <- c("forestsearch", "list")
    return(out)
  }

  # ===========================================================================
  # SECTION 3-GRF: GRF-SELECTION MODE (subgroup_method = "grf")
  # ===========================================================================
  # Delegate subgroup *selection* to the GRF subgroup-identification routine
  # (grf.subg.harm.survival / grf.subg.harm.glm), bypassing LASSO/DINA
  # screening and the consistency search.  The GRF policy-tree harm leaf is
  # returned as sg.harm and consumed by the existing (label / treat.recommend-
  # based) estimation and bootstrap machinery.  The GRF run reuses this call's
  # existing GRF arguments (frac.tau, dmin.grf, grf_depth, is.RCT, seedit, and
  # the GLM extras); there is no separate grf_args list, mirroring how GRF
  # screening is configured.
  if (subgroup_method == "grf") {
    grf_selection <- match.arg(grf_selection)
    grf_select_statistic <- match.arg(grf_select_statistic)
    # In frontier mode the selection rule comes from sg_focus (which the tree
    # path ignores).  All five sg_focus values map to the aligned frontier rule;
    # an unrecognized value falls back to the robust default (effMaxSG).
    .sgf <- tryCatch(.normalize_sg_focus(sg_focus), error = function(e) sg_focus)
    frontier_rule <- switch(as.character(.sgf),
                            hr = "eff", hrMaxSG = "effMaxSG", hrMinSG = "effMinSG",
                            maxSG = "maxSG", minSG = "minSG",
                            "effMaxSG")
    gsel <- .forestsearch_grf_select(
      df = df, df.predict = df.predict, df.test = df.test,
      confounders.name = confounders.name, outcome.name = outcome.name,
      event.name = event.name, treat.name = treat.name, id.name = id.name,
      outcome_type = outcome_type, n.min = n.min,
      frac.tau = frac.tau, dmin.grf = dmin.grf, grf_depth = grf_depth,
      is.RCT = is.RCT, adverse_outcome = adverse_outcome,
      offset.name = offset.name, overdispersion = overdispersion,
      grf_count_transform = grf_count_transform,
      grf_res = grf_res, seedit = seedit,
      grf_selection = grf_selection, frontier_rule = frontier_rule,
      effect_neighborhood = effect_neighborhood, details = details,
      grf_select_statistic = grf_select_statistic,
      effect_measure = effect_measure, adjust_covariates = adjust_covariates)

    t.min_all <- (proc.time()[3] - t.start_all) / 60

    if (details) {
      if (isTRUE(gsel$found)) {
        cat("GRF-selected subgroup:",
            paste(gsel$sg.harm, collapse = " & "), "\n")
      } else {
        cat("GRF selection: no subgroup met the criteria\n")
      }
    }

    out <- list(
      grp.consistency       = gsel$grp.consistency,
      find.grps             = NULL,
      confounders.candidate = confounders.name,
      confounders.evaluated = confounders.name,
      df.est                = gsel$df.est,
      df.predict            = gsel$df.predict,
      df.test               = gsel$df.test,
      minutes_all           = t.min_all,
      grf_res               = gsel$grf_res,
      sg_focus              = sg_focus,
      sg.harm               = gsel$sg.harm,
      grf_cuts              = gsel$sg.harm,
      dina_res              = NULL,
      dina_cuts             = NULL,
      prop_maxk             = NA_real_,
      max_sg_est            = NA_real_,
      grf_plot              = NULL,
      args_call_all         = args_call_all,
      consistency_algorithm = "grf",
      outcome_type          = outcome_type,
      effect_measure        = effect_measure,
      threshold_config      = if (exists("threshold_config")) threshold_config
                              else NULL,
      subgroup_method       = "grf"
    )

    # SECTION 9B (GRF): Tier-2 de-biased gate over GRF's OWN DR-candidate family
    # (approach (i)).  Frontier ranks this family directly; the policy tree uses
    # it as the cut-defined universe.  Re-selection matches sg_focus.  GRF cut
    # operator for non-"left" is ">" (matches .grf_sg_def_from_candidate).
    # Defensive: any failure leaves out$debias_gate NULL; GRF output unaffected.
    out$debias_gate <- NULL
    if (isTRUE(debias_gate) && !is.null(gsel$grp.consistency) &&
        !is.null(gsel$grp.consistency$sg.harm.id)) {
      .dg_df   <- gsel$df.est
      # Mirror SECTION 9B (see DINA branch): survival uses "HR"/log(hr.threshold);
      # GLM uses effect_measure/effect_threshold.  effect_measure is NULL for
      # survival, so never test it directly.
      .dg_em   <- if (identical(outcome_type, "survival")) "HR" else effect_measure
      .dg_cscr <- if (identical(outcome_type, "survival")) log(hr.threshold)
                  else effect_threshold
      .dg_spec <- .fs_dg_spec(
        outcome_type, .dg_em, treat.name, outcome.name, event.name,
        offset.name, adjust_covariates,
        adverse_outcome = if (identical(outcome_type, "survival")) TRUE
                          else adverse_outcome,
        df = .dg_df)
      # Restrict the gate family to candidates competitive on the SAME statistic
      # that selected the winner.  With grf_select_statistic = "effect" (frontier
      # mode), restrict on the per-candidate effect (`sel_effect`, link scale,
      # attached by the select wrapper); otherwise restrict on GRF's native
      # doubly-robust score (additive scale).  Default: match effect_neighborhood
      # (the band GRF's frontier selects within).  Set
      # debias_gate_args$family_native_neighborhood >= 1 to disable (full family).
      .dg_nbhd <- debias_gate_args$family_native_neighborhood
      if (is.null(.dg_nbhd)) .dg_nbhd <- effect_neighborhood
      .dg_eff_sel <- identical(gsel$select_statistic, "effect") &&
                     "sel_effect" %in% names(gsel$candidates)
      .dg_stat    <- if (.dg_eff_sel) gsel$candidates$sel_effect
                     else gsel$candidates$effect
      .dg_tab  <- .fs_dg_restrict_native(
        gsel$candidates, .dg_stat,
        neighborhood = .dg_nbhd,
        log_scale = if (.dg_eff_sel) .dg_em %in% c("HR", "OR", "IRR") else FALSE)
      .dg_fam  <- .fs_dg_family_from_table(.dg_df, .dg_tab,
                                           op_right = ">", n_min = n.min)
      out$debias_gate <- .fs_apply_debias_gate(
        df = .dg_df, candidates = .dg_fam,
        selected_members = which(gsel$grp.consistency$sg.harm.id == 1L),
        spec = .dg_spec, c_screen = .dg_cscr, c_consistency = 0,
        p_star = pconsistency.threshold,
        effect_neighborhood = effect_neighborhood,
        reselection_default = .fs_dg_reselection_from_focus(sg_focus, engine = "effect"),
        selection_rule_default = selection_rule,
        debias_gate_args = debias_gate_args, seedit = seedit)
    }

    out$harm_flag_debiased <- if (!is.null(out$debias_gate))
      isTRUE(out$debias_gate$harm_flag) else NA
    class(out) <- c("forestsearch", "list")
    return(out)
  }

  # ===========================================================================
  # SECTION 3A: GRF CUT GENERATION (if use_grf = TRUE)
  # ===========================================================================

  # If using grf and cuts not already populated, run GRF subgroup identification
  if (use_grf && (is.null(grf_res) || is.null(grf_res$tree.cuts))) {

    grf_res <- tryCatch({
      # Argument-list construction is delegated to the internal builders
      # in R/grf_args.R so the public extract_grf_args() and the call
      # site here cannot drift.
      if (outcome_type == "survival") {
        surv_args <- .build_grf_survival_args(
          data                      = df.analysis,
          confounders.name          = confounders.name,
          outcome.name              = outcome.name,
          event.name                = event.name,
          id.name                   = id.name,
          treat.name                = treat.name,
          frac.tau                  = frac.tau,
          n.min                     = n.min,
          dmin.grf                  = dmin.grf,
          is.RCT                    = is.RCT,
          grf_depth                 = grf_depth,
          seedit                    = seedit,
          return_selected_cuts_only = return_selected_cuts_only
        )
        do.call(grf.subg.harm.survival, surv_args)
      } else {
        glm_args <- .build_grf_glm_args(
          data                      = df.analysis,
          confounders.name          = confounders.name,
          outcome.name              = outcome.name,
          treat.name                = treat.name,
          id.name                   = id.name,
          outcome_type              = outcome_type,
          n.min                     = n.min,
          dmin.grf                  = dmin.grf,
          is.RCT                    = is.RCT,
          grf_depth                 = grf_depth,
          seedit                    = seedit,
          return_selected_cuts_only = return_selected_cuts_only,
          adverse_outcome           = adverse_outcome,
          offset.name               = offset.name,
          overdispersion            = overdispersion,
          grf_count_transform       = grf_count_transform
        )
        do.call(grf.subg.harm.glm, glm_args)
      }
    },
      error = function(e) {
        warning("GRF analysis failed: ", e$message)
        return(NULL)
      }
    )

    # Extract GRF cuts if successful
    if (!is.null(grf_res) && !inherits(grf_res, "try-error")) {
      grf_cuts <- grf_res$tree.cuts

      # GLM GRF returns tree.cuts as a named list of numeric values
      # (e.g., list(bm1 = 1, age = c(70, 35))), but get_FSdata expects
      # a character vector of cut expressions ("bm1 <= 1", "age <= 70").
      # Convert list format to character expression format.
      if (is.list(grf_cuts) && !is.null(names(grf_cuts))) {
        grf_cuts <- unlist(lapply(names(grf_cuts), function(nm) {
          paste0(nm, " <= ", grf_cuts[[nm]])
        }))
      }

      if (details) {
        # Concise GRF summary: subgroup found + cuts.  tidy_cut_display() is
        # display-only (collapses padding, rounds thresholds to integer to match
        # collapse_cuts); the underlying grf_cuts / sg.harm.id are unchanged.
        grf_sg <- grf_res$sg.harm.id
        if (!is.null(grf_sg) && length(grf_sg) > 0 &&
            !all(is.na(grf_sg)) && !all(grf_sg == "")) {
          cat("GRF subgroup:",
              tidy_cut_display(paste(grf_sg, collapse = " & ")), "\n")
        } else {
          cat("GRF: no subgroup identified\n")
        }
        cat("GRF cuts identified:", length(grf_cuts), "\n")
        if (length(grf_cuts) > 0) {
          cat("  Cuts:", paste(tidy_cut_display(grf_cuts), collapse = ", "), "\n")
        }
      }

      # Generate GRF plot if requested
      if (plot.grf && !is.null(grf_res$tree)) {
        grf_plot <- tryCatch({
          plot(grf_res$tree)
        }, error = function(e) {
          warning("GRF plot failed: ", e$message)
          NULL
        })
      }
    } else {
      grf_cuts <- NULL
      if (details) {
        cat("GRF analysis did not produce cuts\n")
      }
    }
  }

  # ===========================================================================
  # SECTION 3B: DINA CUT GENERATION (if use_dina = TRUE)
  # ===========================================================================
  # Mirrors SECTION 3A: fit a DINA model on the analysis data and extract a
  # per-covariate Pareto frontier as additional screening-stage candidate
  # cuts ("x1 <= 0.5"), merged into the candidate pool by get_FSdata()
  # exactly as the GRF cuts are.  These are candidates, not forced
  # selections; the consistency stage remains the gatekeeper.

  if (use_dina && is.null(dina_cuts)) {

    dina_cuts <- tryCatch({
      da <- .resolve_dina_args(dina_args, outcome_type,
                               n_min_default = n.min,
                               seed_default  = seedit)

      if (is.null(dina_res)) {
        # Cox requires a status column; GLM families do not.
        status_arg <- if (identical(da$fit$family, "cox")) event.name else NULL
        fit_call <- c(
          list(df = df.analysis, outcome = outcome.name,
               treatment = treat.name, covariates = confounders.name,
               status = status_arg),
          da$fit
        )
        dina_res <- do.call(dina, fit_call)
      }

      fr <- do.call(dina_frontier,
                    c(list(fit = dina_res, df = df.analysis,
                           covariates = confounders.name),
                      da$frontier))

      if (isTRUE(da$selected_only)) {
        # Selected-cut screening: contribute the SINGLE cut dina_subgroup()
        # selects, using forestsearch's own selection criteria (so this
        # matches what subgroup_method = "dina" would choose).  m_diff is the
        # harm floor on the link scale: log(hr.threshold) for ratio families
        # (cox/binomial/poisson), identity for gaussian.
        m_diff_sel <- if (identical(da$fit$family, "gaussian")) hr.threshold
                      else log(hr.threshold)
        sgsel <- dina_subgroup(
          fit                 = dina_res,
          df                  = df.analysis,
          covariates          = confounders.name,
          m_diff              = m_diff_sel,
          n_min               = n.min,
          max_depth           = da$select$max_depth,
          grid_probs          = da$select$grid_probs,
          sg_focus            = sg_focus,
          selection_rule      = selection_rule,
          effect_neighborhood = effect_neighborhood
        )
        if (isTRUE(sgsel$found)) {
          # Vectorized over the 1-or-2 selected cuts: a depth-2 conjunction
          # contributes both component cuts to the screening pool, which the
          # consistency search composes as usual.
          op_sel <- ifelse(sgsel$direction == "left", "<=", ">=")
          paste0(sgsel$covariate, " ", op_sel,
                 " ", signif(sgsel$threshold, da$frontier$digits))
        } else {
          character(0)
        }
      } else {
        fr$cut_expr
      }
    },
      error = function(e) {
        warning("DINA analysis failed: ", e$message)
        return(NULL)
      }
    )

    if (details) {
      da_mode <- .resolve_dina_args(dina_args, outcome_type,
                                    n_min_default = n.min, seed_default = seedit)
      mode_lbl <- if (isTRUE(da_mode$selected_only))
                    "selected cut (dina_subgroup, forestsearch criteria)"
                  else "frontier candidates"
      n_dc <- length(dina_cuts)
      cat("DINA screening mode:", mode_lbl, "\n")
      cat("DINA cuts identified:", n_dc, "\n")
      if (n_dc > 0) {
        cat("  Cuts:", paste(dina_cuts, collapse = ", "), "\n")
      }
    }
  }

  # ===========================================================================
  # SECTION 4: DATA PREPARATION (get_FSdata)
  # ===========================================================================

  # For GLM outcomes, get_FSdata validates that event.name is numeric.
  # Map it to something valid so the factor-construction machinery works.
  if (outcome_type == "binary") {
    # Binary outcome IS the event indicator
    if (!event.name %in% names(df.analysis)) {
      event.name <- outcome.name
    }
  } else if (outcome_type %in% c("continuous", "count")) {
    # Continuous / count outcomes have no censoring event; create a placeholder
    if (!event.name %in% names(df.analysis)) {
      df.analysis[[".event_placeholder"]] <- 1L
      event.name <- ".event_placeholder"
      df <- df.analysis
    }
  }

  # Update args_call_all with resolved event.name for downstream calls
  args_call_all <- .sync_args_call_all(args_call_all, environment(),
                                       "event.name")

  FSdata <- tryCatch(
    do.call(
      get_FSdata,
      filter_call_args(
        args_call_all,
        get_FSdata,
        list(df.analysis = df.analysis, grf_cuts = grf_cuts,
             dina_cuts = dina_cuts,
             details = FALSE)
      )
    ),
    error = function(e) {
      warning("Error in get_FSdata: ", e$message)
      return(NULL)
    }
  )

  if (inherits(FSdata, "try-error") || is.null(FSdata)) {
    warning("FSdata failure - returning NULL result")
    return(list(sg.harm       = NULL,
                args_call_all = args_call_all,
                error_log     = list(stage   = "get_FSdata",
                                     message = "FSdata returned NULL or error")))
  }

  # Extract FSdata components
  lassoomit <- FSdata$lassoomit
  lassokeep <- FSdata$lassokeep
  df <- FSdata$df

  # Concise LASSO summary
  if (details && use_lasso) {
    n_kept <- length(lassokeep)
    n_total <- n_kept + length(lassoomit)
    cat("Cox-LASSO selected:", n_kept, "of", n_total, "candidate factors\n")
    if (length(lassoomit) > 0) {
      cat("  Omitted:", paste(lassoomit, collapse = ", "), "\n")
    }
  }

  Y <- df[, outcome.name]
  if (outcome_type == "survival") {
    Event <- df[, event.name]
  } else {
    # For GLM outcomes, Event is not meaningful but some downstream
    # code references it.  Use a vector of 1s as placeholder.
    Event <- rep(1L, nrow(df))
  }
  Treat <- df[, treat.name]

  FSconfounders.name <- FSdata$confs_names
  confs_labels <- FSdata$confs

  if (details) {
    cat("Candidate factors:", length(confs_labels), "\n")
    print(confs_labels)
  }

  if (is.null(df.predict)) df.predict <- df

  # ===========================================================================
  # SECTION 5: GRF VARIABLE IMPORTANCE SCREENING
  # ===========================================================================

  if (!is.null(vi.grf.min)) {
    X <- apply(df[, FSconfounders.name, drop = FALSE], 2, as.numeric)

    if (outcome_type == "survival") {
      # Survival path: causal_survival_forest (unchanged)
      # Guard against zero events in either arm: max(numeric(0)) is -Inf
      # (with a warning), which would propagate to horizon = frac.tau * -Inf
      # and crash causal_survival_forest.  Skip the GRF VI step instead.
      y_evt_t <- Y[Treat == 1 & Event == 1]
      y_evt_c <- Y[Treat == 0 & Event == 1]

      if (length(y_evt_t) == 0L || length(y_evt_c) == 0L) {
        if (!quiet) {
          message("[forestsearch] Skipping GRF variable-importance screening: ",
                  "zero events in ",
                  if (length(y_evt_t) == 0L) "treatment" else "control",
                  " arm.")
        }
        cs.forest <- structure("zero-events", class = "try-error")
      } else {
        tau.rmst <- min(max(y_evt_t), max(y_evt_c))
        tune_arg <- if (tune_grf) "all" else "none"

        if (!is.RCT) {
          cs.forest <- try(suppressWarnings(
            grf::causal_survival_forest(X, Y, Treat, Event,
                                         horizon = frac.tau * tau.rmst, seed = seedit,
                                         tune.parameters = tune_arg)
          ), TRUE)
        } else {
          cs.forest <- try(suppressWarnings(
            grf::causal_survival_forest(X, Y, Treat, Event, W.hat = 0.5,
                                         horizon = frac.tau * tau.rmst, seed = seedit,
                                         tune.parameters = tune_arg)
          ), TRUE)
        }
      }
    } else {
      # GLM path: causal_forest (no event/horizon)
      if (outcome_type == "count") {
        # Count data: variance-stabilising log transform (same as
        # grf.subg.harm.glm with grf_count_transform = "log")
        Y_vi <- if (grf_count_transform == "log") log(Y + 0.5) else Y
      } else if (outcome_type == "binary") {
        # Binary: complement (1-Y) so positive CATE = treatment helps
        Y_vi <- if (adverse_outcome) 1L - Y else Y
      } else {
        # Continuous: negate (-Y) so positive CATE = treatment helps
        Y_vi <- if (adverse_outcome) -Y else Y
      }
      cs.forest <- try(suppressWarnings(
        fit_causal_forest_glm(X, Y_vi, Treat, is.RCT, seedit = seedit,
                              tune_grf = tune_grf)
      ), TRUE)
    }

    if (!inherits(cs.forest, "try-error")) {
      vi.cs <- round(grf::variable_importance(cs.forest), 4)
      vi.cs2 <- data.frame(confs_labels, FSconfounders.name, vi.cs)
      vi.order <- order(vi.cs, decreasing = TRUE)
      vi.cs2 <- vi.cs2[vi.order, ]

      conf.screen <- vi.cs2[, 2]
      vi_max <- max(vi.cs2[, 3])
      if (vi_max > 0) {
        vi_ratio <- vi.cs2[, 3] / vi_max
        selected.vars <- which(vi_ratio > vi.grf.min)
        conf.screen <- conf.screen[selected.vars]
        conf.screen <- conf.screen[seq_len(min(length(conf.screen), max_n_confounders))]
      }
      # Fall back to full set if GRF screening eliminated everything
      if (length(conf.screen) == 0L) {
        conf.screen <- FSconfounders.name
      }
    } else {
      # GRF fit failed: skip screening, use all confounders
      conf.screen <- FSconfounders.name
    }
  } else {
    conf.screen <- FSconfounders.name
  }

  # ===========================================================================
  # SECTION 6: SUBGROUP SEARCH
  # ===========================================================================

  # Prepare data for subgroup search
  # NOTE: df[, conf.screen] columns are factors with levels {0, 1}.
  # dummy() expands each 2-level factor into two indicator columns
  # (e.g., q1 -> q1.0, q1.1), so the search explores BOTH directions
  # of each cut (subgroup and complement).

  df.confounders <- df[, conf.screen, drop = FALSE]
  df.confounders <- dummy(df.confounders)


  id <- df[, c(id.name)]
  df.fs <- data.frame(Y, Event, Treat, id, df.confounders)
  Z <- as.matrix(df.confounders)
  colnames(Z) <- names(df.confounders)

  # For GLM outcomes: the estimator closure expects columns with the original
  # names (e.g., "treat", "response"), not the renamed Y/Event/Treat.
  # Attach the original columns so the closure can find them.
  if (!is.null(estimator_fn)) {
    if (!outcome.name %in% names(df.fs) && outcome.name %in% names(df)) {
      df.fs[[outcome.name]] <- df[[outcome.name]]
    }
    if (!treat.name %in% names(df.fs) && treat.name %in% names(df)) {
      df.fs[[treat.name]] <- df[[treat.name]]
    }
    if (!is.null(offset.name) &&
        !offset.name %in% names(df.fs) && offset.name %in% names(df)) {
      df.fs[[offset.name]] <- df[[offset.name]]
    }
  }

  # Covariate adjustment (Cox / survival path): carry the *raw* columns
  # referenced by adjust_covariates onto df.fs so the adjusted scoring
  # formula resolves inside the consistency engine.  df.fs otherwise keeps
  # only Y/Event/Treat/id plus dummied candidate factors, and strata()/
  # pspline() terms need the original (un-dichotomised) column.  The split
  # frames built downstream are row subsets of df.fs and inherit these.
  if (!is.null(adjust_covariates)) {
    adj_vars <- .fs_adjust_vars(adjust_covariates)
    reserved <- c("Y", "Event", "Treat", "id")
    bad_reserved <- intersect(adj_vars, reserved)
    if (length(bad_reserved) > 0L) {
      stop("adjust_covariates may not reference reserved column name(s): ",
           paste(bad_reserved, collapse = ", "),
           call. = FALSE)
    }
    missing_adj <- setdiff(adj_vars, names(df))
    if (length(missing_adj) > 0L) {
      stop("adjust_covariates references column(s) not in the analysis data: ",
           paste(missing_adj, collapse = ", "),
           call. = FALSE)
    }
    for (.av in adj_vars) {
      df.fs[[.av]] <- df[[.av]]
    }
    if (details) {
      .adj_model <- if (identical(outcome_type, "survival")) "Cox" else "GLM"
      cat("  Covariate adjustment (", .adj_model, "): ",
          paste(.fs_adjust_terms(adjust_covariates), collapse = " + "), "\n",
          sep = "")
    }
  }


  search_overrides <- list(
        Y = Y,
        Event = Event,
        Treat = Treat,
        Z = Z,
        d0.min = d0.min,
        d1.min = d1.min,
        n.min = n.min,
        hr.threshold = if (!is.null(effect_threshold)) effect_threshold else hr.threshold,
        max.minutes = max.minutes,
        minp = minp,
        maxk = maxk,
        parallel_workers = 1,  # Ensure sequential execution if inside parallel context
        details = details       # Optionally suppress details
  )

  # Only pass GLM params when active (same rationale as consistency_overrides)
  if (!is.null(estimator_fn)) {
    search_overrides$estimator_fn <- estimator_fn
    search_overrides$df_analysis <- df
    search_overrides$effect_threshold <- effect_threshold
    search_overrides$effect_measure <- effect_measure
    search_overrides$outcome_type <- outcome_type
  }

  # Survival covariate adjustment: pass the adjustment terms and a frame that
  # carries the raw covariate columns (df.fs, row-aligned with Y/Event/Treat/Z)
  # so the candidate-search Cox scorer can adjust consistently with the
  # consistency engine.  NULL on the GLM path (handled above).
  if (is.null(estimator_fn) && !is.null(adjust_covariates)) {
    search_overrides$adjust_covariates <- adjust_covariates
    search_overrides$df_analysis <- df.fs
  }

  # Merge and filter arguments
  search_args <- modifyList(args_call_all, search_overrides)
  # Optionally, filter only valid arguments for subgroup.search
  valid_args <- names(formals(subgroup.search))
  search_args <- search_args[names(search_args) %in% valid_args]

  find.grps <- tryCatch(
    do.call(subgroup.search, search_args),
    error = function(e) {
      warning("Error in subgroup.search: ", e$message)
      return(NULL)
    }
  )
  if (is.null(find.grps) || inherits(find.grps, "try-error")) {
    warning("Subgroup search failed")
    return(list(sg.harm = NULL, args_call_all = args_call_all,
                error_log = list(stage = "subgroup.search",
                                 message = "Search returned NULL or error")))
  }


  # ===========================================================================
  # SECTION 7: INITIALIZE OUTPUT VARIABLES
  # ===========================================================================

  sg.harm <- NULL
  df.est_out <- NULL
  df.predict_out <- NULL
  df.test_out <- NULL
  grp.consistency <- NULL

  max_sg_est <- find.grps$max_sg_est
  prop_maxk <- find.grps$prop_max_count

  t.end_all <- proc.time()[3]
  t.min_all <- (t.end_all - t.start_all) / 60

  # ===========================================================================
  # SECTION 8: CHECK FOR VALID SUBGROUPS
  # ===========================================================================

  has_subgroups <- FALSE

  if (!is.null(find.grps) &&
      !inherits(find.grps, "try-error") &&
      !is.null(find.grps$out.found) &&
      !is.null(find.grps$out.found$hr.subgroups)) {

    hr_values <- find.grps$out.found$hr.subgroups$HR
    # For survival, check HR > hr.consistency (natural scale)
    # For GLM, check effect > consistency_threshold (already on correct scale)
    check_threshold <- if (!is.null(consistency_threshold)) {
      consistency_threshold
    } else {
      hr.consistency
    }
    has_subgroups <- any(hr_values > check_threshold, na.rm = TRUE)
  }

  # ===========================================================================
  # SECTION 9: CONSISTENCY EVALUATION (WITH TWO-STAGE OPTION)
  # ===========================================================================

  if (has_subgroups) {
    # Set plotting parameter if needed
    if (plot.sg && is.null(by.risk)) {
      by.risk <- round(max(Y) / 12, 0)
    }

    if (details) {
      n_candidates <- nrow(find.grps$out.found$hr.subgroups)
      cat("# of candidate subgroups (meeting all criteria) =", n_candidates, "\n")
    }

    # -------------------------------------------------------------------------
    # Build arguments for subgroup.consistency
    # -------------------------------------------------------------------------

    # Base override arguments
    consistency_overrides <- list(
      df = df.fs,
      hr.subgroups = find.grps$out.found$hr.subgroups,
      Lsg = find.grps$L,
      confs_labels = confs_labels,
      n.splits = fs.splits,
      stop_Kgroups = max_subgroups_search,
      # NEW: Pass two-stage parameters
      use_twostage = use_twostage,
      twostage_args = twostage_args,
      consistency_method = consistency_method,
      # Covariate adjustment for Cox subgroup scoring (NULL on GLM path)
      adjust_covariates = if (is.null(estimator_fn)) adjust_covariates else NULL
    )

    # Only pass GLM closure params when active (non-NULL).
    # This ensures the survival path never sends unknown args to
    # subgroup.consistency -- critical for compatibility with
    # the CRAN-installed version.
    if (!is.null(estimator_fn)) {
      consistency_overrides$estimator_fn <- estimator_fn
      consistency_overrides$consistency_threshold <- consistency_threshold
      # The HR column in hr.subgroups contains log-scale values for GLM
      # outcomes; hr.threshold must be on the same scale so the re-filtering
      # inside subgroup.consistency (line: found.hrs$HR >= hr.threshold) works.
      consistency_overrides$hr.threshold <- effect_threshold
      # Display labels: use the resolved effect_measure name and flag
      # log-scale so the candidate table exponentiates for display.
      consistency_overrides$effect_label <- effect_measure
      consistency_overrides$effect_log_scale <-
        effect_measure %in% c("OR", "RR", "IRR")
      # Resampling-approximation spec for consistency_method = "resample".
      # consistency_threshold is already on the comparison scale (log for
      # ratio measures, identity for identity measures), so it is passed as
      # comparison_threshold to avoid a second log-transform. adjust_covariates
      # is carried here because the GLM adjustment lives inside estimator_fn and
      # is otherwise NULL on the consistency path.
      consistency_overrides$glm_resample_spec <- list(
        outcome_type         = outcome_type,
        effect_measure       = effect_measure,
        treat.name           = treat.name,
        outcome.name         = outcome.name,
        offset.name          = offset.name,
        adjust_covariates    = adjust_covariates,
        adverse_outcome      = adverse_outcome,
        comparison_threshold = consistency_threshold
      )
    }

    # Run subgroup consistency analysis with error handling
    sc_filtered_args <- filter_call_args(
      args_call_all,
      subgroup.consistency,
      consistency_overrides
    )

    # Diagnostic: detect zero-length arguments that would cause
    # "argument is of length zero" errors in if() statements
    # Exclude twostage_args and parallel_args: empty list() is valid (use defaults)
    sc_zero_len <- names(sc_filtered_args)[
      vapply(sc_filtered_args, function(x) length(x) == 0 && !is.null(x), logical(1))
    ]
    sc_zero_len <- setdiff(sc_zero_len, c("twostage_args", "parallel_args",
                                          "estimator_fn", "consistency_threshold"))
    if (length(sc_zero_len) > 0) {
      warning("Zero-length arguments passed to subgroup.consistency: ",
              paste(sc_zero_len, collapse = ", "))
      if (details) {
        cat("*** DIAGNOSTIC: Zero-length args for subgroup.consistency:\n")
        for (zarg in sc_zero_len) {
          cat("    ", zarg, " = ", deparse(sc_filtered_args[[zarg]]), "\n")
        }
      }
    }

    grp.consistency <- tryCatch({
      do.call(subgroup.consistency, sc_filtered_args)
    }, error = function(e) {
      warning("Error in subgroup.consistency: ", e$message)
      if (details) {
        cat("*** subgroup.consistency error traceback:\n")
        cat("    Message: ", e$message, "\n")
        cat("    Arg names: ", paste(names(sc_filtered_args), collapse = ", "), "\n")
        # Report lengths of key arguments
        sc_key_args <- c("df", "hr.subgroups", "Lsg", "confs_labels",
                         "sg_focus", "stop_Kgroups", "parallel_args")
        for (ka in intersect(sc_key_args, names(sc_filtered_args))) {
          val <- sc_filtered_args[[ka]]
          if (is.data.frame(val)) {
            cat("    ", ka, ": data.frame [", nrow(val), " x ", ncol(val), "]\n")
          } else {
            cat("    ", ka, ": length=", length(val), " class=", class(val)[1], "\n")
          }
        }
      }
      return(NULL)
    })

    # Handle errors gracefully
    if (is.null(grp.consistency) || inherits(grp.consistency, "try-error")) {
      if (details) {
        cat("Consistency analysis failed - proceeding without results\n")
      }
      grp.consistency <- NULL
    }

    # Update timing
    t.end_all <- proc.time()[3]
    t.min_all <- (t.end_all - t.start_all) / 60

    if (details) {
      cat("Seconds and minutes forestsearch overall =", round(c(t.min_all*60,t.min_all), 4), "\n")
      if (!is.null(grp.consistency$algorithm)) {
        cat("Consistency algorithm used:", grp.consistency$algorithm, "\n")
      }
    }

    # -------------------------------------------------------------------------
    # Process results if consistency analysis succeeded
    # -------------------------------------------------------------------------

    if (!is.null(grp.consistency) && !is.null(grp.consistency$sg.harm)) {
      sg.harm <- grp.consistency$sg.harm

      if (details) {
        cat("Subgroup identified:", paste(sg.harm, collapse = " & "), "\n")
      }

      # Extract prediction datasets
      temp <- grp.consistency$df_flag

      # Merge to analysis data
      # Guard: temp$id is the join key (built by extract_subgroup() with the
      # literal column name "id" by convention).  If id.name != "id" and df
      # also carries an unrelated column named "id", the merge succeeds but
      # df$id survives in the output and can confuse downstream consumers.
      if (id.name != "id" && "id" %in% names(df)) {
        stop("forestsearch(): `df.analysis` contains a column named 'id' ",
             "that is distinct from `id.name = \"", id.name, "\"`.  ",
             "Rename or drop it before calling forestsearch() to avoid an ",
             "ambiguous merge with internal subgroup-flag tables.",
             call. = FALSE)
      }
      df.est_out <- merge(df, temp, by.x = id.name, by.y = "id", all.x = TRUE)

      # Return df.predict
      if (!is.null(df.predict)) {
        if (id.name != "id" && "id" %in% names(df.predict)) {
          stop("forestsearch(): `df.predict` contains a column named 'id' ",
               "that is distinct from `id.name = \"", id.name, "\"`.  ",
               "Rename or drop it before calling forestsearch() to avoid an ",
               "ambiguous merge with internal subgroup-flag tables.",
               call. = FALSE)
        }
        df.predict_out <- merge(df.predict, temp, by.x = id.name, by.y = "id", all.x = TRUE)
      }

      # Return df.test
      if (!is.null(df.test)) {
        df.test_out <- get_dfpred(
          df.predict = df.test,
          sg.harm = grp.consistency$sg.harm,
          version = 2
        )
      }
    }
  } # End has_subgroups

  # ===========================================================================
  # SECTION 9B: TIER-2 DE-BIASED GATE (optional)
  # ===========================================================================
  # Fast multiplier-bootstrap approximation of the bootstrap bias-corrected
  # effect for the selected subgroup, with a near-null "consistent with harm"
  # flag.  The selection-bias term must mirror the FULL search optimism, so the
  # candidate family here is the SAME full enumerated space the search optimized
  # over (all <= maxk cut combinations of Z), re-applying the screen on each
  # multiplier draw -- NOT just the post-screening survivors in hr.subgroups
  # (which would capture only selection-stage, not screening-stage, optimism).
  # Per-candidate dfbeta is re-derived via the consistency engine's own pieces
  # (~seconds; a future ledger from the search makes this zero-refit).
  #
  # Two paths are supported: the GLM resample-consistency path (estimator_fn
  # non-NULL, consistency_method = "resample"), and the survival/Cox path
  # (estimator_fn NULL, outcome_type = "survival").  On the Cox path the gate is
  # a fast standalone approximation -- it re-derives Cox dfbeta and re-selects
  # over the family itself, so it does not require resample consistency and runs
  # under the default split-consistency search.  Default off.
  debias_gate_out <- NULL
  .dg_glm_ok <- consistency_method == "resample" && !is.null(estimator_fn)
  .dg_cox_ok <- outcome_type == "survival" && is.null(estimator_fn)
  if (isTRUE(debias_gate) && !is.null(sg.harm) &&
      (.dg_glm_ok || .dg_cox_ok) &&
      !is.null(grp.consistency) && !is.null(grp.consistency$sg.harm.id)) {

    .g_dg <- function(a, b) if (is.null(a)) b else a

    debias_gate_out <- tryCatch({
      # Full candidate space: enumerate all <= maxk combinations of the cut
      # matrix Z with the search's OWN helpers, so the family is identical to
      # the space subgroup.search() ranked over (faithful to the Tier-1 search).
      L_dg     <- ncol(Z)
      combo_dg <- generate_combination_indices(L_dg, maxk)
      tot_dg   <- calculate_max_combinations(L_dg, maxk)
      fam <- list()
      for (kk in seq_len(tot_dg)) {
        covs.in <- get_covs_in(
          kk, maxk, L_dg,
          combo_dg$counts_1, combo_dg$indices_1,
          combo_dg$counts_2, combo_dg$indices_2,
          combo_dg$counts_3, combo_dg$indices_3)
        k_sel <- sum(covs.in)
        if (k_sel < 1L || k_sel > maxk) next
        mem <- which(get_subgroup_membership(Z, covs.in))
        if (length(mem) >= n.min)
          fam[[paste(colnames(Z)[covs.in == 1], collapse = " & ")]] <- mem
      }

      # Spec + comparison-scale thresholds differ by path.  Survival uses the
      # internal df.fs column names (Y/Event/Treat) and log-HR thresholds; GLM
      # uses the user-facing names and the already-log effect/consistency
      # thresholds.  Cox dispatch happens inside the engine via outcome_type.
      if (.dg_cox_ok) {
        adj_dg <- intersect(adjust_covariates, names(df.fs))  # only if present
        gspec <- list(outcome_type = "survival", effect_measure = "HR",
                      treat.name = "Treat", outcome.name = "Y",
                      event.name = "Event", offset.name = NULL,
                      adjust_covariates = if (length(adj_dg)) adj_dg else NULL,
                      adverse_outcome = TRUE)
        c_screen_dg      <- log(hr.threshold)
        c_consistency_dg <- log(max(hr.consistency, 1e-3))
      } else {
        gspec <- list(outcome_type = outcome_type, effect_measure = effect_measure,
                      treat.name = treat.name, outcome.name = outcome.name,
                      event.name = event.name, offset.name = offset.name,
                      adjust_covariates = adjust_covariates,
                      adverse_outcome = adverse_outcome)
        c_screen_dg      <- effect_threshold      # comparison scale (log for ratio)
        c_consistency_dg <- consistency_threshold # comparison scale
      }

      fs_debias_gate(
        df = df.fs, candidates = fam, spec = gspec,
        selected_members = which(grp.consistency$sg.harm.id == 1),
        c_screen      = c_screen_dg,
        c_consistency = c_consistency_dg,
        p_star        = pconsistency.threshold,
        t_gate        = debias_gate_args$t_gate,            # NULL -> near-null default
        gate          = .g_dg(debias_gate_args$gate,        "point"),
        reselection   = .g_dg(debias_gate_args$reselection,
                              .fs_dg_reselection_from_focus(sg_focus,
                                                            engine = "consistency")),
        effect_neighborhood = effect_neighborhood,
        selection_rule = .g_dg(debias_gate_args$selection_rule, selection_rule),
        draws         = .g_dg(debias_gate_args$draws,       2000L),
        multiplier    = .g_dg(debias_gate_args$multiplier,  "poisson"),
        include_complement = .g_dg(debias_gate_args$include_complement, TRUE),
        ci_method     = .g_dg(debias_gate_args$ci_method,   "ij"),
        seed          = .g_dg(debias_gate_args$seed,        seedit))
    }, error = function(e) {
      warning("debias_gate failed: ", conditionMessage(e)); NULL
    })

    if (!quiet && !is.null(debias_gate_out) &&
        !is.na(debias_gate_out$harm_flag)) {
      cat(sprintf("De-biased gate: %s = %.3f (gate %s %.2f) -> %s\n",
                  .g_dg(debias_gate_out$measure, "effect"),
                  debias_gate_out$debiased$est, debias_gate_out$gate$type,
                  debias_gate_out$gate$t_gate,
                  if (debias_gate_out$harm_flag) "consistent with harm"
                  else "not flagged"))
    }
  }

  # ===========================================================================
  # SECTION 10: COMPILE AND RETURN OUTPUT
  # ===========================================================================

  out <- list(
    grp.consistency = grp.consistency,
    find.grps = find.grps,
    confounders.candidate = FSconfounders.name,
    confounders.evaluated = confs_labels,
    df.est = df.est_out,
    df.predict = df.predict_out,
    df.test = df.test_out,
    minutes_all = t.min_all,
    grf_res = grf_res,
    sg_focus = sg_focus,
    sg.harm = sg.harm,
    # Tier-2 de-biased gate (NULL unless debias_gate = TRUE)
    debias_gate = debias_gate_out,
    harm_flag_debiased = if (!is.null(debias_gate_out)) debias_gate_out$harm_flag
                         else NA,
    grf_cuts = grf_cuts,
    dina_res = dina_res,
    dina_cuts = dina_cuts,
    prop_maxk = prop_maxk,
    max_sg_est = max_sg_est,
    grf_plot = grf_plot,
    args_call_all = args_call_all,
    # NEW: Include algorithm information
    consistency_algorithm = if (!is.null(grp.consistency)) grp.consistency$algorithm else NA,
    # NEW: GLM outcome information
    outcome_type = outcome_type,
    effect_measure = effect_measure,
    # NEW: Resolved threshold configuration
    threshold_config = if (exists("threshold_config")) threshold_config
                       else NULL
  )

  class(out) <- c("forestsearch", "list")
  return(out)
}


