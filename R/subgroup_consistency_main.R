# =============================================================================
# Subgroup Consistency Evaluation with Early Stopping
# =============================================================================
#
# Main function for evaluating subgroup consistency using split-sample
# validation. Supports both fixed-sample and two-stage adaptive algorithms.
#
# Key features:
#   - Split-sample consistency evaluation
#   - Two-stage adaptive algorithm option
#   - Parallel processing support
#   - Early stopping when stop_threshold is met
#   - Configurable batch size for parallel execution
#
# =============================================================================

#' Evaluate Subgroup Consistency
#'
#' Evaluates candidate subgroups using split-sample consistency validation.
#' For each candidate, repeatedly splits the data and checks whether the
#' treatment effect direction is consistent across splits.
#'
#' @param df Data frame containing the analysis dataset. Must include columns
#'   for outcome (Y), event indicator (Event), and treatment (Treat).
#' @param hr.subgroups Data.table of candidate subgroups from subgroup search,
#'   containing columns: HR, n, E, K, d0, d1, m0, m1, grp, and factor indicators.
#' @param hr.threshold Numeric. Minimum hazard ratio threshold for candidates.
#'   Default: 1.0
#' @param hr.consistency Numeric. Minimum HR required in each split for
#'   consistency. Default: 1.0
#' @param pconsistency.threshold Numeric. Minimum proportion of splits that
#'   must be consistent. Default: 0.9
#' @param m1.threshold Numeric. Maximum m1 threshold for filtering. Default: Inf
#' @param n.splits Integer. Number of splits for consistency evaluation.
#'   Default: 100
#' @param details Logical. Print progress details. Default: FALSE
#' @param by.risk Numeric. Risk interval for KM plots. Default: 12
#' @param plot.sg Logical. Generate subgroup plots. Default: FALSE
#' @param maxk Integer. Maximum number of factors in subgroup. Default: 7
#' @param Lsg List of subgroup parameters.
#' @param confs_labels Character vector mapping factor names to labels.
#' @param sg_focus Character. Subgroup selection criterion. One of:
#'   \describe{
#'     \item{\code{"hr"} (or \code{"eff"})}{Sort by
#'       \eqn{(-Pcons, -hr, K)}; pick top.}
#'     \item{\code{"maxSG"}}{Sort by \eqn{(-N, -Pcons, K)}; pick top.}
#'     \item{\code{"minSG"}}{Sort by \eqn{(N, -Pcons, K)}; pick top.}
#'     \item{\code{"hrMaxSG"} (or \code{"effMaxSG"})}{Among candidates
#'       with effect size within \code{effect_neighborhood} of the
#'       maximum, pick the one with the largest \eqn{N}.}
#'     \item{\code{"hrMinSG"} (or \code{"effMinSG"})}{Among candidates
#'       with effect size within \code{effect_neighborhood} of the
#'       maximum, pick the one with the smallest \eqn{N}.}
#'   }
#'   Default: \code{"hr"}.  The \code{"eff*"} forms are aliases for
#'   the \code{"hr*"} forms and read more naturally in GLM contexts
#'   (continuous MD, binary OR/RR/RD, count IRR).  Both vocabularies
#'   produce identical results.
#' @param selection_rule Character. Rule defining the candidate
#'   inclusion set for \code{"hrMaxSG"} / \code{"hrMinSG"}.  One of
#'   \code{"neighborhood"} (default; current behaviour),
#'   \code{"pareto"} (2D Pareto-non-dominated set in (effect, N)),
#'   or \code{"both"} (intersection of the two).  Forwarded to
#'   \code{\link{sort_subgroups_preview}} (pre-consistency triage)
#'   and \code{\link{sort_subgroups}} (post-consistency selection).
#'   Must be \code{"neighborhood"} for other \code{sg_focus} values.
#' @param effect_neighborhood Numeric in \code{[0, 1)}.  Relative
#'   tolerance for \code{"hrMaxSG"} and \code{"hrMinSG"}.  A candidate
#'   is in the \emph{effect-size neighborhood} iff its (natural-scale)
#'   effect is at least \code{(1 - effect_neighborhood) * max(effect)}.
#'   Default \code{0.10} (i.e., within 10\% of the strongest effect).
#'   For GLM ratio measures (OR, IRR), the neighborhood test is applied
#'   on the natural scale (after exponentiation of \code{hr}).  Ignored
#'   when \code{sg_focus} is \code{"hr"}, \code{"maxSG"}, or
#'   \code{"minSG"}.  When \code{selection_rule = "pareto"}, must be
#'   left at its default (the rule does not consult it).
#' @param stop_Kgroups Integer. Maximum number of candidates to evaluate.
#'   Default: 10
#' @param stop_threshold Numeric in \code{[0, 1]} or \code{NULL}.
#'   When a candidate subgroup's consistency probability (\code{Pcons})
#'   meets or exceeds this threshold, evaluation stops early -- remaining
#'   candidates are skipped. Set to \code{NULL} to disable early stopping
#'   and evaluate all candidates up to \code{stop_Kgroups}. Default: \code{NULL}.
#'
#'   \strong{Note:} Values > 1.0 are not permitted. To disable early
#'   stopping, use \code{stop_threshold = NULL}, not a value above 1.
#'
#'   \strong{Interaction with \code{sg_focus}:}
#'   \describe{
#'     \item{\code{"hr"}, \code{"maxSG"}, \code{"minSG"}}{Early stopping
#'       is valid because candidates are sorted by a single criterion.
#'       The first candidate passing the threshold is optimal under that
#'       criterion.}
#'     \item{\code{"hrMaxSG"}, \code{"hrMinSG"}}{Should generally be
#'       \code{NULL}, because these compound criteria require comparing
#'       HR \emph{and} size across all candidates. \code{forestsearch()}
#'       automatically resets to \code{NULL} with a warning for these.}
#'   }
#'
#'   For parallel execution, early stopping is checked after each batch
#'   completes, so some additional candidates beyond the first meeting the
#'   threshold may be evaluated. Use a smaller \code{batch_size} in
#'   \code{parallel_args} for finer-grained early stopping.
#' @param show_candidate_summary Logical. If \code{TRUE}, prints two
#'   diagnostic tables: a \strong{pre-consistency} candidate preview
#'   (all top-\code{stop_Kgroups} candidates entering consistency
#'   evaluation, with Frontier and InBand flags) and a
#'   \strong{post-consistency} summary (passing candidates only, with
#'   Pcons, Frontier, InBand, and Selected flags).  The two views
#'   together make the rule's filter visible end-to-end.  Pcons is
#'   shown only for passing candidates (the consistency evaluator
#'   returns NULL for non-passing ones).  Default: FALSE.
#' @param pconsistency.digits Integer. Decimal places for consistency
#'   proportion. Default: 2
#' @param seed Integer. Random seed for reproducible consistency splits.
#'   Default: 8316951. Set to NULL for non-reproducible random splits.
#'   The seed is used both for sequential execution (via set.seed()) and
#'   parallel execution (via future.seed).
#' @param checking Logical. Enable additional validation checks. Default: FALSE
#' @param estimator_fn Closure or \code{NULL}. Effect-estimator closure from
#'   \code{\link{make_effect_estimator}} for GLM outcomes.  Default \code{NULL}
#'   (uses Cox model).
#' @param consistency_threshold Numeric or \code{NULL}. Threshold for GLM
#'   consistency evaluation.  Default \code{NULL} (uses \code{hr.consistency}).
#' @param adjust_covariates Character vector or \code{NULL}. Additional Cox
#'   model terms used when scoring survival subgroups during consistency
#'   evaluation, e.g. \code{"strata(x1)"} for a stratified model or
#'   \code{"x2"} for linear adjustment.  Referenced columns must be present
#'   in \code{df}.  Ignored on the GLM path.  Default \code{NULL}
#'   (treatment-only, unadjusted).
#' @param consistency_method Character. \code{"split"} (default) for the literal
#'   repeated 50/50 split-and-refit consistency calculation, or
#'   \code{"resample"} for the multiplier (influence-function) approximation via
#'   \code{\link{consistency_resample}}.  On the GLM path \code{"resample"}
#'   requires a \code{glm_resample_spec}; without one (or when
#'   \code{consistency_resample} is unavailable) it falls back to
#'   \code{"split"}.
#' @param glm_resample_spec List or \code{NULL}. GLM resampling specification
#'   threaded from \code{\link{forestsearch}} when
#'   \code{consistency_method = "resample"} on a GLM outcome: components
#'   \code{outcome_type}, \code{effect_measure}, \code{treat.name},
#'   \code{outcome.name}, \code{offset.name}, \code{adjust_covariates},
#'   \code{adverse_outcome}, and \code{comparison_threshold} (the consistency
#'   threshold already on the comparison scale, i.e. \code{log(threshold)} for
#'   ratio measures).  \code{NULL} on the survival path.
#' @param effect_label Character. Column label for the effect measure in
#'   diagnostic output (e.g., candidate subgroup table).  Default \code{"HR"}.
#'   Set automatically by \code{\link{forestsearch}} to the resolved
#'   \code{effect_measure} (e.g., \code{"IRR"}, \code{"OR"}) for GLM outcomes.
#' @param effect_log_scale Logical. If \code{TRUE}, the \code{HR} column in
#'   \code{hr.subgroups} stores log-scale values (log-OR, log-IRR) and will
#'   be exponentiated for display in diagnostic output.  Default \code{FALSE}
#'   (Cox stores natural-scale HR).
#' @param use_twostage Logical. Use two-stage adaptive algorithm. Default: FALSE
#' @param twostage_args List. Parameters for two-stage algorithm:
#'   \describe{
#'     \item{n.splits.screen}{Splits for Stage 1 screening. Default: 30}
#'     \item{screen.threshold}{Consistency threshold for Stage 1. Default: auto}
#'     \item{batch.size}{Splits per batch in Stage 2. Default: 20}
#'     \item{conf.level}{Confidence level for early stopping. Default: 0.95}
#'     \item{min.valid.screen}{Minimum valid Stage 1 splits. Default: 10}
#'   }
#' @param parallel_args List. Parallel processing configuration:
#'   \describe{
#'     \item{plan}{Future plan: "multisession", "multicore", or "sequential"}
#'     \item{workers}{Number of parallel workers}
#'     \item{batch_size}{Number of candidates to evaluate per batch. Smaller
#'       values provide finer-grained early stopping but may increase overhead.
#'       Default: When stop_threshold is set and sg_focus is "hr" or "minSG",
#'       defaults to 1 (stop immediately when first candidate passes). For other
#'       sg_focus values with stop_threshold, defaults to min(workers, n_candidates/4).
#'       When stop_threshold is NULL, defaults to workers*2 for efficiency.}
#'     \item{show_message}{Print parallel config messages}
#'   }
#'
#' @return A list containing:
#'   \describe{
#'     \item{out_sg}{Selected subgroup results.  When non-\code{NULL},
#'       contains \code{result} (sorted candidate table; top row is the
#'       selected subgroup) and \code{pareto_frontier} (data.table of
#'       non-dominated candidates on (effect, N), both maximized -- a
#'       post-hoc diagnostic, not used for selection).}
#'     \item{sg_focus}{Selection criterion used}
#'     \item{df_flag}{Data frame with treatment recommendations}
#'     \item{sg.harm}{\strong{Subgroup definition labels} -- character
#'       vector of factor-level cut names (e.g., \code{c("z1.1", "z2.1")}),
#'       of length equal to the number of cuts defining the subgroup.
#'       \code{NULL} if no subgroup was identified.}
#'     \item{sg.harm.id}{\strong{Per-subject subgroup-membership
#'       indicator} -- integer vector of length \code{nrow(df)} with
#'       \code{1} if subject \eqn{i} is in the identified subgroup and
#'       \code{0} otherwise.  \code{NULL} if no subgroup was identified.
#'       \strong{Not} a character vector of cut expressions; see the
#'       \emph{Field naming collision} section below.}
#'     \item{algorithm}{"twostage" or "fixed"}
#'     \item{n_candidates_evaluated}{Number of candidates actually evaluated}
#'     \item{n_candidates_total}{Total candidates available}
#'     \item{n_passed}{Number meeting consistency threshold}
#'     \item{early_stop_triggered}{Logical indicating if early stop occurred}
#'     \item{early_stop_candidate}{Index of candidate triggering early stop}
#'     \item{stop_threshold}{Threshold used for early stopping}
#'     \item{seed}{Random seed used for reproducibility (NULL if not set)}
#'   }
#'
#' @section Field naming collision with GRF results:
#' The field name \code{sg.harm.id} has \strong{different semantics} on
#' this object versus on the result objects returned by
#' \code{\link{grf.subg.harm.glm}} and \code{\link{grf.subg.harm.survival}}:
#'
#' \tabular{lll}{
#'   \strong{Object} \tab \strong{\code{sg.harm.id} contains} \tab \strong{Length / type} \cr
#'   \code{subgroup.consistency()} result (this function) \tab per-subject 0/1 membership indicator \tab integer, length \code{nrow(df)} \cr
#'   \code{grf.subg.harm.glm()} result \tab character vector of cut expressions \tab character, length = depth of selected tree \cr
#'   \code{grf.subg.harm.survival()} result \tab character vector of cut expressions \tab character, length = depth of selected tree \cr
#' }
#'
#' \strong{Practical consequence.}  Do not paste \code{sg.harm.id} with
#' \code{paste(..., collapse = " & ")} to print "the identified subgroup"
#' without first checking object class.  For subgroup labels from a
#' \code{\link{forestsearch}} or \code{subgroup.consistency()} result,
#' use \code{sg.harm} (character vector of cut names) -- the FS main
#' result exposes \code{sg.harm} at the top level.
#'
#' This naming collision is a documented CRAN-stable API for v0.1.x
#' and v0.2.x; it is expected to be resolved via a deprecation cycle
#' in a future minor release.
#'
#' @examples
#' \dontrun{
#' # Standard evaluation
#' result <- subgroup.consistency(
#'   df = trial_data,
#'   hr.subgroups = candidates,
#'   sg_focus = "hr",
#'   n.splits = 400,
#'   parallel_args = list(plan = "multisession", workers = 6)
#' )
#'
#' # Show top 10 candidates before evaluation
#' result <- subgroup.consistency(
#'   df = trial_data,
#'   hr.subgroups = candidates,
#'   sg_focus = "hr",
#'   show_candidate_summary = TRUE,  # Post-consistency summary
#'   n.splits = 400
#' )
#'
#' # With early stopping and custom batch size
#' result <- subgroup.consistency(
#'   df = trial_data,
#'   hr.subgroups = candidates,
#'   sg_focus = "hr",
#'   stop_threshold = 0.95,
#'   show_candidate_summary = TRUE,
#'   parallel_args = list(
#'     plan = "multisession",
#'     workers = 6,
#'     batch_size = 2  # Check early stopping after every 2 candidates
#'   )
#' )
#' }
#'
#' @importFrom data.table copy as.data.table is.data.table
#' @importFrom survival coxph Surv
#' @importFrom future.apply future_lapply
#' @importFrom utils modifyList
#' @export

subgroup.consistency <- function(df,
                                 hr.subgroups,
                                 hr.threshold = 1.0,
                                 hr.consistency = 1.0,
                                 pconsistency.threshold = 0.9,
                                 m1.threshold = Inf,
                                 n.splits = 100,
                                 details = FALSE,
                                 by.risk = 12,
                                 plot.sg = FALSE,
                                 maxk = 7,
                                 Lsg,
                                 confs_labels,
                                 sg_focus = "hr",
                                 selection_rule = "neighborhood",
                                 effect_neighborhood = 0.10,
                                 stop_Kgroups = 10,
                                 stop_threshold = NULL,
                                 show_candidate_summary = FALSE,
                                 pconsistency.digits = 2,
                                 seed = 8316951,
                                 checking = FALSE,
                                 use_twostage = FALSE,
                                 twostage_args = list(),
                                 parallel_args = list(),
                                 estimator_fn = NULL,
                                 consistency_threshold = NULL,
                                 effect_label = "HR",
                                 effect_log_scale = FALSE,
                                 adjust_covariates = NULL,
                                 consistency_method = c("split", "resample"),
                                 glm_resample_spec = NULL) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION
  # ===========================================================================

  # Accept the GLM-natural vocabulary (effMaxSG / effMinSG / eff) alongside
  # the canonical Cox-flavored names (hrMaxSG / hrMinSG / hr).  Normalize
  # at the entry point so all downstream code -- including the early-stop
  # batching at Section 8, sort_subgroups_preview() / sort_subgroups() in
  # the helpers, and the result-building in Section 10 -- sees the
  # canonical form.
  sg_focus <- .normalize_sg_focus(sg_focus)

  # Consistency calculation method: literal splitting vs. resampling
  # approximation. 'resample' covers the Cox (survival) path and the GLM path
  # (when forestsearch supplies a glm_resample_spec for a supported measure);
  # a missing consistency_resample() or an unspecified GLM path fall back to
  # 'split'.
  consistency_method <- match.arg(consistency_method)
  if (identical(consistency_method, "resample")) {
    if (!exists("consistency_resample", mode = "function")) {
      warning("consistency_resample() not found; using 'split'. Install or ",
              "source consistency_resample.R to enable the approximation.")
      consistency_method <- "split"
    } else if (!is.null(estimator_fn) && is.null(glm_resample_spec)) {
      warning("consistency_method = 'resample' on the GLM (estimator_fn) path ",
              "requires a glm_resample_spec; using 'split'.")
      consistency_method <- "split"
    }
  }

  if (!is.data.frame(df) || nrow(df) == 0) {
    stop("df must be a non-empty data frame")
  }

  # For survival (default), require Y/Event/Treat; for GLM, the estimator

  # closure handles column resolution internally
  if (is.null(estimator_fn)) {
    required_cols <- c("Y", "Event", "Treat")
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop("df missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }

  if (!data.table::is.data.table(hr.subgroups)) {
    hr.subgroups <- data.table::as.data.table(hr.subgroups)
  }

  if (nrow(hr.subgroups) == 0) {
    warning("No candidate subgroups provided")
    return(list(
      out_sg = NULL, sg_focus = sg_focus,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = 0, n_candidates_total = 0, n_passed = 0,
      early_stop_triggered = FALSE, early_stop_candidate = NA_integer_,
      stop_threshold = stop_threshold,
      seed = seed
    ))
  }

  # Validate stop_threshold
  if (!is.null(stop_threshold)) {
    if (!is.numeric(stop_threshold) || length(stop_threshold) != 1 ||
        stop_threshold < 0 || stop_threshold > 1) {
      stop("stop_threshold must be NULL or a numeric value between 0 and 1")
    }
    if (stop_threshold < pconsistency.threshold) {
      warning("stop_threshold (", stop_threshold, ") is less than ",
              "pconsistency.threshold (", pconsistency.threshold, "). ",
              "Early stopping may never trigger for passing candidates.")
    }
  }

  # Validate batch_size if provided
  if (!is.null(parallel_args$batch_size)) {
    if (!is.numeric(parallel_args$batch_size) ||
        parallel_args$batch_size < 1) {
      stop("parallel_args$batch_size must be a positive integer")
    }
  }

  # ===========================================================================
  # SECTION 1b: SET RANDOM SEED FOR REPRODUCIBILITY
  # ===========================================================================

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ===========================================================================
  # SECTION 2: EXTRACT FACTOR NAMES
  # ===========================================================================

  names.Z <- names(hr.subgroups)[
    !names(hr.subgroups) %in% c("K", "n", "E", "d0", "d1", "m0", "m1",
                                 "HR", "L(HR)", "U(HR)", "grp")
  ]

  if (length(names.Z) == 0) {
    stop("No factor columns found in hr.subgroups")
  }

  # ===========================================================================
  # SECTION 3: SETUP TWO-STAGE PARAMETERS
  # ===========================================================================

  ts_defaults <- list(
    n.splits.screen = 30,
    screen.threshold = NULL,
    batch.size = 20,
    conf.level = 0.95,
    min.valid.screen = 10
  )
  ts_params <- modifyList(ts_defaults, twostage_args)

  if (is.null(ts_params$screen.threshold)) {
    se_screen <- sqrt(pconsistency.threshold * (1 - pconsistency.threshold) /
                        ts_params$n.splits.screen)
    ts_params$screen.threshold <- max(0.5, pconsistency.threshold - 2.5 * se_screen)
  }

  # ===========================================================================
  # SECTION 4: FILTER CANDIDATES BY HR THRESHOLD
  # ===========================================================================

  if (nrow(hr.subgroups) == 0 || !"HR" %in% names(hr.subgroups)) {
    warning("No valid hr.subgroups")
    return(list(
      out_sg = NULL, sg_focus = sg_focus,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = 0, n_candidates_total = 0, n_passed = 0,
      early_stop_triggered = FALSE, early_stop_candidate = NA_integer_,
      stop_threshold = stop_threshold,
      seed = seed
    ))
  }

  if (is.finite(m1.threshold)) {
    hr.subgroups <- hr.subgroups[!is.na(hr.subgroups$m1), ]
    if (nrow(hr.subgroups) == 0) {
      warning("All subgroups removed after filtering NA m1 values")
      return(list(
        out_sg = NULL, sg_focus = sg_focus,
        df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
        algorithm = ifelse(use_twostage, "twostage", "fixed"),
        n_candidates_evaluated = 0, n_candidates_total = 0, n_passed = 0,
        early_stop_triggered = FALSE, early_stop_candidate = NA_integer_,
        stop_threshold = stop_threshold,
      seed = seed
      ))
    }
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold &
                                hr.subgroups$m1 <= m1.threshold, ]
  } else {
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]
  }

  if (nrow(found.hrs) == 0) {
    if (details) {
      cat("No subgroups meet criteria (", effect_label, " >=",
          if (effect_log_scale) exp(hr.threshold) else hr.threshold)
      if (is.finite(m1.threshold)) cat(" and m1 <=", m1.threshold)
      cat(")\n")
    }
    return(list(
      out_sg = NULL, sg_focus = sg_focus,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = 0, n_candidates_total = 0, n_passed = 0,
      early_stop_triggered = FALSE, early_stop_candidate = NA_integer_,
      stop_threshold = stop_threshold,
      seed = seed
    ))
  }

  # ===========================================================================
  # SECTION 5: REMOVE DUPLICATES AND SORT BY SG_FOCUS
  # ===========================================================================

  if (nrow(found.hrs) > 1) {
    n_before <- nrow(found.hrs)

    tryCatch({
      found.hrs <- remove_near_duplicate_subgroups(found.hrs, details = details)
    }, error = function(e) {
      warning("Error removing duplicates: ", e$message,
              ". Proceeding with original subgroups.")
    })

    if (nrow(found.hrs) == 0) {
      stop("All subgroups removed during duplicate removal.")
    }
  }


  found.hrs <- sort_subgroups_preview(found.hrs, sg_focus,
                                      selection_rule      = selection_rule,
                                      effect_neighborhood = effect_neighborhood,
                                      effect_log_scale    = effect_log_scale)


  # Extract index matrix
  index.Z <- found.hrs[, names.Z, with = FALSE]

  if (details) {
    cat("# of unique initial candidates:", nrow(found.hrs), "\n")
  }

  # Limit to top stop_Kgroups candidates
  maxsgs <- min(nrow(found.hrs), stop_Kgroups)
  found.hrs <- found.hrs[seq_len(maxsgs), ]
  index.Z <- index.Z[seq_len(maxsgs), ]

  n_candidates <- nrow(found.hrs)

  if (details) {
    cat("# Restricting to top stop_Kgroups =", stop_Kgroups, "\n")
    cat("# of candidates to evaluate:", n_candidates, "\n")
    if (!is.null(stop_threshold)) {
      cat("# Early stop threshold:", stop_threshold, "\n")
    }
  }

  # ===========================================================================
  # SECTION 5b: PRE-CONSISTENCY CANDIDATE PREVIEW
  # ===========================================================================
  # Shows the candidates that will be evaluated for consistency, in their
  # preview-sort order.  Companion to Section 10b (post-consistency
  # summary).  Both views are gated on the same flag.
  if (isTRUE(show_candidate_summary) && details) {
    print_candidate_preview(
      found.hrs           = found.hrs,
      index.Z             = index.Z,
      names.Z             = names.Z,
      confs_labels        = confs_labels,
      sg_focus            = sg_focus,
      selection_rule      = selection_rule,
      effect_neighborhood = effect_neighborhood,
      effect_log_scale    = effect_log_scale,
      effect_label        = effect_label
    )
  }

  # ===========================================================================
  # SECTION 6: VALIDATE PARALLEL CONFIGURATION
  # ===========================================================================

  use_parallel <- length(parallel_args) > 0 && !is.null(parallel_args[[1]])

  if (use_parallel) {
    required_parallel <- c("plan", "workers")
    if (!all(required_parallel %in% names(parallel_args))) {
      warning("parallel_args missing required elements. Using sequential.")
      use_parallel <- FALSE
    }

    valid_plans <- c("multisession", "multicore", "callr", "sequential")
    if (use_parallel && !parallel_args$plan %in% valid_plans) {
      warning("Invalid parallel plan. Using sequential.")
      use_parallel <- FALSE
    }

    if (use_parallel && (!is.numeric(parallel_args$workers) ||
                         parallel_args$workers < 1)) {
      warning("Invalid workers value. Using sequential.")
      use_parallel <- FALSE
    }
  }

  # ===========================================================================
  # SECTION 7: INITIALIZE TRACKING VARIABLES
  # ===========================================================================

  early_stop_triggered <- FALSE
  early_stop_candidate <- NA_integer_
  n_evaluated <- 0L
  results_list <- vector("list", n_candidates)

  # ===========================================================================
  # SECTION 8: EVALUATE CANDIDATES
  # ===========================================================================

  if (!use_parallel) {
    # -------------------------------------------------------------------------
    # SEQUENTIAL EXECUTION
    # -------------------------------------------------------------------------

    for (m in seq_len(n_candidates)) {

      n_evaluated <- m

      if (use_twostage) {
        results_list[[m]] <- evaluate_consistency_twostage(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = details,
          n.splits.screen = ts_params$n.splits.screen,
          screen.threshold = ts_params$screen.threshold,
          n.splits.max = n.splits,
          batch.size = ts_params$batch.size,
          conf.level = ts_params$conf.level,
          min.valid.screen = ts_params$min.valid.screen,
          estimator_fn = estimator_fn,
          consistency_threshold = consistency_threshold,
          adjust_covariates = adjust_covariates,
          consistency_method = consistency_method,
          glm_resample_spec = glm_resample_spec
        )
      } else {
        results_list[[m]] <- evaluate_subgroup_consistency(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          n.splits = n.splits,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = details,
          estimator_fn = estimator_fn,
          consistency_threshold = consistency_threshold,
          adjust_covariates = adjust_covariates,
          consistency_method = consistency_method,
          glm_resample_spec = glm_resample_spec
        )
      }

      # Check early stopping condition
      if (!is.null(stop_threshold) && !is.null(results_list[[m]])) {
        pcons_m <- as.numeric(results_list[[m]]["Pcons"])

        if (!is.na(pcons_m) && pcons_m >= stop_threshold) {
          early_stop_triggered <- TRUE
          early_stop_candidate <- m

          if (details) {
            cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
            cat("EARLY STOP TRIGGERED\n")
            cat("  Candidate:", m, "of", n_candidates, "\n")
            cat("  Pcons:", round(pcons_m, 4), ">=", stop_threshold, "\n")
            cat("  HR:", round(as.numeric(results_list[[m]]["hr"]), 3), "\n")
            cat(paste(rep("=", 50), collapse = ""), "\n\n", sep = "")
          }
          break
        }
      }
    }

    if (details) {
      cat("Evaluated", n_evaluated, "of", n_candidates, "candidates",
          if (early_stop_triggered) "(early stop)" else "(complete)", "\n")
    }

  } else {
    # -------------------------------------------------------------------------
    # PARALLEL EXECUTION WITH BATCHED EARLY STOPPING
    # -------------------------------------------------------------------------

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    # Suppress package version warnings during parallel setup
    suppressWarnings({
      setup_parallel_SGcons(parallel_args)
    })

    # Determine batch size for parallel execution
    n_workers <- parallel_args$workers

    if (!is.null(parallel_args$batch_size)) {
      # User-specified batch size
      batch_size_parallel <- min(as.integer(parallel_args$batch_size), n_candidates)
    } else if (!is.null(stop_threshold)) {
      # Early stopping enabled - batch size depends on sg_focus
      if (sg_focus %in% c("hr", "minSG")) {
        # For "hr" or "minSG", first passing candidate is optimal
        # Use batch_size = 1 to stop immediately
        batch_size_parallel <- 1L
      } else {
        # For other sg_focus values, use smaller batches for granularity
        batch_size_parallel <- min(n_workers, max(1L, n_candidates %/% 4), n_candidates)
      }
    } else {
      # No early stopping - larger batches for efficiency
      batch_size_parallel <- min(n_workers * 2L, n_candidates)
    }

    if (details) {
      cat("Parallel config: workers =", n_workers,
          ", batch_size =", batch_size_parallel, "\n")
    }

    # Create evaluation function closure
    if (use_twostage) {
      eval_fun <- function(m) {
        evaluate_consistency_twostage(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = FALSE,
          n.splits.screen = ts_params$n.splits.screen,
          screen.threshold = ts_params$screen.threshold,
          n.splits.max = n.splits,
          batch.size = ts_params$batch.size,
          conf.level = ts_params$conf.level,
          min.valid.screen = ts_params$min.valid.screen,
          estimator_fn = estimator_fn,
          consistency_threshold = consistency_threshold,
          adjust_covariates = adjust_covariates,
          consistency_method = consistency_method,
          glm_resample_spec = glm_resample_spec
        )
      }
    } else {
      eval_fun <- function(m) {
        evaluate_subgroup_consistency(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          n.splits = n.splits,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = FALSE,
          estimator_fn = estimator_fn,
          consistency_threshold = consistency_threshold,
          adjust_covariates = adjust_covariates,
          consistency_method = consistency_method,
          glm_resample_spec = glm_resample_spec
        )
      }
    }

    # Process in batches
    n_batches <- ceiling(n_candidates / batch_size_parallel)

    for (batch_num in seq_len(n_batches)) {

      start_idx <- (batch_num - 1L) * batch_size_parallel + 1L
      end_idx <- min(batch_num * batch_size_parallel, n_candidates)
      batch_indices <- seq.int(start_idx, end_idx)

      if (details) {
        cat("Batch", batch_num, "/", n_batches,
            ": candidates", start_idx, "-", end_idx, "\n")
      }

      # Parallel evaluation of batch (suppress package version warnings)
      # Use seed for reproducible parallel RNG
      batch_results <- suppressWarnings({
        future.apply::future_lapply(
          batch_indices,
          eval_fun,
          future.seed = if (!is.null(seed)) seed else TRUE
        )
      })

      # Store results
      for (i in seq_along(batch_indices)) {
        results_list[[batch_indices[i]]] <- batch_results[[i]]
      }

      n_evaluated <- end_idx

      # Check early stopping (process in order to respect HR sorting)
      if (!is.null(stop_threshold)) {
        for (i in seq_along(batch_indices)) {
          result_i <- batch_results[[i]]

          if (!is.null(result_i)) {
            pcons_i <- as.numeric(result_i["Pcons"])

            if (!is.na(pcons_i) && pcons_i >= stop_threshold) {
              early_stop_triggered <- TRUE
              early_stop_candidate <- batch_indices[i]

              if (details) {
                cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
                cat("EARLY STOP TRIGGERED (batch", batch_num, ")\n")
                cat("  Candidate:", batch_indices[i], "of", n_candidates, "\n")
                cat("  Pcons:", round(pcons_i, 4), ">=", stop_threshold, "\n")
                cat(paste(rep("=", 50), collapse = ""), "\n\n", sep = "")
              }
              break
            }
          }
        }
      }

      if (early_stop_triggered) break
    }

    if (details) {
      cat("Evaluated", n_evaluated, "of", n_candidates, "candidates",
          if (early_stop_triggered) "(early stop)" else "(complete)", "\n")
    }
  }

  # ===========================================================================
  # SECTION 9: COMBINE AND FILTER RESULTS
  # ===========================================================================

  # Filter to non-NULL results
  results_list_valid <- Filter(Negate(is.null), results_list)

  # Define any.found
  any.found <- length(results_list_valid)

  if (any.found == 0) {
    if (details) cat("No subgroups found meeting consistency threshold\n")

    # Emit the SUMMARY banner even in the zero-passed case so the
    # downstream display (e.g., comparison vignette) shows a real
    # post-consistency summary block rather than "(not available)".
    # print_candidate_summary() handles out_sg = NULL: it emits the
    # banner, the Evaluated/Passed=0 line, and the closing bar.
    if (isTRUE(show_candidate_summary) && details) {
      print_candidate_summary(
        out_sg              = NULL,
        n_evaluated         = n_evaluated,
        sg_focus            = sg_focus,
        selection_rule      = selection_rule,
        effect_neighborhood = effect_neighborhood,
        effect_log_scale    = effect_log_scale,
        effect_label        = effect_label
      )
    }

    return(list(
      out_sg = NULL,
      sg_focus = sg_focus,
      df_flag = NULL,
      sg.harm = NULL,
      sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = n_evaluated,
      n_candidates_total = n_candidates,
      n_passed = 0L,
      early_stop_triggered = early_stop_triggered,
      early_stop_candidate = early_stop_candidate,
      stop_threshold = stop_threshold,
      seed = seed
    ))
  }

  # Convert to data.table
  res <- data.table::as.data.table(do.call(rbind, results_list_valid))

  if (details) {
    cat(any.found, "subgroups passed consistency threshold\n")
  }

  # Convert columns to numeric
  cols_to_numeric <- c("Pcons", "hr", "N", "E", "K")
  missing_cols <- setdiff(cols_to_numeric, names(res))
  if (length(missing_cols) > 0) {
    stop("Result data.table missing expected columns: ",
         paste(missing_cols, collapse = ", "))
  }

  tryCatch({
    res[, (cols_to_numeric) := lapply(.SD, as.numeric), .SDcols = cols_to_numeric]
  }, error = function(e) {
    stop("Error converting result columns to numeric: ", e$message)
  })

  # ===========================================================================
  # SECTION 10: GENERATE OUTPUT
  # ===========================================================================

  out_sg <- NULL
  df_flag <- sg.harm <- sg.harm.id <- NULL

  if (any.found > 0) {
    result_new <- data.table::copy(res)

    # Validate sg_focus.  Note: hrMaxSG and hrMinSG are passed through
    # unchanged (no longer collapsed to maxSG/minSG), so that the
    # neighborhood-based selection in sort_subgroups() is applied.
    valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG")
    if (!sg_focus %in% valid_sg_focus) {
      stop(sprintf("Unknown sg_focus value: %s", sg_focus))
    }

    sgdetails <- ifelse(plot.sg, TRUE, FALSE)

    out_sg <- tryCatch({
      sg_consistency_out(
        df = df,
        result_new = result_new,
        sg_focus = sg_focus,
        details = sgdetails,
        plot.sg = sgdetails,
        index.Z = index.Z,
        names.Z = names.Z,
        by.risk = by.risk,
        confs_labels = confs_labels,
        is_glm = !is.null(estimator_fn),
        selection_rule      = selection_rule,
        effect_neighborhood = effect_neighborhood,
        effect_log_scale    = effect_log_scale
      )
    }, error = function(e) {
      warning("Error in sg_consistency_out for '", sg_focus, "': ", e$message)
      NULL
    })

    if (!is.null(out_sg)) {
      df_flag <- out_sg$df_flag
      sg.harm <- out_sg$sg.harm_label
      sg.harm.id <- out_sg$sg.harm.id
    }

    if (details) cat("SG focus =", sg_focus, "\n")
  }

  # ===========================================================================
  # SECTION 10b: PRINT POST-CONSISTENCY CANDIDATE SUMMARY
  # ===========================================================================
  # Replaces the legacy pre-consistency "TOP 10" preview (removed above).
  # Shows up to max_subgroups_search (= n_candidates) evaluated candidates,
  # but only those that passed consistency have rows in out_sg$result.
  # Frontier / band / selected flags are derived from out_sg.
  if (isTRUE(show_candidate_summary) && details) {
    print_candidate_summary(
      out_sg              = out_sg,
      n_evaluated         = n_evaluated,
      sg_focus            = sg_focus,
      selection_rule      = selection_rule,
      effect_neighborhood = effect_neighborhood,
      effect_log_scale    = effect_log_scale,
      effect_label        = effect_label
    )
  }

  # ===========================================================================
  # SECTION 11: RETURN OUTPUT
  # ===========================================================================

  return(list(
    out_sg = out_sg,
    sg_focus = sg_focus,
    df_flag = df_flag,
    sg.harm = sg.harm,
    sg.harm.id = sg.harm.id,
    algorithm = ifelse(use_twostage, "twostage", "fixed"),
    n_candidates_evaluated = n_evaluated,
    n_candidates_total = n_candidates,
    n_passed = any.found,
    early_stop_triggered = early_stop_triggered,
    early_stop_candidate = early_stop_candidate,
    stop_threshold = stop_threshold,
    seed = seed
  ))
}
