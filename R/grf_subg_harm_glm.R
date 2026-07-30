# =============================================================================
# grf.subg.harm.glm.R
# GRF-based subgroup identification for GLM outcomes
#
# Extends the forestsearch GLM extension to handle:
#   Binary outcomes   : log-OR, log-RR, RD
#   Continuous outcomes: mean difference
#   Count / rate outcomes: log-IRR (Poisson + offset, quasi-Poisson, negbin)
#
# GRF screening strategy for count outcomes (grf_count_transform):
#   "log"      [default] causal_forest() applied to log(Y + 0.5)  --
#              variance-stabilising transformation; avoids log(0) for zeros.
#   "identity" causal_forest() applied to Y directly (raw counts).
#              Appropriate when counts are large and overdispersion moderate.
#   This step is used ONLY for candidate factor selection and policy-tree
#   cut generation. All downstream effect estimation uses the correctly
#   specified GLM (Poisson / quasi-Poisson / negative binomial).
# =============================================================================

#' GRF Subgroup Identification for GLM Outcomes
#'
#' Identifies treatment effect subgroups for binary, continuous, and count
#' outcomes using Generalized Random Forests (GRF) for candidate factor
#' screening, followed by exhaustive subgroup enumeration with GLM-based
#' effect estimation.
#'
#' For count / rate outcomes, the GRF screening step applies a
#' variance-stabilising \eqn{\log(Y + 0.5)} transformation so that
#' \code{grf::causal_forest()} can be applied without dedicated
#' count-specific infrastructure. All subgroup effect estimates are
#' computed from correctly specified Poisson, quasi-Poisson, or
#' negative-binomial GLMs via \code{\link{create_glm_row}}.
#'
#' @param data Data frame containing the analysis data.
#' @param confounders.name Character vector of potential effect modifier
#'   (confounder) names. These are the candidate split variables for GRF
#'   and the subgroup enumeration.
#' @param outcome.name Character. Name of the outcome column.
#' @param treat.name Character. Name of the binary treatment indicator (0/1).
#'   Default: \code{"treat"}.
#' @param id.name Character. Name of the subject ID column.
#'   Default: \code{"id"}.
#' @param outcome_type Character. Type of outcome. One of \code{"binary"}
#'   (outcome is 0/1; use \code{effect_measure} to select log-OR, log-RR,
#'   or RD), \code{"continuous"} (estimates mean difference), or
#'   \code{"count"} (count or rate outcome; use \code{effect_measure = "log_IRR"}
#'   with \code{offset.name} for incidence rate ratios).
#'   Default: \code{"binary"}.
#' @param effect_measure Character. Effect measure to estimate. One of
#'   \code{"log_OR"} (log odds ratio, binary/logistic),
#'   \code{"log_RR"} (log relative risk, binary/modified Poisson),
#'   \code{"RD"} (risk difference, binary/identity link),
#'   \code{"MD"} (mean difference, continuous/Gaussian), or
#'   \code{"log_IRR"} (log incidence rate ratio, count/Poisson with offset).
#'   Default is \code{"log_OR"} for binary, \code{"MD"} for continuous,
#'   and \code{"log_IRR"} for count.
#' @param offset.name Character or \code{NULL}. Name of the person-time
#'   (exposure) column. When supplied, \code{log(offset)} enters the linear
#'   predictor of the GLM as an offset term. Required when
#'   \code{outcome_type = "count"} and \code{effect_measure = "log_IRR"}.
#'   Ignored for binary and continuous outcomes. Default: \code{NULL}.
#' @param overdispersion Character. Overdispersion correction for count models.
#'   One of \code{"none"} (standard Poisson, model-based SEs),
#'   \code{"quasi"} (quasi-Poisson, dispersion-corrected SEs), or
#'   \code{"negbin"} (negative binomial via \code{MASS::glm.nb()}).
#'   Ignored for binary and continuous outcomes. Default: \code{"none"}.
#' @param grf_count_transform Character. Transformation applied to the count
#'   outcome \emph{before} passing it to \code{grf::causal_forest()} for
#'   factor screening. Has no effect on binary or continuous outcomes.
#'   \code{"log"} (default) applies \eqn{Y^* = \log(Y + 0.5)}, which
#'   prevents \eqn{\log(0)} when zero counts are present and approximates
#'   variance stabilisation for Poisson data. \code{"identity"} passes
#'   raw counts \eqn{Y^* = Y} directly, appropriate when counts are large
#'   and rarely zero. Both options affect only GRF factor screening; all
#'   subgroup effect estimates use the correctly specified GLM regardless.
#'   Default: \code{"log"}.
#' @param n.min Integer. Minimum subgroup sample size for a valid split.
#'   Default: \code{60}.
#' @param dmin.grf Numeric. Minimum absolute CATE magnitude (on the
#'   doubly robust score scale) required for the identified subgroup to
#'   be declared valid.  The score equals \code{-CATE(H)}, so a positive
#'   threshold requires the subgroup CATE to be sufficiently negative
#'   (treatment hurts).  Default: \code{0}.
#' @param frac.tau Numeric. Fraction of the sample used for the GRF horizon
#'   (time-horizon analogue for non-survival outcomes). Default: \code{0.5}.
#' @param maxdepth Integer. Maximum depth for GRF policy trees. Default: \code{2}.
#' @param RCT Logical. Whether data come from a randomised trial. When
#'   \code{TRUE}, propensity scores are fixed at 0.5. Default: \code{TRUE}.
#' @param sg.criterion Character. Subgroup selection criterion. One of
#'   \code{"mDiff"} (maximum score difference) or \code{"Nsg"} (largest
#'   valid subgroup). Default: \code{"mDiff"}.
#' @param conf.level Numeric. Confidence level for subgroup effect estimates.
#'   Default: \code{0.95}.
#' @param seedit Integer. Random seed for GRF. Default: \code{8316951}.
#' @param return_selected_cuts_only Logical. If \code{TRUE}, returns only
#'   cuts from the tree depth that identified the selected subgroup. If
#'   \code{FALSE} (default), returns all cuts from all fitted trees.
#'   Matches the interface of \code{\link{grf.subg.harm.survival}}.
#'   Default: \code{FALSE}.
#' @param adverse_outcome Logical. If \code{TRUE}, a higher outcome value
#'   indicates a worse result (e.g., adverse event count), so the "harm"
#'   subgroup is where treatment increases the outcome. If \code{FALSE}
#'   (default), higher values are better. Default: \code{FALSE}.
#' @param tune_grf Logical. If \code{TRUE}, enables cross-validated
#'   hyperparameter tuning for the causal forest via
#'   \code{tune.parameters = "all"}, which tunes \code{min.node.size},
#'   \code{mtry}, \code{sample.fraction}, \code{alpha}, and
#'   \code{imbalance.penalty}. If \code{FALSE} (default), all GRF
#'   parameters remain at their defaults, preserving existing behavior.
#'   Most impactful for observational data where nuisance models benefit
#'   from tuning; defaults are near-optimal for RCTs (Dandl et al., 2024).
#' @param grf_selection Character, one of \code{"tree"} (default) or
#'   \code{"frontier"}. \code{"tree"} uses the GRF policy tree (standard
#'   behavior). \code{"frontier"} is an \strong{experimental} alternative that,
#'   on the same doubly-robust scores, enumerates single-covariate thresholds
#'   (and depth-2 covariate-pair conjunctions when \code{maxdepth >= 2}), takes
#'   the Pareto frontier of (harm-effect, size), and selects one subgroup under
#'   \code{frontier_rule}. The selected subgroup is always a single conjunction
#'   (never the disjunctive union the tree path can return). Provided for
#'   comparison and exploration, not as a recommended default.
#' @param frontier_rule Character, one of \code{"effMaxSG"} (default),
#'   \code{"eff"}, \code{"maxSG"}, \code{"minSG"}, or \code{"effMinSG"}; the rule
#'   applied to the frontier when \code{grf_selection = "frontier"}. See
#'   \code{\link{grf.subg.harm.survival}}.
#' @param effect_neighborhood Numeric in (0, 1); relative neighborhood for the
#'   \code{"effMaxSG"} / \code{"effMinSG"} rules. Default 0.10. Used only when
#'   \code{grf_selection = "frontier"}.
#' @param details Logical. Print GRF diagnostic information. Default: \code{FALSE}.
#' @param verbose Logical. Print progress messages. Default: \code{FALSE}.
#'
#' @return A list of class \code{"grf_glm_result"} containing:
#'   \describe{
#'     \item{sg.harm.id}{\strong{Character vector of cut expressions}
#'       defining the identified subgroup (e.g., \code{c("v1=1", "v3=1")}),
#'       length equal to the depth of the selected policy tree.
#'       \code{NULL} if no subgroup was found.  \strong{Not} a per-subject
#'       membership indicator; see the \emph{Field naming collision}
#'       section below.}
#'     \item{data}{The input data frame with a \code{treat.recommend} column
#'       added (\code{0} = in harm/questionable subgroup, \code{1} =
#'       complement).}
#'     \item{tree.cuts}{Named list of GRF policy-tree split points.}
#'     \item{grf_varimp}{Named numeric vector of GRF variable importances.}
#'     \item{outcome.name}{Character: name of the outcome column.}
#'     \item{treat.name}{Character: name of the treatment column.}
#'     \item{id.name}{Character: name of the subject-identifier column.}
#'     \item{confounders.name}{Character vector of covariate column names.}
#'     \item{effect_measure}{Character: the effect measure used.}
#'     \item{outcome_type}{Character: the outcome type.}
#'     \item{overdispersion}{Character: the overdispersion correction applied.}
#'     \item{offset.name}{Character or \code{NULL}: the offset variable name.}
#'     \item{grf_count_transform}{Character: the transform argument as matched
#'       (\code{"log"} or \code{"identity"}).}
#'     \item{grf_y_transform}{Character: transformation actually applied to Y
#'       for GRF screening (\code{"log(Y + 0.5)"}, \code{"identity"}, or
#'       \code{"none"}).}
#'   }
#'
#' @section Field naming collision with forestsearch / subgroup.consistency:
#' The field name \code{sg.harm.id} has \strong{different semantics} on
#' this object versus on the \code{grp.consistency} list returned by
#' \code{\link{subgroup.consistency}} (and nested inside the
#' \code{\link{forestsearch}} result):
#'
#' \tabular{lll}{
#'   \strong{Object} \tab \strong{\code{sg.harm.id} contains} \tab \strong{Length / type} \cr
#'   \code{grf.subg.harm.glm()} result (this function) \tab character vector of cut expressions \tab character, length = depth of selected tree \cr
#'   \code{subgroup.consistency()} result \tab per-subject 0/1 membership indicator \tab integer, length \code{nrow(df)} \cr
#' }
#'
#' \strong{Practical consequence.}  On this object, pasting
#' \code{paste(obj$sg.harm.id, collapse = " & ")} correctly renders
#' the identified subgroup.  On a \code{grp.consistency} list the same
#' expression silently concatenates a long 0/1 indicator and produces
#' output like \code{"0 & 0 & 1 & 0 & ..."}.  When writing code that
#' must handle both object types, dispatch on class first, or use the
#' top-level \code{sg.harm} field that the \code{\link{forestsearch}}
#' result exposes.
#'
#' This naming collision is a documented CRAN-stable API for v0.1.x
#' and v0.2.x; it is expected to be resolved via a deprecation cycle
#' in a future minor release.
#'
#' @details
#' ## GRF Screening for Count Outcomes
#'
#' \code{grf::causal_forest()} is designed for continuous outcomes. For count
#' data, this function transforms \eqn{Y} before passing it to GRF according
#' to the \code{grf_count_transform} argument:
#'
#' \describe{
#'   \item{\code{"log"} (default)}{Applies \code{Y_star = log(Y + 0.5)}.
#'     The +0.5 shift avoids \code{log(0)} for zero-count observations;
#'     the log scale approximates variance stabilisation for Poisson data and
#'     is recommended when any counts are zero or small.}
#'   \item{\code{"identity"}}{Passes raw counts \code{Y_star = Y} directly.
#'     Appropriate when counts are large, rarely zero, and the raw scale is
#'     near-continuous (e.g., daily hospitalisation volumes in a large centre).}
#' }
#'
#' In both cases the transformation applies \emph{only} to GRF factor
#' screening and policy-tree cut generation -- it determines \emph{which
#' variables} and \emph{which cut-points} are candidate subgroup definitions.
#' The actual effect estimates in each subgroup are always computed from a
#' correctly specified Poisson (or quasi-Poisson / negative-binomial) GLM via
#' \code{\link{create_glm_row}}, regardless of \code{grf_count_transform}.
#'
#' Alternatively, setting \code{use_grf = FALSE} in the parent
#' \code{\link{forestsearch}} call routes to LASSO-based factor selection,
#' which supports \code{family = "poisson"} natively via \pkg{glmnet}.
#'
#' ## Offset and Person-Time
#'
#' When \code{offset.name} is supplied, \eqn{\log(\text{exposure}_i)} enters
#' every GLM fit in the subgroup enumeration loop, so that subgroup effects
#' are on the incidence-rate-ratio scale adjusted for differential follow-up.
#'
#' @examples
#' \dontrun{
#' # Simulate count data with heterogeneous treatment effect
#' set.seed(123)
#' n <- 600
#' z1 <- rbinom(n, 1, 0.4)   # binary subgroup factor
#' df_count <- data.frame(
#'   id          = seq_len(n),
#'   treat       = rbinom(n, 1, 0.5),
#'   v1          = as.factor(z1),
#'   v2          = as.factor(rbinom(n, 1, 0.5)),
#'   person_time = runif(n, 0.5, 2.0),
#'   age         = rnorm(n, 60, 10)
#' )
#' # Generate outcome: higher event rate for treated in z1 = 1 subgroup
#' log_rate <- -1.5 +
#'   0.8 * df_count$treat * (df_count$v1 == 1) -
#'   0.3 * df_count$treat * (df_count$v1 == 0) +
#'   log(df_count$person_time)
#' df_count$events <- rpois(n, lambda = exp(log_rate))
#'
#' # Run GRF subgroup identification for count outcome (log transform, default)
#' grf_glm <- grf.subg.harm.glm(
#'   data                = df_count,
#'   confounders.name    = c("v1", "v2", "age"),
#'   outcome.name        = "events",
#'   treat.name          = "treat",
#'   id.name             = "id",
#'   outcome_type        = "count",
#'   effect_measure      = "log_IRR",
#'   offset.name         = "person_time",
#'   overdispersion      = "quasi",
#'   grf_count_transform = "log",      # default: log(Y + 0.5)
#'   n.min               = 60,
#'   verbose             = TRUE
#' )
#' print(grf_glm$sg.harm.id)
#'
#' # Alternative: identity transform (large counts, no zeros)
#' grf_glm_id <- grf.subg.harm.glm(
#'   data                = df_count,
#'   confounders.name    = c("v1", "v2", "age"),
#'   outcome.name        = "events",
#'   treat.name          = "treat",
#'   id.name             = "id",
#'   outcome_type        = "count",
#'   effect_measure      = "log_IRR",
#'   offset.name         = "person_time",
#'   overdispersion      = "none",
#'   grf_count_transform = "identity",  # raw counts passed to GRF
#'   n.min               = 60
#' )
#' }
#'
#' @seealso
#'   \code{\link{grf.subg.harm.survival}} for time-to-event outcomes,
#'   \code{\link{create_glm_row}} for per-subgroup GLM estimation,
#'   \code{\link{glm_effect_profile}} for continuous biomarker profiles.
#'
#' @importFrom grf causal_forest variable_importance
#' @importFrom policytree policy_tree double_robust_scores
#' @importFrom stats complete.cases
#' @export
grf.subg.harm.glm <- function(
    data,
    confounders.name,
    outcome.name,
    treat.name          = "treat",
    id.name             = "id",
    outcome_type        = c("binary", "continuous", "count"),
    effect_measure      = NULL,
    offset.name         = NULL,
    overdispersion      = c("none", "quasi", "negbin"),
    grf_count_transform = c("log", "identity"),
    n.min               = 60L,
    dmin.grf            = 0,
    frac.tau            = 0.5,
    maxdepth            = 2L,
    RCT                 = TRUE,
    sg.criterion        = c("mDiff", "Nsg"),
    conf.level          = 0.95,
    seedit              = 8316951L,
    return_selected_cuts_only = FALSE,
    adverse_outcome     = FALSE,
    tune_grf            = FALSE,
    grf_selection       = c("tree", "frontier"),
    frontier_rule       = c("effMaxSG", "eff", "maxSG", "minSG", "effMinSG"),
    effect_neighborhood = 0.10,
    details             = FALSE,
    verbose             = FALSE
) {

  # ---------------------------------------------------------------------------
  # Argument matching
  # ---------------------------------------------------------------------------
  outcome_type        <- match.arg(outcome_type)
  overdispersion      <- match.arg(overdispersion)
  sg.criterion        <- match.arg(sg.criterion)
  grf_count_transform <- match.arg(grf_count_transform)
  grf_selection       <- match.arg(grf_selection)
  frontier_rule       <- match.arg(frontier_rule)

  # Default effect_measure by outcome type
  if (is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      binary     = "log_OR",
      continuous = "MD",
      count      = "log_IRR"
    )
  }

  # ---------------------------------------------------------------------------
  # Validate count-specific arguments
  # ---------------------------------------------------------------------------
  if (outcome_type == "count") {
    if (effect_measure == "log_IRR" && is.null(offset.name)) {
      warning(
        "outcome_type = 'count' with effect_measure = 'log_IRR' but ",
        "offset.name = NULL. Fitting Poisson GLM without offset -- ",
        "effects will be on the log-rate (not log-IRR) scale. ",
        "Supply offset.name = '<person_time_column>' for incidence rate ratios."
      )
    }
    if (!effect_measure %in% c("log_IRR", "log_RR")) {
      stop(
        "For outcome_type = 'count', effect_measure must be 'log_IRR' or ",
        "'log_RR'. Got: '", effect_measure, "'."
      )
    }
  }

  if (!is.null(offset.name) && !offset.name %in% names(data)) {
    stop("offset.name '", offset.name, "' not found in data.")
  }

  # ---------------------------------------------------------------------------
  # Required columns check
  # ---------------------------------------------------------------------------
  req_cols <- c(outcome.name, treat.name, id.name, confounders.name)
  if (!is.null(offset.name)) req_cols <- c(req_cols, offset.name)
  missing_cols <- setdiff(req_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Columns not found in data: ", paste(missing_cols, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Complete-case subset
  # ---------------------------------------------------------------------------
  cc_cols <- c(outcome.name, treat.name, id.name, confounders.name)
  if (!is.null(offset.name)) cc_cols <- c(cc_cols, offset.name)
  keep <- stats::complete.cases(data[, cc_cols])
  if (sum(keep) < n.min) {
    stop(sprintf(
      "Only %d complete observations; need at least n.min = %d.",
      sum(keep), n.min
    ))
  }
  data <- data[keep, , drop = FALSE]
  n    <- nrow(data)

  if (verbose) {
    message(sprintf(
      "[grf.subg.harm.glm] n = %d | outcome_type = '%s' | effect_measure = '%s'",
      n, outcome_type, effect_measure
    ))
    if (outcome_type == "count") {
      message(sprintf(
        "[grf.subg.harm.glm] offset = '%s' | overdispersion = '%s'",
        if (is.null(offset.name)) "none" else offset.name, overdispersion
      ))
    }
  }

  # ---------------------------------------------------------------------------
  # Build covariate matrix X for GRF
  # ---------------------------------------------------------------------------
  X <- .build_grf_X(data, confounders.name)

  Y_raw   <- data[[outcome.name]]
  W       <- as.numeric(data[[treat.name]])

  # ---------------------------------------------------------------------------
  # GRF screening: outcome transformation for count data
  # ---------------------------------------------------------------------------
  grf_y_transform <- "none"

  if (outcome_type == "count") {
    if (grf_count_transform == "log") {
      Y_grf           <- log(Y_raw + 0.5)
      grf_y_transform <- "log(Y + 0.5)"
      if (verbose) {
        message(
          "[grf.subg.harm.glm] Count outcome: applying log(Y + 0.5) for ",
          "GRF factor screening (grf_count_transform = 'log'). ",
          "Effect estimation uses Poisson/quasi GLM."
        )
      }
    } else {
      # "identity"  -- pass raw counts directly
      Y_grf           <- Y_raw
      grf_y_transform <- "identity"
      if (verbose) {
        message(
          "[grf.subg.harm.glm] Count outcome: passing raw counts to GRF ",
          "(grf_count_transform = 'identity'). Suitable when counts are ",
          "large and rarely zero."
        )
      }
    }
  } else {
    Y_grf <- Y_raw
  }

  # ---------------------------------------------------------------------------
  # Adverse outcome flip: align CATE direction with ForestSearch convention.
  # When adverse_outcome = TRUE, Y is flipped (binary: 1-Y, continuous: -Y)
  # so that positive CATE = treatment increases the BAD outcome = harm.
  # This ensures the most-negative-CATE subgroup is the true "harm" group,
  # matching forestsearch()'s internal GRF pre-screening behavior.
  # ---------------------------------------------------------------------------
  if (adverse_outcome) {
    if (outcome_type == "binary") {
      Y_grf <- 1L - Y_grf
      if (verbose) message("[grf.subg.harm.glm] adverse_outcome=TRUE: using 1-Y for causal_forest()")
    } else if (outcome_type == "continuous") {
      Y_grf <- -Y_grf
      if (verbose) message("[grf.subg.harm.glm] adverse_outcome=TRUE: using -Y for causal_forest()")
    }
    # count: log transform already applied above; flip not needed
  }

  # ---------------------------------------------------------------------------
  # Fit causal forest
  # ---------------------------------------------------------------------------
  W_hat <- if (RCT) rep(0.5, n) else NULL
  tune_arg <- if (tune_grf) "all" else "none"

  set.seed(seedit)
  cs_forest <- tryCatch({
    if (RCT) {
      grf::causal_forest(X, Y_grf, W, W.hat = W_hat, seed = seedit,
                         tune.parameters = tune_arg)
    } else {
      grf::causal_forest(X, Y_grf, W, seed = seedit,
                         tune.parameters = tune_arg)
    }
  }, error = function(e) {
    stop("grf::causal_forest() failed: ", e$message)
  })

  # Variable importance
  vi <- grf::variable_importance(cs_forest)
  names(vi) <- confounders.name
  if (verbose) {
    vi_sorted <- sort(vi, decreasing = TRUE)
    message("[grf.subg.harm.glm] GRF variable importance (top 5):")
    message(paste(sprintf("  %s: %.4f", names(vi_sorted)[seq_len(min(5, length(vi_sorted)))],
                          vi_sorted[seq_len(min(5, length(vi_sorted)))]),
                  collapse = "\n"))
  }

  # ---------------------------------------------------------------------------
  # Doubly-robust scores and policy trees
  # ---------------------------------------------------------------------------
  dr_scores <- tryCatch(
    policytree::double_robust_scores(cs_forest),
    error = function(e) {
      stop("double_robust_scores() failed: ", e$message)
    }
  )

  # DR-score candidate family (depth-1 + depth-2) in (v1,d1,c1,v2,d2,c2) form,
  # exposed for multiplier resampling (MR): it re-selects over the SAME family
  # GRF ranks (frontier) or, under the policy tree, the cut-defined universe.
  .mr_candidates <- tryCatch({
    cc <- .grf_dr_candidates(X, dr_scores, n_min = n.min)
    if (maxdepth >= 2L) {
      cc2 <- .grf_dr_candidates_d2(X, dr_scores, n_min = n.min)
      cc  <- if (is.null(cc)) cc2 else if (is.null(cc2)) cc else rbind(cc, cc2)
    }
    cc
  }, error = function(e) NULL)

  # ---------------------------------------------------------------------------
  # FRONTIER SELECTION (grf_selection = "frontier")
  # ---------------------------------------------------------------------------
  # Alternative to the policy tree: enumerate threshold (and depth-2 pair)
  # candidates on the SAME DR scores, take the Pareto frontier, select under
  # frontier_rule.  dr_scores columns are [control (Gamma_0), treated (Gamma_1)];
  # the frontier's harm-effect = mean(control) - mean(treated) matches the tree
  # path's score = -CATE convention.  When adverse_outcome = TRUE the Y-flip was
  # applied before causal_forest(), so dr_scores are already harm-oriented and
  # need no special handling here.  grf_selection = "tree" skips this entirely.
  if (identical(grf_selection, "frontier")) {
    cand <- .grf_dr_candidates(X, dr_scores, n_min = n.min)
    if (maxdepth >= 2L) {
      cand2 <- .grf_dr_candidates_d2(X, dr_scores, n_min = n.min)
      cand  <- if (is.null(cand)) cand2
               else if (is.null(cand2)) cand else rbind(cand, cand2)
    }
    sel <- .grf_frontier_select(cand, dmin = dmin.grf, rule = frontier_rule,
                                nbhd = effect_neighborhood)
    if (is.null(sel)) {
      if (verbose || details)
        message("[grf.subg.harm.glm] frontier: no eligible subgroup.")
      data$treat.recommend <- 1L
      return(structure(
        list(sg.harm.id = NULL,
             sg_def = list(conjunctions = list(), labels = NULL,
                           definition = NA_character_, is_disjunction = FALSE),
             data = data, tree.cuts = NULL, grf_varimp = vi,
             outcome.name = outcome.name, treat.name = treat.name,
             id.name = id.name, confounders.name = confounders.name,
             effect_measure = effect_measure, outcome_type = outcome_type,
             overdispersion = overdispersion, offset.name = offset.name,
             grf_count_transform = grf_count_transform,
             grf_y_transform = grf_y_transform,
             return_selected_cuts_only = return_selected_cuts_only,
             adverse_outcome = adverse_outcome, selection = "frontier"),
        class = "grf_glm_result"))
    }
    sg_def     <- .grf_sg_def_from_candidate(sel)
    sg_harm_id <- if (!is.null(sg_def$labels))
                    gsub("^\\{|\\}$", "", sg_def$labels) else sg_def$definition
    data$treat.recommend <- .grf_evaluate_subgroup(sg_def, data)
    if (verbose || details) {
      message(sprintf(paste0("[grf.subg.harm.glm] frontier rule=%s dmin=%.3g ",
                             "selected: %s (effect=%.4f, size=%d)"),
                      frontier_rule, dmin.grf, sg_def$definition,
                      sel$effect, sel$size))
    }
    return(structure(
      list(sg.harm.id = sg_harm_id, sg_def = sg_def, data = data,
           tree.cuts = NULL, grf_varimp = vi,
           outcome.name = outcome.name, treat.name = treat.name,
           id.name = id.name, confounders.name = confounders.name,
           effect_measure = effect_measure, outcome_type = outcome_type,
           overdispersion = overdispersion, offset.name = offset.name,
           grf_count_transform = grf_count_transform,
           grf_y_transform = grf_y_transform,
           return_selected_cuts_only = return_selected_cuts_only,
           adverse_outcome = adverse_outcome,
           trees = NULL, dr_scores = dr_scores, cs_forest = cs_forest,
           cate_sg = -sel$effect, cate_sgc = NA_real_,
           best_depth = if (is.na(sel$v2)) 1L else 2L,
           selection = "frontier",
           frontier = .grf_mark_frontier(cand),
           candidates = .mr_candidates),
      class = "grf_glm_result"))
  }

  trees  <- vector("list", maxdepth)
  values <- vector("list", maxdepth)
  n_max  <- n  # exclude full population

  for (d in seq_len(maxdepth)) {
    trees[[d]] <- tryCatch(
      policytree::policy_tree(X, dr_scores, depth = d),
      error = function(e) NULL
    )
    if (!is.null(trees[[d]])) {
      leaf_action <- predict(trees[[d]], X)
      sg_idx      <- which(leaf_action == 1L)  # action 1 = recommend control
      sg_n        <- length(sg_idx)
      if (sg_n >= n.min && sg_n < n_max) {
        # Treatment effect contrast in the identified subgroup.
        # DR scores: column 1 = Gamma_0 (control reward),
        #            column 2 = Gamma_1 (treated reward).
        # CATE = E[Y(1) - Y(0)] = Gamma_1 - Gamma_0.
        # score = -CATE_sg: positive when treatment hurts the subgroup
        #   (Gamma_0 > Gamma_1), matching the survival GRF convention
        #   where diff = control_mean - treated_mean.
        cate_sg  <- mean(dr_scores[sg_idx, 2L] - dr_scores[sg_idx, 1L])
        cate_sgc <- mean(dr_scores[-sg_idx, 2L] - dr_scores[-sg_idx, 1L])
        values[[d]] <- list(
          depth    = d,
          n_sg     = sg_n,
          score    = -cate_sg,  # positive = treatment hurts subgroup
          cate_sg  = cate_sg,
          cate_sgc = cate_sgc,
          sg_idx   = sg_idx,
          leaf_action = leaf_action
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Select best subgroup
  # ---------------------------------------------------------------------------
  best <- .select_best_sg_glm(values, sg.criterion, dmin.grf, n, n.min)

  if (is.null(best)) {
    if (verbose) {
      message("[grf.subg.harm.glm] No valid subgroup identified.")
    }
    data$treat.recommend <- 1L
    return(structure(
      list(
        sg.harm.id          = NULL,
        data                = data,
        tree.cuts           = NULL,
        grf_varimp          = vi,
        outcome.name        = outcome.name,
        treat.name          = treat.name,
        id.name             = id.name,
        confounders.name    = confounders.name,
        effect_measure      = effect_measure,
        outcome_type        = outcome_type,
        overdispersion      = overdispersion,
        offset.name         = offset.name,
        grf_count_transform = grf_count_transform,
        grf_y_transform     = grf_y_transform,
        return_selected_cuts_only = return_selected_cuts_only,
        adverse_outcome     = adverse_outcome
      ),
      class = "grf_glm_result"
    ))
  }

  # ---------------------------------------------------------------------------
  # Extract tree cuts and subgroup definition
  # ---------------------------------------------------------------------------
  # The GLM subgroup is the set of action-1 (recommend-control) leaves, which
  # can span MORE THAN ONE leaf -- a union of axis-aligned boxes.  Build the
  # exact root-to-leaf rule for each such leaf, with correct directions, via
  # .grf_build_subgroup_definition(); the result is a disjunction of
  # conjunctions when the subgroup spans multiple leaves.  The earlier
  # .build_sg_harm_id() emitted every internal-node split as "<=", ignoring
  # path and direction, so its definition did not match the subgroup GRF
  # actually selected.
  best_tree    <- trees[[best$depth]]
  best_node_id <- predict(best_tree, X, type = "node.id")
  harm_leaves  <- sort(unique(best_node_id[best$sg_idx]))
  tree_cuts    <- .extract_tree_cuts(best_tree, X, confounders.name)
  sg_def       <- .grf_build_subgroup_definition(best_tree, harm_leaves)
  # sg.harm.id: component-label vector for a single conjunction (back-compatible
  # shape), or the disjunctive definition string for a multi-leaf union.
  sg_harm_id   <- if (!is.null(sg_def$labels))
                    gsub("^\\{|\\}$", "", sg_def$labels) else sg_def$definition

  data$treat.recommend <- ifelse(seq_len(n) %in% best$sg_idx, 0L, 1L)

  if (verbose) {
    n_harm <- sum(data$treat.recommend == 0L)
    message(sprintf(
      "[grf.subg.harm.glm] Subgroup found: n(H) = %d (%.1f%%) | definition: %s",
      n_harm, 100 * n_harm / n,
      paste(sg_harm_id, collapse = " & ")
    ))
    message(sprintf(
      "[grf.subg.harm.glm] CATE(H) = %.4f, CATE(Hc) = %.4f, score = %.4f (depth %d)",
      best$cate_sg, best$cate_sgc, best$score, best$depth
    ))
  }

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  structure(
    list(
      sg.harm.id          = sg_harm_id,
      sg_def              = sg_def,
      data                = data,
      tree.cuts           = tree_cuts,
      grf_varimp          = vi,
      outcome.name        = outcome.name,
      treat.name          = treat.name,
      id.name             = id.name,
      confounders.name    = confounders.name,
      effect_measure      = effect_measure,
      outcome_type        = outcome_type,
      overdispersion      = overdispersion,
      offset.name         = offset.name,
      grf_count_transform = grf_count_transform,
      grf_y_transform     = grf_y_transform,
      return_selected_cuts_only = return_selected_cuts_only,
      adverse_outcome     = adverse_outcome,
      trees               = trees,
      dr_scores           = dr_scores,
      cs_forest           = cs_forest,
      cate_sg             = best$cate_sg,
      cate_sgc            = best$cate_sgc,
      best_depth          = best$depth,
      candidates          = .mr_candidates   # DR-candidate family for MR
    ),
    class = "grf_glm_result"
  )
}


# =============================================================================
# create_glm_row()   -- per-subgroup GLM effect estimate
# =============================================================================

#' Compute GLM Effect Estimate for a Single Subgroup
#'
#' Fits a GLM within the specified subgroup and returns a one-row data frame
#' containing the effect estimate, confidence interval, sample size, and event
#' count. Supports binary, continuous, and count (with offset) outcomes.
#'
#' @param df Data frame for the subgroup (already subset).
#' @param outcome.name Character. Name of the outcome column.
#' @param treat.name Character. Name of the treatment indicator.
#' @param effect_measure Character. Effect measure; see
#'   \code{\link{grf.subg.harm.glm}} for options.
#' @param offset.name Character or \code{NULL}. Name of the person-time column.
#'   When supplied, \code{log(offset)} is added as an offset term.
#'   Default: \code{NULL}.
#' @param overdispersion Character. One of \code{"none"}, \code{"quasi"},
#'   \code{"negbin"}. Default: \code{"none"}.
#' @param conf.level Numeric. Confidence level. Default: \code{0.95}.
#' @param min_arm_n Integer. Minimum observations per arm; returns \code{NA}
#'   row if either arm is below this threshold. Default: \code{5L}.
#' @param min_arm_events Integer. For binary/count outcomes, minimum events
#'   per arm; returns \code{NA} row if check fails. Default: \code{3L}.
#' @param verbose Logical. Default: \code{FALSE}.
#'
#' @return A one-row data frame with columns:
#'   \code{N}, \code{n_treat}, \code{n_control},
#'   \code{events_treat}, \code{events_control} (for binary/count),
#'   \code{est}, \code{lower}, \code{upper}, \code{se}, \code{effect_measure}.
#'
#' @seealso \code{\link{grf.subg.harm.glm}}
#'
#' @importFrom stats glm binomial poisson quasipoisson gaussian coef vcov qnorm complete.cases as.formula
#' @export
create_glm_row <- function(
    df,
    outcome.name,
    treat.name      = "treat",
    effect_measure  = "log_OR",
    offset.name     = NULL,
    overdispersion  = "none",
    conf.level      = 0.95,
    min_arm_n       = 5L,
    min_arm_events  = 3L,
    verbose         = FALSE
) {

  z_mult <- stats::qnorm(1 - (1 - conf.level) / 2)

  # NA row template
  na_row <- data.frame(
    N               = nrow(df),
    n_treat         = sum(df[[treat.name]] == 1L, na.rm = TRUE),
    n_control       = sum(df[[treat.name]] == 0L, na.rm = TRUE),
    events_treat    = NA_integer_,
    events_control  = NA_integer_,
    est             = NA_real_,
    lower           = NA_real_,
    upper           = NA_real_,
    se              = NA_real_,
    effect_measure  = effect_measure,
    stringsAsFactors = FALSE
  )

  # Basic sparsity check
  n_treat   <- sum(df[[treat.name]] == 1L, na.rm = TRUE)
  n_control <- sum(df[[treat.name]] == 0L, na.rm = TRUE)
  if (n_treat < min_arm_n || n_control < min_arm_n) return(na_row)

  Y <- df[[outcome.name]]
  W <- df[[treat.name]]

  # Event counts
  if (effect_measure %in% c("log_OR", "log_RR", "RD", "log_IRR")) {
    ev_treat   <- sum(Y[W == 1L], na.rm = TRUE)
    ev_control <- sum(Y[W == 0L], na.rm = TRUE)
    na_row$events_treat   <- as.integer(round(ev_treat))
    na_row$events_control <- as.integer(round(ev_control))
    if (ev_treat < min_arm_events || ev_control < min_arm_events) {
      return(na_row)
    }
  } else {
    na_row$events_treat   <- NA_integer_
    na_row$events_control <- NA_integer_
  }

  # Build GLM formula with optional offset
  if (!is.null(offset.name) && offset.name %in% names(df)) {
    log_off_col <- ".log_offset_cr"
    df[[log_off_col]] <- log(df[[offset.name]])
    fml <- stats::as.formula(
      paste(outcome.name, "~", treat.name, "+ offset(", log_off_col, ")")
    )
  } else {
    log_off_col <- NULL
    fml <- stats::as.formula(paste(outcome.name, "~", treat.name))
  }

  # Select family
  family_obj <- .cr_resolve_family(effect_measure, overdispersion)
  is_negbin  <- identical(family_obj, "negbin")

  fit <- tryCatch({
    if (is_negbin) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("MASS required for negbin")
      }
      MASS::glm.nb(fml, data = df)
    } else {
      stats::glm(fml, data = df, family = family_obj)
    }
  }, error = function(e) {
    if (verbose) message("[create_glm_row] GLM failed: ", e$message)
    NULL
  })

  if (is.null(fit)) return(na_row)

  # Extract treat coefficient
  coef_names <- names(stats::coef(fit))
  treat_idx  <- grep(paste0("^", treat.name), coef_names)
  if (length(treat_idx) == 0L) return(na_row)

  est    <- stats::coef(fit)[treat_idx]
  vcov_m <- .cr_get_vcov(fit, effect_measure, overdispersion)
  se     <- sqrt(max(vcov_m[treat_idx, treat_idx], 0))

  data.frame(
    N              = nrow(df),
    n_treat        = n_treat,
    n_control      = n_control,
    events_treat   = na_row$events_treat,
    events_control = na_row$events_control,
    est            = as.numeric(est),
    lower          = as.numeric(est - z_mult * se),
    upper          = as.numeric(est + z_mult * se),
    se             = as.numeric(se),
    effect_measure = effect_measure,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

# Build numeric covariate matrix for GRF from data frame
#' @noRd
.build_grf_X <- function(data, confounders.name) {
  X <- data[, confounders.name, drop = FALSE]

  # Convert factors/characters to numeric.

  # IMPORTANT: Must use the same conversion as evaluate_comparison()

  # in forestsearch_helpers.R.  For factors with all-numeric levels

  # (e.g., levels "0"/"1"), use as.numeric(as.character()) to get

  # the actual level values {0, 1}.  For non-numeric levels (e.g.,

  # "male"/"female"), fall back to integer codes {1, 2, ...}.
  #
  # If these scales diverge, GRF split_values become inconsistent
  # with downstream cut evaluation, causing valid cuts to be
  # silently dropped (e.g., "z1 <= 1" is trivially TRUE on
  # the {0, 1} scale but meaningful on the {1, 2} scale).
  for (v in names(X)) {
    if (is.factor(X[[v]]) || is.character(X[[v]])) {
      lvls <- if (is.factor(X[[v]])) levels(X[[v]]) else unique(X[[v]])
      if (!anyNA(suppressWarnings(as.numeric(lvls)))) {
        # All-numeric levels: preserve original values
        X[[v]] <- as.numeric(as.character(X[[v]]))
      } else {
        # Non-numeric levels: integer codes
        X[[v]] <- as.integer(as.factor(X[[v]]))
      }
    }
  }
  as.matrix(X)
}


# Select best subgroup from policy tree values
#' @noRd
.select_best_sg_glm <- function(values, sg.criterion, dmin.grf, n, n.min) {
  valid <- Filter(Negate(is.null), values)
  if (length(valid) == 0L) return(NULL)

  if (sg.criterion == "mDiff") {
    # Pick the tree with the largest absolute harm signal (-CATE in
    # the subgroup), subject to dmin.grf threshold.  This matches
    # the survival GRF convention where diff = control_mean - treated_mean.
    scores <- vapply(valid, `[[`, numeric(1), "score")
    best_v <- valid[[which.max(scores)]]
    if (max(scores) < dmin.grf) return(NULL)
    return(best_v)
  }

  if (sg.criterion == "Nsg") {
    # Pick the valid tree with the largest subgroup size
    ns <- vapply(valid, `[[`, numeric(1), "n_sg")
    return(valid[[which.max(ns)]])
  }


  NULL
}


# Extract variable cut-points from a policy tree
#' @noRd
.extract_tree_cuts <- function(tree, X, confounders.name) {
  if (is.null(tree)) return(NULL)
  tryCatch({
    nodes <- tree$nodes
    cuts  <- list()
    for (nd in nodes) {
      if (!is.null(nd$split_variable) && !is.null(nd$split_value)) {
        var_name <- confounders.name[nd$split_variable]
        cuts[[var_name]] <- c(cuts[[var_name]], nd$split_value)
      }
    }
    lapply(cuts, unique)
  }, error = function(e) NULL)
}


# Superseded by .grf_build_subgroup_definition() in grf_subgroup_labels.R,
# which builds the correct root-to-leaf rule(s) with proper directions
# (the former .build_sg_harm_id() emitted every node split as "<=", ignoring
# path and direction).  Removed in this release.


# Resolve GLM family for create_glm_row()
#' @noRd
.cr_resolve_family <- function(effect_measure, overdispersion) {
  switch(effect_measure,
    log_OR  = stats::binomial(link = "logit"),
    RD      = stats::binomial(link = "identity"),
    MD      = stats::gaussian(link = "identity"),
    log_RR  = stats::poisson(link = "log"),
    log_IRR = {
      if (overdispersion == "quasi")  return(stats::quasipoisson(link = "log"))
      if (overdispersion == "negbin") return("negbin")
      stats::poisson(link = "log")
    },
    stop("Unknown effect_measure: ", effect_measure)
  )
}


# Get vcov matrix, with sandwich correction when appropriate
#' @noRd
.cr_get_vcov <- function(fit, effect_measure, overdispersion) {
  use_sandwich <- (effect_measure == "log_RR" && overdispersion == "none")
  if (use_sandwich && requireNamespace("sandwich", quietly = TRUE)) {
    return(sandwich::vcovHC(fit, type = "HC0"))
  }
  stats::vcov(fit)
}
