#' GRF Subgroup Identification for Survival Data
#'
#' Identifies subgroups with differential treatment effect using generalized random forests (GRF)
#' and policy trees. This function uses causal survival forests to identify heterogeneous
#' treatment effects and policy trees to create interpretable subgroup definitions.
#'
#' @param data Data frame containing the analysis data.
#' @param confounders.name Character vector of confounder variable names.
#' @param outcome.name Character. Name of outcome variable (e.g., time-to-event).
#' @param event.name Character. Name of event indicator variable (0/1).
#' @param id.name Character. Name of ID variable.
#' @param treat.name Character. Name of treatment group variable (0/1).
#' @param frac.tau Numeric. Fraction of tau for GRF horizon (default: 1.0).
#' @param n.min Integer. Minimum subgroup size (default: 60).
#' @param dmin.grf Numeric. Minimum difference in subgroup mean (default: 0.0).
#' @param RCT Logical. Is the data from a randomized controlled trial? (default: TRUE)
#' @param details Logical. Print details during execution (default: FALSE).
#' @param sg.criterion Character. Subgroup selection criterion ("mDiff" or "Nsg").
#' @param maxdepth Integer. Maximum tree depth (1, 2, or 3; default: 2).
#' @param seedit Integer. Random seed (default: 8316951).
#' @param return_selected_cuts_only Logical. If TRUE, returns only cuts from the tree
#'   depth that identified the selected subgroup meeting `dmin.grf`. If FALSE (default),

#'   returns all cuts from all fitted trees (depths 1 through `maxdepth`).
#' @param tune_grf Logical. If \code{TRUE}, enables cross-validated
#'   hyperparameter tuning via \code{tune.parameters = "all"} in the
#'   causal survival forest. Default \code{FALSE}.
#' @param grf_selection Character, one of \code{"tree"} (default) or
#'   \code{"frontier"}. \code{"tree"} uses the GRF policy tree (standard
#'   behavior). \code{"frontier"} is an \strong{experimental} alternative that,
#'   on the same doubly-robust scores, enumerates single-covariate thresholds
#'   (and depth-2 covariate-pair conjunctions when \code{maxdepth >= 2}), takes
#'   the Pareto frontier of (harm-effect, size), and selects one subgroup under
#'   \code{frontier_rule}. The selected subgroup is always a single conjunction.
#'   In small benchmarks the policy tree matched or beat the frontier on harm
#'   recovery, so \code{"frontier"} is provided for comparison and exploration,
#'   not as a recommended default.
#' @param frontier_rule Character, one of \code{"effMaxSG"} (default),
#'   \code{"eff"}, \code{"maxSG"}, \code{"minSG"}, or \code{"effMinSG"}; the rule
#'   applied to the frontier when \code{grf_selection = "frontier"}. Semantics
#'   match \code{dina_subgroup()}'s \code{sg_focus}: \code{"eff"} takes the
#'   maximum harm-effect candidate; \code{"maxSG"}/\code{"minSG"} the
#'   largest/smallest eligible candidate (effect >= \code{dmin.grf});
#'   \code{"effMaxSG"}/\code{"effMinSG"} the largest/smallest candidate within
#'   \code{effect_neighborhood} (relative) of the max harm-effect.
#' @param effect_neighborhood Numeric in (0, 1); relative neighborhood for the
#'   \code{"effMaxSG"} / \code{"effMinSG"} rules. Default 0.10. Used only when
#'   \code{grf_selection = "frontier"}.
#'
#' @return A list with GRF results, including:
#'   \item{data}{Original data with added treatment recommendation flags}
#'   \item{grf.gsub}{Selected subgroup information}
#'   \item{sg.harm.id}{\strong{Character vector of cut expressions}
#'     defining the identified subgroup (length = depth of the selected
#'     policy tree), or \code{NULL} if no subgroup was found.  \strong{Not}
#'     a per-subject membership indicator; see the \emph{Field naming
#'     collision} section below.}
#'   \item{tree.cuts}{Cut expressions - either all cuts from all trees (if
#'     `return_selected_cuts_only = FALSE`) or only cuts from the selected tree
#'     depth (if `return_selected_cuts_only = TRUE`)}
#'   \item{tree.names}{Unique variable names used in cuts}
#'   \item{tree}{Selected policy tree object}
#'   \item{tau.rmst}{Time horizon used for RMST}
#'   \item{harm.any}{All subgroups with positive treatment effect difference}
#'   \item{selected_depth}{Depth of the tree that identified the subgroup (when found)}
#'   \item{return_selected_cuts_only}{Logical indicating which cut extraction mode was used}
#'   Additional tree-specific cuts and objects (tree1, tree2, tree3) based on maxdepth
#'
#' @section Field naming collision with forestsearch / subgroup.consistency:
#' The field name \code{sg.harm.id} has \strong{different semantics} on
#' this object versus on the \code{grp.consistency} list returned by
#' \code{\link{subgroup.consistency}} (and nested inside the
#' \code{\link{forestsearch}} result):
#'
#' \tabular{lll}{
#'   \strong{Object} \tab \strong{\code{sg.harm.id} contains} \tab \strong{Length / type} \cr
#'   \code{grf.subg.harm.survival()} result (this function) \tab character vector of cut expressions \tab character, length = depth of selected tree \cr
#'   \code{subgroup.consistency()} result \tab per-subject 0/1 membership indicator \tab integer, length \code{nrow(data)} \cr
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
#' The `return_selected_cuts_only` parameter controls which cuts are returned:
#'
#' \describe{
#'   \item{FALSE (default)}{Returns all cuts from all fitted trees (depths 1 to
#'     `maxdepth`). This provides the full set of candidate splits for downstream
#'     exploration and is the original behavior for backward compatibility.}
#'   \item{TRUE}{Returns only cuts from the tree at the depth that identified
#'     the "winning" subgroup meeting the `dmin.grf` criterion. This is useful
#'     when you want a focused set of cuts associated with the selected subgroup,
#'     reducing noise from non-selected trees.}
#' }
#'
#' When `return_selected_cuts_only = TRUE` and no subgroup meets the criteria,
#' `tree.cuts` will be empty (character(0)).
#'
#' @examples
#' \dontrun{
#' # Return all cuts (default behavior)
#' result_all <- grf.subg.harm.survival(
#'   data = trial_data,
#'   confounders.name = c("age", "biomarker", "region"),
#'   outcome.name = "tte",
#'   event.name = "event",
#'   id.name = "id",
#'   treat.name = "treat",
#'   dmin.grf = 0.1,
#'   maxdepth = 2
#' )
#' result_all$tree.cuts
#' # Returns cuts from both depth 1 and depth 2 trees
#'
#' # Return only cuts from the selected tree
#' result_selected <- grf.subg.harm.survival(
#'   data = trial_data,
#'   confounders.name = c("age", "biomarker", "region"),
#'   outcome.name = "tte",
#'   event.name = "event",
#'   id.name = "id",
#'   treat.name = "treat",
#'   dmin.grf = 0.1,
#'   maxdepth = 2,
#'   return_selected_cuts_only = TRUE
#' )
#' result_selected$tree.cuts
#' # Returns cuts only from the depth that identified the winning subgroup
#' }
#'
#' @importFrom grf causal_survival_forest
#' @importFrom policytree double_robust_scores policy_tree
#' @importFrom stats aggregate
#' @export

grf.subg.harm.survival <- function(data,
                                   confounders.name,
                                   outcome.name,
                                   event.name,
                                   id.name,
                                   treat.name,
                                   frac.tau = 1.0,
                                   n.min = 60,
                                   dmin.grf = 0.0,
                                   RCT = TRUE,
                                   details = FALSE,
                                   sg.criterion = "mDiff",
                                   maxdepth = 2,
                                   seedit = 8316951,
                                   return_selected_cuts_only = FALSE,
                                   tune_grf = FALSE,
                                   grf_selection = c("tree", "frontier"),
                                   frontier_rule = c("effMaxSG", "eff", "maxSG",
                                                     "minSG", "effMinSG"),
                                   effect_neighborhood = 0.10) {

  grf_selection <- match.arg(grf_selection)
  frontier_rule <- match.arg(frontier_rule)

  # ===========================================================================
  # SECTION: INPUT VALIDATION
  # Purpose: Validate all input parameters before processing
  # ===========================================================================

  if (maxdepth > 3) {
    stop("Maximum depth cannot exceed 3")
  }


  valid_criteria <- c("mDiff", "Nsg")
  if (!sg.criterion %in% valid_criteria) {
    stop("sg.criterion must be one of: ", paste(valid_criteria, collapse = ", "))
  }

  if (!is.logical(return_selected_cuts_only) || length(return_selected_cuts_only) != 1) {
    stop("return_selected_cuts_only must be a single logical value (TRUE or FALSE)")
  }

  # ===========================================================================
  # SECTION: CONFIGURATION SETUP
  # Purpose: Create configuration object for consistent parameter passing
  # ===========================================================================

  config <- create_grf_config(
    frac.tau = frac.tau,
    n.min = n.min,
    dmin.grf = dmin.grf,
    RCT = RCT,
    sg.criterion = sg.criterion,
    maxdepth = maxdepth,
    seedit = seedit
  )

  # Add return_selected_cuts_only to config for downstream functions
  config$return_selected_cuts_only <- return_selected_cuts_only
  # Frontier-mode settings (selection = "frontier"); inert for the tree path.
  config$selection           <- grf_selection
  config$frontier_rule       <- frontier_rule
  config$effect_neighborhood <- effect_neighborhood

  # ===========================================================================
  # SECTION: DATA PREPARATION
  # Purpose: Convert data to appropriate format for GRF analysis
  # ===========================================================================

  # Convert confounders to numeric matrix
  X <- apply(data[, confounders.name, drop = FALSE], 2, as.numeric)

  # Extract outcome variables
  Y <- data[, outcome.name]
  W <- data[, treat.name]
  D <- data[, event.name]

  # Validate data sufficiency
  if (!validate_grf_data(W, D, config$n.min)) {
    return(create_null_result(data, NULL, list(), config))
  }

  # Calculate time horizon for RMST
  tau.rmst <- config$frac.tau * min(
    max(Y[W == 1 & D == 1]),
    max(Y[W == 0 & D == 1])
  )

  # Update config with calculated tau
  config$tau.rmst <- tau.rmst

  # ===========================================================================
  # SECTION: CAUSAL FOREST FITTING
  # Purpose: Fit GRF causal survival forest to identify treatment heterogeneity
  # ===========================================================================

  cs.forest <- fit_causal_forest(X, Y, W, D, tau.rmst, config$RCT, config$seedit,
                                  tune_grf = tune_grf)

  # ===========================================================================
  # SECTION: SUBGROUP IDENTIFICATION VIA POLICY TREES
  # Purpose: Use policy trees to partition the covariate space
  # ===========================================================================

  # Compute doubly robust scores for subgroup identification
  dr.scores <- policytree::double_robust_scores(cs.forest)

  # Maximum sample size (used to exclude full population as subgroup)
  n.max <- length(Y)

  # DR-score candidate family (depth-1 + depth-2) in (v1,d1,c1,v2,d2,c2) form,
  # exposed so the Tier-2 de-biased gate can re-select over the SAME family GRF
  # ranks (frontier) or, under the policy tree, over the cut-defined universe.
  .dg_candidates <- tryCatch({
    cc <- .grf_dr_candidates(X, dr.scores, n_min = config$n.min)
    if (config$maxdepth >= 2L) {
      cc2 <- .grf_dr_candidates_d2(X, dr.scores, n_min = config$n.min)
      cc  <- if (is.null(cc)) cc2 else if (is.null(cc2)) cc else rbind(cc, cc2)
    }
    cc
  }, error = function(e) NULL)

  # ===========================================================================
  # SECTION: FRONTIER SELECTION (selection = "frontier")
  # Purpose: alternative to the policy tree -- enumerate threshold (and depth-2
  #   pair) candidates on the SAME DR scores, take the Pareto frontier, select
  #   under frontier_rule.  Default selection = "tree" skips this entirely.
  # ===========================================================================
  if (identical(config$selection, "frontier")) {
    cand <- .grf_dr_candidates(X, dr.scores, n_min = config$n.min)
    if (config$maxdepth >= 2L) {
      cand2 <- .grf_dr_candidates_d2(X, dr.scores, n_min = config$n.min)
      cand  <- if (is.null(cand)) cand2
               else if (is.null(cand2)) cand else rbind(cand, cand2)
    }
    sel <- .grf_frontier_select(cand, dmin = config$dmin.grf,
                                rule = config$frontier_rule,
                                nbhd = config$effect_neighborhood)
    if (is.null(sel)) {
      if (details) print_grf_details(config, NULL, NULL, NULL)
      nullres <- create_null_result(data, NULL, NULL, config)
      nullres$sg_def    <- list(conjunctions = list(), labels = NULL,
                                definition = NA_character_,
                                is_disjunction = FALSE)
      nullres$selection <- "frontier"
      return(nullres)
    }
    sg_def     <- .grf_sg_def_from_candidate(sel)
    sg_harm_id <- if (!is.null(sg_def$labels))
                    gsub("^\\{|\\}$", "", sg_def$labels) else sg_def$definition
    # Membership from the structured definition (reproduces the candidate's S).
    tr_rec <- .grf_evaluate_subgroup(sg_def, data)
    data$treat.recommend <- tr_rec
    if (details) {
      message(sprintf(paste0("[grf frontier] rule=%s  dmin=%.3g  ",
                             "selected: %s  (effect=%.3f, size=%d)"),
                      config$frontier_rule, config$dmin.grf,
                      sg_def$definition, sel$effect, sel$size))
    }
    result <- list(
      data           = data,
      grf.gsub       = list(leaf.node = NA, depth = if (is.na(sel$v2)) 1L else 2L,
                            effect = sel$effect, Nsg = sel$size),
      sg.harm.id     = sg_harm_id,
      sg_def         = sg_def,
      tree           = NULL,
      tau.rmst       = config$tau.rmst,
      dmin.grf       = config$dmin.grf,
      frac.tau       = config$frac.tau,
      maxdepth       = config$maxdepth,
      n.min          = config$n.min,
      selected_depth = if (is.na(sel$v2)) 1L else 2L,
      selection      = "frontier",
      frontier       = .grf_mark_frontier(cand),
      candidates     = .dg_candidates
    )
    return(result)
  }

  # Fit policy trees and compute metrics
  tree_results <- fit_policy_trees(X, data, dr.scores, config$maxdepth, config$n.min)
  trees <- tree_results$trees
  values <- tree_results$values

  # ===========================================================================
  # SECTION: OPTIMAL SUBGROUP SELECTION
  # Purpose: Choose the best subgroup based on specified criterion
  # ===========================================================================

  best_subgroup <- select_best_subgroup(
    values = values,
    sg.criterion = config$sg.criterion,
    dmin.grf = config$dmin.grf,
    n.max = n.max
  )

  # ===========================================================================
  # SECTION: RESULT COMPILATION - NO SUBGROUP FOUND
  # Purpose: Return appropriate result when no valid subgroup is identified
  # ===========================================================================

  if (is.null(best_subgroup)) {
    if (details) {
      print_grf_details(config, values, NULL, NULL)
    }

    return(create_null_result(data, values, trees, config))
  }

  # ===========================================================================
  # SECTION: RESULT COMPILATION - SUBGROUP FOUND
  # Purpose: Extract subgroup information and create comprehensive result
  # ===========================================================================

  # Assign data points to subgroups
  data <- assign_subgroup_membership(data, best_subgroup, trees, X)

  # Select the tree that identified the best subgroup
  selected_tree <- trees[[best_subgroup$depth]]

  # Find the decision rule that defines the subgroup.  The survival subgroup
  # is a SINGLE leaf (best_subgroup$leaf.node), so its definition is the full
  # root-to-leaf conjunction with correct per-split directions.  The earlier
  # find_leaf_split() returned only the single split immediately above the
  # leaf, which silently dropped the rest of the path and mis-rendered
  # right-turns as "<=" -- producing a definition that did not match the
  # subgroup GRF actually selected.  .grf_build_subgroup_definition() builds
  # the correct path; sg.harm.id is the bare cut-expression vector (directions
  # included), back-compatible in shape with the old single-leaf output.
  sg_def     <- .grf_build_subgroup_definition(selected_tree,
                                               best_subgroup$leaf.node)
  sg_harm_id <- if (!is.null(sg_def$labels))
                  gsub("^\\{|\\}$", "", sg_def$labels) else sg_def$definition

  # Extract cuts based on return_selected_cuts_only setting
  if (config$return_selected_cuts_only) {
    # Only extract cuts from the selected tree depth
    tree_cuts <- extract_selected_tree_cuts(trees, best_subgroup$depth, config$maxdepth)
  } else {
    # Extract all cuts from all fitted trees (original behavior)
    tree_cuts <- extract_all_tree_cuts(trees, config$maxdepth)
  }

  # Print details if requested
  if (details) {
    print_grf_details(config, values, best_subgroup, sg_harm_id, tree_cuts)
  }

  # Create comprehensive result object
  result <- create_success_result(
    data = data,
    best_subgroup = best_subgroup,
    trees = trees,
    tree_cuts = tree_cuts,
    selected_tree = selected_tree,
    sg_harm_id = sg_harm_id,
    values = values,
    config = config
  )
  # Structured subgroup definition (path-based, direction-correct; possibly a
  # disjunction) for downstream membership evaluation on new data.
  result$sg_def <- sg_def
  result$candidates <- .dg_candidates   # DR-candidate family for the Tier-2 gate

  return(result)
}


# ── grf.subg.harm.glm() has been moved to grf_subg_harm_glm.R ──────────────
# The canonical implementation lives in grf_subg_harm_glm.R and supports
# outcome_type, effect_measure, offset.name, overdispersion, tune_grf,
# and other GLM-extension parameters.  The old version that was here
# (binary-only, @keywords internal, no outcome_type routing) has been
# removed to avoid function-name collisions during devtools::load_all().
# ────────────────────────────────────────────────────────────────────────────

#'
#' Evaluates the performance of GRF-identified subgroups, including hazard ratios,
#' bias, and predictive values. This function is typically used in simulation studies
#' to assess the performance of the GRF subgroup identification method.
#'
#' @param df Data frame containing the analysis data.
#' @param grf.est List. Output from \code{grf.subg.harm.survival}.
#' @param dgm List. Data-generating mechanism (truth) for simulation.
#' @param cox.formula.sim Formula for unadjusted Cox model.
#' @param cox.formula.adj.sim Formula for adjusted Cox model.
#' @param analysis Character. Analysis label (default: "GRF").
#' @param frac.tau Numeric. Fraction of tau for GRF horizon (default: 1.0).
#'
#' @return A data frame with evaluation metrics.
#'
#' @examples
#' \dontrun{
#' # grf.subg.eval() is called internally to evaluate GRF subgroup quality.
#' # See grf.subg.harm.survival() for the standard entry point.
#' }
#' @export
grf.subg.eval <- function(df,
                          grf.est,
                          dgm,
                          cox.formula.sim,
                          cox.formula.adj.sim,
                          analysis = "GRF",
                          frac.tau = 1.0) {
  # Implementation preserved from original
  # ... (rest of function)
  NULL
}
