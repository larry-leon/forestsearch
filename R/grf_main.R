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
#'
#' @return A list with GRF results, including:
#'   \item{data}{Original data with added treatment recommendation flags}
#'   \item{grf.gsub}{Selected subgroup information}
#'   \item{sg.harm.id}{Expression defining the identified subgroup}
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
                                   return_selected_cuts_only = FALSE) {

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

  # ===========================================================================
  # SECTION: DATA PREPARATION
  # Purpose: Convert data to appropriate format for GRF analysis
  # ===========================================================================

  # Convert confounders to numeric matrix
  temp_matrix <- as.matrix(data[, confounders.name])
  X <- apply(temp_matrix, 2, as.numeric)

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

  cs.forest <- fit_causal_forest(X, Y, W, D, tau.rmst, config$RCT, config$seedit)

  # ===========================================================================
  # SECTION: SUBGROUP IDENTIFICATION VIA POLICY TREES
  # Purpose: Use policy trees to partition the covariate space
  # ===========================================================================

  # Compute doubly robust scores for subgroup identification
  dr.scores <- policytree::double_robust_scores(cs.forest)

  # Maximum sample size (used to exclude full population as subgroup)
  n.max <- length(Y)

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

  # Find the specific split that defines the subgroup
  sg_harm_id <- find_leaf_split(selected_tree, best_subgroup$leaf.node)

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

  return(result)
}


#' GRF Subgroup Identification for GLM Outcomes
#'
#' Identifies subgroups with differential treatment effect using generalized
#' random forests (GRF) and policy trees for binary or continuous outcomes.
#' This is the GLM counterpart to [grf.subg.harm.survival()], using
#' [grf::causal_forest()] instead of [grf::causal_survival_forest()].
#'
#' @section Sign Convention and `adverse_outcome`:
#'
#' [grf::causal_forest()] estimates the conditional average treatment effect
#' (CATE) as `E[Y(1) - Y(0) | X]`, where higher Y is implicitly "better."
#' The downstream policy tree machinery selects subgroups where the control
#' action has higher doubly-robust scores (i.e., treatment is harmful).
#'
#' This works correctly when **higher Y = better** (e.g., survival time,
#' quality of life score).  However, in clinical trials the outcome is
#' typically an adverse event (Y = 1 is death, relapse, etc.), meaning
#' **higher Y = worse**.  In this case:
#'
#' \itemize{
#'   \item The harmful subgroup has `diff = control_score - treated_score < 0`
#'     because treatment increases the bad outcome
#'   \item The `diff >= dmin.grf` gate selects the **complement** instead
#' }
#'
#' Setting `adverse_outcome = TRUE` (the default for binary outcomes) flips
#' the outcome to `1 - Y` before fitting the causal forest, so that higher
#' values = no event = better.  This aligns the sign convention with
#' [grf::causal_survival_forest()], where higher RMST = better.
#'
#' The flip is internal to GRF candidate selection only.  It does **not**
#' affect the downstream search, consistency, or bootstrap pipeline, which
#' operate on the original outcome scale.
#'
#' @param data Data frame containing the analysis data.
#' @param confounders.name Character vector of confounder variable names.
#' @param outcome.name Character. Name of outcome variable.
#' @param treat.name Character. Name of treatment group variable (0/1).
#' @param id.name Character. Name of ID variable.
#' @param n.min Integer. Minimum subgroup size (default: 60).
#' @param dmin.grf Numeric. Minimum difference in subgroup mean DR scores
#'   for a leaf node to qualify as a candidate subgroup (default: 0.0).
#'   Units are on the probability scale for binary outcomes (e.g.,
#'   `dmin.grf = 0.05` requires a 5 percentage-point difference).
#' @param RCT Logical. Is the data from a randomized controlled trial?
#'   (default: TRUE)
#' @param details Logical. Print details during execution (default: FALSE).
#' @param sg.criterion Character. Subgroup selection criterion
#'   ("mDiff" or "Nsg").
#' @param maxdepth Integer. Maximum tree depth (1, 2, or 3; default: 2).
#' @param seedit Integer. Random seed (default: 8316951).
#' @param return_selected_cuts_only Logical. If TRUE, returns only cuts from
#'   the selected tree depth (default: FALSE).
#' @param adverse_outcome Logical. If TRUE (the default), the outcome Y = 1
#'   is treated as an adverse event.  The outcome is flipped to `1 - Y` before
#'   fitting the causal forest so that the harm-detection sign convention
#'   aligns with [grf::causal_survival_forest()].  Set to FALSE when
#'   higher Y = better (e.g., continuous quality-of-life scores).
#'
#' @return A list with GRF results (same structure as
#'   [grf.subg.harm.survival()]).
#'
#' @examples
#' \dontrun{
#' # Binary adverse event (default: adverse_outcome = TRUE)
#' result <- grf.subg.harm.glm(
#'   data = trial_data,
#'   confounders.name = c("age", "biomarker"),
#'   outcome.name = "response",
#'   treat.name = "treat",
#'   id.name = "id"
#' )
#'
#' # Continuous outcome where higher = better (e.g., QoL)
#' result_qol <- grf.subg.harm.glm(
#'   data = trial_data,
#'   confounders.name = c("age", "biomarker"),
#'   outcome.name = "qol_score",
#'   treat.name = "treat",
#'   id.name = "id",
#'   adverse_outcome = FALSE
#' )
#' }
#'
#' @keywords internal
grf.subg.harm.glm <- function(data,
                               confounders.name,
                               outcome.name,
                               treat.name,
                               id.name,
                               n.min = 60,
                               dmin.grf = 0.0,
                               RCT = TRUE,
                               details = FALSE,
                               sg.criterion = "mDiff",
                               maxdepth = 2,
                               seedit = 8316951,
                               return_selected_cuts_only = FALSE,
                               adverse_outcome = TRUE) {

  if (maxdepth > 3) stop("Maximum depth cannot exceed 3")

  valid_criteria <- c("mDiff", "Nsg")
  if (!sg.criterion %in% valid_criteria) {
    stop("sg.criterion must be one of: ", paste(valid_criteria, collapse = ", "))
  }

  config <- create_grf_config(
    frac.tau = 1.0,
    n.min = n.min,
    dmin.grf = dmin.grf,
    RCT = RCT,
    sg.criterion = sg.criterion,
    maxdepth = maxdepth,
    seedit = seedit
  )
  config$return_selected_cuts_only <- return_selected_cuts_only

  # Data preparation — no event variable needed
  temp_matrix <- as.matrix(data[, confounders.name])
  X <- apply(temp_matrix, 2, as.numeric)
  Y <- data[, outcome.name]
  W <- data[, treat.name]

  # -----------------------------------------------------------------------
  # Adverse outcome flip
  #
  # causal_forest estimates CATE = E[Y(1) - Y(0) | X] and the policy tree
  # selects leaves where control has higher DR scores (diff > 0).  When
  # Y = adverse event (higher = worse), the harmful subgroup has diff < 0,
  # so the gate diff >= dmin.grf selects the wrong group.
  #
  # Flipping to 1 - Y makes higher = no event = better, aligning with
  # causal_survival_forest's convention (higher RMST = better).
  # -----------------------------------------------------------------------
  if (adverse_outcome) {
    if (details) {
      cat("GRF GLM: flipping outcome (adverse_outcome = TRUE)\n")
      cat("  Original Y range: [", min(Y), ",", max(Y), "]\n")
    }
    Y <- 1L - Y
    if (details) {
      cat("  Flipped  Y range: [", min(Y), ",", max(Y), "]\n")
    }
  }

  # Validate minimum arm sizes
  n_treat <- sum(W == 1)
  n_ctrl  <- sum(W == 0)
  if (n_treat < config$n.min || n_ctrl < config$n.min) {
    return(create_null_result(data, NULL, list(), config))
  }

  # Fit causal forest (GLM — no horizon/tau needed)
  cs.forest <- fit_causal_forest_glm(X, Y, W, config$RCT, config$seedit)

  # Everything downstream is identical to the survival path
  dr.scores <- policytree::double_robust_scores(cs.forest)
  n.max <- length(Y)

  tree_results <- fit_policy_trees(X, data, dr.scores, config$maxdepth, config$n.min)
  trees <- tree_results$trees
  values <- tree_results$values

  best_subgroup <- select_best_subgroup(
    values = values,
    sg.criterion = config$sg.criterion,
    dmin.grf = config$dmin.grf,
    n.max = n.max
  )

  if (is.null(best_subgroup)) {
    if (details) print_grf_details(config, values, NULL, NULL)
    return(create_null_result(data, values, trees, config))
  }

  data <- assign_subgroup_membership(data, best_subgroup, trees, X)
  selected_tree <- trees[[best_subgroup$depth]]
  sg_harm_id <- find_leaf_split(selected_tree, best_subgroup$leaf.node)

  if (config$return_selected_cuts_only) {
    tree_cuts <- extract_selected_tree_cuts(trees, best_subgroup$depth, config$maxdepth)
  } else {
    tree_cuts <- extract_all_tree_cuts(trees, config$maxdepth)
  }

  if (details) print_grf_details(config, values, best_subgroup, sg_harm_id, tree_cuts)

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

  # Attach the fitted causal forest for downstream PS extraction
  result$forest <- cs.forest

  return(result)
}
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
