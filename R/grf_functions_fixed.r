#' GRF Subgroup Identification for Survival Data
#'
#' Identifies subgroups with differential treatment effect using generalized random forests (GRF)
#' and policy trees. The function fits causal survival forests at multiple tree depths and selects
#' the optimal subgroup based on either maximum difference in treatment effect or largest subgroup size.
#'
#' @param data Data frame containing the analysis data.
#' @param confounders.name Character vector of confounder variable names.
#' @param outcome.name Character. Name of outcome variable (e.g., time-to-event).
#' @param event.name Character. Name of event indicator variable (0/1).
#' @param id.name Character. Name of ID variable.
#' @param treat.name Character. Name of treatment group variable (0/1).
#' @param frac.tau Numeric. Fraction of maximum follow-up time to use as horizon (0-1). Default: 1.0.
#' @param n.min Integer. Minimum subgroup size required. Default: 60.
#' @param dmin.grf Numeric. Minimum treatment effect difference (on RMST scale) required. Default: 0.0.
#' @param RCT Logical. Is the data from a randomized controlled trial? If TRUE, uses W.hat=0.5. Default: TRUE.
#' @param details Logical. Print diagnostic information during execution? Default: FALSE.
#' @param sg.criterion Character. Subgroup selection criterion:
#'   \itemize{
#'     \item "mDiff": Select subgroup with maximum treatment effect difference (default)
#'     \item "Nsg": Select largest subgroup meeting dmin.grf threshold
#'   }
#' @param maxdepth Integer. Maximum tree depth to evaluate (1, 2, or 3). Default: 2.
#' @param seedit Integer. Random seed for reproducibility. Default: 8316951.
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{Original data with added treatment recommendation and node assignments}
#'   \item{grf.gsub}{Optimal subgroup information (effect difference, size, etc.)}
#'   \item{sg.harm.id}{Subgroup definition as a string (e.g., "age <= 65")}
#'   \item{tree.cuts}{All cut points from the selected tree}
#'   \item{tree.names}{Variable names used in the selected tree}
#'   \item{tree1.cuts, tree2.cuts, tree3.cuts}{Cut points from each depth}
#'   \item{tree1.names, tree2.names, tree3.names}{Variable names from each depth}
#'   \item{harm.any}{All candidate subgroups with positive treatment effects}
#'   \item{tree}{Selected policy tree object}
#'   \item{tau.rmst}{Restricted mean survival time horizon used}
#'   \item{tree1, tree2, tree3}{Policy tree objects for each depth}
#' }
#'
#' @details
#' The function works in three main steps:
#' \enumerate{
#'   \item Fits a causal survival forest to estimate treatment effects
#'   \item Fits policy trees at depths 1-3 to identify subgroups
#'   \item Evaluates subgroups based on effect size and sample size criteria
#' }
#'
#' A subgroup is only returned if:
#' \itemize{
#'   \item Treatment effect difference >= dmin.grf
#'   \item Subgroup size >= n.min
#'   \item Subgroup size < total sample size (not the entire population)
#' }
#'
#' @importFrom grf causal_survival_forest
#' @importFrom policytree double_robust_scores policy_tree
#' @importFrom stats aggregate
#' @export
#'
#' @examples
#' \dontrun{
#' # Identify subgroups in RCT data
#' result <- grf.subg.harm.survival(
#'   data = trial_data,
#'   confounders.name = c("age", "sex", "baseline_score"),
#'   outcome.name = "time",
#'   event.name = "event",
#'   id.name = "patient_id",
#'   treat.name = "treatment",
#'   n.min = 50,
#'   dmin.grf = 0.1
#' )
#' }
grf.subg.harm.survival <- function(
    data,
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
    seedit = 8316951
) {
  
  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================
  
  # Validate data frame
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (nrow(data) < 10) {
    stop("'data' must have at least 10 rows")
  }
  
  # Validate character inputs
  if (!is.character(confounders.name) || length(confounders.name) < 1) {
    stop("'confounders.name' must be a character vector with at least 1 element")
  }
  if (!is.character(outcome.name) || length(outcome.name) != 1) {
    stop("'outcome.name' must be a single character value")
  }
  if (!is.character(event.name) || length(event.name) != 1) {
    stop("'event.name' must be a single character value")
  }
  if (!is.character(id.name) || length(id.name) != 1) {
    stop("'id.name' must be a single character value")
  }
  if (!is.character(treat.name) || length(treat.name) != 1) {
    stop("'treat.name' must be a single character value")
  }
  
  # Validate columns exist
  required_cols <- c(confounders.name, outcome.name, event.name, id.name, treat.name)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in 'data': ", paste(missing_cols, collapse = ", "))
  }
  
  # Validate numeric parameters
  if (!is.numeric(frac.tau) || length(frac.tau) != 1 || 
      frac.tau <= 0 || frac.tau > 1) {
    stop("'frac.tau' must be a single numeric value between 0 and 1")
  }
  if (!is.numeric(n.min) || length(n.min) != 1 || n.min < 1) {
    stop("'n.min' must be a positive integer")
  }
  if (!is.numeric(dmin.grf) || length(dmin.grf) != 1) {
    stop("'dmin.grf' must be a single numeric value")
  }
  if (!is.numeric(seedit) || length(seedit) != 1) {
    stop("'seedit' must be a single numeric value")
  }
  
  # Validate logical parameters
  if (!is.logical(RCT) || length(RCT) != 1 || is.na(RCT)) {
    stop("'RCT' must be a single logical value (TRUE or FALSE)")
  }
  if (!is.logical(details) || length(details) != 1 || is.na(details)) {
    stop("'details' must be a single logical value (TRUE or FALSE)")
  }
  
  # Validate sg.criterion
  if (!sg.criterion %in% c("mDiff", "Nsg")) {
    stop("'sg.criterion' must be either 'mDiff' or 'Nsg'")
  }
  
  # Validate maxdepth
  if (!maxdepth %in% c(1, 2, 3)) {
    stop("'maxdepth' must be 1, 2, or 3")
  }
  
  # ============================================================================
  # STEP 1: Prepare data matrices
  # ============================================================================
  
  # Convert confounders to numeric matrix (required by GRF)
  X_temp <- as.matrix(data[, confounders.name, drop = FALSE])
  X <- apply(X_temp, 2, as.numeric)
  
  # Validate X conversion
  if (any(is.na(X))) {
    stop("NA values detected in confounder matrix after numeric conversion")
  }
  
  # Extract outcome variables
  Y <- data[, outcome.name]
  W <- data[, treat.name]
  D <- data[, event.name]
  
  # Validate outcome variables
  if (!is.numeric(Y) || any(is.na(Y))) {
    stop("'", outcome.name, "' must be numeric without NA values")
  }
  if (!is.numeric(W) || any(is.na(W)) || !all(W %in% c(0, 1))) {
    stop("'", treat.name, "' must be numeric 0/1 without NA values")
  }
  if (!is.numeric(D) || any(is.na(D)) || !all(D %in% c(0, 1))) {
    stop("'", event.name, "' must be numeric 0/1 without NA values")
  }
  
  # ============================================================================
  # STEP 2: Calculate time horizon
  # ============================================================================
  
  # Use minimum of maximum observed event times in each arm
  events_treat <- Y[W == 1 & D == 1]
  events_control <- Y[W == 0 & D == 1]
  
  if (length(events_treat) == 0 || length(events_control) == 0) {
    stop("Insufficient events in one or both treatment arms")
  }
  
  max_time_treat <- max(events_treat)
  max_time_control <- max(events_control)
  tau.rmst <- frac.tau * min(max_time_treat, max_time_control)
  
  if (tau.rmst <= 0) {
    stop("Calculated tau.rmst is <= 0, check your data")
  }
  
  # ============================================================================
  # STEP 3: Fit causal survival forest
  # ============================================================================
  
  if (RCT) {
    # For RCT: fix propensity score at 0.5
    cs.forest <- try(
      suppressWarnings(
        grf::causal_survival_forest(
          X, Y, W, D,
          W.hat = 0.5,
          horizon = tau.rmst,
          seed = seedit
        )
      ),
      silent = TRUE
    )
  } else {
    # For observational study: estimate propensity scores
    cs.forest <- try(
      suppressWarnings(
        grf::causal_survival_forest(
          X, Y, W, D,
          horizon = tau.rmst,
          seed = seedit
        )
      ),
      silent = TRUE
    )
  }
  
  if (inherits(cs.forest, "try-error")) {
    stop("Failed to fit causal survival forest: ", 
         attr(cs.forest, "condition")$message)
  }
  
  # ============================================================================
  # STEP 4: Compute doubly robust scores
  # ============================================================================
  
  # These scores estimate treatment effect for each individual
  dr.scores <- try(policytree::double_robust_scores(cs.forest), silent = TRUE)
  
  if (inherits(dr.scores, "try-error")) {
    stop("Failed to compute doubly robust scores: ", 
         attr(dr.scores, "condition")$message)
  }
  
  # ============================================================================
  # STEP 5: Fit policy trees at multiple depths
  # ============================================================================
  
  # Depth 1 tree
  tree1 <- try(policytree::policy_tree(X, dr.scores, depth = 1), silent = TRUE)
  if (inherits(tree1, "try-error")) {
    stop("Failed to fit policy tree at depth 1: ", 
         attr(tree1, "condition")$message)
  }
  
  data$predict1.node <- predict(tree1, X, type = "node.id")
  
  values1 <- aggregate_policy_tree_results(
    data = data,
    dr.scores = dr.scores,
    node_col = "predict1.node",
    n.min = n.min,
    depth = 1
  )
  
  values <- values1
  
  # Depth 2 tree (if requested)
  tree2 <- NULL
  if (maxdepth >= 2) {
    tree2 <- try(policytree::policy_tree(X, dr.scores, depth = 2), silent = TRUE)
    if (!inherits(tree2, "try-error")) {
      data$predict2.node <- predict(tree2, X, type = "node.id")
      
      values2 <- aggregate_policy_tree_results(
        data = data,
        dr.scores = dr.scores,
        node_col = "predict2.node",
        n.min = n.min,
        depth = 2
      )
      
      values <- rbind(values1, values2)
    } else if (details) {
      cat("Warning: Failed to fit policy tree at depth 2\n")
    }
  }
  
  # Depth 3 tree (if requested)
  tree3 <- NULL
  if (maxdepth == 3) {
    tree3 <- try(policytree::policy_tree(X, dr.scores, depth = 3), silent = TRUE)
    if (!inherits(tree3, "try-error")) {
      data$predict3.node <- predict(tree3, X, type = "node.id")
      
      values3 <- aggregate_policy_tree_results(
        data = data,
        dr.scores = dr.scores,
        node_col = "predict3.node",
        n.min = n.min,
        depth = 3
      )
      
      values <- rbind(values1, values2, values3)
    } else if (details) {
      cat("Warning: Failed to fit policy tree at depth 3\n")
    }
  }
  
  # ============================================================================
  # STEP 6: Select optimal subgroup
  # ============================================================================
  
  max.diff <- select_optimal_subgroup(
    values = values,
    sg.criterion = sg.criterion,
    dmin.grf = dmin.grf
  )
  
  # Assign the appropriate tree and node predictions based on selected depth
  if (max.diff$depth == 1) {
    data$predict.node <- data$predict1.node
    tree <- tree1
  } else if (maxdepth >= 2 && max.diff$depth == 2) {
    data$predict.node <- data$predict2.node
    tree <- tree2
  } else if (maxdepth >= 3 && max.diff$depth == 3) {
    data$predict.node <- data$predict3.node
    tree <- tree3
  } else {
    stop("No matching depth found for tree assignment")
  }
  
  # Print diagnostics if requested
  if (details) {
    cat("Tau (time horizon), maxdepth:", c(tau.rmst, maxdepth), "\n")
    temp <- values[, c("leaf.node", "control", "depth")]
    print(round(temp, 2))
    temp <- max.diff[, c("leaf.node", "control", "depth")]
    print(round(temp, 2))
  }
  
  # ============================================================================
  # STEP 7: Extract tree cuts and create results
  # ============================================================================
  
  n.max <- length(Y)  # Total sample size
  
  # Check if valid subgroup was found
  subgroup_is_valid <- max.diff[, "diff"] >= dmin.grf && max.diff[, "Nsg"] < n.max
  
  if (subgroup_is_valid) {
    # Extract subgroup information from selected tree
    sg_info <- extract_subgroup_from_tree(
      tree = tree,
      sg_node = max.diff$leaf.node,
      data = data
    )
    
    # Extract cuts from all tree depths
    tree_cuts_all <- extract_all_tree_cuts(
      tree1 = tree1,
      tree2 = tree2,
      tree3 = tree3,
      maxdepth = maxdepth
    )
    
    if (details) {
      cat("GRF subgroup found\n")
      cat("All splits:\n")
      print(sg_info$all_cuts)
      cat("Terminating node at max.diff (sg.harm.id):\n")
      print(sg_info$sg_definition)
    }
    
    result <- list(
      data = data,
      grf.gsub = max.diff,
      sg.harm.id = sg_info$sg_definition,
      tree.cuts = sg_info$all_cuts,
      tree.names = unique(sg_info$variable_names),
      tree1.cuts = tree_cuts_all$tree1_cuts,
      tree1.names = unique(tree_cuts_all$tree1_names),
      tree2.cuts = tree_cuts_all$tree2_cuts,
      tree2.names = unique(tree_cuts_all$tree2_names),
      tree3.cuts = tree_cuts_all$tree3_cuts,
      tree3.names = unique(tree_cuts_all$tree3_names),
      harm.any = values[values$diff > 0, ],
      tree = tree,
      tau.rmst = tau.rmst,
      tree1 = tree1,
      tree2 = tree2,
      tree3 = tree3
    )
    
  } else {
    # No valid subgroup found
    if (details) {
      cat("GRF subgroup NOT found\n")
    }
    
    harm.any <- NULL
    cand.sgs <- which(values$diff > 0)
    if (length(cand.sgs) > 0) {
      harm.any <- values[cand.sgs, ]
    }
    
    result <- list(
      data = data,
      grf.gsub = NULL,
      sg.harm.id = NULL,
      harm.any = harm.any,
      tree = tree,
      tau.rmst = tau.rmst,
      dmin.grf = dmin.grf,
      frac.tau = frac.tau,
      maxdepth = maxdepth,
      n.min = n.min
    )
  }
  
  return(result)
}
