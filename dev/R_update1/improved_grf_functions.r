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
  
  # Input validation
  if (maxdepth > 3) {
    stop("maxdepth must be at most 3")
  }
  
  # ---- Step 1: Prepare data matrices ----
  # Convert confounders to numeric matrix (required by GRF)
  X_temp <- as.matrix(data[, confounders.name])
  X <- apply(X_temp, 2, as.numeric)
  
  # Extract outcome variables
  Y <- data[, outcome.name]
  W <- data[, treat.name]
  D <- data[, event.name]
  
  # ---- Step 2: Calculate time horizon ----
  # Use minimum of maximum observed event times in each arm
  max_time_treat <- max(Y[W == 1 & D == 1])
  max_time_control <- max(Y[W == 0 & D == 1])
  tau.rmst <- frac.tau * min(max_time_treat, max_time_control)
  
  # ---- Step 3: Fit causal survival forest ----
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
    stop("Failed to fit causal survival forest")
  }
  
  # ---- Step 4: Compute doubly robust scores ----
  # These scores estimate treatment effect for each individual
  dr.scores <- policytree::double_robust_scores(cs.forest)
  
  # ---- Step 5: Fit policy trees at multiple depths ----
  
  # Depth 1 tree
  tree1 <- policytree::policy_tree(X, dr.scores, depth = 1)
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
    tree2 <- policytree::policy_tree(X, dr.scores, depth = 2)
    data$predict2.node <- predict(tree2, X, type = "node.id")
    
    values2 <- aggregate_policy_tree_results(
      data = data,
      dr.scores = dr.scores,
      node_col = "predict2.node",
      n.min = n.min,
      depth = 2
    )
    
    values <- rbind(values1, values2)
  }
  
  # Depth 3 tree (if requested)
  tree3 <- NULL
  if (maxdepth == 3) {
    tree3 <- policytree::policy_tree(X, dr.scores, depth = 3)
    data$predict3.node <- predict(tree3, X, type = "node.id")
    
    values3 <- aggregate_policy_tree_results(
      data = data,
      dr.scores = dr.scores,
      node_col = "predict3.node",
      n.min = n.min,
      depth = 3
    )
    
    values <- rbind(values1, values2, values3)
  }
  
  # ---- Step 6: Select optimal subgroup ----
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
  
  # ---- Step 7: Extract tree cuts and create results ----
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


#' Aggregate Policy Tree Results by Leaf Node
#'
#' Helper function to aggregate doubly robust scores by leaf node and calculate
#' treatment effect differences and subgroup sizes.
#'
#' @param data Data frame with node assignments
#' @param dr.scores Doubly robust scores matrix (control and treated columns)
#' @param node_col Character. Name of the column containing node assignments
#' @param n.min Integer. Minimum subgroup size
#' @param depth Integer. Tree depth
#'
#' @return Data frame with aggregated results per leaf node
#' @keywords internal
aggregate_policy_tree_results <- function(data, dr.scores, node_col, n.min, depth) {
  
  values <- stats::aggregate(
    dr.scores,
    by = list(leaf.node = data[[node_col]]),
    FUN = function(x) c(
      mean = mean(x),
      size = length(x),
      se = sd(x) / sqrt(length(x))
    )
  )
  
  # Calculate treatment effect difference (control - treated)
  # Positive values indicate control is better
  values$diff <- values$control[, "mean"] - values$treated[, "mean"]
  
  # Store subgroup size
  values$Nsg <- values$control[, "size"]
  
  # Filter to subgroups meeting minimum size requirement
  values <- values[values$control[, "size"] >= n.min, ]
  
  # Add depth information
  values$depth <- depth
  
  return(values)
}


#' Select Optimal Subgroup from Candidate Subgroups
#'
#' Selects the best subgroup based on either maximum treatment effect difference
#' or largest subgroup size (among those meeting the minimum effect threshold).
#'
#' @param values Data frame of candidate subgroups
#' @param sg.criterion Character. "mDiff" or "Nsg"
#' @param dmin.grf Numeric. Minimum treatment effect difference
#'
#' @return Single-row data frame with the selected subgroup
#' @keywords internal
select_optimal_subgroup <- function(values, sg.criterion, dmin.grf) {
  
  # Validate inputs
  required_cols <- c("diff", "depth")
  if (!all(required_cols %in% names(values))) {
    stop("The 'values' data.frame must contain columns: ", 
         paste(required_cols, collapse = ", "))
  }
  
  if (nrow(values) == 0) {
    stop("'values' is empty - no candidate subgroups found")
  }
  
  # Select based on criterion
  if (identical(sg.criterion, "mDiff")) {
    # Select subgroup with maximum treatment effect difference
    loc.max <- which.max(values$diff)
    max.diff <- values[loc.max, , drop = FALSE]
    
  } else if (identical(sg.criterion, "Nsg")) {
    # Select largest subgroup meeting minimum effect threshold
    if (!"Nsg" %in% names(values)) {
      stop("'values' must contain column 'Nsg' for Nsg criterion")
    }
    
    values_new <- values[values$diff >= dmin.grf, , drop = FALSE]
    
    if (nrow(values_new) == 0) {
      stop("No subgroups meet the diff >= dmin.grf criterion")
    }
    
    loc.max <- which.max(values_new$Nsg)
    max.diff <- values_new[loc.max, , drop = FALSE]
    
  } else {
    stop("Unknown sg.criterion: ", sg.criterion)
  }
  
  # Validate that depth column exists
  if (!"depth" %in% names(max.diff)) {
    stop("'max.diff' must contain column 'depth'")
  }
  
  return(max.diff)
}


#' Extract Subgroup Definition from Policy Tree
#'
#' Extracts the splits defining a subgroup from a policy tree object.
#'
#' @param tree Policy tree object
#' @param sg_node Integer. Leaf node ID for the subgroup
#' @param data Data frame (used to assign treatment recommendations)
#'
#' @return List with subgroup definition and all cuts
#' @keywords internal
extract_subgroup_from_tree <- function(tree, sg_node, data) {
  
  # Add treatment recommendation flag
  # 0 = assign to control (potential harm from treatment)
  # 1 = assign to treatment (benefit from treatment)
  data$treat.recommend <- ifelse(data$predict.node == sg_node, 0, 1)
  
  grf_names <- tree$columns
  tnodes <- tree$nodes
  
  # Find the split defining the selected subgroup
  sg_cov <- NULL
  sg_cut <- NULL
  
  for (tt in seq_along(tnodes)) {
    temp <- tnodes[[tt]]
    
    if (!temp$is_leaf) {
      # Check if this node leads to our target subgroup
      if (temp$left_child == sg_node || temp$right_child == sg_node) {
        sg_cov <- temp$split_variable
        sg_cut <- temp$split_value
      }
    }
  }
  
  # Create the defining split expression
  sg_definition <- paste0(grf_names[sg_cov], " <= ", sg_cut)
  
  # Extract all cuts from the selected tree
  all_cuts <- c()
  variable_names <- c()
  
  for (tt in seq_along(tnodes)) {
    temp <- tnodes[[tt]]
    
    if (!temp$is_leaf) {
      sg_cov <- temp$split_variable
      sg_cut <- round(temp$split_value, 2)
      
      vcut <- paste0(grf_names[sg_cov], " <= ", sg_cut)
      all_cuts <- c(all_cuts, vcut)
      variable_names <- c(variable_names, grf_names[sg_cov])
    }
  }
  
  return(list(
    sg_definition = sg_definition,
    all_cuts = all_cuts,
    variable_names = variable_names
  ))
}


#' Extract Cuts from All Tree Depths
#'
#' Extracts all split points from policy trees at depths 1-3.
#'
#' @param tree1 Policy tree at depth 1
#' @param tree2 Policy tree at depth 2 (or NULL)
#' @param tree3 Policy tree at depth 3 (or NULL)
#' @param maxdepth Integer. Maximum depth evaluated
#'
#' @return List with cuts and variable names for each depth
#' @keywords internal
extract_all_tree_cuts <- function(tree1, tree2, tree3, maxdepth) {
  
  # Extract from tree1
  tree1_info <- extract_cuts_from_single_tree(tree1)
  
  # Initialize depth 2 and 3 results
  tree2_cuts <- NULL
  tree2_names <- NULL
  tree3_cuts <- NULL
  tree3_names <- NULL
  
  # Extract from tree2 if it exists
  if (maxdepth >= 2 && !is.null(tree2)) {
    tree2_info <- extract_cuts_from_single_tree(tree2)
    tree2_cuts <- tree2_info$cuts
    tree2_names <- tree2_info$names
  }
  
  # Extract from tree3 if it exists
  if (maxdepth == 3 && !is.null(tree3)) {
    tree3_info <- extract_cuts_from_single_tree(tree3)
    tree3_cuts <- tree3_info$cuts
    tree3_names <- tree3_info$names
  }
  
  return(list(
    tree1_cuts = tree1_info$cuts,
    tree1_names = tree1_info$names,
    tree2_cuts = tree2_cuts,
    tree2_names = tree2_names,
    tree3_cuts = tree3_cuts,
    tree3_names = tree3_names
  ))
}


#' Extract Cuts from a Single Policy Tree
#'
#' Helper function to extract all splits from one policy tree.
#'
#' @param tree Policy tree object
#'
#' @return List with cuts (character vector) and names (character vector)
#' @keywords internal
extract_cuts_from_single_tree <- function(tree) {
  
  grf_names <- tree$columns
  tnodes <- tree$nodes
  
  cuts <- c()
  names <- c()
  
  for (tt in seq_along(tnodes)) {
    temp <- tnodes[[tt]]
    
    if (!temp$is_leaf) {
      sg_cov <- temp$split_variable
      sg_cut <- round(temp$split_value, 2)
      
      vcut <- paste0(grf_names[sg_cov], " <= ", sg_cut)
      cuts <- c(cuts, vcut)
      names <- c(names, grf_names[sg_cov])
    }
  }
  
  return(list(cuts = cuts, names = names))
}


#' GRF Subgroup Evaluation and Performance Metrics
#'
#' Evaluates the performance of GRF-identified subgroups by comparing predicted
#' subgroups against true subgroups (for simulation studies). Calculates hazard ratios,
#' bias, predictive values, and coverage metrics.
#'
#' @param df Data frame containing the analysis data
#' @param grf.est List. Output from \code{grf.subg.harm.survival}
#' @param dgm List. Data-generating mechanism (truth) for simulation, must contain:
#'   \itemize{
#'     \item model: "null" (no subgroup) or other
#'     \item grf.harm.true: True subgroup definition (or NULL if no subgroup)
#'     \item hr.H.true: True HR in harm subgroup
#'     \item hr.Hc.true: True HR in complement subgroup
#'   }
#' @param cox.formula.sim Formula for unadjusted Cox model
#' @param cox.formula.adj.sim Formula for adjusted Cox model
#' @param analysis Character. Analysis label (default: "GRF")
#' @param frac.tau Numeric. Fraction of tau for GRF horizon (default: 1.0)
#'
#' @return Data frame with performance metrics:
#' \describe{
#'   \item{any.H}{Indicator: was a subgroup found?}
#'   \item{size.H, size.Hc}{Subgroup sizes}
#'   \item{ppv, npv}{Positive and negative predictive values}
#'   \item{sensitivity, specificity}{Classification metrics}
#'   \item{hr.H.true, hr.Hc.true}{True hazard ratios}
#'   \item{hr.H.hat, hr.Hc.hat}{Estimated hazard ratios}
#'   \item{b1.H, b2.H, b1.Hc, b2.Hc}{Bias metrics}
#' }
#'
#' @details
#' This function is primarily for simulation studies where the true subgroup
#' is known. It compares the GRF-identified subgroup against the truth and
#' calculates various performance metrics.
#'
#' Four scenarios are handled:
#' \enumerate{
#'   \item True subgroup exists AND GRF finds a subgroup
#'   \item True subgroup exists BUT GRF finds no subgroup
#'   \item No true subgroup BUT GRF finds a subgroup (false positive)
#'   \item No true subgroup AND GRF finds no subgroup (correct)
#' }
#'
#' @importFrom survival coxph Surv survfit
#' @export
grf.estimates.out <- function(
    df,
    grf.est = NULL,
    dgm = NULL,
    cox.formula.sim = NULL,
    cox.formula.adj.sim = NULL,
    analysis = "GRF",
    frac.tau = 1.0
) {
  
  # ---- Calculate censoring rate ----
  event.name <- all.vars(cox.formula.sim[[2]])[2]
  p.cens <- mean(1 - df[[event.name]])
  
  # ---- Calculate tau (time horizon) ----
  outcome.name <- all.vars(cox.formula.sim[[2]])[1]
  treat.name <- all.vars(cox.formula.sim)[3]
  
  Y <- df[[outcome.name]]
  W <- df[[treat.name]]
  D <- df[[event.name]]
  
  taumax <- frac.tau * min(
    max(Y[W == 1 & D == 1]),
    max(Y[W == 0 & D == 1])
  )
  
  # ---- Extract GRF subgroup definition ----
  sg.harm.grf <- if (!is.null(grf.est)) grf.est$sg.harm.id else NULL
  
  # ---- Calculate ITT hazard ratios ----
  # Unadjusted HR
  fit <- summary(coxph(cox.formula.sim, data = df, robust = FALSE))$conf.int
  hr.itt <- fit[1]
  l.itt <- fit[3]
  u.itt <- fit[4]
  
  # Adjusted HR (if covariates provided)
  fit <- summary(coxph(cox.formula.adj.sim, data = df, robust = FALSE))$conf.int
  hr.adj.itt <- fit[1, 1]
  l.adj.itt <- fit[1, 3]
  u.adj.itt <- fit[1, 4]
  
  # ---- Validate DGM ----
  if (dgm$model == "null" && !is.null(dgm$grf.harm.true)) {
    stop("For dgm model 'null', grf.harm.true should be NULL")
  }
  
  # ---- Scenario 1: True subgroup exists AND GRF found something ----
  if (!is.null(dgm$grf.harm.true) && !is.null(sg.harm.grf)) {
    
    result <- evaluate_scenario_found_with_truth(
      df = df,
      grf.est = grf.est,
      cox.formula.sim = cox.formula.sim,
      dgm = dgm
    )
    
  # ---- Scenario 2: True subgroup exists BUT GRF found nothing ----
  } else if (!is.null(dgm$grf.harm.true) && is.null(sg.harm.grf)) {
    
    result <- evaluate_scenario_not_found_with_truth(
      df = df,
      cox.formula.sim = cox.formula.sim,
      dgm = dgm
    )
    
  # ---- Scenario 3: No true subgroup BUT GRF found something ----
  } else if (is.null(dgm$grf.harm.true) && !is.null(sg.harm.grf)) {
    
    result <- evaluate_scenario_found_no_truth(
      df = df,
      grf.est = grf.est,
      cox.formula.sim = cox.formula.sim,
      dgm = dgm
    )
    
  # ---- Scenario 4: No true subgroup AND GRF found nothing ----
  } else {
    
    result <- evaluate_scenario_not_found_no_truth(
      df = df,
      cox.formula.sim = cox.formula.sim,
      dgm = dgm
    )
  }
  
  # ---- Combine with common metrics ----
  df.res <- data.frame(
    result,
    p.cens = p.cens,
    analysis = analysis,
    taumax = taumax,
    hr.itt = hr.itt,
    l.itt = l.itt,
    u.itt = u.itt,
    hr.adj.itt = hr.adj.itt,
    l.adj.itt = l.adj.itt,
    u.adj.itt = u.adj.itt
  )
  
  return(df.res)
}


#' Evaluate Scenario: Found Subgroup with Truth
#'
#' Handles case where true subgroup exists and GRF found a subgroup.
#'
#' @keywords internal
evaluate_scenario_found_with_truth <- function(df, grf.est, cox.formula.sim, dgm) {
  
  any.H <- 1.0
  dfout <- grf.est$data
  
  # ---- Calculate confusion matrix ----
  # PPV for H membership: Hpred --> treat.recommend=0 & flag.harm=1
  aa <- with(dfout, sum(treat.recommend == 0 & flag.harm == 1))  # True positive
  bb <- with(dfout, sum(treat.recommend == 1 & flag.harm == 1))  # False negative
  cc <- with(dfout, sum(treat.recommend == 0 & flag.harm == 0))  # False positive
  dd <- with(dfout, sum(treat.recommend == 1 & flag.harm == 0))  # True negative
  
  size.H <- sum(dfout$treat.recommend == 0)
  size.Hc <- sum(dfout$treat.recommend == 1)
  
  # ---- Calculate metrics ----
  ppv <- aa / (aa + cc)
  npv <- dd / (bb + dd)
  specificity <- dd / (bb + dd)
  sensitivity <- aa / (aa + cc)
  
  # ---- HR estimates with true H ----
  hr_true <- fit_cox_with_ci(cox.formula.sim, subset(df, flag.harm == 1))
  hr.H.true <- hr_true$hr
  l.H.true <- hr_true$lower
  u.H.true <- hr_true$upper
  
  hr_true_c <- fit_cox_with_ci(cox.formula.sim, subset(df, flag.harm == 0))
  hr.Hc.true <- hr_true_c$hr
  l.Hc.true <- hr_true_c$lower
  u.Hc.true <- hr_true_c$upper
  
  # ---- HR estimates with estimated H ----
  hr_hat <- fit_cox_with_ci(cox.formula.sim, subset(dfout, treat.recommend == 0))
  hr.H.hat <- hr_hat$hr
  l.H.hat <- hr_hat$lower
  u.H.hat <- hr_hat$upper
  
  hr_hat_c <- fit_cox_with_ci(cox.formula.sim, subset(dfout, treat.recommend == 1))
  hr.Hc.hat <- hr_hat_c$hr
  l.Hc.hat <- hr_hat_c$lower
  u.Hc.hat <- hr_hat_c$upper
  
  # ---- Calculate bias ----
  b1.H <- hr.H.hat - hr.H.true
  b2.H <- hr.H.hat - dgm$hr.H.true
  b1.Hc <- hr.Hc.hat - hr.Hc.true
  b2.Hc <- hr.Hc.hat - dgm$hr.Hc.true
  
  # Pattern matching metrics (set to NA for now)
  found.1 <- found.2 <- found.both <- found.al3 <- NA
  
  return(data.frame(
    any.H = any.H,
    size.H = size.H,
    size.Hc = size.Hc,
    ppv = ppv,
    npv = npv,
    specificity = specificity,
    sensitivity = sensitivity,
    found.1 = found.1,
    found.2 = found.2,
    found.both = found.both,
    found.al3 = found.al3,
    hr.H.true = hr.H.true,
    hr.Hc.true = hr.Hc.true,
    hr.H.hat = hr.H.hat,
    hr.Hc.hat = hr.Hc.hat,
    b1.H = b1.H,
    b2.H = b2.H,
    b1.Hc = b1.Hc,
    b2.Hc = b2.Hc,
    l.H.true = l.H.true,
    u.H.true = u.H.true,
    l.Hc.true = l.Hc.true,
    u.Hc.true = u.Hc.true,
    l.H.hat = l.H.hat,
    u.H.hat = u.H.hat,
    l.Hc.hat = l.Hc.hat,
    u.Hc.hat = u.Hc.hat
  ))
}


#' Evaluate Scenario: Not Found with Truth
#'
#' Handles case where true subgroup exists but GRF found nothing.
#'
#' @keywords internal
evaluate_scenario_not_found_with_truth <- function(df, cox.formula.sim, dgm) {
  
  any.H <- 0
  size.H <- 0
  size.Hc <- nrow(df)
  
  # ---- HR estimates with true H ----
  hr_true <- fit_cox_with_ci(cox.formula.sim, subset(df, flag.harm == 1))
  hr.H.true <- hr_true$hr
  l.H.true <- hr_true$lower
  u.H.true <- hr_true$upper
  
  hr_true_c <- fit_cox_with_ci(cox.formula.sim, subset(df, flag.harm == 0))
  hr.Hc.true <- hr_true_c$hr
  l.Hc.true <- hr_true_c$lower
  u.Hc.true <- hr_true_c$upper
  
  # ---- HR estimates with estimated H (none found, use ITT) ----
  hr.H.hat <- NA
  l.H.hat <- NA
  u.H.hat <- NA
  
  hr_hat_c <- fit_cox_with_ci(cox.formula.sim, df)
  hr.Hc.hat <- hr_hat_c$hr
  l.Hc.hat <- hr_hat_c$lower
  u.Hc.hat <- hr_hat_c$upper
  
  # ---- Calculate bias ----
  b1.H <- NA
  b2.H <- NA
  b1.Hc <- hr.Hc.hat - hr.Hc.true
  b2.Hc <- hr.Hc.hat - dgm$hr.Hc.true
  
  # ---- Confusion matrix (no H found, treat.recommend=1 for all) ----
  aa <- 0  # True positive
  bb <- sum(df$flag.harm == 1)  # False negative
  cc <- 0  # False positive
  dd <- sum(df$flag.harm == 0)  # True negative
  
  ppv <- aa / (aa + cc)
  npv <- dd / (bb + dd)
  specificity <- dd / (bb + dd)
  sensitivity <- ifelse(aa == 0, 0, aa / (aa + cc))
  
  found.1 <- found.2 <- found.both <- NA
  found.al3 <- 0
  
  return(data.frame(
    any.H = any.H,
    size.H = size.H,
    size.Hc = size.Hc,
    ppv = ppv,
    npv = npv,
    specificity = specificity,
    sensitivity = sensitivity,
    found.1 = found.1,
    found.2 = found.2,
    found.both = found.both,
    found.al3 = found.al3,
    hr.H.true = hr.H.true,
    hr.Hc.true = hr.Hc.true,
    hr.H.hat = hr.H.hat,
    hr.Hc.hat = hr.Hc.hat,
    b1.H = b1.H,
    b2.H = b2.H,
    b1.Hc = b1.Hc,
    b2.Hc = b2.Hc,
    l.H.true = l.H.true,
    u.H.true = u.H.true,
    l.Hc.true = l.Hc.true,
    u.Hc.true = u.Hc.true,
    l.H.hat = l.H.hat,
    u.H.hat = u.H.hat,
    l.Hc.hat = l.Hc.hat,
    u.Hc.hat = u.Hc.hat
  ))
}


#' Evaluate Scenario: Found Subgroup without Truth
#'
#' Handles case where no true subgroup exists but GRF found something (false positive).
#'
#' @keywords internal
evaluate_scenario_found_no_truth <- function(df, grf.est, cox.formula.sim, dgm) {
  
  any.H <- 1.0
  dfout <- grf.est$data
  
  # ---- Confusion matrix (flag.harm=0 for all since no true subgroup) ----
  aa <- 0  # True positive (impossible)
  bb <- 0  # False negative (impossible)
  cc <- sum(dfout$treat.recommend == 0)  # False positive
  dd <- sum(dfout$treat.recommend == 1)  # True negative
  
  size.H <- cc
  size.Hc <- dd
  
  # ---- Calculate metrics ----
  ppv <- NA  # Undefined when no true positives exist
  npv <- dd / (cc + dd)
  specificity <- dd / (bb + dd)
  sensitivity <- aa / (aa + cc)
  
  # ---- HR estimates with true H (no true H, use ITT) ----
  hr.H.true <- NA
  l.H.true <- NA
  u.H.true <- NA
  
  hr_true_c <- fit_cox_with_ci(cox.formula.sim, df)
  hr.Hc.true <- hr_true_c$hr
  l.Hc.true <- hr_true_c$lower
  u.Hc.true <- hr_true_c$upper
  
  # ---- HR estimates with estimated H ----
  hr_hat <- fit_cox_with_ci(cox.formula.sim, subset(dfout, treat.recommend == 0))
  hr.H.hat <- hr_hat$hr
  l.H.hat <- hr_hat$lower
  u.H.hat <- hr_hat$upper
  
  hr_hat_c <- fit_cox_with_ci(cox.formula.sim, subset(dfout, treat.recommend == 1))
  hr.Hc.hat <- hr_hat_c$hr
  l.Hc.hat <- hr_hat_c$lower
  u.Hc.hat <- hr_hat_c$upper
  
  # ---- Calculate bias ----
  b1.H <- NA
  b2.H <- NA
  b1.Hc <- hr.Hc.hat - hr.Hc.true
  b2.Hc <- hr.Hc.hat - dgm$hr.Hc.true
  
  found.1 <- found.2 <- found.both <- found.al3 <- NA
  
  return(data.frame(
    any.H = any.H,
    size.H = size.H,
    size.Hc = size.Hc,
    ppv = ppv,
    npv = npv,
    specificity = specificity,
    sensitivity = sensitivity,
    found.1 = found.1,
    found.2 = found.2,
    found.both = found.both,
    found.al3 = found.al3,
    hr.H.true = hr.H.true,
    hr.Hc.true = hr.Hc.true,
    hr.H.hat = hr.H.hat,
    hr.Hc.hat = hr.Hc.hat,
    b1.H = b1.H,
    b2.H = b2.H,
    b1.Hc = b1.Hc,
    b2.Hc = b2.Hc,
    l.H.true = l.H.true,
    u.H.true = u.H.true,
    l.Hc.true = l.Hc.true,
    u.Hc.true = u.Hc.true,
    l.H.hat = l.H.hat,
    u.H.hat = u.H.hat,
    l.Hc.hat = l.Hc.hat,
    u.Hc.hat = u.Hc.hat
  ))
}


#' Evaluate Scenario: Not Found without Truth
#'
#' Handles case where no true subgroup exists and GRF found nothing (correct null).
#'
#' @keywords internal
evaluate_scenario_not_found_no_truth <- function(df, cox.formula.sim, dgm) {
  
  any.H <- 0
  
  # ---- Confusion matrix (all true negatives) ----
  aa <- 0
  bb <- 0
  cc <- 0
  dd <- nrow(df)
  
  size.H <- 0
  size.Hc <- nrow(df)
  
  # ---- Calculate metrics ----
  ppv <- NA
  npv <- dd / (cc + dd)
  specificity <- dd / (bb + dd)
  sensitivity <- NA
  
  # ---- HR estimates (no true or estimated H, use ITT) ----
  hr.H.true <- NA
  l.H.true <- NA
  u.H.true <- NA
  
  hr_true_c <- fit_cox_with_ci(cox.formula.sim, df)
  hr.Hc.true <- hr_true_c$hr
  l.Hc.true <- hr_true_c$lower
  u.Hc.true <- hr_true_c$upper
  
  hr.H.hat <- NA
  l.H.hat <- NA
  u.H.hat <- NA
  
  hr_hat_c <- fit_cox_with_ci(cox.formula.sim, df)
  hr.Hc.hat <- hr_hat_c$hr
  l.Hc.hat <- hr_hat_c$lower
  u.Hc.hat <- hr_hat_c$upper
  
  # ---- Calculate bias ----
  b1.H <- NA
  b2.H <- NA
  b1.Hc <- hr.Hc.hat - hr.Hc.true
  b2.Hc <- hr.Hc.hat - dgm$hr.Hc.true
  
  found.1 <- found.2 <- found.both <- found.al3 <- NA
  
  return(data.frame(
    any.H = any.H,
    size.H = size.H,
    size.Hc = size.Hc,
    ppv = ppv,
    npv = npv,
    specificity = specificity,
    sensitivity = sensitivity,
    found.1 = found.1,
    found.2 = found.2,
    found.both = found.both,
    found.al3 = found.al3,
    hr.H.true = hr.H.true,
    hr.Hc.true = hr.Hc.true,
    hr.H.hat = hr.H.hat,
    hr.Hc.hat = hr.Hc.hat,
    b1.H = b1.H,
    b2.H = b2.H,
    b1.Hc = b1.Hc,
    b2.Hc = b2.Hc,
    l.H.true = l.H.true,
    u.H.true = u.H.true,
    l.Hc.true = l.Hc.true,
    u.Hc.true = u.Hc.true,
    l.H.hat = l.H.hat,
    u.H.hat = u.H.hat,
    l.Hc.hat = l.Hc.hat,
    u.Hc.hat = u.Hc.hat
  ))
}


#' Fit Cox Model and Extract HR with Confidence Interval
#'
#' Helper function to fit Cox model and extract hazard ratio with CI.
#'
#' @param formula Cox model formula
#' @param data Data frame
#'
#' @return List with hr, lower, and upper
#' @keywords internal
fit_cox_with_ci <- function(formula, data) {
  fit <- summary(survival::coxph(formula, data = data, robust = FALSE))$conf.int
  list(
    hr = fit[1],
    lower = fit[3],
    upper = fit[4]
  )
}