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
  # Input validation
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!is.matrix(dr.scores)) stop("'dr.scores' must be a matrix")
  if (!is.character(node_col) || !node_col %in% names(data)) {
    stop("'node_col' must be a column name in 'data'")
  }
  if (!is.numeric(n.min) || n.min < 1) stop("'n.min' must be a positive integer")
  if (!is.numeric(depth) || !depth %in% c(1, 2, 3)) {
    stop("'depth' must be 1, 2, or 3")
  }
  
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
  # Input validation
  if (is.null(tree)) stop("'tree' cannot be NULL")
  if (!is.numeric(sg_node) || length(sg_node) != 1) {
    stop("'sg_node' must be a single numeric value")
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!"predict.node" %in% names(data)) {
    stop("'data' must contain column 'predict.node'")
  }
  
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
  
  # Validate that we found the split
  if (is.null(sg_cov) || is.null(sg_cut)) {
    stop("Could not find split defining subgroup node ", sg_node)
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
  # Input validation
  if (is.null(tree1)) stop("'tree1' cannot be NULL")
  if (!maxdepth %in% c(1, 2, 3)) stop("'maxdepth' must be 1, 2, or 3")
  
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
  # Input validation
  if (is.null(tree)) stop("'tree' cannot be NULL")
  if (!is.list(tree) || !"nodes" %in% names(tree)) {
    stop("'tree' must be a policy tree object with 'nodes' element")
  }
  
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
  # Input validation
  if (!inherits(formula, "formula")) stop("'formula' must be a formula")
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (nrow(data) < 2) stop("'data' must have at least 2 rows")
  
  # Fit model with error handling
  fit_result <- try(
    survival::coxph(formula, data = data, robust = FALSE),
    silent = TRUE
  )
  
  if (inherits(fit_result, "try-error")) {
    stop("Failed to fit Cox model: ", attr(fit_result, "condition")$message)
  }
  
  fit <- summary(fit_result)$conf.int
  
  # Validate output
  if (is.null(fit) || length(fit) < 4) {
    stop("Cox model fit did not return expected results")
  }
  
  list(
    hr = fit[1],
    lower = fit[3],
    upper = fit[4]
  )
}
