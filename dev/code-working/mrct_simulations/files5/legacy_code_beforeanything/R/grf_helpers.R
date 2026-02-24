#' Helper Functions for GRF Subgroup Analysis
#'
#' This file contains helper functions used by grf.subg.harm.survival()
#' to improve readability and modularity.

#' Create GRF configuration object
#'
#' Creates a configuration object to organize GRF parameters
#'
#' @param frac.tau Numeric. Fraction of tau for GRF horizon
#' @param n.min Integer. Minimum subgroup size
#' @param dmin.grf Numeric. Minimum difference in subgroup mean
#' @param RCT Logical. Is the data from a randomized controlled trial?
#' @param sg.criterion Character. Subgroup selection criterion
#' @param maxdepth Integer. Maximum tree depth
#' @param seedit Integer. Random seed
#' @return List with configuration parameters
#' @export

create_grf_config <- function(frac.tau, n.min, dmin.grf, RCT,
                              sg.criterion, maxdepth, seedit) {
  list(
    frac.tau = frac.tau,
    n.min = n.min,
    dmin.grf = dmin.grf,
    RCT = RCT,
    sg.criterion = sg.criterion,
    maxdepth = maxdepth,
    seedit = seedit,
    valid_criteria = c("mDiff", "Nsg"),
    max_tree_depth = 3,
    return_selected_cuts_only = FALSE
  )
}

#' Fit causal survival forest
#'
#' Wrapper function to fit GRF causal survival forest with appropriate settings
#'
#' @param X Matrix. Covariate matrix
#' @param Y Numeric vector. Outcome variable
#' @param W Numeric vector. Treatment indicator
#' @param D Numeric vector. Event indicator
#' @param tau.rmst Numeric. Time horizon for RMST
#' @param RCT Logical. Is this RCT data?
#' @param seedit Integer. Random seed
#' @return Causal survival forest object
#' @export

fit_causal_forest <- function(X, Y, W, D, tau.rmst, RCT, seedit) {
  if (RCT) {
    cs.forest <- try(
      suppressWarnings(
        grf::causal_survival_forest(X, Y, W, D, W.hat = 0.5,
                                    horizon = tau.rmst, seed = seedit)
      ),
      TRUE
    )
  } else {
    cs.forest <- try(
      suppressWarnings(
        grf::causal_survival_forest(X, Y, W, D,
                                    horizon = tau.rmst, seed = seedit)
      ),
      TRUE
    )
  }

  if (inherits(cs.forest, "try-error")) {
    stop("Failed to fit causal survival forest: ", attr(cs.forest, "condition")$message)
  }

  return(cs.forest)
}

#' Compute node metrics for a policy tree
#'
#' Aggregates scores by leaf node and calculates treatment effect differences
#'
#' @param data Data frame. Original data
#' @param dr.scores Matrix. Doubly robust scores
#' @param tree Policy tree object
#' @param X Matrix. Covariate matrix
#' @param n.min Integer. Minimum subgroup size
#' @return Data frame with node metrics
#' @importFrom stats aggregate
#' @export

compute_node_metrics <- function(data, dr.scores, tree, X, n.min) {
  data$predict.node <- predict(tree, X, type = "node.id")

  values <- stats::aggregate(
    dr.scores,
    by = list(leaf.node = data$predict.node),
    FUN = function(x) c(
      mean = mean(x),
      size = length(x),
      se = sd(x) / sqrt(length(x))
    )
  )

  values$diff <- values$control[, "mean"] - values$treated[, "mean"]
  values$Nsg <- values$control[, "size"]

  # Filter by minimum size requirement
  values <- values[values$control[, "size"] >= n.min, ]

  return(values)
}

#' Fit policy trees up to specified depth
#'
#' Fits policy trees of depths 1 through maxdepth and computes metrics
#'
#' @param X Matrix. Covariate matrix
#' @param data Data frame. Original data
#' @param dr.scores Matrix. Doubly robust scores
#' @param maxdepth Integer. Maximum tree depth (1-3)
#' @param n.min Integer. Minimum subgroup size
#' @return List with trees and combined values
#' @importFrom policytree policy_tree
#' @export

fit_policy_trees <- function(X, data, dr.scores, maxdepth, n.min) {
  trees <- list()
  all_values <- list()

  for (depth in 1:maxdepth) {
    trees[[depth]] <- policytree::policy_tree(X, dr.scores, depth = depth)

    values <- compute_node_metrics(data, dr.scores, trees[[depth]], X, n.min)
    values$depth <- depth
    all_values[[depth]] <- values
  }

  # Combine all values across depths
  combined_values <- do.call(rbind, all_values)

  return(list(
    trees = trees,
    values = combined_values
  ))
}

#' Select best subgroup based on criterion
#'
#' Identifies the optimal subgroup according to the specified criterion
#'
#' @param values Data frame. Node metrics from policy trees
#' @param sg.criterion Character. "mDiff" for maximum difference, "Nsg" for largest size
#' @param dmin.grf Numeric. Minimum difference threshold
#' @param n.max Integer. Maximum allowed subgroup size (total sample size)
#' @return Data frame row with best subgroup or NULL if none found
#' @export

select_best_subgroup <- function(values, sg.criterion, dmin.grf, n.max) {
  if (nrow(values) == 0) {
    return(NULL)
  }

  if (sg.criterion == "mDiff") {
    # Select subgroup with maximum difference
    loc.max <- which.max(values$diff)
    max.diff <- values[loc.max, , drop = FALSE]

  } else if (sg.criterion == "Nsg") {
    # Select largest subgroup meeting minimum difference threshold
    eligible <- values[values$diff >= dmin.grf, , drop = FALSE]

    if (nrow(eligible) == 0) {
      return(NULL)
    }

    loc.max <- which.max(eligible$Nsg)
    max.diff <- eligible[loc.max, , drop = FALSE]

  } else {
    stop("Unknown sg.criterion: ", sg.criterion)
  }

  # Validate subgroup: must meet difference threshold and not be entire population
  if (max.diff$diff < dmin.grf || max.diff$Nsg >= n.max) {
    return(NULL)
  }

  return(max.diff)
}

#' Extract cut information from a policy tree
#'
#' Extracts all split points and variables from a policy tree
#'
#' @param tree Policy tree object
#' @return List with cuts (expressions) and names (unique variables)
#' @keywords internal

extract_tree_cuts <- function(tree) {
  grf_names <- tree$columns
  tnodes <- tree$nodes
  cuts <- character()
  names <- character()

  for (tt in seq_along(tnodes)) {
    node <- tnodes[[tt]]
    if (!node$is_leaf) {
      split_variable_index <- node$split_variable
      split_value <- round(node$split_value, 2)

      variable_name <- grf_names[split_variable_index]
      cut_expression <- paste0(variable_name, " <= ", split_value)

      cuts <- c(cuts, cut_expression)
      names <- c(names, variable_name)
    }
  }

  return(list(
    cuts = cuts,
    names = unique(names)
  ))
}

#' Find the split that leads to a specific leaf node
#'
#' Identifies the split point that creates a given leaf node
#'
#' @param tree Policy tree object
#' @param leaf_node Integer. Leaf node identifier
#' @return Character string with split expression or NULL
#' @export

find_leaf_split <- function(tree, leaf_node) {
  grf_names <- tree$columns
  tnodes <- tree$nodes

  for (tt in seq_along(tnodes)) {
    node <- tnodes[[tt]]
    if (!node$is_leaf) {
      if (node$left_child == leaf_node || node$right_child == leaf_node) {
        split_variable_index <- node$split_variable
        split_value <- node$split_value
        variable_name <- grf_names[split_variable_index]

        return(paste0(variable_name, " <= ", split_value))
      }
    }
  }

  return(NULL)
}

#' Extract all cuts from fitted trees
#'
#' Consolidates cut information from all fitted policy trees. This is the
#' default behavior that returns cuts from all trees regardless of which
#' tree identified the selected subgroup.
#'
#' @param trees List. Policy trees (indexed by depth)
#' @param maxdepth Integer. Maximum tree depth
#' @return List with cuts and names for each tree and combined
#' @keywords internal

extract_all_tree_cuts <- function(trees, maxdepth) {
  result <- list()
  all_cuts <- character()
  all_names <- character()

  for (depth in 1:maxdepth) {
    if (!is.null(trees[[depth]])) {
      tree_info <- extract_tree_cuts(trees[[depth]])

      # Store tree-specific cuts
      result[[paste0("tree", depth)]] <- tree_info$cuts
      result[[paste0("names", depth)]] <- tree_info$names

      # Accumulate all cuts
      all_cuts <- c(all_cuts, tree_info$cuts)
      all_names <- c(all_names, tree_info$names)
    }
  }

  result$all <- unique(all_cuts)
  result$all_names <- unique(all_names)

  return(result)
}

#' Extract cuts from selected tree only
#'
#' Extracts cut information only from the tree at the specified selected depth.
#' This provides a focused set of cuts from the tree that identified the
#' subgroup meeting the `dmin.grf` criterion, rather than cuts from all trees.
#'
#' @param trees List. Policy trees (indexed by depth)
#' @param selected_depth Integer. Depth of the selected tree (from best_subgroup$depth)
#' @param maxdepth Integer. Maximum tree depth (for populating tree-specific slots)
#' @return List with cuts and names, structured similarly to extract_all_tree_cuts
#'   but with only the selected tree's cuts in the `all` field
#'
#' @details
#' This function is used when `return_selected_cuts_only = TRUE` in
#' `grf.subg.harm.survival()`. It returns:
#' \itemize{
#'   \item `tree1`, `tree2`, `tree3`: Individual tree cuts (still populated for reference)
#'   \item `names1`, `names2`, `names3`: Individual tree variable names
#'   \item `all`: Cuts from the SELECTED tree only (not union of all trees)
#'   \item `all_names`: Variable names from the SELECTED tree only
#'   \item `selected_depth`: The depth that was selected
#' }
#'
#' @examples
#' \dontrun{
#' # When best_subgroup$depth = 2, only tree2 cuts are returned in $all
#' tree_cuts <- extract_selected_tree_cuts(trees, selected_depth = 2, maxdepth = 2)
#' tree_cuts$all
#' # Returns only cuts from depth 2 tree
#' }
#'
#' @keywords internal

extract_selected_tree_cuts <- function(trees, selected_depth, maxdepth) {
  result <- list()

  # Still extract cuts from all trees for reference
  for (depth in 1:maxdepth) {
    if (!is.null(trees[[depth]])) {
      tree_info <- extract_tree_cuts(trees[[depth]])

      # Store tree-specific cuts
      result[[paste0("tree", depth)]] <- tree_info$cuts
      result[[paste0("names", depth)]] <- tree_info$names
    }
  }

  # For the combined 'all' field, only use the selected tree
  if (!is.null(trees[[selected_depth]])) {
    selected_tree_info <- extract_tree_cuts(trees[[selected_depth]])
    result$all <- unique(selected_tree_info$cuts)
    result$all_names <- unique(selected_tree_info$names)
  } else {
    result$all <- character()
    result$all_names <- character()
  }

  # Record which depth was selected

result$selected_depth <- selected_depth

  return(result)
}

#' Assign data to subgroups based on selected node
#'
#' Creates treatment recommendation flags based on identified subgroup
#'
#' @param data Data frame. Original data
#' @param best_subgroup Data frame row. Selected subgroup information
#' @param trees List. Policy trees
#' @param X Matrix. Covariate matrix
#' @return Data frame with added predict.node and treat.recommend columns
#' @keywords internal

assign_subgroup_membership <- function(data, best_subgroup, trees, X) {
  depth <- best_subgroup$depth
  selected_tree <- trees[[depth]]

  # Assign nodes based on selected tree
  data$predict.node <- predict(selected_tree, X, type = "node.id")

  # Create treatment recommendation flag
  # 0 = belongs to identified harm subgroup
  # 1 = belongs to complement (recommend treatment)
  data$treat.recommend <- ifelse(
    data$predict.node == best_subgroup$leaf.node,
    0,
    1
  )

  return(data)
}

#' Create result object for successful subgroup identification
#'
#' Builds comprehensive result object when a subgroup is found
#'
#' @param data Data frame. Original data with subgroup assignments
#' @param best_subgroup Data frame row. Selected subgroup information
#' @param trees List. All fitted policy trees
#' @param tree_cuts List. Cut information from trees
#' @param selected_tree Policy tree. The tree that identified the subgroup
#' @param sg_harm_id Character. Expression defining the subgroup
#' @param values Data frame. All node metrics
#' @param config List. GRF configuration
#' @return List with complete GRF results
#' @keywords internal

create_success_result <- function(data, best_subgroup, trees, tree_cuts,
                                  selected_tree, sg_harm_id, values, config) {

  # Find all subgroups with positive difference
  harm_any <- values[values$diff > 0, ]

  result <- list(
    # Core results
    data = data,
    grf.gsub = best_subgroup,
    sg.harm.id = sg_harm_id,

    # Tree cuts
    tree.cuts = tree_cuts$all,
    tree.names = tree_cuts$all_names,

    # Tree-specific cuts
    tree1.cuts = tree_cuts$tree1,
    tree1.names = tree_cuts$names1,

    # Additional subgroups
    harm.any = harm_any,

    # Tree objects
    tree = selected_tree,
    tree1 = trees[[1]],

    # Configuration
    tau.rmst = config$tau.rmst,
    dmin.grf = config$dmin.grf,
    frac.tau = config$frac.tau,
    maxdepth = config$maxdepth,
    n.min = config$n.min,

    # Track selected depth and cut extraction mode
    selected_depth = best_subgroup$depth,
    return_selected_cuts_only = config$return_selected_cuts_only
  )

  # Add tree2 cuts if maxdepth >= 2
  if (config$maxdepth >= 2) {
    result$tree2.cuts <- tree_cuts$tree2
    result$tree2.names <- tree_cuts$names2
    result$tree2 <- trees[[2]]
  }

  # Add tree3 cuts if maxdepth >= 3
  if (config$maxdepth >= 3) {
    result$tree3.cuts <- tree_cuts$tree3
    result$tree3.names <- tree_cuts$names3
    result$tree3 <- trees[[3]]
  }

  return(result)
}

#' Create result object when no subgroup is found
#'
#' Builds result object for cases where no valid subgroup is identified
#'
#' @param data Data frame. Original data
#' @param values Data frame. Node metrics (may be empty)
#' @param trees List. Fitted policy trees
#' @param config List. GRF configuration
#' @return List with limited GRF results
#' @keywords internal

create_null_result <- function(data, values, trees, config) {

  # Find any subgroups with positive difference (for diagnostic purposes)
  harm_any <- NULL
  if (!is.null(values) && nrow(values) > 0) {
    cand.sgs <- which(values$diff > 0)
    if (length(cand.sgs) > 0) {
      harm_any <- values[cand.sgs, ]
    }
  }

  # Use the first tree as default (or the deepest available)
  selected_tree <- trees[[1]]
  if (config$maxdepth > 1 && !is.null(trees[[config$maxdepth]])) {
    selected_tree <- trees[[config$maxdepth]]
  }

  result <- list(
    data = data,
    grf.gsub = NULL,
    sg.harm.id = NULL,
    tree.cuts = character(),
    tree.names = character(),
    harm.any = harm_any,
    tree = selected_tree,
    tau.rmst = config$tau.rmst,
    dmin.grf = config$dmin.grf,
    frac.tau = config$frac.tau,
    maxdepth = config$maxdepth,
    n.min = config$n.min,
    selected_depth = NULL,
    return_selected_cuts_only = config$return_selected_cuts_only
  )

  # Add tree objects if they exist
  for (i in 1:config$maxdepth) {
    if (!is.null(trees[[i]])) {
      result[[paste0("tree", i)]] <- trees[[i]]
    }
  }

  return(result)
}

#' Print detailed output for debugging
#'
#' Displays detailed information about the GRF analysis
#'
#' @param config List. GRF configuration
#' @param values Data frame. Node metrics
#' @param best_subgroup Data frame row. Selected subgroup (or NULL)
#' @param sg_harm_id Character. Subgroup definition (or NULL)
#' @param tree_cuts List. Cut information
#' @export

print_grf_details <- function(config, values, best_subgroup, sg_harm_id, tree_cuts = NULL) {
  cat("tau, maxdepth =", c(config$tau.rmst, config$maxdepth), "\n")

  if (!is.null(values) && nrow(values) > 0) {
    temp <- values[, c("leaf.node", "control", "depth")]
    print(round(temp, 2))
  }

  if (!is.null(best_subgroup)) {
    temp <- best_subgroup[, c("leaf.node", "control", "depth")]
    cat("\nSelected subgroup:\n")
    print(round(temp, 2))

    if (!is.null(sg_harm_id)) {
      cat("\nGRF subgroup found\n")
      cat("Terminating node at max.diff (sg.harm.id):\n")
      print(sg_harm_id)

      if (!is.null(tree_cuts)) {
        if (isTRUE(config$return_selected_cuts_only)) {
          cat("\nCuts from selected tree (depth =", best_subgroup$depth, "):\n")
        } else {
          cat("\nAll splits (from all trees):\n")
        }
        print(tree_cuts$all)
      }
    }
  } else {
    cat("GRF subgroup NOT found\n")
  }
}

#' Validate input data for GRF analysis
#'
#' Checks that input data meets requirements for GRF analysis
#'
#' @param W Numeric vector. Treatment indicator
#' @param D Numeric vector. Event indicator
#' @param n.min Integer. Minimum subgroup size
#' @return Logical. TRUE if data is valid, FALSE with warning otherwise
#' @export

validate_grf_data <- function(W, D, n.min) {
  # Check treatment variation
  if (length(unique(W)) < 2) {
    warning("Treatment variable has only one unique value. Cannot identify subgroups.")
    return(FALSE)
  }

  # Check event variation
  if (length(unique(D)) < 2) {
    warning("Event variable has only one unique value. Cannot identify subgroups.")
    return(FALSE)
  }

  # Check minimum sample size per arm
  n_treat <- sum(W == 1)
  n_control <- sum(W == 0)

  if (n_treat < n.min || n_control < n.min) {
    warning(sprintf(
      "Insufficient sample size: treatment=%d, control=%d, required=%d per arm.",
      n_treat, n_control, n.min
    ))
    return(FALSE)
  }

  return(TRUE)
}
