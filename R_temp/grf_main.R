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
#'
#' @return A list with GRF results, including:
#'   \item{data}{Original data with added treatment recommendation flags}
#'   \item{grf.gsub}{Selected subgroup information}
#'   \item{sg.harm.id}{Expression defining the identified subgroup}
#'   \item{tree.cuts}{All cut expressions from the selected tree}
#'   \item{tree.names}{Unique variable names used in cuts}
#'   \item{tree}{Selected policy tree object}
#'   \item{tau.rmst}{Time horizon used for RMST}
#'   \item{harm.any}{All subgroups with positive treatment effect difference}
#'   Additional tree-specific cuts and objects (tree1, tree2, tree3) based on maxdepth
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
                                   seedit = 8316951) {

  # =========================================================================
  # SECTION: INPUT VALIDATION
  # Purpose: Validate all input parameters before processing
  # =========================================================================

  if (maxdepth > 3) {
    stop("Maximum depth cannot exceed 3")
  }

  valid_criteria <- c("mDiff", "Nsg")
  if (!sg.criterion %in% valid_criteria) {
    stop("sg.criterion must be one of: ", paste(valid_criteria, collapse = ", "))
  }

  # =========================================================================
  # SECTION: CONFIGURATION SETUP
  # Purpose: Create configuration object for consistent parameter passing
  # =========================================================================

  config <- create_grf_config(
    frac.tau = frac.tau,
    n.min = n.min,
    dmin.grf = dmin.grf,
    RCT = RCT,
    sg.criterion = sg.criterion,
    maxdepth = maxdepth,
    seedit = seedit
  )

  # =========================================================================
  # SECTION: DATA PREPARATION
  # Purpose: Convert data to appropriate format for GRF analysis
  # =========================================================================

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

  # =========================================================================
  # SECTION: CAUSAL FOREST FITTING
  # Purpose: Fit GRF causal survival forest to identify treatment heterogeneity
  # =========================================================================

  cs.forest <- fit_causal_forest(X, Y, W, D, tau.rmst, config$RCT, config$seedit)

  # =========================================================================
  # SECTION: SUBGROUP IDENTIFICATION VIA POLICY TREES
  # Purpose: Use policy trees to partition the covariate space
  # =========================================================================

  # Compute doubly robust scores for subgroup identification
  dr.scores <- policytree::double_robust_scores(cs.forest)

  # Maximum sample size (used to exclude full population as subgroup)
  n.max <- length(Y)

  # Fit policy trees and compute metrics
  tree_results <- fit_policy_trees(X, data, dr.scores, config$maxdepth, config$n.min)
  trees <- tree_results$trees
  values <- tree_results$values

  # =========================================================================
  # SECTION: OPTIMAL SUBGROUP SELECTION
  # Purpose: Choose the best subgroup based on specified criterion
  # =========================================================================

  best_subgroup <- select_best_subgroup(
    values = values,
    sg.criterion = config$sg.criterion,
    dmin.grf = config$dmin.grf,
    n.max = n.max
  )

  # =========================================================================
  # SECTION: RESULT COMPILATION - NO SUBGROUP FOUND
  # Purpose: Return appropriate result when no valid subgroup is identified
  # =========================================================================

  if (is.null(best_subgroup)) {
    if (details) {
      print_grf_details(config, values, NULL, NULL)
    }

    return(create_null_result(data, values, trees, config))
  }

  # =========================================================================
  # SECTION: RESULT COMPILATION - SUBGROUP FOUND
  # Purpose: Extract subgroup information and create comprehensive result
  # =========================================================================

  # Assign data points to subgroups
  data <- assign_subgroup_membership(data, best_subgroup, trees, X)

  # Select the tree that identified the best subgroup
  selected_tree <- trees[[best_subgroup$depth]]

  # Find the specific split that defines the subgroup
  sg_harm_id <- find_leaf_split(selected_tree, best_subgroup$leaf.node)

  # Extract all cuts from fitted trees
  tree_cuts <- extract_all_tree_cuts(trees, config$maxdepth)

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


#' GRF Subgroup Evaluation and Performance Metrics
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
#' @return Data frame with performance metrics including:
#'   \item{any.H}{Indicator for whether a subgroup was found}
#'   \item{size.H}{Size of identified harm subgroup}
#'   \item{size.Hc}{Size of complement subgroup}
#'   \item{ppv}{Positive predictive value}
#'   \item{npv}{Negative predictive value}
#'   \item{specificity}{Specificity of subgroup identification}
#'   \item{sensitivity}{Sensitivity of subgroup identification}
#'   \item{hr.H.true}{True hazard ratio in harm subgroup}
#'   \item{hr.H.hat}{Estimated hazard ratio in identified subgroup}
#'   \item{b1.H}{Bias type 1 for harm subgroup}
#'   \item{b2.H}{Bias type 2 for harm subgroup}
#'   Additional metrics for complement subgroup and confidence intervals
#'
#' @importFrom survival coxph Surv
#' @export

grf.estimates.out <- function(df,
                              grf.est = NULL,
                              dgm = NULL,
                              cox.formula.sim = NULL,
                              cox.formula.adj.sim = NULL,
                              analysis = "GRF",
                              frac.tau = 1.0) {

  # Calculate censoring proportion
  p.cens <- mean(1 - df[, event.name])

  # Extract outcome variables
  Y <- df[, outcome.name]
  W <- df[, treat.name]
  D <- df[, event.name]

  # Calculate tau maximum
  taumax <- frac.tau * min(
    max(Y[W == 1 & D == 1]),
    max(Y[W == 0 & D == 1])
  )

  # Extract subgroup identifier
  sg.harm.grf <- if (!is.null(grf.est)) grf.est$sg.harm.id else NULL

  # =========================================================================
  # Calculate ITT estimates
  # =========================================================================

  fit <- summary(coxph(cox.formula.sim, data = df, robust = FALSE))$conf.int
  hr.itt <- fit[1]
  l.itt <- fit[3]
  u.itt <- fit[4]

  fit <- summary(coxph(cox.formula.adj.sim, data = df, robust = FALSE))$conf.int
  hr.adj.itt <- fit[1, 1]
  l.adj.itt <- fit[1, 3]
  u.adj.itt <- fit[1, 4]

  # =========================================================================
  # Validate DGM consistency
  # =========================================================================

  if (dgm$model == "null" && !is.null(dgm$grf.harm.true)) {
    stop("For dgm model null, grf.harm.true should be null")
  }

  # =========================================================================
  # Case 1: True subgroup exists AND subgroup was found
  # =========================================================================

  if (!is.null(dgm$grf.harm.true) && !is.null(sg.harm.grf)) {
    any.H <- 1.0
    dfout <- grf.est$data

    # Calculate confusion matrix
    aa <- sum(dfout$treat.recommend == 0 & dfout$flag.harm == 1)  # True positive
    bb <- sum(dfout$treat.recommend == 1 & dfout$flag.harm == 1)  # False negative
    cc <- sum(dfout$treat.recommend == 0 & dfout$flag.harm == 0)  # False positive
    dd <- sum(dfout$treat.recommend == 1 & dfout$flag.harm == 0)  # True negative

    size.H <- sum(dfout$treat.recommend == 0)
    size.Hc <- sum(dfout$treat.recommend == 1)

    # Calculate HR for true subgroups
    fit <- summary(coxph(cox.formula.sim, data = subset(df, flag.harm == 1), robust = FALSE))$conf.int
    hr.H.true <- fit[1]
    l.H.true <- fit[3]
    u.H.true <- fit[4]

    fit <- summary(coxph(cox.formula.sim, data = subset(df, flag.harm == 0), robust = FALSE))$conf.int
    hr.Hc.true <- fit[1]
    l.Hc.true <- fit[3]
    u.Hc.true <- fit[4]

    # Calculate HR for identified subgroups
    fit <- summary(coxph(cox.formula.sim, data = subset(dfout, treat.recommend == 0), robust = FALSE))$conf.int
    hr.H.hat <- fit[1]
    l.H.hat <- fit[3]
    u.H.hat <- fit[4]

    fit <- summary(coxph(cox.formula.sim, data = subset(dfout, treat.recommend == 1), robust = FALSE))$conf.int
    hr.Hc.hat <- fit[1]
    l.Hc.hat <- fit[3]
    u.Hc.hat <- fit[4]

    # Calculate bias metrics
    b1.H <- hr.H.hat - hr.H.true
    b2.H <- hr.H.hat - dgm$hr.H.true
    b1.Hc <- hr.Hc.hat - hr.Hc.true
    b2.Hc <- hr.Hc.hat - dgm$hr.Hc.true

    # Calculate performance metrics
    ppv <- aa / (aa + bb)
    npv <- dd / (cc + dd)
    specificity <- dd / (bb + dd)
    sensitivity <- aa / (aa + cc)

    found.1 <- found.2 <- found.both <- found.al3 <- NA

    # =========================================================================
    # Case 2: True subgroup exists BUT no subgroup was found
    # =========================================================================

  } else if (!is.null(dgm$grf.harm.true) && is.null(sg.harm.grf)) {
    any.H <- 0
    size.H <- 0
    size.Hc <- nrow(df)

    # Calculate HR for true subgroups
    fit <- summary(coxph(cox.formula.sim, data = subset(df, flag.harm == 1), robust = FALSE))$conf.int
    hr.H.true <- fit[1]
    l.H.true <- fit[3]
    u.H.true <- fit[4]

    fit <- summary(coxph(cox.formula.sim, data = subset(df, flag.harm == 0), robust = FALSE))$conf.int
    hr.Hc.true <- fit[1]
    l.Hc.true <- fit[3]
    u.Hc.true <- fit[4]

    # No identified subgroup for H
    hr.H.hat <- l.H.hat <- u.H.hat <- NA

    # Hc is entire population (ITT)
    fit <- summary(coxph(cox.formula.sim, data = df, robust = FALSE))$conf.int
    hr.Hc.hat <- fit[1]
    l.Hc.hat <- fit[3]
    u.Hc.hat <- fit[4]

    # Calculate bias metrics
    b1.H <- b2.H <- NA
    b1.Hc <- hr.Hc.hat - hr.Hc.true
    b2.Hc <- hr.Hc.hat - dgm$hr.Hc.true

    # Performance metrics when no subgroup found
    aa <- 0  # No true positives
    bb <- sum(df$flag.harm == 1)  # All harm group misclassified
    cc <- 0  # No false positives
    dd <- sum(df$flag.harm == 0)  # All non-harm correctly classified

    ppv <- 0
    npv <- dd / (cc + dd)
    specificity <- dd / (bb + dd)
    sensitivity <- 0

    found.1 <- found.2 <- found.both <- NA
    found.al3 <- 0

    # =========================================================================
    # Case 3: No true subgroup BUT subgroup was found (Type I error)
    # =========================================================================

  } else if (is.null(dgm$grf.harm.true) && !is.null(sg.harm.grf)) {
    any.H <- 1.0
    dfout <- grf.est$data

    # All data is truly non-harm (flag.harm = 0 for all)
    aa <- 0  # No true positives possible
    bb <- 0  # No false negatives possible
    cc <- sum(dfout$treat.recommend == 0)  # All identified as harm are false positives
    dd <- sum(dfout$treat.recommend == 1)  # All identified as non-harm are true negatives

    size.H <- sum(dfout$treat.recommend == 0)
    size.Hc <- sum(dfout$treat.recommend == 1)

    # True values (ITT since no real subgroup)
    hr.H.true <- l.H.true <- u.H.true <- NA

    fit <- summary(coxph(cox.formula.sim, data = df, robust = FALSE))$conf.int
    hr.Hc.true <- fit[1]
    l.Hc.true <- fit[3]
    u.Hc.true <- fit[4]

    # Estimated values
    fit <- summary(coxph(cox.formula.sim, data = subset(dfout, treat.recommend == 0), robust = FALSE))$conf.int
    hr.H.hat <- fit[1]
    l.H.hat <- fit[3]
    u.H.hat <- fit[4]

    fit <- summary(coxph(cox.formula.sim, data = subset(dfout, treat.recommend == 1), robust = FALSE))$conf.int
    hr.Hc.hat <- fit[1]
    l.Hc.hat <- fit[3]
    u.Hc.hat <- fit[4]

    # Bias metrics
    b1.H <- b2.H <- NA
    b1.Hc <- hr.Hc.hat - hr.Hc.true
    b2.Hc <- hr.Hc.hat - dgm$hr.Hc.true

    # Performance metrics
    ppv <- NA
    npv <- dd / (cc + dd)
    specificity <- dd / (bb + dd)
    sensitivity <- 0

    found.1 <- found.2 <- found.both <- found.al3 <- NA

    # =========================================================================
    # Case 4: No true subgroup AND no subgroup found (Correct null)
    # =========================================================================

  } else {
    any.H <- 0

    # All classified as non-harm (treat.recommend = 1 for all)
    aa <- bb <- cc <- 0
    dd <- nrow(df)

    size.H <- 0
    size.Hc <- nrow(df)

    # All estimates are ITT
    hr.H.true <- l.H.true <- u.H.true <- NA
    hr.H.hat <- l.H.hat <- u.H.hat <- NA

    fit <- summary(coxph(cox.formula.sim, data = df, robust = FALSE))$conf.int
    hr.Hc.true <- hr.Hc.hat <- fit[1]
    l.Hc.true <- l.Hc.hat <- fit[3]
    u.Hc.true <- u.Hc.hat <- fit[4]

    # Bias metrics
    b1.H <- b2.H <- NA
    b1.Hc <- 0  # No bias when correctly identifying no subgroup
    b2.Hc <- hr.Hc.hat - dgm$hr.Hc.true

    # Performance metrics
    ppv <- sensitivity <- NA
    npv <- specificity <- 1.0

    found.1 <- found.2 <- found.both <- found.al3 <- NA
  }

  # =========================================================================
  # Compile results
  # =========================================================================

  df.res <- data.frame(
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
    p.cens = p.cens,
    analysis = analysis,
    taumax = taumax,
    hr.itt = hr.itt,
    l.itt = l.itt,
    u.itt = u.itt,
    hr.adj.itt = hr.adj.itt,
    l.adj.itt = l.adj.itt,
    u.adj.itt = u.adj.itt,
    l.H.true = l.H.true,
    u.H.true = u.H.true,
    l.Hc.true = l.Hc.true,
    u.Hc.true = u.Hc.true,
    l.H.hat = l.H.hat,
    u.H.hat = u.H.hat,
    l.Hc.hat = l.Hc.hat,
    u.Hc.hat = u.Hc.hat
  )

  return(df.res)
}
