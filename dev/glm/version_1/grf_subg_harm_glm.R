# =============================================================================
# grf_subg_harm_glm.R
#
# GRF-based candidate factor selection for GLM outcomes (binary, continuous).
# Mirrors grf.subg.harm.survival() but uses causal_forest() instead of
# causal_survival_forest(), and a tau-threshold on the risk-difference scale
# for binary outcomes, or a mean-difference scale for continuous outcomes.
# =============================================================================


#' GRF Candidate Factor Selection for GLM Outcomes
#'
#' Uses a Generalized Random Forest causal forest (`grf::causal_forest()`) to
#' estimate conditional average treatment effects (CATEs) and then fits a
#' policy tree to identify binary splitting rules indicative of a detrimental
#' (or beneficial) subgroup.  Designed for binary (0/1) or continuous numeric
#' outcomes.
#'
#' For binary outcomes the CATE is on the **risk-difference** scale.  For
#' continuous outcomes it is on the **mean-difference** scale.  The `tau.threshold`
#' argument specifies the CATE magnitude below which a leaf is considered to
#' favour the control arm (i.e., potential harm from treatment).
#'
#' @param data          A `data.frame` containing all analysis variables.
#' @param confounders.name Character vector of covariate names used as GRF
#'   predictors.
#' @param outcome.name  Character. Name of the outcome column (binary 0/1 or
#'   numeric).
#' @param treat.name    Character. Name of the binary treatment indicator (0/1).
#' @param id.name       Character or `NULL`. Row-identifier column name.
#' @param outcome_type  Character. `"binary"` or `"continuous"`.
#' @param maxdepth      Integer. Maximum depth of the policy tree (default 2).
#' @param n.min         Integer. Minimum leaf size in the policy tree (default 60).
#' @param tau.threshold Numeric. CATE threshold below which a leaf is flagged as
#'   potentially harmful.  For binary outcomes this is a risk difference (e.g.,
#'   `0.05`); for continuous outcomes it is a mean difference (e.g., `2.0`).
#'   Default `0` selects leaves with any excess risk for treatment.
#' @param num.trees     Integer. Number of trees in the causal forest (default 2000).
#' @param details       Logical. If `TRUE`, return the full GRF fit and policy
#'   trees in addition to the candidate-cut summary (default `FALSE`).
#' @param return_selected_cuts_only Logical. If `TRUE`, return only the factor
#'   cuts meeting the `tau.threshold` criterion (default `TRUE`).
#' @param seed          Integer or `NULL`. Random seed for reproducibility.
#' @param ...           Additional arguments passed to `grf::causal_forest()`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`selected_cuts`}{A `data.frame` of selected binary factor definitions
#'     (variable name, cut direction, cut value).}
#'   \item{`tau_hat`}{Numeric vector of individual CATE estimates (length `nrow(data)`).}
#'   \item{`tree1`}{Policy tree of depth 1 (if `details = TRUE`).}
#'   \item{`tree2`}{Policy tree of depth `maxdepth` (if `details = TRUE`).}
#'   \item{`grf_fit`}{The `causal_forest` object (if `details = TRUE`).}
#' }
#'
#' @seealso [grf.subg.harm.survival()] for the survival-endpoint counterpart.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   n <- 400
#'   df <- data.frame(
#'     id    = seq_len(n),
#'     treat = rbinom(n, 1, 0.5),
#'     x1    = rnorm(n),
#'     x2    = rbinom(n, 1, 0.5),
#'     resp  = rbinom(n, 1, 0.3)
#'   )
#'   grf_res <- grf.subg.harm.glm(
#'     data             = df,
#'     confounders.name = c("x1", "x2"),
#'     outcome.name     = "resp",
#'     treat.name       = "treat",
#'     id.name          = "id",
#'     outcome_type     = "binary",
#'     tau.threshold    = 0.0,
#'     num.trees        = 500,
#'     details          = TRUE
#'   )
#' }
#'
#' @export
grf.subg.harm.glm <- function(
    data,
    confounders.name,
    outcome.name,
    treat.name,
    id.name          = NULL,
    outcome_type     = c("binary", "continuous"),
    maxdepth         = 2L,
    n.min            = 60L,
    tau.threshold    = 0,
    num.trees        = 2000L,
    details          = FALSE,
    return_selected_cuts_only = TRUE,
    seed             = NULL,
    ...
) {
  outcome_type <- match.arg(outcome_type)

  # ---- Input validation ----------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  required_cols <- c(confounders.name, outcome.name, treat.name)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "The following columns are missing from `data`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  treat_vals <- unique(data[[treat.name]])
  if (!all(treat_vals %in% c(0L, 1L, 0, 1))) {
    stop(
      "`treat.name` column must contain only 0 and 1 values.",
      call. = FALSE
    )
  }
  if (outcome_type == "binary") {
    outcome_vals <- unique(data[[outcome.name]])
    if (!all(outcome_vals %in% c(0L, 1L, 0, 1, NA))) {
      stop(
        "For `outcome_type = 'binary'`, the outcome column must contain only 0, 1, or NA.",
        call. = FALSE
      )
    }
  }

  # ---- Prepare GRF inputs --------------------------------------------------
  complete_idx <- stats::complete.cases(
    data[, c(confounders.name, outcome.name, treat.name), drop = FALSE]
  )
  data_cc <- data[complete_idx, , drop = FALSE]

  X <- as.matrix(data_cc[, confounders.name, drop = FALSE])
  Y <- as.numeric(data_cc[[outcome.name]])
  W <- as.numeric(data_cc[[treat.name]])

  # ---- Fit causal forest ---------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  cf <- grf::causal_forest(
    X         = X,
    Y         = Y,
    W         = W,
    num.trees = as.integer(num.trees),
    ...
  )

  tau_hat <- predict(cf)$predictions   # risk diff (binary) or mean diff (continuous)

  # ---- Policy tree depth 1 -------------------------------------------------
  tree1 <- policytree::policy_tree(
    X      = X,
    Gamma  = policytree::double_robust_scores(cf),
    depth  = 1L,
    min.node.size = as.integer(n.min)
  )

  # ---- Policy tree depth maxdepth ------------------------------------------
  tree2 <- if (maxdepth >= 2L) {
    policytree::policy_tree(
      X      = X,
      Gamma  = policytree::double_robust_scores(cf),
      depth  = as.integer(maxdepth),
      min.node.size = as.integer(n.min)
    )
  } else {
    tree1
  }

  # ---- Extract candidate cuts from trees -----------------------------------
  cuts_list <- lapply(list(tree1, tree2), function(tree) {
    .extract_grf_cuts_glm(
      tree             = tree,
      X                = X,
      tau_hat          = tau_hat,
      confounders.name = confounders.name,
      tau.threshold    = tau.threshold
    )
  })
  cuts_all <- do.call(rbind, cuts_list)
  # De-duplicate
  if (!is.null(cuts_all) && nrow(cuts_all) > 0L) {
    cuts_all <- unique(cuts_all)
  }

  if (return_selected_cuts_only && (!is.null(cuts_all)) && nrow(cuts_all) > 0L) {
    selected_cuts <- cuts_all[cuts_all$selected, , drop = FALSE]
  } else {
    selected_cuts <- cuts_all
  }

  # ---- Return --------------------------------------------------------------
  out <- list(
    selected_cuts = selected_cuts,
    tau_hat       = tau_hat
  )
  if (details) {
    out$tree1    <- tree1
    out$tree2    <- tree2
    out$grf_fit  <- cf
  }
  out
}


# -----------------------------------------------------------------------------
# Internal helper: extract binary cut definitions from a policy tree
# -----------------------------------------------------------------------------

#' @noRd
.extract_grf_cuts_glm <- function(
    tree,
    X,
    tau_hat,
    confounders.name,
    tau.threshold
) {
  # policytree stores the tree as a list; extract internal nodes
  nodes <- tree$nodes
  if (is.null(nodes)) return(NULL)

  cut_rows <- lapply(seq_along(nodes), function(i) {
    node <- nodes[[i]]
    if (is.null(node$split_variable) || node$is_leaf) return(NULL)

    var_idx   <- node$split_variable        # 1-based index into X columns
    cut_val   <- node$split_value
    var_name  <- confounders.name[var_idx]

    # Mean CATE in the "harm" leaf (treatment effect above threshold)
    in_leaf   <- X[, var_idx] <= cut_val
    tau_leaf  <- mean(tau_hat[in_leaf], na.rm = TRUE)
    selected  <- tau_leaf >= tau.threshold

    data.frame(
      variable  = var_name,
      direction = "<=",
      cut_value = cut_val,
      mean_tau  = tau_leaf,
      selected  = selected,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), cut_rows))
}
