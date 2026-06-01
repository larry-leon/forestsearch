# =============================================================================
# grf_subgroup_labels.R - Correct subgroup-definition construction for GRF
# =============================================================================
#
# GRF policy trees define the identified subgroup as a set of leaves.  The
# textual subgroup definition must be the EXACT root-to-leaf decision path(s)
# for those leaves, with correct per-split directions:
#   * left  child  =>  "var <= value"
#   * right child  =>  "var >  value"
# A subgroup spanning multiple leaves is therefore a DISJUNCTION of
# conjunctions (a union of axis-aligned boxes), which a single AND-composed
# label vector cannot represent.
#
# The helpers here build that definition directly from the policy-tree
# topology (node `left_child`/`right_child` give the path; `split_variable`
# indexes `tree$columns`; `predict(tree, X, type = "node.id")` returns the
# 1-based position in `tree$nodes`, so leaf targets are node positions).
#
# Two earlier builders were structurally unable to do this and produced
# subgroup definitions that did not match GRF's own membership whenever a
# path turned right or the subgroup spanned more than one leaf:
#   * find_leaf_split()    (survival): returned only the single split
#                          immediately above one leaf.
#   * .build_sg_harm_id()  (GLM):      emitted every internal-node split as
#                          "<=", ignoring path and direction.
# Both now delegate here.
# =============================================================================


# ---------------------------------------------------------------------------
# Internal: parent/direction map for a policy_tree
# ---------------------------------------------------------------------------
# Returns, for each node position, its parent position and whether it is the
# parent's left ("<=") or right (">") child.  Root has parent 0L.
#' @noRd
.grf_parent_map <- function(tree) {
  nodes <- tree$nodes
  k <- length(nodes)
  parent <- integer(k); side <- character(k)
  parent[] <- 0L; side[] <- ""
  for (i in seq_len(k)) {
    nd <- nodes[[i]]
    if (isTRUE(nd$is_leaf)) next
    lc <- nd$left_child; rc <- nd$right_child
    if (!is.null(lc) && lc >= 1L && lc <= k) { parent[lc] <- i; side[lc] <- "left" }
    if (!is.null(rc) && rc >= 1L && rc <= k) { parent[rc] <- i; side[rc] <- "right" }
  }
  list(parent = parent, side = side)
}


# ---------------------------------------------------------------------------
# Internal: root-to-leaf conjunction for one leaf
# ---------------------------------------------------------------------------
# Walks from `leaf` up to the root, recording each ancestor split as a
# (variable, op, value) triple with the direction taken to reach the leaf.
# Returns a data.frame of components in root-to-leaf order.
#' @noRd
.grf_leaf_conjunction <- function(tree, leaf, pmap) {
  nodes   <- tree$nodes
  columns <- tree$columns
  comps   <- list()
  cur     <- leaf
  guard   <- 0L
  while (cur >= 1L && pmap$parent[cur] != 0L) {
    par   <- pmap$parent[cur]
    sd    <- pmap$side[cur]
    pnode <- nodes[[par]]
    var   <- columns[pnode$split_variable]
    op    <- if (identical(sd, "left")) "<=" else ">"
    comps[[length(comps) + 1L]] <- data.frame(
      variable = var, op = op, value = as.numeric(pnode$split_value),
      stringsAsFactors = FALSE)
    cur   <- par
    guard <- guard + 1L
    if (guard > length(nodes)) break  # cycle guard; trees are acyclic
  }
  if (length(comps) == 0L)
    return(data.frame(variable = character(0), op = character(0),
                      value = numeric(0), stringsAsFactors = FALSE))
  out <- do.call(rbind, rev(comps))  # root-to-leaf order
  .grf_simplify_conjunction(out)
}


# ---------------------------------------------------------------------------
# Internal: drop dominated same-variable, same-direction cuts
# ---------------------------------------------------------------------------
# Within one conjunction (AND), repeated cuts on the same variable in the same
# direction collapse to the binding one:
#   multiple "x <= a"  -> smallest a   (tightest upper bound)
#   multiple "x >  b"  -> largest  b   (tightest lower bound)
# Opposite-direction cuts on the same variable (an interval) are both kept.
#' @noRd
.grf_simplify_conjunction <- function(cj) {
  if (nrow(cj) <= 1L) return(cj)
  keep <- rep(TRUE, nrow(cj))
  for (v in unique(cj$variable)) {
    for (o in c("<=", ">")) {
      idx <- which(cj$variable == v & cj$op == o)
      if (length(idx) > 1L) {
        bind <- if (o == "<=") idx[which.min(cj$value[idx])]
                else            idx[which.max(cj$value[idx])]
        drop <- setdiff(idx, bind)
        keep[drop] <- FALSE
      }
    }
  }
  cj[keep, , drop = FALSE]
}


# ---------------------------------------------------------------------------
# Internal: format one conjunction as a brace-label vector
# ---------------------------------------------------------------------------
# Each component becomes "{var op value}".  digits controls value formatting.
#' @noRd
.grf_conjunction_labels <- function(cj, digits = 17L) {
  if (nrow(cj) == 0L) return(character(0))
  sprintf("{%s %s %s}", cj$variable, cj$op,
          formatC(cj$value, format = "g", digits = digits))
}


# ---------------------------------------------------------------------------
# Build the GRF subgroup definition from the tree + target leaves
# ---------------------------------------------------------------------------
#' Construct a correct subgroup definition from a GRF policy tree
#'
#' Returns the exact root-to-leaf decision rule(s) for `target_leaves`,
#' with correct directions, as a structured definition consumed by
#' \code{\link{evaluate_grf_subgroup}} and rendered for display.
#'
#' @param tree A \code{policy_tree} object (has \code{$nodes}, \code{$columns}).
#' @param target_leaves Integer vector of leaf node positions (as returned by
#'   \code{predict(tree, X, type = "node.id")}) that constitute the subgroup.
#' @param digits Significant digits for threshold formatting in labels.
#'
#' @return A list with:
#'   \describe{
#'     \item{conjunctions}{list of data.frames, one per target leaf, each with
#'       columns variable/op/value (simplified).}
#'     \item{labels}{when there is exactly ONE conjunction, a character vector
#'       of \code{"{var op value}"} component labels (back-compatible with
#'       \code{\link{get_dfpred}}); \code{NULL} for a union.}
#'     \item{definition}{a single human-readable string: one
#'       parenthesised conjunction, or several joined by \code{" | "} for a
#'       union.}
#'     \item{is_disjunction}{logical; \code{TRUE} when the subgroup spans more
#'       than one leaf (cannot be a single AND-label vector).}
#'   }
#' @keywords internal
#' @noRd
.grf_build_subgroup_definition <- function(tree, target_leaves, digits = 17L) {
  pmap  <- .grf_parent_map(tree)
  leaves <- sort(unique(as.integer(target_leaves)))
  conjunctions <- lapply(leaves, function(lf)
    .grf_leaf_conjunction(tree, lf, pmap))

  # Drop any empty conjunctions (a degenerate root-only "tree" selecting all).
  nonempty <- vapply(conjunctions, nrow, integer(1)) > 0L
  conjunctions <- conjunctions[nonempty]

  is_disj <- length(conjunctions) > 1L

  if (length(conjunctions) == 0L) {
    return(list(conjunctions = list(), labels = NULL,
                definition = NA_character_, is_disjunction = FALSE))
  }

  conj_strings <- vapply(conjunctions, function(cj) {
    parts <- sprintf("%s %s %s", cj$variable, cj$op,
                     formatC(cj$value, format = "g", digits = digits))
    paste(parts, collapse = " & ")
  }, character(1))

  if (is_disj) {
    definition <- paste0("(", conj_strings, ")", collapse = " | ")
    labels <- NULL
  } else {
    definition <- conj_strings[1]
    labels <- .grf_conjunction_labels(conjunctions[[1]], digits = digits)
  }

  list(conjunctions = conjunctions, labels = labels,
       definition = definition, is_disjunction = is_disj)
}


# ---------------------------------------------------------------------------
# Evaluate a GRF subgroup definition on new data (AND of comps, OR of leaves)
# ---------------------------------------------------------------------------
#' Evaluate a GRF subgroup definition and return membership
#'
#' Computes subgroup membership for an arbitrary data frame from the structured
#' definition built by \code{.grf_build_subgroup_definition()}.  Membership is
#' the union (OR) over conjunctions, each the intersection (AND) over its
#' components -- so it reproduces \code{predict(tree, X)} on the fit data and
#' applies cleanly to \code{df.predict} / \code{df.test}.
#'
#' @param def A definition list from \code{.grf_build_subgroup_definition()}.
#' @param df A data frame containing the split variables.
#'
#' @return Integer vector \code{treat.recommend}: 0 for subgroup members (harm),
#'   1 for the complement.
#' @keywords internal
#' @noRd
.grf_evaluate_subgroup <- function(def, df) {
  n <- nrow(df)
  conjunctions <- def$conjunctions
  if (length(conjunctions) == 0L) {
    # No subgroup -> everyone in complement.
    return(rep(1L, n))
  }
  in_any <- rep(FALSE, n)
  for (cj in conjunctions) {
    in_cj <- rep(TRUE, n)
    for (r in seq_len(nrow(cj))) {
      v  <- cj$variable[r]; op <- cj$op[r]; val <- cj$value[r]
      if (!v %in% names(df))
        stop("Subgroup variable '", v, "' not found in data.", call. = FALSE)
      x <- df[[v]]
      member <- switch(op,
                       "<=" = x <= val,
                       ">"  = x >  val,
                       "<"  = x <  val,
                       ">=" = x >= val,
                       stop("Unsupported operator '", op, "'.", call. = FALSE))
      in_cj <- in_cj & member
    }
    in_any <- in_any | in_cj
  }
  ifelse(in_any, 0L, 1L)
}
