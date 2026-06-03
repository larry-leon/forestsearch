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


# =============================================================================
# GRF doubly-robust-score Pareto frontier (selection = "frontier")
# =============================================================================
# An alternative to the policy tree that consumes the SAME doubly-robust score
# matrix.  It enumerates single-covariate thresholds (depth 1) and, when
# grf_depth >= 2, covariate-pair conjunctions (depth 2), scores each candidate
# subgroup S by its DR harm-effect (mean(control[S]) - mean(treated[S])) and
# size, takes the Pareto frontier, and selects one subgroup under a rule
# (effMaxSG / eff / maxSG).  The selected subgroup is always a single
# conjunction, so it maps directly onto the standard sg_def structure.
#
# Sign convention matches the tree path: dr columns are per-action value
# estimates (higher = better outcome); diff = control - treated > 0 means
# treatment harms S.  dmin is the harm-eligibility threshold (== dmin.grf).
# =============================================================================

#' Enumerate single-covariate DR-score candidates
#' @noRd
.grf_dr_candidates <- function(X, dr, grid_probs = seq(0.1, 0.9, 0.1),
                               n_min = 60L) {
  ctrl <- dr[, 1L]; trt <- dr[, 2L]; n <- nrow(X); cn <- colnames(X)
  out <- list(); k <- 0L
  for (j in seq_len(ncol(X))) {
    xj <- X[, j]
    cuts <- unique(stats::quantile(xj, probs = grid_probs, names = FALSE,
                                   type = 7))
    for (cc in cuts) for (dir in c("left", "right")) {
      S <- if (dir == "left") xj <= cc else xj > cc
      nS <- sum(S)
      if (nS < n_min || nS > n - 1L) next
      k <- k + 1L
      out[[k]] <- data.frame(
        v1 = cn[j], d1 = dir, c1 = as.numeric(cc),
        v2 = NA_character_, d2 = NA_character_, c2 = NA_real_,
        effect = mean(ctrl[S]) - mean(trt[S]), size = nS,
        stringsAsFactors = FALSE)
    }
  }
  if (!k) return(NULL)
  do.call(rbind, out)
}

#' Enumerate depth-2 covariate-pair DR-score candidates
#' @noRd
.grf_dr_candidates_d2 <- function(X, dr, grid_probs = seq(0.2, 0.8, 0.2),
                                  n_min = 60L) {
  ctrl <- dr[, 1L]; trt <- dr[, 2L]; n <- nrow(X); p <- ncol(X); cn <- colnames(X)
  conds <- vector("list", p)
  for (j in seq_len(p)) {
    xj <- X[, j]
    cuts <- unique(stats::quantile(xj, probs = grid_probs, names = FALSE,
                                   type = 7))
    cl <- list()
    for (cc in cuts) for (dir in c("left", "right"))
      cl[[length(cl) + 1L]] <- list(cov = cn[j], dir = dir, cut = cc,
                                    v = if (dir == "left") xj <= cc else xj > cc)
    conds[[j]] <- cl
  }
  out <- list(); k <- 0L
  if (p >= 2L) for (j in seq_len(p - 1L)) for (jj in (j + 1L):p) {
    for (a in conds[[j]]) for (b in conds[[jj]]) {
      S <- a$v & b$v; nS <- sum(S)
      if (nS < n_min || nS > n - 1L) next
      k <- k + 1L
      out[[k]] <- data.frame(
        v1 = a$cov, d1 = a$dir, c1 = as.numeric(a$cut),
        v2 = b$cov, d2 = b$dir, c2 = as.numeric(b$cut),
        effect = mean(ctrl[S]) - mean(trt[S]), size = nS,
        stringsAsFactors = FALSE)
    }
  }
  if (!k) return(NULL)
  do.call(rbind, out)
}

#' Mark the Pareto frontier on (effect maximize, size maximize)
#' @noRd
.grf_mark_frontier <- function(cand) {
  ord <- order(-cand$effect, -cand$size)
  cand <- cand[ord, , drop = FALSE]
  best_size <- -Inf; on <- logical(nrow(cand))
  for (i in seq_len(nrow(cand))) {
    if (cand$size[i] > best_size) { on[i] <- TRUE; best_size <- cand$size[i] }
  }
  cand$on_frontier <- on
  cand
}

#' Select one frontier subgroup under a rule (harm side: effect >= dmin)
#'
#' @param rule one of \code{"effMaxSG"} (default), \code{"eff"}, \code{"maxSG"},
#'   \code{"minSG"}, \code{"effMinSG"}; semantics match \code{dina_subgroup()}'s
#'   \code{sg_focus}.  \code{"eff"} = max harm-effect; \code{"maxSG"}/\code{"minSG"}
#'   = largest/smallest eligible (effect >= dmin); \code{"effMaxSG"}/
#'   \code{"effMinSG"} = largest/smallest within a relative neighborhood of the
#'   max harm-effect.
#' @param nbhd relative neighborhood for the \code{eff*SG} rules: keep candidates
#'   with effect >= (1 - nbhd) * max-eligible-effect, then take the
#'   largest/smallest.
#' @return one-row data.frame (the selected candidate), or NULL.
#' @noRd
.grf_frontier_select <- function(cand, dmin, rule = "effMaxSG", nbhd = 0.10) {
  if (is.null(cand) || !nrow(cand)) return(NULL)
  elig <- cand[cand$effect >= dmin, , drop = FALSE]
  if (!nrow(elig)) return(NULL)
  emax <- max(elig$effect)
  if (rule == "maxSG") {
    elig[which.max(elig$size), , drop = FALSE]
  } else if (rule == "minSG") {
    elig[which.min(elig$size), , drop = FALSE]
  } else if (rule == "eff") {
    elig[which.max(elig$effect), , drop = FALSE]
  } else { # effMaxSG / effMinSG: restrict to the near-max-effect band first
    band <- elig[elig$effect >= emax * (1 - nbhd), , drop = FALSE]
    if (rule == "effMinSG") {
      band[which.min(band$size), , drop = FALSE]
    } else { # effMaxSG (default)
      band[which.max(band$size), , drop = FALSE]
    }
  }
}

#' Build the standard sg_def structure from a selected frontier candidate
#'
#' A frontier subgroup is a single conjunction (one or two cuts), so this
#' returns the same shape as \code{.grf_build_subgroup_definition()} for a
#' one-leaf subgroup: a labels vector, a definition string, conjunctions list,
#' and is_disjunction = FALSE.
#' @noRd
.grf_sg_def_from_candidate <- function(cand_row, digits = 17L) {
  comps <- list(data.frame(
    variable = cand_row$v1,
    op       = if (cand_row$d1 == "left") "<=" else ">",
    value    = as.numeric(cand_row$c1), stringsAsFactors = FALSE))
  if (!is.na(cand_row$v2)) {
    comps[[2L]] <- data.frame(
      variable = cand_row$v2,
      op       = if (cand_row$d2 == "left") "<=" else ">",
      value    = as.numeric(cand_row$c2), stringsAsFactors = FALSE)
  }
  cj <- do.call(rbind, comps)
  cj <- .grf_simplify_conjunction(cj)
  parts <- sprintf("%s %s %s", cj$variable, cj$op,
                   formatC(cj$value, format = "g", digits = digits))
  list(conjunctions   = list(cj),
       labels         = .grf_conjunction_labels(cj, digits = digits),
       definition     = paste(parts, collapse = " & "),
       is_disjunction = FALSE)
}
