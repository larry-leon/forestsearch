# =============================================================================
# Tier-2 de-biased gate: method-agnostic application layer
# -----------------------------------------------------------------------------
# fs_debias_gate() itself is already method-agnostic: it takes a candidate
# family (named list of row-index vectors), the selected members, a Cox/GLM
# `spec`, and a `reselection` rule, and de-biases the selected subgroup's
# Cox/GLM effect with an IJ interval.  These helpers let the consistency, DINA,
# and GRF branches of forestsearch() all feed it THEIR OWN candidate family
# (approach (i)): each method is de-biased against the family it ranked over,
# with `reselection` matching the method's selection criterion
# (maxeff / maxSG / minSG / effMaxSG / effMinSG).
#
# Candidates are reconstructed uniformly from a table with columns
# (v1, d1, c1, v2, d2, c2) -- variable name, direction ("left"/"right"), cut
# value; slot 2 is NA for singletons -- via the proven .grf_evaluate_subgroup()
# membership engine (0 = in-subgroup/harm, 1 = complement).  DINA's "right"
# means ">=" while GRF's non-"left" means ">", so the operator is parameterised.
# =============================================================================

#' Row-index membership for one candidate conjunction.
#'
#' @param df Analysis data frame.
#' @param cj Conjunction data frame with columns `variable`, `op`, `value`.
#' @return Integer row indices of subjects in the (harm) subgroup.
#' @keywords internal
.fs_dg_members_from_conj <- function(df, cj) {
  def <- list(conjunctions = list(cj), labels = NULL,
              definition = NA_character_, is_disjunction = FALSE)
  which(.grf_evaluate_subgroup(def, df) == 0L)
}

#' Build a gate candidate family (named list of row-index vectors) from a
#' uniform candidate table.
#'
#' @param df Analysis data frame the memberships are evaluated on.
#' @param tab Data frame with columns `v1`, `d1`, `c1`, `v2`, `d2`, `c2`.
#' @param op_right Operator used when direction is not `"left"`: `">="` (DINA)
#'   or `">"` (GRF).
#' @param n_min Minimum subgroup size; smaller candidates are dropped.
#' @param digits Threshold label formatting (matches the package convention).
#' @return Named list of integer row-index vectors (one per candidate).
#' @keywords internal
.fs_dg_family_from_table <- function(df, tab, op_right = ">=",
                                     n_min = 1L, digits = 17L) {
  if (is.null(tab) || !nrow(tab)) return(list())
  fam <- list()
  for (i in seq_len(nrow(tab))) {
    op1 <- if (identical(as.character(tab$d1[i]), "left")) "<=" else op_right
    comps <- list(data.frame(variable = as.character(tab$v1[i]), op = op1,
                             value = as.numeric(tab$c1[i]),
                             stringsAsFactors = FALSE))
    if (!is.na(tab$v2[i])) {
      op2 <- if (identical(as.character(tab$d2[i]), "left")) "<=" else op_right
      comps[[2L]] <- data.frame(variable = as.character(tab$v2[i]), op = op2,
                                value = as.numeric(tab$c2[i]),
                                stringsAsFactors = FALSE)
    }
    cj  <- do.call(rbind, comps)
    mem <- tryCatch(.fs_dg_members_from_conj(df, cj),
                    error = function(e) integer(0))
    if (length(mem) >= n_min) {
      lab <- paste(sprintf("%s %s %s", cj$variable, cj$op,
                           formatC(cj$value, format = "g", digits = digits,
                                   width = 1L)),
                   collapse = " & ")
      fam[[lab]] <- mem
    }
  }
  fam
}

#' Map `sg_focus` to the gate's re-selection rule, faithfully per engine.
#'
#' The gate must re-select under the same rule the search used.  The only rule
#' that differs by engine is the `"hr"`/`"eff"` focus: the consistency search
#' ranks by consistency rate, whose gate analog is `"maxcons"`; DINA/GRF rank by
#' effect, whose analog is `"maxeff"`.  Every other focus maps identically
#' (`maxSG`/`minSG` pass through; `hrMaxSG`/`hrMinSG` -> the size-within-effect-
#' neighborhood rules `effMaxSG`/`effMinSG`).  This makes `sg_focus` the single
#' source of truth, so callers never need to specify `reselection` by hand.
#'
#' @param sg_focus The search focus (`"eff"`, `"maxSG"`, `"effMaxSG"`, ...).
#' @param engine `"consistency"` (so `"hr"` -> `"maxcons"`) or `"effect"`
#'   (DINA/GRF; `"hr"` -> `"maxeff"`).
#' @keywords internal
.fs_dg_reselection_from_focus <- function(sg_focus,
                                          engine = c("effect", "consistency")) {
  engine  <- match.arg(engine)
  hr_rule <- if (identical(engine, "consistency")) "maxcons" else "maxeff"
  sgf <- tryCatch(.normalize_sg_focus(sg_focus), error = function(e) sg_focus)
  switch(as.character(sgf),
         hr       = hr_rule,
         hrMaxSG  = "effMaxSG",
         hrMinSG  = "effMinSG",
         maxSG    = "maxSG",
         minSG    = "minSG",
         maxcons  = "maxcons",
         maxeff   = "maxeff",
         effMaxSG = "effMaxSG",
         effMinSG = "effMinSG",
         hr_rule)
}

#' Apply the Tier-2 de-biased gate with forestsearch's defaulting rules.
#'
#' Thin wrapper around [fs_debias_gate()] shared by the consistency, DINA, and
#' GRF branches so the call (and `debias_gate_args` overrides) stays identical
#' across methods.  Wrapped in tryCatch so a reconstruction failure yields
#' `NULL` (the branch's normal output is unaffected) rather than aborting.
#'
#' @param reselection_default Method-appropriate default re-selection rule
#'   (`"maxcons"` for consistency; [.fs_dg_reselection_from_focus()] otherwise).
#' @keywords internal
.fs_apply_debias_gate <- function(df, candidates, selected_members, spec,
                                  c_screen, c_consistency, p_star,
                                  effect_neighborhood, reselection_default,
                                  selection_rule_default = "neighborhood",
                                  debias_gate_args = list(), seedit = NULL) {
  .g <- function(a, b) if (is.null(a)) b else a
  if (is.null(debias_gate_args)) debias_gate_args <- list()
  tryCatch(
    fs_debias_gate(
      df               = df,
      candidates       = candidates,
      spec             = spec,
      selected_members = selected_members,
      c_screen         = c_screen,
      c_consistency    = c_consistency,
      p_star           = p_star,
      t_gate           = debias_gate_args$t_gate,            # NULL -> near-null
      gate             = .g(debias_gate_args$gate,        "point"),
      reselection      = .g(debias_gate_args$reselection, reselection_default),
      effect_neighborhood = effect_neighborhood,
      selection_rule   = .g(debias_gate_args$selection_rule, selection_rule_default),
      draws            = .g(debias_gate_args$draws,       2000L),
      multiplier       = .g(debias_gate_args$multiplier,  "poisson"),
      include_complement = .g(debias_gate_args$include_complement, TRUE),
      ci_method        = .g(debias_gate_args$ci_method,   "ij"),
      seed             = .g(debias_gate_args$seed,        seedit)),
    error = function(e) {
      warning("debias_gate (", spec$outcome_type, ") failed: ",
              conditionMessage(e), call. = FALSE)
      NULL
    })
}

#' Build the Cox/GLM `spec` for the gate from forestsearch()'s resolved columns.
#'
#' DINA/GRF operate on the user-facing analysis frame `df` (not the standardized
#' `df.fs`), so survival uses the user column names here -- unlike the
#' consistency hook, which uses the internal Y/Event/Treat names on `df.fs`.
#' @keywords internal
.fs_dg_spec <- function(outcome_type, effect_measure, treat.name, outcome.name,
                        event.name, offset.name, adjust_covariates,
                        adverse_outcome, df) {
  adj <- if (length(adjust_covariates))
    intersect(adjust_covariates, names(df)) else NULL
  list(outcome_type    = outcome_type,
       effect_measure  = effect_measure,
       treat.name      = treat.name,
       outcome.name    = outcome.name,
       event.name      = event.name,
       offset.name     = offset.name,
       adjust_covariates = if (length(adj)) adj else NULL,
       adverse_outcome = adverse_outcome)
}
