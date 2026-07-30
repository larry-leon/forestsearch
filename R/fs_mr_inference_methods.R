# =============================================================================
# Tier-2 de-biased gate: method-agnostic application layer
# -----------------------------------------------------------------------------
# fs_mr_inference() itself is already method-agnostic: it takes a candidate
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
.fs_mr_members_from_conj <- function(df, cj) {
  def <- list(conjunctions = list(cj), labels = NULL,
              definition = NA_character_, is_disjunction = FALSE)
  which(.grf_evaluate_subgroup(def, df) == 0L)
}

#' Restrict a candidate table to rows competitive on a method's native ranking
#' statistic (DINA tau-hat, GRF doubly-robust score)
#'
#' The Tier-2 gate re-selects winners on each multiplier draw by perturbed
#' Cox/GLM effect.  DINA and GRF, however, select on their own statistic
#' (subgroup-mean tau-hat / mean DR score), so over a large family the
#' Cox-effect re-selection explores candidates the native rule would never
#' reach, inflating the selection-bias term relative to the full bootstrap
#' (which re-runs the native selection).  Limiting the family to candidates
#' within a band of the best native statistic makes the gate's re-selection
#' neighbourhood mirror the one the bootstrap actually explores.
#'
#' The DINA and GRF gate branches default the band to `effect_neighborhood` (the
#' band those engines select within); a value `>= 1` (or `NULL`) disables the
#' restriction and returns `tab` unchanged.  The consistency method does not call
#' this -- its gate re-selection (`maxcons`) already matches the search statistic.
#'
#' @param tab Candidate table (rows = candidates).
#' @param stat Numeric native statistic aligned to `tab` rows (harm = high).
#' @param neighborhood Multiplicative band in `[0, 1)` on the natural-effect
#'   scale: keep `e >= (1 - neighborhood) * max(e)`, with `e = exp(stat)` when
#'   `log_scale = TRUE`.  `NULL` or a value `>= 1` disables the restriction.
#' @param log_scale Logical; exponentiate `stat` before banding (ratio scales,
#'   e.g. DINA log-HR tau-hat).  GRF DR scores are additive (`FALSE`).
#' @return The retained subset of `tab` (always at least the best row); `tab`
#'   unchanged when disabled, mis-sized, or when no positive harm-side maximum
#'   exists.
#' @keywords internal
#' @noRd
.fs_mr_restrict_native <- function(tab, stat, neighborhood = NULL,
                                   log_scale = FALSE) {
  if (is.null(neighborhood) || !is.finite(neighborhood) ||
      neighborhood < 0 || neighborhood >= 1 ||
      is.null(stat) || is.null(tab) || !nrow(tab)) return(tab)
  s <- suppressWarnings(as.numeric(stat))
  if (length(s) != nrow(tab) || !any(is.finite(s))) return(tab)
  e  <- if (isTRUE(log_scale)) exp(s) else s
  mx <- suppressWarnings(max(e[is.finite(e)]))
  if (!is.finite(mx) || mx <= 0) return(tab)            # no harm-side band
  keep <- is.finite(e) & e >= (1 - neighborhood) * mx
  keep[which.max(replace(e, !is.finite(e), -Inf))] <- TRUE  # always keep best
  out <- tab[keep, , drop = FALSE]
  if (nrow(out)) out else tab
}

#' Build candidate-membership lists from a (v1,d1,c1,v2,d2,c2) table
#'
#' @param df Analysis data frame (for membership evaluation).
#' @param tab Data frame with columns `v1`, `d1`, `c1`, `v2`, `d2`, `c2`.
#' @param op_right Operator used when direction is not `"left"`: `">="` (DINA)
#'   or `">"` (GRF).
#' @param n_min Minimum subgroup size; smaller candidates are dropped.
#' @param digits Threshold label formatting (matches the package convention).
#' @return Named list of integer row-index vectors (one per candidate).
#' @keywords internal
.fs_mr_family_from_table <- function(df, tab, op_right = ">=",
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
    mem <- tryCatch(.fs_mr_members_from_conj(df, cj),
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
.fs_mr_reselection_from_focus <- function(sg_focus,
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
         # maxeffCons -> the gate's "maxeff", which is
         # passers[which.max(beta[passers])]: argmax effect among PASSERS, and
         # `passers` is the gate's consistency-qualifying set (driven by
         # p_star = pconsistency.threshold).  So the gate rule named "maxeff"
         # is the maxeffCons rule, NOT sg_focus = "maxeff" (which is ungated).
         maxeffCons = "maxeff",
         effMaxSG = "effMaxSG",
         effMinSG = "effMinSG",
         hr_rule)
}

#' Apply the Tier-2 de-biased gate with forestsearch's defaulting rules.
#'
#' Thin wrapper around [fs_mr_inference()] shared by the consistency, DINA, and
#' GRF branches so the call (and `mr_inference_args` overrides) stays identical
#' across methods.  Wrapped in tryCatch so a reconstruction failure yields
#' `NULL` (the branch's normal output is unaffected) rather than aborting.
#'
#' @param reselection_default Method-appropriate default re-selection rule
#'   (`"maxcons"` for consistency; [.fs_mr_reselection_from_focus()] otherwise).
#' @keywords internal
.fs_apply_mr <- function(df, candidates, selected_members, spec,
                                  c_screen, c_consistency, p_star,
                                  effect_neighborhood, reselection_default,
                                  selection_rule_default = "neighborhood",
                                  mr_inference_args = list(), seedit = NULL) {
  .g <- function(a, b) if (is.null(a)) b else a
  if (is.null(mr_inference_args)) mr_inference_args <- list()
  tryCatch(
    fs_mr_inference(
      df               = df,
      candidates       = candidates,
      spec             = spec,
      selected_members = selected_members,
      c_screen         = c_screen,
      c_consistency    = c_consistency,
      p_star           = p_star,
      t_confirm        = mr_inference_args$t_confirm,          # NULL -> near-null
      confirm_rule     = .g(mr_inference_args$confirm_rule, "point"),
      reselection      = .g(mr_inference_args$reselection, reselection_default),
      effect_neighborhood = effect_neighborhood,
      selection_rule   = .g(mr_inference_args$selection_rule, selection_rule_default),
      draws            = .g(mr_inference_args$draws,       2000L),
      multiplier       = .g(mr_inference_args$multiplier,  "poisson"),
      include_complement = .g(mr_inference_args$include_complement, TRUE),
      ci_method        = .g(mr_inference_args$ci_method,   "ij"),
      seed             = .g(mr_inference_args$seed,        seedit)),
    error = function(e) {
      warning("mr_inference (", spec$outcome_type, ") failed: ",
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
.fs_mr_spec <- function(outcome_type, effect_measure, treat.name, outcome.name,
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
