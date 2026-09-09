# =============================================================================
# forestsearch_methods.R - S3 Methods for forestsearch Objects
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helpers (not exported)
# -----------------------------------------------------------------------------

#' Safely extract a value from nested list slots
#'
#' Tries multiple accessor paths and returns the first non-NULL value.
#' Eliminates fragile slot-specific code scattered across methods.
#'
#' @param x A forestsearch object.
#' @param ... Character vectors, each a path of nested names to try.
#' @return The first non-NULL value found, or NULL.
#' @keywords internal
#' @noRd
.fs_get <- function(x, ...) {
  paths <- list(...)
  for (path in paths) {
    val <- x
    ok <- TRUE
    for (key in path) {
      if (is.null(val) || !is.list(val) || !key %in% names(val)) {
        ok <- FALSE
        break
      }
      val <- val[[key]]
    }
    if (ok && !is.null(val)) return(val)
  }
  NULL
}


#' Extract human-readable subgroup labels
#'
#' Prioritises \code{grp.consistency$out_sg$sg.harm_label} (readable labels
#' like \code{"er <= 0"}) over bare factor names in \code{sg.harm}.
#'
#' @param x A forestsearch object.
#' @return Character vector of labels, or NULL.
#' @keywords internal
#' @noRd
.fs_sg_labels <- function(x) {
  labels <- .fs_get(
    x,
    c("grp.consistency", "out_sg", "sg.harm_label")
  )
  if (!is.null(labels)) {
    labels <- labels[!is.na(labels) & labels != ""]
    if (length(labels) > 0) return(labels)
  }
  x$sg.harm
}



#' Assemble the certified post-selection products from an MR result
#'
#' Reads \code{x$mr_inference} -- the \code{\link{fs_mr_inference}} return that
#' \code{\link{forestsearch}} attaches under \code{mr_inference = TRUE} -- and
#' returns the quantities \code{print}/\code{summary} report, or \code{NULL}
#' when MR did not run or produced no field block.  Every element is read, none
#' is recomputed.  The complement's certified bound is \strong{field-s}
#' (\code{upper_1s_s}, present under the default
#' \code{field_scale_complement = "selected"}); the unscaled \code{upper_1s} is
#' the documented fallback and is labelled as such.  The joint pair is taken
#' from \code{field$joint_s} when field-s is what is being reported, so the
#' pair and the marginal bound come from the same field.
#'
#' Note on names: the joint element's members are \code{bonf_lower_H} /
#' \code{bonf_upper_Hc} / \code{bonf_gamma} (\code{fs_mr_inference.R}); the
#' \code{bonf_loH} / \code{bonf_upHc} spellings belong to the simulation
#' template's recorder columns, not to the object.
#'
#' @param x A forestsearch object.
#' @return A list, or NULL when there is nothing certified to report.
#' @keywords internal
#' @noRd
.fs_mr_products <- function(x) {
  g <- x$mr_inference
  if (is.null(g) || !is.list(g)) return(NULL)
  f <- g$field
  if (!is.list(f) || !is.null(f$note)) return(NULL)
  fc <- f$complement
  scaled <- is.list(fc) && is.null(fc$note) && !is.null(fc$upper_1s_s)
  comp_up <- if (scaled) fc$upper_1s_s else if (is.list(fc)) fc$upper_1s else NULL
  jt <- if (scaled && is.list(f$joint_s)) f$joint_s else f$joint
  ph <- if (is.list(g$reselection)) g$reselection$p_hat else NULL
  lab <- g$selected_label
  p_hat_H <- if (!is.null(ph) && !is.null(lab) && lab %in% names(ph))
    unname(ph[[lab]]) else NA_real_
  top <- if (!is.null(ph) && length(ph))
    ph[order(-ph)][seq_len(min(3L, length(ph)))] else NULL
  list(
    measure   = g$measure %||% "HR",
    label     = lab,
    harm_lo   = f$lower_1s,
    comp_up   = comp_up,
    comp_tag  = if (scaled) "field-s" else "field (unscaled; field-s absent)",
    joint_lo  = if (is.list(jt)) jt$bonf_lower_H  else NULL,
    joint_up  = if (is.list(jt)) jt$bonf_upper_Hc else NULL,
    joint_gam = if (is.list(jt)) jt$bonf_gamma    else NULL,
    joint_tag = if (scaled && is.list(f$joint_s)) "field-s" else "field",
    ij_lo     = g$debiased$lower,
    ij_hi     = g$debiased$upper,
    se_ij     = g$debiased$se_ij,
    se_field  = f$se_field,
    se_comp_s = if (scaled) fc$se_field_s else if (is.list(fc)) fc$se_field else NULL,
    se_comp_naive = g$complement$debiased$se_wald,
    p_hat_H   = p_hat_H,
    p_hat_sum = if (!is.null(ph)) sum(ph, na.rm = TRUE) else NA_real_,
    p_hat_top = top,
    n_family  = g$n_family,
    ci_method = g$ci_method)
}


#' The caveat lines, worded from dev/notes/NOTE_survival_products_2026-09-09.md
#'
#' Each string is a faithful compression of a sentence of that NOTE and nothing
#' else; no claim here is generated from the analysis at hand.  The p-hat
#' threshold is the NOTE's own "crossing zero near p-hat ~ 0.5" and is
#' \strong{descriptive, not calibrated}: it marks which side of the bias
#' crossing the analysis sits on, not a decision rule.
#' @keywords internal
#' @noRd
.fs_mr_caveats <- function(p) {
  out <- character(0)
  # NOTE, "Two-sided intervals are not certified.": "... its harm-block
  # coverage falls to 0.913-0.917 at 12.4% prevalence with n >= 1000
  # (0.971-0.981 at 31%).  Read two-sided statements at low prevalence and
  # large n with that caveat."
  if (!is.null(p$ij_lo) && is.finite(p$ij_lo))
    out <- c(out, paste(
      "Two-sided intervals are not certified: the IJ two-term interval's",
      "harm-block coverage falls to 0.913-0.917 at 12.4% prevalence with",
      "n >= 1000 (0.971-0.981 at 31%).  Read two-sided statements at low",
      "prevalence and large n with that caveat."))
  # NOTE, "Analysis-time diagnostic.": bias "crossing zero near p-hat ~ 0.5:
  # ... under-correction at high p-hat (the stable-pick regime, +0.02)" and
  # "its p-hat flag is directional, not calibrated".
  if (is.finite(p$p_hat_H) && p$p_hat_H >= 0.5)
    out <- c(out, paste(
      "p-hat(H) >= 0.5 is the stable-pick regime, where the harm-block",
      "correction is under-corrected (+0.02 log units).  The flag is",
      "directional, not calibrated."))
  out
}


#' Wrap and print the caveat lines under a hanging indent
#' @keywords internal
#' @noRd
.fs_cat_caveat <- function(txt, indent = "  ") {
  for (s in txt)
    cat(paste0(indent, strwrap(s, width = 76, prefix = ""), collapse = "\n"), "\n")
}


#' The shared certified-products block used by print() and summary()
#'
#' @param p The list from \code{.fs_mr_products()}.
#' @param long Logical; \code{TRUE} adds the SEs, the re-selection mass and the
#'   certified/not-certified paragraph (the \code{summary} form).
#' @keywords internal
#' @noRd
.fs_print_mr_products <- function(p, long = FALSE) {
  m <- p$measure
  f3 <- function(v) if (is.null(v) || !is.finite(v)) "NA" else formatC(v, format = "f", digits = 3)
  cat("\nPost-selection inference (certified products):\n")
  cat(sprintf("  Harm subgroup H:        one-sided 95%% lower bound on %s   %s\n", m, f3(p$harm_lo)))
  cat(sprintf("  Complement Hc:          one-sided 95%% upper bound on %s   %s   [%s]\n",
              m, f3(p$comp_up), p$comp_tag))
  if (!is.null(p$joint_lo))
    cat(sprintf("  Joint (Bonferroni):     H lower %s, Hc upper %s  (gamma %s each side) [%s]\n",
                f3(p$joint_lo), f3(p$joint_up), f3(p$joint_gam), p$joint_tag))
  cat(sprintf("  Two-sided (IJ, secondary): H (%s, %s)\n", f3(p$ij_lo), f3(p$ij_hi)))
  if (long) {
    cat("\n  Standard errors (log scale):\n")
    cat(sprintf("    field (H)                 se_field    %s\n", f3(p$se_field)))
    cat(sprintf("    field-s (Hc)              se_field_s  %s   (naive complement SE %s)\n",
                f3(p$se_comp_s), f3(p$se_comp_naive)))
    cat(sprintf("    IJ two-term (H)           se_ij       %s\n", f3(p$se_ij)))
  }
  cat(sprintf("  Re-selection frequency  p-hat(H) = %s\n", f3(p$p_hat_H)))
  if (long) {
    if (!is.null(p$p_hat_top)) {
      cat(sprintf("    top-3 re-selection mass:  %s\n",
                  paste(sprintf("%s %s", names(p$p_hat_top),
                                vapply(unname(p$p_hat_top), f3, "")), collapse = " | ")))
    }
    cat(sprintf("    p_hat_sum = %s over a family of %s candidates\n",
                f3(p$p_hat_sum), if (is.null(p$n_family)) "NA" else p$n_family))
    cat("\n  Certified: the one-sided lower bound on H (field) and the one-sided\n")
    cat("  upper bound on Hc (field-s), and the Bonferroni joint pair at gamma =\n")
    cat("  0.025 each side.  Not certified: any two-sided interval.  p-hat(H) is a\n")
    cat("  recorded diagnostic -- no construction reads it.  Read every bound by\n")
    cat("  location against a clinically meaningful effect size, never as\n")
    cat("  significance at the null.  Source: dev/notes/NOTE_survival_products_2026-09-09.md.\n")
  }
  cv <- .fs_mr_caveats(p)
  if (length(cv)) { cat("\n"); .fs_cat_caveat(cv) }
  invisible(NULL)
}

# =============================================================================
# print.forestsearch
# =============================================================================

#' Print Method for forestsearch Objects
#'
#' Displays a concise summary of ForestSearch results including the
#' identified subgroup definition, consistency metrics, algorithm details,
#' and computation time.
#'
#' @section Post-selection inference:
#' When the fit carries multiplier-resampling results (`mr_inference = TRUE`,
#' which attaches the [fs_mr_inference()] return as `x$mr_inference`) and the
#' field block ran, a further block reports the **certified** survival
#' post-selection products of
#' `dev/notes/NOTE_survival_products_2026-09-09.md`, in that document's order
#' of standing: the one-sided 95% **lower** bound on the harm subgroup
#' \eqn{\beta(\widehat H)} (the field) and the one-sided 95% **upper** bound
#' on the complement \eqn{\beta(\widehat H^c)} (**field-s**, the studentized
#' complement field) first; the Bonferroni **joint** pair at
#' \eqn{\gamma = 0.025} each side second; the IJ two-term **two-sided**
#' interval third, explicitly labelled secondary because no two-sided interval
#' is certified; and the re-selection frequency \eqn{\hat p(\widehat H)}
#' last, as a recorded diagnostic that no construction reads.
#'
#' The field and field-s bounds exist only under `ci_method = "field"`, which
#' is the default; `ci_method = "ij"` omits the field block and with it this
#' entire section apart from the two-sided line.  When MR results are absent
#' the printed output is unchanged from a build without this section.
#'
#' Caveat lines are printed only when they apply and are compressions of the
#' NOTE, not of the analysis: a two-sided caveat whenever a two-sided interval
#' is shown, and a stable-pick note when
#' \eqn{\hat p(\widehat H) \ge 0.5}.  That threshold is the NOTE's own
#' description of where the harm-block bias crosses zero and is
#' **descriptive, not calibrated** -- it says which side of the crossing the
#' analysis sits on, it is not a decision rule.  Bounds are reported by
#' location and must be read against a clinically meaningful effect size,
#' never as significance at the null.
#'
#' @param x A \code{forestsearch} object returned by
#'   \code{\link{forestsearch}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \dontrun{
#' fs <- forestsearch(df.analysis = mydata, ...)
#' print(fs)
#' # or simply:
#' fs
#' }
#'
#' @seealso \code{\link{summary.forestsearch}} for detailed output,
#'   \code{\link{plot.forestsearch}} for visualization.
#' @export
print.forestsearch <- function(x, ...) {
  cat("ForestSearch Results\n")
  cat("====================\n\n")

  if (is.null(x$sg.harm)) {
    cat("No subgroup identified.\n")
    return(invisible(x))
  }

  # --- Subgroup definition ---
  labels <- .fs_sg_labels(x)
  cat("Selected Subgroup:\n")
  cat("  Definition:", paste(labels, collapse = " & "), "\n")

  # --- sg_focus (try multiple locations) ---
  sg_focus <- .fs_get(
    x,
    c("grp.consistency", "sg_focus"),
    c("args_call_all", "sg_focus"),
    c("sg_focus")
  )
  if (!is.null(sg_focus)) {
    cat("  sg_focus:", sg_focus, "\n")
  }

  # --- Top-ranked consistency result ---
  top_result <- .fs_get(x, c("grp.consistency", "out_sg", "result"))

  if (!is.null(top_result) && nrow(top_result) > 0) {
    top <- top_result[1, ]
    if ("N" %in% names(top))
      cat("  N:", top$N, "\n")
    if ("hr" %in% names(top))
      cat("  HR:", round(as.numeric(top$hr), 3), "\n")
    if ("Pcons" %in% names(top))
      cat("  Pcons:", round(as.numeric(top$Pcons), 3), "\n")
  }

  # --- Algorithm info ---
  algorithm <- .fs_get(x, c("grp.consistency", "algorithm"))
  if (!is.null(algorithm)) {
    cat("  Algorithm:", algorithm, "\n")
  }

  # --- Candidate-family status (descriptive; see .fs_family_status()) ---
  if (!is.null(x$family_status)) {
    cat("  Candidate family:", x$family_status, "\n")
  }

  # --- Admission set MR re-selects over (see .fs_resolve_admission()) ---
  if (!is.null(x$admission)) {
    cat("  Admission set:", .fs_format_admission(x$admission), "\n")
  }

  n_eval <- .fs_get(x, c("grp.consistency", "n_candidates_evaluated"))
  n_pass <- .fs_get(x, c("grp.consistency", "n_passed"))
  if (!is.null(n_eval)) {
    cat("  Candidates evaluated:", n_eval, "\n")
    cat("  Candidates passed:", n_pass, "\n")
  }

  # --- Post-selection inference (MR / field), when present ---------------
  # Absent MR results (mr_inference = FALSE, or a fit whose field block did
  # not run) leave every line above untouched: .fs_mr_products() returns NULL
  # and nothing is printed, so output is byte-identical to a pre-extension
  # build.  See dev/notes/NOTE_survival_products_2026-09-09.md.
  mrp <- .fs_mr_products(x)
  if (!is.null(mrp)) .fs_print_mr_products(mrp, long = FALSE)

  # --- Timing ---
  if (!is.null(x$minutes_all)) {
    cat("\nComputation time:", round(x$minutes_all, 2), "minutes\n")
  }

  invisible(x)
}


# =============================================================================
# summary.forestsearch
# =============================================================================

#' Summary Method for forestsearch Objects
#'
#' Provides a detailed summary of a ForestSearch analysis including input
#' parameters, variable selection results, consistency evaluation, and
#' the selected subgroup with key metrics.
#'
#' @section Post-selection inference:
#' Carries the same certified-products block as [print.forestsearch()] -- see
#' its documentation for the ordering, the `ci_method = "field"` requirement,
#' the by-location reading convention and the descriptive (not calibrated)
#' \eqn{\hat p} threshold -- and adds, on the log scale, the field SE for the
#' harm block, the field-s SE for the complement with the naive complement SE
#' beside it, and the IJ two-term SE; the three largest re-selection
#' frequencies with their labels; `p_hat_sum` over the re-selection family;
#' and one paragraph naming what is certified and what is not, sourced from
#' `dev/notes/NOTE_survival_products_2026-09-09.md`.  Absent MR results, the
#' output is unchanged from a build without this section.
#'
#' @param object A \code{forestsearch} object returned by
#'   \code{\link{forestsearch}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{object}.
#'
#' @examples
#' \dontrun{
#' fs <- forestsearch(df.analysis = mydata, ...)
#' summary(fs)
#' }
#'
#' @export
summary.forestsearch <- function(object, ...) {
  cat("ForestSearch Summary\n")
  cat("====================\n\n")

  # -------------------------------------------------------------------------
  # Analysis Parameters
  # -------------------------------------------------------------------------
  params <- object$args_call_all
  if (!is.null(params)) {
    cat("Analysis Parameters:\n")
    .print_param <- function(label, value) {
      if (!is.null(value)) cat("  ", label, ": ", value, "\n", sep = "")
    }
    .print_param("sg_focus",                params$sg_focus)
    .print_param("hr.threshold",            params$hr.threshold)
    .print_param("hr.consistency",          params$hr.consistency)
    .print_param("pconsistency.threshold",  params$pconsistency.threshold)
    .print_param("n.min",                   params$n.min)
    .print_param("fs.splits",               params$fs.splits)
    .print_param("maxk",                    params$maxk)
    .print_param("use_twostage",            params$use_twostage)
    .print_param("use_lasso",               params$use_lasso)
    .print_param("use_grf",                 params$use_grf)
    .print_param("use_dina",                params$use_dina)
    # Candidate-family status: "no-front-end" (no fitted model shapes the
    # family on the observed data -- weaker than a Section 2.1 fixed family,
    # which also needs resample-invariant cuts), "conditional-removable" (a
    # front end is on; turning it off reaches "no-front-end"), or
    # "conditional-inherent" (DINA/GRF, family generated by fitting a model).
    .print_param("candidate family",        object$family_status)
    .print_param("admission set",           .fs_format_admission(object$admission))
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Variable Selection
  # -------------------------------------------------------------------------
  n_candidate <- length(object$confounders.candidate)
  n_evaluated <- length(object$confounders.evaluated)
  if (n_candidate > 0 || n_evaluated > 0) {
    cat("Variable Selection:\n")
    cat("  Candidate confounders:", n_candidate, "\n")
    cat("  Confounders evaluated:", n_evaluated, "\n")

    # GRF screening info
    if (!is.null(object$grf_res) && !inherits(object$grf_res, "try-error")) {
      n_grf_cuts <- length(object$grf_cuts)
      if (n_grf_cuts > 0) {
        cat("  GRF cuts applied:", n_grf_cuts, "\n")
      }
    }
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Search Space
  # -------------------------------------------------------------------------
  if (!is.null(object$prop_maxk) || !is.null(object$max_sg_est)) {
    cat("Search Space:\n")
    if (!is.null(object$prop_maxk)) {
      cat("  Proportion of max combinations searched:",
          round(object$prop_maxk, 4), "\n")
    }
    if (!is.null(object$max_sg_est)) {
      eff_lbl <- if (!is.null(object$effect_measure)) {
        object$effect_measure
      } else {
        "HR"
      }
      # For GLM ratio measures, max_sg_est is on log scale
      max_val <- if (!is.null(object$effect_measure) &&
                     object$effect_measure %in% c("OR", "RR", "IRR")) {
        exp(object$max_sg_est)
      } else {
        object$max_sg_est
      }
      cat(sprintf("  Maximum subgroup %s estimate: %.3f\n", eff_lbl, max_val))
    }
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Consistency Evaluation
  # -------------------------------------------------------------------------
  gc <- object$grp.consistency
  if (!is.null(gc)) {
    cat("Consistency Evaluation:\n")
    algorithm <- gc$algorithm %||% "fixed"
    cat("  Algorithm:", algorithm, "\n")

    if (!is.null(gc$n_candidates_evaluated)) {
      cat("  Candidates evaluated:", gc$n_candidates_evaluated, "\n")
      cat("  Candidates passed:", gc$n_passed, "\n")
    }

    # Two-stage specific info
    if (identical(algorithm, "twostage")) {
      ts <- params$twostage_args
      if (!is.null(ts)) {
        cat("  Stage 1 screening splits:", ts$n.splits.screen, "\n")
        cat("  Stage 2 batch size:", ts$batch.size, "\n")
      }
    }
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Selected Subgroup
  # -------------------------------------------------------------------------
  if (!is.null(object$sg.harm)) {
    labels <- .fs_sg_labels(object)
    cat("Selected Subgroup:\n")
    cat("  Definition:", paste(labels, collapse = " & "), "\n")

    # Factor-level names (technical)
    if (!identical(labels, object$sg.harm)) {
      cat("  Factor names:", paste(object$sg.harm, collapse = " & "), "\n")
    }

    top_result <- .fs_get(object, c("grp.consistency", "out_sg", "result"))
    if (!is.null(top_result) && nrow(top_result) > 0) {
      top <- top_result[1, ]
      if ("N" %in% names(top))
        cat("  Sample size:", top$N, "\n")
      if ("hr" %in% names(top))
        cat("  Hazard ratio:", round(as.numeric(top$hr), 3), "\n")
      if ("Pcons" %in% names(top))
        cat("  Consistency:",
            round(as.numeric(top$Pcons) * 100, 1), "%\n")
    }

    # Dataset sizes
    if (!is.null(object$df.est)) {
      n_est <- nrow(object$df.est)
      n_harm <- sum(object$df.est$treat.recommend == 0, na.rm = TRUE)
      cat("  Estimation data: n =", n_est,
          "(harm =", n_harm,
          ", complement =", n_est - n_harm, ")\n")
    }
  } else {
    cat("No subgroup identified.\n")
  }

  # -------------------------------------------------------------------------
  # Post-selection inference (MR / field), when present
  # -------------------------------------------------------------------------
  # As in print(): NULL when MR did not run, so the output above is unchanged.
  mrp <- .fs_mr_products(object)
  if (!is.null(mrp)) .fs_print_mr_products(mrp, long = TRUE)

  # -------------------------------------------------------------------------
  # Timing
  # -------------------------------------------------------------------------
  if (!is.null(object$minutes_all)) {
    cat("\nComputation time:", round(object$minutes_all, 2), "minutes\n")
  }

  invisible(object)
}


# =============================================================================
# plot.forestsearch
# =============================================================================

#' Plot ForestSearch Results
#'
#' Dispatches to \code{\link{plot_sg_results}} for Kaplan-Meier curves,
#' hazard-ratio forest plots, or combined panels.
#'
#' @param x A \code{forestsearch} object returned by
#'   \code{\link{forestsearch}}.
#' @param type Character. Type of plot:
#'   \describe{
#'     \item{\code{"combined"}}{KM curves + forest plot (default)}
#'     \item{\code{"km"}}{Kaplan-Meier survival curves only}
#'     \item{\code{"forest"}}{Hazard-ratio forest plot only}
#'     \item{\code{"summary"}}{Summary statistics panel}
#'   }
#' @param outcome.name Character. Name of time-to-event column.
#'   Default: \code{"Y"}.
#' @param event.name Character. Name of event indicator column.
#'   Default: \code{"Event"}.
#' @param treat.name Character. Name of treatment column.
#'   Default: \code{"Treat"}.
#' @param ... Additional arguments passed to \code{\link{plot_sg_results}},
#'   such as \code{by.risk}, \code{conf.level}, \code{est.scale},
#'   \code{sg0_name}, \code{sg1_name}, \code{treat_labels}, \code{colors},
#'   \code{title}, \code{show_events}, \code{show_ci}, \code{show_logrank},
#'   \code{show_hr}.
#'
#' @return Invisibly returns the plot result from
#'   \code{\link{plot_sg_results}}.
#'
#' @seealso \code{\link{plot_sg_results}} for full control over appearance,
#'   \code{\link{plot_sg_weighted_km}} for weighted KM curves,
#'   \code{\link{plot_subgroup_results_forestplot}} for publication-ready
#'   forest plots.
#'
#' @examples
#' \dontrun{
#' fs <- forestsearch(df.analysis = mydata, ...)
#'
#' # Combined KM + forest plot (default)
#' plot(fs)
#'
#' # KM curves only
#' plot(fs, type = "km")
#'
#' # Forest plot only
#' plot(fs, type = "forest")
#'
#' # With non-standard column names
#' plot(fs, type = "km",
#'      outcome.name = "os_months",
#'      event.name = "os_event",
#'      treat.name = "treatment")
#'
#' # With custom labels
#' plot(fs, sg0_name = "High Risk", sg1_name = "Standard Risk",
#'      treat_labels = c("0" = "Placebo", "1" = "Active Drug"))
#' }
#'
#' @export
plot.forestsearch <- function(x,
                              type = c("combined", "km",
                                       "forest", "summary"),
                              outcome.name = "Y",
                              event.name = "Event",
                              treat.name = "Treat",
                              ...) {

  type <- match.arg(type)

  # No-subgroup guard.  The CONSISTENCY path returns df.est = NULL when no
  # subgroup is identified, but the DINA/GRF selection paths return a
  # populated df.est WITHOUT a treat.recommend column, so keying on df.est
  # alone lets a no-subgroup DINA/GRF fit fall through to plot_sg_results(),
  # which then hard-errors on the missing treat.recommend column instead of
  # the intended graceful no-op.  Key on sg.harm (the actual no-subgroup
  # contract) as well.  A subgroup-found fit has sg.harm non-NULL, so the
  # condition reduces to is.null(x$df.est) -- unchanged for found fits.
  if (is.null(x$sg.harm) || is.null(x$df.est)) {
    message("No subgroup identified -- nothing to plot.")
    return(invisible(x))
  }

  result <- plot_sg_results(
    fs.est       = x,
    plot_type    = type,
    outcome.name = outcome.name,
    event.name   = event.name,
    treat.name   = treat.name,
    ...
  )

  invisible(result)
}
