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
    # Membership-agreement diagnostics, present only under
    # fs_mr_inference(field_recovery = TRUE) (TASK_field_recovery_2026-09-09).
    # NULL when absent, and summary() prints the line only when it is not.
    recov_sens_H = if (is.list(f$recovery)) f$recovery$sens_H else NULL,
    recov_n_used = if (is.list(f$recovery)) f$recovery$n_used else NULL,
    n_family  = g$n_family,
    ci_method = g$ci_method)
}


#' The fixed note printed beneath the post-selection block
#'
#' One note, printed whenever the block prints (print() and summary()); it
#' does not depend on the analysis at hand.  \code{p} is accepted for
#' interface stability and is not read.
#' @keywords internal
#' @noRd
.fs_mr_caveats <- function(p) {
  c("One-sided bounds are the primary products; the two-sided IJ interval is",
    "secondary. p-hat(H) is a descriptive diagnostic. Read each bound by its",
    "location relative to a clinically meaningful effect size. Methods and",
    "simulation operating characteristics: Leon & Anderson (2026),",
    "arXiv:2609.38361.")
}


#' Wrap and print the caveat lines under a hanging indent
#' @keywords internal
#' @noRd
.fs_cat_caveat <- function(txt, indent = "  ") {
  for (s in txt)
    cat(paste0(indent, strwrap(s, width = 76, prefix = ""), collapse = "\n"), "\n")
}


#' The shared post-selection block used by print() and summary()
#'
#' @param p The list from \code{.fs_mr_products()}.
#' @param long Logical; \code{TRUE} adds the SEs and the re-selection mass
#'   (the \code{summary} form).
#' @keywords internal
#' @noRd
.fs_print_mr_products <- function(p, long = FALSE) {
  m <- p$measure
  f3 <- function(v) if (is.null(v) || !is.finite(v)) "NA" else formatC(v, format = "f", digits = 3)
  cat("\nPost-selection inference (selection-adjusted; Leon & Anderson 2026):\n")
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
    # Membership agreement beside p-hat: p-hat counts EXACT re-selections of
    # Hhat, sens_H measures how much of Hhat the re-selections retain, so a low
    # p-hat with a high sens_H says the re-selections were near-twins rather
    # than unrelated regions.  Printed only when field_recovery = TRUE put it
    # there; summary() only, and no construction reads it.
    if (!is.null(p$recov_sens_H) && is.finite(p$recov_sens_H))
      cat(sprintf(
        "    membership recovery: sens_H = %s (mean share of Hhat retained over %s draws)\n",
        f3(p$recov_sens_H),
        if (is.null(p$recov_n_used)) "NA" else format(p$recov_n_used)))
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
#' attached as `x$mr_inference`) and the field block ran
#' (`ci_method = "field"`, the default), a further block reports: the
#' one-sided 95% lower bound for the identified subgroup (the field); the
#' one-sided 95% upper bound for its complement (field-s, the studentized
#' complement field); the Bonferroni pair for a joint two-subgroup
#' statement; the infinitesimal-jackknife (IJ) two-term two-sided interval,
#' as a secondary summary; and the re-selection frequency
#' \eqn{\hat p(\widehat H)} as a descriptive diagnostic.  Bounds are read by
#' their location relative to a clinically meaningful effect size, not as
#' tests at the null.  With `ci_method = "ij"` only the two-sided interval is
#' shown; without MR results the output is unchanged.  The constructions and
#' their simulation operating characteristics are described in León and
#' Anderson (2026), [arXiv:2609.38361](https://arxiv.org/abs/2609.38361).
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
#' Carries the same block as [print.forestsearch()] (see there) and adds the
#' log-scale standard errors (field, field-s with the naive complement SE
#' beside it, and IJ two-term), the three largest re-selection frequencies
#' with their labels, and `p_hat_sum` over the re-selection family.  Absent MR
#' results, the output is unchanged.
#'
#' When the fit carries the membership-agreement diagnostics -- that is, under
#' `mr_inference_args = list(field_recovery = TRUE)` -- one further line
#' reports `sens_H`, the mean share of the identified subgroup that the field's
#' re-selections retain, beside \eqn{\hat p}.  The pair separates two very
#' different low-\eqn{\hat p} situations: \eqn{\hat p} counts *exact*
#' re-selections of \eqn{\hat H}, so a low \eqn{\hat p} with a high `sens_H`
#' says the draws re-selected near-twins of \eqn{\hat H}, while a low
#' \eqn{\hat p} with a low `sens_H` says they re-selected unrelated regions.
#' Like \eqn{\hat p} it is descriptive and no construction reads it, and it
#' answers a narrower question than a bootstrap or cross-validation recovery
#' rate: re-selection within the fixed candidate family under perturbation, not
#' re-discovery from resampled data.  The line is absent -- and the rest of the
#' output byte-identical -- when the diagnostics were not requested.
#' [print.forestsearch()] does not report it.
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
    # c2 / p* are echoed whatever the path selected on; annotate them where
    # the path returns before the consistency stage and never reads them.
    .c2_note <- .fs_c2_inert_note(params$subgroup_method)
    .print_param("sg_focus",                params$sg_focus)
    .print_param("hr.threshold",            params$hr.threshold)
    .print_param("hr.consistency",
                 if (is.null(params$hr.consistency)) NULL
                 else paste0(params$hr.consistency, .c2_note))
    .print_param("pconsistency.threshold",
                 if (is.null(params$pconsistency.threshold)) NULL
                 else paste0(params$pconsistency.threshold, .c2_note))
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
