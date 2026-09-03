# plot_subgroup_sims.R ----------------------------------------------------
# [delivery sentinel: p44r1-3d0fe35f]
# Figure construction for the extreme-subgroups study: a plot() method
# assembling the gg_forest() calls of the committed vignettes (commit
# 002f4f37) parameter-for-parameter, and forest_height() converting
# panel row counts to the fig.height values the vignette chunks use.
#
# Every default below -- axis limits, tick positions, reference-line
# styling, annotation labels and rounding, footnote wording -- is lifted
# from the corresponding vignette chunk for that (metric, panel) pair,
# so the ported documents render identical figures.  All of it is
# overridable through the arguments or `...` (which is merged over the
# defaults before the gg_forest() call).
#
# Phase 4.3 (scale-aware plots): summaries carrying effect metadata
# (Phase 4.2) take a generic branch -- identity scale (xlog = FALSE,
# null line at 0) for MD-style measures, ratio scale (null at 1)
# otherwise, with annotation columns and footnotes built from the
# resolved labels/thresholds and NA-threshold columns dropped.  Legacy
# summaries take the verbatim pre-4.3 branch: every gg_forest()
# argument is byte-identical.
#
# Phase 4.4 axis policy: metadata-carrying results delegate limits and
# ticks to gg_forest()'s data-driven defaults on BOTH scales -- the
# fixed HR constants can hide binary medians (median UB(OR) at small N
# sits at or beyond 9.0, and clip arrows mark whiskers only).
# Metadata-less generic summaries (survival results summarised with
# explicit threshold overrides) keep the HR panel constants, and the
# legacy branch is untouched, so every pre-4.4 path renders unchanged.
# The constants remain one argument away via xlim= / ticks_at=.

#' Recommended figure height for a forest panel
#'
#' Converts a panel's row count to the `fig.height` (inches) used by the
#' vignette chunks: `rows * row_height + overhead`, rounded to 0.1.
#'
#' @param x A `"subgroup_sims_summary"` object.
#' @param panel `"single"`, `"combo"`, or `"highrisk"`.
#' @param row_height Inches per row (default 0.45).
#' @param overhead Fixed inches for title, x-axis, and footnote
#'   (default 1.5).
#'
#' @return A single numeric, rounded to one decimal.
#' @seealso [plot.subgroup_sims_summary()]
#' @export
#' @examples
#' \dontrun{
#' fh_single <- forest_height(S, "single")
#' }
forest_height <- function(x, panel = c("single", "combo", "highrisk"),
                          row_height = 0.45, overhead = 1.5) {
  stopifnot(inherits(x, "subgroup_sims_summary"))
  panel <- match.arg(panel)
  rows <- switch(panel,
                 single   = x$n_single,
                 combo    = x$n_combo,
                 highrisk = x$highrisk$n)
  round(rows * row_height + overhead, 1)
}

# Internal: footnote text per (metric, panel), exactly as in the
# vignette chunks.  `np` (per-trial N) appears only in the hr/single
# footnote; when unknown (NA) that clause is dropped.
.subgroup_sims_footnote <- function(metric, panel, n_sims, hr_true,
                                    np, thr) {
  ht <- sprintf("%.2f", hr_true)
  if (metric == "hr") {
    if (panel == "single" && !is.na(np)) {
      paste0("Median HR (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials, N = ", np, " per trial  |  ",
             "True HR = ", ht, " (red dashed)")
    } else if (panel %in% c("single", "combo")) {
      paste0("Median HR (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials  |  ",
             "True HR = ", ht, " (red dashed)")
    } else {
      paste0("Same subgroups and ordering as the UB(HR) panel above  |  ",
             n_sims, " simulated trials  |  ",
             "Arrowheads mark whiskers extending beyond the axis  |  ",
             "Dashed: true HR = ", ht, ";  Dotted: HR = 1.0")
    }
  } else {
    if (panel %in% c("single", "combo")) {
      paste0("Median UB(HR) (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials  |  ",
             "Dotted line: UB = 1.0;  Dashed line: true HR = ", ht)
    } else {
      paste0("Filter: Pr(UB(HR) >= 2.0) >= ", round(100 * thr), "%  |  ",
             "ITT included for comparison  |  ",
             n_sims, " simulated trials  |  ",
             "Arrowheads mark whiskers extending beyond the axis  |  ",
             "Dotted: UB = 1.0;  Dashed: true HR = ", ht)
    }
  }
}

# Internal: effect-aware footnote text (Phase 4.3), parameterized by the
# resolved labels, the null-line display string ("1.0" ratio / "0"
# identity), and the first UB-threshold string (NA when that tail is
# disabled, which also disables the high-risk filter clause).  The
# legacy builder above is untouched and byte-frozen.
.subgroup_sims_footnote_effect <- function(metric, panel, n_sims, hr_true,
                                           np, thr, est_label, ub_label,
                                           null_str, ub_t1) {
  ht <- sprintf("%.2f", hr_true)
  if (metric == "hr") {
    if (panel == "single" && !is.na(np)) {
      paste0("Median ", est_label,
             " (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials, N = ", np, " per trial  |  ",
             "True ", est_label, " = ", ht, " (red dashed)")
    } else if (panel %in% c("single", "combo")) {
      paste0("Median ", est_label,
             " (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials  |  ",
             "True ", est_label, " = ", ht, " (red dashed)")
    } else {
      paste0("Same subgroups and ordering as the ", ub_label,
             " panel above  |  ",
             n_sims, " simulated trials  |  ",
             "Arrowheads mark whiskers extending beyond the axis  |  ",
             "Dashed: true ", est_label, " = ", ht, ";  Dotted: ",
             est_label, " = ", null_str)
    }
  } else {
    if (panel %in% c("single", "combo")) {
      paste0("Median ", ub_label,
             " (point) + 1st-99th percentile ECI  |  ",
             n_sims, " simulated trials  |  ",
             "Dotted line: UB = ", null_str, ";  Dashed line: true ",
             est_label, " = ", ht)
    } else {
      filt <- if (is.na(ub_t1)) {
        "High-risk filter disabled (no UB threshold); ITT anchor shown"
      } else {
        paste0("Filter: Pr(", ub_label, " >= ", ub_t1, ") >= ",
               round(100 * thr), "%  |  ITT included for comparison")
      }
      paste0(filt, "  |  ",
             n_sims, " simulated trials  |  ",
             "Arrowheads mark whiskers extending beyond the axis  |  ",
             "Dotted: UB = ", null_str, ";  Dashed: true ",
             est_label, " = ", ht)
    }
  }
}

#' Forest plot of an extreme-subgroups summary
#'
#' Draws one of the vignettes' six forest panels from a
#' [summary.subgroup_sims()] object via [gg_forest()]: the empirical-CI
#' distribution of the HR point estimate or of the upper 95% bound
#' UB(HR), for the single-variable panel, the combination + random-
#' benchmark panel, or the high-risk panel (Pr(UB>=2) filter with the
#' ITT row anchored first).
#'
#' Defaults reproduce the committed vignette figures exactly:
#' metric-specific axis limits, tick positions, reference/vertical line
#' styling, the four annotation columns (N plus tail probabilities and
#' the unconditional median), and the panel-specific footnotes. Any
#' default can be overridden through the named arguments or `...`,
#' which is merged over the assembled argument list before the
#' [gg_forest()] call (so e.g. `ref_col = "black"` or
#' `widths = ...` in `...` win).
#'
#' Summaries carrying effect metadata (the generic path of
#' [summary.subgroup_sims()], e.g. from [subgroup_glm()] fits) take a
#' scale-aware branch instead: identity-scale measures (MD) plot with
#' `xlog = FALSE` and the null line at 0, ratio-scale measures (OR)
#' with `xlog = TRUE` and the null line at 1, and both delegate
#' limits/ticks to [gg_forest()]'s data-driven defaults (the HR panel
#' constants can hide binary medians and remain available via
#' `xlim=`/`ticks_at=`); annotation
#' columns and footnotes are built from the resolved labels and
#' threshold strings (legacy literals preserved when a pair equals the
#' legacy pair), and columns for disabled (`NA`) thresholds are dropped
#' (the default `widths` shrink to match; an explicit `widths` is used
#' as given). `hr_true` is interpreted on the fitter's estimate scale
#' throughout. Legacy summaries produce byte-identical [gg_forest()]
#' calls.
#'
#' @param x A `"subgroup_sims_summary"` object.
#' @param metric `"hr"` (HR point estimate) or `"ub"` (upper 95% bound).
#' @param panel `"single"`, `"combo"`, or `"highrisk"`.
#' @param hr_true True hazard ratio for the reference/vertical line and
#'   footnotes; defaults to the value stored on `x` and must be
#'   non-`NULL`.
#' @param base_size,point_size,line_size Passed to [gg_forest()];
#'   defaults match the vignettes' `fp_*` values.
#' @param widths Column width vector passed to [gg_forest()].
#' @param clip_marker Passed to [gg_forest()]; `"arrow"` marks whiskers
#'   clipped at the axis limits.
#' @param xlim,ticks_at,xlab,footnote `NULL` (default) selects the
#'   vignette value for the chosen `metric`/`panel`.
#' @param ... Further arguments merged over the defaults and passed to
#'   [gg_forest()].
#'
#' @return The [gg_forest()] plot object.
#' @seealso [forest_height()], [summary.subgroup_sims()]
#' @export
#' @examples
#' \dontrun{
#' S <- summary(sims_uniform, hr_true = 0.70)
#' plot(S, metric = "ub", panel = "combo")
#' }
plot.subgroup_sims_summary <- function(x,
                                       metric = c("hr", "ub"),
                                       panel = c("single", "combo",
                                                 "highrisk"),
                                       hr_true = x$hr_true,
                                       base_size = 14,
                                       point_size = 2.2,
                                       line_size = 0.75,
                                       widths = c(3.5, 5, 1.1, 1.2, 1.2,
                                                  1.1),
                                       clip_marker = "arrow",
                                       xlim = NULL,
                                       ticks_at = NULL,
                                       xlab = NULL,
                                       footnote = NULL,
                                       ...) {
  metric <- match.arg(metric)
  panel  <- match.arg(panel)
  if (is.null(hr_true)) {
    stop("`hr_true` is required for the reference lines and footnotes; ",
         "pass it here or to summary().")
  }

  q <- if (metric == "hr") x$hr_q else x$ub_q

  # Row selection per panel, exactly as in the vignette chunks
  if (panel == "highrisk") {
    sel <- x$highrisk$include
    ord <- x$highrisk$ord
  } else {
    sel <- if (metric == "hr") x$ok else x$ok_ub
    ord <- if (metric == "hr") {
      if (panel == "single") x$idx_hr_single else x$idx_hr_combo
    } else {
      if (panel == "single") x$idx_ub_single else x$idx_ub_combo
    }
  }
  pick <- function(v) v[sel][ord]

  # Effect-aware branch selector (Phase 4.3): summaries from Phase 4.2's
  # generic path carry `labels`/`thresholds` (and usually `effect`).
  # Legacy summaries take the verbatim pre-4.3 assembly below.
  generic <- !is.null(x$labels)

  if (!generic) {
  # Annotation block per metric (labels, rounding as in the chunks)
  annot <- if (metric == "hr") {
    list(
      "N"      = as.character(pick(x$mean_n)),
      "HR<0.5" = paste0(round(100 * pick(x$pr_hr_lt050), 1), "%"),
      "HR>1.0" = paste0(round(100 * pick(x$pr_hr_gt1),   1), "%"),
      "mHR"    = pick(x$mhr_uncond)
    )
  } else {
    list(
      "N"     = as.character(pick(x$mean_n)),
      "UB>=2" = paste0(round(100 * pick(x$pr_ub_ge2), 1), "%"),
      "UB>=3" = paste0(round(100 * pick(x$pr_ub_ge3), 1), "%"),
      "mUB"   = pick(x$mub_uncond)
    )
  }

  args <- list(
    subgroups   = pick(x$sg_names),
    est         = q[sel, 2][ord],
    lo          = q[sel, 1][ord],
    hi          = q[sel, 3][ord],
    cat_vec     = pick(x$cat),
    cat_colours = x$cat_cols,
    annot       = annot,
    point_size  = point_size,
    line_size   = line_size,
    base_size   = base_size,
    widths      = widths,
    clip_marker = clip_marker,
    xlog        = TRUE
  )

  if (metric == "hr") {
    args$ref_line   <- hr_true
    args$vert_lines <- 1.00
    args$ref_col    <- "firebrick"
    args$vert_col   <- "grey55"
    args$vert_lty   <- "dotted"
    if (panel == "highrisk") args$ref_lty <- "dashed"
    args$xlim     <- if (is.null(xlim)) c(0.15, 3.5) else xlim
    args$ticks_at <- if (is.null(ticks_at)) {
      c(0.20, 0.35, 0.50, 0.70, 1.00, 1.50, 2.50)
    } else ticks_at
    args$xlab <- if (is.null(xlab)) {
      "Hazard Ratio (Cox, stratified by grade)"
    } else xlab
  } else {
    args$ref_line   <- 1.00
    args$vert_lines <- hr_true
    args$ref_col    <- "grey55"
    args$ref_lty    <- "dotted"
    args$vert_col   <- "firebrick"
    args$vert_lty   <- "dashed"
    args$xlim     <- if (is.null(xlim)) c(0.30, 9.0) else xlim
    args$ticks_at <- if (is.null(ticks_at)) {
      c(0.40, 0.70, 1.00, 1.50, 2.50, 4.00, 8.00)
    } else ticks_at
    args$xlab <- if (is.null(xlab)) {
      "Upper Bound of 95% CI  [UB(HR)]"
    } else xlab
  }

  args$footnote <- if (is.null(footnote)) {
    .subgroup_sims_footnote(metric, panel, x$n_sims, hr_true,
                            np  = if (is.null(x$n_per_trial)) NA_integer_
                                  else x$n_per_trial,
                            thr = x$highrisk$threshold)
  } else footnote

  } else {
  # ----- generic (effect-aware) assembly ---------------------------------
  eff        <- x$effect
  log_scale  <- if (!is.null(eff)) isTRUE(eff$log_scale) else TRUE
  null_value <- if (!is.null(eff)) eff$null_value else 1
  est_label  <- x$labels$est
  ub_label   <- x$labels$ub
  tE <- x$thresholds$est
  tU <- x$thresholds$ub
  thr_strs <- function(t, legacy_vals, legacy_strs, fallbacks) {
    if (identical(as.numeric(t), legacy_vals)) return(legacy_strs)
    vapply(seq_along(t), function(i) {
      if (is.na(t[i])) fallbacks[i] else
        format(t[i], trim = TRUE, drop0trailing = TRUE)
    }, character(1L))
  }
  est_s <- thr_strs(tE, c(0.5, 1.0), c("0.5", "1.0"), c("lo", "hi"))
  ub_s  <- thr_strs(tU, c(2, 3),     c("2", "3"),     c("t1", "t2"))
  fmt_col <- function(v) ifelse(is.na(v), "-",
                                paste0(round(100 * v, 1), "%"))

  # Annotation columns: NA-threshold columns are dropped entirely.
  annot <- list("N" = as.character(pick(x$mean_n)))
  if (metric == "hr") {
    if (!is.na(tE[1])) {
      annot[[sprintf("%s<%s", est_label, est_s[1])]] <-
        fmt_col(pick(x$pr_hr_lt050))
    }
    if (!is.na(tE[2])) {
      annot[[sprintf("%s>%s", est_label, est_s[2])]] <-
        fmt_col(pick(x$pr_hr_gt1))
    }
    annot[[paste0("m", est_label)]] <- pick(x$mhr_uncond)
  } else {
    if (!is.na(tU[1])) {
      annot[[sprintf("UB>=%s", ub_s[1])]] <- fmt_col(pick(x$pr_ub_ge2))
    }
    if (!is.na(tU[2])) {
      annot[[sprintf("UB>=%s", ub_s[2])]] <- fmt_col(pick(x$pr_ub_ge3))
    }
    annot[["mUB"]] <- pick(x$mub_uncond)
  }
  w <- widths
  if (missing(widths)) {
    k <- length(annot)
    w <- c(3.5, 5, 1.1, rep(1.2, max(0L, k - 2L)),
           if (k >= 2L) 1.1)
  }

  args <- list(
    subgroups   = pick(x$sg_names),
    est         = q[sel, 2][ord],
    lo          = q[sel, 1][ord],
    hi          = q[sel, 3][ord],
    cat_vec     = pick(x$cat),
    cat_colours = x$cat_cols,
    annot       = annot,
    point_size  = point_size,
    line_size   = line_size,
    base_size   = base_size,
    widths      = w,
    clip_marker = clip_marker,
    xlog        = log_scale
  )
  null_str <- if (log_scale) "1.0" else "0"

  if (metric == "hr") {
    args$ref_line   <- hr_true
    args$vert_lines <- null_value
    args$ref_col    <- "firebrick"
    args$vert_col   <- "grey55"
    args$vert_lty   <- "dotted"
    if (panel == "highrisk") args$ref_lty <- "dashed"
    # Axis policy (Phase 4.4): metadata-carrying results delegate
    # limits/ticks to gg_forest()'s data-driven defaults (NULL) on both
    # scales; metadata-less ratio summaries (survival with explicit
    # threshold overrides) keep the HR panel constants. See the file
    # header for the rationale and the reachable-path table.
    args$xlim <- if (is.null(xlim)) {
      if (log_scale && is.null(eff)) c(0.15, 3.5) else NULL
    } else xlim
    args$ticks_at <- if (is.null(ticks_at)) {
      if (log_scale && is.null(eff))
        c(0.20, 0.35, 0.50, 0.70, 1.00, 1.50, 2.50) else NULL
    } else ticks_at
    args$xlab <- if (is.null(xlab)) est_label else xlab
  } else {
    args$ref_line   <- null_value
    args$vert_lines <- hr_true
    args$ref_col    <- "grey55"
    args$ref_lty    <- "dotted"
    args$vert_col   <- "firebrick"
    args$vert_lty   <- "dashed"
    args$xlim <- if (is.null(xlim)) {
      if (log_scale && is.null(eff)) c(0.30, 9.0) else NULL
    } else xlim
    args$ticks_at <- if (is.null(ticks_at)) {
      if (log_scale && is.null(eff))
        c(0.40, 0.70, 1.00, 1.50, 2.50, 4.00, 8.00) else NULL
    } else ticks_at
    args$xlab <- if (is.null(xlab)) ub_label else xlab
  }

  args$footnote <- if (is.null(footnote)) {
    .subgroup_sims_footnote_effect(
      metric, panel, x$n_sims, hr_true,
      np  = if (is.null(x$n_per_trial)) NA_integer_ else x$n_per_trial,
      thr = x$highrisk$threshold,
      est_label = est_label, ub_label = ub_label,
      null_str = null_str,
      ub_t1 = if (is.na(tU[1])) NA_character_ else ub_s[1])
  } else footnote
  }

  dots <- list(...)
  if (length(dots)) args[names(dots)] <- dots

  do.call(gg_forest, args)
}
