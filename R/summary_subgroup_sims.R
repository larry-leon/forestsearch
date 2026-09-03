# summary_subgroup_sims.R -------------------------------------------------
# Single source of truth for the extreme-subgroups summary statistics.
#
# Definitions are lifted verbatim from the committed vignettes (commit
# 002f4f37): chunk `uniform-results-table` (quantiles, tail
# probabilities, convergence counts, the formatted results table) plus
# the panel machinery of `forest-theme-params` and
# `ub-forest-highrisk-setup` (category mapping, colour map, ok masks,
# single/combo index vectors, high-risk inclusion and ordering).  The
# design-comparison memo's per-design statistics are these same
# definitions; consolidating them here removes the three-way duplication
# that previously had to be kept in sync by hand.
#
# Phase 4.2 (effect-aware summaries): tail thresholds and display labels
# are now parameterized, resolving explicit arguments -> the result's
# `effect` metadata (stamped by run_subgroup_sims() from subgroup_glm())
# -> the HR legacy literals.  A survival/legacy call with no overrides
# takes a short-circuit branch whose expressions are byte-for-byte the
# pre-4.2 code, so legacy summaries -- fields, printed table, everything
# -- are identical (dev/accept_phase42_summary_compare.R).

#' Summarize an extreme-subgroups simulation study
#'
#' Computes, from the raw matrices of a [run_subgroup_sims()] result (or
#' a loaded RDS payload coerced with `class<-`), every statistic used by
#' the vignettes' results table and forest panels: 1st/50th/99th
#' percentile matrices for the estimate and its upper bound, convergence
#' counts, mean subgroup N, tail probabilities, unconditional median
#' display strings, the formatted `results_tbl`, category labels and
#' colours, validity masks, the single/combination panel index vectors,
#' and the high-risk panel (threshold filter, ITT anchor, ordering).
#'
#' Field names are structural and retained across outcome types
#' (`hr_q`, `pr_hr_lt050`, `pr_ub_ge2`, ...), exactly as `sim_hrs` and
#' `hr.threshold` serve generic duty elsewhere in the package: on a GLM
#' result they hold the fitter's estimate-scale statistics. Estimates
#' follow the estimator layer's harm-positive normalization, so the
#' upper tail (`Pr(est > null)`, `Pr(UB >= t)`) always points toward
#' harm whatever the raw outcome's direction.
#'
#' Threshold and label resolution: explicit arguments win; otherwise the
#' result's `effect` metadata supplies them (for [subgroup_glm()] fits:
#' `est_thresholds = c(NA, null_value)`, `ub_thresholds = c(NA, NA)`,
#' labels from the effect measure); otherwise the HR legacy values
#' `c(0.5, 1.0)` / `c(2, 3)` with the vignettes' literal headers. An
#' `NA` threshold degrades gracefully: its probability field is `NA`
#' and its table column renders `"-"` under a placeholder header.
#'
#' Structurally empty subgroups (all-`NA` columns) yield `NA` medians
#' and `NaN` means/proportions, exactly as in the vignettes; the
#' formatted table renders them as `"-"`.
#'
#' @param object A `"subgroup_sims"` object.
#' @param hr_true True effect on the fitter's estimate scale, for
#'   reference lines (a hazard ratio for Cox fits, a mean difference
#'   for [subgroup_glm()] fits); defaults to the value stored on
#'   `object`.
#' @param probs Quantile probabilities for the ECI matrices; the default
#'   `c(0.01, 0.50, 0.99)` matches the vignettes and downstream code
#'   indexes column 2 as the median.
#' @param est_thresholds Length-2 numeric `c(low, high)`: the estimate
#'   tails are `Pr(est < low)` and `Pr(est > high)`. `NULL` resolves
#'   via the `effect` metadata, then the HR legacy `c(0.5, 1.0)`. `NA`
#'   entries disable that tail (see above).
#' @param ub_thresholds Length-2 numeric: the upper-bound tails are
#'   `Pr(UB >= ub_thresholds[1])` and `Pr(UB >= ub_thresholds[2])`;
#'   the first also drives the high-risk panel. Resolution as for
#'   `est_thresholds` (HR legacy `c(2, 3)`).
#' @param est_label,ub_label Display labels for the estimate and its
#'   upper bound in the table headers (legacy `"HR"` / `"UB"`).
#' @param single_cats,combo_cats Category labels routed to the
#'   single-variable and combination forest panels.
#' @param pr_thresh_highrisk High-risk panel inclusion threshold on
#'   the first upper-bound tail probability, default `0.10`.
#' @param itt_name Name of the ITT row, always included in and anchored
#'   at the top of the high-risk panel.
#' @param cat_cols Named colour map for categories.
#' @param ... Unused.
#'
#' @return An object of class `"subgroup_sims_summary"`: a list with
#'   fields named as in the vignettes (`hr_q`, `ub_q`, `n_valid`,
#'   `pct_valid`, `mean_n`, `pr_hr_lt050`, `pr_hr_gt1`, `pr_ub_ge2`,
#'   `pr_ub_ge3`, `mhr_uncond`, `mub_uncond`, `results_tbl`,
#'   `sg_names`, `cat`, `cat_cols`, `ok`, `ok_ub`, `idx_hr_single`,
#'   `idx_hr_combo`, `idx_ub_single`, `idx_ub_combo`, `n_single`,
#'   `n_combo`, `highrisk`, `n_sims`, `design`, `hr_true`). Legacy
#'   calls return exactly this set; effect-aware calls append `effect`
#'   (the result's metadata, possibly `NULL` when only explicit
#'   overrides were given), `thresholds` (`list(est =, ub =)` as
#'   resolved), and `labels` (`list(est =, ub =)`).
#' @seealso [run_subgroup_sims()], [subgroup_glm()]
#' @export
summary.subgroup_sims <- function(object,
                                  hr_true = object$hr_true,
                                  probs = c(0.01, 0.50, 0.99),
                                  est_thresholds = NULL,
                                  ub_thresholds = NULL,
                                  est_label = NULL,
                                  ub_label = NULL,
                                  single_cats = c("ITT", "Clinical",
                                                  "Continuous"),
                                  combo_cats = c("Interaction", "3-way",
                                                 "Random"),
                                  pr_thresh_highrisk = 0.10,
                                  itt_name = "All Patients",
                                  cat_cols = c(
                                    "ITT"         = "black",
                                    "Clinical"    = "steelblue",
                                    "Continuous"  = "darkgreen",
                                    "Interaction" = "purple4",
                                    "3-way"       = "violetred3",
                                    "Random"      = "darkorange"),
                                  ...) {

  sim_hrs <- object$sim_hrs
  sim_ubs <- object$sim_ubs
  sim_ns  <- object$sim_ns
  n_sims  <- object$n_sims
  sg_names <- colnames(sim_hrs)

  # Median and 1st-99th percentile ECI across simulations
  hr_q <- t(apply(sim_hrs, 2, stats::quantile, probs = probs, na.rm = TRUE))
  ub_q <- t(apply(sim_ubs, 2, stats::quantile, probs = probs, na.rm = TRUE))

  # Convergence: number of simulations where estimation succeeded
  n_valid   <- apply(sim_hrs, 2, function(x) sum(!is.na(x)))
  pct_valid <- round(100 * n_valid / n_sims, 1)

  mean_n <- round(apply(sim_ns, 2, mean, na.rm = TRUE))

  # Effect-aware resolution (Phase 4.2): explicit arguments -> the
  # result's `effect` metadata -> HR legacy. A legacy call (no metadata,
  # no overrides) takes the short-circuit branches below, whose
  # expressions are byte-for-byte the pre-4.2 code.
  eff <- object$effect
  legacy <- is.null(eff) && is.null(est_thresholds) &&
    is.null(ub_thresholds) && is.null(est_label) && is.null(ub_label)
  if (!legacy) {
    if (is.null(est_thresholds)) est_thresholds <-
      if (!is.null(eff$est_thresholds)) eff$est_thresholds else c(0.5, 1.0)
    if (is.null(ub_thresholds)) ub_thresholds <-
      if (!is.null(eff$ub_thresholds)) eff$ub_thresholds else c(2, 3)
    if (is.null(est_label)) est_label <-
      if (!is.null(eff$est_label)) eff$est_label else "HR"
    if (is.null(ub_label)) ub_label <-
      if (!is.null(eff$ub_label)) eff$ub_label else "UB"
    stopifnot(is.numeric(est_thresholds), length(est_thresholds) == 2L,
              is.numeric(ub_thresholds),  length(ub_thresholds)  == 2L,
              is.character(est_label), length(est_label) == 1L,
              is.character(ub_label),  length(ub_label)  == 1L)
  }

  if (legacy) {
    # Tail probabilities, denominated over converged fits (na.rm = TRUE)
    pr_hr_lt050 <- apply(sim_hrs, 2, function(x) mean(x < 0.5,  na.rm = TRUE))
    pr_hr_gt1   <- apply(sim_hrs, 2, function(x) mean(x > 1.0,  na.rm = TRUE))
    pr_ub_ge2   <- apply(sim_ubs, 2, function(x) mean(x >= 2.0, na.rm = TRUE))
    pr_ub_ge3   <- apply(sim_ubs, 2, function(x) mean(x >= 3.0, na.rm = TRUE))
  } else {
    # Same denominators (converged fits, na.rm = TRUE); an NA threshold
    # disables its tail, yielding an all-NA vector (rendered "-").
    tail_pr <- function(m, thr, op) {
      if (is.na(thr)) {
        return(stats::setNames(rep(NA_real_, ncol(m)), colnames(m)))
      }
      apply(m, 2, function(x) mean(op(x, thr), na.rm = TRUE))
    }
    pr_hr_lt050 <- tail_pr(sim_hrs, est_thresholds[1], `<`)
    pr_hr_gt1   <- tail_pr(sim_hrs, est_thresholds[2], `>`)
    pr_ub_ge2   <- tail_pr(sim_ubs, ub_thresholds[1], `>=`)
    pr_ub_ge3   <- tail_pr(sim_ubs, ub_thresholds[2], `>=`)
  }

  # Unconditional median displays (character, "-" when empty)
  mhr_uncond <- ifelse(is.na(hr_q[, 2]), "-",
                       as.character(round(hr_q[, 2], 2)))
  mub_uncond <- ifelse(is.na(ub_q[, 2]), "-",
                       as.character(round(ub_q[, 2], 2)))

  fmt_pct <- function(x) ifelse(is.nan(x), "-",
                                paste0(round(100 * x, 1), "%"))

  if (legacy) {
    results_tbl <- data.frame(
      Subgroup    = sg_names,
      Mean_N      = ifelse(is.nan(mean_n), "-", as.character(mean_n)),
      N_valid     = paste0(n_valid, " (", pct_valid, "%)"),
      Pr_HR_lt050 = fmt_pct(pr_hr_lt050),
      Pr_HR_gt1   = fmt_pct(pr_hr_gt1),
      mHR         = mhr_uncond,
      Pr_UB_ge2   = fmt_pct(pr_ub_ge2),
      Pr_UB_ge3   = fmt_pct(pr_ub_ge3),
      mUB         = mub_uncond,
      stringsAsFactors = FALSE
    )
    colnames(results_tbl) <- c(
      "Subgroup", "N", "#(% converged)",
      "Pr(HR<0.5)", "Pr(HR>1.0)", "mHR",
      "Pr(UB>=2)", "Pr(UB>=3)", "mUB"
    )
  } else {
    fmt_pct_g <- function(x) ifelse(is.na(x), "-",
                                    paste0(round(100 * x, 1), "%"))
    # When a resolved pair equals its legacy pair, reuse the legacy
    # literal strings so untouched headers never drift (e.g. a UB-only
    # override keeps "Pr(HR>1.0)", not "Pr(HR>1)").
    thr_strs <- function(t, legacy_vals, legacy_strs, fallbacks) {
      if (identical(as.numeric(t), legacy_vals)) return(legacy_strs)
      vapply(seq_along(t), function(i) {
        if (is.na(t[i])) fallbacks[i] else
          format(t[i], trim = TRUE, drop0trailing = TRUE)
      }, character(1L))
    }
    est_s <- thr_strs(est_thresholds, c(0.5, 1.0), c("0.5", "1.0"),
                      c("lo", "hi"))
    ub_s  <- thr_strs(ub_thresholds, c(2, 3), c("2", "3"),
                      c("t1", "t2"))
    results_tbl <- data.frame(
      Subgroup  = sg_names,
      Mean_N    = ifelse(is.nan(mean_n), "-", as.character(mean_n)),
      N_valid   = paste0(n_valid, " (", pct_valid, "%)"),
      Pr_est_lo = fmt_pct_g(pr_hr_lt050),
      Pr_est_hi = fmt_pct_g(pr_hr_gt1),
      m_est     = mhr_uncond,
      Pr_ub_t1  = fmt_pct_g(pr_ub_ge2),
      Pr_ub_t2  = fmt_pct_g(pr_ub_ge3),
      m_ub      = mub_uncond,
      stringsAsFactors = FALSE
    )
    colnames(results_tbl) <- c(
      "Subgroup", "N", "#(% converged)",
      sprintf("Pr(%s<%s)",  est_label, est_s[1]),
      sprintf("Pr(%s>%s)",  est_label, est_s[2]),
      paste0("m", est_label),
      sprintf("Pr(%s>=%s)", ub_label, ub_s[1]),
      sprintf("Pr(%s>=%s)", ub_label, ub_s[2]),
      paste0("m", ub_label)
    )
  }

  # Category mapping -- first-match order identical to the vignettes'
  # dplyr::case_when block, in base R to avoid the dependency.
  grp <- vapply(object$subgroups, `[[`, character(1L), "grp")
  cat_lab <- ifelse(grp == "ITT", "ITT",
             ifelse(grp == "Clinical", "Clinical",
             ifelse(grp == "Continuous", "Continuous",
             ifelse(grepl("^Interaction", grp), "Interaction",
             ifelse(grp == "3-way", "3-way", "Random")))))

  ok    <- !is.na(hr_q[, 2])
  ok_ub <- !is.na(ub_q[, 2])

  idx_hr_single <- which(cat_lab[ok]    %in% single_cats)
  idx_hr_combo  <- which(cat_lab[ok]    %in% combo_cats)
  idx_ub_single <- which(cat_lab[ok_ub] %in% single_cats)
  idx_ub_combo  <- which(cat_lab[ok_ub] %in% combo_cats)

  # High-risk panel: meets the first upper-bound tail threshold OR is
  # the ITT anchor; ordered ITT first, then descending tail probability.
  if (legacy) {
    include_highrisk <- ok_ub &
      (pr_ub_ge2 >= pr_thresh_highrisk | sg_names == itt_name)
    ord_highrisk <- order(
      sg_names[include_highrisk] != itt_name,
      -pr_ub_ge2[include_highrisk]
    )
  } else {
    # NA-safe: a disabled first UB tail admits only the ITT anchor.
    pr_hr2 <- ifelse(is.na(pr_ub_ge2), -Inf, pr_ub_ge2)
    include_highrisk <- ok_ub &
      (pr_hr2 >= pr_thresh_highrisk | sg_names == itt_name)
    ord_highrisk <- order(
      sg_names[include_highrisk] != itt_name,
      -pr_hr2[include_highrisk]
    )
  }

  # Per-trial N for footnotes: from the runner's sim_config when present;
  # legacy payloads (pre-wrapper, or fixed-design hand-built lists) fall
  # back to the ITT mean N, which equals the panel size by construction.
  n_per_trial <- object$sim_config$n_per_trial
  if (is.null(n_per_trial)) {
    n_per_trial <- if (itt_name %in% names(mean_n) &&
                       is.finite(mean_n[[itt_name]])) {
      as.integer(mean_n[[itt_name]])
    } else NA_integer_
  }
  n_per_trial <- as.integer(n_per_trial)

  out <- list(
    hr_q = hr_q, ub_q = ub_q,
    n_per_trial = n_per_trial,
    n_valid = n_valid, pct_valid = pct_valid, mean_n = mean_n,
    pr_hr_lt050 = pr_hr_lt050, pr_hr_gt1 = pr_hr_gt1,
    pr_ub_ge2 = pr_ub_ge2, pr_ub_ge3 = pr_ub_ge3,
    mhr_uncond = mhr_uncond, mub_uncond = mub_uncond,
    results_tbl = results_tbl,
    sg_names = sg_names, cat = cat_lab, cat_cols = cat_cols,
    ok = ok, ok_ub = ok_ub,
    idx_hr_single = idx_hr_single, idx_hr_combo = idx_hr_combo,
    idx_ub_single = idx_ub_single, idx_ub_combo = idx_ub_combo,
    n_single = length(idx_hr_single), n_combo = length(idx_hr_combo),
    highrisk = list(threshold = pr_thresh_highrisk, itt_name = itt_name,
                    include = include_highrisk, ord = ord_highrisk,
                    n = sum(include_highrisk)),
    single_cats = single_cats, combo_cats = combo_cats,
    n_sims = n_sims, design = object$design, hr_true = hr_true
  )
  if (!legacy) {
    out$effect     <- eff
    out$thresholds <- list(est = est_thresholds, ub = ub_thresholds)
    out$labels     <- list(est = est_label, ub = ub_label)
  }
  class(out) <- c("subgroup_sims_summary", "list")
  out
}

#' @export
print.subgroup_sims_summary <- function(x, ...) {
  true_lab <- if (is.null(x$hr_true)) "" else if (!is.null(x$labels)) {
    sprintf(" | true %s: %.2f", x$labels$est, x$hr_true)
  } else {
    sprintf(" | true HR: %.2f", x$hr_true)
  }
  cat(sprintf(
    "<subgroup_sims_summary>  design: %s | n_sims: %s | subgroups: %d%s\n\n",
    x$design, format(x$n_sims, big.mark = ","), length(x$sg_names),
    true_lab))
  print(x$results_tbl, row.names = FALSE)
  invisible(x)
}
