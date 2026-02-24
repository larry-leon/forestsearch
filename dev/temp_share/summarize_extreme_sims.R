# =============================================================================
# summarize_extreme_sims.R
#
# Standalone utility functions for summarising, saving, loading, and plotting
# extreme subgroup simulation results.
#
# WORKFLOW FOR NEW RESULTS:
# ────────────────────────
#   source("summarize_extreme_sims.R")
#
#   # (A) Load pre-computed results
#   res <- load_sim_results("results/extreme_sims_5000.rds")
#
#   # (B) Summarize
#   S <- summarize_extreme_sims(res$sim_hrs, res$sim_ubs,
#                                res$sim_ns,  res$subgroups)
#
#   # (C) Print table
#   print(S$results_tbl, row.names = FALSE)
#
#   # (D) Forest plots — all four panels
#   plot_hr_forest(S, panel = "single")
#   plot_hr_forest(S, panel = "combo")
#   plot_ub_forest(S, panel = "single")
#   plot_ub_forest(S, panel = "combo")
#
#   # (E) Compare two runs
#   compare_extreme_sims(S_old, S_new, metric = "ub")
#
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 1. SAVE / LOAD
# ─────────────────────────────────────────────────────────────────────────────

#' Save simulation matrices to RDS
#'
#' Saves the raw simulation output (three matrices plus subgroup definitions
#' and metadata) to a single RDS file that can be reloaded for summarisation
#' and plotting without re-running the simulation loop.
#'
#' @param sim_hrs   Matrix [n_sims x n_subgroups] of HR point estimates.
#' @param sim_ubs   Matrix [n_sims x n_subgroups] of upper 95% CI bounds.
#' @param sim_ns    Matrix [n_sims x n_subgroups] of subgroup sample sizes.
#' @param subgroups List of subgroup definitions (each with $name and $grp).
#' @param file      Character path for the output RDS file.
#' @param hr_true   Scalar true HR used in the DGM (default 0.70).
#' @param metadata  Optional named list of additional metadata
#'   (e.g. \code{list(dgm_seed = 99, cens_adjust = 0.42)}).
#'
#' @return Invisibly returns the file path.
#' @export
save_sim_results <- function(sim_hrs, sim_ubs, sim_ns, subgroups,
                             file, hr_true = 0.70, metadata = NULL) {
  out <- list(
    sim_hrs    = sim_hrs,
    sim_ubs    = sim_ubs,
    sim_ns     = sim_ns,
    subgroups  = subgroups,
    hr_true    = hr_true,
    n_sims     = nrow(sim_hrs),
    n_sg       = ncol(sim_hrs),
    saved_at   = Sys.time(),
    metadata   = metadata
  )
  if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE)
  saveRDS(out, file)
  message("Saved simulation results to: ", file,
          "  (", nrow(sim_hrs), " sims x ", ncol(sim_hrs), " subgroups)")
  invisible(file)
}


#' Load simulation results from RDS
#'
#' @param file Character path to the RDS file saved by \code{save_sim_results()}.
#' @return A list with components: sim_hrs, sim_ubs, sim_ns, subgroups,
#'   hr_true, n_sims, n_sg, saved_at, metadata.
#' @export
load_sim_results <- function(file) {
  stopifnot(file.exists(file))
  res <- readRDS(file)
  message("Loaded: ", file, "  (", res$n_sims, " sims x ", res$n_sg,
          " subgroups, saved ", format(res$saved_at, "%Y-%m-%d %H:%M"), ")")
  res
}


# ─────────────────────────────────────────────────────────────────────────────
# 2. SUMMARIZE
# ─────────────────────────────────────────────────────────────────────────────

#' Summarise Extreme Subgroup Simulation Results
#'
#' Computes all summary statistics needed for tables and forest plots from
#' the raw simulation matrices.  This is the central function: call it once
#' on any set of results (fresh or loaded from disk) and all downstream
#' tables and plots are driven by the returned object.
#'
#' @param sim_hrs   Matrix [n_sims x n_subgroups] of HR point estimates.
#' @param sim_ubs   Matrix [n_sims x n_subgroups] of upper 95% CI bounds.
#' @param sim_ns    Matrix [n_sims x n_subgroups] of subgroup sample sizes.
#' @param subgroups List of subgroup definitions (each with \code{$name}
#'   and \code{$grp}).
#' @param hr_true   Scalar true HR (default 0.70).
#' @param eci_probs Quantile probabilities for the empirical CI
#'   (default \code{c(0.01, 0.50, 0.99)}).
#' @param hr_threshold HR threshold for "meaningful benefit" (default 0.80).
#'
#' @return A list (class \code{"extreme_sim_summary"}) containing:
#' \describe{
#'   \item{sg_names}{Character vector of subgroup names.}
#'   \item{cat_ok}{Category labels per subgroup (ITT, Clinical, etc.).}
#'   \item{hr_q}{Matrix [n_sg x 3] of HR quantiles (lo, median, hi).}
#'   \item{ub_q}{Matrix [n_sg x 3] of UB quantiles (lo, median, hi).}
#'   \item{mean_n}{Mean subgroup N across simulations.}
#'   \item{n_valid}{Count of simulations with valid (non-NA) estimates.}
#'   \item{pct_valid}{Percentage of valid simulations.}
#'   \item{pr_hr_lt_thresh}{Pr(HR < hr_threshold).}
#'   \item{pr_ub_lt1}{Pr(UB < 1.0) --- power to exclude no effect.}
#'   \item{pr_ub_ge1}{Pr(UB >= 1.0) --- failure rate.}
#'   \item{cond_ub}{Matrix [n_sg x 3] of conditional UB quantiles given UB >= 1
#'     (median, 75th, 99th percentiles).}
#'   \item{results_tbl}{Formatted data.frame ready for printing.}
#'   \item{ok, ok_ub}{Logical vectors: subgroups with valid HR/UB medians.}
#'   \item{idx_hr_single, idx_hr_combo}{Index vectors for HR forest panels.}
#'   \item{idx_ub_single, idx_ub_combo}{Index vectors for UB forest panels.}
#'   \item{cat_cols}{Named colour map for subgroup categories.}
#'   \item{n_sims}{Total number of simulations.}
#'   \item{hr_true}{The true HR used.}
#'   \item{hr_threshold}{The HR threshold used.}
#'   \item{lab_cmed_ub}{Formatted conditional median UB labels for annotations.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- load_sim_results("results/extreme_sims_5000.rds")
#' S   <- summarize_extreme_sims(res$sim_hrs, res$sim_ubs,
#'                                res$sim_ns, res$subgroups)
#' print(S$results_tbl, row.names = FALSE)
#' }
#' @export
summarize_extreme_sims <- function(
    sim_hrs,
    sim_ubs,
    sim_ns,
    subgroups,
    hr_true      = 0.70,
    eci_probs    = c(0.01, 0.50, 0.99),
    hr_threshold = 0.80
) {

  sg_names <- sapply(subgroups, `[[`, "name")
  n_sims   <- nrow(sim_hrs)
  n_sg     <- ncol(sim_hrs)

  # --- Quantiles ---
  hr_q <- t(apply(sim_hrs, 2, quantile, probs = eci_probs, na.rm = TRUE))
  ub_q <- t(apply(sim_ubs, 2, quantile, probs = eci_probs, na.rm = TRUE))

  # --- Counts and proportions ---
  n_valid     <- apply(sim_hrs, 2, function(x) sum(!is.na(x)))
  pct_valid   <- round(100 * n_valid / n_sims, 1)
  mean_n      <- round(apply(sim_ns, 2, mean, na.rm = TRUE))

  pr_hr_lt_thresh <- apply(sim_hrs, 2, function(x) {
    mean(x < hr_threshold, na.rm = TRUE)
  })
  pr_ub_lt1 <- apply(sim_ubs, 2, function(x) mean(x < 1.0, na.rm = TRUE))
  pr_ub_ge1 <- 1 - pr_ub_lt1

  # --- Conditional UB distribution given UB >= 1 ---
  cond_ub <- t(apply(sim_ubs, 2, function(x) {
    xx <- x[!is.na(x) & x >= 1.0]
    if (length(xx) < 2) return(c(NA_real_, NA_real_, NA_real_))
    quantile(xx, probs = c(0.50, 0.75, 0.99))
  }))

  # --- Formatted table ---
  fmt_pct <- function(x) ifelse(is.nan(x), "-", paste0(round(100 * x, 1), "%"))
  fmt_num <- function(x) ifelse(is.na(x),  "-", as.character(round(x, 2)))

  results_tbl <- data.frame(
    Subgroup      = sg_names,
    Mean_N        = ifelse(is.nan(mean_n), "-", as.character(mean_n)),
    N_valid       = paste0(n_valid, " (", pct_valid, "%)"),
    HR_ECI        = ifelse(is.na(hr_q[,2]), "-",
                      paste0(round(hr_q[,2], 2),
                             " (", round(hr_q[,1], 2), ", ",
                             round(hr_q[,3], 2), ")")),
    Pr_HR_lt      = fmt_pct(pr_hr_lt_thresh),
    Pr_UB_lt1     = fmt_pct(pr_ub_lt1),
    Pr_UB_ge1     = fmt_pct(pr_ub_ge1),
    Med_UB_cond   = fmt_num(cond_ub[, 1]),
    P99_UB_cond   = fmt_num(cond_ub[, 3]),
    stringsAsFactors = FALSE
  )
  colnames(results_tbl) <- c(
    "Subgroup", "N", "#(% converged)",
    "Median HR (1st, 99th ECI)",
    paste0("Pr(HR<", hr_threshold, ")"),
    "Pr(UB<1.0)", "Pr(UB>=1.0)",
    "Med UB|UB>=1", "P99 UB|UB>=1"
  )

  # --- Category classification ---
  sg_grps <- sapply(subgroups, `[[`, "grp")
  cat_ok <- dplyr::case_when(
    sg_grps == "ITT"               ~ "ITT",
    sg_grps == "Clinical"          ~ "Clinical",
    sg_grps == "Continuous"        ~ "Continuous",
    grepl("^Interaction", sg_grps) ~ "Interaction",
    sg_grps == "3-way"             ~ "3-way",
    TRUE                           ~ "Random"
  )

  cat_cols <- c(
    "ITT"         = "black",
    "Clinical"    = "steelblue",
    "Continuous"  = "darkgreen",
    "Interaction" = "purple4",
    "3-way"       = "violetred3",
    "Random"      = "darkorange"
  )

  # --- Index vectors for forest panels ---
  ok    <- !is.na(hr_q[, 2])
  ok_ub <- !is.na(ub_q[, 2])

  single_cats <- c("ITT", "Clinical", "Continuous")
  combo_cats  <- c("Interaction", "3-way", "Random")

  idx_hr_single <- which(cat_ok[ok]    %in% single_cats)
  idx_hr_combo  <- which(cat_ok[ok]    %in% combo_cats)
  idx_ub_single <- which(cat_ok[ok_ub] %in% single_cats)
  idx_ub_combo  <- which(cat_ok[ok_ub] %in% combo_cats)

  # --- Pre-format conditional-median labels for UB annotations ---
  lab_cmed_ub <- ifelse(
    is.na(cond_ub[ok_ub, 1]), "-",
    as.character(round(cond_ub[ok_ub, 1], 2))
  )

  # --- Assemble output ---
  out <- list(
    sg_names        = sg_names,
    cat_ok          = cat_ok,
    hr_q            = hr_q,
    ub_q            = ub_q,
    mean_n          = mean_n,
    n_valid         = n_valid,
    pct_valid       = pct_valid,
    pr_hr_lt_thresh = pr_hr_lt_thresh,
    pr_ub_lt1       = pr_ub_lt1,
    pr_ub_ge1       = pr_ub_ge1,
    cond_ub         = cond_ub,
    results_tbl     = results_tbl,
    ok              = ok,
    ok_ub           = ok_ub,
    idx_hr_single   = idx_hr_single,
    idx_hr_combo    = idx_hr_combo,
    idx_ub_single   = idx_ub_single,
    idx_ub_combo    = idx_ub_combo,
    cat_cols        = cat_cols,
    n_sims          = n_sims,
    hr_true         = hr_true,
    hr_threshold    = hr_threshold,
    lab_cmed_ub     = lab_cmed_ub
  )
  class(out) <- c("extreme_sim_summary", "list")
  out
}


#' Print method for extreme_sim_summary
#' @export
print.extreme_sim_summary <- function(x, ...) {
  cat("Extreme-subgroup simulation summary\n")
  cat("  Simulations :", x$n_sims, "\n")
  cat("  Subgroups   :", length(x$sg_names), "\n")
  cat("  True HR     :", x$hr_true, "\n\n")
  print(x$results_tbl, row.names = FALSE, ...)
  invisible(x)
}


# ─────────────────────────────────────────────────────────────────────────────
# 3. FOREST-PLOT WRAPPERS
# ─────────────────────────────────────────────────────────────────────────────

#' Plot HR forest panel from a summary object
#'
#' Wrapper around \code{gg_forest()} that extracts the correct index vectors,
#' labels, and annotation columns from a \code{summarize_extreme_sims()} result.
#'
#' @param S      An \code{extreme_sim_summary} object from
#'   \code{summarize_extreme_sims()}.
#' @param panel  \code{"single"} (ITT + Clinical + Continuous) or
#'   \code{"combo"} (Interaction + 3-way + Random benchmarks).
#' @param base_size   ggplot base font size (default 9 for beamer).
#' @param point_size  Size of median point symbol (default 1.8).
#' @param line_size   CI whisker width (default 0.55).
#' @param xlim        x-axis limits (default c(0.15, 3.5)).
#' @param widths      Patchwork column widths (default c(3.0, 5, 0.8, 1.0)).
#' @param footnote_extra Optional extra text appended to the footnote.
#'
#' @return A patchwork plot object.
#' @export
plot_hr_forest <- function(
    S,
    panel          = c("single", "combo"),
    base_size      = 9,
    point_size     = 1.8,
    line_size      = 0.55,
    xlim           = c(0.15, 3.5),
    widths         = c(3.0, 5, 0.8, 1.0),
    footnote_extra = NULL
) {

  panel <- match.arg(panel)
  ok <- S$ok

  if (panel == "single") {
    idx <- S$idx_hr_single
  } else {
    idx <- S$idx_hr_combo
  }

  sg   <- S$sg_names[ok][idx]
  cats <- S$cat_ok[ok][idx]

  fn <- paste0("Median HR + 1st-99th ECI  |  ",
               S$n_sims, " sims  |  True HR=", S$hr_true, " (red)")
  if (!is.null(footnote_extra)) fn <- paste0(fn, "  |  ", footnote_extra)

  gg_forest(
    subgroups   = sg,
    est         = S$hr_q[ok, 2][idx],
    lo          = S$hr_q[ok, 1][idx],
    hi          = S$hr_q[ok, 3][idx],
    cat_vec     = cats,
    cat_colours = S$cat_cols,
    annot       = list(
      "N"          = as.character(S$mean_n[ok][idx]),
      "Pr(UB<1.0)" = paste0(round(100 * S$pr_ub_lt1[ok][idx], 1), "%")
    ),
    ref_line    = S$hr_true,
    vert_lines  = 1.00,
    ref_col     = "firebrick",
    vert_col    = "grey55",
    vert_lty    = "dotted",
    xlim        = xlim,
    ticks_at    = c(0.20, 0.35, 0.50, 0.70, 1.00, 1.50, 2.50),
    xlog        = TRUE,
    xlab        = "Hazard Ratio (Cox, stratified)",
    footnote    = fn,
    point_size  = point_size,
    line_size   = line_size,
    base_size   = base_size,
    widths      = widths
  )
}


#' Plot UB(HR) forest panel from a summary object
#'
#' @inheritParams plot_hr_forest
#' @param xlim  x-axis limits (default c(0.30, 9.0)).
#' @param widths Patchwork column widths (default c(3.0, 5, 0.8, 0.8, 1.2)).
#' @export
plot_ub_forest <- function(
    S,
    panel          = c("single", "combo"),
    base_size      = 9,
    point_size     = 1.8,
    line_size      = 0.55,
    xlim           = c(0.30, 9.0),
    widths         = c(3.0, 5, 0.8, 0.8, 1.2),
    footnote_extra = NULL
) {

  panel <- match.arg(panel)
  ok_ub <- S$ok_ub

  if (panel == "single") {
    idx <- S$idx_ub_single
  } else {
    idx <- S$idx_ub_combo
  }

  sg   <- S$sg_names[ok_ub][idx]
  cats <- S$cat_ok[ok_ub][idx]

  fn <- paste0("Median UB(HR) + 1st-99th ECI  |  ",
               S$n_sims, " sims  |  ",
               "Dotted: UB=1.0;  Dashed: HR=", S$hr_true)
  if (!is.null(footnote_extra)) fn <- paste0(fn, "  |  ", footnote_extra)

  gg_forest(
    subgroups   = sg,
    est         = S$ub_q[ok_ub, 2][idx],
    lo          = S$ub_q[ok_ub, 1][idx],
    hi          = S$ub_q[ok_ub, 3][idx],
    cat_vec     = cats,
    cat_colours = S$cat_cols,
    annot       = list(
      "N"           = as.character(S$mean_n[ok_ub][idx]),
      "UB>=1"       = paste0(round(100 * S$pr_ub_ge1[ok_ub][idx], 1), "%"),
      "mUB | UB>=1" = S$lab_cmed_ub[idx]
    ),
    ref_line    = 1.00,
    vert_lines  = S$hr_true,
    ref_col     = "grey55",
    ref_lty     = "dotted",
    vert_col    = "firebrick",
    vert_lty    = "dashed",
    xlim        = xlim,
    ticks_at    = c(0.40, 0.70, 1.00, 1.50, 2.50, 4.00, 8.00),
    xlog        = TRUE,
    xlab        = "Upper Bound of 95% CI  [UB(HR)]",
    footnote    = fn,
    point_size  = point_size,
    line_size   = line_size,
    base_size   = base_size,
    widths      = widths
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# 4. COMPARISON UTILITY
# ─────────────────────────────────────────────────────────────────────────────

#' Compare two extreme-subgroup simulation summaries
#'
#' Produces a side-by-side data.frame comparing key metrics between two runs
#' (e.g. different sample sizes, different n_sims, or different DGMs).
#' Matching is done by subgroup name.
#'
#' @param S1    First \code{extreme_sim_summary} object.
#' @param S2    Second \code{extreme_sim_summary} object.
#' @param metric  \code{"hr"} for HR point estimates, \code{"ub"} for upper
#'   bounds, or \code{"both"}.
#' @param label1  Label for the first run (default \code{"Run 1"}).
#' @param label2  Label for the second run (default \code{"Run 2"}).
#'
#' @return A data.frame with columns:
#'   Subgroup, N_1, N_2, Median_HR_1, Median_HR_2, Delta_HR,
#'   Pr_UB_ge1_1, Pr_UB_ge1_2, Delta_Pr, Med_UB_cond_1, Med_UB_cond_2.
#' @export
compare_extreme_sims <- function(
    S1, S2,
    metric = c("both", "hr", "ub"),
    label1 = "Run 1",
    label2 = "Run 2"
) {
  metric <- match.arg(metric)

  # Match subgroups by name
  common <- intersect(S1$sg_names, S2$sg_names)
  if (length(common) == 0) stop("No common subgroup names between the two runs.")

  i1 <- match(common, S1$sg_names)
  i2 <- match(common, S2$sg_names)

  out <- data.frame(
    Subgroup = common,
    stringsAsFactors = FALSE
  )

  out[[paste0("N_", label1)]]  <- S1$mean_n[i1]
  out[[paste0("N_", label2)]]  <- S2$mean_n[i2]

  if (metric %in% c("both", "hr")) {
    out[[paste0("MedHR_", label1)]] <- round(S1$hr_q[i1, 2], 3)
    out[[paste0("MedHR_", label2)]] <- round(S2$hr_q[i2, 2], 3)
    out$Delta_MedHR <- round(S2$hr_q[i2, 2] - S1$hr_q[i1, 2], 4)
  }

  if (metric %in% c("both", "ub")) {
    out[[paste0("PrUBge1_", label1)]] <- round(S1$pr_ub_ge1[i1], 3)
    out[[paste0("PrUBge1_", label2)]] <- round(S2$pr_ub_ge1[i2], 3)
    out$Delta_PrUBge1 <- round(S2$pr_ub_ge1[i2] - S1$pr_ub_ge1[i1], 4)

    out[[paste0("MedUBcond_", label1)]] <- round(S1$cond_ub[i1, 1], 2)
    out[[paste0("MedUBcond_", label2)]] <- round(S2$cond_ub[i2, 1], 2)
  }

  out
}


# ─────────────────────────────────────────────────────────────────────────────
# 5. FIGURE-HEIGHT HELPER
# ─────────────────────────────────────────────────────────────────────────────

#' Compute recommended figure height for a forest panel
#'
#' @param S     An \code{extreme_sim_summary} object.
#' @param panel \code{"single"} or \code{"combo"}.
#' @param type  \code{"hr"} or \code{"ub"} (determines which index vector).
#' @param row_height Inches per row (default 0.32 for beamer, 0.45 for vignettes).
#' @param overhead   Fixed overhead inches (default 1.2).
#' @return Numeric figure height in inches.
#' @export
forest_fig_height <- function(
    S,
    panel      = c("single", "combo"),
    type       = c("hr", "ub"),
    row_height = 0.32,
    overhead   = 1.2
) {
  panel <- match.arg(panel)
  type  <- match.arg(type)

  idx_name <- paste0("idx_", type, "_", panel)
  n_rows <- length(S[[idx_name]])
  round(n_rows * row_height + overhead, 1)
}
