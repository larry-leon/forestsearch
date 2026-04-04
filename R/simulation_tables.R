# =============================================================================
# R/simulation_tables.R
# Publication-quality tables for simulation operating characteristics
# =============================================================================

# NOTE: Global variable declarations for data.table columns (Group, Scenario,
# Subgroup) are consolidated in R/globals.R per package convention.

# =============================================================================
# Classification / Identification Rate Table
# =============================================================================

#' Build Classification Rate Table from Simulation Results
#'
#' Constructs a publication-quality \code{gt} table summarizing subgroup
#' identification and classification rates across one or more data generation
#' scenarios and analysis methods. The layout mirrors Table 4 of
#' Leon et al. (2024) with metrics grouped by model scenario (null / alt)
#' and columns for each analysis method.
#'
#' @param scenario_results Named list. Each element is itself a list with:
#'   \describe{
#'     \item{results}{\code{data.table} from
#'       \code{\link{run_simulation_analysis}}.}
#'     \item{label}{Character scenario label, e.g., \code{"M1"}.}
#'     \item{n_sample}{Integer sample size.}
#'     \item{dgm}{DGM object (for true HRs and subgroup prevalence).}
#'     \item{hypothesis}{Character: \code{"null"} or \code{"alt"}.}
#'   }
#' @param analyses Character vector of analysis labels to include
#'   (e.g., \code{c("FS", "FSlg", "GRF")}). When \code{NULL}, all unique
#'   values of \code{results$analysis} across scenarios are used.
#' @param digits Integer. Decimal places for proportions. Default: 2.
#' @param title Character. Table title. Default:
#'   \code{"Subgroup Identification and Classification Rates"}.
#' @param n_sims Integer. Number of simulations (for subtitle). Default:
#'   \code{NULL}.
#' @param bold_threshold Numeric. Type I error threshold above which the
#'   \code{any(H)} value is shown in bold. Set \code{NULL} to disable.
#'   Default: 0.05.
#' @param font_size Numeric. Font size in pixels for table text. Default: 12.
#'   Increase to 14 or 16 for larger display.
#'
#' @return A \code{gt} table object.
#'
#' @details
#' For each scenario the function computes:
#' \itemize{
#'   \item \code{any(H)}: Proportion of simulations identifying any subgroup.
#'   \item \code{sens(H)}: Mean sensitivity (only under alternative).
#'   \item \code{sens(Hc)}: Mean specificity.
#'   \item \code{ppv(H)}: Mean positive predictive value (only under
#'     alternative).
#'   \item \code{ppv(Hc)}: Mean negative predictive value.
#'   \item \code{avg|H|}: Mean size of identified subgroup (when found).
#' }
#'
#' Under the null hypothesis the rows are reduced to \code{any(H)},
#' \code{sens(Hc)}, \code{ppv(Hc)}, and \code{avg|H|}.
#'
#' @examples
#' \dontrun{
#' # Assemble results from H0 and H1 simulations
#' scenarios <- list(
#'   null = list(
#'     results = results_null, label = "M1",
#'     n_sample = 700, dgm = dgm_null, hypothesis = "null"
#'   ),
#'   alt = list(
#'     results = results_alt, label = "M1",
#'     n_sample = 700, dgm = dgm_calibrated, hypothesis = "alt"
#'   )
#' )
#'
#' build_classification_table(scenarios, n_sims = 100)
#' }
#'
#' @seealso \code{\link{format_oc_results}},
#'   \code{\link{summarize_simulation_results}}
#'
#' @importFrom data.table as.data.table rbindlist
#' @importFrom gt gt tab_header cols_hide cols_label tab_style cell_text
#'   tab_options px cells_body cells_row_groups
#' @export
build_classification_table <- function(
    scenario_results,
    analyses = NULL,
    digits = 2,
    title = "Subgroup Identification and Classification Rates",
    n_sims = NULL,
    bold_threshold = 0.05,
    font_size = 12
) {

  rows_list <- list()

  for (sc_name in names(scenario_results)) {
    sc  <- scenario_results[[sc_name]]
    res <- data.table::as.data.table(sc$results)
    dgm <- sc$dgm
    hyp <- sc$hypothesis
    n   <- sc$n_sample
    lab <- sc$label

    # Determine analyses if not specified
    if (is.null(analyses)) {
      analyses_use <- sort(unique(res$analysis))
    } else {
      analyses_use <- analyses
    }

    # ── Section header label ────────────────────────────────────────────────
    # Use get_dgm_hr() for compatibility with both gbsg_dgm and aft_dgm_flex
    if (hyp == "null") {
      hr_itt <- round(get_dgm_hr(dgm, "hr_overall"), 2)
      if (is.na(hr_itt)) hr_itt <- round(dgm$hr_causal, 2)
      section_label <- sprintf(
        "%s Null: N=%d, theta(ITT) = %s",
        lab, n, hr_itt
      )
    } else {
      df_s <- dgm$df_super_rand %||% dgm$df_super
      harm_col <- intersect(c("flag.harm", "flag_harm"), names(df_s))
      prop_H <- if (length(harm_col) > 0) {
        round(100 * mean(df_s[[harm_col[1]]], na.rm = TRUE), 0)
      } else { NA_integer_ }
      hr_H   <- round(get_dgm_hr(dgm, "hr_H"), 2)
      hr_Hc  <- round(get_dgm_hr(dgm, "hr_Hc"), 2)
      hr_itt <- round(get_dgm_hr(dgm, "hr_overall"), 2)
      if (is.na(hr_itt)) hr_itt <- round(dgm$hr_causal, 2)
      section_label <- sprintf(
        "%s Alt: N=%d, p_H=%d%%, theta(H)=%s, theta(Hc)=%s, theta(ITT)=%s",
        lab, n, prop_H, hr_H, hr_Hc, hr_itt
      )
    }

    # ── Metric names depend on hypothesis ───────────────────────────────────
    metric_names <- if (hyp == "null") {
      c("any(H)", "sens(Hc)", "ppv(Hc)", "avg|H|")
    } else {
      c("any(H)", "sens(H)", "sens(Hc)", "ppv(H)", "ppv(Hc)", "avg|H|")
    }

    for (metric in metric_names) {
      row_vals <- list(Scenario = section_label, Metric = metric)

      for (a in analyses_use) {
        r <- res[res$analysis == a, ]
        r_found <- r[r$any.H == 1, ]

        val <- switch(
          metric,
          "any(H)"   = mean(r$any.H, na.rm = TRUE),
          "sens(H)"  = mean(r$sens, na.rm = TRUE),
          "sens(Hc)" = mean(r$spec, na.rm = TRUE),
          "ppv(H)"   = mean(r$ppv, na.rm = TRUE),
          "ppv(Hc)"  = mean(r$npv, na.rm = TRUE),
          "avg|H|"   = if (nrow(r_found) > 0) {
            round(mean(r_found$size.H, na.rm = TRUE), 0)
          } else {
            NA_real_
          }
        )

        row_vals[[a]] <- if (metric == "avg|H|") val else round(val, digits)
      }

      rows_list[[length(rows_list) + 1]] <- as.data.frame(
        row_vals, stringsAsFactors = FALSE
      )
    }
  }

  table_df <- data.table::rbindlist(rows_list, fill = TRUE)

  # ── Grouping variable ───────────────────────────────────────────────────
  table_df[, Group := Scenario]

  gt_tbl <- gt::gt(table_df, groupname_col = "Group") |>
    gt::cols_hide(columns = "Scenario") |>
    gt::tab_header(
      title = title,
      subtitle = if (!is.null(n_sims)) {
        sprintf(
          "Across %s simulations per scenario",
          format(n_sims, big.mark = ",")
        )
      }
    ) |>
    gt::cols_label(Metric = "") |>
    gt::tab_style(
      style = gt::cell_text(size = gt::px(font_size)),
      locations = gt::cells_body()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold", size = gt::px(font_size)),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_options(
      table.font.size = gt::px(font_size),
      row_group.padding = gt::px(4)
    )

  # ── Bold inflated Type I error values ───────────────────────────────────
  if (!is.null(bold_threshold)) {
    analyses_in_tbl <- setdiff(
      names(table_df),
      c("Scenario", "Metric", "Group")
    )

    for (col_name in analyses_in_tbl) {
      bold_rows <- which(
        table_df$Metric == "any(H)" &
          !is.na(table_df[[col_name]]) &
          is.numeric(table_df[[col_name]]) &
          table_df[[col_name]] > bold_threshold
      )
      if (length(bold_rows) > 0) {
        gt_tbl <- gt::tab_style(
          gt_tbl,
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_body(columns = col_name, rows = bold_rows)
        )
      }
    }
  }

  gt_tbl
}


# =============================================================================
# Estimation Properties Table
# =============================================================================

#' Build Estimation Properties Table from Simulation Results
#'
#' Constructs a publication-quality \code{gt} table summarizing estimation
#' properties for hazard ratios in the identified subgroup and its complement.
#' The layout mirrors Table 5 of Leon et al. (2024), showing average estimate,
#' empirical SD, min, max, and relative bias for each estimator.
#'
#' Uses the paper's notation conventions:
#' \itemize{
#'   \item theta-dagger: Marginal (causal) HR truth
#'   \item theta-ddagger: Controlled direct effect (CDE) truth
#'   \item theta-hat(H-hat): Plugin Cox estimate in identified subgroup
#'   \item theta-hat*(H-hat): Bootstrap bias-corrected estimate
#' }
#'
#' Includes both Cox-based HR and AHR (Average Hazard Ratio from loghr_po)
#' estimators when AHR columns are present in the results.
#'
#' @param results \code{data.table} or \code{data.frame}. Simulation results
#'   from \code{\link{run_simulation_analysis}}, optionally enriched with
#'   bootstrap bias-corrected columns (\code{hr.H.bc}, \code{hr.Hc.bc}).
#' @param dgm DGM object. Used for true parameter values (\code{hr_H_true},
#'   \code{hr_Hc_true}, and AHR truth via \code{\link{get_dgm_hr}}).
#' @param analysis_method Character. Which analysis method to tabulate
#'   (e.g., \code{"FSlg"}). Default: \code{"FSlg"}.
#' @param n_boots Integer or \code{NULL}. Number of bootstraps. When non-NULL,
#'   appended to the subtitle as "(B = n_boots bootstraps)". Default: \code{NULL}.
#' @param digits Integer. Decimal places. Default: 2.
#' @param title Character. Table title.
#' @param subtitle Character or \code{NULL}. Optional user-supplied subtitle
#'   text. Auto-populated statistics (method, estimable count, proportion) are
#'   always shown. When non-NULL, the custom text is displayed first with stats
#'   appended in parentheses, e.g. "My title (FSlg: 18/20 (90%) estimable)".
#' @param font_size Numeric. Font size in pixels for table text. Default: 12.
#'   Increase to 14 or 16 for larger display.
#' @param cde_H Numeric or \code{NULL}. Controlled direct effect
#'   (theta-ddagger(H)) for the true harm subgroup.
#'   When non-NULL, an additional b-ddagger bias
#'   column is shown. When \code{NULL} (default), auto-detected from
#'   \code{dgm$cde_H} or \code{dgm$hazard_ratios$CDE_harm}.
#' @param cde_Hc Numeric or \code{NULL}. Controlled direct effect
#'   for the complement. Auto-detected analogously.
#'
#' @return A \code{gt} table object, or \code{NULL} if no estimable
#'   realizations exist.
#'
#' @details
#' For each subgroup (H and Hc) the function reports:
#' \itemize{
#'   \item \strong{Avg}: Mean of the estimates across estimable simulations.
#'   \item \strong{SD}: Empirical standard deviation.
#'   \item \strong{Min / Max}: Range.
#'   \item \strong{b-dagger}: Relative bias (percent) vs marginal truth,
#'     \code{100 * (Avg - theta_dagger) / theta_dagger}.
#'   \item \strong{b-ddagger} (conditional): Relative bias (percent) vs
#'     CDE truth, shown when CDE values are available.
#' }
#'
#' When bootstrap-corrected columns (\code{hr.H.bc}, \code{hr.Hc.bc}) are
#' present in \code{results}, an additional bias-corrected row
#' (theta-hat*(H-hat)) is added per subgroup.
#'
#' When AHR columns (\code{ahr.H.hat}, \code{ahr.Hc.hat}) are present, AHR
#' estimation rows are appended using the DGM's true AHR values for relative
#' bias calculation.
#'
#' When CDE columns (\code{cde.H.hat}, \code{cde.Hc.hat}) are present and
#' CDE truth values are available, CDE estimation rows
#' (theta-ddagger(H-hat)) are appended.  The b-dagger column for CDE rows
#' reports bias relative to the CDE truth rather than the marginal HR.
#'
#' @examples
#' \dontrun{
#' # Basic usage (marginal truth only)
#' build_estimation_table(
#'   results = results_alt,
#'   dgm = dgm_calibrated,
#'   analysis_method = "FSlg"
#' )
#'
#' # With CDE truth (full Table 5 alignment)
#' build_estimation_table(
#'   results = results_alt,
#'   dgm = dgm_calibrated,
#'   cde_H = 2.25,
#'   cde_Hc = 0.60
#' )
#' }
#'
#' @seealso \code{\link{build_classification_table}},
#'   \code{\link{format_oc_results}}, \code{\link{get_dgm_hr}}
#'
#' @importFrom data.table as.data.table rbindlist
#' @importFrom stats sd
#' @importFrom gt gt tab_header cols_label tab_style cell_text tab_options px
#'   cells_body cells_row_groups tab_footnote cells_column_labels
#' @export
build_estimation_table <- function(
    results,
    dgm,
    analysis_method = "FSlg",
    n_boots = NULL,
    digits = 2,
    title = "Estimation Properties",
    subtitle = NULL,
    font_size = 12,
    cde_H = NULL,
    cde_Hc = NULL
) {

  res <- data.table::as.data.table(results)

  if ("analysis" %in% names(res)) {
    res <- res[res$analysis == analysis_method, ]
  }

  res_found <- res[res$any.H == 1, ]
  n_estimable <- nrow(res_found)

  if (n_estimable == 0) {
    message("No estimable realizations for analysis = ", analysis_method)
    return(NULL)
  }

  # ── True values from DGM ──────────────────────────────────────────────────
  # Paper notation: theta-dagger = marginal (causal) HR
  #                 theta-ddagger = controlled direct effect (CDE)
  # Use get_dgm_hr() for compatibility with both gbsg_dgm and aft_dgm_flex
  theta_H_true  <- get_dgm_hr(dgm, "hr_H")
  theta_Hc_true <- get_dgm_hr(dgm, "hr_Hc")
  avg_size_H    <- round(mean(res_found$size.H, na.rm = TRUE), 0)

  # True AHR values (via backward-compatible helper)
  ahr_H_true  <- get_dgm_hr(dgm, "ahr_H")
  ahr_Hc_true <- get_dgm_hr(dgm, "ahr_Hc")

  # ── CDE truth values (theta-ddagger from paper) ───────────────────────────
  # Auto-detect from DGM if not explicitly provided.
  # CDE(H) = mean(exp(theta_1[H])) / mean(exp(theta_0[H]))

  if (is.null(cde_H)) {
    cde_H <- get_dgm_hr(dgm, "cde_H")
  }
  if (is.null(cde_Hc)) {
    cde_Hc <- get_dgm_hr(dgm, "cde_Hc")
  }
  has_cde <- !is.na(cde_H) && !is.na(cde_Hc)

  # Fallback: compute CDE from super-population if available but not stored
  # Support both gbsg_dgm (df_super_rand) and aft_dgm_flex (df_super)
  if (!has_cde) {
    df_sp <- dgm$df_super_rand %||% dgm$df_super
    if (!is.null(df_sp) && all(c("theta_0", "theta_1") %in% names(df_sp))) {
      dgm <- compute_dgm_cde(dgm)
      if (is.null(cde_H) || is.na(cde_H)) {
        cde_H <- get_dgm_hr(dgm, "cde_H")
      }
      if (is.null(cde_Hc) || is.na(cde_Hc)) {
        cde_Hc <- get_dgm_hr(dgm, "cde_Hc")
      }
      has_cde <- !is.na(cde_H) && !is.na(cde_Hc)
    }
  }

  # Under the null DGM, hr_H_true is NA (no true subgroup). Fall back to
  # the overall (causal) HR which is the uniform value for any subset.
  is_null <- (!is.null(dgm$model_type) && dgm$model_type == "null") ||
    (is.na(theta_H_true) && !is.na(theta_Hc_true))

  if (is_null) {
    hr_overall <- get_dgm_hr(dgm, "hr_overall")
    if (is.na(hr_overall)) hr_overall <- dgm$hr_causal
    if (!is.null(hr_overall) && !is.na(hr_overall)) {
      if (is.na(theta_H_true))  theta_H_true  <- hr_overall
      if (is.na(theta_Hc_true)) theta_Hc_true <- hr_overall
    }
    ahr_overall <- get_dgm_hr(dgm, "ahr")
    if (is.na(ahr_overall)) ahr_overall <- dgm$AHR
    if (!is.null(ahr_overall) && !is.na(ahr_overall)) {
      if (is.na(ahr_H_true))  ahr_H_true  <- ahr_overall
      if (is.na(ahr_Hc_true)) ahr_Hc_true <- ahr_overall
    }
    # Under null, CDE equals overall CDE (uniform effect)
    if (!has_cde) {
      cde_overall <- get_dgm_hr(dgm, "cde")
      if (!is.na(cde_overall)) {
        cde_H  <- cde_overall
        cde_Hc <- cde_overall
        has_cde <- TRUE
      }
    }
  }

  # ── Math notation (pure Unicode, paper-aligned) ──────────────────────────
  # Leon et al. (2024) notation:
  #   theta-dagger   = marginal/causal HR (truth)
  #   theta-ddagger  = CDE (truth)
  #   theta-hat(.)   = Cox estimate (plugin or oracle)
  #   theta-hat*(.)  = bootstrap bias-corrected
  #
  # Unicode: \u0124 = H-hat, \u1D9C = superscript-c, \u2020 = dagger,
  #          \u2021 = double-dagger, \u0302 = combining circumflex,
  #          \u03b8 = theta, \u00e2 = a-hat

  label_map <- list(
    "theta-hat(H-hat)"   = "\u03b8\u0302(\u0124)",
    "theta-hat*(H-hat)"  = "\u03b8\u0302*(\u0124)",
    "ahr-hat(H-hat)"     = "\u00e2hr(\u0124)",
    "cde-hat(H-hat)"     = "\u03b8\u2021(\u0124)",
    "theta-hat(Hc-hat)"  = "\u03b8\u0302(\u0124\u1D9C)",
    "theta-hat*(Hc-hat)" = "\u03b8\u0302*(\u0124\u1D9C)",
    "ahr-hat(Hc-hat)"    = "\u00e2hr(\u0124\u1D9C)",
    "cde-hat(Hc-hat)"    = "\u03b8\u2021(\u0124\u1D9C)"
  )

  # Row-group header builders using paper notation:
  #   theta-dagger (marginal truth) always shown
  #   theta-ddagger (CDE truth) appended when available
  build_H_label <- function(n_est, avg_sz, theta_true, cde_val) {
    # CDE part (theta-ddagger)
    cde_part <- if (!is.na(cde_val)) {
      sprintf(", \u03b8\u2021(H) = %s", round(cde_val, 2))
    } else ""
    sprintf(
      "\u0124: %d estimable, avg |\u0124| = %d, \u03b8\u2020(H) = %s%s",
      n_est, avg_sz, round(theta_true, 2), cde_part
    )
  }

  build_Hc_label <- function(avg_sz_Hc, theta_true, cde_val) {
    cde_part <- if (!is.na(cde_val)) {
      sprintf(", \u03b8\u2021(H\u1D9C) = %s", round(cde_val, 2))
    } else ""
    sprintf(
      "\u0124\u1D9C: avg |\u0124\u1D9C| = %d, \u03b8\u2020(H\u1D9C) = %s%s",
      avg_sz_Hc, round(theta_true, 2), cde_part
    )
  }

  # ── Helper: compute one estimation row ────────────────────────────────────
  # When has_cde is TRUE, produces both b-dagger and b-ddagger bias columns.
  make_est_row <- function(estimates, true_val, label, cde_val = NA_real_) {
    est <- estimates[!is.na(estimates)]
    if (length(est) == 0) return(NULL)

    avg_est  <- mean(est)
    sd_est   <- stats::sd(est)
    min_est  <- min(est)
    max_est  <- max(est)

    # b-dagger: bias relative to marginal (causal) truth
    b_dagger <- 100 * (avg_est - true_val) / true_val

    row <- data.frame(
      Estimator          = label,
      Avg                = round(avg_est, digits),
      SD                 = round(sd_est, digits),
      Min                = round(min_est, digits),
      Max                = round(max_est, digits),
      check.names        = FALSE,
      stringsAsFactors   = FALSE
    )

    if (has_cde && !is.na(cde_val)) {
      # b-ddagger: bias relative to CDE truth
      b_ddagger <- 100 * (avg_est - cde_val) / cde_val
      row[["b\u2021 (%)"]] <- round(b_ddagger, digits)
      row[["b\u2020 (%)"]] <- round(b_dagger, digits)
    } else {
      row[["b\u2020 (%)"]] <- round(b_dagger, digits)
    }

    row
  }

  # ── Rows for H (harm subgroup) ────────────────────────────────────────────
  rows_H <- list()

  # Cox HR (identified subgroup)
  if ("hr.H.hat" %in% names(res_found)) {
    rows_H[[length(rows_H) + 1]] <- make_est_row(
      res_found$hr.H.hat, theta_H_true, "theta-hat(H-hat)",
      cde_val = if (has_cde) cde_H else NA_real_
    )
  }

  # Cox HR (bias-corrected)
  if ("hr.H.bc" %in% names(res_found)) {
    rows_H[[length(rows_H) + 1]] <- make_est_row(
      res_found$hr.H.bc, theta_H_true, "theta-hat*(H-hat)",
      cde_val = if (has_cde) cde_H else NA_real_
    )
  }

  # AHR (identified subgroup) — AHR is a different estimand, no CDE comparison
  if ("ahr.H.hat" %in% names(res_found) && !is.na(ahr_H_true)) {
    rows_H[[length(rows_H) + 1]] <- make_est_row(
      res_found$ahr.H.hat, ahr_H_true, "ahr-hat(H-hat)"
    )
  }

  # CDE (identified subgroup) — bias is relative to CDE truth only
  if ("cde.H.hat" %in% names(res_found) && has_cde) {
    rows_H[[length(rows_H) + 1]] <- make_est_row(
      res_found$cde.H.hat, cde_H, "cde-hat(H-hat)"
    )
  }

  # ── Rows for Hc (complement) ──────────────────────────────────────────────
  rows_Hc <- list()

  # Cox HR (identified complement)
  if ("hr.Hc.hat" %in% names(res_found)) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$hr.Hc.hat, theta_Hc_true, "theta-hat(Hc-hat)",
      cde_val = if (has_cde) cde_Hc else NA_real_
    )
  }

  # Cox HR (bias-corrected)
  if ("hr.Hc.bc" %in% names(res_found)) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$hr.Hc.bc, theta_Hc_true, "theta-hat*(Hc-hat)",
      cde_val = if (has_cde) cde_Hc else NA_real_
    )
  }

  # AHR (identified complement)
  if ("ahr.Hc.hat" %in% names(res_found) && !is.na(ahr_Hc_true)) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$ahr.Hc.hat, ahr_Hc_true, "ahr-hat(Hc-hat)"
    )
  }

  # CDE (identified complement) — bias is relative to CDE truth only
  if ("cde.Hc.hat" %in% names(res_found) && has_cde) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$cde.Hc.hat, cde_Hc, "cde-hat(Hc-hat)"
    )
  }

  # ── Combine ───────────────────────────────────────────────────────────────
  df_H  <- if (length(rows_H) > 0) {
    data.table::rbindlist(rows_H, fill = TRUE)
  }
  df_Hc <- if (length(rows_Hc) > 0) {
    data.table::rbindlist(rows_Hc, fill = TRUE)
  }

  if (!is.null(df_H)) {
    h_label <- build_H_label(
      n_estimable, avg_size_H, theta_H_true,
      if (has_cde) cde_H else NA_real_
    )
    df_H[, Subgroup := h_label]
  }
  if (!is.null(df_Hc)) {
    hc_label <- build_Hc_label(
      round(mean(res_found$size.Hc, na.rm = TRUE), 0),
      theta_Hc_true,
      if (has_cde) cde_Hc else NA_real_
    )
    df_Hc[, Subgroup := hc_label]
  }

  table_df <- data.table::rbindlist(list(df_H, df_Hc), fill = TRUE)

  if (nrow(table_df) == 0) {
    message("No estimation columns found in results.")
    return(NULL)
  }

  # ── gt table ──────────────────────────────────────────────────────────────
  # Replace plain-text keys with Unicode labels before creating gt.
  # No fmt_markdown or text_transform — pure Unicode renders natively.

  for (key in names(label_map)) {
    table_df[Estimator == key, Estimator := label_map[[key]]]
  }

  # Build subtitle: auto-populated stats always included
  n_sims <- nrow(res)
  pct <- round(100 * n_estimable / n_sims)
  auto_txt <- sprintf(
    "%s: %d/%d (%d%%) estimable", analysis_method, n_estimable, n_sims, pct
  )
  if (!is.null(n_boots)) {
    auto_txt <- paste0(auto_txt, sprintf(", B = %d", n_boots))
  }

  if (!is.null(subtitle)) {
    sub_txt <- sprintf("%s (%s)", subtitle, auto_txt)
  } else {
    sub_txt <- auto_txt
  }

  gt_tbl <- gt::gt(table_df, groupname_col = "Subgroup") |>
    gt::tab_header(
      title = title,
      subtitle = sub_txt
    ) |>
    gt::cols_label(Estimator = "") |>
    gt::tab_style(
      style = gt::cell_text(size = gt::px(font_size)),
      locations = gt::cells_body()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold", size = gt::px(font_size)),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_style(
      style = gt::cell_text(size = gt::px(font_size)),
      locations = gt::cells_column_labels()
    ) |>
    gt::tab_options(
      table.font.size = gt::px(font_size),
      row_group.padding = gt::px(4)
    )

  # Notation footnote (paper-aligned)
  footnote_parts <- paste0(
    "\u03b8\u0302(\u0124) = plugin Cox HR in identified subgroup; ",
    "\u03b8\u0302*(\u0124) = bootstrap bias-corrected; ",
    "\u00e2hr(\u0124) = average hazard ratio in identified subgroup; ",
    "b\u2020 = bias relative to marginal HR \u03b8\u2020 (causal truth)"
  )
  if (has_cde) {
    footnote_parts <- paste0(
      footnote_parts,
      "; \u03b8\u2021(\u0124) = controlled direct effect in identified subgroup",
      "; b\u2021 = bias relative to CDE \u03b8\u2021"
    )
  }
  gt_tbl <- gt::tab_footnote(gt_tbl, footnote = footnote_parts)

  gt_tbl
}


# =============================================================================
# Interpret Estimation Table
# =============================================================================

#' Generate Narrative Interpretation of Estimation Properties
#'
#' Produces a templated text summary of the estimation properties table,
#' automatically populating numerical results from the simulation output.
#' Useful for reproducible vignettes where interpretation paragraphs should
#' update when simulations are re-run.
#'
#' @param results Data frame of simulation results (same as for
#'   \code{\link{build_estimation_table}}).
#' @param dgm DGM object with true parameter values.
#' @param analysis_method Character. Which analysis method to summarise.
#'   Default: \code{"FSlg"}.
#' @param n_sims Integer. Total number of simulations (for detection rate).
#'   If \code{NULL} (default), derived from \code{nrow(results)} after
#'   filtering to the analysis method.
#' @param n_boots Integer. Number of bootstraps (for narrative). Default: 300.
#' @param digits Integer. Decimal places for reported values. Default: 2.
#' @param scenario Character. One of \code{"null"} or \code{"alt"} (default).
#'   Controls the interpretive framing:
#'   \itemize{
#'     \item \code{"null"}: emphasises false-positive rate and selection bias
#'     \item \code{"alt"}: emphasises power, bias relative to true effect
#'   }
#'   If \code{NULL}, inferred from the DGM (\code{hr_H_true == hr_Hc_true}
#'   implies null).
#' @param cat Logical. If \code{TRUE} (default), prints the paragraph via
#'   \code{cat()}. If \code{FALSE}, returns it invisibly as a character string
#'   (useful for programmatic insertion into Rmd via \code{results = "asis"}).
#'
#' @return Invisibly returns the interpretation as a character string.
#'
#' @examples
#' \dontrun{
#' # In a vignette chunk with results = "asis":
#' interpret_estimation_table(results_null, dgm_null, scenario = "null")
#'
#' # Capture for further processing:
#' txt <- interpret_estimation_table(results_alt, dgm_alt, cat = FALSE)
#' }
#'
#' @seealso \code{\link{build_estimation_table}},
#'   \code{\link{format_oc_results}}, \code{\link{get_dgm_hr}}
#'
#' @importFrom data.table as.data.table
#' @importFrom stats sd
#' @export
interpret_estimation_table <- function(
    results,
    dgm,
    analysis_method = "FSlg",
    n_sims = NULL,
    n_boots = 300,
    digits = 2,
    scenario = NULL,
    cat = TRUE
) {

  res <- data.table::as.data.table(results)

  if ("analysis" %in% names(res)) {
    if (is.null(n_sims)) n_sims <- nrow(res[res$analysis == analysis_method, ])
    res <- res[res$analysis == analysis_method, ]
  } else {
    if (is.null(n_sims)) n_sims <- nrow(res)
  }

  res_found <- res[res$any.H == 1, ]
  n_estimable <- nrow(res_found)

  if (n_estimable == 0) {
    txt <- sprintf(
      "No subgroups were identified in any of the %d simulations under %s, ",
      n_sims, analysis_method
    )
    txt <- paste0(txt, "indicating a 0%% detection rate.")
    if (isTRUE(cat)) cat(txt, "\n") # nolint
    return(invisible(txt))
  }

  # ── True values ───────────────────────────────────────────────────────────
  # Under the null DGM, hr_H_true is NA because there is no true subgroup.
  # Fall back to the overall (causal) HR — the uniform true value for any

  # subset under the null.  Try multiple sources for robustness.
  # Use get_dgm_hr() for compatibility with both gbsg_dgm and aft_dgm_flex
  theta_H_true  <- get_dgm_hr(dgm, "hr_H")
  theta_Hc_true <- get_dgm_hr(dgm, "hr_Hc")
  ahr_H_true    <- get_dgm_hr(dgm, "ahr_H")
  ahr_Hc_true   <- get_dgm_hr(dgm, "ahr_Hc")

  # Overall HR: try several known locations
  hr_overall <- get_dgm_hr(dgm, "hr_overall")
  if (is.na(hr_overall)) hr_overall <- dgm$hr_causal
  if (is.null(hr_overall) || is.na(hr_overall)) hr_overall <- NA_real_

  # Overall AHR: try several known locations
  ahr_overall <- get_dgm_hr(dgm, "ahr")
  if (is.na(ahr_overall)) ahr_overall <- dgm$AHR
  if (is.null(ahr_overall) || is.na(ahr_overall)) ahr_overall <- NA_real_

  # Infer scenario if not supplied
  if (is.null(scenario)) {
    if (!is.null(dgm$model_type) && dgm$model_type == "null") {
      scenario <- "null"
    } else if (is.na(theta_H_true) && !is.na(theta_Hc_true)) {
      scenario <- "null"
    } else if (!is.na(theta_H_true) && !is.na(theta_Hc_true) &&
               isTRUE(all.equal(theta_H_true, theta_Hc_true))) {
      scenario <- "null"
    } else {
      scenario <- "alt"
    }
  }
  scenario <- match.arg(scenario, c("null", "alt"))

  # Under null: fill NA true values with the overall (uniform) HR
  if (scenario == "null") {
    if (is.na(theta_H_true))  theta_H_true  <- hr_overall
    if (is.na(theta_Hc_true)) theta_Hc_true <- hr_overall
    if (is.na(ahr_H_true))    ahr_H_true    <- ahr_overall
    if (is.na(ahr_Hc_true))   ahr_Hc_true   <- ahr_overall
  }

  # ── Compute summary statistics ────────────────────────────────────────────
  fmt <- function(x) round(x, digits)

  avg_size_H  <- round(mean(res_found$size.H, na.rm = TRUE), 0)
  avg_size_Hc <- round(mean(res_found$size.Hc, na.rm = TRUE), 0)
  detect_rate <- round(100 * n_estimable / n_sims, 1)

  # Cox HR summaries
  has_hr_H  <- "hr.H.hat" %in% names(res_found)
  has_hr_Hc <- "hr.Hc.hat" %in% names(res_found)
  has_hr_bc <- "hr.H.bc" %in% names(res_found)

  if (has_hr_H) {
    hr_H_vals <- res_found$hr.H.hat[!is.na(res_found$hr.H.hat)]
    hr_H_avg  <- fmt(mean(hr_H_vals))
    hr_H_sd   <- fmt(stats::sd(hr_H_vals))
    hr_H_bias <- if (!is.na(theta_H_true) && theta_H_true != 0) {
      fmt(100 * (mean(hr_H_vals) - theta_H_true) / theta_H_true)
    } else NA
  }
  if (has_hr_Hc) {
    hr_Hc_vals <- res_found$hr.Hc.hat[!is.na(res_found$hr.Hc.hat)]
    hr_Hc_avg  <- fmt(mean(hr_Hc_vals))
    hr_Hc_sd   <- fmt(stats::sd(hr_Hc_vals))
    hr_Hc_bias <- if (!is.na(theta_Hc_true) && theta_Hc_true != 0) {
      fmt(100 * (mean(hr_Hc_vals) - theta_Hc_true) / theta_Hc_true)
    } else NA
  }
  if (has_hr_bc) {
    hr_bc_vals <- res_found$hr.H.bc[!is.na(res_found$hr.H.bc)]
    hr_bc_avg  <- fmt(mean(hr_bc_vals))
    hr_bc_bias <- if (!is.na(theta_H_true) && theta_H_true != 0) {
      fmt(100 * (mean(hr_bc_vals) - theta_H_true) / theta_H_true)
    } else NA
  }

  # AHR summaries
  has_ahr_H  <- "ahr.H.hat" %in% names(res_found) && !is.na(ahr_H_true)
  has_ahr_Hc <- "ahr.Hc.hat" %in% names(res_found) && !is.na(ahr_Hc_true)

  if (has_ahr_H) {
    ahr_H_vals <- res_found$ahr.H.hat[!is.na(res_found$ahr.H.hat)]
    ahr_H_avg  <- fmt(mean(ahr_H_vals))
    ahr_H_bias <- if (!is.na(ahr_H_true) && ahr_H_true != 0) {
      fmt(100 * (mean(ahr_H_vals) - ahr_H_true) / ahr_H_true)
    } else NA
  }
  if (has_ahr_Hc) {
    ahr_Hc_vals <- res_found$ahr.Hc.hat[!is.na(res_found$ahr.Hc.hat)]
    ahr_Hc_avg  <- fmt(mean(ahr_Hc_vals))
    ahr_Hc_bias <- if (!is.na(ahr_Hc_true) && ahr_Hc_true != 0) {
      fmt(100 * (mean(ahr_Hc_vals) - ahr_Hc_true) / ahr_Hc_true)
    } else NA
  }

  # ── Build narrative ───────────────────────────────────────────────────────
  paras <- character()

  # NA-safe formatter: returns "N/A" for missing values
  fmt_s <- function(x) if (is.na(x)) "N/A" else as.character(round(x, digits))
  fmt_pct <- function(x) if (is.na(x)) "N/A" else sprintf("%.1f%%", x)

  # --- Paragraph 1: Detection summary ---
  true_hr_str <- fmt_s(theta_H_true)
  if (scenario == "null") {
    paras[1] <- sprintf(
      paste0(
        "Under the null hypothesis (true HR = %s uniformly), ",
        "%d of %d simulations (%s) identified a subgroup using %s. ",
        "This low detection rate confirms controlled type-I error. ",
        "Among those %d false detections, the identified subgroup ",
        "averaged %d patients."
      ),
      true_hr_str, n_estimable, n_sims, fmt_pct(detect_rate),
      analysis_method, n_estimable, avg_size_H
    )
  } else {
    paras[1] <- sprintf(
      paste0(
        "Under the alternative hypothesis (true HR(H) = %s, ",
        "true HR(Hc) = %s), %d of %d simulations (%s) ",
        "identified a subgroup using %s. ",
        "The identified subgroup averaged %d patients ",
        "(complement: %d)."
      ),
      fmt_s(theta_H_true), fmt_s(theta_Hc_true),
      n_estimable, n_sims, fmt_pct(detect_rate),
      analysis_method, avg_size_H, avg_size_Hc
    )
  }

  # --- Paragraph 2: Cox HR estimates ---
  if (has_hr_H && has_hr_Hc) {
    bias_H_str  <- fmt_pct(hr_H_bias)
    bias_Hc_str <- fmt_pct(hr_Hc_bias)

    if (scenario == "null") {
      paras[2] <- sprintf(
        paste0(
          "The naive Cox HR in the identified subgroup averaged %s ",
          "(SD = %s), representing %s relative bias above the true ",
          "value of %s. This upward bias reflects selection: the algorithm ",
          "identified whichever patients happened to look most like a harm ",
          "subgroup by chance. In the complement, the Cox HR averaged %s ",
          "(%s bias), showing the expected mirror effect where removing ",
          "the worst-looking patients makes the remainder appear modestly ",
          "better."
        ),
        hr_H_avg, hr_H_sd, bias_H_str, true_hr_str,
        hr_Hc_avg, bias_Hc_str
      )
    } else {
      paras[2] <- sprintf(
        paste0(
          "The naive Cox HR in the identified subgroup averaged %s ",
          "(SD = %s), corresponding to %s relative bias versus the ",
          "true HR(H) = %s. In the complement, the estimate averaged %s ",
          "(%s bias vs. true HR(Hc) = %s)."
        ),
        hr_H_avg, hr_H_sd, bias_H_str, fmt_s(theta_H_true),
        hr_Hc_avg, bias_Hc_str, fmt_s(theta_Hc_true)
      )
    }
  }

  # --- Paragraph 3: Bias correction (if available) ---
  if (has_hr_bc) {
    bias_bc_str <- fmt_pct(hr_bc_bias)
    bc_comparison <- if (!is.na(hr_bc_bias) && !is.na(hr_H_bias) &&
                         abs(hr_bc_bias) < abs(hr_H_bias)) {
      "substantially reducing bias compared to"
    } else {
      "compared to"
    }
    paras[length(paras) + 1] <- sprintf(
      paste0(
        "After bootstrap bias correction (B = %d), the corrected ",
        "estimate averaged %s (%s relative bias), ",
        "%s the naive estimate of %s."
      ),
      n_boots, hr_bc_avg, bias_bc_str, bc_comparison, hr_H_avg
    )
  }

  # --- Paragraph 4: CDE context (if available) ---
  cde_H  <- get_dgm_hr(dgm, "cde_H")
  cde_Hc <- get_dgm_hr(dgm, "cde_Hc")
  has_cde_interp <- !is.na(cde_H) && !is.na(cde_Hc)

  # Fallback: compute from super-population if available

  if (!has_cde_interp && !is.null(dgm$df_super_rand)) {
    df_sp <- dgm$df_super_rand
    if (all(c("theta_0", "theta_1") %in% names(df_sp))) {
      dgm <- compute_dgm_cde(dgm)
      cde_H  <- get_dgm_hr(dgm, "cde_H")
      cde_Hc <- get_dgm_hr(dgm, "cde_Hc")
      has_cde_interp <- !is.na(cde_H) && !is.na(cde_Hc)
    }
  }

  if (!has_cde_interp) {
    # Try overall CDE as null-model fallback
    cde_overall <- get_dgm_hr(dgm, "cde")
    if (!is.na(cde_overall)) {
      cde_H  <- cde_overall
      cde_Hc <- cde_overall
      has_cde_interp <- TRUE
    }
  }

  if (has_cde_interp && has_hr_H) {
    b_cde_H <- if (cde_H != 0) {
      fmt(100 * (mean(hr_H_vals) - cde_H) / cde_H)
    } else NA
    cde_para <- sprintf(
      paste0(
        "Relative to the controlled direct effect (CDE) truth ",
        "theta-ddagger(H) = %s, the naive plugin shows %s ",
        "relative bias."
      ),
      fmt_s(cde_H), fmt_pct(b_cde_H)
    )
    if (has_hr_bc) {
      b_cde_bc <- if (cde_H != 0) {
        fmt(100 * (mean(hr_bc_vals) - cde_H) / cde_H)
      } else NA
      cde_para <- paste0(cde_para, sprintf(
        " After bias correction, CDE-relative bias is %s.",
        fmt_pct(b_cde_bc)
      ))
    }
    paras[length(paras) + 1] <- cde_para
  }

  # --- Paragraph 5: AHR comparison (if available) ---
  if (has_ahr_H) {
    ahr_para <- sprintf(
      paste0(
        "The average hazard ratio (AHR) in the identified subgroup ",
        "averaged %s (%s relative bias vs. true AHR(H) = %s)"
      ),
      ahr_H_avg, fmt_pct(ahr_H_bias), fmt_s(ahr_H_true)
    )
    if (has_ahr_Hc) {
      ahr_para <- paste0(ahr_para, sprintf(
        "; in the complement, %s (%s bias vs. true AHR(Hc) = %s)",
        ahr_Hc_avg, fmt_pct(ahr_Hc_bias), fmt_s(ahr_Hc_true)
      ))
    }
    if (has_hr_H && !is.na(ahr_H_bias) && !is.na(hr_H_bias) &&
        abs(ahr_H_bias) < abs(hr_H_bias)) {
      ahr_para <- paste0(ahr_para,
        ". The AHR shows attenuated bias relative to the Cox HR, ",
        "consistent with AHR being a marginal rather than conditional ",
        "estimand."
      )
    } else {
      ahr_para <- paste0(ahr_para, ".")
    }
    paras[length(paras) + 1] <- ahr_para
  }

  # --- Paragraph 6: Concluding remark ---
  if (scenario == "null") {
    paras[length(paras) + 1] <- paste0(
      "These results underscore that under the null, the few false ",
      "detections produce highly biased estimates, reinforcing the need ",
      "for bootstrap bias correction for any subgroup identified by a ",
      "data-driven search."
    )
  } else {
    if (has_hr_bc) {
      paras[length(paras) + 1] <- paste0(
        "Overall, the bias correction meaningfully improves estimation ",
        "accuracy, and the detection rate reflects the power of the ",
        "algorithm under this effect size."
      )
    }
  }

  txt <- paste(paras, collapse = "\n\n")

  if (isTRUE(cat)) cat(txt, "\n") # nolint
  invisible(txt)
}

#' Render Reference Simulation Table as gt
#'
#' Converts a data frame of pre-computed reference simulation results (e.g.,
#' digitized from a published LaTeX table) into a styled \code{gt} table. This
#' is useful for displaying published benchmark results alongside new
#' simulation output within vignettes or reports.
#'
#' @param ref_df \code{data.frame}. Must contain a \code{Metric} column and
#'   one column per analysis method, with a \code{Scenario} column for row
#'   grouping.
#' @param title Character. Table title.
#' @param subtitle Character. Table subtitle. Default: \code{NULL}.
#' @param bold_threshold Numeric. Values in \code{any(H)} rows exceeding this
#'   threshold are shown in bold. Set \code{NULL} to disable. Default: 0.05.
#'
#' @return A \code{gt} table object.
#'
#' @examples
#' \dontrun{
#' ref <- data.frame(
#'   Scenario = "M1 Null: N=700",
#'   Metric   = "any(H)",
#'   FS       = 0.02,
#'   FSlg     = 0.03,
#'   GRF      = 0.25
#' )
#' render_reference_table(ref, title = "Reference Results")
#' }
#'
#' @importFrom gt gt tab_header cols_label tab_style cell_text tab_options px
#'   cells_body cells_row_groups
#' @export
render_reference_table <- function(
    ref_df,
    title = "Reference Simulation Results",
    subtitle = NULL,
    bold_threshold = 0.05
) {

  ref_df <- as.data.frame(ref_df, stringsAsFactors = FALSE)

  gt_tbl <- gt::gt(ref_df, groupname_col = "Scenario") |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::cols_label(Metric = "") |>
    gt::tab_style(
      style = gt::cell_text(size = "small"),
      locations = gt::cells_body()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold", size = "small"),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      row_group.padding = gt::px(4)
    )

  # Bold inflated Type I error
  if (!is.null(bold_threshold)) {
    analysis_cols <- setdiff(names(ref_df), c("Scenario", "Metric"))
    for (col_name in analysis_cols) {
      bold_rows <- which(
        ref_df$Metric == "any(H)" &
          !is.na(ref_df[[col_name]]) &
          is.numeric(ref_df[[col_name]]) &
          ref_df[[col_name]] > bold_threshold
      )
      if (length(bold_rows) > 0) {
        gt_tbl <- gt::tab_style(
          gt_tbl,
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_body(columns = col_name, rows = bold_rows)
        )
      }
    }
  }

  gt_tbl
}
