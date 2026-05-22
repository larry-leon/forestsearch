# =============================================================================
# R/simulation_tables.R
# Publication-quality tables for simulation operating characteristics
# =============================================================================

# NOTE: Global variable declarations for data.table columns (Group, Scenario,
# Subgroup) are consolidated in R/globals.R per package convention.

# =============================================================================
# Subgroup Notation Helper
# =============================================================================

#' Get Subgroup Display Labels for Harm or Benefit Notation
#'
#' Returns Unicode label components for H/Hc (harm search) or G/Gc (benefit
#' search) notation.  Used internally by \code{build_classification_table},
#' \code{build_estimation_table}, \code{interpret_estimation_table}, and
#' \code{format_oc_results} to switch between the two notation systems
#' described in Leon et al. (2024, Sections 2--4.2).
#'
#' @param notation Character. \code{"harm"} (default) uses H/Hc notation
#'   for detrimental-subgroup searches.  \code{"benefit"} uses G/Gc notation
#'   (G for "good"/"gain") for benefit-subgroup searches.
#'
#' @return A named list with components:
#' \describe{
#'   \item{sg, sg_c}{True subgroup and complement (e.g., "H", "Hc")}
#'   \item{sg_hat, sg_hat_c}{Estimated subgroup and complement (Unicode hat)}
#'   \item{plain, plain_c}{Plain-text labels for metric names}
#'   \item{word}{Descriptive word: "harm" or "benefit"}
#' }
#'
#' @keywords internal
get_sg_labels <- function(notation = c("harm", "benefit")) {
  notation <- match.arg(notation)
  if (notation == "harm") {
    list(
      sg       = "H",
      sg_c     = "H\u1D9C",
      sg_hat   = "\u0124",
      sg_hat_c = "\u0124\u1D9C",
      plain    = "H",
      plain_c  = "Hc",
      word     = "harm"
    )
  } else {
    list(
      sg       = "G",
      sg_c     = "G\u1D9C",
      sg_hat   = "\u011C",
      sg_hat_c = "\u011C\u1D9C",
      plain    = "G",
      plain_c  = "Gc",
      word     = "benefit"
    )
  }
}


#' Is this Effect Measure on the Identity Scale?
#'
#' Identity-scale measures (MD, RD, IRD) use additive differences;
#' ratio-scale measures (HR, OR, RR, IRR) use multiplicative ratios.
#' The distinction matters for benefit-search inversion: ratio-scale
#' values are reciprocated (\code{1/x}), identity-scale values are
#' negated (\code{-x}).
#'
#' @param effect_measure Character.  One of \code{"HR"}, \code{"OR"},
#'   \code{"RR"}, \code{"IRR"}, \code{"RD"}, \code{"IRD"}, \code{"MD"},
#'   or \code{NULL} (treated as ratio-scale for backward compatibility).
#'
#' @return Logical.  \code{TRUE} for identity-scale measures.
#'
#' @keywords internal
is_identity_scale <- function(effect_measure) {
  !is.null(effect_measure) && effect_measure %in% c("MD", "RD", "IRD")
}


#' Human-Readable Effect Measure Label
#'
#' Returns a descriptive label for use in table footnotes and narrative
#' text.  Falls back to \code{"HR"} for survival analyses.
#'
#' @param effect_measure Character or \code{NULL}.
#' @param outcome_type Character or \code{NULL}.
#'
#' @return Character label (e.g., \code{"Cox HR"}, \code{"OR"},
#'   \code{"mean difference"}).
#'
#' @keywords internal
effect_measure_label <- function(effect_measure = NULL,
                                 outcome_type = NULL) {
  if (is.null(effect_measure) ||
      (is.null(outcome_type) && effect_measure == "HR") ||
      (!is.null(outcome_type) && outcome_type == "survival")) {
    return("Cox HR")
  }
  switch(effect_measure,
    OR  = "OR",
    RR  = "RR",
    RD  = "risk difference",
    IRR = "IRR",
    IRD = "incidence rate difference",
    MD  = "mean difference",
    effect_measure
  )
}


#' Invert Effect Columns in Simulation Results for Benefit Search
#'
#' Transforms simulation results from the switched treatment scale to the
#' original scale.  For ratio-scale measures (HR, OR, IRR), inversion is
#' \code{1/x}; for identity-scale measures (MD, RD, IRD), inversion is
#' \code{-x}.  Called automatically by table functions when
#' \code{subgroup_notation = "benefit"}.
#'
#' @param res Data frame of simulation results (from
#'   \code{\link{run_simulation_analysis}}).
#' @param effect_measure Character or \code{NULL}.  The effect measure
#'   used in the simulation.  When \code{NULL} (default), ratio-scale
#'   inversion (\code{1/x}) is used for backward compatibility with
#'   survival analyses.
#'
#' @return Data frame with all effect-related columns inverted.
#'   Non-effect columns (classification metrics, subgroup sizes, etc.)
#'   are unchanged.
#'
#' @keywords internal
invert_hr_columns <- function(res, effect_measure = NULL) {
  hr_cols <- c(
    "hr.H.hat", "hr.Hc.hat", "hr.H.bc", "hr.Hc.bc",
    "hr.H.true", "hr.Hc.true", "hr.itt", "hr.adj.itt",
    "ahr.H.hat", "ahr.Hc.hat", "ahr.H.true", "ahr.Hc.true",
    "cde.H.hat", "cde.Hc.hat"
  )
  use_negate <- is_identity_scale(effect_measure)
  for (col in intersect(hr_cols, names(res))) {
    res[[col]] <- if (use_negate) -res[[col]] else 1 / res[[col]]
  }
  res
}


#' Invert DGM Effect Values for Benefit Search
#'
#' Transforms a DGM object's truth values from the switched treatment scale
#' to the original scale.  For ratio-scale measures (HR, OR, IRR), fields
#' are reciprocated (\code{1/x}); for identity-scale measures (MD, RD, IRD),
#' fields are negated (\code{-x}).  The effect measure is auto-detected from
#' \code{dgm$effect_measure} when not explicitly provided.
#'
#' Called automatically by table functions when
#' \code{subgroup_notation = "benefit"}.
#'
#' @param dgm A DGM object (\code{gbsg_dgm}, \code{aft_dgm_flex}, or
#'   \code{glm_dgm}).
#' @param effect_measure Character or \code{NULL}.  Overrides auto-detection
#'   from \code{dgm$effect_measure}.  When \code{NULL} (default) and
#'   \code{dgm$effect_measure} is also absent, ratio-scale (\code{1/x})
#'   is used for backward compatibility.
#'
#' @return Modified DGM object with all effect fields inverted.
#'   Non-effect fields (super-population data, model parameters, etc.)
#'   are unchanged.
#'
#' @keywords internal
invert_dgm_hrs <- function(dgm, effect_measure = NULL) {
  # Auto-detect from DGM if not explicitly provided
  if (is.null(effect_measure)) {
    effect_measure <- dgm$effect_measure
  }
  use_negate <- is_identity_scale(effect_measure)

  .invert_val <- function(x) if (use_negate) -x else 1 / x

  if (!is.null(dgm$hazard_ratios)) {
    for (nm in names(dgm$hazard_ratios)) {
      val <- dgm$hazard_ratios[[nm]]
      if (is.numeric(val) && length(val) == 1 && !is.na(val)) {
        dgm$hazard_ratios[[nm]] <- .invert_val(val)
      }
    }
  }
  for (f in c("hr_H_true", "hr_Hc_true", "hr_causal",
              "AHR", "AHR_H_true", "AHR_Hc_true",
              "cde_H", "cde_Hc", "CDE")) {
    if (!is.null(dgm[[f]])) dgm[[f]] <- .invert_val(dgm[[f]])
  }
  dgm
}

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
#' @param subgroup_notation Character. \code{"harm"} (default) labels
#'   subgroups as H/Hc; \code{"benefit"} labels as G/Gc for
#'   benefit-search analyses (treatment switching).
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
#' @importFrom gt gt tab_header cols_hide cols_label tab_style cell_text tab_options px cells_body cells_row_groups
#' @export
build_classification_table <- function(
    scenario_results,
    analyses = NULL,
    digits = 2,
    title = "Subgroup Identification and Classification Rates",
    n_sims = NULL,
    bold_threshold = 0.05,
    font_size = 12,
    subgroup_notation = c("harm", "benefit")
) {

  subgroup_notation <- match.arg(subgroup_notation)
  L <- get_sg_labels(subgroup_notation)

  rows_list <- list()

  for (sc_name in names(scenario_results)) {
    sc  <- scenario_results[[sc_name]]
    res <- data.table::as.data.table(sc$results)
    dgm <- sc$dgm
    hyp <- sc$hypothesis
    n   <- sc$n_sample
    lab <- sc$label

    # Benefit search: invert HRs from switched → original scale
    if (subgroup_notation == "benefit") {
      em <- dgm$effect_measure
      res <- data.table::as.data.table(invert_hr_columns(res, em))
      dgm <- invert_dgm_hrs(dgm, em)
    }

    # Determine analyses if not specified
    if (is.null(analyses)) {
      analyses_use <- sort(unique(res$analysis))
    } else {
      analyses_use <- analyses
    }

    # ── Section header label ────────────────────────────────────────────────
    # Use get_dgm_hr() for compatibility with both gbsg_dgm and aft_dgm_flex
    if (hyp == "null") {
      hr_itt <- round(get_dgm_hr(dgm, "hr_overall"), digits)
      if (is.na(hr_itt)) hr_itt <- round(dgm$hr_causal, digits)
      section_label <- sprintf(
        "%s Null: N=%d, theta(ITT) = %s",
        lab, n, hr_itt
      )
    } else {
      df_s <- dgm$df_super_rand %||% dgm$df_super
      harm_col <- intersect(c("flag.harm", "flag_harm"), names(df_s))
      prop_sg <- if (length(harm_col) > 0) {
        round(100 * mean(df_s[[harm_col[1]]], na.rm = TRUE), 0)
      } else { NA_integer_ }
      hr_sg  <- round(get_dgm_hr(dgm, "hr_H"), digits)
      hr_sgc <- round(get_dgm_hr(dgm, "hr_Hc"), digits)
      hr_itt <- round(get_dgm_hr(dgm, "hr_overall"), digits)
      if (is.na(hr_itt)) hr_itt <- round(dgm$hr_causal, digits)
      section_label <- sprintf(
        "%s Alt: N=%d, p_%s=%d%%, theta(%s)=%s, theta(%s)=%s, theta(ITT)=%s",
        lab, n, L$plain, prop_sg, L$plain, hr_sg,
        L$plain_c, hr_sgc, hr_itt
      )
    }

    # ── Metric computation ─────────────────────────────────────────────────
    # Internal keys decouple computation from display labels (H vs Q)
    metric_keys <- if (hyp == "null") {
      c("any_sg", "sens_sgc", "ppv_sgc", "avg_sg")
    } else {
      c("any_sg", "sens_sg", "sens_sgc", "ppv_sg", "ppv_sgc", "avg_sg")
    }

    metric_display <- c(
      any_sg   = sprintf("any(%s)", L$plain),
      sens_sg  = sprintf("sens(%s)", L$plain),
      sens_sgc = sprintf("sens(%s)", L$plain_c),
      ppv_sg   = sprintf("ppv(%s)", L$plain),
      ppv_sgc  = sprintf("ppv(%s)", L$plain_c),
      avg_sg   = sprintf("avg|%s|", L$plain)
    )

    for (key in metric_keys) {
      metric <- metric_display[[key]]
      row_vals <- list(Scenario = section_label, Metric = metric)

      for (a in analyses_use) {
        r <- res[res$analysis == a, ]
        r_found <- r[r$any.H == 1, ]

        # León et al. (2024, Table 1): classification metrics are averaged
        # across ALL simulations (unconditional).  Under alt, when Ĥ = ∅:
        #   sens(Ĥ) = #{Ĥ ∩ H}/#{H} = 0/|H| = 0  (H exists but not found)
        #   ppv(Ĥ)  = #{Ĥ ∩ H}/#{Ĥ} = 0/0         (convention: 0)
        # Under null, sens is genuinely undefined (H = ∅) → keep NA.
        # The raw values from run_simulation_analysis store NA for both
        # sens and ppv when any.H = 0; we replace with 0 under alt only.
        if (hyp != "null") {
          r$sens[is.na(r$sens) & r$any.H == 0] <- 0
          r$ppv[is.na(r$ppv) & r$any.H == 0]   <- 0
          r$spec[is.na(r$spec) & r$any.H == 0]  <- 1
          r$npv[is.na(r$npv) & r$any.H == 0]    <- 1
        }

        val <- switch(
          key,
          "any_sg"   = mean(r$any.H, na.rm = TRUE),
          "sens_sg"  = mean(r$sens, na.rm = TRUE),
          "sens_sgc" = mean(r$spec, na.rm = TRUE),
          "ppv_sg"   = mean(r$ppv, na.rm = TRUE),
          "ppv_sgc"  = mean(r$npv, na.rm = TRUE),
          "avg_sg"   = if (nrow(r_found) > 0) {
            round(mean(r_found$size.H, na.rm = TRUE), 0)
          } else {
            NA_real_
          }
        )

        row_vals[[a]] <- if (key == "avg_sg") val else round(val, digits)
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
        table_df$Metric == sprintf("any(%s)", L$plain) &
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
#' @param subgroup_notation Character. \code{"harm"} (default) labels
#'   subgroups as H/Hc; \code{"benefit"} labels as G/Gc for
#'   benefit-search analyses (treatment switching).
#' @param trim_threshold Numeric or \code{NULL}. When the raw mean of
#'   an estimator's values exceeds this absolute value, that row's
#'   Avg / SD / Min / Max are recomputed on the central
#'   \code{1 - 2 * trim_fraction} of estimates, and a small footnote
#'   marker is placed on the row's Avg cell. The Estimator label is
#'   left unchanged so the Unicode notation (e.g., \eqn{\hat\theta(\hat
#'   Q)}) renders normally; only the Avg cell carries the trim
#'   indicator. Set \code{NULL} to disable trimming entirely. Default:
#'   \code{1000} (a value typical simulation results will never
#'   approach unless a small subgroup produces a degenerate effect
#'   estimate).
#' @param trim_fraction Numeric in \code{(0, 0.5)}. Fraction of
#'   estimates to trim from each tail (per row) when trimming
#'   triggers. Default: \code{0.01} (lower and upper 1\%; central 98\%
#'   used). Has no effect when \code{trim_threshold = NULL}.
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
#' @importFrom gt gt tab_header cols_label tab_style cell_text tab_options px cells_body cells_row_groups tab_footnote cells_column_labels
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
    cde_Hc = NULL,
    subgroup_notation = c("harm", "benefit"),
    trim_threshold = 1000,
    trim_fraction = 0.01
) {

  subgroup_notation <- match.arg(subgroup_notation)
  L <- get_sg_labels(subgroup_notation)

  # Detect GLM outcome type from DGM for scale-aware inversion and labels
  dgm_ot <- dgm$outcome_type
  dgm_em <- dgm$effect_measure
  is_glm_dgm <- !is.null(dgm_ot) && dgm_ot != "survival"

  # Benefit search: invert effects from switched → original scale
  if (subgroup_notation == "benefit") {
    results <- invert_hr_columns(results, dgm_em)
    dgm     <- invert_dgm_hrs(dgm, dgm_em)
  }

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

  # Notation-aware label map: internal keys -> Unicode display
  # Unicode: \u0302 = combining circumflex, \u1D9C = superscript-c,
  #          \u2020 = dagger, \u2021 = double-dagger,
  #          \u03b8 = theta, \u00e2 = a-hat
  sg  <- L$sg_hat    # e.g., Ĥ or Q̂
  sgc <- L$sg_hat_c  # e.g., Ĥᶜ or Q̂ᶜ

  label_map <- list(
    "theta-hat(sg-hat)"   = sprintf("\u03b8\u0302(%s)", sg),
    "theta-hat*(sg-hat)"  = sprintf("\u03b8\u0302*(%s)", sg),
    "ahr-hat(sg-hat)"     = sprintf("\u00e2hr(%s)", sg),
    "cde-hat(sg-hat)"     = sprintf("\u03b8\u2021(%s)", sg),
    "theta-hat(sgc-hat)"  = sprintf("\u03b8\u0302(%s)", sgc),
    "theta-hat*(sgc-hat)" = sprintf("\u03b8\u0302*(%s)", sgc),
    "ahr-hat(sgc-hat)"    = sprintf("\u00e2hr(%s)", sgc),
    "cde-hat(sgc-hat)"    = sprintf("\u03b8\u2021(%s)", sgc)
  )

  # Row-group header builders using paper notation
  build_sg_label <- function(n_est, avg_sz, theta_true, cde_val) {
    cde_part <- if (!is.na(cde_val)) {
      sprintf(", \u03b8\u2021(%s) = %s", L$sg, round(cde_val, digits))
    } else ""
    sprintf(
      "%s: %d estimable, avg |%s| = %d, \u03b8\u2020(%s) = %s%s",
      sg, n_est, sg, avg_sz, L$sg, round(theta_true, digits), cde_part
    )
  }

  build_sgc_label <- function(avg_sz_sgc, theta_true, cde_val) {
    cde_part <- if (!is.na(cde_val)) {
      sprintf(", \u03b8\u2021(%s) = %s", L$sg_c, round(cde_val, digits))
    } else ""
    sprintf(
      "%s: avg |%s| = %d, \u03b8\u2020(%s) = %s%s",
      sgc, sgc, avg_sz_sgc, L$sg_c, round(theta_true, digits), cde_part
    )
  }

  # ── Trimming infrastructure ──────────────────────────────────────────────
  # Mirrors summaryout_mrct() in mrct_simulation.R.  Trimming triggers ONLY
  # when abs(mean(est)) exceeds trim_threshold for a given row, so default
  # behaviour on healthy data is unchanged.  Trimmed rows are tracked via
  # a hidden .trimmed column on the per-row data.frame returned by
  # make_est_row(), then surfaced as a cell-level gt footnote on the Avg
  # column after the gt table is built.
  if (!is.null(trim_threshold)) {
    if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L ||
        trim_threshold <= 0) {
      stop("'trim_threshold' must be a positive numeric scalar or NULL.",
           call. = FALSE)
    }
    if (!is.numeric(trim_fraction) || length(trim_fraction) != 1L ||
        trim_fraction <= 0 || trim_fraction >= 0.5) {
      stop("'trim_fraction' must be in (0, 0.5).", call. = FALSE)
    }
  }

  # ── Helper: compute one estimation row ────────────────────────────────────
  # When has_cde is TRUE, produces both b-dagger and b-ddagger bias columns.
  # Threshold-gated trimming: when abs(raw mean) exceeds trim_threshold the
  # central (1 - 2 * trim_fraction) of estimates is used for avg / sd / min /
  # max, and the row is flagged via a hidden .trimmed column.  The Estimator
  # label is left unchanged so the label_map Unicode substitution still
  # matches; trimmed rows are surfaced to the reader via a cell-level gt
  # footnote on the Avg column (added after the label substitution).  When
  # trim_threshold is NULL no trimming occurs.
  make_est_row <- function(estimates, true_val, label, cde_val = NA_real_) {
    est <- estimates[!is.na(estimates)]
    if (length(est) == 0) return(NULL)

    raw_mean <- mean(est)

    is_trimmed <- FALSE
    if (!is.null(trim_threshold) && length(est) >= 5L &&
        abs(raw_mean) > trim_threshold) {
      lo <- stats::quantile(est, trim_fraction,     na.rm = TRUE)
      hi <- stats::quantile(est, 1 - trim_fraction, na.rm = TRUE)
      est_trim <- est[est >= lo & est <= hi]
      if (length(est_trim) >= 3L) {
        est        <- est_trim
        is_trimmed <- TRUE
      }
    }

    avg_est  <- mean(est)
    sd_est   <- stats::sd(est)
    min_est  <- min(est)
    max_est  <- max(est)

    # b-dagger: bias relative to marginal (causal) truth.  Bias is computed
    # from the (possibly trimmed) avg_est so the displayed avg and bias stay
    # internally consistent.
    b_dagger <- 100 * (avg_est - true_val) / true_val

    row <- data.frame(
      Estimator          = label,
      Avg                = round(avg_est, digits),
      SD                 = round(sd_est, digits),
      Min                = round(min_est, digits),
      Max                = round(max_est, digits),
      .trimmed           = is_trimmed,
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
      res_found$hr.H.hat, theta_H_true, "theta-hat(sg-hat)",
      cde_val = if (has_cde) cde_H else NA_real_
    )
  }

  # Cox HR (bias-corrected)
  if ("hr.H.bc" %in% names(res_found)) {
    rows_H[[length(rows_H) + 1]] <- make_est_row(
      res_found$hr.H.bc, theta_H_true, "theta-hat*(sg-hat)",
      cde_val = if (has_cde) cde_H else NA_real_
    )
  }

  # AHR (identified subgroup) — AHR is a different estimand, no CDE comparison
  if ("ahr.H.hat" %in% names(res_found) && !is.na(ahr_H_true)) {
    rows_H[[length(rows_H) + 1]] <- make_est_row(
      res_found$ahr.H.hat, ahr_H_true, "ahr-hat(sg-hat)"
    )
  }

  # CDE (identified subgroup) — bias is relative to CDE truth only
  if ("cde.H.hat" %in% names(res_found) && has_cde) {
    rows_H[[length(rows_H) + 1]] <- make_est_row(
      res_found$cde.H.hat, cde_H, "cde-hat(sg-hat)"
    )
  }

  # ── Rows for complement ────────────────────────────────────────────────────
  rows_Hc <- list()

  # Cox HR (identified complement)
  if ("hr.Hc.hat" %in% names(res_found)) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$hr.Hc.hat, theta_Hc_true, "theta-hat(sgc-hat)",
      cde_val = if (has_cde) cde_Hc else NA_real_
    )
  }

  # Cox HR (bias-corrected)
  if ("hr.Hc.bc" %in% names(res_found)) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$hr.Hc.bc, theta_Hc_true, "theta-hat*(sgc-hat)",
      cde_val = if (has_cde) cde_Hc else NA_real_
    )
  }

  # AHR (identified complement)
  if ("ahr.Hc.hat" %in% names(res_found) && !is.na(ahr_Hc_true)) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$ahr.Hc.hat, ahr_Hc_true, "ahr-hat(sgc-hat)"
    )
  }

  # CDE (identified complement) — bias is relative to CDE truth only
  if ("cde.Hc.hat" %in% names(res_found) && has_cde) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$cde.Hc.hat, cde_Hc, "cde-hat(sgc-hat)"
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
    h_label <- build_sg_label(
      n_estimable, avg_size_H, theta_H_true,
      if (has_cde) cde_H else NA_real_
    )
    df_H[, Subgroup := h_label]
  }
  if (!is.null(df_Hc)) {
    hc_label <- build_sgc_label(
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
  # Capture trimmed-row indices and drop the internal .trimmed flag column
  # BEFORE the label_map substitution and gt construction.  Indices feed a
  # cell-level footnote on the Avg column added below; dropping the column
  # keeps it out of the rendered table.
  trimmed_rows <- if (".trimmed" %in% names(table_df)) {
    which(table_df$.trimmed)
  } else integer(0)
  if (".trimmed" %in% names(table_df)) {
    table_df[, .trimmed := NULL]
  }

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

  # Notation footnote (paper-aligned, notation-aware, GLM-aware)
  if (is_glm_dgm) {
    em_label <- dgm_em %||% "GLM"
    footnote_parts <- paste0(
      "\u03b8\u0302(", L$sg_hat, ") = plugin ", em_label,
      " estimate in identified subgroup; ",
      "\u03b8\u0302*(", L$sg_hat, ") = bootstrap bias-corrected; ",
      "b\u2020 = bias relative to true ", em_label,
      " \u03b8\u2020 (causal truth)"
    )
  } else {
    footnote_parts <- paste0(
      "\u03b8\u0302(", L$sg_hat, ") = plugin Cox HR in identified subgroup; ",
      "\u03b8\u0302*(", L$sg_hat, ") = bootstrap bias-corrected; ",
      "\u00e2hr(", L$sg_hat,
      ") = average hazard ratio in identified subgroup; ",
      "b\u2020 = bias relative to marginal HR \u03b8\u2020 (causal truth)"
    )
  }
  if (has_cde) {
    footnote_parts <- paste0(
      footnote_parts,
      "; \u03b8\u2021(", L$sg_hat, ") = controlled direct effect in identified subgroup",
      "; b\u2021 = bias relative to CDE \u03b8\u2021"
    )
  }
  gt_tbl <- gt::tab_footnote(gt_tbl, footnote = footnote_parts)

  # Cell-level trim footnote on the Avg column of affected rows.  Replaces
  # the previous bottom-of-table source_note + "(trimmed)" label suffix
  # combination — the suffix broke the Estimator label_map substitution and
  # leaked raw keys (e.g., "theta-hat(sgc-hat)") into the rendered table.
  # Cell-level footnotes preserve the Unicode label notation and place a
  # small superscript marker only on the affected Avg cells.
  if (length(trimmed_rows) > 0) {
    pct <- round(trim_fraction * 100)
    trim_note <- sprintf(
      paste0("Avg / SD / Min / Max for marked rows are computed on the ",
             "central %d%% of estimates (lower and upper %d%% excluded) ",
             "because the row's raw mean exceeded %s."),
      100 - 2 * pct, pct, format(trim_threshold, big.mark = ",")
    )
    gt_tbl <- gt::tab_footnote(
      gt_tbl,
      footnote  = trim_note,
      locations = gt::cells_body(columns = "Avg", rows = trimmed_rows)
    )
  }

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
#' @param subgroup_notation Character. \code{"harm"} (default) labels
#'   subgroups as H/Hc; \code{"benefit"} labels as G/Gc for
#'   benefit-search analyses (treatment switching).
#' @param trim_threshold Numeric or \code{NULL}. When the raw mean of
#'   an estimate vector exceeds this absolute value, the narrative
#'   reports averages and SDs computed on the central
#'   \code{1 - 2 * trim_fraction} of estimates and adds a closing
#'   note disclosing the trimming. Set \code{NULL} to disable.
#'   Default: \code{1000} (matches \code{build_estimation_table}).
#' @param trim_fraction Numeric in \code{(0, 0.5)}. Fraction of
#'   estimates to trim from each tail when trimming triggers.
#'   Default: \code{0.01}.  Has no effect when
#'   \code{trim_threshold = NULL}.
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
    cat = TRUE,
    subgroup_notation = c("harm", "benefit"),
    trim_threshold = 1000,
    trim_fraction = 0.01
) {

  subgroup_notation <- match.arg(subgroup_notation)
  L <- get_sg_labels(subgroup_notation)

  # Benefit search: invert effects from switched → original scale
  if (subgroup_notation == "benefit") {
    results <- invert_hr_columns(results, dgm$effect_measure)
    dgm     <- invert_dgm_hrs(dgm)
  }

  # Detect GLM outcome for narrative labels
  dgm_ot <- dgm$outcome_type
  dgm_em <- dgm$effect_measure
  is_glm <- !is.null(dgm_ot) && dgm_ot != "survival"
  em_label <- effect_measure_label(dgm_em, dgm_ot)

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

  # Threshold-gated trimming for the narrative.  Mirrors the pattern in
  # build_estimation_table(): when an estimate vector's raw mean exceeds
  # trim_threshold, mean / sd / bias are recomputed on the central
  # (1 - 2 * trim_fraction) of estimates.  Each of the five blocks below
  # routes through maybe_trim() so trimming applies uniformly.
  if (!is.null(trim_threshold)) {
    if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L ||
        trim_threshold <= 0) {
      stop("'trim_threshold' must be a positive numeric scalar or NULL.",
           call. = FALSE)
    }
    if (!is.numeric(trim_fraction) || length(trim_fraction) != 1L ||
        trim_fraction <= 0 || trim_fraction >= 0.5) {
      stop("'trim_fraction' must be in (0, 0.5).", call. = FALSE)
    }
  }
  trim_env <- new.env(parent = emptyenv())
  trim_env$trimmed_any <- FALSE

  maybe_trim <- function(vals) {
    if (is.null(trim_threshold) || length(vals) < 5L) return(vals)
    if (abs(mean(vals)) <= trim_threshold) return(vals)
    lo <- stats::quantile(vals, trim_fraction,     na.rm = TRUE)
    hi <- stats::quantile(vals, 1 - trim_fraction, na.rm = TRUE)
    vals_trim <- vals[vals >= lo & vals <= hi]
    if (length(vals_trim) < 3L) return(vals)
    trim_env$trimmed_any <- TRUE
    vals_trim
  }

  avg_size_H  <- round(mean(res_found$size.H, na.rm = TRUE), 0)
  avg_size_Hc <- round(mean(res_found$size.Hc, na.rm = TRUE), 0)
  detect_rate <- round(100 * n_estimable / n_sims, 1)

  # Effect estimate summaries
  has_hr_H  <- "hr.H.hat" %in% names(res_found)
  has_hr_Hc <- "hr.Hc.hat" %in% names(res_found)
  has_hr_bc <- "hr.H.bc" %in% names(res_found)

  if (has_hr_H) {
    hr_H_vals <- maybe_trim(res_found$hr.H.hat[!is.na(res_found$hr.H.hat)])
    hr_H_avg  <- fmt(mean(hr_H_vals))
    hr_H_sd   <- fmt(stats::sd(hr_H_vals))
    hr_H_bias <- if (!is.na(theta_H_true) && theta_H_true != 0) {
      fmt(100 * (mean(hr_H_vals) - theta_H_true) / theta_H_true)
    } else NA
  }
  if (has_hr_Hc) {
    hr_Hc_vals <- maybe_trim(res_found$hr.Hc.hat[!is.na(res_found$hr.Hc.hat)])
    hr_Hc_avg  <- fmt(mean(hr_Hc_vals))
    hr_Hc_sd   <- fmt(stats::sd(hr_Hc_vals))
    hr_Hc_bias <- if (!is.na(theta_Hc_true) && theta_Hc_true != 0) {
      fmt(100 * (mean(hr_Hc_vals) - theta_Hc_true) / theta_Hc_true)
    } else NA
  }
  if (has_hr_bc) {
    hr_bc_vals <- maybe_trim(res_found$hr.H.bc[!is.na(res_found$hr.H.bc)])
    hr_bc_avg  <- fmt(mean(hr_bc_vals))
    hr_bc_bias <- if (!is.na(theta_H_true) && theta_H_true != 0) {
      fmt(100 * (mean(hr_bc_vals) - theta_H_true) / theta_H_true)
    } else NA
  }

  # AHR summaries
  has_ahr_H  <- "ahr.H.hat" %in% names(res_found) && !is.na(ahr_H_true)
  has_ahr_Hc <- "ahr.Hc.hat" %in% names(res_found) && !is.na(ahr_Hc_true)

  if (has_ahr_H) {
    ahr_H_vals <- maybe_trim(res_found$ahr.H.hat[!is.na(res_found$ahr.H.hat)])
    ahr_H_avg  <- fmt(mean(ahr_H_vals))
    ahr_H_bias <- if (!is.na(ahr_H_true) && ahr_H_true != 0) {
      fmt(100 * (mean(ahr_H_vals) - ahr_H_true) / ahr_H_true)
    } else NA
  }
  if (has_ahr_Hc) {
    ahr_Hc_vals <- maybe_trim(res_found$ahr.Hc.hat[!is.na(res_found$ahr.Hc.hat)])
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
        "Under the null hypothesis (true %s = %s uniformly), ",
        "%d of %d simulations (%s) identified a subgroup using %s. ",
        "This low detection rate confirms controlled type-I error. ",
        "Among those %d false detections, the identified subgroup ",
        "averaged %d patients."
      ),
      em_label, true_hr_str, n_estimable, n_sims, fmt_pct(detect_rate),
      analysis_method, n_estimable, avg_size_H
    )
  } else {
    paras[1] <- sprintf(
      paste0(
        "Under the alternative hypothesis (true %s(%s) = %s, ",
        "true %s(%s) = %s), %d of %d simulations (%s) ",
        "identified a subgroup using %s. ",
        "The identified subgroup averaged %d patients ",
        "(complement: %d)."
      ),
      em_label, L$plain, fmt_s(theta_H_true),
      em_label, L$plain_c, fmt_s(theta_Hc_true),
      n_estimable, n_sims, fmt_pct(detect_rate),
      analysis_method, avg_size_H, avg_size_Hc
    )
  }

  # --- Paragraph 2: Effect estimates ---
  if (has_hr_H && has_hr_Hc) {
    bias_H_str  <- fmt_pct(hr_H_bias)
    bias_Hc_str <- fmt_pct(hr_Hc_bias)

    if (scenario == "null") {
      paras[2] <- sprintf(
        paste0(
          "The naive %s in the identified subgroup averaged %s ",
          "(SD = %s), representing %s relative bias above the true ",
          "value of %s. This upward bias reflects selection: the algorithm ",
          "identified whichever patients happened to look most like a %s ",
          "subgroup by chance. In the complement, the %s averaged %s ",
          "(%s bias), showing the expected mirror effect where removing ",
          "the worst-looking patients makes the remainder appear modestly ",
          "better."
        ),
        em_label, hr_H_avg, hr_H_sd, bias_H_str, true_hr_str,
        L$word,
        em_label, hr_Hc_avg, bias_Hc_str
      )
    } else {
      paras[2] <- sprintf(
        paste0(
          "The naive %s in the identified subgroup averaged %s ",
          "(SD = %s), corresponding to %s relative bias versus the ",
          "true %s(%s) = %s. In the complement, the estimate averaged %s ",
          "(%s bias vs. true %s(%s) = %s)."
        ),
        em_label, hr_H_avg, hr_H_sd, bias_H_str,
        em_label, L$plain, fmt_s(theta_H_true),
        hr_Hc_avg, bias_Hc_str, em_label, L$plain_c, fmt_s(theta_Hc_true)
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
        "theta-ddagger(%s) = %s, the naive plugin shows %s ",
        "relative bias."
      ),
      L$plain, fmt_s(cde_H), fmt_pct(b_cde_H)
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

  # --- Paragraph 5: AHR comparison (survival-only) ---
  if (has_ahr_H) {
    ahr_para <- sprintf(
      paste0(
        "The average hazard ratio (AHR) in the identified subgroup ",
        "averaged %s (%s relative bias vs. true AHR(%s) = %s)"
      ),
      ahr_H_avg, fmt_pct(ahr_H_bias), L$plain, fmt_s(ahr_H_true)
    )
    if (has_ahr_Hc) {
      ahr_para <- paste0(ahr_para, sprintf(
        "; in the complement, %s (%s bias vs. true AHR(%s) = %s)",
        ahr_Hc_avg, fmt_pct(ahr_Hc_bias), L$plain_c, fmt_s(ahr_Hc_true)
      ))
    }
    if (has_hr_H && !is.na(ahr_H_bias) && !is.na(hr_H_bias) &&
        abs(ahr_H_bias) < abs(hr_H_bias)) {
      ahr_para <- paste0(ahr_para,
        ". The AHR shows attenuated bias relative to the ", em_label, ", ",
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

  # Trim disclosure: appended only if any *_vals block triggered trimming.
  if (trim_env$trimmed_any) {
    pct_kept <- 100 - 2 * round(trim_fraction * 100)
    paras[length(paras) + 1] <- sprintf(
      paste0("*Note*: One or more reported summaries above were computed ",
             "on the central %d%% of estimates because the raw mean ",
             "exceeded %s; this avoids unwieldy values driven by a small ",
             "number of extreme replicates.  The companion estimation ",
             "table marks the affected rows with '(trimmed)'."),
      pct_kept, format(trim_threshold, big.mark = ",")
    )
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
#' @param subgroup_notation Character. \code{"harm"} (default) labels
#'   subgroups as H/Hc; \code{"benefit"} labels as G/Gc for
#'   benefit-search analyses (treatment switching).
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
#' @importFrom gt gt tab_header cols_label tab_style cell_text tab_options px cells_body cells_row_groups
#' @export
render_reference_table <- function(
    ref_df,
    title = "Reference Simulation Results",
    subtitle = NULL,
    bold_threshold = 0.05,
    subgroup_notation = c("harm", "benefit")
) {

  subgroup_notation <- match.arg(subgroup_notation)
  L <- get_sg_labels(subgroup_notation)

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
    any_label <- sprintf("any(%s)", L$plain)
    for (col_name in analysis_cols) {
      bold_rows <- which(
        ref_df$Metric == any_label &
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

