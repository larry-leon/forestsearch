# Required packages for parallel bootstrap workers
BOOTSTRAP_REQUIRED_PACKAGES <- c(
  "data.table",    # Data manipulation in parallel workers
  "foreach",       # Foreach iteration framework
  "doFuture",      # Future-based parallel backend
  "doRNG",         # Reproducible random numbers in parallel
  "survival"       # Cox model fitting for bias correction
)

#' Resolve parallel processing arguments for bootstrap
#'
#' If parallel_args not provided, falls back to forestsearch call's
#' parallel configuration. Always reports configuration to user.
#'
#' @param parallel_args List or empty list
#' @param forestsearch_call_args List from original forestsearch call
#' @return List with plan, workers, show_message
#' @keywords internal

resolve_bootstrap_parallel_args <- function(parallel_args, forestsearch_call_args) {
  # Use provided args if non-empty, otherwise inherit from forestsearch call
  if (length(parallel_args) == 0) {
    resolved_args <- as.list(forestsearch_call_args$parallel_args)
    message("Bootstrap parallel config: Using 'observed' data analysis forestsearch settings")
  } else {
    resolved_args <- parallel_args
    message("Bootstrap parallel config: Using user-provided settings")
  }

  # Report configuration to user
  max_cores <- parallel::detectCores()
  message("System max cores available: ", max_cores)
  message("Bootstrap will use: ", resolved_args$workers, " workers with '",
          resolved_args$plan, "' plan")

  resolved_args
}


#' Functions required in parallel bootstrap environment
#'
#' Organized by functional category for maintainability
BOOTSTRAP_REQUIRED_FUNCTIONS <- list(
  statistics = c(
    "calc_cov",
    "calculate_counts",
    "analyze_subgroups",
    "calculate_potential_hr",
    "ci.est",
    "count.id",
    "CV_sgs",
    "get_targetEst",
    "get_dfRes",
    "getCIs",
    "getci_Cox",
    "SummaryStat",
    "var_summary",
    "format_results",
    "format_CI",
    "qlow", "qhigh"
  ),
  subgroup_analysis = c(
    "cox_summary",
    "km_summary",
    "n_pcnt",
    "plot_subgroup",
    "plot_weighted_km",
    "prepare_subgroup_data",
    "quiet",
    "rmst_calculation",
    "sg_tables",
    "sort_subgroups",
    "remove_redundant_subgroups",
    "sg_consistency_out",
    "get_split_hr"
  ),
  data_prep = c(
    "get_FSdata",
    "dummy", "dummy2",
    "acm.disjctif", "acm.util.df2", "acm.util.df",
    "ztrail", "one.zero",
    "prepare_data",
    "clean_data"
  ),
  forestsearch_core = c(
    "forestsearch",
    "forestsearch_bootstrap_dofuture",
    "run_bootstrap",
    "run_grf",
    "evaluate_subgroups",
    "summarize_results",
    "SG_tab_estimates",
    "get_dfpred",
    "get_combinations_info",
    "get_subgroup_membership",
    "get_covs_in",
    "extract_idx_flagredundancy",
    "get_cut_name",
    "cut_var",
    "thiscut",
    "FS_labels"
  ),
  grf_policy = c(
    "grf.subg.harm.survival",
    "grf.estimates.out",
    "policy_tree",
    "causal_survival_forest"
  ),
  search_consistency = c(
    "subgroup.search",
    "subgroup.consistency",
    "extract_subgroup",
    "filter_by_lassokeep",
    "is.continuous",
    "process_conf_force_expr",
    "is_flag_continuous",
    "is_flag_drop",
    "lasso_selection",
    "get_Cox_sg",
    "get_conf_force"
  ),
  bootstrap_parallel = c(
    "bootstrap_results",
    "bootstrap_ystar",
    "ensure_packages",
    "fit_cox_models",
    "build_cox_formula",
    "cox.formula.boot",
    "setup_parallel_SGcons"
  )
)

#' Flatten required functions list to character vector
#'
#' @return Character vector of all function names needed in parallel workers
#' @keywords internal
get_bootstrap_exports <- function() {
  unlist(BOOTSTRAP_REQUIRED_FUNCTIONS, use.names = FALSE)
}

#' Find integer pairs (x, y) such that x * y = z and y >= x
#'
#' Given an integer z, this function finds all integer pairs (x, y) such that x * y = z and y >= x.
#' Optionally, you can return only the pair with the largest value of x or y.
#'
#' @param z Integer. The target product.
#' @param return_largest Character. If \"x\", returns the pair with the largest x. If \"y\", returns the pair with the largest y. If NULL, returns all pairs.
#'
#' @return A matrix of integer pairs (x, y) satisfying the conditions, or a single pair if return_largest is specified.
#' @examples
#' find_xy_given_z(12)
#' find_xy_given_z(12, return_largest = \"x\")
#' find_xy_given_z(12, return_largest = \"y\")
#' @keywords internal

find_xy_given_z <- function(z, return_largest = NULL) {
  pairs <- list()
  for (x in 1:z) {
    if (z %% x == 0) {
      y <- z / x
      if (y >= x && y %% 1 == 0) {
        pairs[[length(pairs) + 1]] <- c(x, y)
      }
    }
  }
  result <- do.call(rbind, pairs)
  if (!is.null(return_largest)) {
    if (return_largest == "x") {
      idx <- which.max(result[,1])
      return(result[idx, , drop = FALSE])
    } else if (return_largest == "y") {
      idx <- which.max(result[,2])
      return(result[idx, , drop = FALSE])
    }
  }
  return(result)
}


#' Ensure Required Packages Are Installed and Loaded
#'
#' Installs and loads required packages if not already available.
#'
#' @param pkgs Character vector of package names.
#' @return None. Packages are loaded into the session.
#' @export

ensure_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

#' Format Confidence Interval for Estimates
#'
#' Formats confidence interval for estimates.
#'
#' @param estimates Data frame or data.table of estimates.
#' @param col_names Character vector of column names for estimate, lower, upper.
#' @return Character string formatted as \"estimate (lower, upper)\".
#' @export

format_CI <- function(estimates, col_names) {
  resH <- estimates[, ..col_names]
  Hstat <- round(unlist(resH[1, ]), 2)
  paste0(Hstat[1], " (", Hstat[2], ",", Hstat[3], ")")
}


#' Bootstrap Results for ForestSearch (legacy, in case new one doesn't work out)
#'
#' Runs bootstrap analysis for ForestSearch, fitting Cox models and bias correction.
#'
#' @param fs.est ForestSearch results object.
#' @param df_boot_analysis Data frame for bootstrap analysis.
#' @param cox.formula.boot Cox model formula for bootstrap.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @param show_three Logical. Show details for first three bootstraps.
#' @param H_obs Numeric. Observed HR for subgroup H.
#' @param Hc_obs Numeric. Observed HR for subgroup Hc.
#' @param reset_parallel Logical. Reset parallel plan for bootstrap.
#' @param boot_workers Integer. Number of parallel workers.
#' @return Data.table with bias-adjusted estimates and search metrics.
#' @importFrom foreach foreach
#' @importFrom data.table data.table
#' @importFrom doFuture %dofuture%
#' @export

bootstrap_results_legacy <- function(fs.est, df_boot_analysis, cox.formula.boot, nb_boots, show_three, H_obs, Hc_obs) {

  NN <- nrow(df_boot_analysis)
  id0 <- seq_len(NN)
  foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE,
                           add = get_bootstrap_exports()
    ),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {
    show3 <- FALSE
    if (show_three) show3 <- (boot <= 3)
    set.seed(8316951 + boot * 100)
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df_boot_analysis[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))

    # Bootstrap data evaluated at H: H_star
    fitH_star <- get_Cox_sg(df_sg = subset(df_boot, treat.recommend == 0), cox.formula = cox.formula.boot, est.loghr = TRUE)
    H_star <- fitH_star$est_obs
    fitHc_star <- get_Cox_sg(df_sg = subset(df_boot, treat.recommend == 1), cox.formula = cox.formula.boot, est.loghr = TRUE)
    Hc_star <- fitHc_star$est_obs

    # Bias corrections
    H_biasadj_1 <- H_biasadj_2 <- NA
    Hc_biasadj_1 <- Hc_biasadj_2 <- NA
    tmins_search <- NA
    max_sg_est <- NA
    prop_maxk <- NA
    L <- NA
    max_count <- NA

    # Drop initial confounders
    drop.vars <- c(fs.est$confounders.candidate, "treat.recommend")
    dfnew <- df_boot_analysis[, !(names(df_boot_analysis) %in% drop.vars)]
    dfnew_boot <- df_boot[, !(names(df_boot) %in% drop.vars)]
    # Extract arguments in forestsearch (observed) data analysis
    args_FS_boot <- fs.est$args_call_all
    args_FS_boot$df.analysis <- dfnew_boot
    args_FS_boot$df.predict <- dfnew
    # CATEGORY 1: OUTPUT SUPPRESSION (parallelization efficiency)
    args_FS_boot$details <- show3                # Only show first 3 for debugging
    args_FS_boot$showten_subgroups <- FALSE      # Suppress large output
    args_FS_boot$plot.sg <- FALSE                # Can't plot from parallel worker
    args_FS_boot$plot.grf <- FALSE               # Can't plot from parallel worker
    # CATEGORY 2: VARIABLE RE-SELECTION (bootstrap variability)
    # Force fresh GRF run on bootstrap (oracle uses predictions from bootstrap)
    args_FS_boot$grf_res <- NULL
    args_FS_boot$grf_cuts <- NULL
    # CATEGORY 3: SEQUENTIAL EXECUTION (prevent nested parallelization)
    # Each bootstrap is already running in a parallel worker
    # Nested parallelization causes resource contention and deadlocks
    args_FS_boot$parallel_args$plan <- "sequential"
    args_FS_boot$parallel_args$workers <- 1L
    args_FS_boot$parallel_args$show_message <- FALSE

    run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)

    if (inherits(run_bootstrap, "try-error")) {
      warning("Bootstrap ", boot, " failed: ", as.character(run_bootstrap))
    }


    if (!inherits(run_bootstrap, "try-error") && !is.null(run_bootstrap$sg.harm)) {
      df_PredBoot <- run_bootstrap$df.predict
      dfboot_PredBoot <- run_bootstrap$df.est
      max_sg_est <- as.numeric(run_bootstrap$find.grps$max_sg_est)
      tmins_search <- as.numeric(run_bootstrap$find.grps$time_search)
      prop_maxk <- as.numeric(run_bootstrap$prop_maxk)
      max_count <- run_bootstrap$find.grps$max_count
      L <- run_bootstrap$find.grps$L
      fitHstar_obs <- get_Cox_sg(df_sg = subset(df_PredBoot, treat.recommend == 0), cox.formula = cox.formula.boot, est.loghr = TRUE)
      Hstar_obs <- fitHstar_obs$est_obs
      fitHstar_star <- get_Cox_sg(df_sg = subset(dfboot_PredBoot, treat.recommend == 0), cox.formula = cox.formula.boot, est.loghr = TRUE)
      Hstar_star <- fitHstar_star$est_obs
      rm(fitHstar_star)
      H_biasadj_1 <- H_obs - (Hstar_star - Hstar_obs)
      H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)
      fitHcstar_obs <- get_Cox_sg(df_sg = subset(df_PredBoot, treat.recommend == 1), cox.formula = cox.formula.boot, est.loghr = TRUE)
      Hcstar_obs <- fitHcstar_obs$est_obs
      rm(fitHcstar_obs)
      fitHcstar_star <- get_Cox_sg(df_sg = subset(dfboot_PredBoot, treat.recommend == 1), cox.formula = cox.formula.boot, est.loghr = TRUE)
      Hcstar_star <- fitHcstar_star$est_obs
      Hc_biasadj_1 <- Hc_obs - (Hcstar_star - Hcstar_obs)
      Hc_biasadj_2 <- 2 * Hc_obs - (Hc_star + Hcstar_star - Hcstar_obs)
    }
    dfres <- data.table::data.table(H_biasadj_1, H_biasadj_2,
                                    Hc_biasadj_1, Hc_biasadj_2,
                                    tmins_search, max_sg_est, prop_maxk, L, max_count)
    return(dfres)
  }

}

#' Format Bootstrap Results Table with gt
#'
#' Creates a publication-ready table from ForestSearch bootstrap results,
#' with bias-corrected confidence intervals and informative formatting.
#'
#' @param FSsg_tab Data frame or matrix from forestsearch_bootstrap_dofuture()$FSsg_tab
#' @param nb_boots Integer. Number of bootstrap iterations performed
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#' @param boot_success_rate Numeric. Proportion of bootstraps that found subgroups
#' @param title Character. Custom title (optional)
#' @param subtitle Character. Custom subtitle (optional)
#'
#' @return A gt table object
#' @importFrom gt gt tab_header tab_spanner tab_footnote md cols_label
#' @export

format_bootstrap_table <- function(FSsg_tab, nb_boots, est.scale = "hr",
                                   boot_success_rate = NULL,
                                   title = NULL, subtitle = NULL) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table formatting. Install with: install.packages('gt')")
  }

  # CRITICAL FIX: Convert matrix to data frame if needed
  if (is.matrix(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  # Ensure it's a data frame
  if (!is.data.frame(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  # Default title and subtitle
  if (is.null(title)) {
    title <- "Treatment Effect by Subgroup"
  }

  if (is.null(subtitle)) {
    effect_label <- ifelse(est.scale == "hr", "Hazard Ratio", "Inverse Hazard Ratio (1/HR)")
    subtitle <- sprintf("Bootstrap bias-corrected estimates (%d iterations)", nb_boots)
  }

  # Get column names (handle different possible names)
  col_names <- colnames(FSsg_tab)

  # Create base labels list
  labels_list <- list(
    Subgroup = "Subgroup"
  )

  # Add labels for columns that exist
  if ("n" %in% col_names) labels_list$n <- "N"
  if ("n1" %in% col_names) labels_list$n1 <- gt::md("N<sub>treat</sub>")
  if ("events" %in% col_names) labels_list$events <- "Events"
  if ("m1" %in% col_names) labels_list$m1 <- gt::md("Med<sub>treat</sub>")
  if ("m0" %in% col_names) labels_list$m0 <- gt::md("Med<sub>ctrl</sub>")
  if ("RMST" %in% col_names) labels_list$RMST <- gt::md("RMST<sub>diff</sub>")

  # Handle HR column (might be "HR (95% CI)" or similar)
  hr_col <- grep("HR.*CI", col_names, value = TRUE)[1]
  if (!is.na(hr_col) && length(hr_col) > 0) {
    labels_list[[hr_col]] <- gt::md("HR<br/>(95% CI)<sup>†</sup>")
  }

  # Handle adjusted HR column (might be "HR*" or similar)
  hr_adj_col <- grep("HR\\*", col_names, value = TRUE)[1]
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    labels_list[[hr_adj_col]] <- gt::md("HR<sup>‡</sup><br/>(95% CI)")
  }

  # Create the gt table
  tbl <- FSsg_tab |>
    gt::gt() |>

    # Title and subtitle
    gt::tab_header(
      title = gt::md(paste0("**", title, "**")),
      subtitle = subtitle
    ) |>

    # Column labels with better formatting
    gt::cols_label(.list = labels_list)

  # Add spanners if columns exist
  sample_size_cols <- intersect(c("n", "n1", "events"), col_names)
  if (length(sample_size_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Sample Size",
        columns = sample_size_cols
      )
  }

  survival_cols <- intersect(c("m1", "m0", "RMST"), col_names)
  if (length(survival_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Survival",
        columns = survival_cols
      )
  }

  hr_cols <- grep("HR", col_names, value = TRUE)
  if (length(hr_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Treatment Effect",
        columns = gt::starts_with("HR")
      )
  }

  # Add footnotes
  if (!is.na(hr_col) && length(hr_col) > 0) {
    # Use dplyr::all_of() which comes with gt, or just use the column directly
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md("**Unadjusted HR**: Standard Cox regression hazard ratio with robust standard errors"),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_col))
      )
  }

  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md(sprintf(
          "**Bias-corrected HR**: Bootstrap-adjusted estimate using infinitesimal jacknife method (%d iterations). Corrects for optimism in subgroup selection.",
          nb_boots
        )),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_adj_col))
      )
  }

  # Source note with technical details
  source_note_text <- "*Note*: Med = Median survival time (months). RMST<sub>diff</sub> = Restricted mean survival time difference."
  if (!is.null(boot_success_rate)) {
    source_note_text <- paste0(
      source_note_text,
      sprintf(" Subgroup identified in %.1f%% of bootstrap samples.", boot_success_rate * 100)
    )
  }

  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(source_note_text))

  # Styling
  tbl <- tbl |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f0f0f0"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_column_labels()
    )

  # Center align numeric columns if they exist
  numeric_cols <- intersect(c("n", "n1", "events", "m1", "m0", "RMST"), col_names)
  if (length(numeric_cols) > 0) {
    tbl <- tbl |>
      gt::tab_style(
        style = gt::cell_text(align = "center"),
        locations = gt::cells_body(columns = dplyr::all_of(numeric_cols))
      )
  }

  # Highlight the bias-corrected estimates if column exists
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    tbl <- tbl |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(columns = dplyr::all_of(hr_adj_col))
      )
  }

  # Add subtle borders
  tbl <- tbl |>
    gt::tab_options(
      table.border.top.style = "solid",
      table.border.top.width = gt::px(3),
      table.border.top.color = "#333333",
      heading.border.bottom.style = "solid",
      heading.border.bottom.width = gt::px(2),
      heading.border.bottom.color = "#333333",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = gt::px(2),
      column_labels.border.bottom.color = "#333333",
      table.font.size = gt::px(14)
    )

  return(tbl)
}


#' Enhanced Bootstrap Results Summary
#'
#' Creates comprehensive output including formatted table, diagnostic plots,
#' and bootstrap quality metrics.
#'
#' @param boot_results List. Output from forestsearch_bootstrap_dofuture()
#' @param create_plots Logical. Generate diagnostic plots (default: TRUE)
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#'
#' @return List with formatted table and diagnostics
#' @export

summarize_bootstrap_results <- function(boot_results, create_plots = TRUE,
                                        est.scale = "hr") {

  # Extract components
  FSsg_tab <- boot_results$FSsg_tab
  results <- boot_results$results
  H_estimates <- boot_results$H_estimates
  Hc_estimates <- boot_results$Hc_estimates

  # Calculate bootstrap success rate
  boot_success_rate <- mean(!is.na(results$H_biasadj_2))
  nb_boots <- nrow(results)

  # Create formatted table
  formatted_table <- format_bootstrap_table(
    FSsg_tab = FSsg_tab,
    nb_boots = nb_boots,
    est.scale = est.scale,
    boot_success_rate = boot_success_rate
  )

  # Bootstrap diagnostics
  diagnostics <- list(
    n_boots = nb_boots,
    success_rate = boot_success_rate,
    n_successful = sum(!is.na(results$H_biasadj_2)),
    n_failed = sum(is.na(results$H_biasadj_2)),
    median_search_time = median(results$tmins_search, na.rm = TRUE),
    total_search_time = sum(results$tmins_search, na.rm = TRUE)
  )

  # Print summary
  cat("\n=== Bootstrap Analysis Summary ===\n")
  cat(sprintf("Total bootstrap iterations: %d\n", diagnostics$n_boots))
  cat(sprintf("Successful subgroup identification: %d (%.1f%%)\n",
              diagnostics$n_successful, diagnostics$success_rate * 100))
  cat(sprintf("Failed to find subgroup: %d (%.1f%%)\n",
              diagnostics$n_failed, (1 - diagnostics$success_rate) * 100))
  cat(sprintf("Median search time per bootstrap: %.2f minutes\n",
              diagnostics$median_search_time))
  cat(sprintf("Total computation time: %.2f minutes\n",
              diagnostics$total_search_time))
  cat("\n")

  # Create diagnostic plots if requested
  plots <- NULL
  if (create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- create_bootstrap_diagnostic_plots(results, H_estimates, Hc_estimates)

    # Add combined plot if patchwork is available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      plots$combined <- (plots$H_distribution | plots$Hc_distribution) +
        patchwork::plot_annotation(
          title = "Bootstrap Distributions",
          subtitle = sprintf("%d iterations", nrow(boot_results$results)),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold")
                        )
        )
    }
  }
  list(
    table = formatted_table,
    diagnostics = diagnostics,
    plots = plots
  )
}


#' Create Bootstrap Diagnostic Plots
#'
#' Generates plots showing bootstrap distribution and bias correction
#'
#' @param results Data.table of bootstrap results
#' @param H_estimates Data.table of H subgroup estimates
#' @param Hc_estimates Data.table of Hc subgroup estimates
#'
#' @return List of ggplot objects
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline
#' @keywords internal

create_bootstrap_diagnostic_plots <- function(results, H_estimates, Hc_estimates) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }

  # Plot 1: Bootstrap distribution of bias-corrected estimates H
  p1 <- ggplot2::ggplot(results, ggplot2::aes(x = H_biasadj_2)) +
    ggplot2::geom_histogram(bins = 30, fill = "#4472C4", alpha = 0.7, color = "white") +
    ggplot2::geom_vline(xintercept = c(log(H_estimates$H2),log(H_estimates$H0)), color = c("green","red"),
                        linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "Bootstrap Distribution of Bias-Corrected
      log(hazard ratio): Subgroup H",
      subtitle = "Red line shows final bias-corrected estimate",
      x = "Log Hazard Ratio (bias-corrected)",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
    )

  # Plot 2: Bootstrap distribution of bias-corrected estimates Hc
  p2 <- ggplot2::ggplot(results, ggplot2::aes(x = Hc_biasadj_2)) +
    ggplot2::geom_histogram(bins = 30, fill = "#4472C4", alpha = 0.7, color = "white") +
    ggplot2::geom_vline(xintercept = c(log(Hc_estimates$H2),log(Hc_estimates$H0)), color = c("green","red"),
                        linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "Bootstrap Distribution of Bias-Corrected
      log(hazard ratio): Subgroup Hc",
      subtitle = "Red line shows final bias-corrected estimate",
      x = "Log Hazard Ratio (bias-corrected)",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray40")
    )



  # Plot 3: Bias correction impact
  bias_data <- data.frame(
    iteration = 1:nrow(results),
    unadjusted = H_estimates$H0,
    adjusted = results$H_biasadj_2
  )

  p3 <- ggplot2::ggplot(results, ggplot2::aes(x = H_biasadj_1, y = H_biasadj_2)) +
    ggplot2::geom_point(alpha = 0.5, color = "#4472C4") +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Comparison of H Bias Correction Methods",
      x = "Method 1 (Simple optimism correction)",
      y = "Method 2 (Double bootstrap correction)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  p4 <- ggplot2::ggplot(results, ggplot2::aes(x = Hc_biasadj_1, y = Hc_biasadj_2)) +
    ggplot2::geom_point(alpha = 0.5, color = "#4472C4") +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Comparison of Hc Bias Correction Methods",
      x = "Method 1 (Simple optimism correction)",
      y = "Method 2 (Double bootstrap correction)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )



  list(
    H_distribution = p1,
    Hc_distribution = p2,
    H_methods_comparison = p3,
    Hc_methods_comparison = p4
  )
}


#' Simple Caption Generator for Bootstrap Table
#'
#' Generates an informative caption based on bootstrap results
#'
#' @param nb_boots Integer. Number of bootstrap iterations
#' @param boot_success_rate Numeric. Proportion of successful bootstraps
#' @param est.scale Character. "hr" or "1/hr"
#'
#' @return Character string with formatted caption
#' @export

generate_bootstrap_caption <- function(nb_boots, boot_success_rate, est.scale = "hr") {

  effect_name <- ifelse(est.scale == "hr", "hazard ratio", "inverse hazard ratio (1/HR)")

  caption <- sprintf(
    "Treatment effect estimates by identified subgroups. The unadjusted %s (HR) is from standard Cox regression with robust standard errors. The bias-corrected estimate (HR*) adjusts for optimism in subgroup selection using bootstrap methods (%d iterations, %.1f%% successful). Values less than 1.0 indicate benefit from treatment; bias correction typically shifts estimates toward the null, reflecting the optimism inherent in data-driven subgroup identification.",
    effect_name,
    nb_boots,
    boot_success_rate * 100
  )

  return(caption)
}


#' Summarize Bootstrap Event Counts
#'
#' Provides summary statistics for event counts across bootstrap iterations
#'
#' @param boot_results Data.table from bootstrap_results()
#' @param threshold Integer. Minimum event threshold (default: 5)
#' @return List with summary statistics
#' @export

summarize_bootstrap_events <- function(boot_results, threshold = 5) {

  results <- boot_results$results
  nb_boots <- nrow(results)

  # Check low events for ORIGINAL subgroup on BOOTSTRAP samples
  low_H_0 <- sum(results$events_H_0 < threshold, na.rm = TRUE)
  low_H_1 <- sum(results$events_H_1 < threshold, na.rm = TRUE)
  low_H_either <- sum(results$events_H_0 < threshold | results$events_H_1 < threshold, na.rm = TRUE)

  low_Hc_0 <- sum(results$events_Hc_0 < threshold, na.rm = TRUE)
  low_Hc_1 <- sum(results$events_Hc_1 < threshold, na.rm = TRUE)
  low_Hc_either <- sum(results$events_Hc_0 < threshold | results$events_Hc_1 < threshold, na.rm = TRUE)

  # Check low events for NEW subgroup on ORIGINAL data
  low_Hstar_0 <- sum(results$events_Hstar_0 < threshold, na.rm = TRUE)
  low_Hstar_1 <- sum(results$events_Hstar_1 < threshold, na.rm = TRUE)
  low_Hstar_either <- sum(results$events_Hstar_0 < threshold | results$events_Hstar_1 < threshold, na.rm = TRUE)

  low_Hcstar_0 <- sum(results$events_Hcstar_0 < threshold, na.rm = TRUE)
  low_Hcstar_1 <- sum(results$events_Hcstar_1 < threshold, na.rm = TRUE)
  low_Hcstar_either <- sum(results$events_Hcstar_0 < threshold | results$events_Hcstar_1 < threshold, na.rm = TRUE)

  # Summary statistics
  cat("\n=== Bootstrap Event Count Summary ===\n")
  cat(sprintf("Total bootstrap iterations: %d\n", nb_boots))
  cat(sprintf("Event threshold: <%d events\n\n", threshold))

  cat("ORIGINAL Subgroup H on BOOTSTRAP samples:\n")
  cat(sprintf("  Control arm <%d events: %d (%.1f%%)\n",
              threshold, low_H_0, 100*low_H_0/nb_boots))
  cat(sprintf("  Treatment arm <%d events: %d (%.1f%%)\n",
              threshold, low_H_1, 100*low_H_1/nb_boots))
  cat(sprintf("  Either arm <%d events: %d (%.1f%%)\n\n",
              threshold, low_H_either, 100*low_H_either/nb_boots))

  cat("ORIGINAL Subgroup Hc on BOOTSTRAP samples:\n")
  cat(sprintf("  Control arm <%d events: %d (%.1f%%)\n",
              threshold, low_Hc_0, 100*low_Hc_0/nb_boots))
  cat(sprintf("  Treatment arm <%d events: %d (%.1f%%)\n",
              threshold, low_Hc_1, 100*low_Hc_1/nb_boots))
  cat(sprintf("  Either arm <%d events: %d (%.1f%%)\n\n",
              threshold, low_Hc_either, 100*low_Hc_either/nb_boots))

  # Count successful bootstraps (found new subgroup)
  n_successful <- sum(!is.na(results$events_Hstar_0))

  if (n_successful > 0) {
    cat(sprintf("NEW Subgroups found: %d (%.1f%%)\n\n",
                n_successful, 100*n_successful/nb_boots))

    cat("NEW Subgroup H* on ORIGINAL data:\n")
    cat(sprintf("  Control arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hstar_0, 100*low_Hstar_0/n_successful))
    cat(sprintf("  Treatment arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hstar_1, 100*low_Hstar_1/n_successful))
    cat(sprintf("  Either arm <%d events: %d (%.1f%% of successful)\n\n",
                threshold, low_Hstar_either, 100*low_Hstar_either/n_successful))

    cat("NEW Subgroup Hc* on ORIGINAL data:\n")
    cat(sprintf("  Control arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_0, 100*low_Hcstar_0/n_successful))
    cat(sprintf("  Treatment arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_1, 100*low_Hcstar_1/n_successful))
    cat(sprintf("  Either arm <%d events: %d (%.1f%% of successful)\n",
                threshold, low_Hcstar_either, 100*low_Hcstar_either/n_successful))
  } else {
    cat("No new subgroups found in bootstrap samples\n")
  }

  cat("\n")

  # Return invisible list with all counts
  invisible(list(
    threshold = threshold,
    nb_boots = nb_boots,
    n_successful = n_successful,
    original_H = list(
      low_control = low_H_0,
      low_treat = low_H_1,
      low_either = low_H_either
    ),
    original_Hc = list(
      low_control = low_Hc_0,
      low_treat = low_Hc_1,
      low_either = low_Hc_either
    ),
    new_Hstar = list(
      low_control = low_Hstar_0,
      low_treat = low_Hstar_1,
      low_either = low_Hstar_either
    ),
    new_Hcstar = list(
      low_control = low_Hcstar_0,
      low_treat = low_Hcstar_1,
      low_either = low_Hcstar_either
    )
  ))
}
