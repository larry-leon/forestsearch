# =============================================================================
# L_eff Calibration I/O Helpers
# =============================================================================
#
# Standardized save/load for L_eff calibration results.
# Calibration documents save to quarto/_data/.
# Selection and validation documents load from the same location.
#
# =============================================================================


#' Save L_eff Calibration Results
#'
#' Writes a standardized \code{.rds} bundle containing the fitted
#' L_eff parameters and all supporting data.  Called at the end of
#' calibration documents (e.g., \code{calibration_binary_leff.qmd}).
#'
#' @param C Numeric. Fitted scale parameter.
#' @param alpha Numeric. Fitted power-law exponent.
#' @param n_min Integer. Reference minimum subgroup size.
#' @param outcome_type Character. One of "binary", "survival",
#'   "count", "continuous".
#' @param cal_data Data frame with columns N, sim_fpr, P1, L_eff
#'   (the raw calibration points).
#' @param dgm_description Character. Free-text description of the
#'   DGM used for calibration.
#' @param n_sims_per_N Integer. Simulations per sample size.
#' @param output_dir Character. Directory for .rds files.
#'   Default: \code{"_data"}.
#' @param extra List. Any additional items to store.
#'
#' @return Invisible path to the saved file.
#'
#' @examples
#' \dontrun{
#' save_leff_calibration(
#'   C = 0.220, alpha = 1.298, n_min = 60,
#'   outcome_type = "binary",
#'   cal_data = cal_df,
#'   dgm_description = "Binary DGM, 4 confounders (bm1, bm2, age, ecog)",
#'   n_sims_per_N = 5000
#' )
#' }
#'
#' @export
save_leff_calibration <- function(C, alpha, n_min,
                                   outcome_type,
                                   cal_data = NULL,
                                   dgm_description = "",
                                   n_sims_per_N = NA_integer_,
                                   output_dir = NULL,
                                   extra = list()) {

  # ── Input validation ───────────────────────────────────────────────────────

  stopifnot(
    "C must be a positive scalar" = is.numeric(C) && length(C) == 1 && C > 0,
    "alpha must be a finite numeric scalar" =
      is.numeric(alpha) && length(alpha) == 1 && is.finite(alpha),
    "n_min must be a positive integer" =
      is.numeric(n_min) && length(n_min) == 1 && n_min > 0,
    "outcome_type must be a non-empty character string" =
      is.character(outcome_type) && nzchar(outcome_type)
  )

  if (is.null(output_dir)) {
    output_dir <- "_data"
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  bundle <- list(
    # Core parameters (what downstream documents need)
    C = C,
    alpha = alpha,
    n_min = n_min,
    outcome_type = outcome_type,

    # Provenance
    dgm_description = dgm_description,
    n_sims_per_N = n_sims_per_N,
    timestamp = Sys.time(),
    R_version = R.version.string,
    forestsearch_version = tryCatch(
      as.character(utils::packageVersion("forestsearch")),
      error = function(e) "unknown"
    ),

    # Raw data for re-analysis / plotting
    cal_data = cal_data,

    # User extras
    extra = extra
  )

  filename <- sprintf("calibration_%s_leff.rds", outcome_type)
  path <- file.path(output_dir, filename)
  saveRDS(bundle, path)

  message(sprintf("L_eff calibration saved: %s", path))
  message(sprintf("  outcome: %s | C = %.4f | alpha = %.4f | n_min = %d",
                  outcome_type, C, alpha, n_min))
  message(sprintf("  Size: %.1f KB", file.size(path) / 1024))

  invisible(path)
}


#' Load L_eff Calibration Results
#'
#' Reads a calibration bundle saved by \code{save_leff_calibration()}.
#' Prints a summary and returns the full bundle.
#'
#' @param outcome_type Character. One of "binary", "survival",
#'   "count", "continuous".
#' @param output_dir Character. Directory containing .rds files.
#'   Default: \code{"_data"}.
#' @param verbose Logical. Print summary on load. Default: TRUE.
#'
#' @return A list with components C, alpha, n_min, outcome_type,
#'   cal_data, dgm_description, timestamp, etc.
#'
#' @examples
#' \dontrun{
#' cal <- load_leff_calibration("binary")
#' C_prior <- cal$C
#' alpha_prior <- cal$alpha
#' }
#'
#' @export
load_leff_calibration <- function(outcome_type,
                                   output_dir = NULL,
                                   verbose = TRUE) {

  if (is.null(output_dir)) {
    output_dir <- "_data"
  }

  filename <- sprintf("calibration_%s_leff.rds", outcome_type)
  path <- file.path(output_dir, filename)

  if (!file.exists(path)) {
    stop(sprintf(
      "Calibration file not found: %s\n  Run the calibration document first, or check output_dir.",
      path))
  }

  bundle <- readRDS(path)

  if (verbose) {
    message(sprintf("=== L_eff Calibration: %s ===", outcome_type))
    message(sprintf("  C = %.4f, alpha = %.4f, n_min = %d",
                    bundle$C, bundle$alpha, bundle$n_min))
    message(sprintf("  DGM: %s", bundle$dgm_description))
    message(sprintf("  Sims per N: %s",
                    as.character(bundle$n_sims_per_N)))
    message(sprintf("  Calibrated: %s",
                    format(bundle$timestamp, "%Y-%m-%d %H:%M")))
    if (!is.null(bundle$cal_data)) {
      message(sprintf("  N values: %s",
                      paste(bundle$cal_data$N, collapse = ", ")))
    }
  }

  bundle
}


#' Get L_eff at a Given Sample Size
#'
#' Convenience function: loads calibration and computes
#' \eqn{L_{\text{eff}}(N) = C \cdot (N / n_{\min})^\alpha}.
#'
#' @param N Integer. Sample size.
#' @param outcome_type Character. Outcome type to load.
#' @param output_dir Character. Default: NULL (auto).
#'
#' @return Numeric L_eff value.
#'
#' @export
get_leff <- function(N, outcome_type, output_dir = NULL) {
  stopifnot(
    "N must be a positive numeric value" =
      is.numeric(N) && all(N > 0)
  )
  cal <- load_leff_calibration(outcome_type, output_dir,
                                verbose = FALSE)
  cal$C * (N / cal$n_min)^cal$alpha
}


#' Find the Most Recent .rds File Matching a Pattern
#'
#' Scans \code{output_dir} for files matching the given prefix
#' pattern and returns the path to the most recently modified one.
#' Date-stamped filenames (e.g., \code{results_20260402.rds}) sort
#' naturally, but this function uses filesystem modification time
#' for robustness.
#'
#' @param prefix Character. Filename prefix to match
#'   (e.g., \code{"selection_survival_frontier"}).
#' @param output_dir Character. Directory to scan.
#'   Default: \code{"_data"}.
#' @param verbose Logical. Print the matched file. Default: TRUE.
#'
#' @return Full path to the most recent matching file, or NULL
#'   if no match is found.
#'
#' @examples
#' \dontrun{
#' path <- latest_rds("selection_survival_frontier")
#' res  <- readRDS(path)
#'
#' path <- latest_rds("calibration_binary_leff_grid")
#' }
#'
#' @export
latest_rds <- function(prefix, output_dir = "_data",
                        verbose = TRUE) {
  if (!dir.exists(output_dir)) {
    if (verbose) {
      message(sprintf("Directory '%s' does not exist.", output_dir))
    }
    return(NULL)
  }

  pattern <- sprintf("^%s.*\\.rds$", prefix)
  files <- list.files(output_dir, pattern = pattern,
                       full.names = TRUE)

  if (length(files) == 0) {
    if (verbose) {
      message(sprintf("No .rds files matching '%s' in %s/",
                      prefix, output_dir))
    }
    return(NULL)
  }

  # Prefer date-stamp sort when filenames contain _YYYYMMDD;

  # fall back to mtime (which can be unreliable on NFS / cloud mounts).
  date_stamps <- regmatches(
    basename(files),
    regexpr("\\d{8}", basename(files))
  )
  has_dates <- length(date_stamps) == length(files) &&
    all(nchar(date_stamps) == 8)

  if (has_dates) {
    idx <- order(date_stamps, decreasing = TRUE)
    sort_label <- "date stamp"
  } else {
    mtimes <- file.mtime(files)
    idx <- order(mtimes, decreasing = TRUE)
    sort_label <- "modification time"
  }
  best <- files[idx[1]]

  if (verbose) {
    message(sprintf(
      "Found %d file(s) matching '%s' (sorted by %s):",
      length(files), prefix, sort_label
    ))
    for (i in seq_along(idx)) {
      flag <- if (i == 1) " [latest]" else ""
      mt <- format(file.mtime(files[idx[i]]), "%Y-%m-%d %H:%M")
      message(sprintf("  %s  (%s)%s",
                      basename(files[idx[i]]), mt, flag))
    }
  }

  best
}
