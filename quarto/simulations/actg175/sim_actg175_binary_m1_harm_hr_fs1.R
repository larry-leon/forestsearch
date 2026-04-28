# ============================================================================
# sim_actg175_binary_m1_harm_hrMaxSG_fs1.R
#
# Simulation core for the configuration:
#   study     : ACTG175
#   outcome   : binary (adverse: failure to improve)
#   DGM       : m1   — calibrated GLM interaction (target OR(Q) = 2.0)
#   objective : harm
#   sg_focus  : hrMaxSG
#   FS bundle : fs1
#
# Companion files (same stem):
#   actg175_binary_m1_harm_hrMaxSG_fs1.qmd  — pkgdown article (vignette)
#   actg175_binary_m1_harm_hrMaxSG_fs1.R    — background runner script
#   actg175_binary_m1_harm_hrMaxSG_fs1.rds  — saved bundle (output)
#
# This file defines a single function, run_sim_actg175_binary_m1_harm_hrMaxSG_fs1(),
# which orchestrates the H1 and H0 foreach loops and returns a structured bundle.
# Caller is responsible for:
#   - library(forestsearch), library(doFuture), library(foreach)
#   - future::plan(...) configuration
#   - data preparation and DGM calibration (results passed in as arguments)
#   - persisting the returned bundle (saveRDS)
#
# The function body MIRRORS the H1/H0 simulation chunks in the companion .qmd.
# Updates to the foreach orchestration must be applied in both places.
# ============================================================================

#' Run ACTG175 binary harm hr fs1 simulations
#'
#' Pure simulation core for configuration
#' \code{actg175_binary_m1_harm_hrMaxSG_fs1}. Runs the H1 (calibrated HTE) and
#' H0 (null) foreach loops via \pkg{doFuture} and returns a self-describing
#' bundle containing metadata, the input DGMs, and the per-replicate results
#' from \code{\link[forestsearch]{collect_results}}.
#'
#' The function does NOT call \code{future::plan()}, write any files, or
#' configure logging. Those are caller responsibilities. Inside each worker,
#' \code{future::plan("sequential")} is set to prevent nested parallelism
#' (matches the vignette pattern).
#'
#' @param dgm_calibrated DGM object from
#'   \code{\link[forestsearch]{calibrate_glm_interaction}} for H1.
#' @param dgm_null DGM object from
#'   \code{\link[forestsearch]{generate_glm_dgm}} with \code{model = "null"} for H0.
#' @param n_alt Integer. Number of H1 replicates. May be 0 if \code{run_alt = FALSE}.
#' @param n_null Integer. Number of H0 replicates. May be 0 if \code{run_null = FALSE}.
#' @param sim_n Integer. Per-trial sample size.
#' @param confounders Character vector of analysis covariates passed to
#'   \code{run_simulation_analysis(confounders_base = ...)}.
#' @param fs_params List of FS parameters passed through to
#'   \code{\link[forestsearch]{forestsearch}} via
#'   \code{\link[forestsearch]{run_simulation_analysis}}. Should be fully
#'   populated, including \code{conf_force} (which depends on data-prep
#'   factor variables and is the caller's responsibility to inject).
#' @param grf_params List of GRF parameters.
#' @param diag_keep Character vector of per-replicate diagnostics to retain.
#'   See \code{\link[forestsearch]{run_simulation_analysis}} for valid values.
#' @param diag_keep_first_n Integer. Number of replicates for which to keep
#'   diagnostics. Set to \code{0L} to disable.
#' @param run_alt Logical. Run the H1 loop. Default \code{TRUE}.
#' @param run_null Logical. Run the H0 loop. Default \code{TRUE}.
#' @param verbose Logical. Print progress and timing messages from this
#'   function (does NOT affect per-replicate verbosity, which is forced
#'   to \code{FALSE} inside the loop).
#' @param config_id Character identifier embedded in the bundle metadata.
#'   Defaults to \code{"actg175_binary_m1_harm_hrMaxSG_fs1"}.
#'
#' @return A list with three top-level elements:
#'   \describe{
#'     \item{\code{meta}}{Named list: \code{config_id}, \code{timestamp},
#'       sample sizes, run flags, \code{confounders}, \code{fs_params},
#'       \code{grf_params}, diagnostic settings, elapsed times in seconds,
#'       failure counts, package version, and \code{sessionInfo()}.}
#'     \item{\code{dgm}}{Named list with elements \code{calibrated} and
#'       \code{null}, holding the input DGMs (preserved for downstream OC
#'       tables and consistency checks). Each element is \code{NULL} if the
#'       corresponding loop was skipped.}
#'     \item{\code{results}}{Named list with elements \code{alt} and
#'       \code{null}, each a data.frame returned by
#'       \code{\link[forestsearch]{collect_results}} (with
#'       \code{n_failed} attribute). Either may be \code{NULL} if the
#'       corresponding loop was skipped.}
#'   }
#'
#' @section Reproducibility:
#'   Both foreach loops use \code{.options.future = list(seed = TRUE)} for
#'   L'Ecuyer-CMRG streams. Set the global seed via \code{set.seed()} before
#'   calling this function for reproducible streams.
#'
#' @section Plan setup:
#'   The caller must call \code{future::plan(...)} before invoking this
#'   function. On Linux, \code{multicore} or \code{multisession} both work.
#'   On macOS, prefer \code{multisession} (multicore silently falls back
#'   in RStudio).
#'
#' @seealso
#'   The companion vignette \code{actg175_binary_m1_harm_hrMaxSG_fs1.qmd}
#'   demonstrates the same loops inline, with mode logic for
#'   \code{import} / \code{demo} / \code{full} runs.
run_sim_actg175_binary_m1_harm_hr_fs1 <- function(
    dgm_calibrated,
    dgm_null,
    n_alt,
    n_null,
    sim_n,
    confounders,
    fs_params,
    grf_params,
    diag_keep         = c("fs_full", "candidate_table"),
    diag_keep_first_n = 5L,
    run_alt           = TRUE,
    run_null          = TRUE,
    verbose           = FALSE,
    config_id         = "actg175_binary_m1_harm_hr_fs1"
) {

  # ── Input validation ──────────────────────────────────────────────────
  stopifnot(
    is.numeric(n_alt),  length(n_alt)  == 1L, n_alt  >= 0,
    is.numeric(n_null), length(n_null) == 1L, n_null >= 0,
    is.numeric(sim_n),  length(sim_n)  == 1L, sim_n  >  0,
    is.character(confounders), length(confounders) > 0L,
    is.list(fs_params),
    is.list(grf_params),
    is.character(diag_keep),
    is.numeric(diag_keep_first_n), length(diag_keep_first_n) == 1L
  )
  if (!run_alt && !run_null) {
    stop("Both run_alt and run_null are FALSE; nothing to do.", call. = FALSE)
  }
  if (run_alt && is.null(dgm_calibrated)) {
    stop("run_alt = TRUE but dgm_calibrated is NULL.", call. = FALSE)
  }
  if (run_null && is.null(dgm_null)) {
    stop("run_null = TRUE but dgm_null is NULL.", call. = FALSE)
  }

  n_alt             <- as.integer(n_alt)
  n_null            <- as.integer(n_null)
  sim_n             <- as.integer(sim_n)
  diag_keep_first_n <- as.integer(diag_keep_first_n)

  worker_packages <- c("forestsearch", "survival", "data.table", "future")

  # ── H1 loop (calibrated HTE) ──────────────────────────────────────────
  results_alt <- NULL
  t_alt       <- NA_real_
  if (run_alt && n_alt > 0L) {
    if (verbose) {
      cat(sprintf("[%s] Running %d H1 replicates...\n", config_id, n_alt))
    }
    t0 <- proc.time()[3]
    raw_alt <- foreach::foreach(
      sim = seq_len(n_alt),
      .errorhandling = "pass",
      .options.future = list(packages = worker_packages, seed = TRUE)
    ) %dofuture% {
      future::plan("sequential")
      run_simulation_analysis(
        sim_id           = sim,
        dgm              = dgm_calibrated,
        n_sample         = sim_n,
        confounders_base = confounders,
        run_fs           = TRUE,
        run_fs_grf       = FALSE,
        run_grf          = TRUE,
        fs_params        = fs_params,
        grf_params       = grf_params,
        keep             = diag_keep,
        keep_first_n     = diag_keep_first_n,
        verbose          = FALSE
      )
    }
    t_alt <- proc.time()[3] - t0
    results_alt <- collect_results(
      raw_alt, "H1 alt",
      keep_diagnostics = diag_keep_first_n > 0L
    )
    if (verbose) {
      cat(sprintf("[%s] H1 done: %.1f min, %d rows, %d failed\n",
                  config_id, t_alt / 60,
                  nrow(results_alt),
                  attr(results_alt, "n_failed")))
    }
  }

  # ── H0 loop (null — homogeneous effect) ───────────────────────────────
  results_null <- NULL
  t_null       <- NA_real_
  if (run_null && n_null > 0L) {
    if (verbose) {
      cat(sprintf("[%s] Running %d H0 replicates...\n", config_id, n_null))
    }
    t0 <- proc.time()[3]
    raw_null <- foreach::foreach(
      sim = seq_len(n_null),
      .errorhandling = "pass",
      .options.future = list(packages = worker_packages, seed = TRUE)
    ) %dofuture% {
      future::plan("sequential")
      run_simulation_analysis(
        sim_id           = sim,
        dgm              = dgm_null,
        n_sample         = sim_n,
        confounders_base = confounders,
        run_fs           = TRUE,
        run_fs_grf       = FALSE,
        run_grf          = TRUE,
        fs_params        = fs_params,
        grf_params       = grf_params,
        verbose          = FALSE
      )
    }
    t_null <- proc.time()[3] - t0
    results_null <- collect_results(raw_null, "H0 null")
    if (verbose) {
      cat(sprintf("[%s] H0 done: %.1f min, %d rows, %d failed\n",
                  config_id, t_null / 60,
                  nrow(results_null),
                  attr(results_null, "n_failed")))
    }
  }

  # ── Build bundle ──────────────────────────────────────────────────────
  list(
    meta = list(
      config_id         = config_id,
      timestamp         = Sys.time(),
      n_alt             = n_alt,
      n_null            = n_null,
      sim_n             = sim_n,
      run_alt           = run_alt,
      run_null          = run_null,
      confounders       = confounders,
      fs_params         = fs_params,
      grf_params        = grf_params,
      diag_keep         = diag_keep,
      diag_keep_first_n = diag_keep_first_n,
      t_alt_sec         = unname(t_alt),
      t_null_sec        = unname(t_null),
      n_failed_alt      = if (!is.null(results_alt))
                            attr(results_alt, "n_failed")  else NA_integer_,
      n_failed_null     = if (!is.null(results_null))
                            attr(results_null, "n_failed") else NA_integer_,
      package_version   = as.character(utils::packageVersion("forestsearch")),
      sessionInfo       = utils::sessionInfo()
    ),
    dgm = list(
      calibrated = if (run_alt)  dgm_calibrated else NULL,
      null       = if (run_null) dgm_null       else NULL
    ),
    results = list(
      alt  = results_alt,
      null = results_null
    )
  )
}
